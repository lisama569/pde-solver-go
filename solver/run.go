package solver

// RunPhases 按 cfg.Plan（或旧的单套配置）从步号 startStep 起推进 steps 步，
// workers 个子域并行，返回最终网格与逐阶段收支对账报告。
//
// 语义：
//   - 每步推进前按当前步号激活对应阶段（阶段切换只发生在步边界）；
//   - 每份报告的 MassStart 在任何边界重置之前测量；阶段激活瞬间 Dirichlet
//     写入造成的逐边总量跳变记入 BoundaryReset，因此跨阶段改边界值的
//     贡献"看得见、对得上"，且整任务各报告 ΔMass 之和恰好等于
//     从头到尾的真实总量变化（MassEnd 末 − 初始 MassStart 首）；
//   - 贴 Dirichlet 边/角的热源只统计实际施加到非固定格上的注入量；
//   - 报告按"本次运行实际执行的步"统计：断点续算时，被截断阶段的报告
//     只覆盖本次执行的步区间，恢复时重新激活当前阶段是幂等的（无跳变）；
//   - 无阶段计划时返回一份 PhaseIndex=-1、覆盖整个运行区间的报告，
//     BoundaryReset 恒为零（旧语义：不主动写入 Dirichlet 值）。
//
// 并行一致性：与 RunParallel 相同的双缓冲行分块，结果与 workers 无关。
func RunPhases(cfg *Config, initial *Grid, workers, startStep, steps int) (*Grid, []PhaseReport, error) {
	if workers < 1 {
		workers = 1
	}
	s, err := NewStepper(cfg, initial)
	if err != nil {
		return nil, nil, err
	}
	s.Step = startStep
	bounds := rowBounds(cfg.NY, workers)

	// nSourcesAt 返回覆盖 step 的阶段配置中的源数目（与激活后的 s.srcs 等长）。
	nSourcesAt := func(step int) int {
		if cfg.Plan != nil {
			if idx := cfg.Plan.PhaseAt(step); idx >= 0 {
				return len(cfg.Plan.Phases[idx].Sources)
			}
			return 0
		}
		if cfg.Src != nil {
			return 1
		}
		return 0
	}

	var reports []PhaseReport
	var cur *PhaseReport
	// openReport 必须在该阶段 Dirichlet 写入之前调用，使 MassStart 反映
	// 重置前的场；activate 随后把重置跳变记入新报告。
	openReport := func(step int) {
		idx := -1
		if cfg.Plan != nil {
			idx = cfg.Plan.PhaseAt(step)
		}
		reports = append(reports, PhaseReport{
			PhaseIndex:   idx,
			FromStep:     step,
			SourceInject: make([]float64, nSourcesAt(step)),
			MassStart:    TotalMass(s.Cur, cfg.DX, cfg.DY),
		})
		cur = &reports[len(reports)-1]
	}
	closeReport := func(step int) {
		cur.ToStep = step
		cur.MassEnd = TotalMass(s.Cur, cfg.DX, cfg.DY)
	}
	addReset := func(j BoundaryBudget) {
		cur.BoundaryReset.Left += j.Left
		cur.BoundaryReset.Right += j.Right
		cur.BoundaryReset.Top += j.Top
		cur.BoundaryReset.Bottom += j.Bottom
	}

	openReport(startStep)
	// 激活起始步所属阶段：任务首跑时把阶段0的 Dirichlet 值写入初始场，
	// 跳变计入首份报告；续算恢复时激活是幂等的，跳变为 0。
	addReset(s.applyPhaseForStep(startStep))

	for k := 0; k < steps; k++ {
		if cfg.Plan != nil && k > 0 && cfg.Plan.PhaseAt(s.Step) != s.phaseIdx {
			// 跨阶段边界：先用旧阶段结束时的场结算上一阶段；
			// 再以"重置前"的场开启新报告（MassStart），随后激活新阶段，
			// 边界值跳变作为 BoundaryReset 进入新报告。
			closeReport(s.Step)
			openReport(s.Step)
			addReset(s.applyPhaseForStep(s.Step))
		}
		// 用更新前的场 u^n 累计瞬时通量与源强（与更新算子严格对应）
		influx := s.BoundaryInflux()
		srcRates := s.SourceRates()
		if err := s.stepParallel(bounds); err != nil {
			return nil, nil, err
		}
		cur.Boundary.Left += influx.Left * cfg.DT
		cur.Boundary.Right += influx.Right * cfg.DT
		cur.Boundary.Top += influx.Top * cfg.DT
		cur.Boundary.Bottom += influx.Bottom * cfg.DT
		for si := range srcRates {
			cur.SourceInject[si] += srcRates[si] * cfg.DT
		}
	}
	closeReport(s.Step)
	return s.Cur, reports, nil
}
