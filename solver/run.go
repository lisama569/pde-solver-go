package solver

// RunPhases 按 cfg.Plan（或旧的单套配置）从步号 startStep 起推进 steps 步，
// workers 个子域并行，返回最终网格与逐阶段收支对账报告。
//
// 语义：
//   - 每步推进前按当前步号激活对应阶段（阶段切换只发生在步边界）；
//   - 报告按"本次运行实际执行的步"统计：断点续算时，被截断阶段的报告
//     只覆盖本次执行的步区间，总量变化与源/边界账目仍严格闭合；
//   - 阶段激活时写入 Dirichlet 固定值造成的总量跳变（重置量）计入该报告
//     对应边界的 Boundary，并单列在 Reset 中——切换那一下的贡献看得见、
//     对得上，各阶段 ΔMass 之和等于任务真实总量变化；
//   - 热源注入只统计实际施加到非固定格的部分（贴边/贴角热源区不误报）；
//   - 无阶段计划时返回一份 PhaseIndex=-1、覆盖整个运行区间的报告。
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

	var reports []PhaseReport
	var cur *PhaseReport
	// openReport 在阶段激活**之前**开启报告：MassStart 取尚未写入本阶段
	// Dirichlet 值的场，使随后发生的重置跳变落在本报告账内。
	// 源条目数直接取阶段配置，不依赖激活顺序。
	openReport := func(step int) {
		idx := -1
		nsrc := len(s.srcs) // 无计划：NewStepper 已从 cfg.Src 装好
		if cfg.Plan != nil {
			idx = cfg.Plan.PhaseAt(step)
			nsrc = len(cfg.Plan.Phases[idx].Sources)
		}
		reports = append(reports, PhaseReport{
			PhaseIndex:   idx,
			FromStep:     step,
			SourceInject: make([]float64, nsrc),
			MassStart:    TotalMass(s.Cur, cfg.DX, cfg.DY),
		})
		cur = &reports[len(reports)-1]
	}
	closeReport := func(step int) {
		cur.ToStep = step
		cur.MassEnd = TotalMass(s.Cur, cfg.DX, cfg.DY)
	}
	// activate 激活当前步所属阶段，并把 Dirichlet 重置量计入当前报告。
	activate := func() {
		reset := s.applyPhaseForStep(s.Step)
		cur.Boundary.Left += reset.Left
		cur.Boundary.Right += reset.Right
		cur.Boundary.Top += reset.Top
		cur.Boundary.Bottom += reset.Bottom
		cur.Reset.Left += reset.Left
		cur.Reset.Right += reset.Right
		cur.Reset.Top += reset.Top
		cur.Reset.Bottom += reset.Bottom
	}

	openReport(startStep)
	activate()
	for k := 0; k < steps; k++ {
		if cfg.Plan != nil && k > 0 && cfg.Plan.PhaseAt(s.Step) != s.phaseIdx {
			// 跨阶段边界：先用旧阶段结束时的场结算上一阶段；
			// 再开新报告（MassStart = 重置前的场），然后激活新阶段，
			// 切换瞬间边界值跳变作为重置量计入新阶段的边界账目。
			closeReport(s.Step)
			openReport(s.Step)
			activate()
		} else {
			activate() // 同阶段幂等，重置量为 0
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
