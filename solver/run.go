package solver

// RunPhases 按 cfg.Plan（或旧的单套配置）从步号 startStep 起推进 steps 步，
// workers 个子域并行，返回最终网格与逐阶段收支对账报告。
//
// 语义：
//   - 每步推进前按当前步号激活对应阶段（阶段切换只发生在步边界）；
//   - 报告按"本次运行实际执行的步"统计：断点续算时，被截断阶段的报告
//     只覆盖本次执行的步区间，总量变化与源/边界账目仍严格闭合；
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

	// 先激活起始步所属阶段，使首份报告的源数目与阶段配置一致
	s.applyPhaseForStep(startStep)

	var reports []PhaseReport
	var cur *PhaseReport
	openReport := func(step int) {
		idx := -1
		if cfg.Plan != nil {
			idx = cfg.Plan.PhaseAt(step)
		}
		reports = append(reports, PhaseReport{
			PhaseIndex:   idx,
			FromStep:     step,
			SourceInject: make([]float64, len(s.srcs)),
			MassStart:    TotalMass(s.Cur, cfg.DX, cfg.DY),
		})
		cur = &reports[len(reports)-1]
	}
	closeReport := func(step int) {
		cur.ToStep = step
		cur.MassEnd = TotalMass(s.Cur, cfg.DX, cfg.DY)
	}

	openReport(startStep)
	for k := 0; k < steps; k++ {
		if cfg.Plan != nil && k > 0 && cfg.Plan.PhaseAt(s.Step) != s.phaseIdx {
			// 跨阶段边界：先用旧阶段结束时的场结算上一阶段，
			// 再激活新阶段（写入新 Dirichlet 值），最后开启新报告。
			// 切换瞬间边界值跳变不计入任一阶段的收支账目。
			closeReport(s.Step)
			s.applyPhaseForStep(s.Step)
			openReport(s.Step)
		} else {
			s.applyPhaseForStep(s.Step)
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
