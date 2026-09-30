package solver

import "sync"

// rowBounds 把 ny 行均分为 workers 个子域，返回 workers+1 个分界：
// 子域 b 负责行区间 [bounds[b], bounds[b+1])。
func rowBounds(ny, workers int) []int {
	bounds := make([]int, workers+1)
	for b := 0; b <= workers; b++ {
		bounds[b] = b * ny / workers
	}
	return bounds
}

// stepParallel 并发推进 s 一步：各子域只写自己的行区间。
// 调用前需已激活当前步的阶段配置（applyPhaseForStep）。
func (s *Stepper) stepParallel(bounds []int) error {
	var wg sync.WaitGroup
	for b := 0; b+1 < len(bounds); b++ {
		j0, j1 := bounds[b], bounds[b+1]
		if j0 == j1 {
			continue
		}
		wg.Add(1)
		go func(j0, j1 int) {
			defer wg.Done()
			s.stepRows(j0, j1)
		}(j0, j1)
	}
	wg.Wait()
	return s.commit()
}

// RunParallel 按行把网格划分为 workers 个子域并发推进 steps 步。
//
// 由于采用双缓冲（读 Cur、写 Next），各子域只写自己的行区间，
// 无数据竞争、无重复更新；每个格点的浮点运算序列与串行完全相同，
// 因此结果与串行 Advance 逐位一致，与并发度无关。
func RunParallel(cfg *Config, initial *Grid, workers, steps int) (*Grid, error) {
	if workers < 1 {
		workers = 1
	}
	s, err := NewStepper(cfg, initial)
	if err != nil {
		return nil, err
	}
	bounds := rowBounds(cfg.NY, workers)
	for k := 0; k < steps; k++ {
		s.applyPhaseForStep(s.Step)
		if err := s.stepParallel(bounds); err != nil {
			return nil, err
		}
	}
	return s.Cur, nil
}
