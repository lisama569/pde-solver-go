package solver

import "sync"

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
	ny := cfg.NY
	// 行分块边界：band b 负责 [b*ny/w, (b+1)*ny/w)
	bounds := make([]int, workers+1)
	for b := 0; b <= workers; b++ {
		bounds[b] = b * ny / workers
	}
	var wg sync.WaitGroup
	for k := 0; k < steps; k++ {
		for b := 0; b < workers; b++ {
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
		if err := s.commit(); err != nil {
			return nil, err
		}
	}
	return s.Cur, nil
}
