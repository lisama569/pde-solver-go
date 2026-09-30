package solver

import "sync"

// rowBounds 返回把 ny 行均分为 workers 个行子域的分界 [b0,...,bw]：
// band b 负责 [bounds[b], bounds[b+1])。
func rowBounds(ny, workers int) []int {
	bounds := make([]int, workers+1)
	for b := 0; b <= workers; b++ {
		bounds[b] = b * ny / workers
	}
	return bounds
}

// stepRowsParallel 按行子域并发执行单步更新；workers=1 时退化为串行。
// 由于双缓冲（读 Cur、写 Next），各子域只写自己的行区间，
// 无数据竞争、无重复更新；每个格点的浮点运算序列与串行完全相同。
func (s *Stepper) stepRowsParallel(workers int) {
	if workers < 1 {
		workers = 1
	}
	bounds := rowBounds(s.Cur.NY, workers)
	var wg sync.WaitGroup
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
}

// RunParallel 按行把网格划分为 workers 个子域并发推进 steps 步。
//
// 由于采用双缓冲（读 Cur、写 Next），各子域只写自己的行区间，
// 无数据竞争、无重复更新；每个格点的浮点运算序列与串行 Advance 逐位一致，
// 因此结果与串行基线逐位一致，与并发度无关。
func RunParallel(cfg *Config, initial *Grid, workers, steps int) (*Grid, error) {
	if workers < 1 {
		workers = 1
	}
	s, err := NewStepper(cfg, initial)
	if err != nil {
		return nil, err
	}
	for k := 0; k < steps; k++ {
		s.stepRowsParallel(workers)
		if err := s.commit(); err != nil {
			return nil, err
		}
	}
	return s.Cur, nil
}
