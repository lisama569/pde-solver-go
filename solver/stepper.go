package solver

import (
	"fmt"
	"math"
)

// Stepper 执行二维热传导方程的显式 FTCS 时间推进：
//
//	u^{n+1} = u^n + dt * (D * ∇²u^n + source)
//
// 采用双缓冲：所有新值只写入 Next，因此并行按行分块与串行结果逐位一致。
//
// 当 Config.Plan 非空时，Stepper 在每一步前按当前步号激活对应阶段的
// 边界与源项配置（见 applyPhaseForStep）；Plan 为空时始终使用
// Config.BC / Config.Src，行为与旧版本完全一致。
type Stepper struct {
	Cfg  *Config
	Cur  *Grid
	Next *Grid
	Step int // 已完成的步数

	bc       Boundaries // 当前生效的边界配置
	srcs     []Source   // 当前生效的源项列表
	phaseIdx int        // 当前生效的阶段索引（-1 = 尚未激活任何阶段）
}

// NewStepper 创建步进器；先校验配置合法性，再校验稳定性条件。
func NewStepper(cfg *Config, initial *Grid) (*Stepper, error) {
	if err := cfg.Validate(); err != nil {
		return nil, err
	}
	if err := cfg.CheckStability(); err != nil {
		return nil, err
	}
	if initial == nil || initial.NX != cfg.NX || initial.NY != cfg.NY {
		return nil, configErr(fmt.Sprintf("config: initial grid size mismatch, want %dx%d", cfg.NX, cfg.NY))
	}
	s := &Stepper{
		Cfg:      cfg,
		Cur:      initial.Clone(),
		Next:     NewGrid(cfg.NX, cfg.NY),
		bc:       cfg.BC,
		phaseIdx: -1,
	}
	if cfg.Src != nil {
		s.srcs = []Source{*cfg.Src}
	}
	return s, nil
}

// applyPhaseForStep 在推进步号 step 之前激活覆盖该步的阶段配置。
// 阶段切换只发生在步边界上，是步号的确定性函数，因此断点续算与
// 不间断运行的切换点完全一致。阶段激活时会把该阶段 Dirichlet 边界的
// 固定值写入 Cur 的边界格（幂等：同阶段重复激活结果相同）。
// 无阶段计划或步号超出计划范围时不做任何事。
//
// 返回本次激活造成的逐边质量跳变（Dirichlet 重置量，正 = 总量增加）；
// 未发生切换或无计划时返回零值。该跳变是切换瞬间的边界重置，不是扩散
// 通量，但必须计入收支账目（见 RunPhases），否则跨阶段改边值时
// 各阶段总量变化之和会对不上任务真实变化。
func (s *Stepper) applyPhaseForStep(step int) BoundaryBudget {
	if s.Cfg.Plan == nil {
		return BoundaryBudget{}
	}
	idx := s.Cfg.Plan.PhaseAt(step)
	if idx < 0 || idx == s.phaseIdx {
		return BoundaryBudget{}
	}
	ph := &s.Cfg.Plan.Phases[idx]
	s.bc = ph.BC
	s.srcs = ph.Sources
	s.phaseIdx = idx
	return s.imposeDirichlet()
}

// imposeDirichlet 把当前生效边界中 Dirichlet 边的固定值写入 Cur 边界格。
// 角点被两条 Dirichlet 边共享时，按 left→right→bottom→top 顺序后者覆盖。
//
// 返回逐边的质量跳变（写入前后差值 × dx·dy 的累计）。角点被两条边先后
// 覆盖时，跳变按写入次序拆分到两条边上（前者记 新值−旧值，后者记
// 最终值−新值），各边之和仍严格等于该格的实际质量变化。
func (s *Stepper) imposeDirichlet() BoundaryBudget {
	var b BoundaryBudget
	cell := s.Cfg.DX * s.Cfg.DY
	nx, ny := s.Cur.NX, s.Cur.NY
	if s.bc.Left.Type == Dirichlet {
		for j := 0; j < ny; j++ {
			old := s.Cur.At(0, j)
			s.Cur.Set(0, j, s.bc.Left.Value)
			b.Left += (s.bc.Left.Value - old) * cell
		}
	}
	if s.bc.Right.Type == Dirichlet {
		for j := 0; j < ny; j++ {
			old := s.Cur.At(nx-1, j)
			s.Cur.Set(nx-1, j, s.bc.Right.Value)
			b.Right += (s.bc.Right.Value - old) * cell
		}
	}
	if s.bc.Bottom.Type == Dirichlet {
		for i := 0; i < nx; i++ {
			old := s.Cur.At(i, 0)
			s.Cur.Set(i, 0, s.bc.Bottom.Value)
			b.Bottom += (s.bc.Bottom.Value - old) * cell
		}
	}
	if s.bc.Top.Type == Dirichlet {
		for i := 0; i < nx; i++ {
			old := s.Cur.At(i, ny-1)
			s.Cur.Set(i, ny-1, s.bc.Top.Value)
			b.Top += (s.bc.Top.Value - old) * cell
		}
	}
	return b
}

// neighbor 按边界条件取 (i+di, j+dj) 处的值。
// di,dj 只有一个非零（五点差分）。
func (s *Stepper) neighbor(g *Grid, i, j, di, dj int) float64 {
	nx, ny := g.NX, g.NY
	ni, nj := i+di, j+dj
	if ni < 0 {
		switch s.bc.Left.Type {
		case Periodic:
			ni = nx - 1
		default: // Neumann：零通量，鬼格 = 边界格自身
			ni = 0
		}
	} else if ni >= nx {
		switch s.bc.Right.Type {
		case Periodic:
			ni = 0
		default:
			ni = nx - 1
		}
	}
	if nj < 0 {
		switch s.bc.Bottom.Type {
		case Periodic:
			nj = ny - 1
		default:
			nj = 0
		}
	} else if nj >= ny {
		switch s.bc.Top.Type {
		case Periodic:
			nj = 0
		default:
			nj = ny - 1
		}
	}
	return g.Data[nj*nx+ni]
}

// fixed 判断 (i,j) 是否位于 Dirichlet 边界上（值固定，不更新）。
func (s *Stepper) fixed(i, j int) bool {
	bc := s.bc
	if i == 0 && bc.Left.Type == Dirichlet {
		return true
	}
	if i == s.Cur.NX-1 && bc.Right.Type == Dirichlet {
		return true
	}
	if j == 0 && bc.Bottom.Type == Dirichlet {
		return true
	}
	if j == s.Cur.NY-1 && bc.Top.Type == Dirichlet {
		return true
	}
	return false
}

// sourceAt 返回 (i,j) 处的源项强度（多个源区域重叠时强度相加）。
func (s *Stepper) sourceAt(i, j int) float64 {
	x := float64(i) / float64(s.Cur.NX-1)
	y := float64(j) / float64(s.Cur.NY-1)
	rate := 0.0
	for k := range s.srcs {
		src := &s.srcs[k]
		if x >= src.X0 && x <= src.X1 && y >= src.Y0 && y <= src.Y1 {
			rate += src.Rate
		}
	}
	return rate
}

// cellUpdate 计算 (i,j) 处的新值。
func (s *Stepper) cellUpdate(i, j int) float64 {
	c := s.Cfg
	g := s.Cur
	u := g.At(i, j)
	lap := (s.neighbor(g, i, j, 1, 0) - 2*u + s.neighbor(g, i, j, -1, 0)) / (c.DX * c.DX)
	lap += (s.neighbor(g, i, j, 0, 1) - 2*u + s.neighbor(g, i, j, 0, -1)) / (c.DY * c.DY)
	return u + c.DT*(c.Diffusivity*lap+s.sourceAt(i, j))
}

// stepRows 计算行区间 [j0,j1) 的新值，写入 Next。可被并行调用。
func (s *Stepper) stepRows(j0, j1 int) {
	for j := j0; j < j1; j++ {
		for i := 0; i < s.Cur.NX; i++ {
			if s.fixed(i, j) {
				s.Next.Set(i, j, s.Cur.At(i, j)) // Dirichlet 边界保持
			} else {
				s.Next.Set(i, j, s.cellUpdate(i, j))
			}
		}
	}
}

// commit 交换缓冲区并检查非物理值。
func (s *Stepper) commit() error {
	s.Cur, s.Next = s.Next, s.Cur
	s.Step++
	for k, v := range s.Cur.Data {
		if math.IsNaN(v) || math.IsInf(v, 0) {
			return stabilityErr(fmt.Sprintf(
				"stability: non-finite value %v at step %d, cell (i=%d,j=%d)",
				v, s.Step, k%s.Cur.NX, k/s.Cur.NX))
		}
	}
	return nil
}

// Advance 串行推进 n 步；出现 NaN/Inf 时返回 ErrKindStability 错误。
func (s *Stepper) Advance(n int) error {
	for k := 0; k < n; k++ {
		s.applyPhaseForStep(s.Step)
		s.stepRows(0, s.Cur.NY)
		if err := s.commit(); err != nil {
			return err
		}
	}
	return nil
}
