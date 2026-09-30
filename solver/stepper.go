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
type Stepper struct {
	Cfg  *Config
	Cur  *Grid
	Next *Grid
	Step int // 已完成的步数

	// 活动配置钩子：为 nil 时使用 Cfg 中的单套边界/源（旧模式，行为不变）。
	// 阶段运行器在每一步推进前设置它们，从而实现按步切换边界与源项。
	activeBC      *Boundaries
	activeSources []PhaseSource
}

// SetActive 设定下一步推进使用的边界条件与热源集合。
// 传 nil（仅指 sources）表示该阶段无源，而不是回退到 Cfg.Src；
// 回退到旧单套配置请调用 ResetActive。
func (s *Stepper) SetActive(bc *Boundaries, sources []PhaseSource) {
	s.activeBC = bc
	if sources == nil {
		// 用非 nil 空切片区分「阶段模式但无源」与「旧模式」。
		s.activeSources = []PhaseSource{}
	} else {
		s.activeSources = sources
	}
}

// ResetActive 回退到 Cfg.BC / Cfg.Src 的旧单套配置（旧代码路径默认即此状态）。
func (s *Stepper) ResetActive() {
	s.activeBC = nil
	s.activeSources = nil
}

// activeBoundaries 返回当前步使用的边界条件。
func (s *Stepper) activeBoundaries() Boundaries {
	if s.activeBC != nil {
		return *s.activeBC
	}
	return s.Cfg.BC
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
	return &Stepper{Cfg: cfg, Cur: initial.Clone(), Next: NewGrid(cfg.NX, cfg.NY)}, nil
}

// neighbor 按边界条件取 (i+di, j+dj) 处的值。
// di,dj 只有一个非零（五点差分）。
func (s *Stepper) neighbor(g *Grid, i, j, di, dj int) float64 {
	nx, ny := g.NX, g.NY
	bc := s.activeBoundaries()
	ni, nj := i+di, j+dj
	if ni < 0 {
		switch bc.Left.Type {
		case Periodic:
			ni = nx - 1
		default: // Neumann：零通量，鬼格 = 边界格自身
			ni = 0
		}
	} else if ni >= nx {
		switch bc.Right.Type {
		case Periodic:
			ni = 0
		default:
			ni = nx - 1
		}
	}
	if nj < 0 {
		switch bc.Bottom.Type {
		case Periodic:
			nj = ny - 1
		default:
			nj = 0
		}
	} else if nj >= ny {
		switch bc.Top.Type {
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
	bc := s.activeBoundaries()
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

// sourceAt 返回 (i,j) 处的源项强度（多源区重叠时强度相加）。
func (s *Stepper) sourceAt(i, j int) float64 {
	if s.activeSources != nil {
		nx := float64(s.Cur.NX - 1)
		ny := float64(s.Cur.NY - 1)
		x := float64(i) / nx
		y := float64(j) / ny
		rate := 0.0
		for k := range s.activeSources {
			src := &s.activeSources[k]
			if x >= src.X0 && x <= src.X1 && y >= src.Y0 && y <= src.Y1 {
				rate += src.Rate
			}
		}
		return rate
	}
	src := s.Cfg.Src
	if src == nil {
		return 0
	}
	x := float64(i) / float64(s.Cur.NX-1)
	y := float64(j) / float64(s.Cur.NY-1)
	if x >= src.X0 && x <= src.X1 && y >= src.Y0 && y <= src.Y1 {
		return src.Rate
	}
	return 0
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
		s.stepRows(0, s.Cur.NY)
		if err := s.commit(); err != nil {
			return err
		}
	}
	return nil
}
