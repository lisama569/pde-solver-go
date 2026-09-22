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
	ni, nj := i+di, j+dj
	if ni < 0 {
		switch s.Cfg.BC.Left.Type {
		case Periodic:
			ni = nx - 1
		default: // Neumann：零通量，鬼格 = 边界格自身
			ni = 0
		}
	} else if ni >= nx {
		switch s.Cfg.BC.Right.Type {
		case Periodic:
			ni = 0
		default:
			ni = nx - 1
		}
	}
	if nj < 0 {
		switch s.Cfg.BC.Bottom.Type {
		case Periodic:
			nj = ny - 1
		default:
			nj = 0
		}
	} else if nj >= ny {
		switch s.Cfg.BC.Top.Type {
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
	bc := s.Cfg.BC
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

// sourceAt 返回 (i,j) 处的源项强度。
func (s *Stepper) sourceAt(i, j int) float64 {
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
