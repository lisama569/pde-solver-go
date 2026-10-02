package solver

import (
	"fmt"
	"math"
)

// BoundaryBudget 记录四边的累计净流入（正 = 流入计算域）。
type BoundaryBudget struct {
	Left, Right, Top, Bottom float64
}

// Total 返回四边净流入之和。
func (b BoundaryBudget) Total() float64 {
	return b.Left + b.Right + b.Top + b.Bottom
}

// BoundaryInflux 按当前场 Cur 计算各边界的瞬时净流入（单位时间）。
//
// 只有 Dirichlet 边有净交换：内部格与固定边界格之间的扩散通量
// D·(u_boundary − u_inside)/d·（垂直边长）。Neumann 边鬼格等于边界格自身，
// 净通量恒为 0；Periodic 边的交换在环上两两相消，净通量恒为 0。
// 角点处若内侧相邻格本身也是固定格（两边均为 Dirichlet），该通量不进入
// 任何被更新的格点，必须跳过——这正是 fixed() 检查的作用。
// 该定义与更新算子严格对应：全场非固定格的更新量之和 = 边界净流入 + 源注入，
// 因此收支对账可闭合到浮点舍入量级。
func (s *Stepper) BoundaryInflux() BoundaryBudget {
	c := s.Cfg
	g := s.Cur
	var b BoundaryBudget
	if s.bc.Left.Type == Dirichlet {
		for j := 0; j < g.NY; j++ {
			if s.fixed(1, j) {
				continue // 内侧格固定（角部），无质量交换
			}
			b.Left += c.Diffusivity * (g.At(0, j) - g.At(1, j)) / c.DX * c.DY
		}
	}
	if s.bc.Right.Type == Dirichlet {
		for j := 0; j < g.NY; j++ {
			if s.fixed(g.NX-2, j) {
				continue
			}
			b.Right += c.Diffusivity * (g.At(g.NX-1, j) - g.At(g.NX-2, j)) / c.DX * c.DY
		}
	}
	if s.bc.Bottom.Type == Dirichlet {
		for i := 0; i < g.NX; i++ {
			if s.fixed(i, 1) {
				continue
			}
			b.Bottom += c.Diffusivity * (g.At(i, 0) - g.At(i, 1)) / c.DY * c.DX
		}
	}
	if s.bc.Top.Type == Dirichlet {
		for i := 0; i < g.NX; i++ {
			if s.fixed(i, g.NY-2) {
				continue
			}
			b.Top += c.Diffusivity * (g.At(i, g.NY-1) - g.At(i, g.NY-2)) / c.DY * c.DX
		}
	}
	return b
}

// SourceRates 返回当前生效的每个源的瞬时注入率（单位时间）：
// rate × 区域内**被实际更新的**格点数 × dx·dy。与 sourceAt 使用同一区域
// 判定，且跳过 Dirichlet 固定格——固定格不被 stepRows 更新，源强从未
// 施加到这些格上。因此贴边/贴角的热源区只按真正进入场内的部分记账，
// 收支对账才能闭合。
func (s *Stepper) SourceRates() []float64 {
	rates := make([]float64, len(s.srcs))
	if len(s.srcs) == 0 {
		return rates
	}
	c := s.Cfg
	nx, ny := s.Cur.NX, s.Cur.NY
	for j := 0; j < ny; j++ {
		y := float64(j) / float64(ny-1)
		for i := 0; i < nx; i++ {
			if s.fixed(i, j) {
				continue // Dirichlet 固定格不更新，源实际未施加
			}
			x := float64(i) / float64(nx-1)
			for k := range s.srcs {
				src := &s.srcs[k]
				if x >= src.X0 && x <= src.X1 && y >= src.Y0 && y <= src.Y1 {
					rates[k] += src.Rate * c.DX * c.DY
				}
			}
		}
	}
	return rates
}

// PhaseReport 是一个阶段的收支对账报告。
type PhaseReport struct {
	PhaseIndex   int            // 阶段索引（无阶段计划时为 -1）
	FromStep     int            // 本报告覆盖的起始步号
	ToStep       int            // 本报告覆盖的结束步号
	SourceInject []float64      // 各热源累计注入量（只计实际施加到非固定格的部分）
	Boundary     BoundaryBudget // 各边界累计净流入（扩散通量 + 阶段激活时的 Dirichlet 重置量）
	Reset        BoundaryBudget // 本报告区间内阶段激活写入 Dirichlet 值造成的逐边质量跳变（已含在 Boundary 中）
	MassStart    float64        // 阶段开始总量（施加本阶段 Dirichlet 重置之前的场）
	MassEnd      float64        // 阶段结束总量
}

// DeltaMass 返回阶段内总量变化。
func (r *PhaseReport) DeltaMass() float64 {
	return r.MassEnd - r.MassStart
}

// SourceTotal 返回全部热源累计注入之和。
func (r *PhaseReport) SourceTotal() float64 {
	sum := 0.0
	for _, v := range r.SourceInject {
		sum += v
	}
	return sum
}

// Closure 返回对账残差：ΔMass − (源注入 + 边界净流入)。
// 边界净流入含阶段激活时的 Dirichlet 重置量（见 Reset 字段），因此
// 跨阶段修改固定边值造成的总量跳变也在账内。理论上恒为 0，实际残差
// 来自浮点求和顺序，量级约为总量的 1e-15 倍。
func (r *PhaseReport) Closure() float64 {
	return r.DeltaMass() - r.SourceTotal() - r.Boundary.Total()
}

// ClosureTol 返回对账容差：1e-9 × 收支规模（各账目绝对值之和，下限 1），
// 相对浮点舍入（~1e-15）留有约 6 个数量级的余量，同时足以发现
// 任何物理量级（~1e-6 相对）的漏账。
func (r *PhaseReport) ClosureTol() float64 {
	scale := 1.0 + math.Abs(r.MassStart) + math.Abs(r.MassEnd)
	for _, v := range r.SourceInject {
		scale += math.Abs(v)
	}
	scale += math.Abs(r.Boundary.Left) + math.Abs(r.Boundary.Right) +
		math.Abs(r.Boundary.Top) + math.Abs(r.Boundary.Bottom)
	return 1e-9 * scale
}

// Balanced 判定对账是否在容差内闭合。
func (r *PhaseReport) Balanced() bool {
	return math.Abs(r.Closure()) <= r.ClosureTol()
}

// TotalDeltaMass 返回各阶段报告的总量变化之和。由于报告区间首尾相接
// 且阶段切换的 Dirichlet 重置量已计入对应报告的 Boundary，该和严格等于
// 任务从头到尾的真实总量变化（末报告 MassEnd − 首报告 MassStart）。
func TotalDeltaMass(reports []PhaseReport) float64 {
	sum := 0.0
	for i := range reports {
		sum += reports[i].DeltaMass()
	}
	return sum
}

// CheckBudgets 校验所有阶段报告都在容差内闭合。不闭合时返回
// ErrKindBudget 错误，信息中包含阶段索引、步区间、残差与容差；
// 具体是哪个热源或哪条边对不上，可查该报告的 SourceInject / Boundary /
// Reset 各分项（逐源、逐边）。全部闭合返回 nil。
func CheckBudgets(reports []PhaseReport) error {
	for i := range reports {
		r := &reports[i]
		if !r.Balanced() {
			return budgetErr(fmt.Sprintf(
				"budget: phase %d steps[%d,%d) does not close: residual %+.6e exceeds tol %.3e "+
					"(deltaMass %+.6e, sources %+.6e, boundary{L:%+.6e R:%+.6e T:%+.6e B:%+.6e} "+
					"incl. reset{L:%+.6e R:%+.6e T:%+.6e B:%+.6e})",
				r.PhaseIndex, r.FromStep, r.ToStep, r.Closure(), r.ClosureTol(),
				r.DeltaMass(), r.SourceTotal(),
				r.Boundary.Left, r.Boundary.Right, r.Boundary.Top, r.Boundary.Bottom,
				r.Reset.Left, r.Reset.Right, r.Reset.Top, r.Reset.Bottom))
		}
	}
	return nil
}
