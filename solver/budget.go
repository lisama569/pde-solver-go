package solver

import (
	"fmt"
	"math"
	"strings"
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
// rate × 区域内实际受激格点数 × dx·dy。与 sourceAt / stepRows 使用同一
// 区域判定与同一套固定格判定——贴边、贴角的源区若压到 Dirichlet 固定格，
// 那些格点在更新算子中不接受源项（值保持边界值），因此也不能计入注入量，
// 否则账目会系统性偏大并误报不平衡。Neumann/Periodic 边不是固定格，
// 贴这些边的源区照常全额计入。
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
				continue // Dirichlet 固定格：源项不会施加到场上
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

// PhaseReport 是一个阶段（或一次运行的一段步区间）的收支对账报告。
type PhaseReport struct {
	PhaseIndex    int            // 阶段索引（无阶段计划时为 -1）
	FromStep      int            // 本报告覆盖的起始步号
	ToStep        int            // 本报告覆盖的结束步号
	SourceInject  []float64      // 各热源累计实际注入量（贴 Dirichlet 边的固定格不计）
	Boundary      BoundaryBudget // 各边界累计扩散净流入（正 = 流入计算域）
	BoundaryReset BoundaryBudget // 阶段激活时边界重置的逐边总量跳变（新值−旧值）
	MassStart     float64        // 报告起点总量（任何边界重置之前）
	MassEnd       float64        // 报告终点总量
}

// DeltaMass 返回阶段内总量变化（含阶段激活时边界重置造成的跳变）。
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

// ResetTotal 返回边界重置跳变之和。
func (r *PhaseReport) ResetTotal() float64 {
	return r.BoundaryReset.Total()
}

// Closure 返回对账残差：ΔMass − (源注入 + 边界扩散净流入 + 边界重置跳变)。
//
// 三项账目口径分别为：
//   - 源注入：只统计实际施加到非固定格上的量（贴 Dirichlet 边/角的固定格不接受源项）；
//   - 边界扩散净流入：更新算子中 Dirichlet 固定格与相邻内部格之间的扩散交换，
//     按更新前的场逐边累计；
//   - 边界重置跳变：阶段激活瞬间把新边界值写入固定格造成的总量跳变，逐边记录
//     （改 Dirichlet 数值、Dirichlet↔Neumann 切换都覆盖；任务首次激活
//     阶段0 的 Dirichlet 值相对初始场的跳变也计入阶段0）。
//
// 理论上恒为 0，实际残差来自浮点求和顺序，量级约为总量的 1e-15 倍。
func (r *PhaseReport) Closure() float64 {
	return r.DeltaMass() - r.SourceTotal() - r.Boundary.Total() - r.ResetTotal()
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
	scale += math.Abs(r.BoundaryReset.Left) + math.Abs(r.BoundaryReset.Right) +
		math.Abs(r.BoundaryReset.Top) + math.Abs(r.BoundaryReset.Bottom)
	return 1e-9 * scale
}

// Balanced 判定对账是否在容差内闭合。
func (r *PhaseReport) Balanced() bool {
	return math.Abs(r.Closure()) <= r.ClosureTol()
}

// VerifyReports 校验一次 RunPhases 输出的全部报告：
//  1. 每份报告内部闭合：ΔMass = 源注入 + 边界扩散净流入 + 边界重置跳变；
//  2. 报告区间首尾相接（后一份 MassStart == 前一份 MassEnd）；
//  3. 各报告 ΔMass 之和 == 末份 MassEnd − 首份 MassStart（跨报告总账）。
//
// 任一不满足返回 ErrKindBudget 错误，信息中带阶段/区间索引、残差、容差，
// 以及逐源、逐边（扩散与重置分列）的全部账目，可直接定位到具体热源或边界。
func VerifyReports(reports []PhaseReport) error {
	if len(reports) == 0 {
		return budgetErr("budget: no reports to verify")
	}
	for k := range reports {
		if err := reports[k].verify(); err != nil {
			return err
		}
		if k > 0 {
			prev, cur := &reports[k-1], &reports[k]
			if prev.ToStep != cur.FromStep {
				return budgetErr(fmt.Sprintf(
					"budget: report %d (phase %d, steps [%d,%d)) and report %d (phase %d, steps [%d,%d)) are not contiguous",
					k-1, prev.PhaseIndex, prev.FromStep, prev.ToStep,
					k, cur.PhaseIndex, cur.FromStep, cur.ToStep))
			}
		}
	}
	sumDelta := 0.0
	for k := range reports {
		sumDelta += reports[k].DeltaMass()
	}
	totalChange := reports[len(reports)-1].MassEnd - reports[0].MassStart
	residual := sumDelta - totalChange
	tol := 1e-9 * (1 + math.Abs(reports[0].MassStart) + math.Abs(reports[len(reports)-1].MassEnd))
	for k := range reports {
		tol += 1e-9 * math.Abs(reports[k].DeltaMass())
	}
	if math.Abs(residual) > tol {
		return budgetErr(fmt.Sprintf(
			"budget: cross-report ledger mismatch: ΣΔMass=%+.6e but true change=%+.6e (residual %+.3e, tol %.3e)",
			sumDelta, totalChange, residual, tol))
	}
	return nil
}

// verify 校验单份报告内部闭合；不闭合时给出可定位的明细错误。
func (r *PhaseReport) verify() error {
	residual := r.Closure()
	if math.Abs(residual) <= r.ClosureTol() {
		return nil
	}
	var sb strings.Builder
	fmt.Fprintf(&sb, "budget: phase %d steps [%d,%d) does not close: residual=%+.6e tol=%.3e\n",
		r.PhaseIndex, r.FromStep, r.ToStep, residual, r.ClosureTol())
	fmt.Fprintf(&sb, "  ΔMass=%+.6e (mass %+.6e -> %+.6e)\n",
		r.DeltaMass(), r.MassStart, r.MassEnd)
	for si, inj := range r.SourceInject {
		fmt.Fprintf(&sb, "  source %d injected=%+.6e\n", si, inj)
	}
	name := []string{"left", "right", "top", "bottom"}
	edges := []float64{r.Boundary.Left, r.Boundary.Right, r.Boundary.Top, r.Boundary.Bottom}
	resets := []float64{r.BoundaryReset.Left, r.BoundaryReset.Right, r.BoundaryReset.Top, r.BoundaryReset.Bottom}
	for i, n := range name {
		fmt.Fprintf(&sb, "  boundary %s influx=%+.6e reset=%+.6e\n", n, edges[i], resets[i])
	}
	fmt.Fprintf(&sb, "  sums: src=%+.6e influx=%+.6e reset=%+.6e",
		r.SourceTotal(), r.Boundary.Total(), r.ResetTotal())
	if hint := locateImbalance(r); hint != "" {
		fmt.Fprintf(&sb, "\n  likely culprit: %s", hint)
	}
	return budgetErr(sb.String())
}

// locateImbalance 在不闭合时给出最可能漏账/错账的单项（源或边）提示。
// 判据：各项对"期望贡献"的偏离量级；残差方向与单项缺额一致时优先点名。
func locateImbalance(r *PhaseReport) string {
	type term struct {
		label string
		value float64
	}
	var terms []term
	for si, inj := range r.SourceInject {
		terms = append(terms, term{fmt.Sprintf("source %d", si), inj})
	}
	name := []string{"left", "right", "top", "bottom"}
	edges := []float64{r.Boundary.Left, r.Boundary.Right, r.Boundary.Top, r.Boundary.Bottom}
	resets := []float64{r.BoundaryReset.Left, r.BoundaryReset.Right, r.BoundaryReset.Top, r.BoundaryReset.Bottom}
	for i, n := range name {
		terms = append(terms, term{fmt.Sprintf("boundary %s influx", n), edges[i]})
		terms = append(terms, term{fmt.Sprintf("boundary %s reset", n), resets[i]})
	}
	// 期望残差为 0：仅当某个账项的量级与残差接近（25% 以内）时才点名，
	// 说明该项疑似整笔漏记/重复记；量级不匹配时宁可不猜，逐源逐边明细
	// 本身已足以人工定位。
	residual := r.Closure()
	best := ""
	bestRatio := 0.25
	for _, t := range terms {
		if t.value == 0 {
			continue
		}
		ratio := math.Abs(math.Abs(residual)-math.Abs(t.value)) / math.Abs(t.value)
		if ratio < bestRatio {
			bestRatio = ratio
			best = t.label
		}
	}
	if best != "" {
		return fmt.Sprintf("%s (magnitude matches the residual within %.0f%%)", best, bestRatio*100)
	}
	return ""
}
