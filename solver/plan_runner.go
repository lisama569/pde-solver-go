package solver

import (
	"fmt"
	"math"
)

// 默认对账容差：闭合判定用「相对 + 绝对」混合容差
//
//	tol = tolAbs + tolRel * max(|StartMass|,|EndMass|,|注入|,|通量|,|重置|)
//
// 离散格式下三者应恒等相消（见 PhaseReport 注释），残差只有浮点舍入量级。
const (
	DefaultBalanceTolAbs = 1e-12
	DefaultBalanceTolRel = 1e-9
)

// EdgeFlow 是单条边在一个阶段内的净通量（正值=流入域内，负值=流出）。
type EdgeFlow struct {
	Edge string  // "left" / "right" / "top" / "bottom"
	Flow float64 // 净流入
}

// SourceInjection 是单个热源在一个阶段内的总注入（已乘 dt·dx·dy）。
type SourceInjection struct {
	Name string  // 未命名时使用 "source#<阶段内序号>"
	Rate float64 // 阶段总注入
}

// PhaseReport 是单个阶段结束时的收支对账报告。
//
// 离散收支闭合关系（正值均为「进入域内」）：
//
//	EndMass - StartMass = Σ 各源注入 + Σ 四边净流入 + BoundaryReset + 舍入
//
// 其中：
//   - Σ源注入：该阶段每步落在源区内、且非 Dirichlet 固定格点上的 dt·dx·dy·Rate 之和；
//   - Σ边通量：Dirichlet 边按五点差分的边界通量（流入为正），Neumann 边为 0，
//     Periodic 边两两相消为 0（逐对计数，不依赖求和顺序）；
//   - BoundaryReset：进入阶段时 Dirichlet 边被重置为该阶段固定值导致的总量跳变
//     （仅阶段计划模式发生；阶段 0 也会把 Dirichlet 边吸附到声明值，跳变计入阶段 0）。
//
// Residual = (Σ注入 + Σ通量 + BoundaryReset) - (EndMass-StartMass)。
type PhaseReport struct {
	PhaseIndex    int
	StartStep     int
	EndStep       int
	StartMass     float64 // 阶段入口总量（重置前）
	EndMass       float64 // 阶段出口总量
	MassChange    float64 // EndMass - StartMass
	Sources       []SourceInjection
	Edges         []EdgeFlow
	BoundaryReset float64
	Residual      float64 // 收入合计 - 总量变化
	Tolerance     float64 // 本次采用的判定容差
	Closed        bool    // |Residual| <= Tolerance
}

// phaseAccum 阶段内累计器。
type phaseAccum struct {
	index     int
	fromStep  int // 入账窗口起点（完整阶段=StartStep；阶段内恢复=恢复步号）
	startMass float64
	sourceInj []float64
	edges     [4]float64 // left,right,top,bottom
	reset     float64
	steps     int
}

var edgeNames = [4]string{"left", "right", "top", "bottom"}

// Runner 按阶段计划推进求解，并在每阶段结束时给出收支对账。
type Runner struct {
	cfg   *Config
	plan  *PhasePlan
	st    *Stepper
	tolA  float64
	tolR  float64
	phase int
	acc   *phaseAccum
}

// NewRunner 从头创建阶段推进器。先校验配置与阶段计划（ErrKindConfig），
// 再校验稳定性（ErrKindStability），二者在推进前都必须通过。
// 目标总步数由 cfg.Steps 给出，阶段计划必须至少覆盖到 cfg.Steps。
func NewRunner(cfg *Config, plan *PhasePlan, initial *Grid) (*Runner, error) {
	if err := cfg.Validate(); err != nil {
		return nil, err
	}
	if plan == nil {
		return nil, configErr("config: NewRunner requires a non-nil phase plan")
	}
	if err := ValidatePlan(plan); err != nil {
		return nil, err
	}
	if plan.End() < cfg.Steps {
		return nil, configErr(fmt.Sprintf(
			"config: phase plan covers steps [0,%d) but cfg.Steps=%d is beyond it", plan.End(), cfg.Steps))
	}
	if err := cfg.CheckStability(); err != nil {
		return nil, err
	}
	if initial == nil || initial.NX != cfg.NX || initial.NY != cfg.NY {
		return nil, configErr(fmt.Sprintf("config: initial grid size mismatch, want %dx%d", cfg.NX, cfg.NY))
	}
	st, err := NewStepper(cfg, initial)
	if err != nil {
		return nil, err
	}
	r := &Runner{cfg: cfg, plan: plan, st: st, tolA: DefaultBalanceTolAbs, tolR: DefaultBalanceTolRel}
	r.startPhase(0) // 阶段 0 入口：吸附 Dirichlet 边并记录入口总量
	return r, nil
}

// ResumeRunner 从已校验的状态恢复阶段推进器。恢复点可位于阶段内部
// （start < EndStep）或刚好在阶段边界（start == 某阶段 StartStep）。
// 恢复前用 State.Verify 严格比对网格、步长、系数与阶段计划指纹，任何不一致
// 都以 ErrKindConfig 拒绝；调用方一般先 LoadState 再调用本函数。
func ResumeRunner(cfg *Config, plan *PhasePlan, state *State) (*Runner, error) {
	if err := cfg.Validate(); err != nil {
		return nil, err
	}
	if plan == nil {
		return nil, configErr("config: ResumeRunner requires a non-nil phase plan")
	}
	if err := ValidatePlan(plan); err != nil {
		return nil, err
	}
	if state == nil {
		return nil, configErr("state: cannot resume from nil state")
	}
	if err := state.Verify(cfg, plan); err != nil {
		return nil, err
	}
	if plan.End() < cfg.Steps {
		return nil, configErr(fmt.Sprintf(
			"config: phase plan covers [0,%d) but cfg.Steps=%d is beyond it", plan.End(), cfg.Steps))
	}
	if state.Step > plan.End() || state.Step > cfg.Steps {
		return nil, configErr(fmt.Sprintf(
			"state: saved step %d beyond plan end %d or cfg.Steps %d", state.Step, plan.End(), cfg.Steps))
	}
	if err := cfg.CheckStability(); err != nil {
		return nil, err
	}
	g := state.Grid()
	st, err := NewStepper(cfg, g)
	if err != nil {
		return nil, err
	}
	st.Step = state.Step
	r := &Runner{cfg: cfg, plan: plan, st: st, tolA: DefaultBalanceTolAbs, tolR: DefaultBalanceTolRel}
	pi := plan.phaseAt(state.Step)
	r.phase = pi
	ph := &plan.Phases[pi]
	atBoundary := state.Step == ph.StartStep
	// 刚好停在阶段边界：与一次性运行相同，进入新阶段时要把 Dirichlet 边吸附到
	// 新固定值，吸附跳变计入该阶段；阶段内部恢复则场已是连续轨迹的一部分，不吸附。
	r.applyActive(pi)
	mBefore := TotalMass(st.Cur, cfg.DX, cfg.DY)
	if atBoundary {
		r.resetDirichlet(ph.BC)
	}
	mAfter := TotalMass(st.Cur, cfg.DX, cfg.DY)
	r.acc = &phaseAccum{
		index:     pi,
		fromStep:  state.Step,
		startMass: mBefore,
		reset:     mAfter - mBefore,
		sourceInj: make([]float64, len(ph.Sources)),
	}
	return r, nil
}

// SetBalanceTolerance 覆盖对账容差（abs 与 rel 必须 >= 0）。
func (r *Runner) SetBalanceTolerance(abs, rel float64) error {
	if math.IsNaN(abs) || math.IsInf(abs, 0) || math.IsNaN(rel) || math.IsInf(rel, 0) || abs < 0 || rel < 0 {
		return configErr(fmt.Sprintf("config: invalid balance tolerance abs=%g rel=%g", abs, rel))
	}
	r.tolA, r.tolR = abs, rel
	return nil
}

// startPhase 进入阶段 pi：设置活动边界/源，并把 Dirichlet 边吸附到该阶段固定值。
// 吸附导致的总量跳变计入该阶段 BoundaryReset。
func (r *Runner) startPhase(pi int) {
	ph := &r.plan.Phases[pi]
	r.phase = pi
	r.applyActive(pi)
	mBefore := TotalMass(r.st.Cur, r.cfg.DX, r.cfg.DY)
	r.resetDirichlet(ph.BC)
	mAfter := TotalMass(r.st.Cur, r.cfg.DX, r.cfg.DY)
	r.acc = &phaseAccum{
		index:     pi,
		fromStep:  ph.StartStep,
		startMass: mBefore,
		reset:     mAfter - mBefore,
		sourceInj: make([]float64, len(ph.Sources)),
	}
}

// applyActive 把阶段边界与源装入步进器活动钩子。
func (r *Runner) applyActive(pi int) {
	ph := &r.plan.Phases[pi]
	r.st.SetActive(&ph.BC, ph.Sources)
}

// resetDirichlet 把当前为 Dirichlet 的边的边界格点设为声明固定值。
// 四角若两边都为 Dirichlet，按 right/top 后写的顺序取该边值。
func (r *Runner) resetDirichlet(bc Boundaries) {
	g := r.st.Cur
	nx, ny := g.NX, g.NY
	if bc.Left.Type == Dirichlet {
		for j := 0; j < ny; j++ {
			g.Set(0, j, bc.Left.Value)
		}
	}
	if bc.Bottom.Type == Dirichlet {
		for i := 0; i < nx; i++ {
			g.Set(i, 0, bc.Bottom.Value)
		}
	}
	if bc.Right.Type == Dirichlet {
		for j := 0; j < ny; j++ {
			g.Set(nx-1, j, bc.Right.Value)
		}
	}
	if bc.Top.Type == Dirichlet {
		for i := 0; i < nx; i++ {
			g.Set(i, ny-1, bc.Top.Value)
		}
	}
}

// isFixed 判断当前活动边界下 (i,j) 是否为 Dirichlet 固定格点。
func (r *Runner) isFixed(i, j int) bool {
	return r.st.fixed(i, j)
}

// accountStep 用推进前的场 g 计算这一步的源注入与边界净通量，并累加到阶段累计器。
// 这些量与实际更新逐格点对应，离散恒等式保证闭合（仅余浮点舍入）。
func (r *Runner) accountStep(g *Grid, acc *phaseAccum) {
	c := r.cfg
	nx, ny := g.NX, g.NY
	dx, dy, dt, dd := c.DX, c.DY, c.DT, c.Diffusivity
	bc := r.st.activeBoundaries()
	area := dx * dy

	// 源注入：仅计入非固定格点（Dirichlet 固定格点不参与更新）。
	for si := range r.st.activeSources {
		src := &r.st.activeSources[si]
		sum := 0.0
		for j := 0; j < ny; j++ {
			y := float64(j) / float64(ny-1)
			if y < src.Y0 || y > src.Y1 {
				continue
			}
			for i := 0; i < nx; i++ {
				if r.isFixed(i, j) {
					continue
				}
				x := float64(i) / float64(nx-1)
				if x >= src.X0 && x <= src.X1 {
					sum += src.Rate
				}
			}
		}
		acc.sourceInj[si] += dt * area * sum
	}

	// 边界净通量（流入为正 = 边值减去内侧邻值），按「非整行固定的行集 /
	// 非整列固定的列集」求和，使角点附近也不漏项：
	// 行集随 bottom/top 是否 Dirichlet 收缩，列集随 left/right 是否 Dirichlet 收缩。
	// 周期边逐对相消，对总量贡献恒为 0，不计入。
	var left, right, top, bottom float64
	if dd != 0 {
		rowLo, rowHi := 0, ny
		if bc.Bottom.Type == Dirichlet {
			rowLo = 1
		}
		if bc.Top.Type == Dirichlet {
			rowHi = ny - 1
		}
		colLo, colHi := 0, nx
		if bc.Left.Type == Dirichlet {
			colLo = 1
		}
		if bc.Right.Type == Dirichlet {
			colHi = nx - 1
		}
		if bc.Left.Type == Dirichlet {
			for j := rowLo; j < rowHi; j++ {
				left += g.At(0, j) - g.At(1, j)
			}
			left *= dd * dt * (dy / dx)
		}
		if bc.Right.Type == Dirichlet {
			for j := rowLo; j < rowHi; j++ {
				right += g.At(nx-1, j) - g.At(nx-2, j)
			}
			right *= dd * dt * (dy / dx)
		}
		if bc.Bottom.Type == Dirichlet {
			for i := colLo; i < colHi; i++ {
				bottom += g.At(i, 0) - g.At(i, 1)
			}
			bottom *= dd * dt * (dx / dy)
		}
		if bc.Top.Type == Dirichlet {
			for i := colLo; i < colHi; i++ {
				top += g.At(i, ny-1) - g.At(i, ny-2)
			}
			top *= dd * dt * (dx / dy)
		}
	}
	acc.edges[0] += left
	acc.edges[1] += right
	acc.edges[2] += top
	acc.edges[3] += bottom
	acc.steps++
}

// finalizePhase 在到达阶段末端步号后生成对账报告。
func (r *Runner) finalizePhase(acc *phaseAccum) PhaseReport {
	ph := &r.plan.Phases[acc.index]
	endMass := TotalMass(r.st.Cur, r.cfg.DX, r.cfg.DY)
	srcs := make([]SourceInjection, len(ph.Sources))
	totSrc := 0.0
	for si, src := range ph.Sources {
		name := src.Name
		if name == "" {
			name = fmt.Sprintf("source#%d", si)
		}
		srcs[si] = SourceInjection{Name: name, Rate: acc.sourceInj[si]}
		totSrc += acc.sourceInj[si]
	}
	edges := make([]EdgeFlow, 4)
	totEdge := 0.0
	for e := 0; e < 4; e++ {
		edges[e] = EdgeFlow{Edge: edgeNames[e], Flow: acc.edges[e]}
		totEdge += acc.edges[e]
	}
	change := endMass - acc.startMass
	in := totSrc + totEdge + acc.reset
	residual := in - change
	scale := math.Max(math.Abs(acc.startMass), math.Abs(endMass))
	scale = math.Max(scale, math.Abs(totSrc))
	scale = math.Max(scale, math.Abs(totEdge))
	scale = math.Max(scale, math.Abs(acc.reset))
	tol := r.tolA + r.tolR*scale
	return PhaseReport{
		PhaseIndex:    acc.index,
		StartStep:     acc.fromStep,
		EndStep:       ph.EndStep,
		StartMass:     acc.startMass,
		EndMass:       endMass,
		MassChange:    change,
		Sources:       srcs,
		Edges:         edges,
		BoundaryReset: acc.reset,
		Residual:      residual,
		Tolerance:     tol,
		Closed:        math.Abs(residual) <= tol,
	}
}

// Advance 从当前绝对步号再推进 steps 步（workers 为按行子域数，1=串行）。
// 每当某阶段在其 EndStep 处走完，就返回一份该阶段的 PhaseReport；
// 跨越多个阶段边界时按时间顺序返回多份。出现 NaN/Inf 时返回 ErrKindStability。
func (r *Runner) Advance(steps, workers int) ([]PhaseReport, error) {
	if steps < 0 {
		return nil, configErr(fmt.Sprintf("config: advance steps must be >= 0, got %d", steps))
	}
	target := r.st.Step + steps
	if target > r.cfg.Steps {
		return nil, configErr(fmt.Sprintf(
			"config: advance to step %d exceeds cfg.Steps=%d", target, r.cfg.Steps))
	}
	if workers < 1 {
		workers = 1
	}
	var reports []PhaseReport
	for r.st.Step < target {
		pi := r.plan.phaseAt(r.st.Step)
		if pi != r.phase {
			// 到达新阶段边界：切活动配置并执行 Dirichlet 吸附。
			r.startPhase(pi)
		}
		ph := &r.plan.Phases[pi]
		before := r.st.Cur.Clone()
		r.accountStep(before, r.acc)
		r.st.stepRowsParallel(workers)
		if err := r.st.commit(); err != nil {
			return reports, err
		}
		if r.st.Step >= ph.EndStep {
			rep := r.finalizePhase(r.acc)
			reports = append(reports, rep)
			// 新阶段入口（Dirichlet 吸附、累计器重置）在下一次迭代顶部完成。
		}
	}
	return reports, nil
}

// Step 返回已完成的绝对步号。
func (r *Runner) Step() int { return r.st.Step }

// Grid 返回当前场的拷贝。
func (r *Runner) Grid() *Grid { return r.st.Cur.Clone() }

// Save 落盘当前状态（含配置签名与阶段计划指纹）。
func (r *Runner) Save(path string) error {
	return SaveState(path, r.st.Step, r.st.Cur, r.cfg, r.plan)
}
