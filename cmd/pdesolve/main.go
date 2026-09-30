// pdesolve 是二维热传导方程的本地求解入口。
// 支持串行/并行推进、时间阶段计划、检查点落盘与断点续算、守恒量与逐阶段收支校验。
package main

import (
	"encoding/json"
	"flag"
	"fmt"
	"math"
	"os"
	"strings"

	"github.com/lisama569/pde-solver-go/solver"
)

// jsonBoundary 是边界条件的 JSON 表示：type ∈ dirichlet|neumann|periodic。
type jsonBoundary struct {
	Type  string  `json:"type"`
	Value float64 `json:"value"`
}

type jsonSource struct {
	Rate float64 `json:"rate"`
	X0   float64 `json:"x0"`
	X1   float64 `json:"x1"`
	Y0   float64 `json:"y0"`
	Y1   float64 `json:"y1"`
}

type jsonPhase struct {
	From    int          `json:"from"`
	To      int          `json:"to"`
	BC      jsonBC       `json:"bc"`
	Sources []jsonSource `json:"sources"`
}

type jsonBC struct {
	Left   jsonBoundary `json:"left"`
	Right  jsonBoundary `json:"right"`
	Top    jsonBoundary `json:"top"`
	Bottom jsonBoundary `json:"bottom"`
}

type jsonPlan struct {
	Phases []jsonPhase `json:"phases"`
}

func parseBoundary(b jsonBoundary) (solver.Boundary, error) {
	switch strings.ToLower(b.Type) {
	case "dirichlet":
		return solver.Boundary{Type: solver.Dirichlet, Value: b.Value}, nil
	case "neumann":
		return solver.Boundary{Type: solver.Neumann}, nil
	case "periodic":
		return solver.Boundary{Type: solver.Periodic}, nil
	}
	return solver.Boundary{}, fmt.Errorf("unknown boundary type %q", b.Type)
}

// loadPlan 读取 JSON 阶段计划文件并转换为 solver.PhasePlan。
func loadPlan(path string) (*solver.PhasePlan, error) {
	raw, err := os.ReadFile(path)
	if err != nil {
		return nil, err
	}
	var jp jsonPlan
	if err := json.Unmarshal(raw, &jp); err != nil {
		return nil, fmt.Errorf("parse %s: %w", path, err)
	}
	plan := &solver.PhasePlan{Phases: make([]solver.Phase, len(jp.Phases))}
	for k, p := range jp.Phases {
		ph := solver.Phase{From: p.From, To: p.To}
		for name, dst := range map[string]struct {
			jb  jsonBoundary
			set func(solver.Boundary)
		}{
			"left":   {p.BC.Left, func(b solver.Boundary) { ph.BC.Left = b }},
			"right":  {p.BC.Right, func(b solver.Boundary) { ph.BC.Right = b }},
			"top":    {p.BC.Top, func(b solver.Boundary) { ph.BC.Top = b }},
			"bottom": {p.BC.Bottom, func(b solver.Boundary) { ph.BC.Bottom = b }},
		} {
			_ = name
			b, err := parseBoundary(dst.jb)
			if err != nil {
				return nil, fmt.Errorf("phase %d: %w", k, err)
			}
			dst.set(b)
		}
		for _, s := range p.Sources {
			ph.Sources = append(ph.Sources, solver.Source{
				Rate: s.Rate, X0: s.X0, X1: s.X1, Y0: s.Y0, Y1: s.Y1,
			})
		}
		plan.Phases[k] = ph
	}
	return plan, nil
}

func main() {
	var (
		nx       = flag.Int("nx", 64, "网格列数（含边界）")
		ny       = flag.Int("ny", 48, "网格行数（含边界）")
		dx       = flag.Float64("dx", 0.01, "空间步长 x")
		dy       = flag.Float64("dy", 0.01, "空间步长 y")
		dt       = flag.Float64("dt", 1e-5, "时间步长")
		diff     = flag.Float64("D", 1.0, "扩散系数")
		steps    = flag.Int("steps", 100, "时间步数（带阶段计划时为任务总步数）")
		until    = flag.Int("until", 0, "只推进到该步号为止（0 = 跑完 steps；用于分段执行阶段计划）")
		workers  = flag.Int("workers", 1, "并行子域数（1 = 串行）")
		srcRate  = flag.Float64("src", 0, "中心矩形源项强度（0 = 无源，被 -plan 覆盖）")
		planPath = flag.String("plan", "", "JSON 阶段计划文件（为空则按单套边界/源项运行）")
		savePath = flag.String("save", "", "结束时保存检查点到该路径")
		resume   = flag.String("resume", "", "从该检查点文件恢复继续计算")
	)
	flag.Parse()

	cfg := &solver.Config{
		NX: *nx, NY: *ny,
		DX: *dx, DY: *dy, DT: *dt,
		Diffusivity: *diff,
		Steps:       *steps,
		BC: solver.Boundaries{
			Left:   solver.Boundary{Type: solver.Dirichlet, Value: 100},
			Right:  solver.Boundary{Type: solver.Neumann},
			Top:    solver.Boundary{Type: solver.Neumann},
			Bottom: solver.Boundary{Type: solver.Dirichlet, Value: 0},
		},
	}
	if *srcRate != 0 {
		cfg.Src = &solver.Source{Rate: *srcRate, X0: 0.4, X1: 0.6, Y0: 0.4, Y1: 0.6}
	}
	if *planPath != "" {
		plan, err := loadPlan(*planPath)
		if err != nil {
			fail(fmt.Errorf("config: %w", err))
		}
		cfg.Plan = plan
	}

	// 初始场：确定性高斯热点，或从检查点恢复
	var initial *solver.Grid
	startStep := 0
	if *resume != "" {
		cp, err := solver.LoadCheckpoint(*resume)
		if err != nil {
			fail(err)
		}
		// 恢复前校验：网格/步长/系数/阶段计划必须与保存时一致
		if err := cp.VerifyConfig(cfg); err != nil {
			fail(err)
		}
		initial = cp.Grid()
		startStep = cp.Step
		fmt.Printf("resumed from %s at step %d (%dx%d)\n", *resume, cp.Step, cp.NX, cp.NY)
	} else {
		initial = solver.NewGrid(cfg.NX, cfg.NY)
		for j := 0; j < cfg.NY; j++ {
			for i := 0; i < cfg.NX; i++ {
				x := float64(i)/float64(cfg.NX-1) - 0.5
				y := float64(j)/float64(cfg.NY-1) - 0.5
				initial.Set(i, j, 100*math.Exp(-(x*x+y*y)/0.02))
			}
		}
	}

	m0 := solver.TotalMass(initial, cfg.DX, cfg.DY)
	fmt.Printf("config: %dx%d dx=%g dy=%g dt=%g D=%g steps=%d workers=%d plan=%v (dt_limit=%g)\n",
		cfg.NX, cfg.NY, cfg.DX, cfg.DY, cfg.DT, cfg.Diffusivity, *steps, *workers,
		cfg.Plan != nil, cfg.StabilityLimit())

	// 带计划时 steps 为任务总步数，本次推进 [startStep, until)；否则推进 steps 步（旧行为）
	runSteps := *steps
	if cfg.Plan != nil {
		end := *until
		if end == 0 {
			end = *steps
		}
		if end > *steps || end < startStep {
			fail(fmt.Errorf("config: invalid run range [%d,%d) within total steps %d", startStep, end, *steps))
		}
		runSteps = end - startStep
	}
	final, reports, err := solver.RunPhases(cfg, initial, *workers, startStep, runSteps)
	if err != nil {
		fail(err)
	}
	endStep := startStep + runSteps

	rep := solver.CheckInvariant(final, cfg, m0, endStep)
	fmt.Printf("done: step=%d mass=%.6f (start %.6f, deviation %+.3e)\n",
		rep.Step, rep.Total, rep.Expected, rep.Deviation)

	// 逐阶段收支对账
	for i := range reports {
		r := &reports[i]
		status := "OK"
		if !r.Balanced() {
			status = "IMBALANCE"
		}
		fmt.Printf("phase %d steps[%d,%d): deltaMass=%+.6e srcTotal=%+.6e "+
			"boundary{L:%+.4e R:%+.4e T:%+.4e B:%+.4e} closure=%+.3e tol=%.1e %s\n",
			r.PhaseIndex, r.FromStep, r.ToStep, r.DeltaMass(), r.SourceTotal(),
			r.Boundary.Left, r.Boundary.Right, r.Boundary.Top, r.Boundary.Bottom,
			r.Closure(), r.ClosureTol(), status)
		for si, inj := range r.SourceInject {
			fmt.Printf("  phase %d source %d injected=%+.6e\n", r.PhaseIndex, si, inj)
		}
		if !r.Balanced() {
			os.Exit(1)
		}
	}

	if *savePath != "" {
		if err := solver.SaveCheckpointWithConfig(*savePath, endStep, final, cfg); err != nil {
			fail(err)
		}
		fmt.Printf("checkpoint saved to %s at step %d\n", *savePath, endStep)
	}
}

// fail 按错误类别打印可区分的失败原因并退出。
func fail(err error) {
	switch {
	case solver.IsConfigError(err):
		fmt.Fprintf(os.Stderr, "CONFIG ERROR: %v\n", err)
	case solver.IsStabilityError(err):
		fmt.Fprintf(os.Stderr, "STABILITY ERROR: %v\n", err)
	default:
		fmt.Fprintf(os.Stderr, "ERROR: %v\n", err)
	}
	os.Exit(1)
}
