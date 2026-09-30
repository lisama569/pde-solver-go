// pdesolve 是二维热传导方程的本地求解入口。
// 支持串行/并行推进、检查点落盘与断点续算、守恒量校验，
// 以及按时间阶段切换源项与四边边界（-plan）。
package main

import (
	"flag"
	"fmt"
	"math"
	"os"

	"github.com/lisama569/pde-solver-go/solver"
)

func main() {
	var (
		nx       = flag.Int("nx", 64, "网格列数（含边界）")
		ny       = flag.Int("ny", 48, "网格行数（含边界）")
		dx       = flag.Float64("dx", 0.01, "空间步长 x")
		dy       = flag.Float64("dy", 0.01, "空间步长 y")
		dt       = flag.Float64("dt", 1e-5, "时间步长")
		diff     = flag.Float64("D", 1.0, "扩散系数")
		steps    = flag.Int("steps", 100, "目标绝对时间步号（要跑到第几步）")
		workers  = flag.Int("workers", 1, "并行子域数（1 = 串行）")
		srcRate  = flag.Float64("src", 0, "中心矩形源项强度（0 = 无源；仅旧单套配置模式）")
		planPath = flag.String("plan", "", "阶段计划 JSON 文件路径；设置后启用按阶段切换边界/源项")
		savePath = flag.String("save", "", "结束时保存状态（v2，含配置签名）到该路径")
		resume   = flag.String("resume", "", "从该状态文件恢复继续计算（旧 v1 检查点请用旧版本流程）")
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

	var plan *solver.PhasePlan
	if *planPath != "" {
		raw, err := os.ReadFile(*planPath)
		if err != nil {
			fail(err)
		}
		plan, err = solver.ParsePhasePlanJSON(raw)
		if err != nil {
			fail(err)
		}
	} else if *srcRate != 0 {
		cfg.Src = &solver.Source{Rate: *srcRate, X0: 0.4, X1: 0.6, Y0: 0.4, Y1: 0.6}
	}

	// 阶段计划模式：带配置签名的状态保存/恢复 + 阶段末对账。
	if plan != nil {
		runPhased(cfg, plan, *workers, *resume, *savePath)
		return
	}

	// 旧单套配置模式（行为与历史版本一致）。
	var initial *solver.Grid
	startStep := 0
	if *resume != "" {
		cp, err := solver.LoadCheckpoint(*resume)
		if err != nil {
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
	fmt.Printf("config: %dx%d dx=%g dy=%g dt=%g D=%g steps=%d workers=%d (dt_limit=%g)\n",
		cfg.NX, cfg.NY, cfg.DX, cfg.DY, cfg.DT, cfg.Diffusivity, *steps, *workers, cfg.StabilityLimit())

	final, err := solver.RunParallel(cfg, initial, *workers, *steps)
	if err != nil {
		fail(err)
	}

	rep := solver.CheckInvariant(final, cfg, m0, startStep+*steps)
	fmt.Printf("done: step=%d mass=%.6f (initial %.6f, deviation %+.3e)\n",
		rep.Step, rep.Total, rep.Expected, rep.Deviation)

	if *savePath != "" {
		if err := solver.SaveCheckpoint(*savePath, startStep+*steps, final); err != nil {
			fail(err)
		}
		fmt.Printf("checkpoint saved to %s at step %d\n", *savePath, startStep+*steps)
	}
}

// runPhased 执行阶段计划模式：恢复时严格校验配置签名，结束时打印每阶段对账。
func runPhased(cfg *solver.Config, plan *solver.PhasePlan, workers int, resumePath, savePath string) {
	var r *solver.Runner
	if resumePath != "" {
		st, err := solver.LoadState(resumePath)
		if err != nil {
			fail(err)
		}
		// 恢复构造内部做签名比对：网格/步长/系数/计划任一不一致都会在推进前拒绝。
		r, err = solver.ResumeRunner(cfg, plan, st)
		if err != nil {
			fail(err)
		}
		fmt.Printf("resumed phased run from %s at step %d\n", resumePath, st.Step)
	} else {
		initial := solver.NewGrid(cfg.NX, cfg.NY)
		for j := 0; j < cfg.NY; j++ {
			for i := 0; i < cfg.NX; i++ {
				x := float64(i)/float64(cfg.NX-1) - 0.5
				y := float64(j)/float64(cfg.NY-1) - 0.5
				initial.Set(i, j, 100*math.Exp(-(x*x+y*y)/0.02))
			}
		}
		var err error
		r, err = solver.NewRunner(cfg, plan, initial)
		if err != nil {
			fail(err)
		}
	}

	remain := cfg.Steps - r.Step()
	fmt.Printf("phased config: %dx%d dx=%g dy=%g dt=%g D=%g targetStep=%d workers=%d phases=%d (dt_limit=%g)\n",
		cfg.NX, cfg.NY, cfg.DX, cfg.DY, cfg.DT, cfg.Diffusivity, cfg.Steps, workers, plan.Len(), cfg.StabilityLimit())
	fmt.Printf("balance tolerance: abs=%g rel=%g (tol=abs+rel*scale)\n",
		solver.DefaultBalanceTolAbs, solver.DefaultBalanceTolRel)

	reps, err := r.Advance(remain, workers)
	if err != nil {
		fail(err)
	}
	printReports(reps)

	fmt.Printf("done: step=%d mass=%.6f\n", r.Step(), solver.TotalMass(r.Grid(), cfg.DX, cfg.DY))
	if savePath != "" {
		if err := r.Save(savePath); err != nil {
			fail(err)
		}
		fmt.Printf("state saved to %s at step %d\n", savePath, r.Step())
	}
}

// printReports 打印每阶段收支对账：注入、四边通量、重置跳变、实测残差与判定。
func printReports(reps []solver.PhaseReport) {
	for _, rep := range reps {
		fmt.Printf("phase %d [%d,%d): mass %.6f -> %.6f (change %+.6e), dirichletReset=%+.3e\n",
			rep.PhaseIndex, rep.StartStep, rep.EndStep, rep.StartMass, rep.EndMass, rep.MassChange, rep.BoundaryReset)
		for _, src := range rep.Sources {
			fmt.Printf("  source[%s] injection=%+.9f\n", src.Name, src.Rate)
		}
		flow := 0.0
		for _, e := range rep.Edges {
			fmt.Printf("  edge[%s] netInflow=%+.9f\n", e.Edge, e.Flow)
			flow += e.Flow
		}
		status := "CLOSED"
		if !rep.Closed {
			status = "NOT-CLOSED"
		}
		fmt.Printf("  balance: residual=%+.3e tolerance=%+.3e -> %s\n", rep.Residual, rep.Tolerance, status)
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
		// 文件缺失/格式不符/JSON 解析等启动前问题归入配置类，绝不与失稳混淆。
		fmt.Fprintf(os.Stderr, "CONFIG ERROR: %v\n", err)
	}
	os.Exit(1)
}
