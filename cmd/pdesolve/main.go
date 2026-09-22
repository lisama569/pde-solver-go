// pdesolve 是二维热传导方程的本地求解入口。
// 支持串行/并行推进、检查点落盘与断点续算、守恒量校验。
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
		steps    = flag.Int("steps", 100, "时间步数")
		workers  = flag.Int("workers", 1, "并行子域数（1 = 串行）")
		srcRate  = flag.Float64("src", 0, "中心矩形源项强度（0 = 无源）")
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

	// 初始场：确定性高斯热点，或从检查点恢复
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
