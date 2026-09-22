// Command pdesim runs a 2D heat-equation simulation with configurable
// boundary conditions, source term, parallel subdomain decomposition and
// checkpoint/resume.
//
// Usage:
//
//	pdesim run    [flags]   # run from t=0, optionally writing a checkpoint
//	pdesim resume [flags]   # resume from --checkpoint
//
// Failures exit non-zero with a typed message: config | stability |
// blowup | checkpoint | invariant.
package main

import (
	"flag"
	"fmt"
	"os"
	"strings"

	"github.com/lisama569/pde-solver-go/pde"
)

type cliOptions struct {
	cfg              pde.Config
	sourceScale      float64
	sourceCx         float64
	sourceCy         float64
	sourceSigma      float64
	bc               [4]string
	bcVal            [4]float64
	workers          int
	checkpoint       string
	checkpointEvery  int
	checkInvariant   bool
	invariantRegions int
	workersRun       int
}

func parseBC(s string) (pde.BoundaryCondition, error) {
	parts := strings.SplitN(s, "=", 2)
	name := strings.ToLower(strings.TrimSpace(parts[0]))
	var v float64
	if len(parts) == 2 {
		if _, err := fmt.Sscanf(strings.TrimSpace(parts[1]), "%g", &v); err != nil {
			return pde.BoundaryCondition{}, fmt.Errorf("cannot parse boundary value %q", parts[1])
		}
	}
	switch name {
	case "d", "dirichlet":
		return pde.BoundaryCondition{Type: pde.BoundaryDirichlet, Value: v}, nil
	case "n", "neumann":
		return pde.BoundaryCondition{Type: pde.BoundaryNeumann, Value: v}, nil
	case "p", "periodic":
		return pde.BoundaryCondition{Type: pde.BoundaryPeriodic}, nil
	default:
		return pde.BoundaryCondition{}, fmt.Errorf("unknown boundary kind %q (want dirichlet/neumann/periodic)", name)
	}
}

func registerFlags(fs *flag.FlagSet, o *cliOptions) {
	c := &o.cfg
	fs.IntVar(&c.Nx, "nx", 32, "interior grid points in x")
	fs.IntVar(&c.Ny, "ny", 32, "interior grid points in y")
	fs.Float64Var(&c.Dx, "dx", 1.0/32, "grid spacing in x")
	fs.Float64Var(&c.Dy, "dy", 1.0/32, "grid spacing in y")
	fs.Float64Var(&c.Alpha, "alpha", 0.1, "thermal diffusivity")
	fs.Float64Var(&c.TotalTime, "total-time", 0.05, "physical end time")
	fs.Float64Var(&c.Dt, "dt", 0.0, "time step (0 => 0.9 * stability limit)")
	fs.Float64Var(&c.InitValue, "init", 1.0, "uniform initial value")
	fs.Float64Var(&c.MaxValue, "max-value", pde.DefaultMaxValue, "finite sanity bound triggering blowup failure")
	fs.Float64Var(&c.ConserveTol, "conserve-tol", pde.DefaultConserveTol, "relative heat-balance tolerance")
	fs.BoolVar(&c.ConserveEnforce, "enforce-conservation", false, "fail the run if the heat balance exceeds --conserve-tol")
	fs.Int64Var(&c.Seed, "seed", 0, "deterministic perturbation seed (0 disables)")
	fs.Float64Var(&c.NoiseAmp, "noise", 0, "deterministic perturbation amplitude")
	fs.Float64Var(&o.sourceScale, "source-scale", 0, "source magnitude (0 disables source)")
	fs.StringVar(&o.bc[0], "x0", "d=0", "x=0 boundary: d=V dirichlet, n=G neumann, p periodic")
	fs.StringVar(&o.bc[1], "x1", "d=0", "x=Lx boundary")
	fs.StringVar(&o.bc[2], "y0", "d=0", "y=0 boundary")
	fs.StringVar(&o.bc[3], "y1", "d=0", "y=Ly boundary")
	fs.Float64Var(&o.sourceCx, "source-cx", 0.5, "gaussian source center x")
	fs.Float64Var(&o.sourceCy, "source-cy", 0.5, "gaussian source center y")
	fs.Float64Var(&o.sourceSigma, "source-sigma", 0.1, "gaussian source width")
	fs.IntVar(&o.workers, "workers", 1, "subdomain worker count (1 serial baseline)")
	fs.IntVar(&o.workersRun, "run-workers", 0, "parallelism at run time, defaults to --workers")
	fs.StringVar(&o.checkpoint, "checkpoint", "", "checkpoint file path")
	fs.IntVar(&o.checkpointEvery, "checkpoint-every", 0, "save checkpoint every N steps")
	fs.BoolVar(&o.checkInvariant, "check-invariant", false, "verify discrete heat balance each step")
	fs.IntVar(&o.invariantRegions, "invariant-regions", 0, "residual localization slabs (default: workers)")
}

func (o *cliOptions) buildConfig() (pde.Config, error) {
	c := o.cfg
	if o.sourceScale != 0 {
		c.Source = pde.SourceConfig{
			Type:  pde.SourceGaussian,
			Scale: o.sourceScale,
			Cx:    o.sourceCx,
			Cy:    o.sourceCy,
			Sigma: o.sourceSigma,
		}
	}
	for k, raw := range o.bc {
		bc, err := parseBC(raw)
		if err != nil {
			return c, err
		}
		switch k {
		case 0:
			c.Boundaries.X0 = bc
		case 1:
			c.Boundaries.X1 = bc
		case 2:
			c.Boundaries.Y0 = bc
		case 3:
			c.Boundaries.Y1 = bc
		}
	}
	// Auto dt: validate geometry/diffusivity first (which do not depend on
	// dt), then choose 0.9 * stability limit, then fully validate.
	if err := c.ValidateGeometryOnly(); err != nil {
		return c, err
	}
	if c.Dt == 0 {
		dt, _, err := c.AutoDt()
		if err != nil {
			return c, err
		}
		c.Dt = dt
	}
	if err := c.Validate(); err != nil {
		return c, err
	}
	if err := c.CheckStability(); err != nil {
		return c, err
	}
	c.Workers = o.workers
	return c, nil
}

func execute(cfg pde.Config, o cliOptions, resume bool) error {
	logf := func(m string) { fmt.Fprintln(os.Stderr, m) }
	logf("parameters: " + cfg.String())
	rx, ry := cfg.DiffusionRatios()
	logf(fmt.Sprintf("stability: rx=%.12g ry=%.12g sum=%.12g (limit <= 0.5, dt_limit=%.12g)",
		rx, ry, rx+ry, cfg.StabilityLimit()))

	opt := pde.RunOptions{
		Workers:          o.workersRun,
		Checkpoint:       o.checkpoint,
		Resume:           resume,
		CheckpointEvery:  o.checkpointEvery,
		CheckInvariant:   o.checkInvariant,
		InvariantRegions: o.invariantRegions,
		Logger:           logf,
	}
	res, err := pde.RunSimulation(cfg, opt)
	if err != nil {
		return err
	}
	logf(fmt.Sprintf("completed: steps=%d time=%.12g max_rel_heat_residual=%.3g",
		res.Steps, res.Time, res.MaxRelDelta))
	if res.LastInvariant != nil {
		logf(fmt.Sprintf("invariant final: heat=%.12g flux=%.12g source=%.12g abs_delta=%.3g rel_delta=%.3g",
			res.LastInvariant.Heat, res.LastInvariant.BoundaryFlux,
			res.LastInvariant.SourceHeat, res.LastInvariant.AbsDelta, res.LastInvariant.RelDelta))
	}
	if o.checkpoint != "" {
		logf("checkpoint: " + o.checkpoint)
	}
	fmt.Printf("steps=%d time=%.12g final_sum=%.12g\n", res.Steps, res.Time, res.Final.SumInterior())
	return nil
}

func run(args []string) error {
	if len(args) == 0 {
		return fmt.Errorf("subcommand required: run | resume")
	}
	cmd, rest := args[0], args[1:]
	switch cmd {
	case "run", "resume":
	default:
		return fmt.Errorf("unknown subcommand %q (want run | resume)", cmd)
	}
	fs := flag.NewFlagSet(cmd, flag.ContinueOnError)
	var o cliOptions
	registerFlags(fs, &o)
	if err := fs.Parse(rest); err != nil {
		return err
	}
	cfg, err := o.buildConfig()
	if err != nil {
		return err
	}
	return execute(cfg, o, cmd == "resume")
}

func main() {
	if err := run(os.Args[1:]); err != nil {
		if f := pde.AsFailure(err); f != nil {
			fmt.Fprintf(os.Stderr, "FAILURE kind=%s step=%d domain=%d: %s\n",
				f.Kind, f.Step, f.Domain, f.Msg)
		} else {
			fmt.Fprintln(os.Stderr, "error:", err)
		}
		os.Exit(1)
	}
}
