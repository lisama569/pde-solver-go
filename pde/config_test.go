package pde

import (
	"math"
	"testing"
)

func validCfg() Config {
	return Config{
		Nx: 12, Ny: 12, Dx: 0.1, Dy: 0.1, Alpha: 0.1,
		TotalTime: 0.1, Dt: 0.02,
		Boundaries: Boundaries{
			X0: BoundaryCondition{Type: BoundaryDirichlet, Value: 1},
			X1: BoundaryCondition{Type: BoundaryDirichlet, Value: 0},
			Y0: BoundaryCondition{Type: BoundaryNeumann, Value: 0},
			Y1: BoundaryCondition{Type: BoundaryNeumann, Value: 0},
		},
		Source:  SourceConfig{Type: SourceConstant, Scale: 1},
		Workers: 1, InitValue: 0.5, MaxValue: 1e6, ConserveTol: 1e-10,
	}
}

func TestValidConfig(t *testing.T) {
	cfg := validCfg()
	if err := cfg.Validate(); err != nil {
		t.Fatalf("expected valid config, got %v", err)
	}
	if err := cfg.CheckStability(); err != nil {
		t.Fatalf("expected stable, got %v", err)
	}
	rx, ry := cfg.diffusionRatios()
	t.Logf("input: %s", cfg)
	t.Logf("stability basis: rx=%.12g ry=%.12g sum=%.12g limit=0.5 dt_limit=%.12g",
		rx, ry, rx+ry, cfg.stabilityLimit())
}

func TestIllegalConfig(t *testing.T) {
	cases := []struct {
		name   string
		mutate func(*Config)
		want   FailureKind
	}{
		{"small grid", func(c *Config) { c.Nx = 1 }, KindConfig},
		{"zero dx", func(c *Config) { c.Dx = 0 }, KindConfig},
		{"negative dy", func(c *Config) { c.Dy = -1 }, KindConfig},
		{"negative alpha", func(c *Config) { c.Alpha = -0.1 }, KindConfig},
		{"nan dt", func(c *Config) { c.Dt = math.NaN() }, KindConfig},
		{"zero total time", func(c *Config) { c.TotalTime = 0 }, KindConfig},
		{"dt does not divide total", func(c *Config) { c.Dt = 0.03; c.TotalTime = 0.1 }, KindConfig},
		{"unmatched periodic x", func(c *Config) {
			c.Boundaries.X0 = BoundaryCondition{Type: BoundaryPeriodic}
		}, KindConfig},
		{"nan dirichlet value", func(c *Config) {
			c.Boundaries.X1 = BoundaryCondition{Type: BoundaryDirichlet, Value: math.NaN()}
		}, KindConfig},
		{"bad gaussian sigma", func(c *Config) {
			c.Source = SourceConfig{Type: SourceGaussian, Scale: 1, Sigma: 0}
		}, KindConfig},
		{"negative workers", func(c *Config) { c.Workers = -2 }, KindConfig},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			cfg := validCfg()
			tc.mutate(&cfg)
			err := NewSolverOrValidate(cfg)
			f := AsFailure(err)
			if f == nil {
				t.Fatalf("input %+v: expected failure, got nil", cfg)
			}
			if f.Kind != tc.want {
				t.Fatalf("input %+v: want kind=%s got kind=%s: %v", cfg, tc.want, f.Kind, f)
			}
			t.Logf("input: Nx=%d Ny=%d dx=%g dy=%g alpha=%g dt=%g T=%g -> rejected kind=%s: %s",
				cfg.Nx, cfg.Ny, cfg.Dx, cfg.Dy, cfg.Alpha, cfg.Dt, cfg.TotalTime, f.Kind, f.Msg)
		})
	}
}

func TestStabilityRejection(t *testing.T) {
	// Stable baseline: rx+ry = 0.4. Increase dt beyond the limit.
	cfg := validCfg()
	rx0, ry0 := cfg.diffusionRatios()
	if got := rx0 + ry0; got > 0.5 {
		t.Fatalf("baseline unexpectedly unstable sum=%g", got)
	}
	cfg.TotalTime = 0.06 // dt divides T (2 steps); instability is the only issue
	cfg.Dt = 0.03        // rx+ry = 0.6 > 0.5
	if err := cfg.Validate(); err != nil {
		t.Fatalf("config should be legal (instability is a distinct failure): %v", err)
	}
	err := cfg.CheckStability()
	f := AsFailure(err)
	if f == nil || f.Kind != KindStability {
		t.Fatalf("want KindStability, got %v", err)
	}
	t.Logf("input: dt=%g rx+ry=%g dt_limit=%.12g -> %s",
		cfg.Dt, 0.6, cfg.stabilityLimit(), f.Msg)
	if _, err := NewSolver(cfg); AsFailure(err).Kind != KindStability {
		t.Fatalf("NewSolver must propagate KindStability, got %v", err)
	}
}

func TestStabilityIsDistinguishableFromConfig(t *testing.T) {
	// A dt that is legal and divides T but violates CFL is a stability
	// failure, not a config failure.
	cfg := validCfg()
	cfg.TotalTime = 0.03
	cfg.Dt = 0.03 // single step, rx+ry=0.6
	if err := cfg.Validate(); err != nil {
		t.Fatalf("Validate should pass: %v", err)
	}
	err := cfg.CheckStability()
	if AsFailure(err).Kind != KindStability {
		t.Fatalf("want stability error, got %v", err)
	}
	t.Logf("distinguishable failures: legal-but-unstable dt=%.3g reported as %s: %s",
		cfg.Dt, AsFailure(err).Kind, err)
}

func TestAutoDtStableAndDivisible(t *testing.T) {
	cfg := validCfg()
	cfg.Dt = 0
	dt, n, err := cfg.AutoDt()
	if err != nil {
		t.Fatal(err)
	}
	cfg.Dt = dt
	if err := cfg.CheckStability(); err != nil {
		t.Fatalf("auto dt unstable: %v", err)
	}
	if math.Abs(dt*float64(n)-cfg.TotalTime) > 1e-12*cfg.TotalTime {
		t.Fatalf("auto dt not divisible: dt=%g n=%d T=%g", dt, n, cfg.TotalTime)
	}
	t.Logf("auto dt=%.15g n=%d (dt_limit=%.12g, rx+ry=%.6g)",
		dt, n, cfg.stabilityLimit(), func() float64 {
			a, b := cfg.diffusionRatios()
			return a + b
		}())
}
