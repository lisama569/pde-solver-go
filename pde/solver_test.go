package pde

import (
	"math"
	"testing"
)

// consistencyCfg is a scenario that exercises every stencil interaction:
// periodic in x (wrap across the subdomain cut), Dirichlet+Neumann in y,
// non-uniform initial field (seeded noise) and a localized source.
func consistencyCfg() Config {
	return Config{
		Nx:        24,
		Ny:        20,
		Dx:        1.0 / 24,
		Dy:        1.0 / 20,
		Alpha:     0.2,
		TotalTime: 0.05,
		Dt:        0.002, // rx = 0.2*0.002*576 = 0.2304, ry = 0.2*0.002*400 = 0.16
		Boundaries: Boundaries{
			X0: BoundaryCondition{Type: BoundaryPeriodic},
			X1: BoundaryCondition{Type: BoundaryPeriodic},
			Y0: BoundaryCondition{Type: BoundaryDirichlet, Value: 1},
			Y1: BoundaryCondition{Type: BoundaryNeumann, Value: 0},
		},
		Source: SourceConfig{
			Type: SourceGaussian, Scale: 3, Cx: 0.5, Cy: 0.5, Sigma: 0.12,
		},
		InitValue: 0.2, Seed: 42, NoiseAmp: 0.05, MaxValue: 1e6, ConserveTol: 1e-10,
	}
}

func fieldMaxAbsDiff(a, b *Field) (maxDiff float64, maxI, maxJ int) {
	for j := 0; j < a.Ny; j++ {
		for i := 0; i < a.Nx; i++ {
			d := math.Abs(a.at(i, j) - b.at(i, j))
			if d > maxDiff {
				maxDiff, maxI, maxJ = d, i, j
			}
		}
	}
	return maxDiff, maxI, maxJ
}

func runToEnd(t *testing.T, cfg Config, workers int) *Field {
	t.Helper()
	cfg.Workers = workers
	s, err := NewSolver(cfg)
	if err != nil {
		t.Fatalf("workers=%d: construction failed: %v", workers, err)
	}
	for s.StepIndex() < s.PlannedSteps() {
		if err := s.Step(workers); err != nil {
			t.Fatalf("workers=%d step %d failed: %v", workers, s.StepIndex(), err)
		}
	}
	return s.cur
}

// TestSerialParallelConsistency verifies that the domain-decomposed parallel
// update is bit-for-bit identical to the serial baseline, and independent of
// worker count and repeated runs. The double-buffered scheme should make the
// maximum difference exactly 0; a loose float64 tolerance is also asserted as
// the documented numerical guarantee.
func TestSerialParallelConsistency(t *testing.T) {
	cfg := consistencyCfg()
	if err := cfg.CheckStability(); err != nil {
		t.Fatalf("test config unstable: %v", err)
	}
	rx, ry := cfg.diffusionRatios()
	t.Logf("input: %s", cfg)
	t.Logf("stability basis: rx=%.6g ry=%.6g sum=%.6g <= 0.5; steps=%d", rx, ry, rx+ry,
		int(math.Round(cfg.TotalTime/cfg.Dt)))

	serial := runToEnd(t, cfg, 1)
	workerCounts := []int{2, 3, 4, 5, 7, 8, 24}
	// Numerical guarantee (documented tolerance).
	const tol = 1e-14
	for _, w := range workerCounts {
		f := runToEnd(t, cfg, w)
		d, i, j := fieldMaxAbsDiff(serial, f)
		t.Logf("workers=%2d vs serial(1): max|delta|=%.3g at (i=%d,j=%d) (tol=%.0e, bitwise-exact=%v)",
			w, d, i, j, tol, d == 0)
		if d > tol {
			t.Fatalf("workers=%d deviates from serial by %.6g > tol %.0e at (%d,%d)", w, d, tol, i, j)
		}
	}
	// Repeatability: same worker count, separate solver instances.
	a := runToEnd(t, cfg, 4)
	b := runToEnd(t, cfg, 4)
	d, i, j := fieldMaxAbsDiff(a, b)
	t.Logf("repeat run workers=4 vs workers=4: max|delta|=%.3g at (%d,%d)", d, i, j)
	if d != 0 {
		t.Fatalf("repeated parallel runs differ by %.6g", d)
	}
}

// TestPartitionRowsDeterministic checks the decomposition properties used to
// argue race freedom and order independence.
func TestPartitionRowsDeterministic(t *testing.T) {
	for _, nx := range []int{1, 2, 3, 7, 24, 25} {
		for _, n := range []int{1, 2, 3, 4, 8, 100} {
			p := partitionRows(nx, n)
			if p[0][0] != 0 || p[len(p)-1][1] != nx {
				t.Fatalf("nx=%d n=%d does not cover domain: %v", nx, n, p)
			}
			for k := 1; k < len(p); k++ {
				if p[k][0] != p[k-1][1] {
					t.Fatalf("nx=%d n=%d gap/overlap at %d: %v", nx, n, k, p)
				}
			}
		}
	}
	p1 := partitionRows(24, 4)
	p2 := partitionRows(24, 4)
	for k := range p1 {
		if p1[k] != p2[k] {
			t.Fatalf("partition not deterministic: %v vs %v", p1, p2)
		}
	}
	t.Logf("partition(24,4)=%v (contiguous, disjoint, deterministic)", p1)
}

// TestRuntimeBlowup verifies that an out-of-range value produced at runtime is
// rejected with KindBlowup (distinct from the a-priori KindStability check),
// naming the step and location. The source is chosen large enough that the
// finite MaxValue bound is crossed in the first few steps despite a stable
// time step.
func TestRuntimeBlowup(t *testing.T) {
	cfg := consistencyCfg()
	cfg.Source = SourceConfig{
		Type: SourceGaussian, Scale: 5e15, Cx: 0.5, Cy: 0.5, Sigma: 0.05,
	}
	cfg.MaxValue = 100
	if err := cfg.CheckStability(); err != nil {
		t.Fatal(err)
	}
	s, err := NewSolver(cfg)
	if err != nil {
		t.Fatal(err)
	}
	var blow *Failure
	for s.StepIndex() < s.PlannedSteps() {
		err := s.Step(3)
		if err != nil {
			blow = AsFailure(err)
			break
		}
	}
	if blow == nil {
		t.Fatal("expected runtime blowup failure, got none")
	}
	if blow.Kind != KindBlowup {
		t.Fatalf("want KindBlowup, got %s: %v", blow.Kind, blow)
	}
	if blow.Step <= 0 {
		t.Fatalf("blowup must name a positive step, got %d", blow.Step)
	}
	t.Logf("input: stable dt=%.3g source_scale=%.0e max_value=%g -> %s",
		cfg.Dt, cfg.Source.Scale, cfg.MaxValue, blow)

	// NaN must never be observable: all current field values stay finite.
	for j := 0; j < cfg.Ny; j++ {
		for i := 0; i < cfg.Nx; i++ {
			if v := s.cur.at(i, j); math.IsNaN(v) || math.IsInf(v, 0) {
				t.Fatalf("non-finite value leaked into field at (%d,%d)=%v", i, j, v)
			}
		}
	}
}

// TestMidRunWorkerSwitch advances one solver switching worker counts every
// few steps (including back to serial), and compares against a purely serial
// run. Results must stay bitwise identical.
func TestMidRunWorkerSwitch(t *testing.T) {
	cfg := consistencyCfg()
	serial := runToEnd(t, cfg, 1)
	s, err := NewSolver(cfg)
	if err != nil {
		t.Fatal(err)
	}
	pattern := []int{4, 1, 7, 2, 1, 3, 8, 1}
	k := 0
	for s.StepIndex() < s.PlannedSteps() {
		w := pattern[k%len(pattern)]
		if err := s.Step(w); err != nil {
			t.Fatal(err)
		}
		k++
	}
	d, i, j := fieldMaxAbsDiff(serial, s.cur)
	t.Logf("worker-switch pattern %v: max|delta vs serial|=%.3g at (%d,%d) over %d steps",
		pattern, d, i, j, s.StepIndex())
	if d != 0 {
		t.Fatalf("mid-run worker switching introduced deviation %.6g", d)
	}
}

// TestNeumannFluxHeatBalance specifically exercises a non-zero outward
// derivative: the predicted boundary influx must match the actual heat change
// to machine precision on a no-source, otherwise-insulated domain.
func TestNeumannFluxHeatBalance(t *testing.T) {
	cfg := Config{
		Nx: 10, Ny: 10, Dx: 0.1, Dy: 0.1, Alpha: 0.1,
		TotalTime: 0.02, Dt: 0.01, // rx+ry=0.2
		Boundaries: Boundaries{
			X0: BoundaryCondition{Type: BoundaryNeumann, Value: 1.0}, // du/dn=+1 outward -> heating
			X1: BoundaryCondition{Type: BoundaryNeumann, Value: 0},
			Y0: BoundaryCondition{Type: BoundaryNeumann, Value: 0},
			Y1: BoundaryCondition{Type: BoundaryNeumann, Value: 0},
		},
		InitValue: 0, MaxValue: 1e9, ConserveTol: 1e-10,
	}
	s, err := NewSolver(cfg)
	if err != nil {
		t.Fatal(err)
	}
	prev := s.cur.clone()
	if err := s.Step(3); err != nil {
		t.Fatal(err)
	}
	rep, err := CheckConservation(cfg, prev, s.cur, 1, 0, cfg.Dt, 3)
	if err != nil {
		t.Fatal(err)
	}
	// First step: interior is flat at 0, x0 ghosts = 0 + g*dx = 0.1.
	// Incoming heat = alpha*dt/dx * 0.1 * Ny * dy
	//               = 0.1*0.01/0.1*0.1 * 10*0.1 = 0.001
	want := 0.1 * 0.01 / 0.1 * 0.1 * float64(cfg.Ny) * cfg.Dy
	t.Logf("input: g_x0=%g -> flux=%.12g source=%.3g dH=%.12g rel_delta=%.3g (analytic flux %.12g)",
		cfg.Boundaries.X0.Value, rep.BoundaryFlux, rep.SourceHeat,
		rep.Heat-rep.PrevHeat, rep.RelDelta, want)
	if math.Abs(rep.BoundaryFlux-want) > 1e-14 {
		t.Fatalf("boundary flux %.12g != analytic %.12g", rep.BoundaryFlux, want)
	}
	if rep.RelDelta > 1e-12 {
		t.Fatalf("conservation residual too large: %.6g", rep.RelDelta)
	}
}
