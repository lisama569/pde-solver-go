package pde

import (
	"math"
	"path/filepath"
	"testing"
)

// TestConservationAllBoundaryKinds verifies the exact discrete heat balance
// H1-H0 = boundary influx + source injection for each boundary-condition
// combination, and that the residual localizes per slab.
func TestConservationAllBoundaryKinds(t *testing.T) {
	cases := []struct {
		name string
		bc   Boundaries
		src  SourceConfig
	}{
		{"all-dirichlet", Boundaries{
			X0: BoundaryCondition{Type: BoundaryDirichlet, Value: 1},
			X1: BoundaryCondition{Type: BoundaryDirichlet, Value: 0},
			Y0: BoundaryCondition{Type: BoundaryDirichlet, Value: 0.5},
			Y1: BoundaryCondition{Type: BoundaryDirichlet, Value: 0},
		}, SourceConfig{Type: SourceConstant, Scale: 1.5}},
		{"all-neumann", Boundaries{
			X0: BoundaryCondition{Type: BoundaryNeumann, Value: 0},
			X1: BoundaryCondition{Type: BoundaryNeumann, Value: 0},
			Y0: BoundaryCondition{Type: BoundaryNeumann, Value: 0},
			Y1: BoundaryCondition{Type: BoundaryNeumann, Value: 0},
		}, SourceConfig{Type: SourceGaussian, Scale: 4, Cx: 0.3, Cy: 0.7, Sigma: 0.1}},
		{"periodic-x-dn-y", Boundaries{
			X0: BoundaryCondition{Type: BoundaryPeriodic},
			X1: BoundaryCondition{Type: BoundaryPeriodic},
			Y0: BoundaryCondition{Type: BoundaryDirichlet, Value: 1},
			Y1: BoundaryCondition{Type: BoundaryNeumann, Value: 0.3},
		}, SourceConfig{Type: SourceConstant, Scale: -2}},
		{"fully-periodic", Boundaries{
			X0: BoundaryCondition{Type: BoundaryPeriodic},
			X1: BoundaryCondition{Type: BoundaryPeriodic},
			Y0: BoundaryCondition{Type: BoundaryPeriodic},
			Y1: BoundaryCondition{Type: BoundaryPeriodic},
		}, SourceConfig{Type: SourceNone}},
	}
	const tol = 1e-12
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			cfg := Config{
				Nx: 16, Ny: 13, Dx: 0.125, Dy: 0.1, Alpha: 0.15,
				TotalTime: 0.018, Dt: 0.003, // rx=0.0288, ry=0.045 -> 0.0738
				Boundaries: tc.bc, Source: tc.src,
				InitValue: 1.0, Seed: 7, NoiseAmp: 0.02,
				MaxValue: 1e8, ConserveTol: tol,
			}
			s, err := NewSolver(cfg)
			if err != nil {
				t.Fatal(err)
			}
			nRegions := 4
			var worst float64
			for s.StepIndex() < s.PlannedSteps() {
				prev := s.cur.clone()
				pt := s.Time()
				if err := s.Step(nRegions); err != nil {
					t.Fatal(err)
				}
				rep, err := CheckConservation(cfg, prev, s.cur, s.StepIndex(), pt, cfg.Dt, nRegions)
				if err != nil {
					t.Fatalf("step %d: %v", s.StepIndex(), err)
				}
				if rep.RelDelta > worst {
					worst = rep.RelDelta
				}
				if len(rep.RegionDeltas) != nRegions {
					t.Fatalf("want %d region deltas, got %d", nRegions, len(rep.RegionDeltas))
				}
				// Per-region residuals must sum to the global residual.
				var sum float64
				for _, r := range rep.RegionDeltas {
					sum += r.Delta
				}
				if math.Abs(sum-(rep.Heat-rep.PrevHeat-rep.SourceHeat-rep.BoundaryFlux)) > tol*math.Max(1, math.Abs(rep.Heat)) {
					t.Fatalf("region residuals do not decompose global residual: sum=%.6g global=%.6g",
						sum, rep.Heat-rep.PrevHeat-rep.SourceHeat-rep.BoundaryFlux)
				}
				if s.StepIndex() == 1 {
					t.Logf("input: %s", cfg)
					t.Logf("step 1: heat=%.10g dH=%.6g flux=%.6g source=%.6g abs_delta=%.3g rel_delta=%.3g",
						rep.Heat, rep.Heat-rep.PrevHeat, rep.BoundaryFlux, rep.SourceHeat,
						rep.AbsDelta, rep.RelDelta)
					for k, r := range rep.RegionDeltas {
						t.Logf("  region %d i=[%d,%d) delta=%.3g", k, r.I0, r.I1, r.Delta)
					}
				}
			}
			t.Logf("%s: worst relative residual over all steps = %.3g (tol %.0e)", tc.name, worst, tol)
		})
	}
}

// TestInvariantViolationLocalized proves an externally corrupted state that
// breaks the balance is reported as KindInvariant with the offending slab
// named, distinctly from blowup/config/stability.
func TestInvariantViolationLocalized(t *testing.T) {
	cfg := Config{
		Nx: 20, Ny: 10, Dx: 0.1, Dy: 0.1, Alpha: 0.1,
		TotalTime: 0.02, Dt: 0.01, // rx+ry = 0.2
		Boundaries: Boundaries{
			X0: BoundaryCondition{Type: BoundaryNeumann, Value: 0},
			X1: BoundaryCondition{Type: BoundaryNeumann, Value: 0},
			Y0: BoundaryCondition{Type: BoundaryNeumann, Value: 0},
			Y1: BoundaryCondition{Type: BoundaryNeumann, Value: 0},
		},
		InitValue: 1, MaxValue: 1e9, ConserveTol: 1e-10,
	}
	s, err := NewSolver(cfg)
	if err != nil {
		t.Fatal(err)
	}
	prev := s.cur.clone()
	if err := s.Step(1); err != nil {
		t.Fatal(err)
	}
	// Sabotage interior cells in slab i in [10,15) by adding a non-physical
	// injection of heat; conservation check must flag that region.
	for j := 0; j < cfg.Ny; j++ {
		for i := 10; i < 15; i++ {
			s.cur.set(i, j, s.cur.at(i, j)+5.0)
		}
	}
	rep, err := CheckConservation(cfg, prev, s.cur, 1, 0, cfg.Dt, 4)
	f := AsFailure(err)
	if f == nil || f.Kind != KindInvariant {
		t.Fatalf("want KindInvariant, got %v", err)
	}
	t.Logf("input: injected +5 at i=[10,15) -> kind=%s step=%d rel_delta=%.3g: %s",
		f.Kind, f.Step, rep.RelDelta, f.Msg)
	// Expected localized delta = 5 * 25 * dx*dy = 1.25 (plus roundoff).
	var regionSum float64
	for _, r := range rep.RegionDeltas {
		if r.I0 == 10 && r.I1 == 15 {
			regionSum = r.Delta
		}
	}
	if math.Abs(regionSum-2.5) > 1e-9 {
		t.Fatalf("offending slab delta=%.6g want 2.5; regions=%v", regionSum, rep.RegionDeltas)
	}
	t.Logf("localized offending slab i=[10,15) delta=%.6g (expected 2.5 = +5 * %d cells * dx*dy=%g)",
		regionSum, 5*cfg.Ny, cfg.Dx*cfg.Dy)
}

// TestRunnerResumeEndToEnd drives RunSimulation twice (split mid-way through
// a checkpoint) and requires bitwise identical final fields and heat.
func TestRunnerResumeEndToEnd(t *testing.T) {
	cfg := resumeCfg()
	cfg.ConserveEnforce = true
	dir := t.TempDir()
	path := filepath.Join(dir, "state.json")

	full, err := RunSimulation(cfg, RunOptions{
		Workers: 2, Checkpoint: path, CheckInvariant: true, InvariantRegions: 2,
		Logger: func(m string) { t.Logf("[full] %s", m) },
	})
	if err != nil {
		t.Fatal(err)
	}

	// Split: run 10 steps via periodic checkpoints, discard solver, resume.
	const split = 10
	s1, err := NewSolver(cfg)
	if err != nil {
		t.Fatal(err)
	}
	for s1.StepIndex() < split {
		if err := s1.Step(1); err != nil {
			t.Fatal(err)
		}
	}
	if err := SaveCheckpoint(path, s1.Snapshot()); err != nil {
		t.Fatal(err)
	}
	resumed, err := RunSimulation(cfg, RunOptions{
		Workers: 3, Checkpoint: path, Resume: true,
		CheckInvariant: true, InvariantRegions: 2,
		Logger: func(m string) { t.Logf("[resume] %s", m) },
	})
	if err != nil {
		t.Fatal(err)
	}
	d, i, j := fieldMaxAbsDiff(full.Final, resumed.Final)
	t.Logf("end-to-end: continuous workers=2 vs split@%d+resume workers=3: max|delta|=%.3g at (%d,%d); max heat rel residual full=%.3g resumed=%.3g",
		split, d, i, j, full.MaxRelDelta, resumed.MaxRelDelta)
	if d != 0 {
		t.Fatalf("resumed trajectory differs: %.6g", d)
	}
	if resumed.Steps != full.Steps {
		t.Fatalf("step count %d != %d", resumed.Steps, full.Steps)
	}
}

// TestRunnerInvariantHardFailure verifies ConserveEnforce-style usage: callers
// can make a conservation mismatch fatal and the failure is typed.
func TestRunnerInvariantHardFailure(t *testing.T) {
	cfg := resumeCfg()
	cfg.ConserveTol = 1e-30 // absurdly tight; roundoff (~1e-16) must exceed it
	_, err := RunSimulation(cfg, RunOptions{Workers: 1, CheckInvariant: true, InvariantRegions: 3})
	f := AsFailure(err)
	if f == nil || f.Kind != KindInvariant {
		t.Fatalf("want KindInvariant under impossible tolerance, got %v", err)
	}
	t.Logf("impossible tol %.0e at step %d -> kind=%s: %s", cfg.ConserveTol, f.Step, f.Kind, f.Msg)
}
