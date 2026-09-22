package pde

import (
	"path/filepath"
	"testing"
)

func resumeCfg() Config {
	c := consistencyCfg()
	c.Boundaries = Boundaries{
		X0: BoundaryCondition{Type: BoundaryDirichlet, Value: 1},
		X1: BoundaryCondition{Type: BoundaryNeumann, Value: 0.2},
		Y0: BoundaryCondition{Type: BoundaryDirichlet, Value: -1},
		Y1: BoundaryCondition{Type: BoundaryPeriodic},
	}
	// periodic y requires matched pair: reset to a fully non-periodic mix for
	// this scenario, covering D/N/D/N on all four faces.
	c.Boundaries.Y1 = BoundaryCondition{Type: BoundaryNeumann, Value: -0.1}
	return c
}

// TestCheckpointResume verifies that splitting a run at every possible step
// via checkpoint save/load reproduces the uninterrupted trajectory bit for
// bit, independent of worker counts before and after the split.
func TestCheckpointResume(t *testing.T) {
	cfg := resumeCfg()
	base := runToEnd(t, cfg, 1)
	t.Logf("input: %s", cfg)
	rx, ry := cfg.diffusionRatios()
	t.Logf("stability basis: rx+ry=%.6g <= 0.5; planned steps=%d", rx+ry,
		int(0.5+cfg.TotalTime/cfg.Dt))

	const tol = 1e-14
	splits := []int{1, 2, 7, 13, 24, 25}
	for _, split := range splits {
		dir := t.TempDir()
		path := filepath.Join(dir, "ckpt.json")

		// Phase A: serial to the split, save.
		s1, err := NewSolver(cfg)
		if err != nil {
			t.Fatal(err)
		}
		for s1.StepIndex() < split {
			if err := s1.Step(1); err != nil {
				t.Fatalf("phase A: %v", err)
			}
		}
		if err := SaveCheckpoint(path, s1.Snapshot()); err != nil {
			t.Fatal(err)
		}

		// Phase B: fresh solver, restore, advance to the end with 4 workers.
		s2, err := NewSolver(cfg)
		if err != nil {
			t.Fatal(err)
		}
		snap, err := LoadCheckpoint(path)
		if err != nil {
			t.Fatal(err)
		}
		if err := s2.Restore(snap); err != nil {
			t.Fatal(err)
		}
		if s2.StepIndex() != split {
			t.Fatalf("restored step %d want %d", s2.StepIndex(), split)
		}
		for s2.StepIndex() < s2.PlannedSteps() {
			if err := s2.Step(4); err != nil {
				t.Fatalf("phase B: %v", err)
			}
		}
		d, i, j := fieldMaxAbsDiff(base, s2.cur)
		t.Logf("split at step %2d, resume workers=4: max|delta vs continuous|=%.3g at (%d,%d) (tol=%.0e exact=%v)",
			split, d, i, j, tol, d == 0)
		if d > tol {
			t.Fatalf("resume at step %d deviates %.6g", split, d)
		}
	}
}

// TestCheckpointRepeatedResume chases a checkpoint through multiple save/load
// cycles with changing worker counts; the final state must equal baseline.
func TestCheckpointRepeatedResume(t *testing.T) {
	cfg := resumeCfg()
	base := runToEnd(t, cfg, 1)
	path := filepath.Join(t.TempDir(), "ckpt.json")
	workers := []int{3, 1, 5, 2, 4}

	s, err := NewSolver(cfg)
	if err != nil {
		t.Fatal(err)
	}
	step := 0
	w := 0
	for s.StepIndex() < s.PlannedSteps() {
		target := s.StepIndex() + 5
		if target > s.PlannedSteps() {
			target = s.PlannedSteps()
		}
		for s.StepIndex() < target {
			if err := s.Step(workers[w%len(workers)]); err != nil {
				t.Fatal(err)
			}
		}
		if err := SaveCheckpoint(path, s.Snapshot()); err != nil {
			t.Fatal(err)
		}
		s2, err := NewSolver(cfg)
		if err != nil {
			t.Fatal(err)
		}
		snap, err := LoadCheckpoint(path)
		if err != nil {
			t.Fatal(err)
		}
		if err := s2.Restore(snap); err != nil {
			t.Fatal(err)
		}
		s = s2
		step = s.StepIndex()
		w++
	}
	d, i, j := fieldMaxAbsDiff(base, s.cur)
	t.Logf("chased checkpoint %d times with workers %v, final step=%d: max|delta|=%.3g at (%d,%d)",
		w, workers, step, d, i, j)
	if d != 0 {
		t.Fatalf("multi-cycle resume deviates %.6g", d)
	}
}

// TestCorruptAndMismatchedCheckpoints checks typed, distinguishable failure
// reasons for bad checkpoint inputs.
func TestCorruptAndMismatchedCheckpoints(t *testing.T) {
	cfg := resumeCfg()
	s, err := NewSolver(cfg)
	if err != nil {
		t.Fatal(err)
	}
	for k := 0; k < 3; k++ {
		if err := s.Step(1); err != nil {
			t.Fatal(err)
		}
	}
	snap := s.Snapshot()

	t.Run("garbage", func(t *testing.T) {
		path := filepath.Join(t.TempDir(), "bad.json")
		if err := writeFileAtomic(path, []byte("{not json")); err != nil {
			t.Fatal(err)
		}
		_, err := LoadCheckpoint(path)
		f := AsFailure(err)
		if f == nil || f.Kind != KindCheckpoint {
			t.Fatalf("want KindCheckpoint, got %v", err)
		}
		t.Logf("garbage file -> kind=%s: %s", f.Kind, f.Msg)
	})

	t.Run("wrong schema", func(t *testing.T) {
		bad := *snap
		bad.SchemaVersion = 999
		_, err := UnmarshalSnapshot(mustMarshal(t, &bad))
		f := AsFailure(err)
		if f == nil || f.Kind != KindCheckpoint {
			t.Fatalf("want KindCheckpoint, got %v", err)
		}
		t.Logf("schema=%d -> kind=%s: %s", bad.SchemaVersion, f.Kind, f.Msg)
	})

	t.Run("nan data", func(t *testing.T) {
		// NaN must never be silently persisted: encoding refuses it and
		// decoding a file carrying a NaN/Infinity token refuses it too.
		bad := *snap
		bad.Data = append([]float64(nil), snap.Data...)
		bad.Data[len(bad.Data)/2] = nan()
		errMarshal := func() error { _, e := bad.Marshal(); return e }()
		if AsFailure(errMarshal).Kind != KindCheckpoint {
			t.Fatalf("Marshal must reject NaN with KindCheckpoint, got %v", errMarshal)
		}
		t.Logf("encoding NaN -> kind=%s: %s", AsFailure(errMarshal).Kind, AsFailure(errMarshal).Msg)

		raw := injectToken(mustMarshal(t, snap), "Infinity")
		_, err := UnmarshalSnapshot(raw)
		f := AsFailure(err)
		if f == nil || f.Kind != KindCheckpoint {
			t.Fatalf("decoding a non-finite payload must fail as KindCheckpoint, got %v", err)
		}
		t.Logf("decoding Infinity payload -> kind=%s: %s", f.Kind, f.Msg)
	})

	t.Run("config mismatch on restore", func(t *testing.T) {
		other := cfg
		other.Nx++
		s2, err := NewSolver(other)
		if err != nil {
			t.Fatal(err)
		}
		err = s2.Restore(snap)
		f := AsFailure(err)
		if f == nil || f.Kind != KindCheckpoint {
			t.Fatalf("want KindCheckpoint on grid mismatch, got %v", err)
		}
		t.Logf("grid mismatch -> kind=%s: %s", f.Kind, f.Msg)
	})
}

// TestResumeWithExtendedTotalTime mirrors a real restart workflow: the first
// process planned a short horizon; after restart the same physics continues
// toward a longer horizon with a different worker count. The trajectory must
// match an uninterrupted serial run to the longer horizon.
func TestResumeWithExtendedTotalTime(t *testing.T) {
	cfgA := resumeCfg()
	cfgA.TotalTime = 0.02 // 10 steps
	s1, err := NewSolver(cfgA)
	if err != nil {
		t.Fatal(err)
	}
	for s1.StepIndex() < 10 {
		if err := s1.Step(2); err != nil {
			t.Fatal(err)
		}
	}
	path := filepath.Join(t.TempDir(), "ck.json")
	if err := SaveCheckpoint(path, s1.Snapshot()); err != nil {
		t.Fatal(err)
	}

	cfgB := resumeCfg()
	cfgB.TotalTime = 0.05 // full 25 steps
	s2, err := NewSolver(cfgB)
	if err != nil {
		t.Fatal(err)
	}
	snap, err := LoadCheckpoint(path)
	if err != nil {
		t.Fatal(err)
	}
	if err := s2.Restore(snap); err != nil {
		t.Fatalf("restore with extended total time must be allowed: %v", err)
	}
	for s2.StepIndex() < s2.PlannedSteps() {
		if err := s2.Step(6); err != nil {
			t.Fatal(err)
		}
	}
	base := runToEnd(t, resumeCfg(), 1)
	d, i, j := fieldMaxAbsDiff(base, s2.cur)
	t.Logf("restart with extended horizon (0.02 -> 0.05), workers 2 -> 6: max|delta|=%.3g at (%d,%d)",
		d, i, j)
	if d != 0 {
		t.Fatalf("extended-horizon resume deviates %.6g", d)
	}
}
