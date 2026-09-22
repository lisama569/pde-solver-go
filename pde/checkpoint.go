package pde

import (
	"encoding/json"
	"fmt"
	"os"
	"path/filepath"
)

// Snapshot is the full resumable state of a simulation.
type Snapshot struct {
	SchemaVersion int       `json:"schema_version"`
	Config        Config    `json:"config"`
	Step          int       `json:"step"`
	Time          float64   `json:"time"`
	Data          []float64 `json:"data"` // interior values in row-major order
}

// Marshal encodes a snapshot as indented, deterministic JSON.
func (s *Snapshot) Marshal() ([]byte, error) {
	b, err := json.MarshalIndent(s, "", "  ")
	if err != nil {
		return nil, failf(KindCheckpoint, s.Step, -1, "cannot encode checkpoint: %v", err)
	}
	return append(b, '\n'), nil
}

// UnmarshalSnapshot decodes and structurally validates a snapshot.
func UnmarshalSnapshot(b []byte) (*Snapshot, error) {
	var s Snapshot
	if err := json.Unmarshal(b, &s); err != nil {
		return nil, failf(KindCheckpoint, -1, -1, "cannot decode checkpoint: %v", err)
	}
	if s.SchemaVersion != SchemaVersion {
		return nil, failf(KindCheckpoint, s.Step, -1,
			"unsupported checkpoint schema version %d (want %d)", s.SchemaVersion, SchemaVersion)
	}
	if err := s.Config.Validate(); err != nil {
		return nil, err
	}
	if err := s.Config.CheckStability(); err != nil {
		return nil, err
	}
	if s.Step < 0 || s.Time < 0 {
		return nil, failf(KindCheckpoint, s.Step, -1, "negative step %d or time %v", s.Step, s.Time)
	}
	if len(s.Data) != s.Config.Nx*s.Config.Ny {
		return nil, failf(KindCheckpoint, s.Step, -1,
			"checkpoint data length %d does not match grid %dx%d", len(s.Data), s.Config.Nx, s.Config.Ny)
	}
	for k, v := range s.Data {
		if !finite(v) {
			return nil, failf(KindCheckpoint, s.Step, -1,
				"non-finite value %v at data index %d (i=%d,j=%d)", v, k, k%s.Config.Nx, k/s.Config.Nx)
		}
	}
	wantTime := float64(s.Step) * s.Config.Dt
	if d := s.Time - wantTime; d < -1e-9*s.Config.Dt || d > 1e-9*s.Config.Dt {
		return nil, failf(KindCheckpoint, s.Step, -1,
			"checkpoint time %.17g inconsistent with step %d * dt = %.17g", s.Time, s.Step, wantTime)
	}
	return &s, nil
}

// trajectoryEqual reports whether the physical trajectory-defining fields of
// two configs coincide. Worker count is an execution strategy, not part of
// trajectory identity (serial and parallel updates are bitwise identical by
// design), so Workers is intentionally excluded; MaxValue and ConserveTol are
// monitoring parameters that do not change the trajectory either.
func trajectoryEqual(a, b Config) string {
	switch {
	case a.Nx != b.Nx || a.Ny != b.Ny:
		return fmt.Sprintf("grid %dx%d vs %dx%d", a.Nx, a.Ny, b.Nx, b.Ny)
	case a.Dx != b.Dx || a.Dy != b.Dy:
		return fmt.Sprintf("spacing dx/dy %g/%g vs %g/%g", a.Dx, a.Dy, b.Dx, b.Dy)
	case a.Alpha != b.Alpha:
		return fmt.Sprintf("alpha %g vs %g", a.Alpha, b.Alpha)
	case a.Dt != b.Dt:
		return fmt.Sprintf("dt %g vs %g", a.Dt, b.Dt)
	case a.Boundaries != b.Boundaries:
		return "boundary conditions differ"
	case a.Source != b.Source:
		return "source configuration differs"
	case a.InitValue != b.InitValue:
		return fmt.Sprintf("init_value %g vs %g", a.InitValue, b.InitValue)
	case a.Seed != b.Seed:
		return fmt.Sprintf("seed %d vs %d", a.Seed, b.Seed)
	case a.NoiseAmp != b.NoiseAmp:
		return fmt.Sprintf("noise_amp %g vs %g", a.NoiseAmp, b.NoiseAmp)
	}
	return ""
}

// checkCompatible checks that a snapshot can resume a solver built from cfg.
// Grid/physics mismatches are rejected rather than silently interpolated.
func (s *Snapshot) checkCompatible(cfg Config) error {
	if s.SchemaVersion != SchemaVersion {
		return failf(KindCheckpoint, s.Step, -1,
			"checkpoint schema %d != current %d", s.SchemaVersion, SchemaVersion)
	}
	if reason := trajectoryEqual(s.Config, cfg); reason != "" {
		return failf(KindCheckpoint, s.Step, -1,
			"checkpoint trajectory does not match solver: %s (worker count may differ; physics must match)", reason)
	}
	if s.Step < 0 {
		return failf(KindCheckpoint, s.Step, -1, "checkpoint step %d must be non-negative", s.Step)
	}
	totalSteps, err := cfg.Steps()
	if err != nil {
		return failf(KindCheckpoint, s.Step, -1, "solver config cannot determine total steps: %v", err)
	}
	if s.Step > totalSteps {
		return failf(KindCheckpoint, s.Step, -1,
			"checkpoint step %d is beyond solver total steps %d", s.Step, totalSteps)
	}
	if len(s.Data) != cfg.Nx*cfg.Ny {
		return failf(KindCheckpoint, s.Step, -1, "data length %d does not match grid", len(s.Data))
	}
	return nil
}

// writeFileAtomic writes data to path via a temp file in the same directory
// and an atomic rename, so a crashed write cannot corrupt a prior checkpoint.
func writeFileAtomic(path string, data []byte) error {
	dir := filepath.Dir(path)
	tmp, err := os.CreateTemp(dir, ".checkpoint-*.tmp")
	if err != nil {
		return failf(KindCheckpoint, -1, -1, "cannot create temp checkpoint in %s: %v", dir, err)
	}
	tmpName := tmp.Name()
	cleanup := func() { _ = os.Remove(tmpName) }
	if _, err := tmp.Write(data); err != nil {
		tmp.Close()
		cleanup()
		return failf(KindCheckpoint, -1, -1, "cannot write checkpoint temp file: %v", err)
	}
	if err := tmp.Sync(); err != nil {
		tmp.Close()
		cleanup()
		return failf(KindCheckpoint, -1, -1, "cannot flush checkpoint: %v", err)
	}
	if err := tmp.Close(); err != nil {
		cleanup()
		return failf(KindCheckpoint, -1, -1, "cannot close checkpoint: %v", err)
	}
	if err := os.Rename(tmpName, path); err != nil {
		cleanup()
		return failf(KindCheckpoint, -1, -1, "cannot rename checkpoint into place: %v", err)
	}
	return nil
}

// SaveCheckpoint atomically writes the snapshot to path.
func SaveCheckpoint(path string, snap *Snapshot) error {
	if snap == nil {
		return failf(KindCheckpoint, -1, -1, "cannot save nil snapshot")
	}
	b, err := snap.Marshal()
	if err != nil {
		return err
	}
	if err := writeFileAtomic(path, b); err != nil {
		return err
	}
	return nil
}

// LoadCheckpoint reads and validates a checkpoint file.
func LoadCheckpoint(path string) (*Snapshot, error) {
	b, err := os.ReadFile(path)
	if err != nil {
		return nil, failf(KindCheckpoint, -1, -1, "cannot read checkpoint %s: %v", path, err)
	}
	snap, err := UnmarshalSnapshot(b)
	if err != nil {
		return nil, fmt.Errorf("%w (%s)", err, path)
	}
	return snap, nil
}
