package pde

import (
	"math"
	"math/rand"
)

// SchemaVersion is the checkpoint format version produced by this package.
const SchemaVersion = 1

// Solver advances a 2D heat equation in time with double-buffered fields:
//
//	u^{n+1}_{i,j} = u^n_{i,j}
//	  + rx (u^n_{i-1,j}+u^n_{i+1,j}-2u^n_{i,j})
//	  + ry (u^n_{i,j-1}+u^n_{i,j+1}-2u^n_{i,j})
//	  + dt * q(x_i, y_j, t^n),
//
// with rx = alpha*dt/dx^2, ry = alpha*dt/dy^2.
type Solver struct {
	cfg      Config
	steps    int
	cur, nxt *Field
	step     int
	time     float64
	rx, ry   float64
	maxValue float64
}

// NewSolver validates the configuration, rejects a-priori unstable
// parameters, allocates fields and builds the deterministic initial state.
func NewSolver(cfg Config) (*Solver, error) {
	if err := cfg.Validate(); err != nil {
		return nil, err
	}
	if err := cfg.CheckStability(); err != nil {
		return nil, err
	}
	steps, err := cfg.Steps()
	if err != nil {
		return nil, err
	}
	rx, ry := cfg.diffusionRatios()
	s := &Solver{
		cfg:      cfg,
		steps:    steps,
		cur:      newField(cfg.Nx, cfg.Ny),
		nxt:      newField(cfg.Nx, cfg.Ny),
		rx:       rx,
		ry:       ry,
		maxValue: cfg.MaxValue,
	}
	if s.maxValue <= 0 {
		s.maxValue = DefaultMaxValue
	}
	s.initState()
	return s, nil
}

// initState writes u = InitValue plus an optional deterministic perturbation
// generated from a seeded local rand source. Serial and parallel runs use the
// exact same (i,j)->value mapping, so initial states are identical.
func (s *Solver) initState() {
	var rng *rand.Rand
	if s.cfg.NoiseAmp != 0 && s.cfg.Seed != 0 {
		rng = rand.New(rand.NewSource(s.cfg.Seed))
	}
	for j := 0; j < s.cfg.Ny; j++ {
		for i := 0; i < s.cfg.Nx; i++ {
			v := s.cfg.InitValue
			if rng != nil {
				v += (rng.Float64()*2 - 1) * s.cfg.NoiseAmp
			}
			s.cur.set(i, j, v)
		}
	}
}

// Step advances exactly one time step with the requested parallelism
// (workers<=1 selects the serial baseline).
func (s *Solver) Step(workers int) error {
	tnew := s.time + s.cfg.Dt
	if err := updateParallel(s.cfg, s.cur, s.nxt, s.rx, s.ry, tnew, s.step+1, workers); err != nil {
		return err
	}
	s.cur, s.nxt = s.nxt, s.cur
	s.step++
	s.time = tnew
	return nil
}

// Run advances n steps (clamped to the remaining planned steps) using the
// solver's configured worker count.
func (s *Solver) Run(n int) error {
	remaining := s.steps - s.step
	if n > remaining {
		n = remaining
	}
	for k := 0; k < n; k++ {
		if err := s.Step(s.cfg.Workers); err != nil {
			return err
		}
	}
	return nil
}

// StepIndex returns the number of completed steps.
func (s *Solver) StepIndex() int { return s.step }

// Time returns the current physical time.
func (s *Solver) Time() float64 { return s.time }

// PlannedSteps returns the total step count implied by TotalTime/Dt.
func (s *Solver) PlannedSteps() int { return s.steps }

// Field returns a snapshot copy of the current field (with ghost layers
// already consistent with the interior state).
func (s *Solver) Field() *Field {
	out := s.cur.clone()
	applyBoundaries(out, s.cfg.Boundaries, s.cfg.Dx, s.cfg.Dy)
	return out
}

// Snapshot captures full resumable state.
func (s *Solver) Snapshot() *Snapshot {
	return &Snapshot{
		SchemaVersion: SchemaVersion,
		Config:        s.cfg,
		Step:          s.step,
		Time:          s.time,
		Data:          s.cur.interiorCopy(),
	}
}

// Restore replaces solver state from a snapshot, which must match the
// solver configuration and have a consistent step/time.
func (s *Solver) Restore(snap *Snapshot) error {
	if err := snap.checkCompatible(s.cfg); err != nil {
		return err
	}
	if err := s.cur.setInterior(snap.Data); err != nil {
		return err
	}
	s.step = snap.Step
	s.time = snap.Time
	return nil
}

// updateSubdomain computes interior points i in [i0,i1), j in [0,Ny) from
// src into dst. It reads only src and writes only its assigned slab of dst,
// so updates over disjoint slabs are race-free.
//
// Every produced value is checked for NaN/Inf and against maxValue; the first
// offending point fails the step with a KindBlowup error naming the step,
// subdomain and grid location. A step is never partially "accepted": dst is a
// scratch buffer and the solver only swaps buffers on full success.
func updateSubdomain(cfg Config, src, dst *Field, rx, ry, tnew float64, step, i0, i1 int) error {
	dt := cfg.Dt
	maxV := cfg.MaxValue
	if maxV <= 0 {
		maxV = DefaultMaxValue
	}
	tn := tnew - dt // q is evaluated at the old time level (forward Euler)
	for j := 0; j < src.Ny; j++ {
		y := cfg.Dy * (float64(j) + 0.5)
		for i := i0; i < i1; i++ {
			uc := src.at(i, j)
			ul := src.ghostAt(i, j, -1, 0)
			ur := src.ghostAt(i, j, 1, 0)
			ud := src.ghostAt(i, j, 0, -1)
			uu := src.ghostAt(i, j, 0, 1)
			x := cfg.Dx * (float64(i) + 0.5)
			v := uc + rx*(ul+ur-2*uc) + ry*(ud+uu-2*uc) + dt*sourceAt(cfg.Source, x, y, tn)
			if math.IsNaN(v) || math.IsInf(v, 0) {
				return failf(KindBlowup, step, -1,
					"non-finite value %v produced at (i=%d,j=%d) t=%.12g", v, i, j, tnew)
			}
			if math.Abs(v) > maxV {
				return failf(KindBlowup, step, -1,
					"|u|=%g exceeds sanity bound %g at (i=%d,j=%d) t=%.12g (numerical blow-up)",
					v, maxV, i, j, tnew)
			}
			dst.set(i, j, v)
		}
	}
	return nil
}
