package pde

// RunOptions controls one (possibly resumed) simulation run. All behavior is
// deterministic: identical Config + initial field + worker count produce
// identical trajectories regardless of how the steps are split across runs.
type RunOptions struct {
	Workers          int    // parallelism; falls back to Config.Workers when <=0
	Checkpoint       string // checkpoint file path ("" disables checkpointing)
	Resume           bool   // start by loading Checkpoint instead of t=0
	CheckpointEvery  int    // save every k steps (0 == save only at end)
	CheckInvariant   bool   // verify discrete heat balance after every step
	InvariantRegions int    // x-slabs for residual localization (default: workers)
	Logger           func(string)
}

// RunResult summarizes a completed run.
type RunResult struct {
	Steps         int
	Time          float64
	Final         *Field
	Snapshot      *Snapshot
	LastInvariant *InvariantReport
	MaxRelDelta   float64
}

// RunSimulation builds a solver for cfg, optionally resumes from a checkpoint,
// advances to the configured end time, optionally checks conservation at each
// step and saves checkpoints at deterministic intervals.
func RunSimulation(cfg Config, opt RunOptions) (*RunResult, error) {
	workers := opt.Workers
	if workers <= 0 {
		workers = cfg.Workers
	}
	solver, err := NewSolver(cfg)
	if err != nil {
		return nil, err
	}
	if opt.Resume {
		if opt.Checkpoint == "" {
			return nil, failf(KindCheckpoint, -1, -1, "resume requested without a checkpoint path")
		}
		snap, err := LoadCheckpoint(opt.Checkpoint)
		if err != nil {
			return nil, err
		}
		if err := solver.Restore(snap); err != nil {
			return nil, err
		}
		if opt.Logger != nil {
			opt.Logger("resumed from " + opt.Checkpoint + " at step " + itoa(snap.Step))
		}
	}

	regions := opt.InvariantRegions
	if regions <= 0 {
		if workers > 1 {
			regions = workers
		} else {
			regions = 1
		}
	}

	res := &RunResult{}
	plan := solver.PlannedSteps()
	checkInvariant := opt.CheckInvariant || cfg.ConserveEnforce
	for solver.StepIndex() < plan {
		prev := solver.cur.clone()
		prevTime := solver.Time()
		if err := solver.Step(workers); err != nil {
			return nil, err
		}
		if checkInvariant {
			rep, err := CheckConservation(cfg, prev, solver.cur, solver.StepIndex(), prevTime, cfg.Dt, regions)
			if err != nil {
				return res, err
			}
			res.LastInvariant = rep
			if rep.RelDelta > res.MaxRelDelta {
				res.MaxRelDelta = rep.RelDelta
			}
		}
		if opt.Checkpoint != "" && opt.CheckpointEvery > 0 && solver.StepIndex()%opt.CheckpointEvery == 0 {
			if err := SaveCheckpoint(opt.Checkpoint, solver.Snapshot()); err != nil {
				return nil, err
			}
		}
	}
	res.Steps = solver.StepIndex()
	res.Time = solver.Time()
	res.Final = solver.cur.clone()
	snap := solver.Snapshot()
	res.Snapshot = snap
	if opt.Checkpoint != "" {
		if err := SaveCheckpoint(opt.Checkpoint, snap); err != nil {
			return nil, err
		}
	}
	return res, nil
}
