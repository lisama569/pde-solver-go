package pde

import "math"

// RegionDelta is the discrete heat-balance residual of one x-column slab.
type RegionDelta struct {
	I0, I1 int
	Delta  float64 // H_new-H_old - source - net_influx for the slab
}

// InvariantReport records a discrete heat-balance check spanning one step.
//
// The discrete scheme obeys, per cell per step,
//
//	h_new-h_old = dt*q + dt*alpha/ds^2 * (neighbor sums),
//
// and neighbor exchange between adjacent cells cancels when summed over a
// slab; only physical-boundary ghost exchanges and slab-interface exchanges
// remain. Interface exchanges cancel between neighboring slabs, so the
// global balance contains only boundary flux plus source injection.
type InvariantReport struct {
	Step         int
	Time         float64
	Heat         float64       // dx*dy * sum of interior u at new time
	PrevHeat     float64       // dx*dy * sum of interior u at old time
	AbsDelta     float64       // |Heat - PrevHeat - SourceHeat - BoundaryFlux|
	RelDelta     float64       // AbsDelta / max(1, |expected new heat|)
	BoundaryFlux float64       // net heat entering through boundaries during the step
	SourceHeat   float64       // heat injected by the source during the step
	RegionDeltas []RegionDelta // per-slab residuals (interface fluxes included)
}

// conservationTolerance returns the effective relative tolerance.
func conservationTolerance(cfg Config) float64 {
	if cfg.ConserveTol > 0 {
		return cfg.ConserveTol
	}
	return DefaultConserveTol
}

// boundaryInfluxX0 sums heat entering across the x=0 face of column 0 during
// one step, using the SAME ghost/interior values the update used.
// Unit-y length per row is dy.
func boundaryInfluxX0(f *Field, cfg Config) float64 {
	var s float64
	for j := 0; j < f.Ny; j++ {
		// incoming = +alpha*dt/dx * (u_ghost - u(0))
		s += cfg.Alpha * cfg.Dt / cfg.Dx * (f.ghostAt(0, j, -1, 0) - f.at(0, j))
	}
	return s * cfg.Dy
}

// boundaryInfluxX1 sums heat entering across the x=Lx face of the last column.
func boundaryInfluxX1(f *Field, cfg Config) float64 {
	var s float64
	for j := 0; j < f.Ny; j++ {
		// incoming = +alpha*dt/dx * (u_ghost - u(Nx-1))
		s += cfg.Alpha * cfg.Dt / cfg.Dx * (f.ghostAt(f.Nx-1, j, 1, 0) - f.at(f.Nx-1, j))
	}
	return s * cfg.Dy
}

// interfaceInflux returns heat entering slab [a,b) through its left interface
// between columns a-1 and a (positive when it warms the slab), and through its
// right interface between b-1 and b. For an interior interface the matching
// term in the neighboring slab has the opposite sign, so terms cancel.
func interfaceInflux(f *Field, cfg Config, a, b int) (left, right float64) {
	for j := 0; j < f.Ny; j++ {
		if a > 0 {
			left += cfg.Alpha * cfg.Dt / cfg.Dx * (f.at(a-1, j) - f.at(a, j))
		}
		if b < f.Nx {
			right += cfg.Alpha * cfg.Dt / cfg.Dx * (f.at(b, j) - f.at(b-1, j))
		}
	}
	return left * cfg.Dy, right * cfg.Dy
}

// CheckConservation verifies the discrete heat balance from prev at time
// prevTime (= step*dt - dt) to cur at step step. nRegions controls how many
// x-slabs the residual is localized into. cur is read only; prev is cloned
// internally before ghost layers are refreshed.
func CheckConservation(cfg Config, prev, cur *Field, step int, prevTime, dt float64, nRegions int) (*InvariantReport, error) {
	if prev.Nx != cur.Nx || prev.Ny != cur.Ny {
		return nil, failf(KindInvariant, step, -1, "field size mismatch in conservation check")
	}
	old := prev.clone()
	applyBoundaries(old, cfg.Boundaries, cfg.Dx, cfg.Dy)

	slabs := partitionRows(cfg.Nx, nRegions)
	regions := make([]RegionDelta, 0, len(slabs))
	cellArea := cfg.Dx * cfg.Dy

	var hNew, hOld, fluxTotal, sourceTotal float64
	for _, sl := range slabs {
		a, b := sl[0], sl[1]
		var hsNew, hsOld, src float64
		for j := 0; j < cfg.Ny; j++ {
			y := cfg.Dy * (float64(j) + 0.5)
			for i := a; i < b; i++ {
				x := cfg.Dx * (float64(i) + 0.5)
				hsNew += cur.at(i, j)
				hsOld += old.at(i, j)
				src += dt * sourceAt(cfg.Source, x, y, prevTime)
			}
		}
		hsNew *= cellArea
		hsOld *= cellArea
		src *= cellArea

		var flux float64
		// Physical x-boundary contributions only for slabs touching sides.
		if a == 0 {
			flux += boundaryInfluxX0(old, cfg)
		}
		if b == cfg.Nx {
			flux += boundaryInfluxX1(old, cfg)
		}
		// Internal x-interface exchanges.
		left, right := interfaceInflux(old, cfg, a, b)
		flux += left + right
		// y-boundary contributions for every column in the slab: incoming
		// across y0/y1 faces, unit-x length dx per column.
		for i := a; i < b; i++ {
			flux += cfg.Alpha * cfg.Dt / cfg.Dy *
				(old.ghostAt(i, 0, 0, -1) - old.at(i, 0) +
					old.ghostAt(i, cfg.Ny-1, 0, 1) - old.at(i, cfg.Ny-1)) * cfg.Dx
		}

		delta := hsNew - hsOld - src - flux
		regions = append(regions, RegionDelta{I0: a, I1: b, Delta: delta})
		hNew += hsNew
		hOld += hsOld
		sourceTotal += src
		fluxTotal += flux
	}

	absDelta := math.Abs(hNew - hOld - sourceTotal - fluxTotal)
	denom := math.Max(1.0, math.Abs(hOld+sourceTotal+fluxTotal))
	rep := &InvariantReport{
		Step:         step,
		Time:         prevTime + dt,
		Heat:         hNew,
		PrevHeat:     hOld,
		AbsDelta:     absDelta,
		RelDelta:     absDelta / denom,
		BoundaryFlux: fluxTotal,
		SourceHeat:   sourceTotal,
		RegionDeltas: regions,
	}
	tol := conservationTolerance(cfg)
	if rep.RelDelta > tol {
		worst := regions[0]
		for _, r := range regions[1:] {
			if math.Abs(r.Delta) > math.Abs(worst.Delta) {
				worst = r
			}
		}
		return rep, failf(KindInvariant, step, -1,
			"heat balance residual relative %.6g exceeds tolerance %.6g (abs=%.6g); worst slab i=[%d,%d) delta=%.6g",
			rep.RelDelta, tol, rep.AbsDelta, worst.I0, worst.I1, worst.Delta)
	}
	return rep, nil
}
