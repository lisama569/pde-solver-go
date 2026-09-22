package pde

import (
	"fmt"
	"math"
)

// BoundaryType enumerates supported boundary conditions.
type BoundaryType int

const (
	BoundaryDirichlet BoundaryType = iota // prescribed value
	BoundaryNeumann                       // prescribed outward normal derivative
	BoundaryPeriodic                      // periodic wrap
)

// String renders a BoundaryType.
func (t BoundaryType) String() string {
	switch t {
	case BoundaryDirichlet:
		return "dirichlet"
	case BoundaryNeumann:
		return "neumann"
	case BoundaryPeriodic:
		return "periodic"
	default:
		return "unknown"
	}
}

// BoundaryCondition describes one side of the rectangular domain.
//
// For Neumann, Value is the outward normal derivative du/dn at that side.
// Value is ignored for periodic sides.
type BoundaryCondition struct {
	Type  BoundaryType
	Value float64
}

// Boundaries holds conditions for the low/high side of each axis.
type Boundaries struct {
	X0 BoundaryCondition
	X1 BoundaryCondition
	Y0 BoundaryCondition
	Y1 BoundaryCondition
}

// SourceType enumerates built-in source terms.
type SourceType int

const (
	SourceNone SourceType = iota
	SourceConstant
	SourceGaussian
)

// SourceConfig configures a source term q(x,y,t).
type SourceConfig struct {
	Type   SourceType
	Scale  float64
	Cx, Cy float64 // Gaussian center
	Sigma  float64 // Gaussian standard deviation
}

// Config fully describes a simulation.
type Config struct {
	Nx, Ny          int     // interior grid points per axis
	Dx, Dy          float64 // grid spacing
	Alpha           float64 // thermal diffusivity (>= 0)
	TotalTime       float64 // physical end time
	Dt              float64 // time step
	Boundaries      Boundaries
	Source          SourceConfig
	Workers         int     // subdomain parallelism (0/1 == serial)
	InitValue       float64 // uniform initial condition
	MaxValue        float64 // finite sanity bound; |u| above it => blowup (<=0 disables)
	ConserveTol     float64 // relative conservation tolerance (<=0 uses default)
	ConserveEnforce bool    // treat conservation violation as hard error
	Seed            int64   // deterministic perturbation seed
	NoiseAmp        float64 // deterministic perturbation amplitude (0 disables)
}

// DefaultConserveTol is used when Config.ConserveTol <= 0.
const DefaultConserveTol = 1e-10

// DefaultMaxValue is used when Config.MaxValue <= 0.
const DefaultMaxValue = 1e12

// finite reports whether v is a usable finite number.
func finite(v float64) bool { return !math.IsNaN(v) && !math.IsInf(v, 0) }

// Steps returns the number of time steps TotalTime/Dt, which must divide the
// total time evenly (relative tolerance 1e-9).
func (c Config) Steps() (int, error) {
	if c.Dt <= 0 || c.TotalTime <= 0 {
		return 0, failf(KindConfig, -1, -1, "cannot compute steps: total_time=%v dt=%v must both be positive", c.TotalTime, c.Dt)
	}
	exact := c.TotalTime / c.Dt
	n := int(math.Round(exact))
	if n < 1 {
		return 0, failf(KindConfig, -1, -1, "total_time %v is shorter than one dt %v", c.TotalTime, c.Dt)
	}
	if math.Abs(float64(n)*c.Dt-c.TotalTime) > 1e-9*c.TotalTime {
		return 0, failf(KindConfig, -1, -1,
			"total_time %v is not an integer multiple of dt %v (got %.12f steps)",
			c.TotalTime, c.Dt, exact)
	}
	return n, nil
}

// validateSide checks one boundary specification.
func validateSide(name string, bc BoundaryCondition) error {
	switch bc.Type {
	case BoundaryDirichlet, BoundaryNeumann:
		if !finite(bc.Value) {
			return failf(KindConfig, -1, -1, "boundary %s: %s value must be finite, got %v", name, bc.Type, bc.Value)
		}
	case BoundaryPeriodic:
		// value ignored
	default:
		return failf(KindConfig, -1, -1, "boundary %s: unknown type %d", name, bc.Type)
	}
	return nil
}

// Validate checks legality of the configuration and returns *Failure values
// of KindConfig. It does not check stability; see CheckStability.
func (c Config) Validate() error {
	if err := c.ValidateGeometryOnly(); err != nil {
		return err
	}
	if !(c.Dt > 0) || !finite(c.Dt) {
		return failf(KindConfig, -1, -1, "dt must be positive and finite, got %v", c.Dt)
	}
	if _, err := c.Steps(); err != nil {
		return err
	}
	return nil
}

// ValidateGeometryOnly validates every field except dt-dependent quantities
// (dt itself and the total-time step divisibility), so callers can compute an
// automatic stable dt from the geometry before the full validation.
func (c Config) ValidateGeometryOnly() error {
	if c.Nx < 2 || c.Ny < 2 {
		return failf(KindConfig, -1, -1, "grid size must be >= 2 in each dimension, got Nx=%d Ny=%d", c.Nx, c.Ny)
	}
	if !(c.Dx > 0) || !finite(c.Dx) {
		return failf(KindConfig, -1, -1, "dx must be a positive finite number, got %v", c.Dx)
	}
	if !(c.Dy > 0) || !finite(c.Dy) {
		return failf(KindConfig, -1, -1, "dy must be a positive finite number, got %v", c.Dy)
	}
	if c.Alpha < 0 || !finite(c.Alpha) {
		return failf(KindConfig, -1, -1, "alpha (diffusivity) must be >= 0 and finite, got %v", c.Alpha)
	}
	if !(c.TotalTime > 0) || !finite(c.TotalTime) {
		return failf(KindConfig, -1, -1, "total_time must be positive and finite, got %v", c.TotalTime)
	}
	if !finite(c.InitValue) {
		return failf(KindConfig, -1, -1, "init_value must be finite, got %v", c.InitValue)
	}
	if c.MaxValue < 0 || !finite(c.MaxValue) {
		return failf(KindConfig, -1, -1, "max_value must be >= 0 and finite, got %v", c.MaxValue)
	}
	if c.ConserveTol < 0 || !finite(c.ConserveTol) {
		return failf(KindConfig, -1, -1, "conserve_tol must be >= 0 and finite, got %v", c.ConserveTol)
	}
	if c.NoiseAmp < 0 || !finite(c.NoiseAmp) {
		return failf(KindConfig, -1, -1, "noise_amp must be >= 0 and finite, got %v", c.NoiseAmp)
	}
	if c.Workers < 0 {
		return failf(KindConfig, -1, -1, "workers must be >= 0, got %d", c.Workers)
	}
	for name, bc := range map[string]BoundaryCondition{
		"x0": c.Boundaries.X0, "x1": c.Boundaries.X1,
		"y0": c.Boundaries.Y0, "y1": c.Boundaries.Y1,
	} {
		if err := validateSide(name, bc); err != nil {
			return err
		}
	}
	// Periodicity must be applied in matched pairs.
	if (c.Boundaries.X0.Type == BoundaryPeriodic) != (c.Boundaries.X1.Type == BoundaryPeriodic) {
		return failf(KindConfig, -1, -1, "periodic boundary must be set on both x sides (x0=%s x1=%s)",
			c.Boundaries.X0.Type, c.Boundaries.X1.Type)
	}
	if (c.Boundaries.Y0.Type == BoundaryPeriodic) != (c.Boundaries.Y1.Type == BoundaryPeriodic) {
		return failf(KindConfig, -1, -1, "periodic boundary must be set on both y sides (y0=%s y1=%s)",
			c.Boundaries.Y0.Type, c.Boundaries.Y1.Type)
	}
	switch c.Source.Type {
	case SourceNone:
	case SourceConstant:
		if !finite(c.Source.Scale) {
			return failf(KindConfig, -1, -1, "constant source scale must be finite, got %v", c.Source.Scale)
		}
	case SourceGaussian:
		if !finite(c.Source.Scale) || !finite(c.Source.Cx) || !finite(c.Source.Cy) {
			return failf(KindConfig, -1, -1, "gaussian source scale/center must be finite")
		}
		if !(c.Source.Sigma > 0) || !finite(c.Source.Sigma) {
			return failf(KindConfig, -1, -1, "gaussian source sigma must be positive finite, got %v", c.Source.Sigma)
		}
	default:
		return failf(KindConfig, -1, -1, "unknown source type %d", c.Source.Type)
	}
	return nil
}

// AutoDt picks a dt that is both stable and an exact divisor of TotalTime:
//
//	nsteps = ceil(T / (0.9*dt_limit)),  dt = T / nsteps.
//
// The 0.9 safety margin leaves headroom below the stability limit.
func (c Config) AutoDt() (float64, int, error) {
	if err := c.ValidateGeometryOnly(); err != nil {
		return 0, 0, err
	}
	margin := 0.9 * c.stabilityLimit()
	n := int(math.Ceil(c.TotalTime/margin - 1e-12))
	if n < 1 {
		n = 1
	}
	return c.TotalTime / float64(n), n, nil
}

// diffusionRatios returns rx = alpha*dt/dx^2 and ry = alpha*dt/dy^2.
func (c Config) diffusionRatios() (rx, ry float64) {
	return c.Alpha * c.Dt / (c.Dx * c.Dx), c.Alpha * c.Dt / (c.Dy * c.Dy)
}

// stabilityLimit returns the FTCS 2D CFL upper bound on dt:
//
//	rx+ry <= 1/2  =>  dt <= 1 / (4*alpha*(1/dx^2 + 1/dy^2)).
//
// It returns +Inf when alpha == 0 (no diffusion means no diffusive
// stability restriction).
func (c Config) stabilityLimit() float64 {
	if c.Alpha == 0 {
		return math.Inf(1)
	}
	return 1.0 / (4.0 * c.Alpha * (1.0/(c.Dx*c.Dx) + 1.0/(c.Dy*c.Dy)))
}

// CheckStability rejects step sizes or coefficients outside the explicit
// FTCS stability range: rx + ry <= 1/2, where rx = alpha*dt/dx^2 and
// ry = alpha*dt/dy^2. Violations produce KindStability errors.
func (c Config) CheckStability() error {
	rx, ry := c.diffusionRatios()
	r := rx + ry
	if math.IsNaN(r) || math.IsInf(r, 0) || r > 0.5 {
		return failf(KindStability, -1, -1,
			"unstable parameters: rx=alpha*dt/dx^2=%.12g ry=alpha*dt/dy^2=%.12g sum=%.12g > 0.5; dt must be <= %.12g",
			rx, ry, r, c.stabilityLimit())
	}
	return nil
}

// String renders a compact deterministic parameter summary used in logs.
func (c Config) String() string {
	return fmt.Sprintf(
		"Nx=%d Ny=%d dx=%g dy=%g alpha=%g dt=%g total_time=%g workers=%d init=%g source=%d scale=%g bcs=[x0:%s(%g) x1:%s(%g) y0:%s(%g) y1:%s(%g)]",
		c.Nx, c.Ny, c.Dx, c.Dy, c.Alpha, c.Dt, c.TotalTime, c.Workers, c.InitValue,
		c.Source.Type, c.Source.Scale,
		c.Boundaries.X0.Type, c.Boundaries.X0.Value, c.Boundaries.X1.Type, c.Boundaries.X1.Value,
		c.Boundaries.Y0.Type, c.Boundaries.Y0.Value, c.Boundaries.Y1.Type, c.Boundaries.Y1.Value)
}

// DiffusionRatios is the exported accessor for rx, ry.
func (c Config) DiffusionRatios() (rx, ry float64) { return c.diffusionRatios() }

// StabilityLimit is the exported accessor for the FTCS dt upper bound.
func (c Config) StabilityLimit() float64 { return c.stabilityLimit() }
