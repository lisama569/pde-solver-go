package pde

import "math"

// sourceAt evaluates the configured source q(x,y,t).
//
//	q = 0                                       (SourceNone)
//	q = Scale                                   (SourceConstant)
//	q = Scale * exp(-((x-Cx)^2+(y-Cy)^2)/(2 Sigma^2))  (SourceGaussian)
//
// All built-in sources are time-independent; t is accepted so callers use a
// uniform signature and custom sources can be added later.
func sourceAt(cfg SourceConfig, x, y, t float64) float64 {
	_ = t
	switch cfg.Type {
	case SourceNone:
		return 0
	case SourceConstant:
		return cfg.Scale
	case SourceGaussian:
		dx := x - cfg.Cx
		dy := y - cfg.Cy
		return cfg.Scale * math.Exp(-(dx*dx+dy*dy)/(2.0*cfg.Sigma*cfg.Sigma))
	default:
		// Unreachable for validated configs; never silently inject 0.
		return math.NaN()
	}
}
