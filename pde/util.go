package pde

import "strconv"

func itoa(n int) string { return strconv.Itoa(n) }

// NewSolverOrValidate runs construction purely to surface validation errors.
func NewSolverOrValidate(cfg Config) error {
	_, err := NewSolver(cfg)
	return err
}
