package pde

// All failures are reported as *Failure with a FailureKind so callers can
// distinguish illegal configuration, a-priori stability rejection, runtime
// numerical blow-up, checkpoint I/O problems and invariant violations.
// Package-level overview lives in doc.go.

import (
	"errors"
	"fmt"
)

// FailureKind classifies a solver failure.
type FailureKind int

const (
	KindUnknown    FailureKind = iota
	KindConfig                 // invalid/illegal configuration
	KindStability              // configuration rejected a-priori by the stability condition
	KindBlowup                 // runtime numerical instability detected (NaN/Inf/out-of-range)
	KindCheckpoint             // checkpoint write/read failure
	KindInvariant              // conservation/invariant violation beyond tolerance
)

// String renders a FailureKind.
func (k FailureKind) String() string {
	switch k {
	case KindConfig:
		return "config"
	case KindStability:
		return "stability"
	case KindBlowup:
		return "blowup"
	case KindCheckpoint:
		return "checkpoint"
	case KindInvariant:
		return "invariant"
	default:
		return "unknown"
	}
}

// Failure is the structured error type returned by the package.
type Failure struct {
	Kind   FailureKind
	Step   int    // time step at which the failure was detected (-1 if N/A)
	Domain int    // subdomain at which the failure was detected (-1 if N/A)
	Msg    string // human-readable explanation
}

func (e *Failure) Error() string {
	if e.Domain >= 0 {
		return fmt.Sprintf("pde %s failure at step %d subdomain %d: %s",
			e.Kind, e.Step, e.Domain, e.Msg)
	}
	if e.Step >= 0 {
		return fmt.Sprintf("pde %s failure at step %d: %s", e.Kind, e.Step, e.Msg)
	}
	return fmt.Sprintf("pde %s failure: %s", e.Kind, e.Msg)
}

// failf constructs *Failure values.
func failf(kind FailureKind, step, domain int, format string, args ...any) *Failure {
	return &Failure{Kind: kind, Step: step, Domain: domain, Msg: fmt.Sprintf(format, args...)}
}

// AsFailure extracts a *Failure from err, or nil if err is not a *Failure.
func AsFailure(err error) *Failure {
	var f *Failure
	if errors.As(err, &f) {
		return f
	}
	return nil
}
