// Package pde is a deterministic, checkpointable 2D heat-equation solver
// with subdomain-parallel execution.
//
// It solves
//
//	∂u/∂t = α(∂²u/∂x² + ∂²u/∂y²) + q(x,y,t)
//
// on a uniform rectangular grid with an explicit FTCS five-point update,
// one ghost layer around the field, and double buffering:
//
//	u^(n+1) = u^n + rx Δx_lap u^n + ry Δy_lap u^n + dt q(t^n),
//	rx = α dt/dx²,  ry = α dt/dy².
//
// Stability: rx+ry <= 1/2 is enforced before the first step
// (Config.CheckStability, KindStability); illegal inputs fail earlier as
// KindConfig. Every produced value is checked at runtime and NaN/Inf or
// values beyond Config.MaxValue fail the step as KindBlowup with step,
// subdomain and grid location; the double buffer means a failed step is
// never committed.
//
// Parallelism: the interior columns are partitioned into deterministic
// contiguous slabs; workers read only the shared old buffer and write
// disjoint slabs of the new buffer, so serial and parallel results are
// bit-for-bit identical and free of races (verified with -race).
//
// Persistence: Snapshot values are serialized to deterministic JSON
// (SaveCheckpoint/LoadCheckpoint) via a temp-file-and-rename write, and
// resume reproduces the uninterrupted trajectory exactly.
//
// Conservation: CheckConservation verifies the per-step discrete heat
// balance H^n+1-H^n = source injection + boundary influx using the same
// ghost values as the update, and localizes residuals to x-slabs. Failures
// beyond tolerance are KindInvariant.
//
// See the repository README for the numerical details, boundary/parallel
// constraints and the verification matrix.
package pde
