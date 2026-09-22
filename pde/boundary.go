package pde

// applyBoundaries fills every edge-adjacent ghost slot of f from its interior
// state and the prescribed boundary conditions. The physical boundary is
// placed half a cell outside the first/last interior point:
//
//	Dirichlet v: boundary sits at the midpoint between the ghost slot and the
//	             first interior cell center, so
//	             v = (u_ghost + u_interior)/2, u_ghost = 2v - u_interior.
//	Neumann  g: g = du/dn is the outward-normal derivative. Along the outward
//	             direction g = (u_ghost - u_interior)/dx (or /dy), hence
//	             u_ghost = u_interior + g*dx  (resp. + g*dy) on EVERY side.
//	Periodic:    the ghost wraps to the opposite interior column/row.
//
// Corner ghost slots are never read by the 5-point stencil and stay zero.
func applyBoundaries(f *Field, b Boundaries, dx, dy float64) {
	nx, ny := f.Nx, f.Ny

	for j := 0; j < ny; j++ {
		u0 := f.at(0, j)
		uN := f.at(nx-1, j)
		switch b.X0.Type {
		case BoundaryDirichlet:
			f.setGhost(0, j, -1, 0, 2.0*b.X0.Value-u0)
		case BoundaryNeumann:
			f.setGhost(0, j, -1, 0, u0+b.X0.Value*dx)
		case BoundaryPeriodic:
			f.setGhost(0, j, -1, 0, uN)
		}
		switch b.X1.Type {
		case BoundaryDirichlet:
			f.setGhost(nx-1, j, 1, 0, 2.0*b.X1.Value-uN)
		case BoundaryNeumann:
			f.setGhost(nx-1, j, 1, 0, uN+b.X1.Value*dx)
		case BoundaryPeriodic:
			f.setGhost(nx-1, j, 1, 0, u0)
		}
	}
	for i := 0; i < nx; i++ {
		u0 := f.at(i, 0)
		uN := f.at(i, ny-1)
		switch b.Y0.Type {
		case BoundaryDirichlet:
			f.setGhost(i, 0, 0, -1, 2.0*b.Y0.Value-u0)
		case BoundaryNeumann:
			f.setGhost(i, 0, 0, -1, u0+b.Y0.Value*dy)
		case BoundaryPeriodic:
			f.setGhost(i, 0, 0, -1, f.at(i, ny-1))
		}
		switch b.Y1.Type {
		case BoundaryDirichlet:
			f.setGhost(i, ny-1, 0, 1, 2.0*b.Y1.Value-uN)
		case BoundaryNeumann:
			f.setGhost(i, ny-1, 0, 1, uN+b.Y1.Value*dy)
		case BoundaryPeriodic:
			f.setGhost(i, ny-1, 0, 1, f.at(i, 0))
		}
	}
}
