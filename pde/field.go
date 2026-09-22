package pde

// Field is a 2D grid with one ghost layer on every side, stored in flat
// row-major order. Interior point (i,j), 0<=i<Nx, 0<=j<Ny, is stored at
//
//	index = (j+1)*stride + (i+1),  stride = Nx+2.
type Field struct {
	Nx, Ny int
	stride int
	data   []float64
}

// newField allocates a zero field.
func newField(nx, ny int) *Field {
	return &Field{
		Nx:     nx,
		Ny:     ny,
		stride: nx + 2,
		data:   make([]float64, (nx+2)*(ny+2)),
	}
}

// idx returns the flat index of interior point (i,j).
func (f *Field) idx(i, j int) int { return (j+1)*f.stride + (i + 1) }

// clone returns an independent deep copy.
func (f *Field) clone() *Field {
	d := make([]float64, len(f.data))
	copy(d, f.data)
	return &Field{Nx: f.Nx, Ny: f.Ny, stride: f.stride, data: d}
}

// at returns the interior value u(i,j).
func (f *Field) at(i, j int) float64 { return f.data[f.idx(i, j)] }

// set writes an interior value.
func (f *Field) set(i, j int, v float64) { f.data[f.idx(i, j)] = v }

// ghostAt returns a ghost-slot value adjacent to interior point (i,j),
// where gx,gy in {-1,0,1} indicate the slot offset. Only edge-adjacent
// slots (one offset non-zero) are ever read by the 5-point stencil.
func (f *Field) ghostAt(i, j, gx, gy int) float64 {
	return f.data[(j+1+gy)*f.stride+(i+1+gx)]
}

// setGhost writes a ghost slot.
func (f *Field) setGhost(i, j, gx, gy int, v float64) {
	f.data[(j+1+gy)*f.stride+(i+1+gx)] = v
}

// sumInterior sums all interior values in fixed row-major order. The fixed
// order makes the result independent of worker count and execution order.
func (f *Field) sumInterior() float64 {
	var s float64
	for j := 0; j < f.Ny; j++ {
		base := (j + 1) * f.stride
		for i := 0; i < f.Nx; i++ {
			s += f.data[base+i+1]
		}
	}
	return s
}

// interiorCopy returns interior values in row-major order.
func (f *Field) interiorCopy() []float64 {
	out := make([]float64, f.Nx*f.Ny)
	for j := 0; j < f.Ny; j++ {
		copy(out[j*f.Nx:(j+1)*f.Nx], f.data[(j+1)*f.stride+1:(j+1)*f.stride+1+f.Nx])
	}
	return out
}

// setInterior loads interior values given in row-major order.
func (f *Field) setInterior(vals []float64) error {
	if len(vals) != f.Nx*f.Ny {
		return failf(KindCheckpoint, -1, -1,
			"checkpoint field length %d does not match grid %dx%d", len(vals), f.Nx, f.Ny)
	}
	for j := 0; j < f.Ny; j++ {
		copy(f.data[(j+1)*f.stride+1:(j+1)*f.stride+1+f.Nx], vals[j*f.Nx:(j+1)*f.Nx])
	}
	return nil
}

// SumInterior is the exported fixed-order sum of interior values.
func (f *Field) SumInterior() float64 { return f.sumInterior() }

// Dims returns the interior grid dimensions.
func (f *Field) Dims() (nx, ny int) { return f.Nx, f.Ny }

// At is the exported interior accessor.
func (f *Field) At(i, j int) float64 { return f.at(i, j) }
