package pde

import "sync"

// partitionRows splits nx interior columns into n contiguous near-equal
// subdomains [start,end). The partition is deterministic and depends only
// on (nx, n), so every run with the same worker count decomposes identically.
// n<=1 yields a single full-width subdomain.
func partitionRows(nx, n int) [][2]int {
	nParts := n
	if nParts <= 1 || nParts > nx {
		nParts = 1
	}
	if nx == 0 {
		return nil
	}
	parts := make([][2]int, nParts)
	base := nx / nParts
	rem := nx % nParts
	start := 0
	for w := 0; w < nParts; w++ {
		width := base
		if w < rem {
			width++
		}
		parts[w] = [2]int{start, start + width}
		start += width
	}
	if start != nx { // defensive; should never happen
		return [][2]int{{0, nx}}
	}
	return parts
}

// runSubdomain executes fn on one worker and converts unexpected panics into
// KindBlowup failures tagged with the subdomain index, so a worker never
// crashes the process and failures stay attributable to a region.
func runSubdomain(fn func() error) (err error) {
	defer func() {
		if r := recover(); r != nil {
			err = failf(KindBlowup, -1, -1, "panic during subdomain update: %v", r)
		}
	}()
	return fn()
}

// updateParallel fills ghost layers of src from its interior, then updates
// all interior points of dst using the given worker count.
//
// Because src and dst are distinct buffers and each worker writes a disjoint
// i-column slab of dst, workers need no halo exchange and cannot race: all
// cross-subdomain values are read from the immutable src, including its
// freshly filled ghost layers (which encode periodic/physical boundaries).
// Results are therefore bit-for-bit independent of worker count and
// scheduling order.
func updateParallel(cfg Config, src, dst *Field, rx, ry, tnew float64, step, workers int) error {
	applyBoundaries(src, cfg.Boundaries, cfg.Dx, cfg.Dy)
	parts := partitionRows(src.Nx, workers)
	if len(parts) == 1 {
		return updateSubdomain(cfg, src, dst, rx, ry, tnew, step, parts[0][0], parts[0][1])
	}

	var wg sync.WaitGroup
	errs := make([]error, len(parts))
	for w, p := range parts {
		wg.Add(1)
		go func(w, i0, i1 int) {
			defer wg.Done()
			errs[w] = runSubdomain(func() error {
				return updateSubdomain(cfg, src, dst, rx, ry, tnew, step, i0, i1)
			})
		}(w, p[0], p[1])
	}
	wg.Wait()
	for w, err := range errs { // fixed order: lowest subdomain reported first
		if err != nil {
			if f := AsFailure(err); f != nil && f.Domain < 0 {
				f.Domain = w
			}
			return err
		}
	}
	return nil
}
