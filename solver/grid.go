package solver

import "math"

// Grid 是行优先存储的二维标量场，行索引 j∈[0,NY)，列索引 i∈[0,NX)。
type Grid struct {
	NX, NY int
	Data   []float64
}

// NewGrid 创建零初始化网格。
func NewGrid(nx, ny int) *Grid {
	return &Grid{NX: nx, NY: ny, Data: make([]float64, nx*ny)}
}

// At 返回 (i,j) 处的值。
func (g *Grid) At(i, j int) float64 { return g.Data[j*g.NX+i] }

// Set 设置 (i,j) 处的值。
func (g *Grid) Set(i, j int, v float64) { g.Data[j*g.NX+i] = v }

// Clone 深拷贝网格。
func (g *Grid) Clone() *Grid {
	c := &Grid{NX: g.NX, NY: g.NY, Data: make([]float64, len(g.Data))}
	copy(c.Data, g.Data)
	return c
}

// MaxAbsDiff 返回两网格最大绝对差；尺寸不一致返回 +Inf。
func (g *Grid) MaxAbsDiff(o *Grid) float64 {
	if g.NX != o.NX || g.NY != o.NY {
		return math.Inf(1)
	}
	max := 0.0
	for k := range g.Data {
		if d := math.Abs(g.Data[k] - o.Data[k]); d > max {
			max = d
		}
	}
	return max
}
