package solver

import (
	"encoding/binary"
	"fmt"
	"math"
	"os"
)

// checkpointMagic 用于识别文件格式，版本号含在低字节。
const checkpointMagic uint64 = 0x504445434B505431 // "PDECKPT1"

// Checkpoint 是可落盘的运行状态。
type Checkpoint struct {
	Step int
	NX   int
	NY   int
	Data []float64
}

// SaveCheckpoint 将状态以二进制格式写入文件（先写临时文件再原子改名）。
func SaveCheckpoint(path string, step int, g *Grid) error {
	buf := make([]byte, 0, 8+8+8+8+8*len(g.Data))
	var hdr [8]byte
	put := func(v uint64) {
		binary.LittleEndian.PutUint64(hdr[:], v)
		buf = append(buf, hdr[:]...)
	}
	put(checkpointMagic)
	put(uint64(step))
	put(uint64(g.NX))
	put(uint64(g.NY))
	for _, v := range g.Data {
		put(math.Float64bits(v))
	}
	tmp := path + ".tmp"
	if err := os.WriteFile(tmp, buf, 0o644); err != nil {
		return fmt.Errorf("checkpoint: write %s: %w", tmp, err)
	}
	if err := os.Rename(tmp, path); err != nil {
		return fmt.Errorf("checkpoint: rename to %s: %w", path, err)
	}
	return nil
}

// LoadCheckpoint 从文件恢复状态；文件损坏时返回描述性错误。
func LoadCheckpoint(path string) (*Checkpoint, error) {
	raw, err := os.ReadFile(path)
	if err != nil {
		return nil, fmt.Errorf("checkpoint: read %s: %w", path, err)
	}
	if len(raw) < 32 {
		return nil, fmt.Errorf("checkpoint: %s too short (%d bytes)", path, len(raw))
	}
	u64 := func(off int) uint64 { return binary.LittleEndian.Uint64(raw[off:]) }
	if u64(0) != checkpointMagic {
		return nil, fmt.Errorf("checkpoint: %s bad magic, not a PDE checkpoint file", path)
	}
	step := int(u64(8))
	nx, ny := int(u64(16)), int(u64(24))
	want := 32 + 8*nx*ny
	if len(raw) != want {
		return nil, fmt.Errorf("checkpoint: %s size mismatch: got %d bytes, want %d for %dx%d grid",
			path, len(raw), want, nx, ny)
	}
	cp := &Checkpoint{Step: step, NX: nx, NY: ny, Data: make([]float64, nx*ny)}
	for k := range cp.Data {
		cp.Data[k] = math.Float64frombits(u64(32 + 8*k))
	}
	return cp, nil
}

// Grid 把检查点内容还原为网格。
func (cp *Checkpoint) Grid() *Grid {
	g := &Grid{NX: cp.NX, NY: cp.NY, Data: make([]float64, len(cp.Data))}
	copy(g.Data, cp.Data)
	return g
}
