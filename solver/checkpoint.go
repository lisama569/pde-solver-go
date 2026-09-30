package solver

import (
	"encoding/binary"
	"fmt"
	"math"
	"os"
)

// 检查点文件格式标识（低字节为版本号）。
const (
	checkpointMagicV1 uint64 = 0x504445434B505431 // "PDECKPT1"：magic + 步号 + 维度 + 数据
	checkpointMagicV2 uint64 = 0x504445434B505432 // "PDECKPT2"：V1 + 网格/步长/系数 + 计划指纹
)

// Checkpoint 是可落盘的运行状态。
//
// V2 格式额外记录保存时的 DX/DY/DT/Diffusivity 与有效阶段计划指纹
// （HasFingerprint 为真），恢复时可用 VerifyConfig 在继续计算前
// 拒绝与保存时不一致的配置。V1 文件不含这些信息，只能校验网格维度。
type Checkpoint struct {
	Step int
	NX   int
	NY   int
	Data []float64

	// 以下为 V2 格式携带的配置指纹（V1 文件为零值且 HasFingerprint=false）
	DX, DY         float64
	DT             float64
	Diffusivity    float64
	PlanHash       uint64
	HasFingerprint bool
}

// SaveCheckpoint 以 V1 格式落盘（不含配置指纹，仅供兼容旧流程）。
// 新流程应使用 SaveCheckpointWithConfig 以获得恢复时的配置校验能力。
func SaveCheckpoint(path string, step int, g *Grid) error {
	return saveCheckpoint(path, step, g, nil)
}

// SaveCheckpointWithConfig 以 V2 格式落盘：除场数据外还记录
// DX/DY/DT/Diffusivity 与 cfg.EffectiveFingerprint()。
func SaveCheckpointWithConfig(path string, step int, g *Grid, cfg *Config) error {
	return saveCheckpoint(path, step, g, cfg)
}

// saveCheckpoint 写检查点（先写临时文件再原子改名）；cfg 非空时写 V2 格式。
func saveCheckpoint(path string, step int, g *Grid, cfg *Config) error {
	buf := make([]byte, 0, 80+8*len(g.Data))
	var hdr [8]byte
	put := func(v uint64) {
		binary.LittleEndian.PutUint64(hdr[:], v)
		buf = append(buf, hdr[:]...)
	}
	if cfg == nil {
		put(checkpointMagicV1)
	} else {
		put(checkpointMagicV2)
	}
	put(uint64(step))
	put(uint64(g.NX))
	put(uint64(g.NY))
	if cfg != nil {
		put(math.Float64bits(cfg.DX))
		put(math.Float64bits(cfg.DY))
		put(math.Float64bits(cfg.DT))
		put(math.Float64bits(cfg.Diffusivity))
		put(cfg.EffectiveFingerprint())
	}
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

// LoadCheckpoint 从文件恢复状态（自动识别 V1/V2）；文件损坏时返回描述性错误。
func LoadCheckpoint(path string) (*Checkpoint, error) {
	raw, err := os.ReadFile(path)
	if err != nil {
		return nil, fmt.Errorf("checkpoint: read %s: %w", path, err)
	}
	if len(raw) < 32 {
		return nil, fmt.Errorf("checkpoint: %s too short (%d bytes)", path, len(raw))
	}
	u64 := func(off int) uint64 { return binary.LittleEndian.Uint64(raw[off:]) }
	magic := u64(0)
	if magic != checkpointMagicV1 && magic != checkpointMagicV2 {
		return nil, fmt.Errorf("checkpoint: %s bad magic, not a PDE checkpoint file", path)
	}
	cp := &Checkpoint{Step: int(u64(8)), NX: int(u64(16)), NY: int(u64(24))}
	off := 32
	if magic == checkpointMagicV2 {
		if len(raw) < 72 {
			return nil, fmt.Errorf("checkpoint: %s too short for v2 header (%d bytes)", path, len(raw))
		}
		cp.DX = math.Float64frombits(u64(32))
		cp.DY = math.Float64frombits(u64(40))
		cp.DT = math.Float64frombits(u64(48))
		cp.Diffusivity = math.Float64frombits(u64(56))
		cp.PlanHash = u64(64)
		cp.HasFingerprint = true
		off = 72
	}
	want := off + 8*cp.NX*cp.NY
	if len(raw) != want {
		return nil, fmt.Errorf("checkpoint: %s size mismatch: got %d bytes, want %d for %dx%d grid",
			path, len(raw), want, cp.NX, cp.NY)
	}
	cp.Data = make([]float64, cp.NX*cp.NY)
	for k := range cp.Data {
		cp.Data[k] = math.Float64frombits(u64(off + 8*k))
	}
	return cp, nil
}

// VerifyConfig 在断点续算前校验配置与保存时是否一致，不一致返回
// ErrKindConfig 错误（绝不按新配置静默继续）。校验项：
//   - 网格维度 NX/NY（始终校验）；
//   - DX/DY/DT/Diffusivity 与有效阶段计划指纹（仅 V2 检查点）。
func (cp *Checkpoint) VerifyConfig(cfg *Config) error {
	if cp.NX != cfg.NX || cp.NY != cfg.NY {
		return configErr(fmt.Sprintf(
			"checkpoint: grid mismatch: saved %dx%d, config %dx%d", cp.NX, cp.NY, cfg.NX, cfg.NY))
	}
	if !cp.HasFingerprint {
		return nil // V1 文件无更多信息可校验
	}
	if cp.DX != cfg.DX || cp.DY != cfg.DY {
		return configErr(fmt.Sprintf(
			"checkpoint: spacing mismatch: saved (%g,%g), config (%g,%g)", cp.DX, cp.DY, cfg.DX, cfg.DY))
	}
	if cp.DT != cfg.DT {
		return configErr(fmt.Sprintf("checkpoint: DT mismatch: saved %g, config %g", cp.DT, cfg.DT))
	}
	if cp.Diffusivity != cfg.Diffusivity {
		return configErr(fmt.Sprintf(
			"checkpoint: diffusivity mismatch: saved %g, config %g", cp.Diffusivity, cfg.Diffusivity))
	}
	if h := cfg.EffectiveFingerprint(); cp.PlanHash != h {
		return configErr(fmt.Sprintf(
			"checkpoint: phase plan / boundary / source mismatch: saved hash %#x, config hash %#x",
			cp.PlanHash, h))
	}
	return nil
}

// Grid 把检查点内容还原为网格。
func (cp *Checkpoint) Grid() *Grid {
	g := &Grid{NX: cp.NX, NY: cp.NY, Data: make([]float64, len(cp.Data))}
	copy(g.Data, cp.Data)
	return g
}
