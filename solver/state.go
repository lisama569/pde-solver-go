package solver

import (
	"encoding/binary"
	"fmt"
	"math"
	"os"
)

// stateMagic 是带配置签名的状态文件（v2）magic，区别于旧版检查点 v1。
const stateMagic uint64 = 0x5044455354415445 // "PDESTATE"

const stateHeaderLen = 80 // magic+step+nx+ny+hasPlan+planHash+DX+DY+DT+D

// StateSig 是决定数值轨迹的配置指纹。恢复时任一字段不一致都必须在继续
// 计算前明确拒绝（返回 ErrKindConfig），不能悄悄按新配置算。
type StateSig struct {
	NX, NY      int
	DX, DY, DT  float64
	Diffusivity float64
	HasPlan     bool
	PlanHash    uint64 // HasPlan=false 时为 0（旧单套配置模式）
}

// State 是可落盘、可校验一致性的运行状态（检查点 v2）。
type State struct {
	Step int
	NX   int
	NY   int
	Data []float64
	Sig  StateSig
}

// stateSignature 提取配置的轨迹决定性字段；plan 为 nil 时为旧模式。
func stateSignature(cfg *Config, plan *PhasePlan) StateSig {
	sig := StateSig{
		NX: cfg.NX, NY: cfg.NY,
		DX: cfg.DX, DY: cfg.DY, DT: cfg.DT,
		Diffusivity: cfg.Diffusivity,
	}
	if plan != nil {
		sig.HasPlan = true
		sig.PlanHash = plan.fingerprint()
	}
	return sig
}

// SaveState 将状态连同配置签名以二进制格式原子落盘（先写临时文件再改名）。
func SaveState(path string, step int, g *Grid, cfg *Config, plan *PhasePlan) error {
	sig := stateSignature(cfg, plan)
	buf := make([]byte, 0, stateHeaderLen+8*len(g.Data))
	var hdr [8]byte
	put := func(v uint64) {
		binary.LittleEndian.PutUint64(hdr[:], v)
		buf = append(buf, hdr[:]...)
	}
	putF := func(v float64) { put(math.Float64bits(v)) }
	put(stateMagic)
	put(uint64(step))
	put(uint64(g.NX))
	put(uint64(g.NY))
	if sig.HasPlan {
		put(1)
	} else {
		put(0)
	}
	put(sig.PlanHash)
	putF(sig.DX)
	putF(sig.DY)
	putF(sig.DT)
	putF(sig.Diffusivity)
	for _, v := range g.Data {
		putF(v)
	}
	tmp := path + ".tmp"
	if err := os.WriteFile(tmp, buf, 0o644); err != nil {
		return fmt.Errorf("state: write %s: %w", tmp, err)
	}
	if err := os.Rename(tmp, path); err != nil {
		return fmt.Errorf("state: rename to %s: %w", path, err)
	}
	return nil
}

// LoadState 从文件恢复带签名的状态；文件损坏或格式不符时返回描述性错误。
func LoadState(path string) (*State, error) {
	raw, err := os.ReadFile(path)
	if err != nil {
		return nil, fmt.Errorf("state: read %s: %w", path, err)
	}
	if len(raw) < stateHeaderLen {
		return nil, fmt.Errorf("state: %s too short (%d bytes)", path, len(raw))
	}
	u64 := func(off int) uint64 { return binary.LittleEndian.Uint64(raw[off:]) }
	f64 := func(off int) float64 { return math.Float64frombits(u64(off)) }
	if u64(0) == checkpointMagic {
		return nil, fmt.Errorf("state: %s is a legacy v1 checkpoint without config signature; phase-plan resume requires a v2 state saved by the -plan workflow", path)
	}
	if u64(0) != stateMagic {
		return nil, fmt.Errorf("state: %s bad magic, not a PDE state file", path)
	}
	step := int(u64(8))
	nx, ny := int(u64(16)), int(u64(24))
	want := stateHeaderLen + 8*nx*ny
	if len(raw) != want {
		return nil, fmt.Errorf("state: %s size mismatch: got %d bytes, want %d for %dx%d grid",
			path, len(raw), want, nx, ny)
	}
	st := &State{
		Step: step, NX: nx, NY: ny,
		Data: make([]float64, nx*ny),
		Sig: StateSig{
			NX: nx, NY: ny,
			HasPlan:     u64(32) == 1,
			PlanHash:    u64(40),
			DX:          f64(48),
			DY:          f64(56),
			DT:          f64(64),
			Diffusivity: f64(72),
		},
	}
	for k := range st.Data {
		st.Data[k] = f64(stateHeaderLen + 8*k)
	}
	return st, nil
}

// Grid 还原状态中的网格（拷贝，调用方可自由修改）。
func (st *State) Grid() *Grid {
	g := &Grid{NX: st.NX, NY: st.NY, Data: make([]float64, len(st.Data))}
	copy(g.Data, st.Data)
	return g
}

// Verify 校验传入配置/计划与保存状态是否一致；任一字段不一致都返回
// ErrKindConfig 并指出具体差异，确保不会悄悄按新配置继续。
func (st *State) Verify(cfg *Config, plan *PhasePlan) error {
	want := stateSignature(cfg, plan)
	mismatch := func(format string, a ...any) error {
		return configErr(fmt.Sprintf("state: resume rejected, config differs from saved state: "+format, a...))
	}
	if st.Sig.NX != want.NX || st.Sig.NY != want.NY {
		return mismatch("grid saved=%dx%d resume=%dx%d", st.Sig.NX, st.Sig.NY, want.NX, want.NY)
	}
	if st.Sig.DX != want.DX || st.Sig.DY != want.DY {
		return mismatch("spacing saved dx=%g dy=%g resume dx=%g dy=%g", st.Sig.DX, st.Sig.DY, want.DX, want.DY)
	}
	if st.Sig.DT != want.DT {
		return mismatch("dt saved=%g resume=%g", st.Sig.DT, want.DT)
	}
	if st.Sig.Diffusivity != want.Diffusivity {
		return mismatch("diffusivity saved=%g resume=%g", st.Sig.Diffusivity, want.Diffusivity)
	}
	if st.Sig.HasPlan != want.HasPlan {
		return mismatch("phase plan presence saved hasPlan=%v resume hasPlan=%v", st.Sig.HasPlan, want.HasPlan)
	}
	if st.Sig.HasPlan && st.Sig.PlanHash != want.PlanHash {
		return mismatch("phase plan content differs (saved hash=%016x resume hash=%016x)",
			st.Sig.PlanHash, want.PlanHash)
	}
	return nil
}
