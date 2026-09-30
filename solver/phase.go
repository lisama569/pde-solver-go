package solver

import (
	"encoding/binary"
	"fmt"
	"hash/fnv"
	"math"
)

// Phase 描述一个时间阶段：步号区间 [From, To) 内生效的边界与源项配置。
type Phase struct {
	From    int        // 起始步号（含）
	To      int        // 结束步号（不含）
	BC      Boundaries // 本阶段四边边界
	Sources []Source   // 本阶段源项列表（可多个，可为空）
}

// PhasePlan 是按时间步顺序排列的连续、不重叠阶段序列。
// 为 nil 时按旧的单套 BC/Src 配置运行，行为不变。
type PhasePlan struct {
	Phases []Phase
}

// validateBoundaries 校验一组四边边界配置（类型合法、周期成对、Dirichlet 值有限）。
func validateBoundaries(bc Boundaries) error {
	for name, b := range map[string]Boundary{
		"left": bc.Left, "right": bc.Right, "top": bc.Top, "bottom": bc.Bottom,
	} {
		if b.Type != Dirichlet && b.Type != Neumann && b.Type != Periodic {
			return configErr(fmt.Sprintf("config: unknown boundary type on %s: %d", name, b.Type))
		}
		if b.Type == Dirichlet && (math.IsNaN(b.Value) || math.IsInf(b.Value, 0)) {
			return configErr(fmt.Sprintf("config: non-finite Dirichlet value on %s: %v", name, b.Value))
		}
	}
	if (bc.Left.Type == Periodic) != (bc.Right.Type == Periodic) {
		return configErr("config: periodic boundary must be paired on left/right")
	}
	if (bc.Top.Type == Periodic) != (bc.Bottom.Type == Periodic) {
		return configErr("config: periodic boundary must be paired on top/bottom")
	}
	return nil
}

// validateSource 校验单个源项（区域在 [0,1] 内且非空、强度有限）。
func validateSource(s *Source) error {
	if math.IsNaN(s.Rate) || math.IsInf(s.Rate, 0) {
		return configErr(fmt.Sprintf("config: non-finite source rate: %v", s.Rate))
	}
	if !(0 <= s.X0 && s.X0 < s.X1 && s.X1 <= 1 && 0 <= s.Y0 && s.Y0 < s.Y1 && s.Y1 <= 1) {
		return configErr(fmt.Sprintf("config: invalid source region: [%g,%g]x[%g,%g]", s.X0, s.X1, s.Y0, s.Y1))
	}
	return nil
}

// Validate 校验阶段计划合法性，返回 ErrKindConfig 错误：
//   - 至少一个阶段（存在有效时间范围）；
//   - 首阶段从 0 开始，末阶段覆盖到 totalSteps（不越界、不欠覆盖）；
//   - 相邻阶段首尾相接（无间隙、无重叠），每阶段长度 > 0；
//   - 每阶段的边界与源项配置合法、数值有限。
func (p *PhasePlan) Validate(totalSteps int) error {
	if p == nil || len(p.Phases) == 0 {
		return configErr("config: phase plan has no phases (no valid time range)")
	}
	prevTo := 0
	for k, ph := range p.Phases {
		if ph.From < 0 || ph.To < 0 || ph.From >= ph.To {
			return configErr(fmt.Sprintf("config: phase %d has invalid range [%d,%d)", k, ph.From, ph.To))
		}
		if ph.From != prevTo {
			return configErr(fmt.Sprintf(
				"config: phase %d starts at %d, want %d (phases must be contiguous and non-overlapping)",
				k, ph.From, prevTo))
		}
		if err := validateBoundaries(ph.BC); err != nil {
			return configErr(fmt.Sprintf("config: phase %d: %v", k, err))
		}
		for si := range ph.Sources {
			if err := validateSource(&ph.Sources[si]); err != nil {
				return configErr(fmt.Sprintf("config: phase %d source %d: %v", k, si, err))
			}
		}
		prevTo = ph.To
	}
	if last := p.Phases[len(p.Phases)-1]; last.To != totalSteps {
		return configErr(fmt.Sprintf(
			"config: phase plan ends at step %d, want %d (plan must cover [0, Steps) exactly)",
			last.To, totalSteps))
	}
	return nil
}

// PhaseAt 返回覆盖步号 step 的阶段索引；无覆盖返回 -1。
func (p *PhasePlan) PhaseAt(step int) int {
	for k := range p.Phases {
		if step >= p.Phases[k].From && step < p.Phases[k].To {
			return k
		}
	}
	return -1
}

// Fingerprint 返回阶段计划内容的 FNV-1a 哈希，用于检查点一致性校验。
func (p *PhasePlan) Fingerprint() uint64 {
	h := fnv.New64a()
	var buf [8]byte
	putU64 := func(v uint64) {
		binary.LittleEndian.PutUint64(buf[:], v)
		h.Write(buf[:])
	}
	putF64 := func(v float64) { putU64(math.Float64bits(v)) }
	putBC := func(bc Boundaries) {
		for _, b := range []Boundary{bc.Left, bc.Right, bc.Top, bc.Bottom} {
			putU64(uint64(b.Type))
			putF64(b.Value)
		}
	}
	putSrc := func(s *Source) {
		putF64(s.Rate)
		putF64(s.X0)
		putF64(s.X1)
		putF64(s.Y0)
		putF64(s.Y1)
	}
	putU64(uint64(len(p.Phases)))
	for k := range p.Phases {
		ph := &p.Phases[k]
		putU64(uint64(ph.From))
		putU64(uint64(ph.To))
		putBC(ph.BC)
		putU64(uint64(len(ph.Sources)))
		for si := range ph.Sources {
			putSrc(&ph.Sources[si])
		}
	}
	return h.Sum64()
}

// EffectiveFingerprint 返回配置的"有效阶段计划"指纹：
// 有计划时为计划指纹；无计划（旧配置）时为单套 BC/Src 的指纹，
// 使得断点续算对两类配置都能校验一致性。
func (c *Config) EffectiveFingerprint() uint64 {
	if c.Plan != nil {
		return c.Plan.Fingerprint()
	}
	legacy := &PhasePlan{Phases: []Phase{{From: 0, To: 1, BC: c.BC}}}
	if c.Src != nil {
		legacy.Phases[0].Sources = []Source{*c.Src}
	}
	return legacy.Fingerprint()
}
