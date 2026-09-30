package solver

import (
	"encoding/json"
	"fmt"
	"hash/fnv"
	"math"
)

// PhaseSource 是单个阶段内的一个矩形热源（归一化坐标 [0,1]）。
//
// 热源出现在某阶段的 Sources 列表中即表示该阶段开启；Rate=0 等价于关闭
// （该区域这一阶段注入为 0）。多个源区允许重叠，重叠格点的注入强度按各源 Rate 相加。
type PhaseSource struct {
	Name string  `json:"name"` // 可选名称，用于对账定位；同阶段内不得重名
	Rate float64 `json:"rate"` // 源强度（单位时间）
	X0   float64 `json:"x0"`   // 归一化 x 区间起点 [0,1)
	X1   float64 `json:"x1"`   // 归一化 x 区间终点 (0,1]
	Y0   float64 `json:"y0"`   // 归一化 y 区间起点 [0,1)
	Y1   float64 `json:"y1"`   // 归一化 y 区间终点 (0,1]
}

// Phase 是一个连续时间阶段，覆盖绝对步号区间 [StartStep, EndStep)：
// 从绝对步 StartStep 推进到 StartStep+1 的那一步起使用本阶段的 BC/Sources，
// 直到到达绝对步 EndStep（该步属于下一阶段）。
type Phase struct {
	StartStep int           `json:"start_step"`
	EndStep   int           `json:"end_step"`
	BC        Boundaries    `json:"bc"`
	Sources   []PhaseSource `json:"sources"`
}

// PhasePlan 是按步号连续、不重叠排列的阶段序列，必须从步 0 开始覆盖到
// 最后一个阶段的 EndStep；不允许空档或重叠。
type PhasePlan struct {
	Phases []Phase `json:"phases"`
}

// ParsePhasePlanJSON 从 JSON 文本解析阶段计划。边界类型可写字符串
// （dirichlet/neumann/periodic）或旧整数（0/1/2）。本函数只做语法解析，
// 结构合法性仍以 ValidatePlan 为准。
func ParsePhasePlanJSON(raw []byte) (*PhasePlan, error) {
	var p PhasePlan
	if err := json.Unmarshal(raw, &p); err != nil {
		return nil, configErr(fmt.Sprintf("config: parse phase plan json: %v", err))
	}
	return &p, nil
}

// ValidatePlan 在推进前校验阶段计划的结构合法性；非法时返回 ErrKindConfig，
// 与运行期的 ErrKindStability 明确区分。
//
// 规则：至少一个阶段；首阶段必须从步 0 开始；阶段首尾相接（StartStep 必须
// 等于上一阶段 EndStep，不允许重叠或空档）；每阶段至少含一步；每阶段边界
// 类型合法且周期边成对；每个源区坐标与强度为有限、有效数值；同阶段源名不重名。
func ValidatePlan(p *PhasePlan) error {
	if p == nil {
		return configErr("config: phase plan is nil")
	}
	if len(p.Phases) == 0 {
		return configErr("config: phase plan has no phases")
	}
	seen := make(map[string]bool)
	for idx := range p.Phases {
		ph := &p.Phases[idx]
		if idx == 0 && ph.StartStep != 0 {
			return configErr(fmt.Sprintf("config: phase %d must start at step 0, got %d", idx, ph.StartStep))
		}
		if idx > 0 {
			prev := p.Phases[idx-1].EndStep
			if ph.StartStep != prev {
				return configErr(fmt.Sprintf(
					"config: phase %d starts at step %d but previous phase ends at %d (phases must be contiguous, no overlap/gap)",
					idx, ph.StartStep, prev))
			}
		}
		if ph.EndStep <= ph.StartStep {
			return configErr(fmt.Sprintf("config: phase %d has empty/invalid range [%d,%d)",
				idx, ph.StartStep, ph.EndStep))
		}
		if ph.StartStep < 0 || ph.EndStep < 0 {
			return configErr(fmt.Sprintf("config: phase %d negative step index", idx))
		}
		if err := validatePhaseBC(idx, ph.BC); err != nil {
			return err
		}
		clear(seen)
		for si, src := range ph.Sources {
			if src.Name != "" {
				if seen[src.Name] {
					return configErr(fmt.Sprintf("config: phase %d duplicate source name %q", idx, src.Name))
				}
				seen[src.Name] = true
			}
			if !(0 <= src.X0 && src.X0 < src.X1 && src.X1 <= 1 &&
				0 <= src.Y0 && src.Y0 < src.Y1 && src.Y1 <= 1) {
				return configErr(fmt.Sprintf("config: phase %d source %d invalid region [%g,%g]x[%g,%g]",
					idx, si, src.X0, src.X1, src.Y0, src.Y1))
			}
			if math.IsNaN(src.Rate) || math.IsInf(src.Rate, 0) {
				return configErr(fmt.Sprintf("config: phase %d source %d rate must be finite, got %g",
					idx, si, src.Rate))
			}
		}
	}
	return nil
}

// validatePhaseBC 校验单阶段四边边界：类型已知、Dirichlet 值有限、周期边成对。
func validatePhaseBC(idx int, bc Boundaries) error {
	type eb struct {
		name string
		b    Boundary
	}
	for _, e := range []eb{
		{"left", bc.Left}, {"right", bc.Right}, {"top", bc.Top}, {"bottom", bc.Bottom},
	} {
		switch e.b.Type {
		case Dirichlet:
			if math.IsNaN(e.b.Value) || math.IsInf(e.b.Value, 0) {
				return configErr(fmt.Sprintf("config: phase %d %s dirichlet value must be finite, got %g",
					idx, e.name, e.b.Value))
			}
		case Neumann, Periodic:
		default:
			return configErr(fmt.Sprintf("config: phase %d unknown boundary type on %s: %d",
				idx, e.name, e.b.Type))
		}
	}
	if (bc.Left.Type == Periodic) != (bc.Right.Type == Periodic) {
		return configErr(fmt.Sprintf("config: phase %d periodic boundary must be paired on left/right", idx))
	}
	if (bc.Top.Type == Periodic) != (bc.Bottom.Type == Periodic) {
		return configErr(fmt.Sprintf("config: phase %d periodic boundary must be paired on top/bottom", idx))
	}
	return nil
}

// phaseAt 返回绝对步号 step 处即将执行的那一步所属阶段索引
// （满足 Phases[i].StartStep <= step < Phases[i].EndStep）。
func (p *PhasePlan) phaseAt(step int) int {
	for i := range p.Phases {
		if step >= p.Phases[i].StartStep && step < p.Phases[i].EndStep {
			return i
		}
	}
	// 越过计划末端时归为最后一个阶段（调用方负责提前校验目标步号）。
	return len(p.Phases) - 1
}

// Len 返回阶段数量。
func (p *PhasePlan) Len() int { return len(p.Phases) }

// End 返回计划覆盖到的绝对步号（最后一个阶段的 EndStep）。
func (p *PhasePlan) End() int {
	if len(p.Phases) == 0 {
		return 0
	}
	return p.Phases[len(p.Phases)-1].EndStep
}

// writeU64 / writeF64 把数值以确定字节序喂给哈希。
func writeU64(h interface{ Write([]byte) (int, error) }, v uint64) {
	var b [8]byte
	for k := 0; k < 8; k++ {
		b[k] = byte(v >> (8 * uint(k)))
	}
	h.Write(b[:])
}

func writeF64(h interface{ Write([]byte) (int, error) }, v float64) {
	writeU64(h, math.Float64bits(v))
}

// fingerprint 返回阶段计划内容的确定性指纹，供检查点恢复时比对计划是否被改动。
func (p *PhasePlan) fingerprint() uint64 {
	h := fnv.New64a()
	writeU64(h, uint64(len(p.Phases)))
	for _, ph := range p.Phases {
		writeU64(h, uint64(ph.StartStep))
		writeU64(h, uint64(ph.EndStep))
		writeU64(h, uint64(ph.BC.Left.Type))
		writeF64(h, ph.BC.Left.Value)
		writeU64(h, uint64(ph.BC.Right.Type))
		writeF64(h, ph.BC.Right.Value)
		writeU64(h, uint64(ph.BC.Top.Type))
		writeF64(h, ph.BC.Top.Value)
		writeU64(h, uint64(ph.BC.Bottom.Type))
		writeF64(h, ph.BC.Bottom.Value)
		writeU64(h, uint64(len(ph.Sources)))
		for _, s := range ph.Sources {
			h.Write([]byte(s.Name))
			writeF64(h, s.Rate)
			writeF64(h, s.X0)
			writeF64(h, s.X1)
			writeF64(h, s.Y0)
			writeF64(h, s.Y1)
		}
	}
	return h.Sum64()
}

// MarshalJSON 使边界类型以可读字符串序列化。
func (t BoundaryType) MarshalJSON() ([]byte, error) {
	switch t {
	case Dirichlet:
		return json.Marshal("dirichlet")
	case Neumann:
		return json.Marshal("neumann")
	case Periodic:
		return json.Marshal("periodic")
	default:
		return nil, fmt.Errorf("unknown boundary type %d", t)
	}
}

// UnmarshalJSON 同时接受字符串（dirichlet/neumann/periodic）与整数（0/1/2）。
func (t *BoundaryType) UnmarshalJSON(data []byte) error {
	var s string
	if err := json.Unmarshal(data, &s); err == nil {
		switch s {
		case "dirichlet":
			*t = Dirichlet
		case "neumann":
			*t = Neumann
		case "periodic":
			*t = Periodic
		default:
			return configErr(fmt.Sprintf("config: unknown boundary type %q", s))
		}
		return nil
	}
	var n int
	if err := json.Unmarshal(data, &n); err == nil {
		switch BoundaryType(n) {
		case Dirichlet, Neumann, Periodic:
			*t = BoundaryType(n)
		default:
			return configErr(fmt.Sprintf("config: unknown boundary type %d", n))
		}
		return nil
	}
	return configErr(fmt.Sprintf("config: boundary type must be string or int, got %s", data))
}
