package solver

import "fmt"

// BoundaryType 边界条件类型。
type BoundaryType int

const (
	Dirichlet BoundaryType = iota // 固定值
	Neumann                       // 零法向通量（绝热）
	Periodic                      // 周期边界
)

// Boundary 单边边界条件。
type Boundary struct {
	Type  BoundaryType
	Value float64 // Dirichlet 时的固定值
}

// Boundaries 四边边界条件。
type Boundaries struct {
	Left, Right, Top, Bottom Boundary
}

// Source 矩形常源项配置（归一化坐标）。
type Source struct {
	Rate float64 // 源强度（单位时间）
	X0   float64 // 源区域归一化起点 [0,1)
	X1   float64 // 源区域归一化终点 (0,1]
	Y0   float64
	Y1   float64
}

// Config 求解配置。
type Config struct {
	NX, NY      int        // 网格点数（含边界）
	DX, DY      float64    // 空间步长
	DT          float64    // 时间步长
	Diffusivity float64    // 扩散系数 D
	Steps       int        // 时间步数
	BC          Boundaries // 边界条件
	Src         *Source    // 源项，nil 表示无源
	Plan        *PhasePlan // 阶段计划，nil 表示按单套 BC/Src 运行（旧行为）
}

// Validate 校验配置合法性，返回 ErrKindConfig 错误。
func (c *Config) Validate() error {
	if c.NX < 3 || c.NY < 3 {
		return configErr(fmt.Sprintf("config: grid too small: NX=%d NY=%d, need >= 3", c.NX, c.NY))
	}
	if c.DX <= 0 || c.DY <= 0 {
		return configErr(fmt.Sprintf("config: spacing must be positive: DX=%g DY=%g", c.DX, c.DY))
	}
	if c.DT <= 0 {
		return configErr(fmt.Sprintf("config: DT must be positive: DT=%g", c.DT))
	}
	if c.Diffusivity < 0 {
		return configErr(fmt.Sprintf("config: diffusivity must be >= 0: D=%g", c.Diffusivity))
	}
	if c.Steps < 0 {
		return configErr(fmt.Sprintf("config: steps must be >= 0: Steps=%d", c.Steps))
	}
	if err := validateBoundaries(c.BC); err != nil {
		return err
	}
	if c.Src != nil {
		if err := validateSource(c.Src); err != nil {
			return err
		}
	}
	if c.Plan != nil {
		if err := c.Plan.Validate(c.Steps); err != nil {
			return err
		}
	}
	return nil
}

// StabilityLimit 返回二维显式 FTCS 格式允许的最大 dt：dt <= 1/(2D(1/dx²+1/dy²))。
func (c *Config) StabilityLimit() float64 {
	if c.Diffusivity == 0 {
		return 1e308 // 无扩散时仅源项，不受 CFL 限制
	}
	return 1.0 / (2.0 * c.Diffusivity * (1.0/(c.DX*c.DX) + 1.0/(c.DY*c.DY)))
}

// CheckStability 校验稳定性条件，违反时返回 ErrKindStability 错误。
func (c *Config) CheckStability() error {
	lim := c.StabilityLimit()
	if c.DT > lim {
		return stabilityErr(fmt.Sprintf(
			"stability: DT=%g exceeds explicit-scheme limit %g (D=%g, DX=%g, DY=%g)",
			c.DT, lim, c.Diffusivity, c.DX, c.DY))
	}
	return nil
}
