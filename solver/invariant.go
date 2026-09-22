package solver

// InvariantReport 守恒量校验报告。
type InvariantReport struct {
	Step      int
	Total     float64 // 当前总量
	Expected  float64 // 理论预期总量
	Deviation float64 // Total - Expected
}

// TotalMass 计算全场积分量 Σ u·dx·dy（全 Neumann + 无源时为守恒量）。
func TotalMass(g *Grid, dx, dy float64) float64 {
	sum := 0.0
	for _, v := range g.Data {
		sum += v
	}
	return sum * dx * dy
}

// CheckInvariant 校验守恒量偏差，报告可定位到具体步号。
func CheckInvariant(g *Grid, cfg *Config, expected float64, step int) InvariantReport {
	total := TotalMass(g, cfg.DX, cfg.DY)
	return InvariantReport{
		Step:      step,
		Total:     total,
		Expected:  expected,
		Deviation: total - expected,
	}
}
