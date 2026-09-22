package solver

import "testing"

// baseConfig 返回一个合法配置作为测试基准。
func baseConfig() *Config {
	return &Config{
		NX: 32, NY: 24,
		DX: 0.01, DY: 0.01, DT: 1e-5,
		Diffusivity: 1.0,
		Steps:       10,
		BC: Boundaries{
			Left:   Boundary{Type: Dirichlet, Value: 100},
			Right:  Boundary{Type: Neumann},
			Top:    Boundary{Type: Neumann},
			Bottom: Boundary{Type: Dirichlet, Value: 0},
		},
		Src: &Source{Rate: 50, X0: 0.4, X1: 0.6, Y0: 0.4, Y1: 0.6},
	}
}

func TestInvalidConfigsAreRejected(t *testing.T) {
	cases := []struct {
		name   string
		mutate func(*Config)
	}{
		{"grid too small", func(c *Config) { c.NX = 2 }},
		{"negative spacing", func(c *Config) { c.DX = -0.01 }},
		{"zero dt", func(c *Config) { c.DT = 0 }},
		{"negative diffusivity", func(c *Config) { c.Diffusivity = -1 }},
		{"negative steps", func(c *Config) { c.Steps = -1 }},
		{"unknown boundary type", func(c *Config) { c.BC.Left.Type = BoundaryType(99) }},
		{"unpaired periodic x", func(c *Config) { c.BC.Left.Type = Periodic }},
		{"unpaired periodic y", func(c *Config) { c.BC.Top.Type = Periodic }},
		{"bad source region", func(c *Config) { c.Src.X0, c.Src.X1 = 0.8, 0.2 }},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			cfg := baseConfig()
			tc.mutate(cfg)
			err := cfg.Validate()
			t.Logf("input=%+v err=%v", *cfg, err)
			if err == nil {
				t.Fatalf("case %q: expected config error, got nil", tc.name)
			}
			if !IsConfigError(err) {
				t.Fatalf("case %q: expected ErrKindConfig, got %v", tc.name, err)
			}
			if IsStabilityError(err) {
				t.Fatalf("case %q: config error misclassified as stability error", tc.name)
			}
		})
	}
}

func TestStabilityLimitEnforced(t *testing.T) {
	cfg := baseConfig()
	lim := cfg.StabilityLimit()
	t.Logf("DX=%g DY=%g D=%g -> dt_limit=%g", cfg.DX, cfg.DY, cfg.Diffusivity, lim)

	// 恰好不超限：允许
	cfg.DT = lim
	if err := cfg.CheckStability(); err != nil {
		t.Fatalf("DT at limit should be accepted, got %v", err)
	}
	// 超限：必须给出 ErrKindStability，且与配置错误可区分
	cfg.DT = lim * 1.5
	err := cfg.CheckStability()
	t.Logf("DT=%g (1.5x limit) err=%v", cfg.DT, err)
	if err == nil {
		t.Fatal("expected stability error for oversized DT")
	}
	if !IsStabilityError(err) {
		t.Fatalf("expected ErrKindStability, got %v", err)
	}
	if IsConfigError(err) {
		t.Fatal("stability violation misclassified as config error")
	}
}

func TestValidConfigAccepted(t *testing.T) {
	cfg := baseConfig()
	if err := cfg.Validate(); err != nil {
		t.Fatalf("valid config rejected: %v", err)
	}
	if err := cfg.CheckStability(); err != nil {
		t.Fatalf("valid config failed stability check: %v", err)
	}
}
