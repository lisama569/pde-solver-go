package solver

import (
	"math"
	"path/filepath"
	"testing"
)

// phaseConfig 返回带三阶段计划的配置（总 60 步）：
//
//	阶段0 [0,20)：中心矩形源加热，左边界 Dirichlet 100，底边 Dirichlet 0
//	阶段1 [20,40)：关源冷却，四边全 Neumann（绝热）
//	阶段2 [40,60)：左边界改为 Dirichlet 200，启用两个新源区
func phaseConfig() *Config {
	cfg := baseConfig()
	cfg.Steps = 60
	cfg.Plan = &PhasePlan{Phases: []Phase{
		{
			From: 0, To: 20,
			BC: Boundaries{
				Left:   Boundary{Type: Dirichlet, Value: 100},
				Right:  Boundary{Type: Neumann},
				Top:    Boundary{Type: Neumann},
				Bottom: Boundary{Type: Dirichlet, Value: 0},
			},
			Sources: []Source{{Rate: 50, X0: 0.4, X1: 0.6, Y0: 0.4, Y1: 0.6}},
		},
		{
			From: 20, To: 40,
			BC: Boundaries{
				Left:   Boundary{Type: Neumann},
				Right:  Boundary{Type: Neumann},
				Top:    Boundary{Type: Neumann},
				Bottom: Boundary{Type: Neumann},
			},
		},
		{
			From: 40, To: 60,
			BC: Boundaries{
				Left:   Boundary{Type: Dirichlet, Value: 200},
				Right:  Boundary{Type: Neumann},
				Top:    Boundary{Type: Neumann},
				Bottom: Boundary{Type: Neumann},
			},
			Sources: []Source{
				{Rate: 30, X0: 0.1, X1: 0.2, Y0: 0.1, Y1: 0.2},
				{Rate: 70, X0: 0.7, X1: 0.9, Y0: 0.7, Y1: 0.9},
			},
		},
	}}
	return cfg
}

// countSourceCells 按 sourceAt 的同一判定统计源区域覆盖的格点数。
func countSourceCells(nx, ny int, s *Source) int {
	n := 0
	for j := 0; j < ny; j++ {
		y := float64(j) / float64(ny-1)
		for i := 0; i < nx; i++ {
			x := float64(i) / float64(nx-1)
			if x >= s.X0 && x <= s.X1 && y >= s.Y0 && y <= s.Y1 {
				n++
			}
		}
	}
	return n
}

// TestPhasePlanValidation 验证合法计划被接受、PhaseAt 定位正确。
func TestPhasePlanValidation(t *testing.T) {
	cfg := phaseConfig()
	if err := cfg.Validate(); err != nil {
		t.Fatalf("valid phase plan rejected: %v", err)
	}
	for step, want := range map[int]int{0: 0, 19: 0, 20: 1, 39: 1, 40: 2, 59: 2, 60: -1, -1: -1} {
		if got := cfg.Plan.PhaseAt(step); got != want {
			t.Fatalf("PhaseAt(%d) = %d, want %d", step, got, want)
		}
	}
	t.Logf("plan: 3 phases covering [0,60), PhaseAt boundary steps verified")
}

// TestInvalidPhasePlans 验证各类非法阶段计划在推进前被拒绝且归类为配置错误。
func TestInvalidPhasePlans(t *testing.T) {
	valid := func() []Phase { return phaseConfig().Plan.Phases }
	cases := []struct {
		name   string
		mutate func(*Config)
	}{
		{"empty plan", func(c *Config) { c.Plan = &PhasePlan{} }},
		{"gap between phases", func(c *Config) {
			c.Plan.Phases[1].From = 25
			c.Plan.Phases[1].To = 45
			c.Plan.Phases[2].From = 45
		}},
		{"overlap between phases", func(c *Config) {
			c.Plan.Phases[0].To = 25
			c.Plan.Phases[2].From = 45
		}},
		{"zero length phase", func(c *Config) {
			c.Plan.Phases[0].To = 0
			c.Plan.Phases[1].From = 0
		}},
		{"negative from", func(c *Config) { c.Plan.Phases[0].From = -5 }},
		{"plan ends early", func(c *Config) { c.Plan.Phases[2].To = 55 }},
		{"plan exceeds steps", func(c *Config) { c.Plan.Phases[2].To = 70 }},
		{"plan not starting at 0", func(c *Config) {
			c.Plan.Phases = valid()[1:] // 首阶段从 20 开始，[0,20) 无覆盖
		}},
		{"bad boundary type in phase", func(c *Config) {
			c.Plan.Phases[1].BC.Left.Type = BoundaryType(99)
		}},
		{"unpaired periodic in phase", func(c *Config) {
			c.Plan.Phases[1].BC.Left.Type = Periodic
		}},
		{"NaN source rate", func(c *Config) {
			c.Plan.Phases[0].Sources[0].Rate = math.NaN()
		}},
		{"Inf Dirichlet value", func(c *Config) {
			c.Plan.Phases[2].BC.Left.Value = math.Inf(1)
		}},
		{"bad source region in phase", func(c *Config) {
			c.Plan.Phases[2].Sources[0].X1 = 1.5
		}},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			cfg := phaseConfig()
			tc.mutate(cfg)
			err := cfg.Validate()
			t.Logf("case=%q err=%v", tc.name, err)
			if err == nil {
				t.Fatalf("case %q: expected config error, got nil", tc.name)
			}
			if !IsConfigError(err) {
				t.Fatalf("case %q: expected ErrKindConfig, got %v", tc.name, err)
			}
			if IsStabilityError(err) {
				t.Fatalf("case %q: plan problem misclassified as stability error", tc.name)
			}
		})
	}
}

// TestPhasePlanStabilityDistinct 验证带合法计划但 dt 超限时归类为失稳而非配置错误。
func TestPhasePlanStabilityDistinct(t *testing.T) {
	cfg := phaseConfig()
	cfg.DT = cfg.StabilityLimit() * 2
	_, _, err := RunPhases(cfg, testInitial(cfg.NX, cfg.NY), 1, 0, cfg.Steps)
	t.Logf("DT=2x limit with valid plan: err=%v", err)
	if err == nil {
		t.Fatal("expected stability error")
	}
	if !IsStabilityError(err) || IsConfigError(err) {
		t.Fatalf("expected pure ErrKindStability, got %v", err)
	}
}

// TestPhaseSwitchingAndBudget 验证阶段切换语义与逐阶段收支对账闭合。
func TestPhaseSwitchingAndBudget(t *testing.T) {
	cfg := phaseConfig()
	initial := testInitial(cfg.NX, cfg.NY)
	final, reports, err := RunPhases(cfg, initial, 1, 0, cfg.Steps)
	if err != nil {
		t.Fatalf("RunPhases: %v", err)
	}
	if len(reports) != 3 {
		t.Fatalf("got %d phase reports, want 3", len(reports))
	}
	for k, r := range reports {
		t.Logf("phase %d steps[%d,%d): mass %.6f -> %.6f (delta %+.6e), srcTotal %+.6e, "+
			"boundary{L:%+.4e R:%+.4e T:%+.4e B:%+.4e}, closure %+.3e tol %.3e balanced=%v",
			r.PhaseIndex, r.FromStep, r.ToStep, r.MassStart, r.MassEnd, r.DeltaMass(),
			r.SourceTotal(), r.Boundary.Left, r.Boundary.Right, r.Boundary.Top, r.Boundary.Bottom,
			r.Closure(), r.ClosureTol(), r.Balanced())
		if r.PhaseIndex != k {
			t.Fatalf("report %d: PhaseIndex = %d", k, r.PhaseIndex)
		}
		if !r.Balanced() {
			t.Fatalf("phase %d budget does not close: residual %g exceeds tol %g",
				k, r.Closure(), r.ClosureTol())
		}
	}
	// 阶段边界检查
	wantRanges := [][2]int{{0, 20}, {20, 40}, {40, 60}}
	for k, r := range reports {
		if r.FromStep != wantRanges[k][0] || r.ToStep != wantRanges[k][1] {
			t.Fatalf("phase %d range = [%d,%d), want %v", k, r.FromStep, r.ToStep, wantRanges[k])
		}
	}
	// 冷却阶段（全 Neumann 无源）：总量应守恒
	cool := reports[1]
	if len(cool.SourceInject) != 0 {
		t.Fatalf("cooling phase should have no sources, got %d", len(cool.SourceInject))
	}
	if d := math.Abs(cool.DeltaMass()); d > cool.ClosureTol() {
		t.Fatalf("cooling phase mass changed by %g (tol %g)", d, cool.ClosureTol())
	}
	t.Logf("cooling phase: |deltaMass|=%.3e (insulated, no source)", math.Abs(cool.DeltaMass()))
	// 加热阶段：源注入为正且总量上升
	if reports[0].SourceTotal() <= 0 || reports[0].DeltaMass() <= 0 {
		t.Fatalf("heating phase: srcTotal=%g deltaMass=%g, want both positive",
			reports[0].SourceTotal(), reports[0].DeltaMass())
	}
	// 阶段2 的左边界 Dirichlet 200 应已写入边界格
	for j := 0; j < cfg.NY; j++ {
		if v := final.At(0, j); v != 200 {
			t.Fatalf("phase 2 left boundary at j=%d = %g, want 200 (Dirichlet switch not applied)", j, v)
		}
	}
	t.Logf("phase 2 Dirichlet switch: left boundary == 200 everywhere")
}

// TestMultiSourceAccounting 验证多源区注入量与解析值一致。
func TestMultiSourceAccounting(t *testing.T) {
	cfg := phaseConfig()
	_, reports, err := RunPhases(cfg, testInitial(cfg.NX, cfg.NY), 1, 0, cfg.Steps)
	if err != nil {
		t.Fatalf("RunPhases: %v", err)
	}
	r := reports[2] // 阶段2 有两个源区
	if len(r.SourceInject) != 2 {
		t.Fatalf("phase 2: got %d source entries, want 2", len(r.SourceInject))
	}
	nSteps := float64(r.ToStep - r.FromStep)
	for si, src := range cfg.Plan.Phases[2].Sources {
		ncells := countSourceCells(cfg.NX, cfg.NY, &src)
		want := src.Rate * float64(ncells) * cfg.DX * cfg.DY * cfg.DT * nSteps
		dev := math.Abs(r.SourceInject[si] - want)
		t.Logf("phase 2 source %d: rate=%g cells=%d injected=%.6e analytic=%.6e dev=%.3e",
			si, src.Rate, ncells, r.SourceInject[si], want, dev)
		if dev > 1e-9*math.Abs(want) {
			t.Fatalf("source %d injection %g deviates from analytic %g by %g",
				si, r.SourceInject[si], want, dev)
		}
	}
}

// TestPhaseParallelConsistency 验证带阶段计划时不同子域数结果逐位一致。
func TestPhaseParallelConsistency(t *testing.T) {
	cfg := phaseConfig()
	initial := testInitial(cfg.NX, cfg.NY)
	base, _, err := RunPhases(cfg, initial, 1, 0, cfg.Steps)
	if err != nil {
		t.Fatalf("baseline: %v", err)
	}
	for _, workers := range []int{2, 3, 5, 8} {
		got, _, err := RunPhases(cfg, initial, workers, 0, cfg.Steps)
		if err != nil {
			t.Fatalf("workers=%d: %v", workers, err)
		}
		diff := base.MaxAbsDiff(got)
		t.Logf("workers=%d steps=%d maxAbsDiff=%.3e tolerance=0 (bit-identical)", workers, cfg.Steps, diff)
		if diff != 0 {
			t.Fatalf("workers=%d: result deviates from single-domain baseline by %g", workers, diff)
		}
	}
}

// runPhaseResume 在 split 步处保存并以 workers 个子域恢复，返回最终网格。
func runPhaseResume(t *testing.T, cfg *Config, split, workers int) *Grid {
	t.Helper()
	dir := t.TempDir()
	path := filepath.Join(dir, "ckpt.bin")
	initial := testInitial(cfg.NX, cfg.NY)

	part, _, err := RunPhases(cfg, initial, 1, 0, split)
	if err != nil {
		t.Fatalf("first segment: %v", err)
	}
	if err := SaveCheckpointWithConfig(path, split, part, cfg); err != nil {
		t.Fatalf("save: %v", err)
	}
	cp, err := LoadCheckpoint(path)
	if err != nil {
		t.Fatalf("load: %v", err)
	}
	if !cp.HasFingerprint {
		t.Fatal("v2 checkpoint should carry config fingerprint")
	}
	if err := cp.VerifyConfig(cfg); err != nil {
		t.Fatalf("verify: %v", err)
	}
	final, _, err := RunPhases(cfg, cp.Grid(), workers, cp.Step, cfg.Steps-cp.Step)
	if err != nil {
		t.Fatalf("resumed segment: %v", err)
	}
	return final
}

// TestPhaseCheckpointResume 验证阶段内部与阶段边界两种位置的保存/恢复
// 都与不间断运行逐位一致，且恢复后改变子域数仍与单子域基线一致。
func TestPhaseCheckpointResume(t *testing.T) {
	cfg := phaseConfig()
	full, _, err := RunPhases(cfg, testInitial(cfg.NX, cfg.NY), 1, 0, cfg.Steps)
	if err != nil {
		t.Fatalf("uninterrupted run: %v", err)
	}
	for _, tc := range []struct {
		name    string
		split   int
		workers int
	}{
		{"inside phase 1", 25, 1},
		{"inside phase 1, parallel resume", 25, 4},
		{"at phase boundary", 40, 1},
		{"at phase boundary, parallel resume", 40, 3},
		{"inside phase 0", 7, 2},
	} {
		t.Run(tc.name, func(t *testing.T) {
			got := runPhaseResume(t, cfg, tc.split, tc.workers)
			diff := full.MaxAbsDiff(got)
			t.Logf("split@%d workers=%d maxAbsDiff=%.3e tolerance=0 (bit-identical)",
				tc.split, tc.workers, diff)
			if diff != 0 {
				t.Fatalf("resumed run deviates from uninterrupted run by %g", diff)
			}
		})
	}
}

// TestResumeConfigMismatch 恢复时配置与保存不一致必须在继续计算前拒绝。
func TestResumeConfigMismatch(t *testing.T) {
	cfg := phaseConfig()
	path := filepath.Join(t.TempDir(), "ckpt.bin")
	part, _, err := RunPhases(cfg, testInitial(cfg.NX, cfg.NY), 1, 0, 25)
	if err != nil {
		t.Fatalf("first segment: %v", err)
	}
	if err := SaveCheckpointWithConfig(path, 25, part, cfg); err != nil {
		t.Fatalf("save: %v", err)
	}
	cp, err := LoadCheckpoint(path)
	if err != nil {
		t.Fatalf("load: %v", err)
	}

	cases := []struct {
		name   string
		mutate func(*Config)
	}{
		{"grid mismatch", func(c *Config) { c.NX = 64 }},
		{"dt mismatch", func(c *Config) { c.DT *= 1.1 }},
		{"dx mismatch", func(c *Config) { c.DX *= 2 }},
		{"diffusivity mismatch", func(c *Config) { c.Diffusivity = 2 }},
		{"plan boundary moved", func(c *Config) {
			c.Plan.Phases[0].To = 15
			c.Plan.Phases[1].From = 15
		}},
		{"source rate changed", func(c *Config) { c.Plan.Phases[0].Sources[0].Rate = 99 }},
		{"boundary value changed", func(c *Config) { c.Plan.Phases[2].BC.Left.Value = 150 }},
		{"plan dropped", func(c *Config) { c.Plan = nil }},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			bad := phaseConfig()
			tc.mutate(bad)
			err := cp.VerifyConfig(bad)
			t.Logf("case=%q err=%v", tc.name, err)
			if err == nil {
				t.Fatalf("case %q: mismatched config accepted", tc.name)
			}
			if !IsConfigError(err) || IsStabilityError(err) {
				t.Fatalf("case %q: expected pure ErrKindConfig, got %v", tc.name, err)
			}
		})
	}
	// 未修改的配置必须通过
	if err := cp.VerifyConfig(phaseConfig()); err != nil {
		t.Fatalf("matching config rejected: %v", err)
	}
}

// TestLegacyRunPhasesCompatible 验证无阶段计划时 RunPhases 与旧入口逐位一致。
func TestLegacyRunPhasesCompatible(t *testing.T) {
	cfg := baseConfig() // 无 Plan
	cfg.Steps = 50
	initial := testInitial(cfg.NX, cfg.NY)

	legacy, err := RunParallel(cfg, initial, 3, cfg.Steps)
	if err != nil {
		t.Fatalf("RunParallel: %v", err)
	}
	got, reports, err := RunPhases(cfg, initial, 3, 0, cfg.Steps)
	if err != nil {
		t.Fatalf("RunPhases: %v", err)
	}
	diff := legacy.MaxAbsDiff(got)
	t.Logf("legacy config: RunPhases vs RunParallel maxAbsDiff=%.3e tolerance=0", diff)
	if diff != 0 {
		t.Fatalf("legacy behavior changed: diff %g", diff)
	}
	if len(reports) != 1 || reports[0].PhaseIndex != -1 {
		t.Fatalf("legacy run should yield one report with PhaseIndex=-1, got %+v", reports)
	}
	if !reports[0].Balanced() {
		t.Fatalf("legacy report does not close: residual %g tol %g",
			reports[0].Closure(), reports[0].ClosureTol())
	}
	t.Logf("legacy report: steps[%d,%d) closure=%.3e tol=%.3e",
		reports[0].FromStep, reports[0].ToStep, reports[0].Closure(), reports[0].ClosureTol())
}
