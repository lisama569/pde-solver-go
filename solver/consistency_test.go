package solver

import (
	"math"
	"path/filepath"
	"testing"
)

// testInitial 构造确定性的初始场（高斯热点），避免随机数。
func testInitial(nx, ny int) *Grid {
	g := NewGrid(nx, ny)
	for j := 0; j < ny; j++ {
		for i := 0; i < nx; i++ {
			x := float64(i)/float64(nx-1) - 0.5
			y := float64(j)/float64(ny-1) - 0.5
			g.Set(i, j, 100*math.Exp(-(x*x+y*y)/0.02))
		}
	}
	return g
}

// TestSerialParallelConsistency 验证不同并发度与串行基线逐位一致（容差 0）。
func TestSerialParallelConsistency(t *testing.T) {
	cfg := baseConfig()
	cfg.Steps = 50
	initial := testInitial(cfg.NX, cfg.NY)

	// 串行基线
	base, err := NewStepper(cfg, initial)
	if err != nil {
		t.Fatalf("baseline stepper: %v", err)
	}
	if err := base.Advance(cfg.Steps); err != nil {
		t.Fatalf("baseline advance: %v", err)
	}

	for _, workers := range []int{1, 2, 3, 4, 7, 16} {
		got, err := RunParallel(cfg, initial, workers, cfg.Steps)
		if err != nil {
			t.Fatalf("workers=%d: %v", workers, err)
		}
		diff := base.Cur.MaxAbsDiff(got)
		t.Logf("workers=%d steps=%d maxAbsDiff=%.3e tolerance=0 (bit-identical)", workers, cfg.Steps, diff)
		if diff != 0 {
			t.Fatalf("workers=%d: parallel result deviates from serial baseline by %g", workers, diff)
		}
	}
}

// TestResumeConsistency 验证从检查点恢复的轨迹与不间断运行逐位一致。
func TestResumeConsistency(t *testing.T) {
	cfg := baseConfig()
	cfg.Steps = 60
	initial := testInitial(cfg.NX, cfg.NY)

	// 不间断运行
	full, err := NewStepper(cfg, initial)
	if err != nil {
		t.Fatalf("full run: %v", err)
	}
	if err := full.Advance(cfg.Steps); err != nil {
		t.Fatalf("full run advance: %v", err)
	}

	// 分段运行：25 步 -> 落盘 -> 恢复 -> 再跑 35 步
	const split = 25
	part, err := NewStepper(cfg, initial)
	if err != nil {
		t.Fatalf("partial run: %v", err)
	}
	if err := part.Advance(split); err != nil {
		t.Fatalf("partial advance: %v", err)
	}
	path := filepath.Join(t.TempDir(), "ckpt.bin")
	if err := SaveCheckpoint(path, part.Step, part.Cur); err != nil {
		t.Fatalf("save checkpoint: %v", err)
	}

	cp, err := LoadCheckpoint(path)
	if err != nil {
		t.Fatalf("load checkpoint: %v", err)
	}
	t.Logf("checkpoint: step=%d grid=%dx%d (saved at step %d of %d)", cp.Step, cp.NX, cp.NY, split, cfg.Steps)
	if cp.Step != split {
		t.Fatalf("checkpoint step = %d, want %d", cp.Step, split)
	}

	resumed, err := NewStepper(cfg, cp.Grid())
	if err != nil {
		t.Fatalf("resume stepper: %v", err)
	}
	resumed.Step = cp.Step
	if err := resumed.Advance(cfg.Steps - split); err != nil {
		t.Fatalf("resumed advance: %v", err)
	}

	diff := full.Cur.MaxAbsDiff(resumed.Cur)
	t.Logf("resume: steps=%d split@%d maxAbsDiff=%.3e tolerance=0 (bit-identical)", cfg.Steps, split, diff)
	if diff != 0 {
		t.Fatalf("resumed trajectory deviates from uninterrupted run by %g", diff)
	}
	if resumed.Step != full.Step {
		t.Fatalf("resumed step count = %d, want %d", resumed.Step, full.Step)
	}
}

// TestInvariantConservation 验证全 Neumann 边界 + 无源时总量守恒。
func TestInvariantConservation(t *testing.T) {
	cfg := baseConfig()
	cfg.BC = Boundaries{
		Left:   Boundary{Type: Neumann},
		Right:  Boundary{Type: Neumann},
		Top:    Boundary{Type: Neumann},
		Bottom: Boundary{Type: Neumann},
	}
	cfg.Src = nil
	cfg.Steps = 40
	initial := testInitial(cfg.NX, cfg.NY)

	m0 := TotalMass(initial, cfg.DX, cfg.DY)
	s, err := NewStepper(cfg, initial)
	if err != nil {
		t.Fatalf("stepper: %v", err)
	}
	maxDev := 0.0
	for k := 0; k < cfg.Steps; k++ {
		if err := s.Advance(1); err != nil {
			t.Fatalf("advance: %v", err)
		}
		rep := CheckInvariant(s.Cur, cfg, m0, s.Step)
		if math.Abs(rep.Deviation) > maxDev {
			maxDev = math.Abs(rep.Deviation)
		}
	}
	const tol = 1e-9
	t.Logf("Neumann+no-source: steps=%d initialMass=%.6f maxDeviation=%.3e tolerance=%.1e",
		cfg.Steps, m0, maxDev, tol)
	if maxDev > tol {
		t.Fatalf("mass not conserved: max deviation %g exceeds tolerance %g", maxDev, tol)
	}
}

// TestNaNGuard 验证绕过前置校验时 NaN 会被运行期防护捕获并归类为失稳。
func TestNaNGuard(t *testing.T) {
	cfg := baseConfig()
	cfg.DT = cfg.StabilityLimit() * 10 // 远超稳定极限
	if err := cfg.CheckStability(); err == nil {
		t.Fatal("expected stability rejection before run")
	}
	// 强行构造步进器（绕过 CheckStability）以触发运行期防护
	s := &Stepper{Cfg: cfg, Cur: testInitial(cfg.NX, cfg.NY), Next: NewGrid(cfg.NX, cfg.NY)}
	err := s.Advance(500)
	t.Logf("DT=10x limit: advance err=%v", err)
	if err == nil {
		t.Fatal("expected non-finite guard to fire")
	}
	if !IsStabilityError(err) {
		t.Fatalf("NaN guard error must be ErrKindStability, got %v", err)
	}
}
