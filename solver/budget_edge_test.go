package solver

import (
	"math"
	"os"
	"path/filepath"
	"strings"
	"testing"
)

// edgeConfig 返回贴近真实实验的两阶段计划（总 60 步）：
//
//	阶段0 [0,30)：左边界 Dirichlet 100、底边界 Dirichlet 0，三个热源区——
//	              src0 贴左边，src1 贴左下角（与 src0 部分重叠），src2 贴右上两边
//	阶段1 [30,60)：左边界改固定值 100→250，底边界改绝热（Dirichlet→Neumann），
//	              只保留一个贴左边的热源
func edgeConfig() *Config {
	cfg := baseConfig()
	cfg.Steps = 60
	cfg.Plan = &PhasePlan{Phases: []Phase{
		{
			From: 0, To: 30,
			BC: Boundaries{
				Left:   Boundary{Type: Dirichlet, Value: 100},
				Right:  Boundary{Type: Neumann},
				Top:    Boundary{Type: Neumann},
				Bottom: Boundary{Type: Dirichlet, Value: 0},
			},
			Sources: []Source{
				{Rate: 50, X0: 0.0, X1: 0.3, Y0: 0.3, Y1: 0.7}, // 贴左边
				{Rate: 40, X0: 0.0, X1: 0.2, Y0: 0.0, Y1: 0.5}, // 贴左下角，与 src0 重叠
				{Rate: 60, X0: 0.6, X1: 1.0, Y0: 0.8, Y1: 1.0}, // 贴右边+上边（均绝热）
			},
		},
		{
			From: 30, To: 60,
			BC: Boundaries{
				Left:   Boundary{Type: Dirichlet, Value: 250}, // 改固定值
				Right:  Boundary{Type: Neumann},
				Top:    Boundary{Type: Neumann},
				Bottom: Boundary{Type: Neumann}, // 固定改绝热
			},
			Sources: []Source{
				{Rate: 25, X0: 0.0, X1: 0.25, Y0: 0.4, Y1: 0.6}, // 贴左边
			},
		},
	}}
	return cfg
}

// fixedCell 与 Stepper.fixed 同一判定：(i,j) 是否为 Dirichlet 固定格。
func fixedCell(nx, ny int, bc Boundaries, i, j int) bool {
	if i == 0 && bc.Left.Type == Dirichlet {
		return true
	}
	if i == nx-1 && bc.Right.Type == Dirichlet {
		return true
	}
	if j == 0 && bc.Bottom.Type == Dirichlet {
		return true
	}
	if j == ny-1 && bc.Top.Type == Dirichlet {
		return true
	}
	return false
}

// countActiveSourceCells 统计源区域内**实际被更新**（非固定）的格点数。
func countActiveSourceCells(nx, ny int, bc Boundaries, s *Source) int {
	n := 0
	for j := 0; j < ny; j++ {
		y := float64(j) / float64(ny-1)
		for i := 0; i < nx; i++ {
			x := float64(i) / float64(nx-1)
			if x >= s.X0 && x <= s.X1 && y >= s.Y0 && y <= s.Y1 && !fixedCell(nx, ny, bc, i, j) {
				n++
			}
		}
	}
	return n
}

// logReport 打印阶段报告的输入参数、容差、判定依据与实测偏差。
func logReport(t *testing.T, r *PhaseReport) {
	t.Helper()
	t.Logf("phase %d steps[%d,%d): mass %.6f -> %.6f (delta %+.6e), srcTotal %+.6e, "+
		"boundary{L:%+.4e R:%+.4e T:%+.4e B:%+.4e}, reset{L:%+.4e R:%+.4e T:%+.4e B:%+.4e}, "+
		"closure %+.3e tol %.3e balanced=%v",
		r.PhaseIndex, r.FromStep, r.ToStep, r.MassStart, r.MassEnd, r.DeltaMass(),
		r.SourceTotal(), r.Boundary.Left, r.Boundary.Right, r.Boundary.Top, r.Boundary.Bottom,
		r.Reset.Left, r.Reset.Right, r.Reset.Top, r.Reset.Bottom,
		r.Closure(), r.ClosureTol(), r.Balanced())
}

// TestEdgeSourceBudget 验证贴边/贴角热源按实际施加量记账，阶段收支闭合，
// 且逐源注入量与"非固定格解析值"一致。
func TestEdgeSourceBudget(t *testing.T) {
	cfg := edgeConfig()
	t.Logf("input: grid=%dx%d dx=%g dy=%g dt=%g D=%g steps=%d phases=%d",
		cfg.NX, cfg.NY, cfg.DX, cfg.DY, cfg.DT, cfg.Diffusivity, cfg.Steps, len(cfg.Plan.Phases))
	_, reports, err := RunPhases(cfg, testInitial(cfg.NX, cfg.NY), 1, 0, cfg.Steps)
	if err != nil {
		t.Fatalf("RunPhases: %v", err)
	}
	if len(reports) != 2 {
		t.Fatalf("got %d reports, want 2", len(reports))
	}
	for k := range reports {
		r := &reports[k]
		logReport(t, r)
		if !r.Balanced() {
			t.Fatalf("phase %d budget does not close: residual %g exceeds tol %g",
				r.PhaseIndex, r.Closure(), r.ClosureTol())
		}
	}
	// 逐源注入量 == rate × 非固定格数 × dx·dy·dt × 步数
	anyExcluded := false
	for k := range reports {
		r := &reports[k]
		ph := &cfg.Plan.Phases[k]
		nSteps := float64(r.ToStep - r.FromStep)
		for si, src := range ph.Sources {
			total := countSourceCells(cfg.NX, cfg.NY, &src)
			active := countActiveSourceCells(cfg.NX, cfg.NY, ph.BC, &src)
			if total > active {
				anyExcluded = true
			}
			want := src.Rate * float64(active) * cfg.DX * cfg.DY * cfg.DT * nSteps
			dev := math.Abs(r.SourceInject[si] - want)
			t.Logf("phase %d source %d: region cells=%d fixed-excluded=%d active=%d "+
				"injected=%.6e analytic=%.6e dev=%.3e tol=%.3e",
				k, si, total, total-active, active, r.SourceInject[si], want, dev, 1e-9*math.Abs(want))
			if dev > 1e-9*math.Abs(want) {
				t.Fatalf("phase %d source %d injection %g deviates from analytic %g by %g",
					k, si, r.SourceInject[si], want, dev)
			}
		}
	}
	if !anyExcluded {
		t.Fatal("test intends edge-touching sources on Dirichlet edges, but no fixed cell was excluded")
	}
}

// TestPhaseSwitchResetAccounting 验证跨阶段修改固定边值 / 固定改绝热时：
// 重置量逐边可见、各阶段闭合、各阶段 ΔMass 之和等于任务真实总量变化。
func TestPhaseSwitchResetAccounting(t *testing.T) {
	cfg := edgeConfig()
	initial := testInitial(cfg.NX, cfg.NY)
	final, reports, err := RunPhases(cfg, initial, 1, 0, cfg.Steps)
	if err != nil {
		t.Fatalf("RunPhases: %v", err)
	}
	for k := range reports {
		logReport(t, &reports[k])
	}

	// 阶段1 把左边界 100→250。解析值用阶段0 结束时的真实场逐格计算：
	// 大部分左列格为 100，但左下角 (0,0) 在阶段0 被底边 Dirichlet 0 覆盖，
	// 因此该格跳变是 0→250 而非 100→250。
	at30, _, err := RunPhases(cfg, initial, 1, 0, 30)
	if err != nil {
		t.Fatalf("phase-0-only run: %v", err)
	}
	wantReset := 0.0
	for j := 0; j < cfg.NY; j++ {
		wantReset += (250.0 - at30.At(0, j)) * cfg.DX * cfg.DY
	}
	gotReset := reports[1].Reset.Left
	t.Logf("phase 1 left reset: got %+.6e analytic %+.6e dev %.3e (corner cell was %.3g, not 100)",
		gotReset, wantReset, math.Abs(gotReset-wantReset), at30.At(0, 0))
	if dev := math.Abs(gotReset - wantReset); dev > 1e-9*math.Abs(wantReset) {
		t.Fatalf("phase 1 left reset %g deviates from analytic %g by %g", gotReset, wantReset, dev)
	}
	// 底边界 固定→绝热：无写入，重置量必须为 0
	if reports[1].Reset.Bottom != 0 {
		t.Fatalf("phase 1 bottom reset = %g, want 0 (Dirichlet->Neumann writes nothing)",
			reports[1].Reset.Bottom)
	}
	t.Logf("phase 1 bottom reset=0 (Dirichlet->Neumann), as expected")

	// 任务级对账：各阶段 ΔMass 之和 == 真实总量变化（末报告 MassEnd − 首报告 MassStart）
	m0 := TotalMass(initial, cfg.DX, cfg.DY)
	m1 := TotalMass(final, cfg.DX, cfg.DY)
	sumDelta := TotalDeltaMass(reports)
	trueDelta := m1 - m0
	spanDelta := reports[len(reports)-1].MassEnd - reports[0].MassStart
	tol := 1e-9 * (1 + math.Abs(m0) + math.Abs(m1))
	t.Logf("task reconcile: sum(phase deltaMass)=%+.6e true(m1-m0)=%+.6e span=%+.6e "+
		"residual=%.3e tol=%.3e", sumDelta, trueDelta, spanDelta, math.Abs(sumDelta-trueDelta), tol)
	if math.Abs(sumDelta-trueDelta) > tol {
		t.Fatalf("sum of phase deltaMass %g != true change %g (residual %g > tol %g)",
			sumDelta, trueDelta, math.Abs(sumDelta-trueDelta), tol)
	}
	if math.Abs(spanDelta-trueDelta) > tol {
		t.Fatalf("report span %g != true change %g", spanDelta, trueDelta)
	}
	// 首报告 MassStart 必须是施加任何重置之前的原始初始场
	if reports[0].MassStart != m0 {
		t.Fatalf("first report MassStart %.9g != initial mass %.9g", reports[0].MassStart, m0)
	}
	if err := CheckBudgets(reports); err != nil {
		t.Fatalf("CheckBudgets: %v", err)
	}
}

// TestEdgeSourceCheckpointResume 验证贴边热源 + 跨阶段改边值场景下，
// 阶段内部与阶段边界两种位置的保存/恢复都与一次性跑到底逐位一致，
// 检查点确实落盘，恢复后改子域数仍与单子域基线一致。
func TestEdgeSourceCheckpointResume(t *testing.T) {
	cfg := edgeConfig()
	initial := testInitial(cfg.NX, cfg.NY)
	full, _, err := RunPhases(cfg, initial, 1, 0, cfg.Steps)
	if err != nil {
		t.Fatalf("uninterrupted run: %v", err)
	}
	for _, tc := range []struct {
		name    string
		split   int
		workers int
	}{
		{"inside phase 0", 12, 1},
		{"inside phase 1, parallel resume", 45, 4},
		{"at phase boundary", 30, 1},
		{"at phase boundary, parallel resume", 30, 3},
	} {
		t.Run(tc.name, func(t *testing.T) {
			dir := t.TempDir()
			path := filepath.Join(dir, "ckpt.bin")
			part, _, err := RunPhases(cfg, initial, 1, 0, tc.split)
			if err != nil {
				t.Fatalf("first segment: %v", err)
			}
			if err := SaveCheckpointWithConfig(path, tc.split, part, cfg); err != nil {
				t.Fatalf("save: %v", err)
			}
			if _, err := os.Stat(path); err != nil {
				t.Fatalf("checkpoint not on disk: %v", err)
			}
			cp, err := LoadCheckpoint(path)
			if err != nil {
				t.Fatalf("load: %v", err)
			}
			if err := cp.VerifyConfig(cfg); err != nil {
				t.Fatalf("verify: %v", err)
			}
			got, resumed, err := RunPhases(cfg, cp.Grid(), tc.workers, cp.Step, cfg.Steps-cp.Step)
			if err != nil {
				t.Fatalf("resumed segment: %v", err)
			}
			diff := full.MaxAbsDiff(got)
			t.Logf("split@%d workers=%d maxAbsDiff=%.3e tolerance=0 (bit-identical)",
				tc.split, tc.workers, diff)
			if diff != 0 {
				t.Fatalf("resumed run deviates from uninterrupted run by %g", diff)
			}
			// 续算段的报告自身也要闭合
			for i := range resumed {
				logReport(t, &resumed[i])
				if !resumed[i].Balanced() {
					t.Fatalf("resumed report %d does not close: residual %g",
						i, resumed[i].Closure())
				}
			}
		})
	}
}

// TestEdgeSourceParallelConsistency 验证贴边热源 + 阶段切换场景下
// 不同子域数结果与单子域基线逐位一致。
func TestEdgeSourceParallelConsistency(t *testing.T) {
	cfg := edgeConfig()
	initial := testInitial(cfg.NX, cfg.NY)
	base, baseReports, err := RunPhases(cfg, initial, 1, 0, cfg.Steps)
	if err != nil {
		t.Fatalf("baseline: %v", err)
	}
	for _, workers := range []int{2, 3, 5, 8} {
		got, reports, err := RunPhases(cfg, initial, workers, 0, cfg.Steps)
		if err != nil {
			t.Fatalf("workers=%d: %v", workers, err)
		}
		diff := base.MaxAbsDiff(got)
		t.Logf("workers=%d steps=%d maxAbsDiff=%.3e tolerance=0 (bit-identical)", workers, cfg.Steps, diff)
		if diff != 0 {
			t.Fatalf("workers=%d: result deviates from single-domain baseline by %g", workers, diff)
		}
		for k := range reports {
			if reports[k].Closure() != baseReports[k].Closure() {
				t.Fatalf("workers=%d phase %d: closure %g differs from baseline %g",
					workers, k, reports[k].Closure(), baseReports[k].Closure())
			}
		}
	}
}

// TestCheckBudgetsLocatesImbalance 验证对账真的不闭合时：
// 错误归类为 ErrKindBudget（与配置/失稳可区分），且信息定位到具体阶段。
func TestCheckBudgetsLocatesImbalance(t *testing.T) {
	balanced := PhaseReport{
		PhaseIndex: 0, FromStep: 0, ToStep: 10,
		SourceInject: []float64{5},
		Boundary:     BoundaryBudget{Left: 2},
		MassStart:    100, MassEnd: 107, // 7 = 5 + 2 闭合
	}
	imbalanced := PhaseReport{
		PhaseIndex: 1, FromStep: 10, ToStep: 20,
		SourceInject: []float64{5, 3},
		Boundary:     BoundaryBudget{Left: 2, Bottom: -1},
		MassStart:    107, MassEnd: 120, // 13 != 8+1 漏账
	}
	if !balanced.Balanced() {
		t.Fatal("fixture: balanced report should close")
	}
	if imbalanced.Balanced() {
		t.Fatal("fixture: imbalanced report should not close")
	}
	if err := CheckBudgets([]PhaseReport{balanced}); err != nil {
		t.Fatalf("balanced reports rejected: %v", err)
	}
	err := CheckBudgets([]PhaseReport{balanced, imbalanced})
	t.Logf("imbalanced reports: err=%v", err)
	if err == nil {
		t.Fatal("expected budget error")
	}
	if !IsBudgetError(err) {
		t.Fatalf("expected ErrKindBudget, got %v", err)
	}
	if IsConfigError(err) || IsStabilityError(err) {
		t.Fatalf("budget error misclassified: %v", err)
	}
	if !strings.Contains(err.Error(), "phase 1") {
		t.Fatalf("error should locate the offending phase, got: %v", err)
	}
	t.Logf("error kind=Budget, locates phase 1; distinguishable from config/stability errors")
}
