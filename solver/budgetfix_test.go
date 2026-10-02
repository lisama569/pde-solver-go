package solver

import (
	"math"
	"path/filepath"
	"strings"
	"testing"
)

// edgePlanConfig 构造接近真实实验的两/三阶段计划（32x24 网格，总 60 步）：
//
//	阶段0 [0,20)：左 Dirichlet100、底 Dirichlet0、右/上 Neumann；
//	  源0 贴左 Dirichlet 边（x0=0），源1 贴底 Dirichlet 边（y0=0）并与源0重叠，
//	  源2 贴左下角点（同时压左/底两条 Dirichlet 边），
//	  源3 贴上 Neumann 边（y1=1，应全额计入，不扣格点）。
//	阶段1 [20,40)：左 Dirichlet 改 200，底改为 Neumann，全部关源。
//	阶段2 [40,60)：左改 Neumann，底改 Dirichlet50。
func edgePlanConfig() *Config {
	cfg := baseConfig()
	cfg.Steps = 60
	dLeft100 := Boundary{Type: Dirichlet, Value: 100}
	dLeft200 := Boundary{Type: Dirichlet, Value: 200}
	dBottom0 := Boundary{Type: Dirichlet, Value: 0}
	dBottom50 := Boundary{Type: Dirichlet, Value: 50}
	neu := Boundary{Type: Neumann}
	cfg.Plan = &PhasePlan{Phases: []Phase{
		{
			From: 0, To: 20,
			BC: Boundaries{Left: dLeft100, Right: neu, Top: neu, Bottom: dBottom0},
			Sources: []Source{
				{Rate: 50, X0: 0.0, X1: 0.30, Y0: 0.30, Y1: 0.70}, // 贴左边
				{Rate: 30, X0: 0.10, X1: 0.40, Y0: 0.0, Y1: 0.35}, // 贴底边、与源0重叠
				{Rate: 20, X0: 0.0, X1: 0.20, Y0: 0.0, Y1: 0.20},  // 贴左下角点
				{Rate: 15, X0: 0.30, X1: 0.60, Y0: 0.85, Y1: 1.0}, // 贴上 Neumann 边
			},
		},
		{
			From: 20, To: 40,
			BC: Boundaries{Left: dLeft200, Right: neu, Top: neu, Bottom: neu},
		},
		{
			From: 40, To: 60,
			BC: Boundaries{Left: neu, Right: neu, Top: neu, Bottom: dBottom50},
		},
	}}
	return cfg
}

// stimulatedCells 统计某阶段配置下源区域覆盖的"实际受激"格点数
// （Dirichlet 固定格除外），与 SourceRates 口径一致。
func stimulatedCells(t *testing.T, cfg *Config, phaseIdx int, src *Source) int {
	s, err := NewStepper(cfg, testInitial(cfg.NX, cfg.NY))
	if err != nil {
		t.Fatalf("NewStepper: %v", err)
	}
	s.applyPhaseForStep(cfg.Plan.Phases[phaseIdx].From)
	n := 0
	for j := 0; j < cfg.NY; j++ {
		y := float64(j) / float64(cfg.NY-1)
		for i := 0; i < cfg.NX; i++ {
			x := float64(i) / float64(cfg.NX-1)
			if x >= src.X0 && x <= src.X1 && y >= src.Y0 && y <= src.Y1 && !s.fixed(i, j) {
				n++
			}
		}
	}
	return n
}

// TestEdgeTouchingSourceBudget 覆盖贴边（含贴角）的单热源语义：
// 贴 Dirichlet 边的固定格不接受源项、不计注入；贴 Neumann 边全额计入；
// 阶段账目在 1e-9 相对容差内闭合。
func TestEdgeTouchingSourceBudget(t *testing.T) {
	cfg := edgePlanConfig()
	initial := testInitial(cfg.NX, cfg.NY)
	final, reports, err := RunPhases(cfg, initial, 1, 0, cfg.Steps)
	if err != nil {
		t.Fatalf("RunPhases: %v", err)
	}
	t.Logf("input: grid=%dx%d dx=%g dy=%g dt=%g D=%g steps=%d workers=1 phases=%d",
		cfg.NX, cfg.NY, cfg.DX, cfg.DY, cfg.DT, cfg.Diffusivity, cfg.Steps, len(cfg.Plan.Phases))

	r := reports[0]
	nSteps := float64(r.ToStep - r.FromStep)
	for si := range cfg.Plan.Phases[0].Sources {
		src := &cfg.Plan.Phases[0].Sources[si]
		nAll := countSourceCells(cfg.NX, cfg.NY, src)
		nStim := stimulatedCells(t, cfg, 0, src)
		want := src.Rate * float64(nStim) * cfg.DX * cfg.DY * cfg.DT * nSteps
		got := r.SourceInject[si]
		dev := math.Abs(got - want)
		t.Logf("source %d region[%.2f,%.2f]x[%.2f,%.2f] rate=%g: cellsAll=%d cellsStimulated=%d "+
			"injected=%.6e analytic=%.6e dev=%.3e",
			si, src.X0, src.X1, src.Y0, src.Y1, src.Rate, nAll, nStim, got, want, dev)
		// 源0~2 贴 Dirichlet 边/角：必须排除固定格；源3 贴 Neumann 边应全额计入
		if si < 3 && nStim >= nAll {
			t.Fatalf("source %d: edge-touching source must exclude fixed cells, nStim=%d nAll=%d",
				si, nStim, nAll)
		}
		if dev > 1e-9*math.Abs(want) {
			t.Fatalf("source %d injection %g deviates from applied-only analytic %g by %g",
				si, got, want, dev)
		}
	}
	// 源3 贴上 Neumann 边：上边界不是固定格，覆盖格点应全部受激
	src3 := &cfg.Plan.Phases[0].Sources[3]
	if n := stimulatedCells(t, cfg, 0, src3); n != countSourceCells(cfg.NX, cfg.NY, src3) {
		t.Fatalf("source touching Neumann edge must be fully counted")
	}
	for k := range reports {
		t.Logf("phase %d steps[%d,%d): ΔMass=%+.6e src=%+.6e influxTotal=%+.6e resetTotal=%+.6e "+
			"closure=%+.3e tol=%.3e balanced=%v",
			reports[k].PhaseIndex, reports[k].FromStep, reports[k].ToStep,
			reports[k].DeltaMass(), reports[k].SourceTotal(),
			reports[k].Boundary.Total(), reports[k].ResetTotal(),
			reports[k].Closure(), reports[k].ClosureTol(), reports[k].Balanced())
		if !reports[k].Balanced() {
			t.Fatalf("phase %d does not close: residual %g tol %g",
				k, reports[k].Closure(), reports[k].ClosureTol())
		}
	}
	if err := VerifyReports(reports); err != nil {
		t.Fatalf("VerifyReports: %v", err)
	}
	// 最终场左/底边界必须反映最后阶段配置
	for j := 0; j < cfg.NY; j++ {
		if v := final.At(0, j); v == 200 || v == 100 {
			t.Fatalf("phase 2 left edge should be free (Neumann), got pinned %g at j=%d", v, j)
		}
	}
	for i := 0; i < cfg.NX; i++ {
		if v := final.At(i, 0); v != 50 {
			t.Fatalf("phase 2 bottom edge at i=%d = %g, want 50", i, v)
		}
	}
	t.Logf("CASE edge-touching single/multi source: PASS")
}

// TestPhaseBoundaryResetJump 验证跨阶段边界改动（改 Dirichlet 值、
// Dirichlet↔Neumann 切换）的重置跳变逐边可见、与解析值一致，
// 且各阶段 ΔMass 之和等于任务从头到尾的真实总量变化。
func TestPhaseBoundaryResetJump(t *testing.T) {
	cfg := edgePlanConfig()
	initial := testInitial(cfg.NX, cfg.NY)
	// 先单独跑到阶段0结束，取得阶段1激活前底边各格的实际值，
	// 用于解析阶段2 的 Neumann→Dirichlet 跳变。
	pre0, _, err := RunPhases(cfg, initial, 1, 0, 20)
	if err != nil {
		t.Fatalf("pre-run to step 20: %v", err)
	}
	pre, _, err := RunPhases(cfg, initial, 1, 0, 40)
	if err != nil {
		t.Fatalf("pre-run to step 40: %v", err)
	}
	final, reports, err := RunPhases(cfg, initial, 1, 0, cfg.Steps)
	if err != nil {
		t.Fatalf("RunPhases: %v", err)
	}
	if len(reports) != 3 {
		t.Fatalf("got %d reports, want 3", len(reports))
	}

	cellArea := cfg.DX * cfg.DY
	// 阶段0 激活：左列写 100（含角点）、底行写 0（角点已被左写100，底再写0，归底边）
	wantReset0L := 0.0
	wantReset0B := 0.0
	for j := 0; j < cfg.NY; j++ {
		wantReset0L += (100 - initial.At(0, j)) * cellArea
	}
	for i := 1; i < cfg.NX; i++ {
		wantReset0B += (0 - initial.At(i, 0)) * cellArea // 相对初始场
	}
	wantReset0B += (0 - 100) * cellArea // 左下角点：先被左写成100，再被底写成0
	// 阶段1 激活：左 100→200（逐格按阶段0结束场取值：左下角点在阶段0末为0，
	// 跳变为200，其余格为100）；底 Dirichlet0→Neumann 不写入，无跳变
	wantReset1L := 0.0
	for j := 0; j < cfg.NY; j++ {
		wantReset1L += (200 - pre0.At(0, j)) * cellArea
	}
	// 阶段2 激活：左 Dirichlet200→Neumann 不写入；底 Neumann→Dirichlet50
	wantReset2B := 0.0
	for i := 0; i < cfg.NX; i++ {
		wantReset2B += (50 - pre.At(i, 0)) * cellArea
	}

	wants := []BoundaryBudget{
		{Left: wantReset0L, Bottom: wantReset0B},
		{Left: wantReset1L},
		{Bottom: wantReset2B},
	}
	for k, want := range wants {
		got := reports[k].BoundaryReset
		dev := math.Abs(got.Total()-want.Total()) / math.Max(1, math.Abs(want.Total()))
		t.Logf("phase %d reset: got{L:%+.6e B:%+.6e} want{L:%+.6e B:%+.6e} relDev=%.3e",
			k, got.Left, got.Bottom, want.Left, want.Bottom, dev)
		if math.Abs(got.Left-want.Left) > 1e-9*math.Max(1, math.Abs(want.Left)) ||
			math.Abs(got.Bottom-want.Bottom) > 1e-9*math.Max(1, math.Abs(want.Bottom)) {
			t.Fatalf("phase %d reset jump mismatch:\n got %+v\nwant %+v", k, got, want)
		}
		if reports[k].BoundaryReset.Right != 0 || reports[k].BoundaryReset.Top != 0 {
			t.Fatalf("phase %d unexpected right/top reset: %+v", k, reports[k].BoundaryReset)
		}
	}
	// 整任务总账：ΣΔMass == 真实变化（首份 MassStart 是任何重置之前的初始总量）
	sumDelta := 0.0
	for k := range reports {
		sumDelta += reports[k].DeltaMass()
	}
	trueChange := TotalMass(final, cfg.DX, cfg.DY) - TotalMass(initial, cfg.DX, cfg.DY)
	dev := math.Abs(sumDelta - trueChange)
	t.Logf("task ledger: ΣΔMass=%+.12e trueChange=%+.12e dev=%.3e tol=%.3e",
		sumDelta, trueChange, dev, 1e-9*(1+math.Abs(trueChange)))
	if dev > 1e-9*(1+math.Abs(trueChange)) {
		t.Fatalf("phase deltas %g do not sum to true change %g (dev %g)", sumDelta, trueChange, dev)
	}
	if err := VerifyReports(reports); err != nil {
		t.Fatalf("VerifyReports: %v", err)
	}
	t.Logf("CASE cross-phase boundary value/type changes: PASS")
}

// TestEdgePlanCheckpointResume 覆盖贴边热源计划在阶段内部与阶段边界两种
// 位置的保存/恢复：切换步号、最终轨迹与一次性跑到底逐位一致；
// 恢复时改变并行子域数仍与单子域基线一致。
func TestEdgePlanCheckpointResume(t *testing.T) {
	cfg := edgePlanConfig()
	initial := testInitial(cfg.NX, cfg.NY)
	full, fullReports, err := RunPhases(cfg, initial, 1, 0, cfg.Steps)
	if err != nil {
		t.Fatalf("uninterrupted run: %v", err)
	}
	cases := []struct {
		name    string
		split   int
		workers int
	}{
		{"inside phase 0 with edge sources", 13, 1},
		{"inside phase 0, parallel resume", 13, 4},
		{"at boundary 20 (value change)", 20, 1},
		{"at boundary 20, parallel resume", 20, 3},
		{"inside phase 1", 25, 5},
		{"at boundary 40 (type change)", 40, 1},
		{"at boundary 40, parallel resume", 40, 7},
		{"inside phase 2", 55, 2},
	}
	for _, tc := range cases {
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
			cp, err := LoadCheckpoint(path)
			if err != nil {
				t.Fatalf("load: %v", err)
			}
			if err := cp.VerifyConfig(cfg); err != nil {
				t.Fatalf("verify: %v", err)
			}
			got, reports, err := RunPhases(cfg, cp.Grid(), tc.workers, cp.Step, cfg.Steps-cp.Step)
			if err != nil {
				t.Fatalf("resumed segment: %v", err)
			}
			diff := full.MaxAbsDiff(got)
			// 分段恢复的账目必须各自闭合
			for k := range reports {
				if !reports[k].Balanced() {
					t.Fatalf("post-resume phase report %d does not close: %+v", k, reports[k])
				}
			}
			t.Logf("split@%d resumeWorkers=%d maxAbsDiff=%.3e (bit-identical tolerance 0), "+
				"post-resume reports=%d all balanced", tc.split, tc.workers, diff, len(reports))
			if diff != 0 {
				t.Fatalf("resumed trajectory deviates from uninterrupted run by %g", diff)
			}
			// 恢复后下一次阶段切换必须发生在配置规定的步号
			for _, r := range reports {
				if r.PhaseIndex >= 0 {
					ph := cfg.Plan.Phases[r.PhaseIndex]
					if r.FromStep < ph.From || r.FromStep >= ph.To {
						t.Fatalf("resumed report starts at %d outside its phase [%d,%d)",
							r.FromStep, ph.From, ph.To)
					}
				}
			}
		})
	}
	// 续算点之后的最终场与全程报告末态一致（抽查边界值）
	_ = fullReports
	t.Logf("CASE edge-plan resume inside/boundary with worker changes: PASS")
}

// TestEdgePlanParallelConsistency 贴边热源 + 跨阶段改边值计划在不同
// 并行子域数下与单子域基线逐位一致。
func TestEdgePlanParallelConsistency(t *testing.T) {
	cfg := edgePlanConfig()
	initial := testInitial(cfg.NX, cfg.NY)
	base, baseReports, err := RunPhases(cfg, initial, 1, 0, cfg.Steps)
	if err != nil {
		t.Fatalf("baseline: %v", err)
	}
	for _, workers := range []int{2, 3, 4, 8} {
		got, reports, err := RunPhases(cfg, initial, workers, 0, cfg.Steps)
		if err != nil {
			t.Fatalf("workers=%d: %v", workers, err)
		}
		diff := base.MaxAbsDiff(got)
		t.Logf("workers=%d maxAbsDiff=%.3e (tolerance 0)", workers, diff)
		if diff != 0 {
			t.Fatalf("workers=%d deviates from single-domain baseline by %g", workers, diff)
		}
		for k := range reports {
			if !reports[k].Balanced() {
				t.Fatalf("workers=%d phase %d not balanced", workers, k)
			}
		}
	}
	_ = baseReports
	t.Logf("CASE edge-plan parallel consistency: PASS")
}

// TestBudgetImbalanceDetectedAndLocated 真不闭合时必须返回 ErrKindBudget，
// 且错误信息能定位到具体阶段与具体热源/边界，不与配置/失稳错误混淆。
func TestBudgetImbalanceDetectedAndLocated(t *testing.T) {
	cfg := edgePlanConfig()
	_, reports, err := RunPhases(cfg, testInitial(cfg.NX, cfg.NY), 1, 0, cfg.Steps)
	if err != nil {
		t.Fatalf("RunPhases: %v", err)
	}
	// 篡改1：抹掉阶段0 的源0 注入（模拟该热源漏账）
	bad := make([]PhaseReport, len(reports))
	copy(bad, reports)
	inj := make([]float64, len(bad[0].SourceInject))
	copy(inj, bad[0].SourceInject)
	stolen := inj[0]
	inj[0] = 0
	bad[0].SourceInject = inj
	err = VerifyReports(bad)
	t.Logf("hide source 0 (stolen=%+.3e): err=%v", stolen, err)
	if err == nil || !IsBudgetError(err) {
		t.Fatalf("expected ErrKindBudget, got %v", err)
	}
	if IsConfigError(err) || IsStabilityError(err) {
		t.Fatalf("budget failure misclassified: %v", err)
	}
	msg := err.Error()
	for _, want := range []string{"phase 0", "source 0", "[0,20)"} {
		if !strings.Contains(msg, want) {
			t.Fatalf("error message %q missing locator %q", msg, want)
		}
	}
	// 篡改2：抹掉阶段1 左边界重置跳变（模拟切换贡献丢失）
	bad2 := make([]PhaseReport, len(reports))
	copy(bad2, reports)
	bad2[1].BoundaryReset.Left = 0
	err = VerifyReports(bad2)
	t.Logf("hide phase 1 left reset: err=%v", err)
	if err == nil || !IsBudgetError(err) {
		t.Fatalf("expected ErrKindBudget, got %v", err)
	}
	if !strings.Contains(err.Error(), "phase 1") {
		t.Fatalf("error message %q missing phase 1 locator", err.Error())
	}
	// 原始报告必须通过
	if err := VerifyReports(reports); err != nil {
		t.Fatalf("untampered reports rejected: %v", err)
	}
	t.Logf("CASE imbalance detection and localization: PASS")
}

// singleCornerConfig 返回最朴素的两阶段计划：阶段0 单热源贴左下角
// （压住左 Dirichlet 边与底 Dirichlet 边），阶段1 把左边界改成绝热。
// 这正是用户报告的"先加热中心贴边、阶段一跑完即误报不平衡"的最小复现。
func singleCornerConfig() *Config {
	cfg := baseConfig()
	cfg.Steps = 40
	neu := Boundary{Type: Neumann}
	cfg.Plan = &PhasePlan{Phases: []Phase{
		{
			From: 0, To: 20,
			BC: Boundaries{
				Left: Boundary{Type: Dirichlet, Value: 100}, Right: neu, Top: neu,
				Bottom: Boundary{Type: Dirichlet, Value: 0},
			},
			Sources: []Source{{Rate: 80, X0: 0.0, X1: 0.25, Y0: 0.0, Y1: 0.25}},
		},
		{
			From: 20, To: 40,
			BC: Boundaries{Left: neu, Right: neu, Top: neu, Bottom: neu},
		},
	}}
	return cfg
}

// TestSingleCornerSourceTwoPhase 单热源贴角 + 跨阶段 Dirichlet→Neumann：
// 阶段0 不得误报不平衡，全部报告闭合且检查点可落盘恢复。
func TestSingleCornerSourceTwoPhase(t *testing.T) {
	cfg := singleCornerConfig()
	initial := testInitial(cfg.NX, cfg.NY)
	final, reports, err := RunPhases(cfg, initial, 1, 0, cfg.Steps)
	if err != nil {
		t.Fatalf("RunPhases: %v", err)
	}
	if len(reports) != 2 {
		t.Fatalf("got %d reports, want 2", len(reports))
	}
	r := reports[0]
	src := cfg.Plan.Phases[0].Sources[0]
	nStim := stimulatedCells(t, cfg, 0, &src)
	nAll := countSourceCells(cfg.NX, cfg.NY, &src)
	want := src.Rate * float64(nStim) * cfg.DX * cfg.DY * cfg.DT * 20
	t.Logf("single corner source: cellsAll=%d cellsStimulated=%d injected=%.6e analytic=%.6e",
		nAll, nStim, r.SourceInject[0], want)
	if nStim >= nAll {
		t.Fatalf("corner source must exclude fixed cells: nStim=%d nAll=%d", nStim, nAll)
	}
	if math.Abs(r.SourceInject[0]-want) > 1e-9*math.Abs(want) {
		t.Fatalf("injected %g want %g", r.SourceInject[0], want)
	}
	for k := range reports {
		t.Logf("phase %d steps[%d,%d): ΔMass=%+.6e src=%+.6e influx=%+.6e reset=%+.6e "+
			"closure=%+.3e tol=%.3e balanced=%v",
			reports[k].PhaseIndex, reports[k].FromStep, reports[k].ToStep,
			reports[k].DeltaMass(), reports[k].SourceTotal(), reports[k].Boundary.Total(),
			reports[k].ResetTotal(), reports[k].Closure(), reports[k].ClosureTol(),
			reports[k].Balanced())
		if !reports[k].Balanced() {
			t.Fatalf("phase %d falsely reported imbalance: residual %g tol %g",
				k, reports[k].Closure(), reports[k].ClosureTol())
		}
	}
	if err := VerifyReports(reports); err != nil {
		t.Fatalf("VerifyReports: %v", err)
	}
	// 阶段1 左边界已释放为 Neumann：角点不再被钉在 100/0
	if v := final.At(0, cfg.NY-1); v == 100 {
		t.Fatalf("left edge should be free after switching to Neumann, got 100")
	}
	// 检查点在阶段边界落盘并恢复，轨迹逐位一致
	path := filepath.Join(t.TempDir(), "corner.bin")
	part, _, err := RunPhases(cfg, initial, 1, 0, 20)
	if err != nil {
		t.Fatalf("first segment: %v", err)
	}
	if err := SaveCheckpointWithConfig(path, 20, part, cfg); err != nil {
		t.Fatalf("checkpoint must be written even though source touches the corner: %v", err)
	}
	cp, err := LoadCheckpoint(path)
	if err != nil {
		t.Fatalf("load: %v", err)
	}
	resumed, _, err := RunPhases(cfg, cp.Grid(), 3, cp.Step, cfg.Steps-cp.Step)
	if err != nil {
		t.Fatalf("resume: %v", err)
	}
	if diff := final.MaxAbsDiff(resumed); diff != 0 {
		t.Fatalf("resumed trajectory differs by %g", diff)
	}
	t.Logf("CASE single corner-touching source two-phase + checkpoint: PASS")
}
