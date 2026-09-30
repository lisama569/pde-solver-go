package solver

import (
	"encoding/json"
	"math"
	"path/filepath"
	"strconv"
	"testing"
)

// --- 测试夹具 ---------------------------------------------------------------

func dirBC(v float64) Boundary { return Boundary{Type: Dirichlet, Value: v} }

var neuBC = Boundary{Type: Neumann}

// heatCoolPlan：0-20 中心加热；20-40 关源冷却（边界不变）；40-60 左墙降温+
// 顶墙改成固定值，并换成一正一负两个源区。
func heatCoolPlan() *PhasePlan {
	return &PhasePlan{Phases: []Phase{
		{StartStep: 0, EndStep: 20,
			BC:      Boundaries{Left: dirBC(100), Right: neuBC, Top: neuBC, Bottom: dirBC(0)},
			Sources: []PhaseSource{{Name: "heater", Rate: 50, X0: .4, X1: .6, Y0: .4, Y1: .6}}},
		{StartStep: 20, EndStep: 40,
			BC:      Boundaries{Left: dirBC(100), Right: neuBC, Top: neuBC, Bottom: dirBC(0)},
			Sources: nil},
		{StartStep: 40, EndStep: 60,
			BC: Boundaries{Left: dirBC(20), Right: neuBC, Top: dirBC(10), Bottom: dirBC(0)},
			Sources: []PhaseSource{
				{Name: "heater", Rate: 10, X0: .2, X1: .35, Y0: .2, Y1: .35},
				{Name: "sink", Rate: -7, X0: .7, X1: .8, Y0: .7, Y1: .8},
			}},
	}}
}

// allDirichletPlan：四边全部 Dirichlet + 多源（含角点覆盖与负强度），
// 用于验证角点附近通量不漏项、阶段切换时四边同时改值。
func allDirichletPlan(steps int) *PhasePlan {
	return &PhasePlan{Phases: []Phase{
		{StartStep: 0, EndStep: steps / 2,
			BC: Boundaries{Left: dirBC(80), Right: dirBC(-10), Top: dirBC(30), Bottom: dirBC(5)},
			Sources: []PhaseSource{
				{Name: "a", Rate: 40, X0: 0, X1: .3, Y0: 0, Y1: .3},
				{Name: "b", Rate: -20, X0: .7, X1: 1, Y0: .7, Y1: 1},
			}},
		{StartStep: steps / 2, EndStep: steps,
			BC:      Boundaries{Left: dirBC(5), Right: dirBC(40), Top: dirBC(-5), Bottom: dirBC(60)},
			Sources: []PhaseSource{{Name: "a", Rate: 5, X0: .1, X1: .9, Y0: .1, Y1: .9}}},
	}}
}

func phaseTestConfig() *Config {
	c := baseConfig()
	c.Steps = 60
	return c
}

// logReport 打印一份阶段报告的输入参数、容差、判定依据与实测偏差。
func logReport(t *testing.T, rep PhaseReport) {
	t.Helper()
	t.Logf("phase %d window [%d,%d) steps=%d: startMass=%.9f endMass=%.9f massChange=%+.6e",
		rep.PhaseIndex, rep.StartStep, rep.EndStep, rep.EndStep-rep.StartStep,
		rep.StartMass, rep.EndMass, rep.MassChange)
	t.Logf("  boundaryReset=%+.6e", rep.BoundaryReset)
	inj := 0.0
	for _, s := range rep.Sources {
		t.Logf("  source[%s] injection=%+.9f", s.Name, s.Rate)
		inj += s.Rate
	}
	flow := 0.0
	for _, e := range rep.Edges {
		t.Logf("  edge[%s] netInflow=%+.9f", e.Edge, e.Flow)
		flow += e.Flow
	}
	t.Logf("  income(inj=%+.6e flow=%+.6e reset=%+.6e) vs massChange=%+.6e -> residual=%+.3e tolerance=%+.3e closed=%v",
		inj, flow, rep.BoundaryReset, rep.MassChange, rep.Residual, rep.Tolerance, rep.Closed)
}

// TestPhaseSwitchingAndBalance 验证多阶段切换正确、每阶段收支在容差内闭合。
func TestPhaseSwitchingAndBalance(t *testing.T) {
	cfg := phaseTestConfig()
	plan := heatCoolPlan()
	init := testInitial(cfg.NX, cfg.NY)
	t.Logf("input: grid=%dx%d dx=%g dy=%g dt=%g D=%g steps=%d dtLimit=%g",
		cfg.NX, cfg.NY, cfg.DX, cfg.DY, cfg.DT, cfg.Diffusivity, cfg.Steps, cfg.StabilityLimit())
	t.Logf("tolerance: abs=%g rel=%g (tol = abs + rel*scale)", DefaultBalanceTolAbs, DefaultBalanceTolRel)

	r, err := NewRunner(cfg, plan, init)
	if err != nil {
		t.Fatalf("NewRunner: %v", err)
	}
	reps, err := r.Advance(cfg.Steps, 1)
	if err != nil {
		t.Fatalf("advance: %v", err)
	}
	if len(reps) != 3 {
		t.Fatalf("got %d phase reports, want 3", len(reps))
	}
	for i, rep := range reps {
		if rep.PhaseIndex != i || rep.StartStep != i*20 || rep.EndStep != (i+1)*20 {
			t.Fatalf("phase %d window mismatch: [%d,%d)", i, rep.StartStep, rep.EndStep)
		}
		logReport(t, rep)
		if !rep.Closed {
			t.Fatalf("phase %d not closed: |residual|=%g exceeds tolerance=%g",
				i, math.Abs(rep.Residual), rep.Tolerance)
		}
		if math.Abs(rep.Residual) > 1e-11 {
			t.Fatalf("phase %d residual %g far above rounding level", i, rep.Residual)
		}
	}
	// 阶段 2 的两个源必须分别可对账（多源区 + 负强度）。
	if len(reps[2].Sources) != 2 {
		t.Fatalf("phase 2 sources = %d, want 2", len(reps[2].Sources))
	}
	if reps[2].Sources[0].Name != "heater" || reps[2].Sources[0].Rate <= 0 {
		t.Fatalf("heater injection expected positive, got %+v", reps[2].Sources[0])
	}
	if reps[2].Sources[1].Name != "sink" || reps[2].Sources[1].Rate >= 0 {
		t.Fatalf("sink injection expected negative, got %+v", reps[2].Sources[1])
	}
	// 阶段 1 无源：注入合计必须为 0，质量变化仅由边界通量解释。
	if len(reps[1].Sources) != 0 {
		t.Fatalf("phase 1 should have no sources, got %d", len(reps[1].Sources))
	}
	// 阶段 2 左墙 100->20、顶墙改为 Dirichlet(10)，应产生负向重置跳变。
	if reps[2].BoundaryReset >= 0 {
		t.Fatalf("phase 2 boundary reset should be negative, got %g", reps[2].BoundaryReset)
	}
}

// TestPhaseBalanceAllDirichlet 四边 Dirichlet + 角点源区，验证角点附近通量
// 不漏项（若求和行/列范围写错，这里会明显不闭合）。
func TestPhaseBalanceAllDirichlet(t *testing.T) {
	cfg := phaseTestConfig()
	plan := allDirichletPlan(60)
	r, err := NewRunner(cfg, plan, testInitial(cfg.NX, cfg.NY))
	if err != nil {
		t.Fatal(err)
	}
	reps, err := r.Advance(60, 1)
	if err != nil {
		t.Fatal(err)
	}
	for _, rep := range reps {
		logReport(t, rep)
		if !rep.Closed {
			t.Fatalf("all-dirichlet phase %d not closed: residual=%g tol=%g",
				rep.PhaseIndex, rep.Residual, rep.Tolerance)
		}
	}
}

// TestInvalidPhasePlansRejected 非法阶段计划必须在推进前按配置问题拒绝，
// 且归类为 ErrKindConfig（而不是失稳）。
func TestInvalidPhasePlansRejected(t *testing.T) {
	validBC := Boundaries{Left: neuBC, Right: neuBC, Top: neuBC, Bottom: neuBC}
	ph := func(start, end int, bc Boundaries, srcs ...PhaseSource) Phase {
		return Phase{StartStep: start, EndStep: end, BC: bc, Sources: srcs}
	}
	plan1 := func(p Phase) *PhasePlan { return &PhasePlan{Phases: []Phase{p}} }
	bad := []struct {
		name string
		plan *PhasePlan
	}{
		{"no phases", &PhasePlan{}},
		{"first phase not from zero", plan1(ph(5, 20, validBC))},
		{"gap between phases", &PhasePlan{Phases: []Phase{ph(0, 10, validBC), ph(12, 20, validBC)}}},
		{"overlapping phases", &PhasePlan{Phases: []Phase{ph(0, 10, validBC), ph(8, 20, validBC)}}},
		{"empty phase range", plan1(ph(0, 0, validBC))},
		{"unknown boundary type", plan1(ph(0, 10, Boundaries{Left: Boundary{Type: BoundaryType(7)}, Right: neuBC, Top: neuBC, Bottom: neuBC}))},
		{"unpaired periodic", plan1(ph(0, 10, Boundaries{Left: Boundary{Type: Periodic}, Right: neuBC, Top: neuBC, Bottom: neuBC}))},
		{"invalid source region", plan1(ph(0, 10, validBC, PhaseSource{Rate: 1, X0: .8, X1: .2, Y0: 0, Y1: 1}))},
		{"nan source rate", plan1(ph(0, 10, validBC, PhaseSource{Rate: math.NaN(), X0: 0, X1: 1, Y0: 0, Y1: 1}))},
		{"inf dirichlet value", plan1(ph(0, 10, Boundaries{Left: dirBC(math.Inf(1)), Right: neuBC, Top: neuBC, Bottom: neuBC}))},
		{"duplicate source name in phase", plan1(ph(0, 10, validBC,
			PhaseSource{Name: "x", Rate: 1, X0: 0, X1: .5, Y0: 0, Y1: 1},
			PhaseSource{Name: "x", Rate: 1, X0: .5, X1: 1, Y0: 0, Y1: 1}))},
	}
	for _, tc := range bad {
		t.Run(tc.name, func(t *testing.T) {
			err := ValidatePlan(tc.plan)
			t.Logf("plan=%+v err=%v", tc.plan, err)
			if err == nil {
				t.Fatalf("expected config rejection for %q", tc.name)
			}
			if !IsConfigError(err) {
				t.Fatalf("%q: want ErrKindConfig, got %v", tc.name, err)
			}
			if IsStabilityError(err) {
				t.Fatalf("%q: config error misclassified as stability", tc.name)
			}
		})
	}

	// 计划本身合法但盖不住 cfg.Steps：同样在推进前拒绝。
	cfg := phaseTestConfig()
	short := plan1(ph(0, 10, validBC))
	if err := ValidatePlan(short); err != nil {
		t.Fatalf("short plan itself is structurally valid: %v", err)
	}
	_, err := NewRunner(cfg, short, testInitial(cfg.NX, cfg.NY))
	if err == nil || !IsConfigError(err) {
		t.Fatalf("plan shorter than Steps: want ErrKindConfig, got %v", err)
	}
	t.Logf("coverage rejection: %v", err)
}

// TestPhaseStabilityStillDetected 合法计划 + 非法 dt 必须报 ErrKindStability，
// 与配置错误明确区分；绕过前置校验时运行期 NaN 防护仍归类为失稳。
func TestPhaseStabilityStillDetected(t *testing.T) {
	cfg := phaseTestConfig()
	cfg.DT = cfg.StabilityLimit() * 1.5
	_, err := NewRunner(cfg, heatCoolPlan(), testInitial(cfg.NX, cfg.NY))
	if err == nil || !IsStabilityError(err) {
		t.Fatalf("oversized dt: want ErrKindStability, got %v", err)
	}
	t.Logf("pre-run stability rejection: %v", err)

	// 合法 dt 启动，但在热源极强 + 不稳定 dt 下运行期出现非有限值（绕过前置校验）。
	cfg2 := phaseTestConfig()
	cfg2.DT = cfg2.StabilityLimit() * 10
	plan := heatCoolPlan()
	st := &Stepper{Cfg: cfg2, Cur: testInitial(cfg2.NX, cfg2.NY), Next: NewGrid(cfg2.NX, cfg2.NY)}
	st.SetActive(&plan.Phases[0].BC, plan.Phases[0].Sources)
	err = st.Advance(500)
	if err == nil || !IsStabilityError(err) {
		t.Fatalf("runtime blowup: want ErrKindStability, got %v", err)
	}
	t.Logf("runtime stability rejection: %v", err)
}

// TestPhaseWorkersConsistency 不同并行子域数与串行基线逐位一致（容差 0）。
func TestPhaseWorkersConsistency(t *testing.T) {
	for pi, plan := range []*PhasePlan{heatCoolPlan(), allDirichletPlan(60)} {
		cfg := phaseTestConfig()
		base, err := NewRunner(cfg, plan, testInitial(cfg.NX, cfg.NY))
		if err != nil {
			t.Fatal(err)
		}
		if _, err := base.Advance(60, 1); err != nil {
			t.Fatal(err)
		}
		for _, w := range []int{1, 2, 3, 4, 7, 16} {
			r, err := NewRunner(phaseTestConfig(), plan, testInitial(cfg.NX, cfg.NY))
			if err != nil {
				t.Fatal(err)
			}
			if _, err := r.Advance(60, w); err != nil {
				t.Fatalf("plan %d workers=%d: %v", pi, w, err)
			}
			d := base.Grid().MaxAbsDiff(r.Grid())
			t.Logf("plan=%d workers=%d maxAbsDiff=%.3e tolerance=0", pi, w, d)
			if d != 0 {
				t.Fatalf("plan %d workers=%d deviates from serial: %g", pi, w, d)
			}
		}
	}
}

// runReference 一次性串行跑到底，返回终场；用于恢复一致性对照。
func runReference(t *testing.T, cfg *Config, plan *PhasePlan) *Grid {
	t.Helper()
	r, err := NewRunner(cfg, plan, testInitial(cfg.NX, cfg.NY))
	if err != nil {
		t.Fatal(err)
	}
	if _, err := r.Advance(cfg.Steps, 1); err != nil {
		t.Fatal(err)
	}
	return r.Grid()
}

// saveAt 在 split 步落盘状态：先用 w1 个子域跑到 split。
func saveAt(t *testing.T, cfg *Config, plan *PhasePlan, split, w1 int, path string) {
	t.Helper()
	r, err := NewRunner(cfg, plan, testInitial(cfg.NX, cfg.NY))
	if err != nil {
		t.Fatal(err)
	}
	if _, err := r.Advance(split, w1); err != nil {
		t.Fatal(err)
	}
	if err := r.Save(path); err != nil {
		t.Fatal(err)
	}
}

// TestPhaseSaveResume 覆盖阶段内部（15、35）与刚好跨阶段边界（20、40）
// 两类保存点；恢复后切换子域数（5），终场必须与不间断串行基线逐位一致。
func TestPhaseSaveResume(t *testing.T) {
	cfg := phaseTestConfig()
	plan := heatCoolPlan()
	reference := runReference(t, cfg, plan)

	for _, split := range []int{15, 20, 35, 40} {
		t.Run(describeSplit(split), func(t *testing.T) {
			path := filepath.Join(t.TempDir(), "state.bin")
			saveAt(t, cfg, plan, split, 3, path)
			st, err := LoadState(path)
			if err != nil {
				t.Fatal(err)
			}
			t.Logf("saved at split=%d: state step=%d grid=%dx%d hasPlan=%v planHash=%016x",
				split, st.Step, st.NX, st.NY, st.Sig.HasPlan, st.Sig.PlanHash)
			if st.Step != split {
				t.Fatalf("state step=%d want %d", st.Step, split)
			}
			r, err := ResumeRunner(cfg, plan, st)
			if err != nil {
				t.Fatalf("resume: %v", err)
			}
			reps, err := r.Advance(cfg.Steps-split, 5) // 恢复后改用 5 个子域
			if err != nil {
				t.Fatalf("advance after resume: %v", err)
			}
			for _, rep := range reps {
				logReport(t, rep)
				if !rep.Closed {
					t.Fatalf("post-resume phase %d not closed: residual=%g tol=%g",
						rep.PhaseIndex, rep.Residual, rep.Tolerance)
				}
			}
			d := reference.MaxAbsDiff(r.Grid())
			t.Logf("split=%d postResumeWorkers=5 finalStep=%d maxAbsDiff=%.3e tolerance=0",
				split, r.Step(), d)
			if d != 0 {
				t.Fatalf("resumed trajectory at split=%d differs: %g", split, d)
			}
			if r.Step() != cfg.Steps {
				t.Fatalf("final step=%d want %d", r.Step(), cfg.Steps)
			}
		})
	}
}

func describeSplit(split int) string {
	switch split {
	case 20, 40:
		return "boundary@" + strconv.Itoa(split)
	default:
		return "midphase@" + strconv.Itoa(split)
	}
}

// TestResumeRejectsChangedConfig 恢复时网格/步长/系数/阶段计划对不上，
// 必须在继续计算前以 ErrKindConfig 明确拒绝。
func TestResumeRejectsChangedConfig(t *testing.T) {
	cfg := phaseTestConfig()
	plan := heatCoolPlan()
	path := filepath.Join(t.TempDir(), "state.bin")
	saveAt(t, cfg, plan, 25, 1, path)
	st, err := LoadState(path)
	if err != nil {
		t.Fatal(err)
	}

	mutate := []struct {
		name string
		cfg  *Config
		plan *PhasePlan
	}{
		{"grid NX changed", mutateConfig(cfg, func(c *Config) { c.NX++ }), clonePlan(plan)},
		{"DT changed", mutateConfig(cfg, func(c *Config) { c.DT *= 1.0001 }), clonePlan(plan)},
		{"DX changed", mutateConfig(cfg, func(c *Config) { c.DX = 0.02 }), clonePlan(plan)},
		{"diffusivity changed", mutateConfig(cfg, func(c *Config) { c.Diffusivity = 1.1 }), clonePlan(plan)},
		{"plan source strength changed", phaseTestConfig(), mutatePlan(plan, func(p *PhasePlan) {
			p.Phases[0].Sources[0].Rate = 99
		})},
		{"plan boundary changed", phaseTestConfig(), mutatePlan(plan, func(p *PhasePlan) {
			p.Phases[0].BC.Left.Value = 50
		})},
		{"plan step boundary changed", phaseTestConfig(), mutatePlan(plan, func(p *PhasePlan) {
			p.Phases[0].EndStep = 21
			p.Phases[1].StartStep = 21
		})},
	}
	for _, tc := range mutate {
		_, err := ResumeRunner(tc.cfg, tc.plan, st)
		t.Logf("%s -> %v", tc.name, err)
		if err == nil {
			t.Fatalf("%s: resume must be rejected", tc.name)
		}
		if !IsConfigError(err) {
			t.Fatalf("%s: want ErrKindConfig, got %v", tc.name, err)
		}
		if IsStabilityError(err) {
			t.Fatalf("%s: mismatch misclassified as stability", tc.name)
		}
	}

	// 保存时有阶段计划，恢复时给 nil 计划也要拒绝。
	if _, err := ResumeRunner(phaseTestConfig(), nil, st); err == nil || !IsConfigError(err) {
		t.Fatalf("nil plan resume: want config error, got %v", err)
	}
}

func mutateConfig(base *Config, f func(*Config)) *Config {
	c := *base
	f(&c)
	return &c
}

func clonePlan(p *PhasePlan) *PhasePlan {
	raw, _ := jsonMarshalPlan(p)
	out, _ := ParsePhasePlanJSON(raw)
	return out
}

func mutatePlan(p *PhasePlan, f func(*PhasePlan)) *PhasePlan {
	q := clonePlan(p)
	f(q)
	return q
}

// TestPhasePlanJSONRoundTrip 字符串/整数两种边界类型写法都能解析，
// 且指纹与原计划一致。
func TestPhasePlanJSONRoundTrip(t *testing.T) {
	plan := heatCoolPlan()
	raw, err := jsonMarshalPlan(plan)
	if err != nil {
		t.Fatal(err)
	}
	t.Logf("json: %s", raw)
	got, err := ParsePhasePlanJSON(raw)
	if err != nil {
		t.Fatal(err)
	}
	if err := ValidatePlan(got); err != nil {
		t.Fatal(err)
	}
	if got.fingerprint() != plan.fingerprint() {
		t.Fatalf("plan fingerprint changed after json round trip: %016x vs %016x",
			got.fingerprint(), plan.fingerprint())
	}
	// 整数写法（兼容旧 BoundaryType 枚举）也必须能解析。
	intForm := []byte(`{"phases":[{"start_step":0,"end_step":10,` +
		`"bc":{"left":{"type":1,"value":0},"right":{"type":1,"value":0},` +
		`"top":{"type":1,"value":0},"bottom":{"type":1,"value":0}},"sources":[]}]}`)
	if _, err := ParsePhasePlanJSON(intForm); err != nil {
		t.Fatalf("integer boundary form rejected: %v", err)
	}
}

func jsonMarshalPlan(p *PhasePlan) ([]byte, error) { return json.Marshal(p) }

// TestPhasePlanLegacyEquivalence 单阶段（边界=旧 BC、热源=旧 Src）的阶段计划
// 必须与旧 Stepper/RunParallel 计算路径逐位一致。
//
// 注意一个有意的语义差异：阶段模式在进入每个阶段（含阶段 0）时会把 Dirichlet
// 边「吸附」到声明固定值并计入 BoundaryReset；旧模式只冻结初始场边值。
// 因此这里用边初值已等于声明值的初始场，使吸附跳变为 0、两路径可逐位比较。
func TestPhasePlanLegacyEquivalence(t *testing.T) {
	cfg := phaseTestConfig()
	init := testInitial(cfg.NX, cfg.NY)
	// 让初始场的 Dirichlet 边已经等于声明值（顺序与 resetDirichlet 一致）。
	for j := 0; j < cfg.NY; j++ {
		init.Set(0, j, cfg.BC.Left.Value)
	}
	for i := 0; i < cfg.NX; i++ {
		init.Set(i, 0, cfg.BC.Bottom.Value)
	}

	// 旧路径
	old, err := NewStepper(cfg, init)
	if err != nil {
		t.Fatal(err)
	}
	if err := old.Advance(60); err != nil {
		t.Fatal(err)
	}

	// 阶段计划：整段一个阶段，使用完全相同的边界与源
	plan := &PhasePlan{Phases: []Phase{{
		StartStep: 0, EndStep: 60, BC: cfg.BC,
		Sources: []PhaseSource{{
			Name: "legacy-src", Rate: cfg.Src.Rate,
			X0: cfg.Src.X0, X1: cfg.Src.X1, Y0: cfg.Src.Y0, Y1: cfg.Src.Y1,
		}},
	}}}
	r, err := NewRunner(cfg, plan, init)
	if err != nil {
		t.Fatal(err)
	}
	reps, err := r.Advance(60, 4)
	if err != nil {
		t.Fatal(err)
	}
	if reps[0].BoundaryReset != 0 {
		t.Fatalf("with pre-set edges, phase-0 reset must be 0, got %g", reps[0].BoundaryReset)
	}
	d := old.Cur.MaxAbsDiff(r.Grid())
	t.Logf("single-phase plan vs legacy stepper: maxAbsDiff=%.3e tolerance=0", d)
	if d != 0 {
		t.Fatalf("single-phase plan changed legacy trajectory: diff=%g", d)
	}
}

// TestPhaseBalancePeriodic 周期边界下逐对相消：边通量应为 0，收支仍闭合。
func TestPhaseBalancePeriodic(t *testing.T) {
	cfg := phaseTestConfig()
	peri := Boundary{Type: Periodic}
	plan := &PhasePlan{Phases: []Phase{
		{StartStep: 0, EndStep: 25,
			BC:      Boundaries{Left: peri, Right: peri, Top: peri, Bottom: peri},
			Sources: []PhaseSource{{Name: "p", Rate: 33, X0: .3, X1: .7, Y0: .3, Y1: .7}}},
		{StartStep: 25, EndStep: 60,
			BC:      Boundaries{Left: peri, Right: peri, Top: peri, Bottom: peri},
			Sources: nil},
	}}
	r, err := NewRunner(cfg, plan, testInitial(cfg.NX, cfg.NY))
	if err != nil {
		t.Fatal(err)
	}
	reps, err := r.Advance(60, 2)
	if err != nil {
		t.Fatal(err)
	}
	for _, rep := range reps {
		logReport(t, rep)
		if !rep.Closed {
			t.Fatalf("periodic phase %d not closed: residual=%g tol=%g", rep.PhaseIndex, rep.Residual, rep.Tolerance)
		}
		for _, e := range rep.Edges {
			if e.Flow != 0 {
				t.Fatalf("periodic edge %s flow must be 0, got %g", e.Edge, e.Flow)
			}
		}
		if rep.BoundaryReset != 0 {
			t.Fatalf("periodic reset must be 0, got %g", rep.BoundaryReset)
		}
	}
}
