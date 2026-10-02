# pde-solver-go

本地可运行的二维热传导（扩散）方程数值求解流程。使用内存网格模拟随时间演化的标量场，
支持串行/按子域并发求解、**按时间阶段切换源项与边界**、检查点落盘与断点续算、
守恒量与逐阶段收支校验，保证不同入口与不同并行度下结果一致、可解释、可复现。

## 数值格式

求解方程：`∂u/∂t = D·∇²u + S(x,y)`（`S` 为可选矩形常源项，可多个叠加）。

- 空间：五点中心差分（二阶精度），行优先一维数组存储。
- 时间：显式 FTCS（向前欧拉），双缓冲推进——每步从 `Cur` 读、写 `Next`，随后整体交换。
- 初始场：确定性高斯热点（无随机数），保证可复现。

## 边界条件

四边可独立配置（`solver.Boundaries`）：

| 类型 | 含义 | 离散处理 |
|---|---|---|
| `Dirichlet` | 固定值 | 边界格不更新，保持设定值 |
| `Neumann` | 零法向通量（绝热） | 鬼格值 = 边界格自身值 |
| `Periodic` | 周期 | 索引回绕；**必须左右、上下成对配置**，否则视为非法配置 |

源项 `Source` 以归一化坐标 `[X0,X1]×[Y0,Y1]` 指定矩形区域与强度 `Rate`。

## 时间阶段计划（Phase Plan）

一次任务可切分为若干**连续、不重叠**的阶段，每阶段独立配置四边边界与任意多个
矩形源区。Go 侧为 `solver.PhasePlan`：

```go
cfg.Plan = &solver.PhasePlan{Phases: []solver.Phase{
    {From: 0, To: 200,                       // 阶段0：加热（步号区间 [0,200)）
     BC:  ...,                               // 本阶段四边边界
     Sources: []solver.Source{{Rate: 50, X0: 0.4, X1: 0.6, Y0: 0.4, Y1: 0.6}}},
    {From: 200, To: 400, BC: 全Neumann},      // 阶段1：关源冷却
    {From: 400, To: 600, BC: 左边界改为200},   // 阶段2：改边界 + 两个新源区
}}
```

CLI 侧用 `-plan plan.json` 传入同构的 JSON（`phases[].from/to/bc/sources`，
边界 `type` 取 `dirichlet|neumann|periodic`）。

**合法性约束**（违反时在推进前拒绝，归类 `ErrKindConfig`）：

- 至少一个阶段；首阶段从 0 开始，末阶段结束于 `Steps`（恰好覆盖 `[0, Steps)`）；
- 相邻阶段首尾相接：无间隙、无重叠；每阶段长度 > 0；
- 每阶段的边界类型合法、周期边界成对、Dirichlet 值有限；源区域在 `[0,1]` 内
  非空、强度有限（NaN/Inf 一律拒绝）。

**切换语义**：

- 阶段切换只发生在**步边界**上：步号 `n ∈ [From, To)` 的更新（`uⁿ→uⁿ⁺¹`）使用该
  阶段的配置。切换是步号的确定性函数，因此断点续算与不间断运行的切换点完全一致。
- 阶段激活时，该阶段 `Dirichlet` 边的 `Value` 会**写入边界格**（幂等）；角点被两条
  Dirichlet 边共享时按 left→right→bottom→top 顺序后者覆盖。切换瞬间边界值的跳变
  属于重置：它**不是扩散通量，但计入收支账目**——按边记入该阶段报告的
  `Boundary`（并单列在 `Reset` 字段），见「逐阶段收支对账」。
- 阶段内多个源区可重叠，强度相加；源区可以贴边、贴角（见下）。

**与旧配置的兼容**：`Plan == nil` 时行为与旧版完全一致——单套 `BC`/`Src` 全程生效，
且旧版不主动写入 Dirichlet 值（边界取初始场的值）的语义保持不变。旧版 V1 检查点
（无配置指纹）仍可读取，恢复时校验网格维度。

## 稳定性条件

显式格式要求（`Config.StabilityLimit` / `CheckStability`）：

```
dt <= 1 / (2·D·(1/dx² + 1/dy²))
```

- 启动前：`CheckStability` 拒绝超限的 `dt`，返回 **ErrKindStability**。
- 运行中：每步提交后扫描全场，一旦出现 NaN/Inf 立即中止，报错含**步号与格点坐标**。
- 非法参数（网格过小、负系数、非法阶段计划等）返回 **ErrKindConfig**。
- 收支对账不闭合返回 **ErrKindBudget**（见「逐阶段收支对账」）。
- 三类错误可通过 `solver.IsConfigError` / `solver.IsStabilityError` /
  `solver.IsBudgetError` 区分；阶段计划问题永远归类为配置错误，与运行期失稳、
  对账失败严格区分。

## 并行约束

`RunParallel` / `RunPhases` 按行把网格均分为 `workers` 个子域并发推进：

- 双缓冲保证各子域只写自己的行区间，**无数据竞争、无重复更新**（`go test -race` 验证）。
- 每个格点的浮点运算序列与串行完全相同，因此结果与串行基线**逐位一致**（容差 0），
  与并发度、执行顺序无关；断点续算后改变子域数结果仍与单子域基线一致。

## 检查点与断点续算

- `SaveCheckpointWithConfig` 以 V2 二进制格式（magic + 步号 + 维度 + DX/DY/DT/D +
  **有效阶段计划指纹** + float64 数据）原子落盘；`SaveCheckpoint` 保留 V1 旧格式。
- 恢复时 `Checkpoint.VerifyConfig` 在**继续计算前**校验：网格维度、空间/时间步长、
  扩散系数、阶段计划（含边界与源项）指纹。任何不一致都返回 `ErrKindConfig` 并拒绝，
  绝不按新配置静默续算。
- 恢复点可以在**阶段内部**或**阶段边界**上，恢复后阶段切换步号与一次性跑到底相同，
  最终轨迹与不间断运行**逐位一致**。

## 逐阶段收支对账

`RunPhases` 为每个阶段输出 `PhaseReport`：各热源累计注入（逐源 `SourceInject`）、
各边界累计净流入（逐边 `Boundary`，正 = 流入）、阶段激活时的 Dirichlet 重置量
（逐边 `Reset`，**已含在 `Boundary` 中**）、阶段首末总量。对账恒等式：

```
ΔMass = Σ源注入 + Σ边界净流入        残差 Closure = ΔMass − 源 − 边界
```

**记账口径**（与更新算子严格对应，这是闭合的前提）：

- **边界通量**：Dirichlet 边计入固定格与内部格的扩散交换（角部内侧格也固定时跳过），
  Neumann 边恒为 0，Periodic 边两两相消恒为 0。
- **贴边/贴角热源**：Dirichlet 固定格不被更新，源强从未施加到这些格上，因此
  注入量只按**实际进入非固定格**的部分统计（`SourceRates` 与 `stepRows` 跳过同一批
  固定格）。源区贴 Neumann/Periodic 边时该区域格点照常更新，全额计入。
- **跨阶段改边值的重置量**：阶段激活写入 Dirichlet 值造成的总量跳变，按边计入
  新阶段报告的 `Boundary`，并单列在 `Reset` 中便于核对；角点被两条边先后覆盖时
  跳变按写入次序拆分到两条边，各边之和严格等于该格实际变化。固定改绝热
  （Dirichlet→Neumann）不写入任何值，重置量为 0。

**容差与判定**：

- 理论残差为 0，实测为浮点求和顺序误差（~1e-15 相对量级）。
- `ClosureTol = 1e-9 × (1 + |MassStart| + |MassEnd| + Σ|源| + Σ|边界|)`，
  比舍入误差高约 6 个数量级，又足以发现任何物理量级的漏账。
- `Balanced()` 判定闭合；`CheckBudgets(reports)` 把首个不闭合阶段报为
  **ErrKindBudget** 错误（与 `ErrKindConfig`/`ErrKindStability` 可区分），
  错误信息含阶段索引、步区间、残差与容差；具体是哪个热源或哪条边对不上，
  查该报告的逐源 `SourceInject` 与逐边 `Boundary`/`Reset` 分项。
- **任务级对账**：报告区间首尾相接且重置量已入账，因此
  `TotalDeltaMass(reports)`（各阶段 ΔMass 之和）== 末报告 `MassEnd` − 首报告
  `MassStart` == 任务从头到尾的真实总量变化；CLI 会打印该对账行。
- 断点续算时报告只覆盖本次实际执行的步区间，账目仍严格闭合；首报告的
  `MassStart` 取施加本阶段 Dirichlet 重置**之前**的场。
- CLI 先落盘检查点、再判定对账：检查点不会因对账结果而被跳过；
  对账不闭合时以 `BUDGET ERROR` 退出（退出码 1）。

## 运行

```sh
# 旧模式：单套边界/源项（行为不变）
go run ./cmd/pdesolve -nx 64 -ny 48 -dt 1e-5 -steps 100 -workers 4 -save ckpt.bin
go run ./cmd/pdesolve -resume ckpt.bin -steps 100 -workers 2

# 阶段计划：总 600 步，先跑前 350 步落盘，再恢复跑完（子域数可变）
go run ./cmd/pdesolve -steps 600 -until 350 -plan plan.json -workers 2 -save seg.bin
go run ./cmd/pdesolve -steps 600 -plan plan.json -workers 3 -resume seg.bin -save final.bin
# 恢复时修改 -dt/-D/-nx/-plan 等会被 CONFIG ERROR 拒绝
```

复现「贴边热源 + 跨阶段改边值」场景的 `plan.json` 示例（源区 `x0=0` 贴左边界、
`y0=0` 贴下边界；阶段 1 把左边界 100→250、底边界改绝热）：

```json
{"phases":[
 {"from":0,"to":300,
  "bc":{"left":{"type":"dirichlet","value":100},"right":{"type":"neumann"},
        "top":{"type":"neumann"},"bottom":{"type":"dirichlet","value":0}},
  "sources":[{"rate":50,"x0":0.0,"x1":0.3,"y0":0.3,"y1":0.7},
             {"rate":40,"x0":0.0,"x1":0.2,"y0":0.0,"y1":0.5}]},
 {"from":300,"to":600,
  "bc":{"left":{"type":"dirichlet","value":250},"right":{"type":"neumann"},
        "top":{"type":"neumann"},"bottom":{"type":"neumann"}},
  "sources":[{"rate":25,"x0":0.0,"x1":0.25,"y0":0.4,"y1":0.6}]}
]}
```

预期输出要点：两阶段 `closure` 都在 `1e-15` 量级（`OK`）；阶段 1 出现
`dirichlet reset{L:+3.7e-01 ...}`（左边界改值的重置量，底边改绝热为 0）；
`task:` 行 `sum(phase deltaMass)` 与 `true deltaMass` 残差为 0；
分段（`-until 350` + `-resume`）与一次性跑到底的最终检查点逐位一致。

## 验证方法

```sh
go test -race -v ./solver/
```

- `TestInvalidConfigsAreRejected` / `TestInvalidPhasePlans`：非法配置与非法阶段计划
  （空计划、间隙、重叠、零长度、越界、非法边界/源数值等）均被拒绝且归类 ErrKindConfig。
- `TestStabilityLimitEnforced` / `TestPhasePlanStabilityDistinct`：超限 dt 归类
  ErrKindStability，与配置错误可区分。
- `TestPhaseSwitchingAndBudget`：三阶段（加热→冷却→改边界）切换语义与逐阶段收支闭合。
- `TestEdgeSourceBudget`：贴边/贴角（含相互重叠、分别贴不同边）热源按实际施加量
  记账，逐源注入量与"非固定格解析值"一致，阶段收支闭合。
- `TestPhaseSwitchResetAccounting`：跨阶段修改固定边值 / 固定改绝热时重置量逐边
  可见（与解析值一致），各阶段 ΔMass 之和等于任务真实总量变化。
- `TestCheckBudgetsLocatesImbalance`：对账真不闭合时报 ErrKindBudget 并定位到
  具体阶段，与配置/失稳错误可区分。
- `TestMultiSourceAccounting`：多源区注入量与解析值一致。
- `TestResumeConsistency` / `TestPhaseCheckpointResume` /
  `TestEdgeSourceCheckpointResume`：阶段内部与阶段边界两种位置的保存/恢复与
  不间断运行逐位一致（含贴边热源场景），恢复后改子域数仍一致，检查点确实落盘。
- `TestResumeConfigMismatch`：网格/步长/系数/阶段计划不一致在续算前被拒绝。
- `TestInvariantConservation`：全 Neumann 无源时总量守恒（偏差 ~1e-15，容差 1e-9）。
- `TestLegacyRunPhasesCompatible`：无阶段计划时新入口与旧入口逐位一致。
- `TestSerialParallelConsistency` / `TestPhaseParallelConsistency` /
  `TestEdgeSourceParallelConsistency`：不同子域数与串行基线逐位一致。
- `TestNaNGuard`：运行期防护捕获非有限值并报告步号与坐标。

测试日志打印输入参数、判定依据（容差）与实测偏差数值。
