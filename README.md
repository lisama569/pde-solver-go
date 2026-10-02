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
- 阶段激活时，该阶段 `Dirichlet` 边的 `Value` 会**写入边界格**（幂等：同阶段重复
  激活结果相同）；角点被两条 Dirichlet 边共享时按 left→right→bottom→top 顺序后者覆盖。
- **切换瞬间边界值的跳变（边界重置）逐边计入新阶段报告的 `BoundaryReset` 账目**，
  口径为 Σ(新值−旧值)·dx·dy（按格面积，与总量 `TotalMass` 一致；角点跳变归实际
  后写入的那条边）。覆盖三种情形：跨阶段修改 Dirichlet 数值、Dirichlet→Neumann
  （释放固定格，不写入、跳变为 0）、Neumann→Dirichlet（把当前值钉到新值，跳变=
  新值−切换前该格实际值）。任务首次激活阶段0 时，Dirichlet 值相对初始场的跳变
  同样计入阶段0，因此各阶段 ΔMass 之和恰好等于任务从头到尾的真实总量变化。
- 阶段内多个源区可重叠，强度相加。

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
- 两类错误可通过 `solver.IsConfigError` / `solver.IsStabilityError` 区分，
  阶段计划问题永远归类为配置错误，与运行期失稳严格区分。

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

`RunPhases` 为每个阶段输出 `PhaseReport`：各热源累计**实际注入**（逐源）、各边界
累计扩散净流入（逐边，正 = 流入）、阶段激活时的**边界重置跳变**（逐边）、阶段首末
总量。对账恒等式：

```
ΔMass = Σ源注入(实际施加) + Σ边界扩散净流入 + Σ边界重置跳变
残差 Closure = ΔMass − 源 − 扩散 − 重置
```

**贴边热源的口径**：源项只施加在更新算子实际更新的格点上。Dirichlet 边界格值被
钉住、不接受源项，因此贴 Dirichlet 边（或同时压住两条 Dirichlet 边的角部）的源区，
其注入量**只统计区域内非固定格**（`rate × 实际受激格点数 × dx·dy × dt`，逐步累计），
不会把压在固定格上的名义注入量计入账目；贴 Neumann / Periodic 边的源区照常全额
计入（这些边的格点不是固定格）。多个源区相互重叠时强度相加，逐源独立记账。
Neumann/Periodic 边恒为 0，Periodic 边两两相消恒为 0。

**跨阶段改边界的口径**：把某条边的 Dirichlet 值改大/改小、或 Dirichlet↔Neumann
切换，写入瞬间造成的全场总量跳变不属于扩散交换，而以 `BoundaryReset` 单独逐边列示。
`solver.VerifyReports` 除逐阶段校验闭合外，还校验：

1. 各报告区间首尾相接；
2. `Σ各阶段 ΔMass == 末份MassEnd − 首份MassStart`（首份 MassStart 在任何重置
   之前测量），即"把各阶段报的总量变化加起来等于任务从头到尾的真实变化"。

**容差判定**：

- `ClosureTol = 1e-9 × (1 + |MassStart| + |MassEnd| + Σ|源注入| + Σ|扩散净流入| + Σ|重置跳变|)`，
  比浮点舍入误差（~1e-15 相对）高约 6 个数量级，又足以发现任何物理量级（~1e-6 相对）
  的漏账；实测正常闭合残差在 1e-15~1e-16 量级。
- 不闭合时返回 **ErrKindBudget**（`solver.IsBudgetError`），错误信息带**具体阶段
  索引、步区间、逐源注入、逐边扩散与重置明细**，可直接定位到具体热源或边界；
  当残差量级与某一账项吻合时还会给出疑似项提示。该类别与 `ErrKindConfig`（非法
  配置，推进前拒绝）、`ErrKindStability`（NaN/Inf 或违反 CFL）严格区分。
- 断点续算时报告只覆盖本次实际执行的步区间，恢复时重新激活当前阶段是幂等的
  （已在边界值上，重置跳变为 0），各段账目仍严格闭合。

## 检查点落盘时机与失败区分

CLI（`cmd/pdesolve`）**先落盘检查点、后做对账校验**：对账真不闭合（或误报）都不会
跳过检查点，恢复时不会出现"找不到文件"。退出原因分四类打印前缀：

| 前缀 | 类别 | 触发时机 |
|---|---|---|
| `CONFIG ERROR` | `ErrKindConfig` | 推进之前：非法网格/系数/阶段计划、续算配置指纹不一致 |
| `STABILITY ERROR` | `ErrKindStability` | 推进中：dt 超 CFL 极限、出现 NaN/Inf（含步号与格点坐标） |
| `BUDGET ERROR` | `ErrKindBudget` | 推进完成后：逐阶段或整任务对账不闭合（含阶段/热源/边界定位） |
| `ERROR` | 其他 | IO 等普通错误（如检查点文件不存在） |

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

### 贴边热源 + 跨阶段改边界的复现计划

下面的两阶段计划复现"先加热（热源贴左/下固定温度边，含贴角）、再改边界继续算"：

```json
{
  "phases": [
    {
      "from": 0, "to": 30,
      "bc": {"left": {"type": "dirichlet", "value": 100}, "right": {"type": "neumann"},
             "top": {"type": "neumann"}, "bottom": {"type": "dirichlet", "value": 0}},
      "sources": [
        {"rate": 50, "x0": 0.0, "x1": 0.3, "y0": 0.4, "y1": 0.6},
        {"rate": 20, "x0": 0.7, "x1": 0.9, "y0": 0.0, "y1": 0.2}
      ]
    },
    {
      "from": 30, "to": 60,
      "bc": {"left": {"type": "dirichlet", "value": 200}, "right": {"type": "neumann"},
             "top": {"type": "neumann"}, "bottom": {"type": "neumann"}},
      "sources": []
    }
  ]
}
```

```sh
# 一次跑到底：每阶段打印 influx / reset 两组逐边账目与 closure/tol，结尾打印整任务总账
go run ./cmd/pdesolve -nx 32 -ny 24 -dt 1e-5 -steps 60 -plan plan-edge.json -workers 2 -save run.bin
# 在阶段边界（30）落盘，恢复时换 4 个子域；轨迹与一次性跑到底逐位一致
go run ./cmd/pdesolve -nx 32 -ny 24 -dt 1e-5 -steps 60 -until 30 -plan plan-edge.json -save seg.bin
go run ./cmd/pdesolve -nx 32 -ny 24 -dt 1e-5 -steps 60 -plan plan-edge.json -workers 4 -resume seg.bin
```

预期日志中：阶段0 两个源的 `injected` 只含非固定格贡献（贴左边源不含左列、贴底边
源不含底行）；阶段1 `reset{L:...}` 恰好等于左列 100→200 的总量跳变；各阶段
`closure` 在 ~1e-15 量级、远小于 `tol`；`task ledger` 行显示 ΣΔMass 与真实变化
相等。把 `-dt` 改大超过 CFL 极限得到 `STABILITY ERROR`；恢复时改任意网格/步长/系数/
计划得到 `CONFIG ERROR`；账目真被破坏时得到带阶段与源/边定位的 `BUDGET ERROR`，
且检查点仍已落盘。

## 验证方法

```sh
go test -race -v ./solver/
```

- `TestInvalidConfigsAreRejected` / `TestInvalidPhasePlans`：非法配置与非法阶段计划
  （空计划、间隙、重叠、零长度、越界、非法边界/源数值等）均被拒绝且归类 ErrKindConfig。
- `TestStabilityLimitEnforced` / `TestPhasePlanStabilityDistinct`：超限 dt 归类
  ErrKindStability，与配置错误可区分。
- `TestPhaseSwitchingAndBudget`：三阶段（加热→冷却→改边界）切换语义与逐阶段收支闭合。
- `TestMultiSourceAccounting`：多源区注入量与解析值一致。
- `TestSerialParallelConsistency` / `TestPhaseParallelConsistency`：不同子域数与
  串行基线逐位一致（含带阶段计划场景）。
- `TestResumeConsistency` / `TestPhaseCheckpointResume`：阶段内部与阶段边界两种位置
  的保存/恢复与不间断运行逐位一致，恢复后改子域数仍一致。
- `TestResumeConfigMismatch`：网格/步长/系数/阶段计划不一致在续算前被拒绝。
- `TestInvariantConservation`：全 Neumann 无源时总量守恒（偏差 ~1e-15，容差 1e-9）。
- `TestLegacyRunPhasesCompatible`：无阶段计划时新入口与旧入口逐位一致。
- `TestNaNGuard`：运行期防护识别为失稳。
- `TestSingleCornerSourceTwoPhase`：**单热源贴角**（同时压左/底两条 Dirichlet 边）
  两阶段计划不误报不平衡、账闭合、检查点正常落盘恢复（用户场景最小复现）。
- `TestEdgeTouchingSourceBudget`：**多热源区**分别贴左/底 Dirichlet 边（含相互重叠）、
  贴 Neumann 边、贴左下角点；注入量只统计实际受激格点、与解析值一致，全部阶段闭合。
- `TestPhaseBoundaryResetJump`：跨阶段**改 Dirichlet 值、Dirichlet↔Neumann 双向切换**
  的逐边重置跳变与解析值一致，且 ΣΔMass 等于真实变化；Neumann→Dirichlet 跳变
  按切换前实际场值计算。
- `TestEdgePlanCheckpointResume`：贴边计划在**阶段内部与阶段边界**两种位置保存/恢复，
  切换步号与最终轨迹和一次性跑到底**逐位一致**，恢复时改并行子域数（2/3/5/7）
  仍与单子域基线一致，各续算段报告均闭合。
- `TestEdgePlanParallelConsistency`：贴边计划在 workers=2/3/4/8 下与单子域基线
  逐位一致（`go test -race` 无数据竞争）。
- `TestBudgetImbalanceDetectedAndLocated`：篡改源注入/重置跳变后 `VerifyReports`
  返回可定位到阶段与具体源（source N）/边界（boundary N influx/reset）的
  `ErrKindBudget`，且不与配置/失稳错误混淆。

测试日志打印输入参数（网格、dx/dy/dt、D、步数、子域数、阶段数）、采用的容差
（closure vs tol）、关键判定依据（格点计数、逐源解析注入、逐边 influx/reset）、
实测偏差与每个业务用例的 `CASE ...: PASS` 通过标记。
