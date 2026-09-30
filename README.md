# pde-solver-go

本地可运行的二维热传导（扩散）方程数值求解流程。使用内存网格模拟随时间演化的标量场，
支持串行/按子域并发求解、检查点落盘与断点续算、守恒量校验，并支持**按时间阶段切换
矩形热源与四边边界条件**，保证不同入口与不同并行度下结果一致、可解释、可复现。

## 数值格式

求解方程：`∂u/∂t = D·∇²u + S(x,y,t)`（`S` 为零个或多个矩形常源项，可随阶段变化）。

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

源项 `Source` / `PhaseSource` 以归一化坐标 `[X0,X1]×[Y0,Y1]` 指定矩形区域与强度 `Rate`。
一个阶段内可配置**多个源区**，允许重叠，重叠格点注入强度按各源 `Rate` 相加；
`Rate` 可正（加热）可负（冷却/热汇）。

## 时间阶段计划（PhasePlan）

一次任务可按**绝对时间步号**切分成连续、不重叠的多个阶段，典型场景：先加热中心区域、
再关源冷却、中途把某条边的温度改掉。

- 阶段 `Phase{StartStep, EndStep, BC, Sources}` 覆盖步号区间 **`[StartStep, EndStep)`**：
  从绝对步 `StartStep` 推进到 `StartStep+1` 的那一步起，使用本阶段的边界与源，直到到达
  绝对步 `EndStep`（那一步属于下一阶段）。
- 阶段序列必须：从步 0 开始、首尾相接（后一个 `StartStep` 必须等于前一个 `EndStep`，
  不允许空档或重叠）、每段至少一步、边界类型合法且周期边成对、源区坐标与强度为有效有限值、
  同阶段源名不重名。任何不合法计划在**开始推进前**以 **ErrKindConfig** 拒绝。
- 阶段计划还必须覆盖到 `cfg.Steps`；计划短于目标步号同样按配置问题拒绝。

### 阶段切换时的语义（明确、可解释）

- **源项**：切换发生在步边界，无重叠、无插值。上一阶段最后一步用旧源，下一阶段第一步
  立即用新源；某阶段 `sources` 为空即该阶段完全无源。
- **Neumann / Periodic 边**：切换类型只改变差分取邻规则，不直接改场。
- **Dirichlet 边**：进入新阶段时，把该边上所有边界格点**一次性吸附（snap）到新阶段的
  固定值**，然后这些格点继续被冻结。吸附导致的全场总量跳变被单独记账为该阶段的
  `BoundaryReset`（见对账公式），因此切换即使瞬时改变能量，账目依然闭合。
  - 四角若两边同为 Dirichlet，按 left → bottom → right → top 的顺序赋值，取后写边的值。
  - 阶段 0 也会执行一次吸附：若你的初始场 Dirichlet 边值本来就等于声明值，跳变为 0。
- 阶段计划通过 JSON 文件传入（`-plan`），边界类型可写字符串
  `dirichlet`/`neumann`/`periodic`，也兼容旧整数枚举 `0/1/2`：

```json
{
  "phases": [
    {"start_step": 0, "end_step": 30,
     "bc": {"left":{"type":"dirichlet","value":100},"right":{"type":"neumann"},
            "top":{"type":"neumann"},"bottom":{"type":"dirichlet","value":0}},
     "sources":[{"name":"heater","rate":50,"x0":0.4,"x1":0.6,"y0":0.4,"y1":0.6}]},
    {"start_step": 30, "end_step": 60,
     "bc": {"left":{"type":"dirichlet","value":100},"right":{"type":"neumann"},
            "top":{"type":"neumann"},"bottom":{"type":"dirichlet","value":0}},
     "sources":[]},
    {"start_step": 60, "end_step": 100,
     "bc": {"left":{"type":"dirichlet","value":20},"right":{"type":"neumann"},
            "top":{"type":"dirichlet","value":10},"bottom":{"type":"dirichlet","value":0}},
     "sources":[{"name":"corner-heat","rate":12,"x0":0.1,"x1":0.3,"y0":0.1,"y1":0.3}]}
  ]
}
```

### 与旧配置的兼容范围

- **不传 `-plan`（`plan == nil`）时，完全走旧的单套 `BC` / `Src` 路径，行为与历史版本
  逐位一致**（既有边界、源项、并行、续算、守恒用例零回归）。
- 唯一新增的语义差异只存在于**阶段计划模式**：进入阶段（含阶段 0）时 Dirichlet 边会被
  吸附到声明值并计入 `BoundaryReset`；旧模式只冻结初始场边值、不做吸附。
- 旧版 v1 检查点（`SaveCheckpoint`，无配置签名）不能用于阶段模式恢复，会被明确拒绝；
  旧模式之间仍可继续用 v1 续算。阶段模式使用带配置签名的 v2 状态（`SaveState`/`LoadState`）。

## 保存 / 恢复（任意阶段内或跨阶段边界）

- v2 状态文件除步号、维度、场外，还保存配置签名：`NX,NY,DX,DY,DT,D` 与**阶段计划指纹**
  （对所有阶段边界类型/值、源区/强度/步号做 FNV-1a 哈希）。
- 恢复时 `ResumeRunner`（内部 `State.Verify`）逐字段比对：网格尺寸、`DX/DY`、`DT`、
  扩散系数、是否有阶段计划、计划指纹。**任一对不上都在继续计算前以 ErrKindConfig
  明确拒绝并指出差异**，绝不悄悄按新配置算。并行子域数不属于配置签名，可以随意更改。
- 恢复点既支持**阶段内部**（`StartStep < step < EndStep`），也支持**刚好停在阶段边界**
  （`step == 某 StartStep`）。边界恢复时会与一次性运行一样补做新阶段的 Dirichlet 吸附，
  因此恢复后的阶段切换步号、吸附时刻与最终轨迹和不间断运行**逐位一致**。

## 阶段收支对账（能量平衡）

每个阶段在其 `EndStep` 处给出一份 `PhaseReport`，可逐阶段、逐热源、逐边对账。离散恒等式：

```
EndMass - StartMass = Σ 各热源注入 + Σ 四边净流入 + BoundaryReset  (+ 浮点舍入)
```

- 质量定义：`Mass = Σ u·dx·dy`。
- 热源注入：阶段内每步落在源区内、且非 Dirichlet 固定格点上的 `dt·dx·dy·Rate` 之和，
  按源名分别列出（未命名自动编号 `source#i`），便于定位到**具体热源**。
- 四边净通量（**流入域内为正**）：Dirichlet 边为 `D·dt·(dy/dx)·Σ(u_边 − u_内侧邻点)`
  （顶/底边用 `dx/dy`）；Neumann 边为 0；Periodic 边逐对相消为 0。求和行/列集会随对边
  是否 Dirichlet 收缩，角点附近不漏项，因此四边混合 Dirichlet 时也能严格闭合，便于定位到
  **具体边**。
- `BoundaryReset`：进入该阶段时 Dirichlet 吸附导致的总量跳变。
- 判定（容差可通过 `Runner.SetBalanceTolerance` 覆盖）：

```
tol = tolAbs + tolRel · scale,  默认 tolAbs = 1e-12, tolRel = 1e-9
scale = max(|StartMass|, |EndMass|, |总注入|, |总通量|, |重置|)
Closed ⇔ |Residual| ≤ tol
```

由于上述三项与实际逐格点更新严格对应，残差只有浮点舍入量级（实测约 1e-16 ~ 1e-15）；
默认容差留了约 6 个数量级余量。对不上时报告直接标明**哪个阶段、哪个源/哪条边**。

## 稳定性条件

显式格式要求（`Config.StabilityLimit` / `CheckStability`）：

```
dt <= 1 / (2·D·(1/dx² + 1/dy²))
```

- 启动前：`CheckStability` 拒绝超限的 `dt`，返回 **ErrKindStability**（阶段计划同样适用）。
- 运行中：每步提交后扫描全场，一旦出现 NaN/Inf 立即中止，报错含**步号与格点坐标**。
- 非法参数 / 非法阶段计划 / 恢复配置不符返回 **ErrKindConfig**；数值失稳返回
  **ErrKindStability**。两者用 `solver.IsConfigError` / `solver.IsStabilityError` 区分，
  CLI 分别打印 `CONFIG ERROR` / `STABILITY ERROR`。

## 并行约束

`RunParallel` 与阶段运行器均按行把网格均分为 `workers` 个子域并发推进：

- 双缓冲保证各子域只写自己的行区间，**无数据竞争、无重复更新**（`go test -race` 验证）。
- 每个格点的浮点运算序列与串行完全相同，因此结果与串行基线**逐位一致**（容差 0），
  与并发度、执行顺序无关；**恢复后改并行子域数结果也与单子域基线一致**。

## 检查点格式

- v1 `SaveCheckpoint`（magic `PDECKPT1`）：magic + 步号 + 维度 + float64 数据；旧模式使用。
- v2 `SaveState`（magic `PDESTATE`）：在 v1 基础上增加配置签名字段
  （hasPlan、planHash、DX、DY、DT、D）；阶段模式使用，支持恢复前一致性校验。
- 两者均先写临时文件再原子改名；读取时校验 magic 与长度，损坏或混用给出描述性错误。

## 运行

```sh
# 旧单套配置（行为不变）
go run ./cmd/pdesolve -nx 64 -ny 48 -dt 1e-5 -steps 100 -workers 4 -save ckpt.bin
go run ./cmd/pdesolve -resume ckpt.bin -steps 100 -workers 2

# 阶段计划模式：加热 -> 冷却 -> 中途改边温（先跑到阶段边界 step 30 并保存）
go run ./cmd/pdesolve -nx 40 -ny 30 -dt 1e-5 -steps 30 -workers 1 -plan plan.json -save s30.bin
# 阶段边界保存后恢复跑到 step 100，并把并行子域数从 1 改成 7
go run ./cmd/pdesolve -nx 40 -ny 30 -dt 1e-5 -steps 100 -workers 7 -plan plan.json -resume s30.bin
# 配置（dt/D/网格/计划）与保存时不一致 -> CONFIG ERROR，拒绝续算
```

## 验证方法

```sh
go test -race -v ./solver/
```

既有用例（不回归）：

- `TestInvalidConfigsAreRejected`：9 类非法配置均被拒绝且归类为 ErrKindConfig。
- `TestStabilityLimitEnforced`：超限 dt 归类为 ErrKindStability，与配置错误可区分。
- `TestSerialParallelConsistency`：workers ∈ {1,2,3,4,7,16} 与串行基线逐位一致。
- `TestResumeConsistency`：分段落盘续算与不间断运行逐位一致。
- `TestInvariantConservation`：全 Neumann 无源时总量守恒（偏差 ~1e-15，容差 1e-9）。
- `TestNaNGuard`：绕过前置校验时运行期防护捕获非有限值并报告步号与坐标。

阶段能力新增用例：

- `TestPhaseSwitchingAndBalance`：加热/关源/改边温三阶段，切换窗口正确，各阶段逐源、
  逐边对账闭合（残差 ~1e-16，容差 1e-12 + 1e-9·scale）。
- `TestPhaseBalanceAllDirichlet`：四边 Dirichlet + 角点源区，验证角点通量不漏项。
- `TestPhaseBalancePeriodic`：周期边界边通量恒为 0、重置为 0，收支闭合。
- `TestInvalidPhasePlansRejected`：空档、重叠、越界、空区间、未知边界、未配对周期、
  非法源区、NaN 强度、Inf 固定值、源重名、计划盖不住步数等均按配置问题拒绝。
- `TestPhaseStabilityStillDetected`：阶段模式下超限 dt 与运行期 NaN 都归类为失稳。
- `TestPhaseWorkersConsistency`：阶段运行 workers ∈ {1,2,3,4,7,16} 与串行逐位一致。
- `TestPhaseSaveResume`：阶段内（15/35）与阶段边界（20/40）保存恢复、恢复后改用 5 个子域，
  终场与不间断串行基线逐位一致，恢复后各阶段仍对账闭合。
- `TestResumeRejectsChangedConfig`：网格/DT/DX/D/源强度/边值/阶段边界改动及计划缺失，
  恢复前全部被 ErrKindConfig 拒绝。
- `TestPhasePlanJSONRoundTrip`：字符串与整数边界写法均可解析，往返后计划指纹不变。
- `TestPhasePlanLegacyEquivalence`：边初值等于声明值时，单阶段计划与旧 Stepper 路径逐位一致。

测试日志打印输入参数（网格、步长、系数、步数、CFL 上限）、采用的容差、关键判定依据
（逐源注入、逐边通量、重置跳变）与实测残差数值。
