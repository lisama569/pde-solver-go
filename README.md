# pde-solver-go

一个在本地运行、内存网格模拟、可断点续算、支持子域并发的二维热传导方程
（含源项）数值求解器。作为 Go 库 `github.com/lisama569/pde-solver-go/pde`
提供能力，并附带命令行入口 `cmd/pdesim`。

- 多种边界条件（Dirichlet / Neumann / 周期）与源项（常数 / 高斯）可自由组合；
- 步长或系数越过稳定范围时**在执行前拒绝**；非法配置与运行时失稳返回
  **可区分类型**的错误；任何情况下都不会静默产生 NaN / 非物理结果；
- 区域分解并发求解，串行基线与任意并行度结果**逐位一致**（无 halo 通信、
  无竞态、与调度顺序无关）；
- 状态原子落盘、可断点续算；恢复后的轨迹与不间断运行一致；
- 每一步可做离散热平衡（守恒量）校验，残差按 x 子域切片定位到具体步骤与区域。

## 1. 控制方程与数值格式

求解二维抛物方程

```
∂u/∂t = α(∂²u/∂x² + ∂²u/∂y²) + q(x,y,t)
```

网格为均匀矩形网格，物理边界位于最外层**内部单元中心外侧半格**处。
时间推进采用显式 FTCS（前向时间、中心空间）五点格式：

```
u^(n+1)_ij = u^n_ij
  + rx (u^n_{i-1,j} + u^n_{i+1,j} - 2 u^n_ij)
  + ry (u^n_ij-1   + u^n_ij+1   - 2 u^n_ij)
  + dt * q(x_i, y_j, t^n)

rx = α dt / dx²,  ry = α dt / dy²
```

源项在前一时刻 `t^n` 取值（前向 Euler）。求解器使用**双缓冲**
（cur/nxt 两块独立网格）：每一步从只读的旧场算出新场，成功后交换缓冲；
单步失败时新缓冲被丢弃，求解器状态停留在上一步，保证不会出现半步状态。

### 边界条件与 ghost 层

网格四周各有一层 ghost。每步更新前由内部状态与边界条件一次性填充：

| 条件 | ghost 取值（以 x0 面为例） | 含义 |
| --- | --- | --- |
| Dirichlet `v` | `u_ghost = 2v - u(0)` | 边界值 v 位于 ghost 与首内部点的中点（半格约定） |
| Neumann `g = ∂u/∂n` | `u_ghost = u(0) + g·dx` | g 为**外法向**导数；四个面符号统一为“加” |
| Periodic | `u_ghost = u(N-1)` / `u(N-1 ghost)=u(0)` | 周期环绕 |

周期边界必须在同一轴的两个面成对出现，否则在配置阶段报错。
角点 ghost 不被五点格式读取，保持零值。

### 源项

- `SourceNone`：无源；
- `SourceConstant`：`q = Scale`；
- `SourceGaussian`：`q = Scale · exp(-((x-Cx)²+(y-Cy)²)/(2σ²))`。

## 2. 稳定性条件（执行前拒绝）

2D 显式 FTCS 的 CFL 型必要（对本格式亦充分）条件：

```
rx + ry ≤ 1/2

即  dt ≤ 1 / (4 α (1/dx² + 1/dy²))
```

`Config.CheckStability()` 在构造求解器时强制检查；不满足时返回
`KindStability` 错误，消息中给出 `rx`、`ry`、二者之和与允许的 `dt` 上界，
**不会开始积分**。α = 0 时无扩散稳定性限制。

CLI 传 `-dt 0`（默认）时自动选择步长：`nsteps = ceil(T/(0.9·dt_lim))`，
`dt = T/nsteps`，既保留 10% 稳定裕量又保证 `dt` 整除总时间。

运行期还会对每个新值检查：非有限值（NaN/±Inf）或 `|u| > MaxValue`
（默认 1e12）立即以 `KindBlowup` 失败，错误信息带**步骤号、子域号、网格
坐标与时刻**；失败步不落盘、不交换缓冲，现场不可见任何 NaN。

### 失败类型（可区分）

| `FailureKind` | 含义 |
| --- | --- |
| `config` | 非法配置：尺寸、间距、α、dt、边界配对、源参数、dt 不整除总时间等 |
| `stability` | 配置合法但 `rx+ry > 1/2`，执行前拒绝 |
| `blowup` | 运行期产生 NaN/Inf 或越过有限值界限（含步骤/子域/坐标） |
| `checkpoint` | 检查点读写、版本、载荷非有限值、网格/物理配置不匹配 |
| `invariant` | 离散热平衡残差超过容差（含步骤与最严重子域） |

用 `pde.AsFailure(err)` 取出结构化错误，读取 `Kind/Step/Domain/Msg`。

## 3. 并发模型与一致性保证

- 按 x 方向把内部列切成 `Workers` 个连续、不交、近均的 slab
  （`partitionRows`，划分只依赖 `(Nx, Workers)`，确定可复现）；
- 每个 goroutine **只读共享的旧场 `src`**（含已填充的 ghost 层），
  **只写自己 slab 的新场 `dst`**。双缓冲 + 不交写集合意味着：
  - 不需要 halo 通信，不存在重复更新或写竞争；
  - 所有浮点运算的操作数和顺序与串行完全一致，
    因此串行与任意并行度（含 slab 切过周期接缝的情形）结果**逐位相同**；
  - 重复运行、改变 goroutine 调度均不影响结果。
- 错误按子域编号固定顺序回收（先报小编号），panic 被 recover 转为
  `KindBlowup`，不会崩溃进程。

数值容差约定：单测以 `1e-14` 作为串并行数值一致性阈值，而实测最大差为
`0`（逐位一致）；判定依据为全场 `max_{i,j} |u_w(i,j)-u_1(i,j)|`。

## 4. 检查点与断点续算

`Snapshot` 以确定性 JSON（schema_version=1）保存：完整配置、步号、物理
时间、行优先内部场。`SaveCheckpoint` 采用“同目录临时文件 + fsync +
原子 rename”，崩溃不会损坏旧检查点。

加载时校验：schema 版本、配置稳定性、步号/时间一致性、数据长度、
**每个数值必须有限**。恢复时只比较决定轨迹身份的物理量（网格、间距、
α、dt、边界、源、初值、种子）；`Workers`、监控阈值等不影响轨迹，
续算时允许不同。计划终点 `TotalTime` 允许在续算时延长（但 dt 必须一致，
且检查点步号不得超过新计划总步数）。

由于格式是确定的（FTCS、确定性初值、定序规约），从第 k 步恢复后继续
积分与一次性积分到终点**逐位一致**。

## 5. 守恒量（离散热平衡）校验

`CheckConservation` 逐步检验离散恒等式

```
H^n+1 − H^n = Q_source + Φ_boundary
H = dx·dy · Σ_ij u_ij
```

其中边界热通量使用与更新**完全相同的 ghost/内部值**按半格约定计算：

```
x 面进入热:  α dt/dx · (u_ghost − u_interior) · dy        （逐行求和）
y 面进入热:  α dt/dy · (u_ghost − u_interior) · dx        （逐列求和）
源项注入:    dx·dy·dt · Σ q(x_i,y_j,t^n)
```

残差还会按若干 x-slab 分解（内部界面通量成对抵消，各 slab 残差之和恒等于
全局残差），从而把偏差**定位到具体步骤与区域**。对本离散格式实测相对残差
为浮点舍入量级（≈1e-16）。默认相对容差 `1e-10`（相对量以
`max(1, |期望热|)` 归一化），可用 `ConserveTol` 调整；`ConserveEnforce`
（CLI: `-enforce-conservation`）使超限直接以 `KindInvariant` 终止。

物理上：全周期或全零通量 Neumann 问题总热只随源项变化；Dirichlet/非零
Neumann 边界存在与上述公式一致的热交换，这是“可解释”的平衡项而非误差。

## 6. 使用

### 库

```go
cfg := pde.Config{
    Nx: 24, Ny: 20, Dx: 1.0/24, Dy: 1.0/20, Alpha: 0.2,
    TotalTime: 0.05, Dt: 0.002,
    Boundaries: pde.Boundaries{ /* X0/X1/Y0/Y1 */ },
    Source:     pde.SourceConfig{Type: pde.SourceGaussian, Scale: 3,
        Cx: 0.5, Cy: 0.5, Sigma: 0.12},
    InitValue: 0.2, Seed: 42, NoiseAmp: 0.05,
    Workers: 4, ConserveEnforce: true,
}
res, err := pde.RunSimulation(cfg, pde.RunOptions{
    Checkpoint: "state.json", CheckpointEvery: 10,
    CheckInvariant: true, InvariantRegions: 4,
})
// 续算：
res2, err := pde.RunSimulation(cfg, pde.RunOptions{
    Checkpoint: "state.json", Resume: true, Workers: 8,
})
```

### 命令行

```
pdesim run    [flags]   # 从 t=0 运行
pdesim resume [flags]   # 从 --checkpoint 恢复

  -nx,-ny,-dx,-dy,-alpha,-dt,-total-time
  -init,-seed,-noise
  -x0,-x1,-y0,-y1   d=V（Dirichlet） | n=G（外法向导数） | p（周期）
  -source-scale,-source-cx,-source-cy,-source-sigma  （scale=0 即无源）
  -workers              子域并行度（1 为串行基线）
  -max-value            运行期有限值界限
  -check-invariant,-enforce-conservation,-conserve-tol,-invariant-regions
  -checkpoint,-checkpoint-every
```

示例：

```
go run ./cmd/pdesim run -nx 20 -ny 20 -dx 0.05 -dy 0.05 -alpha 0.2 \
  -dt 0.002 -total-time 0.024 -x0 d=2 -x1 d=0 -y0 n=0.1 -y1 n=0 \
  -source-scale 5 -workers 3 -check-invariant \
  -checkpoint state.json -checkpoint-every 6

go run ./cmd/pdesim resume ... -total-time 0.048 -workers 5 -checkpoint state.json
```

标准错误流会打印：全部输入参数、`rx/ry/和/dt 上界`稳定性判定依据、
每步守恒残差与最终汇总。

## 7. 验证方法

```
go test ./...          # 全部单测
go test -race ./...    # 竞态检测（已通过，count 可重复）
go vet ./...
```

测试矩阵（日志中均打印输入参数、判定依据与偏差数值）：

1. **非法配置**（`TestIllegalConfig`）：11 类非法输入逐一断言
   `KindConfig`；**合法但不稳定**（`TestStabilityRejection`、
   `TestStabilityIsDistinguishableFromConfig`）断言 `KindStability`，
   验证两类失败可区分；`TestAutoDtStableAndDivisible` 验证自动步长稳定
   且整除。
2. **串并行一致性**（`TestSerialParallelConsistency`）：周期 x 跨越子域
   接缝 + Dirichlet/Neumann y + 非均匀初值 + 局部源，workers ∈
   {1,2,3,4,5,7,8,24}，全场最大差 `0`（阈值 1e-14）；重复运行差为 0。
   `TestPartitionRowsDeterministic` 验证划分连续、不交、确定。
3. **续算一致性**（`TestCheckpointResume`、`TestCheckpointRepeatedResume`、
   `TestResumeWithExtendedTotalTime`、`TestRunnerResumeEndToEnd`、
   CLI `TestCLISplitResume`）：在 1/2/7/13/24/25 步处断点、多轮 save/load、
   续算时更换并行度与延长终点，最终场与连续运行差为 `0`。
4. **守恒量**（`TestConservationAllBoundaryKinds`）：全 D、全 N、
   周期 x+D/N y、全周期 + 各类源，逐步相对残差 ≲ 4e-16（阈值 1e-12），
   并验证子域残差分解之和等于全局残差；
   `TestInvariantViolationLocalized` 人工注入非物理热，断言 `KindInvariant`
   且最严重 slab 被准确指出（注入 2.5 即报 2.5）。
5. **运行期失稳**（`TestRuntimeBlowup`）：稳定 dt + 极大源，验证
   `KindBlowup` 带步骤/子域/坐标，且字段中不存在 NaN/Inf。
6. **检查点损坏**（`TestCorruptAndMismatchedCheckpoints`）：非法 JSON、
   错误 schema、非有限载荷、网格不匹配均返回 `KindCheckpoint`；
   NaN 在编码与解码两侧都被拒绝。

## 8. 范围与限制

- 二维、均匀网格、常系数 α 的热传导/反应-扩散模型；
- 显式格式，刚性问题受 CFL 约束（需要时可在此包内扩展隐式格式）；
- ghost 约定决定了离散边界通量与守恒平衡必须配套使用，
  若修改边界离散需同步修改 `invariant.go`。
