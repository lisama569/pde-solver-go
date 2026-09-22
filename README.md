# pde-solver-go

本地可运行的二维热传导（扩散）方程数值求解流程。使用内存网格模拟随时间演化的标量场，
支持串行/按子域并发求解、检查点落盘与断点续算、守恒量校验，保证不同入口与不同并行度下
结果一致、可解释、可复现。

## 数值格式

求解方程：`∂u/∂t = D·∇²u + S(x,y)`（`S` 为可选矩形常源项）。

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

## 稳定性条件

显式格式要求（`Config.StabilityLimit` / `CheckStability`）：

```
dt <= 1 / (2·D·(1/dx² + 1/dy²))
```

- 启动前：`CheckStability` 拒绝超限的 `dt`，返回 **ErrKindStability**。
- 运行中：每步提交后扫描全场，一旦出现 NaN/Inf 立即中止，报错含**步号与格点坐标**。
- 非法参数（网格过小、负系数、未配对周期边界、非法源区域等）返回 **ErrKindConfig**。
- 两类错误可通过 `solver.IsConfigError` / `solver.IsStabilityError` 区分，绝不静默产生非物理结果。

## 并行约束

`RunParallel` 按行把网格均分为 `workers` 个子域并发推进：

- 双缓冲保证各子域只写自己的行区间，**无数据竞争、无重复更新**（`go test -race` 验证）。
- 每个格点的浮点运算序列与串行完全相同，因此结果与串行基线**逐位一致**（容差 0），
  与并发度、执行顺序无关。

## 检查点与断点续算

- `SaveCheckpoint` 以二进制格式（magic + 步号 + 维度 + float64 数据）原子落盘；
  `LoadCheckpoint` 校验 magic 与长度，损坏文件给出描述性错误。
- 恢复后从检查点步号继续推进，轨迹与不间断运行**逐位一致**。

## 守恒量校验

`TotalMass` = Σu·dx·dy。全 Neumann 边界且无源时为严格守恒量（离散通量两两相消，
边界通量为零），测试容差 1e-9（实测 ~1e-15）。含 Dirichlet 边界或源项时总量按物理
通量/源强变化，偏差可由边界条件解释；`CheckInvariant` 报告步号、当前值、预期值与偏差，
便于定位。

## 运行

```sh
# 串行/并行求解并保存检查点
go run ./cmd/pdesolve -nx 64 -ny 48 -dt 1e-5 -steps 100 -workers 4 -save ckpt.bin
# 断点续算
go run ./cmd/pdesolve -resume ckpt.bin -steps 100 -workers 2
```

## 验证方法

```sh
go test -race -v ./solver/
```

- `TestInvalidConfigsAreRejected`：9 类非法配置均被拒绝且归类为 ErrKindConfig。
- `TestStabilityLimitEnforced`：超限 dt 归类为 ErrKindStability，与配置错误可区分。
- `TestSerialParallelConsistency`：workers ∈ {1,2,3,4,7,16} 与串行基线逐位一致。
- `TestResumeConsistency`：分段落盘续算与不间断运行逐位一致。
- `TestInvariantConservation`：全 Neumann 无源时总量守恒（偏差 ~1e-15，容差 1e-9）。
- `TestNaNGuard`：绕过前置校验时运行期防护捕获非有限值并报告步号与坐标。

测试日志打印输入参数、判定依据（容差）与实测偏差数值。
