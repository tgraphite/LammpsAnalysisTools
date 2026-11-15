# transport — 扩散、Nernst–Einstein 与 Green–Kubo 输运

> 状态：稳定（1.0.0） · 需要拓扑：是 · 需要轨迹：是（≥3 帧，且含速度列）

## 用途

按物种计算扩散系数、Nernst–Einstein 电导与迁移数，并输出 Green–Kubo 运行电导积分。

## 用法

    lat transport -s electrolyte.data -f run.lammpstrj -sel "type 2 5" -temperature 300 -velocity_unit A/ps -o out

## 参数

<!-- lat:params:begin -->
| 参数 | 类型 | 默认 | 说明 |
|---|---|---|---|
| `-sel` | 字符串 | 必填 | charged species to analyse |
| `-temperature` | 浮点数 | — | temperature in K; required for NE/GK |
| `-velocity_unit` | 字符串 | A/ps | trajectory velocity unit: A/ps or A/fs |
| `-frame_range` | 字符串 | — | frame window start:stop[:stride] |
| `-stride` | 整数 | — | sample every N frames (default: 1) |
<!-- lat:params:end -->

公共参数（-s、-f、-o、-frame_range、-stride 等）见 [common.md](common.md)。

## 输出

| 文件 | 内容 |
|---|---|
| <base>.transport.d.out.csv | 各 type 的 D、D_cm2_s、R² |
| <base>.transport.ne.out.csv | Nernst–Einstein 电导与迁移数 |
| <base>.transport.gk.out.csv | Green–Kubo 运行电导 |

所有 CSV 都带统一文件头；标量结果在 # scalars: 行，见 [common.md](common.md#输出约定)。

## 原理

前置条件是轨迹含 vx/vy/vz，且拓扑含电荷。各物种 D_i 由 MSD 斜率除以 6 得到；σ_NE=(F²/RT)Σ z_i²c_iD_i，迁移数 t_i=z_i²c_iD_i/Σ_j z_j²c_jD_j。Green–Kubo 用电流 J(t)=Σ q_i v_i，σ_GK=(1/3Vk_BT)∫⟨J(0)·J(τ)⟩dτ，输出运行积分。速度单位由 -velocity_unit 指定（A/ps 或 A/fs）。D 由 Å²/frame 换算到 cm²/s 或 m²/s 需要 frame_interval。

## 判读

D 应随温度升高、随体系粘度降低；σ_NE 通常高于 σ_GK，因为它忽略离子间关联；迁移数之和应约等于 1。比较不同轨迹时使用相同的拟合窗口和帧间隔。

## 告警与限制

- 没有速度列直接报错。
- 没有电荷、没有 -temperature 或没有 frame_interval 时跳过 NE/GK 并 WARN，只输出 D 表。
- 少于 3 帧直接报错。
- NE 忽略离子关联，通常高估电导；报告时必须声明口径。
- 运行积分是否收敛应通过观察 gk 表尾部平台判断。

## 相关

MSD 细节见 [msd.md](msd.md)。
