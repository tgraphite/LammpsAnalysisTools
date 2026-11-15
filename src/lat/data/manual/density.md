# density — 密度剖面

> 状态：稳定（1.0.0） · 需要拓扑：是 · 需要轨迹：是（≥1 帧）

## 用途

沿一个坐标轴输出一维密度剖面，或在垂直于该轴的平面内输出二维密度网格。支持质量密度与数密度。

## 用法

    lat density -s system.data -f run.lammpstrj -axis z -nbins 200 -quantity mass -o out
lat density -s system.data -f run.lammpstrj -dim 2 -axis z -o out

## 参数

<!-- lat:params:begin -->
| 参数 | 类型 | 默认 | 说明 |
|---|---|---|---|
| `-sel` | 字符串 | — | selection (default: all atoms) |
| `-axis` | 字符串 | z | axis normal to the profile |
| `-nbins` | 整数 | 200 | bins along the axis (bins per side in 2D) |
| `-dim` | 整数 | 1 | 1 = profile along -axis, 2 = plane spanned by the other two axes |
| `-quantity` | 字符串 | mass | mass density (kg/m^3) or number density (1/A^3) |
| `-frame_range` | 字符串 | — | frame window start:stop[:stride] |
| `-stride` | 整数 | — | sample every N frames (default: 1) |
<!-- lat:params:end -->

公共参数（-s、-f、-o、-frame_range、-stride 等）见 [common.md](common.md)。

## 输出

| 文件 | 内容 |
|---|---|
| <base>.density.out.csv | 1D：坐标与密度 |
| <base>.density.2d.out.csv | 2D 长表：x、y、密度 |

所有 CSV 都带统一文件头；标量结果在 # scalars: 行，见 [common.md](common.md#输出约定)。

## 原理

采用 gmx density 的分数坐标分箱：每帧把坐标按该帧盒长归一化到 [0,1)，落入固定的 nbins 个分数薄层，因此 NPT 下 bin 随盒缩放，原子不会掉出边界；最后用平均盒长标定坐标轴。薄层体积为 V_f/nbins。质量密度换算为 kg/m³；数密度单位为 1/Å³。平滑只在 plot1d/plot2d 做，CSV 始终是原始值。

## 判读

均匀液体应给出近似水平的曲线；出现分层、界面峰或吸附层时应与体系预期一致。2D 图适合观察孔道、表面偏析和相分离。比较不同工具时优先看平均值和积分，而非界面处的单点极值。

## 告警与限制

- 质量密度需要每个 type 的质量；没有 Masses 时改用 -quantity number。
- 任一相关轴的帧盒长相对均值偏离超过 ±2% 会 WARN；此时轴按平均盒长标定。
- 三斜或畸变盒的截面归一化是近似，复杂盒请自行确认。

## 相关

绘图见 [plot1d.md](plot1d.md) 与 [plot2d.md](plot2d.md)。
