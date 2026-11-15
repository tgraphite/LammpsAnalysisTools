# sasa — 溶剂可及表面积

> 状态：稳定（1.0.0） · 需要拓扑：是 · 需要轨迹：是（≥1 帧）

## 用途

用 Shrake–Rupley 算法计算溶剂可及表面积，并按选择单独报告某部分原子的贡献。

## 用法

    lat sasa -s system.data -f run.lammpstrj -o out
lat sasa -s system.data -f run.lammpstrj -sel "type 2" -radius 1.445 -sasa_points 92 -o out

## 参数

<!-- lat:params:begin -->
| 参数 | 类型 | 默认 | 说明 |
|---|---|---|---|
| `-radius` | 浮点数 | 1.445 | probe radius in Angstrom |
| `-sasa_points` | 整数 | 92 | sphere sample points per atom |
| `-sel` | 字符串 | — | atoms whose area is reported (surface uses all) |
| `-frame_range` | 字符串 | — | frame window start:stop[:stride] |
| `-stride` | 整数 | — | sample every N frames (default: 1) |
<!-- lat:params:end -->

公共参数（-s、-f、-o、-frame_range、-stride 等）见 [common.md](common.md)。

## 输出

| 文件 | 内容 |
|---|---|
| <base>.sasa.out.csv | frame、总 SASA、选区 SASA |

所有 CSV 都带统一文件头；标量结果在 # scalars: 行，见 [common.md](common.md#输出约定)。

## 原理

每个原子以 r_vdW+r_probe 为半径撒 M 个 Fibonacci 球面点（默认 92），落在任何邻居球内即被遮挡，SASA 为 4πr² 乘以暴露点比例。邻域搜索和遮挡判定全部使用最小镜像，贴面原子与盒中心原子等价。表面始终用全部原子构建；-sel 只决定单独报告哪部分原子的贡献，因此 sasa_sel 是该基团在真实环境中的暴露面积。

## 判读

总 SASA 应随界面/孔隙变化；比较不同探针半径时趋势应一致。与 gmx sasa 对比时，残余差异主要来自球面采样点分布（双方都用 92 点时约 0.6%）。

## 告警与限制

- 需要能推断元素，否则报错。
- 选区为空直接报错。
- 探针半径和采样点数必须与对比工具一致，否则差异不可解释。
- 球面采样是有限离散，增大 -sasa_points 可检查收敛。

## 相关

自由体积与孔径见 [freevol.md](freevol.md)。
