# msd — 均方位移与扩散系数

> 状态：稳定（1.0.0） · 需要拓扑：是 · 需要轨迹：是（≥3 帧）

## 用途

计算均方位移 MSD(τ) 并用 Einstein 关系拟合扩散系数 D。支持原子坐标或分子质心。

## 用法

    lat msd -s system.data -f run.lammpstrj -sel "type 2" -dim 3 -o out

## 参数

<!-- lat:params:begin -->
| 参数 | 类型 | 默认 | 说明 |
|---|---|---|---|
| `-sel` | 字符串 | 必填 | diffusing atoms or molecules |
| `-dim` | 整数 | 3 | 2D or 3D mean squared displacement |
| `-seltype` | 字符串 | atom | atom positions or per-molecule COM |
| `-fit_range` | 字符串 | — | lag fit window start:stop (default: 10%-50% of lags) |
| `-remove_com` | 开关 | true | subtract the selection COM motion (default: on) |
| `-noremove_com` | 开关 | — | keep the COM drift (overrides -remove_com) |
| `-stride` | 整数 | — | sample every N frames (default: 1) |
| `-frame_range` | 字符串 | — | frame window start:stop[:stride] |
<!-- lat:params:end -->

公共参数（-s、-f、-o、-frame_range、-stride 等）见 [common.md](common.md)。

## 输出

| 文件 | 内容 |
|---|---|
| <base>.msd.out.csv | lag、总 MSD 与各分量 MSD |

所有 CSV 都带统一文件头；标量结果在 # scalars: 行，见 [common.md](common.md#输出约定)。

## 原理

MSD(τ)=⟨|r(t+τ)−r(t)|²⟩，对全部时间原点平均，并对每个 lag 取点。优先使用 xu/yu/zu 或 x+i·L 解包；两者都没有时退回 wrapped 坐标并 WARN。默认扣除所选组的整体质心漂移，-noremove_com 关闭。D 由 fit_range 内的 MSD 对 lag 做最小二乘直线得到，斜率除以 6（3D）或 4（2D）；默认拟合窗是 10%–50% 的 lag 范围。D 默认单位 Å²/frame；配置 meta.frame_interval 后同时给出 D_cm2_per_s。

## 判读

长轨迹上 MSD 对 lag 近似线性，R² 接近 1。D 应与温度、粘度趋势一致。短轨迹（几帧）MSD 可能非单调、R² 无意义，这是有限采样而非实现错误。

## 告警与限制

- 少于 3 帧直接报错。
- 轨迹没有 image flags 或解包坐标时可能饱和，结果只作参考。
- 位移接近半盒时 PBC 解包本身有歧义，任何实现都会有几 % 偏差。
- R² < 0.8 会 WARN，此时 D 不可靠。

## 相关

输运系数（电导、迁移数）见 [transport.md](transport.md)。
