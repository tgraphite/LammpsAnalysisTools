# hbond — 氢键

> 状态：稳定（1.0.0） · 需要拓扑：是 · 需要轨迹：是（≥2 帧）

## 用途

按几何判据识别氢键，输出逐帧数量、存在自相关与寿命、以及氢键中点的空间分布。

## 用法

    lat hbond -s system.data -f run.lammpstrj -sel1 "type 5" -sel2 "type 6" -o out

## 参数

<!-- lat:params:begin -->
| 参数 | 类型 | 默认 | 说明 |
|---|---|---|---|
| `-sel1` | 字符串 | — | donor group (default: all donors) |
| `-sel2` | 字符串 | — | acceptor group (default: all acceptors) |
| `-distance` | 浮点数 | 3.5 | D...A (or H...A) cutoff in Angstrom |
| `-distance_mode` | 字符串 | DA | distance measured D...A (gmx default) or H...A |
| `-angle` | 浮点数 | 30 | maximum D-H...A angle in degrees |
| `-nbins` | 整数 | 50 | bins per axis for the xz/yz/xy geometry maps |
| `-lifetime` | 开关 | true | compute the ACF and lifetime tables |
| `-frame_range` | 字符串 | — | frame window start:stop[:stride] |
| `-stride` | 整数 | — | sample every N frames (default: 1) |
<!-- lat:params:end -->

公共参数（-s、-f、-o、-frame_range、-stride 等）见 [common.md](common.md)。

## 输出

| 文件 | 内容 |
|---|---|
| <base>.hbond.out.csv | frame 与氢键数 |
| <base>.hbond.acf.out.csv | lag、间断/连续自相关、寿命标量 |
| <base>.hbond.dist_xz.out.csv | D–A 中点在 xz 平面的分布 |
| <base>.hbond.dist_yz.out.csv | D–A 中点在 yz 平面的分布 |
| <base>.hbond.dist_xy.out.csv | D–A 中点在 xy 平面的分布 |

所有 CSV 都带统一文件头；标量结果在 # scalars: 行，见 [common.md](common.md#输出约定)。

## 原理

判据与 gmx hbond 对齐：D···A 距离 ≤ distance（默认 3.5 Å，-distance_mode DA；HA 则用 H···A），且在给体 D 处的 H–D···A 夹角 ≤ angle（默认 30°）。只有 N/O/F/S/Cl/Br/I 上的 H 才可能是给体；C–H 不是给体。需要明确的键表来识别 D–H，不猜测。自相关按有效起点无偏归一化，寿命同时给积分法 τ_int=Σ_{τ≥1}C(τ) 与指数拟合 τ_fit=−1/斜率，单位是帧。

## 判读

逐帧计数应随温度/浓度合理变化；ACF 应从 1 衰减；τ_int 与 τ_fit 应同量级。与 gmx 比较时注意两点：lat 的 -sel1 是给体组、-sel2 是受体组（有方向），而 gmx 在两组并集内识别给体/受体；lat 的受体集合包含卤素和硫，含 Cl⁻ 的体系计数会高于 gmx，这是化学上更完整的口径。

## 告警与限制

- 没有键信息直接报错。
- sel1 与 sel2 必须相同或完全不重叠，否则报错。
- sel1 内没有电负性给体时直接报错（与 gmx 的 No Donors found 一致）。
- 少于 2 帧直接报错；ACF 在轨迹后半段仍未衰减时会 WARN，寿命被截断。
- 几何分布使用周期中点，跨边界不会丢点。

## 相关

滞留时间见 [residence.md](residence.md)；π–π 堆积见 [pistack.md](pistack.md)。
