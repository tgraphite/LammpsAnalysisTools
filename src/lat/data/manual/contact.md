# contact — 接触数与距离

> 状态：稳定（1.0.0） · 需要拓扑：是 · 需要轨迹：是（≥1 帧）

## 用途

统计两个原子组之间的接触数、最近距离，以及两组质心距离。用于界面接触、吸附和离子对筛查。

## 用法

    lat contact -s system.data -f run.lammpstrj -sel1 "type 2" -sel2 "type 5" -cutoff 5 -o out
lat contact -s system.data -f run.lammpstrj -sel1 "mol 1" -sel2 "mol 2" -cutoff 6 -mode center -o out

## 参数

<!-- lat:params:begin -->
| 参数 | 类型 | 默认 | 说明 |
|---|---|---|---|
| `-sel1` | 字符串 | 必填 | first group |
| `-sel2` | 字符串 | 必填 | second group |
| `-cutoff` | 浮点数 | 必填 | contact cutoff in Angstrom (required; no default) |
| `-mode` | 字符串 | min | count atom pairs (min) or residue pairs by COM (center) |
| `-frame_range` | 字符串 | — | frame window start:stop[:stride] |
| `-stride` | 整数 | — | sample every N frames (default: 1) |
<!-- lat:params:end -->

公共参数（-s、-f、-o、-frame_range、-stride 等）见 [common.md](common.md)。

## 输出

| 文件 | 内容 |
|---|---|
| <base>.contact.out.csv | frame、接触数、最小距离、质心距离 |

所有 CSV 都带统一文件头；标量结果在 # scalars: 行，见 [common.md](common.md#输出约定)。

## 原理

-mode min（默认）在原子对上计数：距离小于 cutoff 的 A–B 原子对数量，并给出这些对中的最小距离。-mode center 先按残基/分子求质心，再在质心对上计数并给出最小质心距离。两种模式都同时输出两组整体质心的周期最小镜像距离 com_dist。所有距离使用最小镜像。

## 判读

接触数随 cutoff 增大应单调增加；min_dist 应大于 0 且与首壳距离一致；com_dist 变化应与整体接近/远离趋势一致。比较不同体系时固定 cutoff 和选择。

## 告警与限制

- cutoff 必填，没有默认值。
- -mode center 的质心是几何中心而非质量加权，与 gmx 的 com of group 可能有约 0.01 Å 的差异。
- 选区跨盒且不是整分子时，质心会 WARN 并给出歧义提示。
- 两个选择完全相同（或重叠）时，中心模式按无序残基对处理；normal 模式按无序原子对处理。

## 相关

滞留时间见 [residence.md](residence.md)，离子聚类见 [agg.md](agg.md)。
