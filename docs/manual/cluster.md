# cluster — 配位数与团簇尺寸分布

> 状态：稳定（1.0.0） · 需要拓扑：是 · 需要轨迹：是（≥1 帧）

## 用途

在残基图上统计中心离子的配位数分布，以及连通团簇的尺寸分布。

## 用法

    lat cluster -s electrolyte.data -f run.lammpstrj -sel1 "type 2" -sel2 "type 5" -cutoff 3.5 -o out
lat cluster -s electrolyte.data -f run.lammpstrj -sel1 "type 2" -sel2 "type 5" -cutoff 3.5 -focus -o out

## 参数

<!-- lat:params:begin -->
| 参数 | 类型 | 默认 | 说明 |
|---|---|---|---|
| `-sel1` | 字符串 | 必填 | center group |
| `-sel2` | 字符串 | 必填 | counter or neighbour group |
| `-cutoff` | 浮点数 | 必填 | edge cutoff in Angstrom (required; no default) |
| `-mode` | 字符串 | min | edge criterion: min atom-atom distance or residue COM distance |
| `-focus` | 开关 | false | only per-center coordination numbers; skip cluster sizes |
| `-frame_range` | 字符串 | — | frame window start:stop[:stride] |
| `-stride` | 整数 | — | sample every N frames (default: 1) |
<!-- lat:params:end -->

公共参数（-s、-f、-o、-frame_range、-stride 等）见 [common.md](common.md)。

## 输出

| 文件 | 内容 |
|---|---|
| <base>.cluster.cn.out.csv | 配位数、计数、概率 |
| <base>.cluster.size.out.csv | 团簇尺寸、计数、概率；-focus 时不输出 |

所有 CSV 都带统一文件头；标量结果在 # scalars: 行，见 [common.md](common.md#输出约定)。

## 原理

与 agg 共用残基图引擎：按 -mode 和 cutoff 建边，对每个中心 A 统计连接的 B 数，并把所有帧的计数累加为分布；非 -focus 时对每个连通分量统计总残基数，得到尺寸分布。概率按全部帧的中心数或团簇数归一。

## 判读

平均配位数应与 rdf 首壳配位数一致；尺寸分布应随浓度/温度移动。对比不同 cutoff 时关注分布整体平移，而不是单个 bin 的噪声。

## 告警与限制

- cutoff 必填。
- -focus 只输出配位数分布，不建团簇尺寸表，速度更快。
- 与 agg 一样需要分子/残基信息。
- 所有帧的计数直接累加，帧数不等时概率仍按总样本归一。

## 相关

形态分类见 [agg.md](agg.md)。
