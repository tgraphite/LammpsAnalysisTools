# agg — 离子聚合形态（SSIP/CIP/AGG）

> 状态：稳定（1.0.0） · 需要拓扑：是 · 需要轨迹：是（≥1 帧）

## 用途

把中心离子与反离子建成残基级二部图，按连通分量分类为 SSIP（溶剂分离离子对）、CIP（接触离子对）和 AGG（聚集体）。

## 用法

    lat agg -s electrolyte.data -f run.lammpstrj -sel1 "type 2" -sel2 "type 5" -cutoff 3.5 -o out
lat agg -s electrolyte.data -f run.lammpstrj -sel1 "type 2" -sel2 "type 5" -cutoff 3.5 -solv "type 1" -o out

## 参数

<!-- lat:params:begin -->
| 参数 | 类型 | 默认 | 说明 |
|---|---|---|---|
| `-sel1` | 字符串 | 必填 | center group, e.g. 'type 2' for Li+ |
| `-sel2` | 字符串 | 必填 | counter-ion group, e.g. 'type 5' for FSI- |
| `-cutoff` | 浮点数 | 必填 | edge cutoff in Angstrom (required; no default) |
| `-mode` | 字符串 | min | edge criterion: min atom-atom distance or residue COM distance |
| `-solv` | 字符串 | — | optional solvent group; adds the A_iB_jS_k table |
| `-frame_range` | 字符串 | — | frame window start:stop[:stride] |
| `-stride` | 整数 | — | sample every N frames (default: 1) |
<!-- lat:params:end -->

公共参数（-s、-f、-o、-frame_range、-stride 等）见 [common.md](common.md)。

## 输出

| 文件 | 内容 |
|---|---|
| <base>.agg.out.csv | 每帧各类 A_iB_j 团簇数量 |
| <base>.agg.labels.out.csv | A/B 两侧的 SSIP/CIP/AGG 分数 |
| <base>.agg.speciation.out.csv | 仅在 -solv 时输出：A_iB_jS_k 形态计数 |

所有 CSV 都带统一文件头；标量结果在 # scalars: 行，见 [common.md](common.md#输出约定)。

## 原理

对中心离子 A 与反离子 B 按 -mode min|center 和 cutoff 建邻接，求残基级连通分量，记录每个分量的 (i,j) 数量。对中心 A 分类：A₁B₀ → SSIP，A₁B₁ → CIP，其余 → AGG。B 侧的标签按同样规则以 B 为中心计数。使用真图而不是贪心最近邻，避免漏掉桥接离子。可选 -solv 统计每个团簇周围接触的溶剂分子数，输出 S_k 形态。

## 判读

同一物种的 SSIP+CIP+AGG 分数应约等于 1；cutoff 应与 rdf 首壳极小一致；增大 cutoff 会单调地把 SSIP 转化为 CIP/AGG。电荷不平衡是 WARN 而不是错误：两侧电荷不必相等，取决于选择。

## 告警与限制

- cutoff 必填。
- 两侧电荷不平衡会 WARN；全部残基连成一个分量时会 WARN（cutoff 可能过大）。
- 需要拓扑中的分子/残基信息；没有时无法建立残基图。
- 输出以中心 A 和被计数 B 的方向声明，不能把 A 侧与 B 侧分数混用。

## 相关

团簇尺寸分布见 [cluster.md](cluster.md)，配位判据见 [rdf.md](rdf.md)。
