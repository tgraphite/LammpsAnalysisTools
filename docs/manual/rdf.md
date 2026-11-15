# rdf — 径向分布函数与配位数

> 状态：稳定（1.0.0） · 需要拓扑：是 · 需要轨迹：是（≥1 帧）

## 用途

计算两个原子组之间的径向分布函数 g(r) 与积分配位数 n(r)。这是判断液体结构、离子配位和首壳距离的入口模块。

## 用法

    lat rdf -s system.data -f run.lammpstrj -sel1 "type 2" -sel2 "type 5" -cutoff 10 -o out

## 参数

<!-- lat:params:begin -->
| 参数 | 类型 | 默认 | 说明 |
|---|---|---|---|
| `-sel1` | 字符串 | 必填 | first selection (A) |
| `-sel2` | 字符串 | 必填 | second selection (B) |
| `-nbins` | 整数 | 200 | number of radial bins |
| `-cutoff` | 浮点数 | 10.0 | maximum radius in Angstrom (keep below half the box edge) |
| `-frame_range` | 字符串 | — | frame window start:stop[:stride] |
| `-stride` | 整数 | — | sample every N frames (default: 1) |
| `-intramol` | 开关 | false | include intramolecular (same-molecule) pairs; default: excluded |
<!-- lat:params:end -->

公共参数（-s、-f、-o、-frame_range、-stride 等）见 [common.md](common.md)。

## 输出

| 文件 | 内容 |
|---|---|
| <base>.rdf.out.csv | r 与 g(r) |
| <base>.rdf.cn.out.csv | r 与积分配位数 n(r) |

所有 CSV 都带统一文件头；标量结果在 # scalars: 行，见 [common.md](common.md#输出约定)。

## 原理

逐帧统计 A–B 原子对的距离直方图。g(r) 按逐帧盒子体积加权，使用精确球壳体积 (4/3)π(r_hi³−r_lo³)；配位数由同一组直方图累加，因此 g 与 n 在恒定体积下严格自洽。sel1 与 sel2 完全相同时按无序对还原有序对（乘 2）。默认排除同一分子内的原子对；-intramol 打开后包含。横坐标是 bin 中心 (i+½)Δr，而 gmx 报告的是左边缘 iΔr，逐点比较时需注意半 bin 偏移。

## 判读

均匀液体在长程回到 g≈1；首峰位置对应典型接触距离；n(r) 单调上升并在首壳极小后进入平台。判断正确性优先看积分量（首壳配位数、平台高度），不要用逐点 max|Δg|：有限采样下逐点涨落可达数个标准误，而配位数的相对误差通常小一个量级。

## 告警与限制

- cutoff 应小于最短盒边的一半，否则会出现 WARN 且 g(r) 失真。
- 拓扑有键但没有分子 ID 时无法排除同分子对，会 WARN；需要分子内 RDF 时显式加 -intramol。
- NPT 下 g 按逐帧 V 加权、n 不含体积因子，两者是同一口径下的近似自洽。
- 少于 2 帧会 WARN，结果只适合冒烟。

## 相关

算法原理与误差来源见 [MODULE_STATUS.md](../MODULE_STATUS.md)；gmx 对拍方法见仓库 tests/regression。
