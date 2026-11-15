# freevol — 自由体积、孔径与连通性

> 状态：稳定（1.0.0） · 需要拓扑：是 · 需要轨迹：是（≥1 帧）

## 用途

用周期体素网格计算自由体积分数（FFV）、孔径分布（PSD）和自由空间连通性。

## 用法

    lat freevol -s system.data -f run.lammpstrj -mode ffv -radius 1.445 -grid_dr 0.5 -o out
lat freevol -s polymer.data -f run.lammpstrj -mode psd -grid_dr 0.5 -o out
lat freevol -s polymer.data -f run.lammpstrj -mode connectivity -connectivity 26 -o out

## 参数

<!-- lat:params:begin -->
| 参数 | 类型 | 默认 | 说明 |
|---|---|---|---|
| `-radius` | 浮点数 | 1.445 | probe radius in Angstrom |
| `-mode` | 字符串 | ffv | ffv = free volume fraction; psd = pore size; connectivity = free-space graph |
| `-grid_dr` | 浮点数 | 0.5 | requested voxel edge in Angstrom (may be coarsened) |
| `-connectivity` | 整数 | 26 | 6 = faces, 26 = faces+edges+corners |
| `-psd_rmax` | 浮点数 | — | upper PSD radius (A); default: data-driven |
| `-components` | 开关 | false | also write the per-component connectivity table |
| `-frame_range` | 字符串 | — | frame window start:stop[:stride] |
| `-stride` | 整数 | — | sample every N frames (default: 1) |
<!-- lat:params:end -->

公共参数（-s、-f、-o、-frame_range、-stride 等）见 [common.md](common.md)。

## 输出

| 文件 | 内容 |
|---|---|
| <base>.freevol.out.csv | ffv 模式：frame 与自由体积分数 |
| <base>.freevol.psd.out.csv | psd 模式：孔半径与 dV/dr |
| <base>.freevol.connectivity.out.csv | connectivity 模式：逐帧分量汇总 |
| <base>.freevol.connectivity.components.out.csv | 加 -components 时的逐分量表 |

所有 CSV 都带统一文件头；标量结果在 # scalars: 行，见 [common.md](common.md#输出约定)。

## 原理

每个格点的占用判据是格点中心落在任一原子的 (r_vdW+r_probe) 球内，各轴按最小镜像回绕。ffv 输出空闲格点比例；connectivity 对自由格点做环面闭合的连通分量，默认 26 邻接，输出分量数、贯通分量数、最大分量体积/占比和贯通体积；psd 对 vdW 占据做周期距离变换，在探针可及区域内取局部极大为孔心，用 26 连通平台去重，再按周期最近孔心做 Voronoi 分配，得到归一化 dV/dr（积分 = 1）。

## 判读

FFV 应与已知孔隙率趋势一致；PSD 的主峰应随链段堆积/温度移动；connectivity 的 frac_largest 接近 1 表示自由空间基本贯通，出现大量孤立小孔时该值下降。网格收敛性应通过减小 grid_dr 检查。

## 告警与限制

- 网格分析要求正交且固定的盒子：三斜盒直接报错；所选帧内盒长漂移 >2% 报错，0.1%–2% 给 WARN 且网格锁定首帧。
- 总格点数超过上限时会自动放大 grid_dr 并 WARN，实际网格可能比请求值粗。
- 需要能推断元素（Masses、dump 的 element 列或 -types），否则报错。
- PSD/connectivity 是体素网格近似，不是 Zeo++ 的 Voronoi 几何算法。

## 相关

表面积见 [sasa.md](sasa.md)。
