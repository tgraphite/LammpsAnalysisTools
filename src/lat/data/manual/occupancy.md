# occupancy — 区域占据、穿越与首次通过

> 状态：稳定（1.0.0） · 需要拓扑：是 · 需要轨迹：是（≥2 帧）

## 用途

把盒子划分为一个或多个笛卡尔区域，统计每帧占据数、区域间穿越、首次到达时间和富集因子。

## 用法

    lat occupancy -s system.data -f run.lammpstrj -sel "type 2" -zone "z:0-10,10-20 x:0-5" -o out

## 参数

<!-- lat:params:begin -->
| 参数 | 类型 | 默认 | 说明 |
|---|---|---|---|
| `-sel` | 字符串 | 必填 | atoms or molecules to track |
| `-zone` | 字符串 | 必填 | axis:min-max[,min-max] ... (half-open intervals; space-separated axes) |
| `-frame_range` | 字符串 | — | frame window start:stop[:stride] |
| `-stride` | 整数 | — | sample every N frames (default: 1) |
<!-- lat:params:end -->

公共参数（-s、-f、-o、-frame_range、-stride 等）见 [common.md](common.md)。

## 输出

| 文件 | 内容 |
|---|---|
| <base>.occupancy.out.csv | frame 与各区域占据数 |
| <base>.occupancy.crossings.out.csv | frame、from_zone、to_zone、count |
| <base>.occupancy.first_passage.out.csv | 首次到达帧与原子数 |
| <base>.occupancy.enrichment.out.csv | 区域体积分数、占据分数、富集因子 |

所有 CSV 都带统一文件头；标量结果在 # scalars: 行，见 [common.md](common.md#输出约定)。

## 原理

区域语法是 axis:min-max[,min-max] ...，多个轴用空格分隔，区间为半开 [min,max)，多轴区间做笛卡尔积。每帧给每个原子打区域号（0 表示区域外）；穿越统计相邻帧区域号变化；首次通过记录每个原子第一次进入任一区域的帧；富集因子 = 区域占据分数 / 区域体积分数，区域体积按区间与盒子的交集计算。

## 判读

均匀分布时富集因子应接近 1；吸附/排斥会使其显著偏离。穿越计数可用来估计交换速率；首次通过分布可判断初始非平衡态。

## 告警与限制

- zone 必填；区间为空、重叠或轴名错误会直接报错。
- 少于 2 帧时穿越/首次通过没有意义，会 WARN。
- 区域可以超出盒子，超出部分按交集计入体积，但占据只会来自盒内坐标。
- 目前不支持按分子体积自动生成 zone。

## 相关

密度剖面见 [density.md](density.md)。
