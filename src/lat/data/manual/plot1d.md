# plot1d — 一维 CSV 绘图

> 状态：稳定（1.0.0） · 需要拓扑：否 · 需要轨迹：否

## 用途

读取 lat 生成的 CSV，把第一列作为 x 轴、其余列作为 y 轴，输出折线图、柱状图或堆积柱状图。计算与绘图严格分离：绘图不会改变 CSV。

## 用法

    lat plot1d -f out.rdf.out.csv -o rdf.png
    lat plot1d -f out.agg.labels.out.csv -columns A_SSIP,A_CIP,A_AGG -o spec.png
    lat plot1d -f out.agg.labels.out.csv -stack -o spec_stack.png

## 参数

| 参数 | 类型 | 默认 | 说明 |
|---|---|---|---|
| -f | 字符串 | 必填 | 输入 CSV（文件头携带 title/xlabel/ylabel/columns） |
| -o | 字符串 | 与输入同名 .png | 输出图片路径 |
| -title | 字符串 | CSV 文件头 | 覆盖标题 |
| -smooth | 字符串 | config plot.smooth（none） | none、gaussian、savgol；只影响绘图 |
| -columns | 列表 | 全部 y 列 | 只画这些列；逗号或空格分隔，x 轴恒为第一列 |
| -bar | 开关 | false | 分组柱状图，每个 y 列一组 |
| -stack | 开关 | false | 归一化堆积柱状图，适合占据/分数列 |

公共参数见 [common.md](common.md)。

## 输出

图片文件（扩展名决定格式，常用 png/pdf/svg）。屏幕只输出 wrote 路径。

## 原理

读取 CSV 文件头的 # columns、# title、# xlabel、# ylabel，按列名解析数据；-columns 支持逗号/空格分隔并保持用户给出的顺序。-bar 与 -stack 互斥，且不做平滑。-stack 会把每个 x 位置的 y 值归一化为分数，适合画 SSIP/CIP/AGG 或区域占据。

## 判读

曲线形状应与对应模块手册的判读一致；柱状/堆积图适合少量离散类别，不适合连续分布。平滑只用于改善视觉效果，定量结论以 CSV 原始值为准。

## 告警与限制

- CSV 缺少文件头或列名不存在会直接报错。
- -bar 与 -stack 互斥；-bar/-stack 不接受 -smooth。
- plot2d 不接受 -columns/-bar/-stack。
- 需要 matplotlib；无图形界面时使用非交互后端即可。

## 相关

二维绘图见 [plot2d.md](plot2d.md)；输出契约见 [common.md](common.md#输出约定)。