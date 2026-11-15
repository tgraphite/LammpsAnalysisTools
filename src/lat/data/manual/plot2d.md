# plot2d — 二维 CSV 绘图

> 状态：稳定（1.0.0） · 需要拓扑：否 · 需要轨迹：否

## 用途

读取 lat 的二维 CSV（长表 x,y,value 或网格），输出热图/等高线图。

## 用法

    lat plot2d -f out.density.2d.out.csv -o density2d.png
    lat plot2d -f out.hbond.dist_xy.out.csv -smooth gaussian -o hbond_xy.png

## 参数

| 参数 | 类型 | 默认 | 说明 |
|---|---|---|---|
| -f | 字符串 | 必填 | 输入 CSV（二维长表） |
| -o | 字符串 | 与输入同名 .png | 输出图片路径 |
| -title | 字符串 | CSV 文件头 | 覆盖标题 |
| -smooth | 字符串 | config plot.smooth（none） | none、gaussian、savgol；只影响绘图 |

公共参数见 [common.md](common.md)。

## 输出

图片文件（png/pdf/svg 等）。

## 原理

从文件头读取列名，把前两列解释为 x/y，第三列解释为值；若 CSV 是长表则重排为网格后绘制。平滑在绘图阶段执行，CSV 保持原始值。

## 判读

二维图用于观察界面、孔道、吸附层和空间分布。颜色范围应结合物理量单位解释；平滑后不要从图上读取精确极值。

## 告警与限制

- CSV 缺少文件头或不是二维长表会报错。
- 不支持 -columns/-bar/-stack。
- 平滑窗口/阶数来自配置，过大窗口会抹掉真实界面。

## 相关

一维绘图见 [plot1d.md](plot1d.md)。