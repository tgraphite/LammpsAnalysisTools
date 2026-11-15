# describe — 模块清单与契约

> 状态：稳定（1.0.0） · 需要拓扑：否 · 需要轨迹：否

## 用途

列出全部模块，或显示一个模块的输入、输出、参数来源和帮助手册路径。

## 用法

    lat describe
lat describe -module rdf

## 参数

<!-- lat:params:begin -->
| 参数 | 类型 | 默认 | 说明 |
|---|---|---|---|
| `-module` | 字符串 | — | module name to describe |
<!-- lat:params:end -->

公共参数（-s、-f、-o、-frame_range、-stride 等）见 [common.md](common.md)。

## 输出

| 文件 | 内容 |
|---|---|
| （屏幕） | 模块清单或单模块契约 |

所有 CSV 都带统一文件头；标量结果在 # scalars: 行，见 [common.md](common.md#输出约定)。

## 原理

直接从模块注册表读取元数据：summary、inputs、outputs、每个参数的 flags/默认来源/help，并在末尾给出手册路径。参数默认值的实际解析顺序是命令行 > 配置文件；describe 只显示定义，不执行计算。

## 判读

当不确定某个模块需要哪些参数时先运行 describe -module <name>，再用 lat <name> -h 查看 argparse 生成的完整帮助和默认值。

## 告警与限制

- describe 不接受 -s/-f，也不需要拓扑。
- 手册路径在源码树中是 docs/manual/，安装后是包内的 data/manual/。

## 相关

模块总览见 [README.md](README.md)。
