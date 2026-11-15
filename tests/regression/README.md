# 回归测试

本目录只保留可随仓库分发的自包含回归。

- `level2_format/`：LAMMPS data/lammpstrj 解析、模块冒烟、PBC 平移不变性。
- `samples/sys1_reduced/`：6600 原子的自包含样例（data/lammpstrj/gro/xtc/mapping/selections/verify）。

单元测试在 `../unit/`。历史 GROMACS 对拍脚本与转换器不随 1.0.0 发布；关键数值结论记录在
`docs/CHANGELOG.md` 与各模块手册的“判读”段落。

运行：

    pytest tests/unit tests/regression/level2_format

已知限制：`.gro` 不携带质量/电荷；MDAnalysis 读 `.lammpstrj` 必须显式 `format="LAMMPSDUMP"`。
