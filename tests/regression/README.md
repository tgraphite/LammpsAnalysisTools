# 回归测试

两级输入回归（详见 `docs/05_regression.md`）：

- `level1_algo/`：原生 GROMACS 输入（tpr/gro/xtc，MDAnalysis 内部直读）→ 与 gmx `.xvg` 比对，验证**算法**。
- `level2_format/`：转成 LAMMPS `data`+`lammpstrj` 后 → 与同一份结果比对，验证**格式解析 + 算法**。
- `convert/`：**已实现**的自研转换器（不依赖 InterMol/VMD）：
  - `convert_gmx_to_lammps.py`：gmx → LAMMPS data + lammpstrj（可缩减分子/帧）。
  - `make_selections.py`：生成元素/分子/模块候选选择。
  - `verify_lammps_roundtrip.py`：坐标/盒/键/元素/选择一致性自检。
- `samples/sys1_reduced/`：**已落盘自包含样例**（data/lammpstrj/gro/xtc/mapping/selections/verify），
  自检全通过（坐标偏差 1e-6 Å，元素与选择一致）。见样例 `README.md`。
- `expected/`：从 `/mnt/f/CRO/40_multi_electrolytes/11-systems` 复制的精简期望值。

覆盖：sys1–sys6（sys4/sys6 含 PMBA 聚合物 + 电解液）。
- rdf/cn、density、msd 有 GT；agg/cluster 做自洽 + 转换一致性。
- hbond 无 GT → 用 `tests/data/` 自造对象，手算点数。
- pistack 沿用 `/mnt/d/dev/83_pipi_stack_code/tests`。

已知限制：`.gro` 不携带质量/电荷（Li 质量会被猜为 0）；权威拓扑见 `*.mapping.csv`。
MDAnalysis 读 `.lammpstrj` 必须显式 `format="LAMMPSDUMP"`。

不建：固体/金属力学、能量、缺陷、SDF、谱。
