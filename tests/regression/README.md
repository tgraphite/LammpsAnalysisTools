# 回归测试

两级输入回归（详见 `docs/05_regression.md`）：

- `level1_algo/`：原生 GROMACS 输入（tpr/gro/xtc，MDAnalysis 内部直读）→ 与 gmx `.xvg` 比对，验证**算法**。
  - `compare_gmx.py`：sys1–sys6 的 rdf/cn/density/msd 对拍。
  - `gmx_ref_lammps.py`：**LAMMPS 输入**的 gmx 原生对拍（gyrate/rms/rmsf/sasa/rdf/hbond），
    需要 `convert/lammps_to_gmx.py` 先转拓扑。见 `level1_algo/README.md`。
- `level2_format/`：转成 LAMMPS `data`+`lammpstrj` 后 → 与同一份结果比对，验证**格式解析 + 算法**。
- `convert/`：**已实现**的自研转换器（不依赖 InterMol/VMD）：
  - `convert_gmx_to_lammps.py`：gmx → LAMMPS data + lammpstrj（可缩减分子/帧）。
  - `lammps_to_gmx.py`：**反向**，LAMMPS data + dump → gmx `top/gro/ndx/tpr/traj.xtc/whole.xtc`，
    专为「用原生 gmx 验证 LAMMPS 路径」而建（按 id 顺序建拓扑、每个分子一个 moleculetype、
    σ 由 lat 的 vdW 表反解）。
  - `make_selections.py`：生成元素/分子/模块候选选择。
  - `verify_lammps_roundtrip.py`：坐标/盒/键/元素/选择一致性自检。
- `samples/sys1_reduced/`：**已落盘自包含样例**（data/lammpstrj/gro/xtc/mapping/selections/verify），
  自检全通过（坐标偏差 1e-6 Å，元素与选择一致）。见样例 `README.md`。
- `expected/`：从 `/mnt/f/CRO/40_multi_electrolytes/11-systems` 复制的精简期望值。

覆盖：sys1–sys6（sys4/sys6 含 PMBA 聚合物 + 电解液）。
- rdf/cn、density、msd 有 GT；agg/cluster 做自洽 + 转换一致性。
- **hbond 已有 GT**：`gmx hbond`（PET 体系 501 帧，逐帧计数差 0.014、最大 1 条键）。
  `sys1_reduced` 只有 C–H、没有电负性给体，gmx 与 lat 都直接拒绝，因此样本级 hbond 测试
  改为「无给体 → 报错」+ 合成 O–H 空帧几何测试。
- pistack 沿用 `/mnt/d/dev/83_pipi_stack_code/tests`。
- **LAMMPS 路径已有 gmx GT**：gyrate/Rg 3e-4 Å、RMSD 1e-4 Å、RMSF 7e-4 Å、CN 0.26%、hbond 0.21%；
  SASA 确认有 14% 的 PBC 缺陷（已记录未修）。

已知限制：`.gro` 不携带质量/电荷（Li 质量会被猜为 0）；权威拓扑见 `*.mapping.csv`。
MDAnalysis 读 `.lammpstrj` 必须显式 `format="LAMMPSDUMP"`。

不建：固体/金属力学、能量、缺陷、SDF、谱。
