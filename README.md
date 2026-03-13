# LNB-MDT

![LNB-MDT Logo](LNB-MDT.jpg)

**LNB-MDT** (Lipid NanoBubble Molecular Dynamics Toolkit) is a comprehensive toolkit designed for molecular dynamics simulations of lipid nanobubbles. It also supports **protein RMSD analysis** (including ligand and complex) for MD trajectories, with optional trajectory alignment and plotting.

## Installation

### Method 1: Install from PyPI (Recommended)

The easiest way to install LNB-MDT is using pip. We recommend using conda to create a virtual environment first:

```bash
# Create conda environment
conda create -n LNB-MDT python=3.11
conda activate LNB-MDT

# Install from PyPI
pip install lnb-mdt
```

### Method 2: Install from Source

If you want to install the latest development version or contribute to the project:

```bash
# Clone the repository
git clone https://github.com/xinyuren-bio/LNB-MDT.git
cd LNB-MDT

# Create conda environment (optional but recommended)
conda create -n LNB-MDT python=3.11
conda activate LNB-MDT

# Install in editable mode
pip install -e .
```

**Note:** For editable installation, you need Python 3.10+ and pip. Using conda is recommended for managing dependencies.

### Verify Installation

After installation, verify that LNB-MDT is correctly installed:

```bash
# Check if command is available
LNB-MDT --help

# Or test in Python
python -c "import LNB_MDT; print('LNB-MDT installed successfully!')"
```

## Quick Start

After installation, you can use LNB-MDT in two ways:

### 1. Command Line Interface (CLI)

```bash
# Launch GUI
LNB-MDT UI

# Run area analysis
LNB-MDT AREA --help

# Run with test data
LNB-MDT AREA -test

# Configure VMD path
LNB-MDT VMD --path /path/to/vmd

# Run protein RMSD analysis (see below for details)
LNB-MDT RMSD --help
```

#### Protein RMSD Analysis

LNB-MDT 支持对蛋白质（及配体/复合物）轨迹进行 RMSD 分析，基于 MDAnalysis 实现轨迹对齐与多组 RMSD 计算。

**基本用法：**

```bash
# 仅蛋白质 RMSD，结果输出到指定目录
LNB-MDT RMSD -input_gro /path/to/complex.pdb -input_xtc /path/to/fit.xtc -output_dir ./results

# 指定蛋白与配体选择，并生成 RMSD 曲线图
LNB-MDT RMSD -input_gro ./complex.pdb -input_xtc ./fit.xtc -protein_name "protein" -lig_name "resname LIG" -plot

# 使用 .gro 拓扑
LNB-MDT RMSD -input_gro ./system.gro -input_xtc ./traj.xtc -output_dir ./results -plot
```

**常用参数：**

| 参数 | 简写 | 说明 | 默认值 |
|------|------|------|--------|
| `-input_gro` | `-g` / `-i` | 拓扑文件路径（.pdb 或 .gro） | 必填 |
| `-input_xtc` | `-x` / `-t` | 轨迹文件路径（如 .xtc） | 必填 |
| `-output_dir` | `-o` | 结果输出目录 | `./results` |
| `-protein_name` | `-p` | 蛋白质选择语句 | `protein` |
| `-lig_name` | `-l` | 配体选择语句（如 `resname LIG`），不填则只算蛋白 | 空 |
| `-ref_frame` | `-ref` | 参考帧索引 | `0` |
| `-plot` | — | 计算完成后绘制 RMSD 曲线并保存 PNG | 否 |
| `-no_align` | — | 不做轨迹对齐，直接计算 RMSD | 否（默认会先对齐） |

结果文件：`<output_dir>/rmsd.csv`（时间与各组 RMSD）；若使用 `-plot`，还会生成 `rmsd.png`。

### 2. Python API

```python
from LNB_MDT.analysis import Area
import MDAnalysis as mda

# Load trajectory
u = mda.Universe("system.gro", "trajectory.xtc")

# Run analysis
area_analysis = Area(u, {'DPPC': ['PO4'], 'CHOL': ['ROH']})
area_analysis.run()
```

## Documentation

For detailed documentation, including installation guide, quick start, user guide, and command line tools, please visit:

📚 **[Read the Docs - LNB-MDT Documentation](https://lnb-mdt.readthedocs.io/en/latest/)**

## File Structure

```
LNB-MDT/
├── main.py                 # Main program entry
├── requirements.txt        # Python dependencies
├── analysis/              # Analysis modules
│   ├── area.py           # Area analysis
│   ├── height.py         # Height analysis
│   ├── cluster.py        # Cluster analysis
│   ├── anisotropy.py     # Anisotropy analysis
│   ├── gyration.py       # Gyration analysis
│   ├── rmsd.py           # Protein/ligand/complex RMSD analysis
│   ├── sz.py             # Sz order parameter analysis
│   └── density.py        # Density analysis (time and radius)
├── preparation/            # Preparation module
└── cases_lnb/             # Example lipid nanobubble data
    ├── lnb.gro           # Example topology file (Martini 3.0, DPPC:DAPC:CHOL=5:3:2)
    ├── lnb.xtc           # Example trajectory file (50-60 ns time window)
    └── README.md         # Example data description
```

## Citation

If you use LNB-MDT in your research, please cite our paper:

```bibtex
@article{ren2025lnb,
   author = {Xinyu Ren and Xubo Lin},
   title = {LNB-MDT: An Integrated Python Toolkit for Preparing and Analyzing Lipid Nanobubble Simulations},
   journal = {Journal of Chemical Information and Modeling},
   volume = {66},
   number = {1},
   year = {2025},
   month = {December},
   doi = {10.1021/acs.jcim.5c02771},
   URL = {https://pubs.acs.org/doi/10.1021/acs.jcim.5c02771},
}
```

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---

**LNB-MDT** - Making lipid nanobubble simulations simpler and more efficient!
