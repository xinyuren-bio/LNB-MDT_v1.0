"""
RMSD 分析模块：使用 MDAnalysis 对蛋白质/配体/复合物轨迹计算均方根偏差（RMSD）。
支持轨迹对齐、多组选择（蛋白、配体、复合物）及可选绘图。
"""
import os
import sys
import argparse
import warnings
warnings.filterwarnings('ignore')

import numpy as np
import pandas as pd
import MDAnalysis as mda
from MDAnalysis.analysis import rms, align

if __name__ == '__main__':
    current_file_path = os.path.abspath(__file__)
    current_dir = os.path.dirname(current_file_path)
    package_root = os.path.abspath(os.path.join(current_dir, '..'))
    if package_root not in sys.path:
        sys.path.insert(0, package_root)

__all__ = ['run_rmsd_analysis']


def run_rmsd_analysis(
    topology_file,
    trajectory_file,
    output_dir,
    selection_protein="protein",
    selection_ligand="",
    selection_backbone="protein and backbone",
    ref_frame=0,
    start_frame=None,
    do_align=True,
    verbose=True,
):
    """
    执行 RMSD 分析：轨迹对齐（可选）、计算 RMSD、写出 CSV。

    Parameters
    ----------
    topology_file : str
        拓扑文件路径（如 .pdb、.gro）。
    trajectory_file : str
        轨迹文件路径（如 .xtc）。
    output_dir : str
        结果输出目录，将在此目录下生成 rmsd.csv。
    selection_protein : str
        蛋白质选择语句，默认 "protein"。
    selection_ligand : str
        配体选择语句，如 "resname LIG"；空串表示不计算配体/复合物。
    selection_backbone : str
        用于对齐的骨架选择，默认 "protein and backbone"。
    ref_frame : int
        参考帧索引（0 表示第一帧）。
    start_frame : int, optional
        从第几帧开始计算；None 表示从 ref_frame 开始。
    do_align : bool
        是否先对轨迹做 backbone 对齐。
    verbose : bool
        是否打印进度。

    Returns
    -------
    pd.DataFrame
        包含 Time (ns) 及各组 RMSD 的 DataFrame。
    """
    os.makedirs(output_dir, exist_ok=True)

    if verbose:
        print("加载轨迹中...")
    u = mda.Universe(topology_file, trajectory_file)

    if do_align:
        if verbose:
            print("执行轨迹对齐拟合（基于 backbone）...")
        align.AlignTraj(u, u, select=selection_backbone, in_memory=True).run()

    if selection_ligand and selection_ligand.strip():
        selection_complex = f"({selection_protein}) or ({selection_ligand})"
        group_sels = [selection_protein, selection_ligand, selection_complex]
        col_names = ["Protein", "Ligand", "Complex"]
    else:
        group_sels = [selection_protein]
        col_names = ["Protein"]

    if verbose:
        print("计算 RMSD...")
    rmsd_analyzer = rms.RMSD(
        u,
        select=selection_backbone,
        groupselections=group_sels,
        ref_frame=ref_frame,
    )
    run_start = start_frame if start_frame is not None else ref_frame
    rmsd_analyzer.run(start=run_start, verbose=verbose)

    # results.rmsd: 列 0=Frame, 1=Time(ps), 2=Backbone, 3+=groupselections
    rmsd_dict = {
        "Time (ns)": rmsd_analyzer.results.rmsd[:, 1] / 1000.0,
        "Backbone_RMSD (Å)": rmsd_analyzer.results.rmsd[:, 2],
    }
    for i, name in enumerate(col_names):
        rmsd_dict[f"{name}_RMSD (Å)"] = rmsd_analyzer.results.rmsd[:, 3 + i]

    df_rmsd = pd.DataFrame(rmsd_dict)
    out_csv = os.path.join(output_dir, "rmsd.csv")
    df_rmsd.to_csv(out_csv, index=False)
    if verbose:
        print(f"RMSD 结果已保存: {out_csv}")
    return df_rmsd


def plot_rmsd(df_rmsd, output_dir, save_path=None, font_family="Arial", font_size=12, dpi=300):
    """
    绘制 RMSD 曲线图并保存。

    Parameters
    ----------
    df_rmsd : pd.DataFrame
        由 run_rmsd_analysis 返回的 DataFrame。
    output_dir : str
        输出目录；若 save_path 未指定，则在此目录下保存 rmsd.png。
    save_path : str, optional
        指定时直接保存到该路径。
    font_family : str
        字体。
    font_size : int
        基础字体大小。
    dpi : int
        图像分辨率。
    """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    plt.rcParams["font.family"] = font_family
    plt.rcParams["font.size"] = font_size
    plt.rcParams["figure.dpi"] = dpi

    time_col = "Time (ns)"
    if time_col not in df_rmsd.columns:
        time_col = df_rmsd.columns[0]

    rmsd_cols = [c for c in df_rmsd.columns if "RMSD" in c and c != time_col]
    if not rmsd_cols:
        print("未找到 RMSD 列，跳过绘图")
        return

    fig, ax = plt.subplots(figsize=(8, 5))
    for col in rmsd_cols:
        ax.plot(df_rmsd[time_col], df_rmsd[col], label=col.replace("_RMSD (Å)", "").replace("_", " "))

    ax.set_xlabel("Time (ns)", fontsize=font_size)
    ax.set_ylabel("RMSD (Å)", fontsize=font_size)
    ax.set_title("RMSD Analysis", fontsize=14)
    ax.legend()
    ax.grid(True, alpha=0.3)
    plt.tight_layout()

    path = save_path if save_path else os.path.join(output_dir, "rmsd.png")
    plt.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close()
    print(f"RMSD 图已保存: {path}")


def parse_args():
    parser = argparse.ArgumentParser(
        description="对分子动力学轨迹进行 RMSD 分析（MDAnalysis），可选绘图。",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例:
  LNB-MDT RMSD -input_gro ./complex.pdb -input_xtc ./fit.xtc -output_dir ./results
  LNB-MDT RMSD -input_gro ./complex.pdb -input_xtc ./fit.xtc -protein_name "protein" -lig_name "resname LIG" -plot
        """,
    )
    parser.add_argument(
        "-input_gro", "-g", "-i",
        type=str,
        required=True,
        dest="input_gro",
        help="拓扑文件路径（.pdb 或 .gro）。",
    )
    parser.add_argument(
        "-input_xtc", "-x", "-t",
        type=str,
        required=True,
        dest="input_xtc",
        help="轨迹文件路径（如 .xtc）。",
    )
    parser.add_argument(
        "-output_dir", "-o",
        type=str,
        default="./results",
        dest="output_dir",
        help="结果输出目录（默认: ./results）。",
    )
    parser.add_argument(
        "-protein_name", "-protein", "-p",
        type=str,
        default="protein",
        dest="protein_name",
        help="蛋白质选择语句（默认: protein）。",
    )
    parser.add_argument(
        "-lig_name", "-lig", "-l",
        type=str,
        default="",
        dest="lig_name",
        help="配体选择语句，如 'resname LIG'；不填则只计算蛋白。",
    )
    parser.add_argument(
        "-backbone",
        type=str,
        default="protein and backbone",
        dest="backbone",
        help="用于对齐的骨架选择（默认: protein and backbone）。",
    )
    parser.add_argument(
        "-ref_frame", "-ref",
        type=int,
        default=0,
        dest="ref_frame",
        help="参考帧索引（默认: 0）。",
    )
    parser.add_argument(
        "-start_frame", "-start", "-s",
        type=int,
        default=None,
        dest="start_frame",
        help="从第几帧开始计算；默认与 ref_frame 一致。",
    )
    parser.add_argument(
        "-no_align",
        action="store_true",
        dest="no_align",
        help="不做轨迹对齐，直接计算 RMSD。",
    )
    parser.add_argument(
        "-plot", "-pfig",
        action="store_true",
        dest="plot",
        help="计算完成后绘制 RMSD 曲线并保存为 PNG。",
    )
    parser.add_argument(
        "-quiet", "-q",
        action="store_true",
        dest="quiet",
        help="减少终端输出。",
    )
    return parser.parse_args()


def main():
    args = parse_args()

    topo = os.path.abspath(os.path.expanduser(args.input_gro))
    traj = os.path.abspath(os.path.expanduser(args.input_xtc))
    out_dir = os.path.abspath(os.path.expanduser(args.output_dir))

    if not os.path.exists(topo):
        print(f"错误: 拓扑文件不存在: {topo}")
        sys.exit(1)
    if not os.path.exists(traj):
        print(f"错误: 轨迹文件不存在: {traj}")
        sys.exit(1)

    df = run_rmsd_analysis(
        topology_file=topo,
        trajectory_file=traj,
        output_dir=out_dir,
        selection_protein=args.protein_name,
        selection_ligand=args.lig_name.strip(),
        selection_backbone=args.backbone,
        ref_frame=args.ref_frame,
        start_frame=args.start_frame,
        do_align=not args.no_align,
        verbose=not args.quiet,
    )

    if args.plot:
        plot_rmsd(df, out_dir)


if __name__ == "__main__":
    main()
