import warnings
warnings.filterwarnings('ignore')

import numpy as np
import scipy.sparse
from MDAnalysis.lib.distances import capped_distance
import matplotlib.pyplot as plt

try:
    from .analysis_base import *
except:
    from analysis_base import *

class ClusterHeight(AnalysisBase):
    def __init__(self, universe, cluster_residues_group, normals_residues_group, N_cutoff=10, cutoff=12, k=14, file_path=None, plot_path=None):
        """
        初始化 ClusterHeight 类，用于分析簇内和簇外脂质的高度。

        参数：
            universe: MDAnalysis.Universe 对象
            cluster_residues_group: dict，指定用于簇分析的脂质及其头部和尾部原子，例如 {'DPPC': (['PO4'], ['C4A', 'C4B']), 'DAPC': (['PO4'], ['C5A', 'C5B'])}
            normals_residues_group: dict，指定用于法向量计算的脂质及其头部原子，例如 {'DPPC': ['PO4'], 'DAPC': ['PO4'], 'DOPC': ['PO4']}
            N_cutoff: int，大于此值的簇被认为是“大簇”，默认 10
            cutoff: float，簇分析的距离截止值（Å），默认 12
            k: int，计算法向量时最近邻数，默认 14
            file_path: str，可选，输出结果的文件路径
            plot_path: str，可选，保存图表的路径
        """
        super().__init__(universe.trajectory)
        self.u = universe
        self.cluster_residues_group = cluster_residues_group
        self.normals_residues_group = normals_residues_group
        self.N_cutoff = N_cutoff
        self.cutoff = cutoff
        self.k = k
        self.file_path = file_path
        self.plot_path = plot_path

        # 选择簇分析的头部原子
        sel_str_cluster = ' or '.join([f'(resname {res} and name {" ".join(atoms[0])})' for res, atoms in cluster_residues_group.items()])
        self.clusterHeadAtoms = self.u.select_atoms(sel_str_cluster, updating=False)
        if self.clusterHeadAtoms.n_atoms != self.clusterHeadAtoms.n_residues:
            raise ValueError('每个簇分析残基必须恰好指定一个头部原子')
        self._n_cluster_residues = self.clusterHeadAtoms.n_residues
        self.cluster_resids = self.clusterHeadAtoms.resids
        self.cluster_resnames = self.clusterHeadAtoms.resnames

        # 选择法向量计算的头部原子
        sel_str_normals = ' or '.join([f'(resname {res} and name {atoms[0]})' for res, atoms in normals_residues_group.items()])
        self.normalsHeadAtoms = self.u.select_atoms(sel_str_normals, updating=False)
        if self.normalsHeadAtoms.n_atoms != self.normalsHeadAtoms.n_residues:
            raise ValueError('每个法向量计算残基必须恰好指定一个头部原子')

        # 选择尾部原子（基于簇分析脂质）
        sel_str_tail = ' or '.join([f'(resname {res} and name {" ".join(atoms[1])})' for res, atoms in cluster_residues_group.items()])
        self.tailAtoms = self.u.select_atoms(sel_str_tail, updating=False)
        if self.tailAtoms.n_residues != self._n_cluster_residues:
            raise ValueError('尾部原子的残基数必须与头部匹配')

    def _prepare(self):
        """准备结果数组，存储每帧簇内、簇外和总体的平均 Height"""
        self.results.height_in_cluster = np.full(self.n_frames, np.nan)
        self.results.height_out_cluster = np.full(self.n_frames, np.nan)
        self.results.height_average = np.full(self.n_frames, np.nan)

    def _single_frame(self):
        """单帧分析：计算簇并分离 Height"""
        # 簇分析
        positions_cluster = self.clusterHeadAtoms.positions
        pairs = capped_distance(positions_cluster, positions_cluster, max_cutoff=self.cutoff,
                               box=self._ts.dimensions, return_distances=False)
        ref, nei = pairs.T
        index_needs = ref != nei
        ref = ref[index_needs]
        nei = nei[index_needs]
        data = np.ones_like(ref)
        neighbours_frame = scipy.sparse.csr_matrix((data, (ref, nei)),
                                                  shape=(self._n_cluster_residues, self._n_cluster_residues))
        _, com_labels = scipy.sparse.csgraph.connected_components(neighbours_frame)
        unique_com_labels, counts = np.unique(com_labels, return_counts=True)
        large_clusters = unique_com_labels[counts > self.N_cutoff]
        in_large_cluster = np.isin(com_labels, large_clusters)

        # 计算法向量
        positions_normals = self.normalsHeadAtoms.positions
        normals = get_normals(self.k, positions_normals)  # 假设 get_normals 已定义

        # 映射法向量到簇分析脂质
        cluster_indices = self.clusterHeadAtoms.residues.ix
        normals_indices = self.normalsHeadAtoms.residues.ix
        mapping = np.searchsorted(normals_indices, cluster_indices)
        cluster_normals = normals[mapping]
        print(f"Frame {self._frame_index}, cluster_normals shape: {cluster_normals.shape}")

        # 计算高度
        center_head = self.clusterHeadAtoms.positions
        center_tail = self.tailAtoms.center_of_geometry(compound='residues')
        head_to_tail = center_head - center_tail
        height = np.abs(np.einsum('ij,ij->i', cluster_normals, head_to_tail)) * 0.1  # Å to nm
        print(f"Frame {self._frame_index}, height shape: {height.shape}")

        # 分离簇内和簇外的 Height 并计算平均值
        self.results.height_in_cluster[self._frame_index] = (np.mean(height[in_large_cluster])
                                                            if np.any(in_large_cluster) else np.nan)
        self.results.height_out_cluster[self._frame_index] = (np.mean(height[~in_large_cluster])
                                                             if np.any(~in_large_cluster) else np.nan)
        self.results.height_average[self._frame_index] = np.mean(height)  # 总体平均高度

    def _conclude(self):
        """结束时将结果写入文件并绘制图表"""
        if self.file_path:
            with open(self.file_path, 'w') as f:
                f.write('frame,height_in_cluster,height_out_cluster,height_average\n')
                for i in range(self.n_frames):
                    f.write(f'{i},{self.results.height_in_cluster[i]},{self.results.height_out_cluster[i]},{self.results.height_average[i]}\n')

        if self.plot_path:
            plt.figure(figsize=(12, 6))
            frames = np.arange(self.n_frames)
            ls = [self.results.height_in_cluster, self.results.height_out_cluster, self.results.height_average]
            txt = ['Height_in_cluster', 'Height_out_cluster', 'Height_average']
            color = ['blue', 'orange', 'green']
            for i in range(3):
                plt.plot(frames, ls[i], label=txt[i], color=color[i], linewidth=2)
                plt.xlabel('Frame', fontsize=12)
                plt.ylabel('Average SZ', fontsize=12)
                plt.title('Average Order Parameter (SZ) Inside, Outside Clusters, and Overall', fontsize=14)
                plt.legend(fontsize=10)
                plt.grid(True, linestyle='--', alpha=0.7)
                plt.tight_layout()
                plt.show()

# 示例使用
if __name__ == "__main__":
    import MDAnalysis as mda
    u = mda.Universe(r"E:\awork\lnb\20250311\lnb.gro", r"E:\awork\lnb\20250311\md3.part0003.xtc")
    cluster_group = {'DPPC': (['PO4'], ['C4A', 'C4B']), 'CHOL': (['ROH'], ['R5'])}
    normals_group = {'DPPC': ['PO4'], 'DUPC': ['PO4'], 'CHOL': ['ROH']}
    analysis = ClusterHeight(u, cluster_group, normals_group,
                             N_cutoff=10, cutoff=12, k=21,
                             file_path='E:/clusterheight.csv', plot_path='E:/clusterheight_plot.png')
    analysis.run(step=3, verbose=True)