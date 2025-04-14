import warnings

warnings.filterwarnings('ignore')

import numpy as np
import scipy.sparse
import scipy.stats
from MDAnalysis.lib.distances import capped_distance
import matplotlib.pyplot as plt

try:
    from .analysis_base import *
except:
    from analysis_base import *

class ClusterSZ(AnalysisBase):
    def __init__(self, universe, cluster_residues_group, normals_residues_group, chain, N_cutoff=10, cutoff=12, k=14,
                 file_path=None, plot_path=None):
        """
        初始化 ClusterSZ 类，用于分析簇内和簇外磷脂的序参数。

        参数：
            universe: MDAnalysis.Universe 对象
            cluster_residues_group: dict，指定用于簇分析的脂质及其头部原子，例如 {'DPPC': ['PO4'], 'DAPC': ['PO4']}
            normals_residues_group: dict，指定用于法向量计算的脂质及其头部原子，例如 {'DPPC': ['PO4'], 'DAPC': ['PO4'], 'DOPC': ['PO4']}
            chain: str，指定分析的尾链，'sn1'、'sn2' 或 'both'
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
        self.chain = chain
        self.N_cutoff = N_cutoff
        self.cutoff = cutoff
        self.k = k
        self.file_path = file_path
        self.plot_path = plot_path

        # 选择簇分析的头部原子
        sel_str_cluster = ' or '.join(
            [f'(resname {res} and name {atoms[0]})' for res, atoms in cluster_residues_group.items()])
        self.clusterHeadAtoms = self.u.select_atoms(sel_str_cluster, updating=False)
        if self.clusterHeadAtoms.n_atoms != self.clusterHeadAtoms.n_residues:
            raise ValueError('每个簇分析残基必须恰好指定一个头部原子')
        self._n_cluster_residues = self.clusterHeadAtoms.n_residues
        self.cluster_resids = self.clusterHeadAtoms.resids
        self.cluster_resnames = self.clusterHeadAtoms.resnames

        # 选择法向量计算的头部原子
        sel_str_normals = ' or '.join(
            [f'(resname {res} and name {atoms[0]})' for res, atoms in normals_residues_group.items()])
        self.normalsHeadAtoms = self.u.select_atoms(sel_str_normals, updating=False)
        if self.normalsHeadAtoms.n_atoms != self.normalsHeadAtoms.n_residues:
            raise ValueError('每个法向量计算残基必须恰好指定一个头部原子')

        # 设置尾部原子模式（基于簇分析脂质）
        tail_patterns = {'sn1': '??A', 'sn2': '??B', 'both': ['??A', '??B']}
        self.tail_chains = ['sn1', 'sn2'] if chain == 'both' else [chain]
        self.tailAtoms = {}
        for c in self.tail_chains:
            pattern = tail_patterns[c]
            sel_str = ' or '.join([f'(resname {res} and name {pattern})' for res in cluster_residues_group])
            self.tailAtoms[c] = self.u.select_atoms(sel_str, updating=False)
            if self.tailAtoms[c].n_atoms < self._n_cluster_residues * 3:
                raise ValueError(
                    f"尾部原子数量不足 for chain {c}: 找到 {self.tailAtoms[c].n_atoms} 个原子，预期至少 {self._n_cluster_residues * 3}")

    def _prepare(self):
        """准备结果数组，存储每帧簇内、簇外和总体的平均 SZ"""
        self.results.sz_in_cluster = np.full(self.n_frames, np.nan)
        self.results.sz_out_cluster = np.full(self.n_frames, np.nan)
        self.results.sz_average = np.full(self.n_frames, np.nan)

    def _single_frame(self):
        """单帧分析：计算簇并分离 SZ"""
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

        # 计算每条链的 SZ
        sz_per_chain = []
        for c in self.tail_chains:
            tail_atoms = self.tailAtoms[c]
            tail_positions = tail_atoms.positions.reshape([self._n_cluster_residues, -1, 3])

            if tail_positions.shape[1] < 3:
                raise ValueError(f"尾部原子数量不足，每残基至少需要3个原子，当前形状：{tail_positions.shape}")

            tail_tail = tail_positions[:, -1, :]
            head_to_tail = positions_cluster - tail_tail
            consistent_normals = ensure_consistent_normals(cluster_normals, head_to_tail)

            vectors = tail_positions[:, :-2, :] - tail_positions[:, 2:, :]
            vectors_norm = np.linalg.norm(vectors, axis=2, keepdims=True)

            dot_products = np.sum(vectors * consistent_normals[:, np.newaxis, :], axis=2) / vectors_norm[:, :, 0]

            theta = np.mean(dot_products, axis=1)
            sz_chain = (3 * theta ** 2 - 1) / 2
            sz_per_chain.append(sz_chain)

        sz = np.mean(sz_per_chain, axis=0) if len(self.tail_chains) > 1 else sz_per_chain[0]

        # 分离簇内和簇外的 SZ 并计算平均值
        self.results.sz_in_cluster[self._frame_index] = (np.mean(sz[in_large_cluster])
                                                         if np.any(in_large_cluster) else np.nan)
        self.results.sz_out_cluster[self._frame_index] = (np.mean(sz[~in_large_cluster])
                                                          if np.any(~in_large_cluster) else np.nan)
        self.results.sz_average[self._frame_index] = np.mean(sz)  # 计算所有脂质的平均 SZ

    def _conclude(self):
        """结束时将结果写入文件并绘制图表"""
        if self.file_path:
            with open(self.file_path, 'w') as f:
                f.write('frame,sz_in_cluster,sz_out_cluster,sz_average\n')
                for i in range(self.n_frames):
                    f.write(
                        f'{i},{self.results.sz_in_cluster[i]},{self.results.sz_out_cluster[i]},{self.results.sz_average[i]}\n')

        if self.plot_path:
            plt.figure(figsize=(12, 6))
            frames = np.arange(self.n_frames)
            ls = [self.results.sz_in_cluster, self.results.sz_out_cluster, self.results.sz_average]
            txt = ['SZ in Clusters', 'SZ out of Clusters', 'SZ Average']
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

            # 保存图片时使用
            # plt.savefig(self.plot_path, dpi=300, format='png')
            # plt.close()

# 示例使用
if __name__ == "__main__":
    import MDAnalysis as mda

    u = mda.Universe(r"E:\awork\lnb\20250311\lnb.gro", r"E:\awork\lnb\20250311\md3.part0003.xtc")
    cluster_group = {'DUPC': ['PO4']}
    normals_group = {'DPPC': ['PO4'], 'DUPC': ['PO4'], 'CHOL': ['ROH']}
    analysis = ClusterSZ(u, cluster_group, normals_group, 'both',
                         N_cutoff=10, cutoff=12, k=21,
                         file_path='E:/clustersz.csv', plot_path='E:/clustersz_plot.png')
    analysis.run(step=3, verbose=True)