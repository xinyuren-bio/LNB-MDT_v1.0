import warnings
warnings.filterwarnings('ignore')

import numpy as np
import scipy.sparse
from MDAnalysis.lib.distances import capped_distance
import matplotlib.pyplot as plt
from scipy.spatial import KDTree
from numpy.linalg import lstsq

try:
    from .analysis_base import *
except:
    from analysis_base import *

class ClusterCurvature(AnalysisBase):
    def __init__(self, universe, cluster_residues_group, normals_residues_group, N_cutoff=10, cutoff=12, k=14,
                 file_path=None, plot_path=None, method='mean'):
        super().__init__(universe.trajectory)
        self.u = universe
        self.cluster_residues_group = cluster_residues_group
        self.normals_residues_group = normals_residues_group
        self.N_cutoff = N_cutoff
        self.cutoff = cutoff
        self.k = k
        self.file_path = file_path
        self.plot_path = plot_path
        self.method = method

        # 选择用于簇分析的头部原子
        sel_str_cluster = ' or '.join([f'(resname {res} and name {atoms[0]})' for res, atoms in cluster_residues_group.items()])
        self.clusterHeadAtoms = self.u.select_atoms(sel_str_cluster, updating=False)
        if self.clusterHeadAtoms.n_atoms != self.clusterHeadAtoms.n_residues:
            raise ValueError('每个簇分析残基必须恰好指定一个头部原子')
        self._n_residues = self.clusterHeadAtoms.n_residues
        self.resids = self.clusterHeadAtoms.resids
        self.resnames = self.clusterHeadAtoms.resnames

        # 选择用于法向量计算的头部原子
        sel_str_normals = ' or '.join([f'(resname {res} and name {atoms[0]})' for res, atoms in normals_residues_group.items()])
        self.normalsHeadAtoms = self.u.select_atoms(sel_str_normals, updating=False)
        if self.normalsHeadAtoms.n_atoms != self.normalsHeadAtoms.n_residues:
            raise ValueError('每个法向量计算残基必须恰好指定一个头部原子')

    def _prepare(self):
        """准备结果数组，存储每帧的簇内、簇外和总体平均曲率"""
        self.results = type('Results', (), {})()
        self.results.curvature_in_cluster = np.full(self.n_frames, np.nan)
        self.results.curvature_out_cluster = np.full(self.n_frames, np.nan)
        self.results.curvature_average = np.full(self.n_frames, np.nan)

    def _single_frame(self):
        """单帧分析：计算簇并分离曲率"""
        # 获取簇分析的头部原子位置
        cluster_positions = self.clusterHeadAtoms.positions

        # 簇分析
        pairs = capped_distance(cluster_positions, cluster_positions, max_cutoff=self.cutoff,
                               box=self._ts.dimensions, return_distances=False)
        ref, nei = pairs.T
        index_needs = ref != nei
        ref = ref[index_needs]
        nei = nei[index_needs]
        data = np.ones_like(ref)
        neighbours_frame = scipy.sparse.csr_matrix((data, (ref, nei)),
                                                  shape=(self._n_residues, self._n_residues))
        _, com_labels = scipy.sparse.csgraph.connected_components(neighbours_frame)
        unique_com_labels, counts = np.unique(com_labels, return_counts=True)
        large_clusters = unique_com_labels[counts > self.N_cutoff]
        in_large_cluster = np.isin(com_labels, large_clusters)

        # 获取法向量计算的头部原子位置
        normals_positions = self.normalsHeadAtoms.positions
        normals = self.get_normals(normals_positions)

        # 将法向量映射到簇分析的原子
        cluster_indices = self.clusterHeadAtoms.residues.ix
        normals_indices = self.normalsHeadAtoms.residues.ix
        mapping = np.searchsorted(normals_indices, cluster_indices)
        cluster_normals = normals[mapping]

        # 计算曲率
        curvature = self.calculate_curvature(cluster_positions, cluster_normals)

        # 分离簇内和簇外的曲率并计算平均值
        self.results.curvature_in_cluster[self._frame_index] = (np.mean(curvature[in_large_cluster])
                                                               if np.any(in_large_cluster) else np.nan)
        self.results.curvature_out_cluster[self._frame_index] = (np.mean(curvature[~in_large_cluster])
                                                                if np.any(~in_large_cluster) else np.nan)
        self.results.curvature_average[self._frame_index] = np.mean(curvature)

    def get_normals(self, positions):
        """计算每个原子的法向量"""
        kd_tree = KDTree(positions)
        normals = np.zeros((positions.shape[0], 3))
        for i, point in enumerate(positions):
            _, idxs = kd_tree.query(point, self.k)
            nearest_neighbors = positions[idxs]
            normal = self.fit_plane_and_get_normal(nearest_neighbors)
            normals[i] = normal
        return normals

    def fit_plane_and_get_normal(self, points):
        """拟合平面并获取法向量"""
        mean = np.mean(points, axis=0)
        centered_points = points - mean
        covariance_matrix = np.cov(centered_points.T)
        eigenvalues, eigenvectors = np.linalg.eig(covariance_matrix)
        normal = eigenvectors[:, np.argmin(eigenvalues)]
        return normal / np.linalg.norm(normal)

    def calculate_curvature(self, positions, normals):
        """计算每个原子的曲率"""
        kd_tree = KDTree(positions)
        curvature = np.zeros(positions.shape[0])
        for i, point in enumerate(positions):
            _, idxs = kd_tree.query(point, self.k)
            nearest_neighbors = positions[idxs]
            mean = np.mean(nearest_neighbors, axis=0)
            centered_points = nearest_neighbors - mean
            # Project z-coordinates onto the normal, keep x and y as is
            z_values = np.dot(centered_points, normals[i])
            local_points = np.column_stack((centered_points[:, 0], centered_points[:, 1], z_values))
            coefficients = self.fit_quadratic_surface(local_points)
            first_style = self.calculate_coefficients(coefficients)
            if self.method == 'mean':
                curvature[i] = self.calculate_mean_curvature(first_style)
            elif self.method == 'gaussian':
                curvature[i] = self.calculate_gaussian_curvature(first_style)
        return curvature

    def fit_quadratic_surface(self, points):
        """拟合二次曲面"""
        A = []
        b = []
        for x, y, z in points:
            A.append([x**2, y**2, x*y, x, y, 1])
            b.append(z)
        A = np.array(A)
        b = np.array(b)
        coefficients, _, _, _ = lstsq(A, b, rcond=None)
        return coefficients

    def calculate_coefficients(self, coefficients):
        """计算曲率系数"""
        E = (coefficients[3]**2) + 1
        F = coefficients[3] * coefficients[4]
        G = (coefficients[4]**2) + 1
        L = 2 * coefficients[0]
        M = coefficients[2]
        N = coefficients[1] * 2
        return np.array([E, F, G, L, M, N])

    def calculate_mean_curvature(self, coefficients):
        """计算平均曲率"""
        E, F, G, L, M, N = coefficients
        EN = E * N
        FM = F * M
        GL = G * L
        EG = E * G
        F2 = F**2
        k_mean = (EN - 2*FM + GL) / (2 * (EG - F2))
        return k_mean

    def calculate_gaussian_curvature(self, coefficients):
        """计算高斯曲率"""
        E, F, G, L, M, N = coefficients
        LN = L * N
        M2 = M**2
        EG = E * G
        F2 = F**2
        k_gaussian = (LN - M2) / (EG - F2)
        return k_gaussian

    def _conclude(self):
        """结束时保存结果并绘图"""
        if self.file_path:
            with open(self.file_path, 'w') as f:
                f.write('frame,curvature_in_cluster,curvature_out_cluster,curvature_average\n')
                for i in range(self.n_frames):
                    f.write(f'{i},{self.results.curvature_in_cluster[i]},{self.results.curvature_out_cluster[i]},{self.results.curvature_average[i]}\n')

        if self.plot_path:
            plt.figure(figsize=(12, 6))
            frames = np.arange(self.n_frames)
            ls = [self.results.curvature_in_cluster, self.results.curvature_out_cluster, self.results.curvature_average]
            txt = ['curvature_in_cluster', 'curvature_out_cluster', 'curvature_average']
            color = ['blue', 'orange', 'green']
            for i in range(3):
                plt.plot(frames, ls[i], label=txt[i], color=color[i], linewidth=2)
                plt.xlabel('Frame', fontsize=12)
                plt.ylabel('Curvature', fontsize=12)
                plt.title('Curvature Profiles', fontsize=14)
                plt.legend(fontsize=10)
                plt.grid(True, linestyle='--', alpha=0.7)
                plt.tight_layout()
                plt.show()
            for i in range(3):
                plt.plot(frames, ls[i], label=txt[i], color=color[i], linewidth=2)
            plt.xlabel('Frame', fontsize=12)
            plt.ylabel('Curvature', fontsize=12)
            plt.title('Curvature Profiles', fontsize=14)
            plt.legend(fontsize=10)
            plt.grid(True, linestyle='--', alpha=0.7)
            plt.tight_layout()
            plt.show()


# 示例使用
if __name__ == "__main__":
    import MDAnalysis as mda
    u = mda.Universe(r"E:\awork\lnb\20250311\lnb.gro", r"E:\awork\lnb\20250311\md3.part0003.xtc")
    cluster_group = {'DUPC': ['PO4']}  # 用于簇分析的脂质
    normals_group = {'DPPC': ['PO4'], 'DUPC': ['PO4'], 'CHOL': ['ROH']}  # 用于法向量计算的脂质
    analysis = ClusterCurvature(u, cluster_group, normals_group, N_cutoff=10, cutoff=12, k=20,
                                file_path='E:/clustercurvature.csv',
                                plot_path='E:/clustercurvature_plot.png',
                                method='mean')
    analysis.run(step=3, verbose=True)