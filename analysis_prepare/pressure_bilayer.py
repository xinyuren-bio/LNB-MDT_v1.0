from tqdm import tqdm
import numpy as np
from scipy.spatial import KDTree
from joblib import Parallel, delayed
from numba import njit


# 使用 Numba 加速周期性边界条件下的距离计算
@njit
def pbc_distance(pos1, pos2, box_size):
    """
    计算考虑周期性边界条件 (PBC) 的距离向量。
    参数:
    - pos1: ndarray, 粒子1的位置 (单位: nm)
    - pos2: ndarray, 粒子2的位置 (单位: nm)
    - box_size: tuple, 模拟盒子尺寸 (单位: nm)
    返回:
    - delta: ndarray, 最小距离向量
    """
    delta = pos2 - pos1
    for i in range(3):
        delta[i] -= round(delta[i] / box_size[i]) * box_size[i]
    return delta


# 并行计算单个粒子对的 Lennard-Jones 力
def compute_force_pair(i, j, positions, box_size, epsilon, sigma):
    r_vec = pbc_distance(positions[i], positions[j], box_size)
    r = np.linalg.norm(r_vec)
    if r == 0:
        return i, j, np.zeros(3), np.zeros(3)
    f_mag = -24 * epsilon * (2 * (sigma / r) ** 12 - (sigma / r) ** 6) / r
    f_ij = f_mag * r_vec
    return i, j, f_ij, -f_ij


def compute_lj_force(positions, box_size, epsilon=1.0, sigma=0.3, n_jobs=6):
    """
    并行计算 Lennard-Jones 力。
    参数:
    - positions: ndarray, shape (N, 3), 所有粒子位置 (单位: nm)
    - box_size: tuple, 模拟盒子尺寸 (单位: nm)
    - epsilon: float, LJ 势能深度 (单位: kcal/mol)
    - sigma: float, LJ 特征长度 (单位: nm)
    - n_jobs: int, 并行任务数，默认 4
    返回:
    - forces: ndarray, shape (N, N, 3), 粒子间力 (单位: kcal/mol/nm)
    """
    N = len(positions)
    tree = KDTree(positions)
    pairs = tree.query_pairs(1.2, output_type='ndarray')

    # 使用 joblib 并行计算所有粒子对的力
    results = Parallel(n_jobs=n_jobs)(
        delayed(compute_force_pair)(i, j, positions, box_size, epsilon, sigma)
        for i, j in pairs
    )

    # 初始化力数组并填充结果
    forces = np.zeros((N, N, 3))
    for i, j, f_ij, f_ji in results:
        forces[i, j] = f_ij
        forces[j, i] = f_ji
    return forces


# 并行计算单个粒子的动能贡献
def compute_kinetic_contribution(i, positions, velocities, masses, slab_thickness, num_slabs):
    pos = positions[i]
    vel = velocities[i]
    mass = masses[i]
    slab_idx = int(pos[2] / slab_thickness)
    if 0 <= slab_idx < num_slabs:
        contrib = np.outer(vel, vel) * mass
        return slab_idx, contrib
    return None


# 并行计算单个粒子的构型贡献
def compute_configurational_contribution(i, positions, forces, box_size, slab_thickness, num_slabs, N_segments):
    N = len(positions)
    contrib = np.zeros((num_slabs, 3, 3))
    for j in range(N):
        if i == j:
            continue
        r_ij = pbc_distance(positions[i], positions[j], box_size)
        f_ij = forces[i, j]
        for lambda_ in np.linspace(0, 1, N_segments + 1):
            l = positions[i] + lambda_ * r_ij
            slab_idx = int(l[2] / slab_thickness)
            if 0 <= slab_idx < num_slabs:
                contrib[slab_idx] += -np.outer(f_ij, r_ij) / N_segments
    return contrib


def calculate_planar_pressure_profile(positions, velocities, masses, box_size, slab_thickness=0.4, N_segments=50,
                                      epsilon=1.0, sigma=0.3, n_jobs=4):
    num_particles = len(positions)
    Lx, Ly, Lz = box_size
    num_slabs = int(Lz / slab_thickness) # 分成了n份

    slab_area = Lx * Ly

    # 计算 Lennard-Jones 力
    forces = compute_lj_force(positions, box_size, epsilon, sigma, n_jobs=n_jobs)
    pressure_field = np.zeros((num_slabs, 3, 3), dtype=np.float64)

    # 并行计算动能部分
    kinetic_results = Parallel(n_jobs=n_jobs)(
        delayed(compute_kinetic_contribution)(i, positions, velocities, masses, slab_thickness, num_slabs)
        for i in tqdm(range(num_particles), desc="计算动能部分")
    )
    for res in kinetic_results:
        if res is not None:
            slab_idx, contrib = res
            pressure_field[slab_idx] += contrib
    print('动能计算完成')

    # 并行计算构型部分
    config_results = Parallel(n_jobs=n_jobs)(
        delayed(compute_configurational_contribution)(
            i, positions, forces, box_size, slab_thickness, num_slabs, N_segments
        ) for i in tqdm(range(num_particles), desc="计算构型部分")
    )
    pressure_field += np.sum(config_results, axis=0)
    print('构型计算完成')

    # 归一化压力张量并转换为 bar
    pressure_field /= (slab_area * slab_thickness)
    pressure_field *= 6.9478e-4  # 转换为 bar

    # 计算 p_L 和 p_zz
    p_L = (pressure_field[:, 0, 0] + pressure_field[:, 1, 1]) / 2
    p_zz = pressure_field[:, 2, 2]
    z_bins = (np.arange(num_slabs) + 0.5) * slab_thickness

    valid = ~np.isnan(p_L) & ~np.isnan(p_zz)
    return z_bins[valid], p_L[valid], p_zz[valid]


# 示例用法保持不变
if __name__ == "__main__":
    import matplotlib.pyplot as plt
    import MDAnalysis as mda

    u = mda.Universe('E:/awork/lnb/pressure/1.pdb')
    N = u.atoms.n_atoms
    box_size = (5.0, 5.0, 8.5)
    positions = u.atoms.positions * 0.1
    T = 300.0  # 目标温度 (K)
    kB = 0.0083144621  # 玻尔兹曼常数 (kJ·mol⁻¹·K⁻¹)
    mass = 72.0  # 所有粒子的质量 (u)
    n_particles = N  # 粒子总数 (根据你的体系调整)
    n_dim = 3  # 每个粒子的速度分量 (x, y, z)

    # 随机种子 (可选，用于可重复性)
    np.random.seed(42)

    # 计算标准差 (统一质量)
    sigma = np.sqrt(kB * T / mass)  # 单位：nm/ps

    # 为所有粒子生成初始速度 (Maxwell-Boltzmann 分布)
    velocities = np.random.normal(0, sigma, (n_particles, n_dim))
    masses = np.full([N], fill_value=72)
    # 动量校正 (确保总动量为零)
    total_momentum = np.sum(mass * velocities, axis=0)
    center_of_mass_velocity = total_momentum / (n_particles * mass)
    velocities -= center_of_mass_velocity

    z_bins, p_L, p_zz = calculate_planar_pressure_profile(
        positions, velocities, masses, box_size, epsilon=1.0, sigma=0.3, n_jobs=4
    )

    plt.plot(z_bins, p_L * 0.1, label='$p_L(z)$', color='blue', linestyle='-', linewidth=2)
    plt.plot(z_bins, p_zz * 0.1, label='$p_{zz}(z)$', color='red', linestyle='--', linewidth=2)
    plt.legend(loc='upper right')
    plt.xlabel('z (nm)')
    plt.ylabel('Pressure (bar)')
    plt.title('Lateral and Normal Pressure Profiles in Planar Bilayer')
    plt.grid(True)
    plt.show()

    print("z坐标 (nm):", z_bins)
    print("p_L (bar):", p_L)
    print("p_zz (bar):", p_zz)