import numpy as np


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
    delta -= np.round(delta / box_size) * box_size
    return delta


def compute_lj_force(positions, box_size, epsilon=1.0, sigma=0.3):
    """
    计算Lennard-Jones力。

    参数:
    - positions: ndarray, shape (N, 3), 所有粒子位置 (单位: nm)
    - box_size: tuple, 模拟盒子尺寸 (单位: nm)
    - epsilon: float, LJ势能深度 (单位: kcal/mol)
    - sigma: float, LJ特征长度 (单位: nm)

    返回:
    - forces: ndarray, shape (N, N, 3), 粒子间力 (单位: kcal/mol/nm)
    """
    N = len(positions)
    forces = np.zeros((N, N, 3))
    for i in range(N):
        for j in range(i + 1, N):
            r_vec = pbc_distance(positions[i], positions[j], box_size)
            r = np.linalg.norm(r_vec)
            if r == 0:
                continue
            f_mag = 24 * epsilon * (2 * (sigma / r) ** 12 - (sigma / r) ** 6) / r ** 2
            f_ij = f_mag * r_vec
            forces[i, j] = f_ij
            forces[j, i] = -f_ij  # 牛顿第三定律
    return forces


def calculate_3d_pressure_field(positions, velocities, masses, box_size, cube_size=0.1, vesicle_center=None,
                                N_segments=100, epsilon=1.0, sigma=0.3):
    """
    计算单层球状磷脂膜的3D压力分布。

    参数:
    - positions: ndarray, shape (N, 3), 粒子位置 (单位: nm)
    - velocities: ndarray, shape (N, 3), 粒子速度 (单位: nm/ps)
    - masses: ndarray, shape (N,), 粒子质量 (单位: g/mol)
    - box_size: tuple, (Lx, Ly, Lz), 模拟盒子尺寸 (单位: nm)
    - cube_size: float, 立方体边长 (单位: nm), 默认0.1 nm
    - vesicle_center: tuple, (x0, y0, z0), 囊泡中心坐标 (单位: nm), 默认盒子中心
    - N_segments: int, IK轮廓分段数, 默认100
    - epsilon: float, LJ势能参数 (单位: kcal/mol)
    - sigma: float, LJ势能参数 (单位: nm)

    返回:
    - r_bins: ndarray, 径向距离分箱
    - p_rr: ndarray, 径向压力分量 p_rr(r)
    - p_T: ndarray, 切向压力分量 p_T(r)
    """

    # 参数初始化
    num_particles = len(positions)
    if vesicle_center is None:
        vesicle_center = np.array(box_size) / 2  # 默认囊泡中心为盒子中心

    # 计算Lennard-Jones力
    forces = compute_lj_force(positions, box_size, epsilon, sigma)

    # 计算立方体数量
    num_cubes = [int(box_size[i] / cube_size) for i in range(3)]
    volume = cube_size ** 3  # 每个立方体的体积 (nm^3)

    # 初始化压力张量场，shape: (nx, ny, nz, 3, 3)
    pressure_field = np.zeros(num_cubes + [3, 3], dtype=np.float64)

    # --- 计算动能部分 ---
    for i in range(num_particles):
        pos = positions[i]
        vel = velocities[i]
        mass = masses[i]
        # 确定粒子所在立方体
        cube_idx = np.floor(pos / cube_size).astype(int)
        # 检查边界
        if np.any(cube_idx < 0) or np.any(cube_idx >= num_cubes):
            continue
        # 添加动能贡献: m * v_alpha * v_beta
        for alpha in range(3):
            for beta in range(3):
                pressure_field[cube_idx[0], cube_idx[1], cube_idx[2], alpha, beta] += (
                        mass * vel[alpha] * vel[beta]
                )

    # --- 计算构型部分 (使用IK轮廓) ---
    for i in range(num_particles):
        for j in range(i + 1, num_particles):
            r_ij = pbc_distance(positions[i], positions[j], box_size)
            f_ij = forces[i, j]
            # 沿IK轮廓（直线）积分
            for lambda_ in np.linspace(0, 1, N_segments + 1):
                l = positions[i] + lambda_ * r_ij
                cube_idx = np.floor(l / cube_size).astype(int)
                if np.any(cube_idx < 0) or np.any(cube_idx >= num_cubes):
                    continue
                # 添加构型贡献: -f_alpha * r_beta / N_segments
                contrib = -f_ij[:, np.newaxis] * r_ij[np.newaxis, :] / N_segments
                pressure_field[cube_idx[0], cube_idx[1], cube_idx[2]] += contrib

    # 归一化压力张量 (单位转换为bar: 1 kcal/mol/nm^3 ≈ 694.76 bar)
    pressure_field /= volume  # 单位: g/mol * nm/ps^2 / nm^3
    pressure_field *= 4184 / 6.022e23 / 1e-24  # 转换为 bar

    # --- 转换为球坐标系 ---
    # 计算每个立方体的中心位置
    cube_centers = np.array([
        (np.arange(num_cubes[i]) + 0.5) * cube_size for i in range(3)
    ])
    X, Y, Z = np.meshgrid(cube_centers[0], cube_centers[1], cube_centers[2], indexing='ij')
    coords = np.stack([X, Y, Z], axis=-1)  # shape: (nx, ny, nz, 3)

    # 计算径向距离 r
    r = np.sqrt(np.sum((coords - vesicle_center) ** 2, axis=-1))  # shape: (nx, ny, nz)

    # 计算球坐标基矢量
    r_hat = (coords - vesicle_center) / r[..., np.newaxis]  # 径向单位向量
    # 计算 p_rr
    p_rr_field = np.einsum('...i,...ij,...j -> ...', r_hat, pressure_field, r_hat)

    # 切向分量 p_T 近似
    trace = np.trace(pressure_field, axis1=-2, axis2=-1)
    p_T_field = (trace - p_rr_field) / 2

    # --- 径向分箱平均 ---
    r_max = np.max(r)
    r_bins = np.linspace(0, r_max, 50)
    dr = r_bins[1] - r_bins[0]
    p_rr = np.zeros(len(r_bins) - 1)
    p_T = np.zeros(len(r_bins) - 1)
    counts = np.zeros(len(r_bins) - 1)

    r_flat = r.flatten()
    p_rr_flat = p_rr_field.flatten()
    p_T_flat = p_T_field.flatten()

    for i in range(len(r_bins) - 1):
        mask = (r_flat >= r_bins[i]) & (r_flat < r_bins[i + 1])
        if np.sum(mask) > 0:
            p_rr[i] = np.mean(p_rr_flat[mask])
            p_T[i] = np.mean(p_T_flat[mask])
            counts[i] = np.sum(mask)

    # 仅保留有效数据
    valid = counts > 0
    r_bins = (r_bins[:-1] + r_bins[1:]) / 2
    return r_bins[valid], p_rr[valid], p_T[valid]


# 示例用法
if __name__ == "__main__":
    # 模拟数据
    np.random.seed(42)
    N = 1000
    box_size = (10.0, 10.0, 10.0)  # nm
    positions = np.random.uniform(0, box_size[0], (N, 3))
    velocities = np.random.normal(0, 1, (N, 3))  # nm/ps
    masses = np.ones(N) * 18  # g/mol

    # 计算压力分布
    r_bins, p_rr, p_T = calculate_3d_pressure_field(
        positions, velocities, masses, box_size, epsilon=1.0, sigma=0.3
    )

    # 输出结果
    print("径向距离 (nm):", r_bins)
    print("p_rr (bar):", p_rr)
    print("p_T (bar):", p_T)