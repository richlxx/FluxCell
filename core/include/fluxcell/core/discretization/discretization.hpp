/// @file discretization.hpp
/// @brief 1D 有限体积法（FVM）空间离散化
#pragma once

#include <fluxcell/core/types.hpp>

#include <vector>

namespace fluxcell::discretization {

/// 1D 网格
class Mesh1D {
public:
    /// 构造均匀网格
    /// @param x_start  起始坐标 [m]
    /// @param x_end    终止坐标 [m]
    /// @param n_cells  控制体数量
    static Mesh1D uniform(real_t x_start, real_t x_end, int n_cells);

    /// 构造非均匀网格（支持边界加密）
    /// @param x_start  起始坐标 [m]
    /// @param x_end    终止坐标 [m]
    /// @param n_cells  控制体数量
    /// @param ratio    边界加密比（>1 时边界网格更密）
    static Mesh1D nonuniform(real_t x_start, real_t x_end, int n_cells, real_t ratio = 1.5);

    int num_cells() const { return n_cells_; }
    real_t x_start() const { return x_start_; }
    real_t x_end() const { return x_end_; }
    real_t length() const { return x_end_ - x_start_; }

    /// 控制体中心坐标
    const std::vector<real_t>& cell_centers() const { return x_c_; }
    real_t cell_center(int i) const { return x_c_[i]; }

    /// 控制体面坐标（n_cells + 1 个面）
    const std::vector<real_t>& faces() const { return x_f_; }
    real_t face(int i) const { return x_f_[i]; }

    /// 控制体宽度
    const std::vector<real_t>& cell_widths() const { return dx_; }
    real_t cell_width(int i) const { return dx_[i]; }

private:
    int n_cells_ = 0;
    real_t x_start_ = 0.0;
    real_t x_end_ = 0.0;
    std::vector<real_t> x_c_;   ///< 控制体中心
    std::vector<real_t> x_f_;   ///< 面位置
    std::vector<real_t> dx_;    ///< 控制体宽度
};

/// 球坐标 1D 网格（用于颗粒内部扩散）
class SphericalMesh {
public:
    /// @param radius  颗粒半径 [m]
    /// @param n_cells 径向控制体数
    static SphericalMesh create(real_t radius, int n_cells);

    int num_cells() const { return mesh_.num_cells(); }
    real_t radius() const { return radius_; }
    const Mesh1D& mesh() const { return mesh_; }

    /// 各控制体体积 V_i = (4/3)π (r_{i+1}^3 - r_i^3)
    const std::vector<real_t>& volumes() const { return volumes_; }

    /// 各面面积 A_i = 4π r_i^2
    const std::vector<real_t>& face_areas() const { return face_areas_; }

private:
    real_t radius_ = 0.0;
    Mesh1D mesh_;
    std::vector<real_t> volumes_;
    std::vector<real_t> face_areas_;
};

// ══════════════════════════════════════════════════════════════
//  FVM 离散化工具函数
// ══════════════════════════════════════════════════════════════

/// 组装 1D 扩散算子 ∇·(D ∇c) 的稀疏矩阵
/// 使用中心差分格式：D_{i+1/2} (c_{i+1} - c_i) / dx_{i+1/2}
/// @param mesh   1D 网格
/// @param D      各控制体的扩散系数 [m²/s]
/// @return 三对角稀疏矩阵
SparseMatrix assemble_diffusion_operator(const Mesh1D& mesh, const std::vector<real_t>& D);

/// 组装 1D 扩散算子（常数扩散系数）
SparseMatrix assemble_diffusion_operator(const Mesh1D& mesh, real_t D);

/// 组装球坐标扩散算子 (1/r²) ∂/∂r (r² D ∂c/∂r)
/// @param mesh   球坐标网格
/// @param D      扩散系数 [m²/s]
/// @return 稀疏矩阵
SparseMatrix assemble_spherical_diffusion(const SphericalMesh& mesh, real_t D);

/// 施加 Neumann 边界条件（通量边界）到扩散矩阵和 RHS
/// @param matrix  扩散矩阵（会被修改）
/// @param rhs     右端向量（会被修改）
/// @param mesh    网格
/// @param D       边界处扩散系数
/// @param flux_left  左边界通量（正=流入）
/// @param flux_right 右边界通量（正=流出）
void apply_neumann_bc(SparseMatrix& matrix, Vector& rhs,
                      const Mesh1D& mesh, real_t D,
                      real_t flux_left, real_t flux_right);

} // namespace fluxcell::discretization
