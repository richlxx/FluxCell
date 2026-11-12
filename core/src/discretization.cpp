/// @file discretization.cpp
/// @brief 1D FVM 空间离散化实现
#include <fluxcell/core/discretization/discretization.hpp>

#include <cmath>
#include <numeric>
#include <stdexcept>

namespace fluxcell::discretization {

// ══════════════════════════════════════════════════════════════
//  Mesh1D
// ══════════════════════════════════════════════════════════════

Mesh1D Mesh1D::uniform(real_t x_start, real_t x_end, int n_cells) {
    if (n_cells <= 0 || x_end <= x_start) {
        throw std::invalid_argument("Mesh1D::uniform: invalid parameters");
    }

    Mesh1D mesh;
    mesh.n_cells_ = n_cells;
    mesh.x_start_ = x_start;
    mesh.x_end_ = x_end;

    real_t dx = (x_end - x_start) / n_cells;

    // 面位置
    mesh.x_f_.resize(n_cells + 1);
    for (int i = 0; i <= n_cells; ++i) {
        mesh.x_f_[i] = x_start + i * dx;
    }

    // 控制体宽度和中心
    mesh.dx_.resize(n_cells, dx);
    mesh.x_c_.resize(n_cells);
    for (int i = 0; i < n_cells; ++i) {
        mesh.x_c_[i] = 0.5 * (mesh.x_f_[i] + mesh.x_f_[i + 1]);
    }

    return mesh;
}

Mesh1D Mesh1D::nonuniform(real_t x_start, real_t x_end, int n_cells, real_t ratio) {
    if (n_cells <= 0 || x_end <= x_start || ratio <= 0.0) {
        throw std::invalid_argument("Mesh1D::nonuniform: invalid parameters");
    }

    Mesh1D mesh;
    mesh.n_cells_ = n_cells;
    mesh.x_start_ = x_start;
    mesh.x_end_ = x_end;

    // 等比数列：dx_i = dx_0 * ratio^i，总和 = L
    real_t L = x_end - x_start;
    real_t dx0;
    if (std::abs(ratio - 1.0) < 1e-12) {
        dx0 = L / n_cells;
    } else {
        dx0 = L * (ratio - 1.0) / (std::pow(ratio, n_cells) - 1.0);
    }

    mesh.x_f_.resize(n_cells + 1);
    mesh.dx_.resize(n_cells);
    mesh.x_c_.resize(n_cells);

    mesh.x_f_[0] = x_start;
    real_t dx_i = dx0;
    for (int i = 0; i < n_cells; ++i) {
        mesh.dx_[i] = dx_i;
        mesh.x_f_[i + 1] = mesh.x_f_[i] + dx_i;
        mesh.x_c_[i] = 0.5 * (mesh.x_f_[i] + mesh.x_f_[i + 1]);
        dx_i *= ratio;
    }

    return mesh;
}

// ══════════════════════════════════════════════════════════════
//  SphericalMesh
// ══════════════════════════════════════════════════════════════

SphericalMesh SphericalMesh::create(real_t radius, int n_cells) {
    SphericalMesh smesh;
    smesh.radius_ = radius;
    smesh.mesh_ = Mesh1D::uniform(0.0, radius, n_cells);

    // 计算各控制体体积和面面积
    const int n = n_cells;
    smesh.volumes_.resize(n);
    smesh.face_areas_.resize(n + 1);

    static constexpr real_t PI = 3.14159265358979323846;

    for (int i = 0; i < n; ++i) {
        real_t r_in = smesh.mesh_.face(i);
        real_t r_out = smesh.mesh_.face(i + 1);
        smesh.volumes_[i] = (4.0 / 3.0) * PI * (r_out * r_out * r_out - r_in * r_in * r_in);
    }

    for (int i = 0; i <= n; ++i) {
        real_t r = smesh.mesh_.face(i);
        smesh.face_areas_[i] = 4.0 * PI * r * r;
    }

    return smesh;
}

// ══════════════════════════════════════════════════════════════
//  FVM 离散化算子
// ══════════════════════════════════════════════════════════════

SparseMatrix assemble_diffusion_operator(const Mesh1D& mesh, const std::vector<real_t>& D) {
    const int n = mesh.num_cells();
    TripletList triplets;
    triplets.reserve(3 * n);

    for (int i = 0; i < n; ++i) {
        real_t diag = 0.0;

        // 左面通量 (i-1/2)
        if (i > 0) {
            real_t D_face = 0.5 * (D[i - 1] + D[i]);  // 算术平均
            real_t dist = 0.5 * (mesh.cell_width(i - 1) + mesh.cell_width(i));
            real_t coeff = D_face / (dist * mesh.cell_width(i));
            triplets.emplace_back(i, i - 1, coeff);
            diag -= coeff;
        }

        // 右面通量 (i+1/2)
        if (i < n - 1) {
            real_t D_face = 0.5 * (D[i] + D[i + 1]);
            real_t dist = 0.5 * (mesh.cell_width(i) + mesh.cell_width(i + 1));
            real_t coeff = D_face / (dist * mesh.cell_width(i));
            triplets.emplace_back(i, i + 1, coeff);
            diag -= coeff;
        }

        triplets.emplace_back(i, i, diag);
    }

    SparseMatrix op(n, n);
    op.setFromTriplets(triplets.begin(), triplets.end());
    return op;
}

SparseMatrix assemble_diffusion_operator(const Mesh1D& mesh, real_t D) {
    std::vector<real_t> D_vec(mesh.num_cells(), D);
    return assemble_diffusion_operator(mesh, D_vec);
}

SparseMatrix assemble_spherical_diffusion(const SphericalMesh& mesh, real_t D) {
    // (1/r²) ∂/∂r (r² D ∂c/∂r)  →  (1/V_i) [A_{i+1/2} D (c_{i+1}-c_i)/Δr - A_{i-1/2} D (c_i-c_{i-1})/Δr]
    const int n = mesh.num_cells();
    const auto& m = mesh.mesh();
    TripletList triplets;
    triplets.reserve(3 * n);

    for (int i = 0; i < n; ++i) {
        real_t diag = 0.0;
        real_t vol = mesh.volumes()[i];

        // 内面 (i-1/2)
        if (i > 0) {
            real_t A = mesh.face_areas()[i];
            real_t dr = 0.5 * (m.cell_width(i - 1) + m.cell_width(i));
            real_t coeff = A * D / (dr * vol);
            triplets.emplace_back(i, i - 1, coeff);
            diag -= coeff;
        }

        // 外面 (i+1/2)
        if (i < n - 1) {
            real_t A = mesh.face_areas()[i + 1];
            real_t dr = 0.5 * (m.cell_width(i) + m.cell_width(i + 1));
            real_t coeff = A * D / (dr * vol);
            triplets.emplace_back(i, i + 1, coeff);
            diag -= coeff;
        }

        triplets.emplace_back(i, i, diag);
    }

    SparseMatrix op(n, n);
    op.setFromTriplets(triplets.begin(), triplets.end());
    return op;
}

void apply_neumann_bc(SparseMatrix& matrix, Vector& rhs,
                      const Mesh1D& mesh, real_t D,
                      real_t flux_left, real_t flux_right) {
    const int n = mesh.num_cells();

    // 左边界：修正第一个控制体
    // 通量 = D * (c[0] - c_ghost) / (dx/2) = flux_left
    // 等效：在对角线加 D/(dx/2 * dx)，RHS 加 flux_left/dx
    real_t dx0 = mesh.cell_width(0);
    rhs(0) += D * flux_left / (0.5 * dx0 * dx0);

    // 右边界：修正最后一个控制体
    real_t dxN = mesh.cell_width(n - 1);
    rhs(n - 1) -= D * flux_right / (0.5 * dxN * dxN);
}

} // namespace fluxcell::discretization
