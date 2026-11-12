/// @file test_discretization.cpp
/// @brief FVM 离散化单元测试
#include <fluxcell/core/discretization/discretization.hpp>

#include <gtest/gtest.h>

#include <numeric>

using namespace fluxcell;
using namespace fluxcell::discretization;

// ══════════════════════════════════════════════════════════════
//  Mesh1D 测试
// ══════════════════════════════════════════════════════════════

TEST(Mesh1D, UniformCreation) {
    auto mesh = Mesh1D::uniform(0.0, 1.0, 10);
    EXPECT_EQ(mesh.num_cells(), 10);
    EXPECT_DOUBLE_EQ(mesh.x_start(), 0.0);
    EXPECT_DOUBLE_EQ(mesh.x_end(), 1.0);
    EXPECT_EQ(mesh.faces().size(), 11u);
    EXPECT_EQ(mesh.cell_centers().size(), 10u);
    EXPECT_EQ(mesh.cell_widths().size(), 10u);
}

TEST(Mesh1D, UniformSpacing) {
    auto mesh = Mesh1D::uniform(0.0, 100e-6, 20);
    real_t expected_dx = 100e-6 / 20.0;
    for (int i = 0; i < 20; ++i) {
        EXPECT_NEAR(mesh.cell_width(i), expected_dx, 1e-15);
    }
}

TEST(Mesh1D, UniformCellCenters) {
    auto mesh = Mesh1D::uniform(0.0, 1.0, 4);
    // 中心应在 0.125, 0.375, 0.625, 0.875
    EXPECT_NEAR(mesh.cell_center(0), 0.125, 1e-12);
    EXPECT_NEAR(mesh.cell_center(1), 0.375, 1e-12);
    EXPECT_NEAR(mesh.cell_center(2), 0.625, 1e-12);
    EXPECT_NEAR(mesh.cell_center(3), 0.875, 1e-12);
}

TEST(Mesh1D, UniformTotalLength) {
    auto mesh = Mesh1D::uniform(0.0, 1.0, 100);
    real_t total = std::accumulate(mesh.cell_widths().begin(), mesh.cell_widths().end(), 0.0);
    EXPECT_NEAR(total, 1.0, 1e-12);
}

TEST(Mesh1D, NonUniformTotalLength) {
    auto mesh = Mesh1D::nonuniform(0.0, 1.0, 20, 2.0);
    real_t total = std::accumulate(mesh.cell_widths().begin(), mesh.cell_widths().end(), 0.0);
    EXPECT_NEAR(total, 1.0, 1e-10);
}

TEST(Mesh1D, InvalidParameters) {
    EXPECT_THROW(Mesh1D::uniform(1.0, 0.0, 10), std::invalid_argument);
    EXPECT_THROW(Mesh1D::uniform(0.0, 1.0, 0), std::invalid_argument);
    EXPECT_THROW(Mesh1D::uniform(0.0, 1.0, -5), std::invalid_argument);
}

// ══════════════════════════════════════════════════════════════
//  SphericalMesh 测试
// ══════════════════════════════════════════════════════════════

TEST(SphericalMesh, TotalVolume) {
    real_t radius = 10e-6;
    auto smesh = SphericalMesh::create(radius, 20);

    // 总控制体体积应等于球体积 (4/3)πR³
    real_t total_vol = std::accumulate(smesh.volumes().begin(), smesh.volumes().end(), 0.0);
    real_t expected = (4.0 / 3.0) * 3.14159265358979323846 * radius * radius * radius;
    EXPECT_NEAR(total_vol, expected, expected * 1e-10);
}

TEST(SphericalMesh, FaceAreas) {
    real_t radius = 1.0;
    auto smesh = SphericalMesh::create(radius, 10);

    // 中心面面积应为 0
    EXPECT_NEAR(smesh.face_areas()[0], 0.0, 1e-15);
    // 外面面积应为 4πR²
    EXPECT_NEAR(smesh.face_areas().back(), 4.0 * 3.14159265358979323846, 1e-10);
}

// ══════════════════════════════════════════════════════════════
//  扩散算子测试
// ══════════════════════════════════════════════════════════════

TEST(DiffusionOperator, SymmetricConstantD) {
    auto mesh = Mesh1D::uniform(0.0, 1.0, 10);
    auto op = assemble_diffusion_operator(mesh, 1.0);

    // 扩散算子应为对称矩阵
    for (int k = 0; k < op.outerSize(); ++k) {
        for (SparseMatrix::InnerIterator it(op, k); it; ++it) {
            EXPECT_NEAR(op.coeff(it.row(), it.col()), op.coeff(it.col(), it.row()), 1e-12);
        }
    }
}

TEST(DiffusionOperator, ConservationProperty) {
    // 常数浓度场的扩散应为零
    auto mesh = Mesh1D::uniform(0.0, 1.0, 10);
    auto op = assemble_diffusion_operator(mesh, 1.0);

    Vector c = Vector::Constant(10, 5.0);
    Vector result = op * c;
    EXPECT_NEAR(result.norm(), 0.0, 1e-12);
}

TEST(SphericalDiffusion, ConservationProperty) {
    auto smesh = SphericalMesh::create(10e-6, 15);
    auto op = assemble_spherical_diffusion(smesh, 3.9e-14);

    Vector c = Vector::Constant(15, 1000.0);
    Vector result = op * c;
    // 内部节点应为零（边界除外）
    for (int i = 1; i < 14; ++i) {
        EXPECT_NEAR(result(i), 0.0, 1e-6);
    }
}
