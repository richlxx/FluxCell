/// @file test_physics.cpp
/// @brief 电化学本构关系单元测试
#include <fluxcell/core/physics/physics.hpp>

#include <gtest/gtest.h>

#include <cmath>

using namespace fluxcell;
using namespace fluxcell::physics;

// ══════════════════════════════════════════════════════════════
//  Butler-Volmer 测试
// ══════════════════════════════════════════════════════════════

TEST(ButlerVolmer, ZeroOverpotential) {
    // η=0 时 j=0
    real_t j = butler_volmer(0.0, 1000.0, 15000.0, 30000.0, 2e-11, 0.5, 0.5, 298.15);
    EXPECT_NEAR(j, 0.0, 1e-15);
}

TEST(ButlerVolmer, PositiveOverpotential) {
    // η>0 时 j>0（阳极方向）
    real_t j = butler_volmer(0.01, 1000.0, 15000.0, 30000.0, 2e-11, 0.5, 0.5, 298.15);
    EXPECT_GT(j, 0.0);
}

TEST(ButlerVolmer, NegativeOverpotential) {
    // η<0 时 j<0（阴极方向）
    real_t j = butler_volmer(-0.01, 1000.0, 15000.0, 30000.0, 2e-11, 0.5, 0.5, 298.15);
    EXPECT_LT(j, 0.0);
}

TEST(ButlerVolmer, SymmetricDerivative) {
    // dj/dη 在 η=0 处应对称
    real_t dj = butler_volmer_deta(0.0, 1000.0, 15000.0, 30000.0, 2e-11, 0.5, 0.5, 298.15);
    EXPECT_GT(dj, 0.0); // 导数应为正
}

// ══════════════════════════════════════════════════════════════
//  OCP 测试
// ══════════════════════════════════════════════════════════════

TEST(GraphiteOCP, RangeCheck) {
    // 石墨 OCP 应在 0~2V 范围内
    for (real_t x = 0.01; x <= 0.99; x += 0.05) {
        real_t u = graphite_ocp(x);
        EXPECT_GE(u, -0.1) << "at x=" << x;
        EXPECT_LE(u, 2.0) << "at x=" << x;
    }
}

TEST(GraphiteOCP, MonotonicDecrease) {
    // 石墨 OCP 总体应随 x 增大而降低
    real_t u_low = graphite_ocp(0.1);
    real_t u_high = graphite_ocp(0.9);
    EXPECT_GT(u_low, u_high);
}

TEST(NMC111OCP, RangeCheck) {
    // NMC111 OCP 在实际工作区间 (y ∈ [0.3, 0.99]) 应在 2.5~4.35V 范围内
    // 注：拟合公式在 y < 0.12 时外推超4.6V，无物理意义（实际电池不运行至此）
    for (real_t y = 0.3; y <= 0.99; y += 0.05) {
        real_t u = nmc111_ocp(y);
        EXPECT_GE(u, 2.5) << "at y=" << y;
        EXPECT_LE(u, 4.35) << "at y=" << y;
    }
}

TEST(LFPOCP, FlatPlateau) {
    // LFP 应在 ~3.4V 附近有平坦平台
    real_t u_30 = lfp_ocp(0.3);
    real_t u_50 = lfp_ocp(0.5);
    real_t u_70 = lfp_ocp(0.7);
    EXPECT_NEAR(u_50 - u_30, 0.0, 0.2);
    EXPECT_NEAR(u_70 - u_50, 0.0, 0.2);
}

// ══════════════════════════════════════════════════════════════
//  OCP 查表插值测试
// ══════════════════════════════════════════════════════════════

TEST(OcpFromLut, LinearInterpolation) {
    const std::vector<real_t> xs = {0.0, 0.5, 1.0};
    const std::vector<real_t> vs = {1.0, 2.0, 4.0};
    EXPECT_NEAR(ocp_from_lut(xs, vs, 0.0), 1.0, 1e-12);   // 下端点
    EXPECT_NEAR(ocp_from_lut(xs, vs, 0.25), 1.5, 1e-12);  // 区间中点
    EXPECT_NEAR(ocp_from_lut(xs, vs, 0.5), 2.0, 1e-12);   // 内部节点
    EXPECT_NEAR(ocp_from_lut(xs, vs, 0.75), 3.0, 1e-12);
    EXPECT_NEAR(ocp_from_lut(xs, vs, 1.0), 4.0, 1e-12);   // 上端点
}

TEST(OcpFromLut, ClampsOutsideRange) {
    const std::vector<real_t> xs = {0.1, 0.9};
    const std::vector<real_t> vs = {1.0, 3.0};
    EXPECT_NEAR(ocp_from_lut(xs, vs, -0.5), 1.0, 1e-12);  // 下界外钳位
    EXPECT_NEAR(ocp_from_lut(xs, vs, 1.5), 3.0, 1e-12);   // 上界外钳位
}

TEST(OcpFromLut, InvalidTableReturnsNaN) {
    const std::vector<real_t> xs = {0.0, 0.5};
    const std::vector<real_t> vs = {1.0};
    EXPECT_TRUE(std::isnan(ocp_from_lut(xs, vs, 0.25)));  // 长度不一致
    const std::vector<real_t> empty;
    const std::vector<real_t> v = {1.0};
    EXPECT_TRUE(std::isnan(ocp_from_lut(empty, v, 0.25))); // 空表
}

// ══════════════════════════════════════════════════════════════
//  电解液传输性质测试
// ══════════════════════════════════════════════════════════════

TEST(Electrolyte, ConductivityPositive) {
    EXPECT_GT(electrolyte_conductivity(1000.0, 298.15), 0.0);
}

TEST(Electrolyte, DiffusivityPositive) {
    EXPECT_GT(electrolyte_diffusivity(1000.0, 298.15), 0.0);
}

TEST(Electrolyte, ThermodynamicFactorPositive) {
    EXPECT_GE(thermodynamic_factor(1000.0, 298.15), 0.5);
}

// ══════════════════════════════════════════════════════════════
//  通用查表插值 + Nyman2008 电解液输运测试
// ══════════════════════════════════════════════════════════════

TEST(LutInterpolate, MatchesOcpLutImplementation) {
    // ocp_from_lut 转调 lut_interpolate：同一表结果一致
    const std::vector<real_t> xs = {0.0, 0.5, 1.0};
    const std::vector<real_t> vs = {1.0, 2.0, 4.0};
    EXPECT_NEAR(lut_interpolate(xs, vs, 0.75), 3.0, 1e-12);
    EXPECT_NEAR(lut_interpolate(xs, vs, 0.75), ocp_from_lut(xs, vs, 0.75), 1e-15);
}

TEST(Nyman2008, ReferenceValuesAt1M) {
    // PyBaMM Chen2020 标定值 @ c_e = 1000 mol/m³
    EXPECT_NEAR(nyman2008_conductivity(1000.0), 0.9487, 1e-12);
    EXPECT_NEAR(nyman2008_diffusivity(1000.0), 1.7694e-10, 1e-16);
}

TEST(Nyman2008, PhysicalShapeOverCalibrationRange) {
    // κ(0) = 0、峰值在 ~0.95 M；D(0) 最大、谷值 ~2.3 M 后回升
    EXPECT_NEAR(nyman2008_conductivity(0.0), 0.0, 1e-15);
    EXPECT_GT(nyman2008_conductivity(950.0), nyman2008_conductivity(1000.0));
    EXPECT_GT(nyman2008_diffusivity(0.0), nyman2008_diffusivity(2300.0));
    // 0~5 M 全域正值（电解液输运表覆盖范围）
    for (real_t c = 0.0; c <= 5000.0; c += 50.0) {
        EXPECT_GE(nyman2008_conductivity(c), 0.0) << "at c=" << c;
        EXPECT_GT(nyman2008_diffusivity(c), 0.0) << "at c=" << c;
    }
}
