// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <gtest/gtest.h>

#include "4C_ale_mesh_quality.hpp"

#include <array>
#include <cmath>
#include <vector>

namespace
{
  using namespace FourC;

  //! unit cube with node numbering of a hex8 element
  std::vector<double> unit_cube()
  {
    return {0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1};
  }

  //! a distorted hexahedron
  std::vector<double> distorted_hex()
  {
    std::vector<double> x = unit_cube();
    const std::vector<double> perturbation{0.1, -0.05, 0.02, 0.3, 0.1, -0.1, -0.1, 0.2, 0.05, 0.0,
        0.1, -0.02, 0.2, 0.0, 0.1, 0.1, -0.1, 0.3, 0.0, 0.05, -0.1, -0.05, 0.1, 0.2};
    for (std::size_t i = 0; i < x.size(); ++i) x[i] += perturbation[i];
    return x;
  }

  TEST(AleMeshQualityTest, IdealElementsHaveUnitQuality)
  {
    EXPECT_NEAR(ALE::mean_ratio_quality(Core::FE::CellType::hex8, unit_cube()), 1.0, 1.0e-14);

    // scaled and rotated cube
    std::vector<double> cube = unit_cube();
    for (std::size_t node = 0; node < 8; ++node)
    {
      const double x = cube[3 * node];
      const double y = cube[3 * node + 1];
      cube[3 * node] = 2.0 * (x - y);
      cube[3 * node + 1] = 2.0 * (x + y);
      cube[3 * node + 2] *= 2.0 * std::sqrt(2.0);
    }
    EXPECT_NEAR(ALE::mean_ratio_quality(Core::FE::CellType::hex8, cube), 1.0, 1.0e-14);

    const double s3 = std::sqrt(3.0);
    const std::vector<double> regular_tet{
        0, 0, 0, 1, 0, 0, 0.5, s3 / 2, 0, 0.5, s3 / 6, std::sqrt(2.0 / 3.0)};
    EXPECT_NEAR(ALE::mean_ratio_quality(Core::FE::CellType::tet4, regular_tet), 1.0, 1.0e-14);
  }

  TEST(AleMeshQualityTest, DistortedAndInvertedElements)
  {
    const double quality = ALE::mean_ratio_quality(Core::FE::CellType::hex8, distorted_hex());
    EXPECT_GT(quality, 0.0);
    EXPECT_LT(quality, 1.0);

    // move node 6 through the opposite face
    std::vector<double> inverted = unit_cube();
    inverted[3 * 6 + 2] = -0.5;
    EXPECT_LT(ALE::mean_ratio_quality(Core::FE::CellType::hex8, inverted), 0.0);
  }

  TEST(AleMeshQualityTest, PotentialVanishesForUndeformedIdealElement)
  {
    const std::vector<double> zero(24, 0.0);
    Core::LinAlg::SerialDenseVector residual;
    const double energy = ALE::evaluate_mesh_quality_potential(
        Core::FE::CellType::hex8, unit_cube(), zero, {}, &residual, nullptr);
    EXPECT_NEAR(energy, 0.0, 1.0e-14);
    for (int i = 0; i < 24; ++i) EXPECT_NEAR(residual(i), 0.0, 1.0e-13);
  }

  TEST(AleMeshQualityTest, DerivativesMatchFiniteDifferences)
  {
    const std::vector<double> reference = distorted_hex();
    std::vector<double> displacement(24);
    for (int i = 0; i < 24; ++i) displacement[i] = 0.01 * std::sin(1.0 + i);

    // the corners grow and shrink slightly, which also checks the size term for growth only
    for (bool growth_only : {false, true})
    {
      const ALE::MeshQualityParameters parameters{.size_weight = 0.7,
          .motion_weight = 0.3,
          .size_growth_only = growth_only,
          .shape_exponent = growth_only ? 2.5 : 1.0};
      Core::LinAlg::SerialDenseVector residual;
      Core::LinAlg::SerialDenseMatrix stiffness;
      ALE::evaluate_mesh_quality_potential(
          Core::FE::CellType::hex8, reference, displacement, parameters, &residual, &stiffness);

      const double h = 1.0e-6;
      for (int j = 0; j < 24; ++j)
      {
        std::vector<double> plus = displacement;
        std::vector<double> minus = displacement;
        plus[j] += h;
        minus[j] -= h;

        Core::LinAlg::SerialDenseVector residual_plus;
        Core::LinAlg::SerialDenseVector residual_minus;
        const double energy_plus = ALE::evaluate_mesh_quality_potential(
            Core::FE::CellType::hex8, reference, plus, parameters, &residual_plus, nullptr);
        const double energy_minus = ALE::evaluate_mesh_quality_potential(
            Core::FE::CellType::hex8, reference, minus, parameters, &residual_minus, nullptr);

        EXPECT_NEAR(residual(j), (energy_plus - energy_minus) / (2.0 * h), 1.0e-7);
        for (int i = 0; i < 24; ++i)
          EXPECT_NEAR(stiffness(i, j), (residual_plus(i) - residual_minus(i)) / (2.0 * h), 1.0e-6);
      }
    }
  }

  TEST(AleMeshQualityTest, TargetConfigurationIsOptimal)
  {
    // a flat element, which is far from the ideal cube
    std::vector<double> flat = unit_cube();
    for (int i = 0; i < 8; ++i) flat[3 * i + 2] *= 0.1;
    const std::vector<double> zero(24, 0.0);

    // relative to the ideal cube, the flat element is distorted
    EXPECT_GT(ALE::evaluate_mesh_quality_potential(
                  Core::FE::CellType::hex8, flat, zero, {}, nullptr, nullptr),
        0.1);

    // with itself as target, it is optimal: no energy and no driving force
    Core::LinAlg::SerialDenseVector residual;
    EXPECT_NEAR(ALE::evaluate_mesh_quality_potential(
                    Core::FE::CellType::hex8, flat, zero, {}, &residual, nullptr, flat),
        0.0, 1.0e-14);
    for (int i = 0; i < 24; ++i) EXPECT_NEAR(residual(i), 0.0, 1.0e-12);

    // a rotated and translated copy of the target is optimal as well
    // displacement of a rotation by 90 degrees about z and a translation
    std::vector<double> moved(24);
    for (int i = 0; i < 8; ++i)
    {
      moved[3 * i + 0] = -flat[3 * i + 1] + 2.0 - flat[3 * i + 0];
      moved[3 * i + 1] = flat[3 * i + 0] + 1.0 - flat[3 * i + 1];
      moved[3 * i + 2] = 0.0;
    }
    EXPECT_NEAR(ALE::evaluate_mesh_quality_potential(
                    Core::FE::CellType::hex8, flat, moved, {}, nullptr, nullptr, flat),
        0.0, 1.0e-12);
  }

  TEST(AleMeshQualityTest, ElementSizeIsTheVolumeOfAffineElements)
  {
    // a parallelepiped spanned by (2, 0, 0), (0.5, 1, 0), (0, 0.3, 0.5)
    std::vector<double> x(24);
    const std::array<std::array<double, 3>, 3> edges{{{2, 0, 0}, {0.5, 1, 0}, {0, 0.3, 0.5}}};
    const std::vector<double> cube = unit_cube();
    for (int n = 0; n < 8; ++n)
      for (int d = 0; d < 3; ++d)
        x[3 * n + d] = cube[3 * n] * edges[0][d] + cube[3 * n + 1] * edges[1][d] +
                       cube[3 * n + 2] * edges[2][d];
    EXPECT_NEAR(ALE::element_size(Core::FE::CellType::hex8, x), 1.0, 1.0e-14);

    const std::vector<double> tet{0, 0, 0, 2, 0, 0, 0, 3, 0, 0, 0, 1};
    EXPECT_NEAR(ALE::element_size(Core::FE::CellType::tet4, tet), 1.0, 1.0e-14);
  }

  TEST(AleMeshQualityTest, MotionTermMeasuresDisplacementInElementSizes)
  {
    // a translation does not change the shape or size, only the motion term contributes
    const ALE::MeshQualityParameters parameters{.motion_weight = 0.4};
    const double d = 0.1;

    // cube with edge length 2, with itself as target and without target
    std::vector<double> cube = unit_cube();
    for (double& x : cube) x *= 2.0;
    std::vector<double> translation(24, 0.0);
    for (int i = 0; i < 8; ++i) translation[3 * i] = d;
    EXPECT_NEAR(ALE::evaluate_mesh_quality_potential(Core::FE::CellType::hex8, cube, translation,
                    parameters, nullptr, nullptr, cube),
        0.4 * std::pow(d / 2.0, 2), 1.0e-14);
    EXPECT_NEAR(ALE::evaluate_mesh_quality_potential(
                    Core::FE::CellType::hex8, cube, translation, parameters, nullptr, nullptr),
        0.4 * std::pow(d / 2.0, 2), 1.0e-14);

    // regular tetrahedron with edge length 2
    const double s3 = std::sqrt(3.0);
    std::vector<double> tet{0, 0, 0, 1, 0, 0, 0.5, s3 / 2, 0, 0.5, s3 / 6, std::sqrt(2.0 / 3.0)};
    for (double& x : tet) x *= 2.0;
    std::vector<double> tet_translation(12, 0.0);
    for (int i = 0; i < 4; ++i) tet_translation[3 * i + 1] = d;
    EXPECT_NEAR(ALE::evaluate_mesh_quality_potential(
                    Core::FE::CellType::tet4, tet, tet_translation, parameters, nullptr, nullptr),
        0.4 * std::pow(d / 2.0, 2), 1.0e-14);
  }

  TEST(AleMeshQualityTest, PotentialOfElementInvertedInTheReferenceIsNaN)
  {
    std::vector<double> reference = unit_cube();
    reference[3 * 6 + 2] = -0.2;

    std::vector<double> displacement(24, 0.0);
    displacement[3 * 6 + 2] = 1.2;
    EXPECT_TRUE(std::isnan(ALE::evaluate_mesh_quality_potential(
        Core::FE::CellType::hex8, reference, displacement, {}, nullptr, nullptr)));
  }

  TEST(AleMeshQualityTest, SizeTermForGrowthOnly)
  {
    // uniform scaling of the unit cube: no distortion, only the size term with v = s^3
    const ALE::MeshQualityParameters two_sided{.size_weight = 2.0};
    const ALE::MeshQualityParameters growth_only{.size_weight = 2.0, .size_growth_only = true};
    for (double s : {0.8, 1.2})
    {
      std::vector<double> displacement = unit_cube();
      for (double& u : displacement) u *= s - 1.0;
      const double v = s * s * s;
      EXPECT_NEAR(ALE::evaluate_mesh_quality_potential(Core::FE::CellType::hex8, unit_cube(),
                      displacement, two_sided, nullptr, nullptr),
          2.0 * (v + 1.0 / v - 2.0), 1.0e-12);
      EXPECT_NEAR(ALE::evaluate_mesh_quality_potential(Core::FE::CellType::hex8, unit_cube(),
                      displacement, growth_only, nullptr, nullptr),
          s > 1.0 ? 2.0 * std::pow(v - 1.0, 3) / (v * v) : 0.0, 1.0e-12);
    }
  }

  TEST(AleMeshQualityTest, ShapeExponentWeightsStrongDistortionsMore)
  {
    const std::vector<double> zero(24, 0.0);
    const ALE::MeshQualityParameters plain{};
    const ALE::MeshQualityParameters squared{.shape_exponent = 2.0};
    // ideal elements have no energy for any exponent
    EXPECT_NEAR(ALE::evaluate_mesh_quality_potential(
                    Core::FE::CellType::hex8, unit_cube(), zero, squared, nullptr, nullptr),
        0.0, 1.0e-14);
    // with mu - 1 = d at all corners of a uniformly sheared cube, (mu^2 - 1) / 2 = d + d^2 / 2
    std::vector<double> shear(24, 0.0);
    for (int node = 0; node < 8; ++node) shear[3 * node] = 0.5 * unit_cube()[3 * node + 1];
    const double d = ALE::evaluate_mesh_quality_potential(
        Core::FE::CellType::hex8, unit_cube(), shear, plain, nullptr, nullptr);
    EXPECT_GT(d, 0.0);
    EXPECT_NEAR(ALE::evaluate_mesh_quality_potential(
                    Core::FE::CellType::hex8, unit_cube(), shear, squared, nullptr, nullptr),
        d + 0.5 * d * d, 1.0e-12);
  }

  TEST(AleMeshQualityTest, PotentialOfInvertedElementIsNaN)
  {
    std::vector<double> displacement(24, 0.0);
    displacement[3 * 6 + 2] = -1.5;
    EXPECT_TRUE(std::isnan(ALE::evaluate_mesh_quality_potential(
        Core::FE::CellType::hex8, unit_cube(), displacement, {}, nullptr, nullptr)));
  }
}  // namespace
