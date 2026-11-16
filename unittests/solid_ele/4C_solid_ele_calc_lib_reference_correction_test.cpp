// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <gtest/gtest.h>

#include "4C_solid_ele_calc_lib_reference_correction.hpp"

#include "4C_linalg_tensor_generators.hpp"
#include "4C_unittest_utils_assertions_test.hpp"

namespace
{
  using namespace FourC;

  constexpr Core::FE::CellType celltype = Core::FE::CellType::hex8;
  constexpr double lambda = 1.3;
  constexpr double mu = 0.7;

  Core::LinAlg::Tensor<double, 3, 3> correction()
  {
    return Core::LinAlg::Tensor<double, 3, 3>{{
        {1.1, 0.2, -0.1},
        {0.05, 0.9, 0.3},
        {-0.2, 0.1, 1.2},
    }};
  }

  Core::LinAlg::SymmetricTensor<double, 3, 3> element_strain()
  {
    return Core::LinAlg::assume_symmetry(Core::LinAlg::Tensor<double, 3, 3>{{
        {0.02, 0.01, -0.03},
        {0.01, -0.04, 0.02},
        {-0.03, 0.02, 0.05},
    }});
  }

  //! St. Venant-Kirchhoff strain energy, quadratic in the strain
  double strain_energy(const Core::LinAlg::SymmetricTensor<double, 3, 3>& E)
  {
    const double trace = Core::LinAlg::trace(E);
    return 0.5 * lambda * trace * trace + mu * Core::LinAlg::ddot(E, E);
  }

  Discret::Elements::Stress<celltype> stress_response(
      const Core::LinAlg::SymmetricTensor<double, 3, 3>& E)
  {
    const auto I = Core::LinAlg::get_full(Core::LinAlg::TensorGenerators::identity<double, 3, 3>);
    Discret::Elements::Stress<celltype> stress;
    stress.pk2_ = Core::LinAlg::assume_symmetry(
        lambda * Core::LinAlg::trace(E) * I + 2.0 * mu * Core::LinAlg::get_full(E));
    stress.cmat_ = Core::LinAlg::assume_symmetry(
        lambda * Core::LinAlg::dyadic(I, I) +
        2.0 * mu *
            Core::LinAlg::get_full(
                Core::LinAlg::TensorGenerators::symmetric_identity<double, 3, 3, 3, 3>));
    return stress;
  }

  //! strain energy per reference volume of the element as a function of the strain of the element
  double element_strain_energy(const Core::LinAlg::SymmetricTensor<double, 3, 3>& element_E)
  {
    return strain_energy(
               Discret::Elements::corrected_green_lagrange_strain(element_E, correction())) /
           Discret::Elements::reference_volume_ratio(correction());
  }

  Discret::Elements::Stress<celltype> element_stress(
      const Core::LinAlg::SymmetricTensor<double, 3, 3>& element_E)
  {
    Discret::Elements::Stress<celltype> stress = stress_response(
        Discret::Elements::corrected_green_lagrange_strain(element_E, correction()));
    Discret::Elements::transform_stress_to_element_reference(stress, correction());
    return stress;
  }

  //! symmetric perturbation of the strain such that the change of a function f of the strain is
  //! h df/dE_ij
  Core::LinAlg::SymmetricTensor<double, 3, 3> perturbation(int i, int j, double h)
  {
    Core::LinAlg::Tensor<double, 3, 3> dE{};
    dE(i, j) += 0.5 * h;
    dE(j, i) += 0.5 * h;
    return Core::LinAlg::assume_symmetry(dE);
  }

  TEST(ReferenceCorrection, StrainIsConsistentWithDeformationGradient)
  {
    const Core::LinAlg::Tensor<double, 3, 3> element_F{{
        {1.05, 0.1, 0.0},
        {-0.02, 0.97, 0.04},
        {0.03, 0.0, 1.02},
    }};
    const auto I = Core::LinAlg::get_full(Core::LinAlg::TensorGenerators::identity<double, 3, 3>);
    const auto element_E =
        Core::LinAlg::assume_symmetry(0.5 * (Core::LinAlg::transpose(element_F) * element_F - I));

    const auto F = Discret::Elements::corrected_deformation_gradient(element_F, correction());
    const auto E = Core::LinAlg::assume_symmetry(0.5 * (Core::LinAlg::transpose(F) * F - I));

    FOUR_C_EXPECT_NEAR(
        Discret::Elements::corrected_green_lagrange_strain(element_E, correction()), E, 1.0e-14);
  }

  TEST(ReferenceCorrection, StressIsDerivativeOfElementStrainEnergy)
  {
    const Discret::Elements::Stress<celltype> stress = element_stress(element_strain());

    // the energy is quadratic in the strain, so central differences are exact up to round-off
    constexpr double h = 1.0e-3;
    for (int i = 0; i < 3; ++i)
    {
      for (int j = 0; j < 3; ++j)
      {
        const double derivative =
            (element_strain_energy(element_strain() + perturbation(i, j, h)) -
                element_strain_energy(element_strain() - perturbation(i, j, h))) /
            (2.0 * h);
        EXPECT_NEAR(stress.pk2_(i, j), derivative, 1.0e-10);
      }
    }
  }

  TEST(ReferenceCorrection, LinearizationIsDerivativeOfStress)
  {
    const Discret::Elements::Stress<celltype> stress = element_stress(element_strain());

    constexpr double h = 1.0e-3;
    for (int k = 0; k < 3; ++k)
    {
      for (int l = 0; l < 3; ++l)
      {
        const auto stress_plus = element_stress(element_strain() + perturbation(k, l, h));
        const auto stress_minus = element_stress(element_strain() - perturbation(k, l, h));
        for (int i = 0; i < 3; ++i)
        {
          for (int j = 0; j < 3; ++j)
          {
            EXPECT_NEAR(stress.cmat_(i, j, k, l),
                (stress_plus.pk2_(i, j) - stress_minus.pk2_(i, j)) / (2.0 * h), 1.0e-10);
          }
        }
      }
    }
  }

  TEST(ReferenceCorrection, IdentityKeepsTheStress)
  {
    Discret::Elements::Stress<celltype> stress = stress_response(element_strain());
    const Discret::Elements::Stress<celltype> reference = stress;
    Discret::Elements::transform_stress_to_element_reference(
        stress, Core::LinAlg::get_full(Core::LinAlg::TensorGenerators::identity<double, 3, 3>));

    FOUR_C_EXPECT_NEAR(stress.pk2_, reference.pk2_, 1.0e-15);
    FOUR_C_EXPECT_NEAR(stress.cmat_, reference.cmat_, 1.0e-15);
  }
}  // namespace
