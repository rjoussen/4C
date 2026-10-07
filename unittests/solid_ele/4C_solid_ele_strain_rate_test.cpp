// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <gtest/gtest.h>

#include "4C_fem_general_cell_type.hpp"
#include "4C_solid_ele.hpp"
#include "4C_solid_ele_calc_displacement_based.hpp"
#include "4C_solid_ele_calc_displacement_based_linear_kinematics.hpp"
#include "4C_solid_ele_calc_fbar.hpp"
#include "4C_solid_ele_calc_lib.hpp"
#include "4C_solid_ele_formulation.hpp"

#include <array>

namespace
{
  using namespace FourC;

  constexpr Core::FE::CellType celltype = Core::FE::CellType::hex8;
  constexpr int num_nodes = Core::FE::num_nodes(celltype);
  constexpr int num_dim = Core::FE::dim<celltype>;
  constexpr int num_dof = num_nodes * num_dim;

  using NodalVector = Core::LinAlg::Matrix<num_dof, 1>;

  //! strain passed to the material at a point of an element with nodal displacements d
  template <typename SolidFormulation>
  struct StrainEvaluator
  {
    Discret::Elements::ElementNodes<celltype> nodes;
    Core::LinAlg::Tensor<double, num_dim> xi;
    Discret::Elements::Solid<3> element{1, 0};

    void set_displacements(const NodalVector& d)
    {
      for (int node = 0; node < Core::FE::num_nodes(celltype); ++node)
      {
        for (int i = 0; i < num_dim; ++i)
        {
          nodes.displacements(i, node) = d(node * num_dim + i);
          nodes.current_coordinates(i, node) =
              nodes.reference_coordinates(i, node) + nodes.displacements(i, node);
        }
      }
    }

    //! call @p action with the strain, the strain rate and the operators of the formulation
    template <typename Action>
    void evaluate(const NodalVector& d, const NodalVector& v, Action action)
    {
      set_displacements(d);
      Discret::Elements::SolidFormulationHistory<SolidFormulation> history_data{};
      const auto preparation_data =
          Discret::Elements::prepare<SolidFormulation>(element, nodes, history_data);
      const auto shape_functions =
          Discret::Elements::evaluate_shape_functions_and_derivs<celltype>(xi, nodes);
      const auto jacobian_mapping =
          Discret::Elements::evaluate_jacobian_mapping(shape_functions, nodes);

      Discret::Elements::evaluate<SolidFormulation>(element, nodes, xi, shape_functions,
          jacobian_mapping, preparation_data, history_data, 0,
          [&](const auto& F, const auto& strain, const auto& linearization)
          {
            const auto strain_rate = Discret::Elements::evaluate_strain_rate<SolidFormulation>(
                jacobian_mapping, F, linearization, v, preparation_data, history_data, 0);
            action(strain, strain_rate,
                [&](const auto& Y, const double timefac_d, const double factor,
                    NodalVector& d_strain_d_displacements,
                    NodalVector& d_strain_rate_d_displacements)
                {
                  Discret::Elements::add_d_strain_d_displacements_contraction<SolidFormulation>(
                      jacobian_mapping, F, linearization, Y, factor, preparation_data, history_data,
                      0, d_strain_d_displacements);
                  Discret::Elements::add_d_strain_rate_d_displacements_contraction<
                      SolidFormulation>(jacobian_mapping, F, linearization, strain_rate, Y,
                      timefac_d, factor, preparation_data, history_data, 0,
                      d_strain_rate_d_displacements);
                });
          });
    }

    Core::LinAlg::SymmetricTensor<double, num_dim, num_dim> strain(const NodalVector& d)
    {
      Core::LinAlg::SymmetricTensor<double, num_dim, num_dim> result{};
      evaluate(d, NodalVector(Core::LinAlg::Initialization::zero),
          [&](const auto& strain, const auto&, const auto&) { result = strain; });
      return result;
    }

    Core::LinAlg::SymmetricTensor<double, num_dim, num_dim> strain_rate(
        const NodalVector& d, const NodalVector& v)
    {
      Core::LinAlg::SymmetricTensor<double, num_dim, num_dim> result{};
      evaluate(d, v, [&](const auto&, const auto& strain_rate, const auto&)
          { result = strain_rate.strain_rate; });
      return result;
    }
  };

  template <typename SolidFormulation>
  class SolidEleStrainRateTest : public ::testing::Test
  {
   protected:
    void SetUp() override
    {
      evaluator_.xi = {{0.1, -0.2, 0.3}};

      // distorted reference configuration of the element
      const std::array<std::array<double, num_dim>, num_nodes> reference_coordinates = {{
          {-1.02, -0.97, -1.04},
          {1.03, -1.01, -0.96},
          {0.98, 1.04, -1.02},
          {-0.96, 0.99, -0.98},
          {-1.01, -1.03, 1.02},
          {0.97, -0.98, 1.03},
          {1.04, 1.02, 0.97},
          {-0.99, 0.96, 1.01},
      }};
      // a finite deformation, i.e., the volume change varies within the element
      const std::array<std::array<double, num_dim>, num_nodes> displacements = {{
          {0.03, -0.02, 0.04},
          {-0.01, 0.05, -0.03},
          {0.04, 0.01, 0.02},
          {-0.05, -0.04, 0.01},
          {0.02, 0.03, -0.05},
          {-0.03, 0.02, 0.04},
          {0.01, -0.05, -0.02},
          {0.05, 0.04, 0.03},
      }};
      const std::array<std::array<double, num_dim>, num_nodes> velocities = {{
          {0.31, -0.12, 0.45},
          {-0.27, 0.08, -0.39},
          {0.14, 0.42, -0.06},
          {-0.48, -0.21, 0.33},
          {0.05, -0.36, 0.19},
          {0.37, 0.24, -0.41},
          {-0.16, 0.49, 0.02},
          {0.28, -0.03, -0.25},
      }};
      for (int node = 0; node < num_nodes; ++node)
      {
        for (int i = 0; i < num_dim; ++i)
        {
          evaluator_.nodes.reference_coordinates(i, node) = reference_coordinates[node][i];
          displacements_(node * num_dim + i) = displacements[node][i];
          velocities_(node * num_dim + i) = velocities[node][i];
        }
      }
    }

    StrainEvaluator<SolidFormulation> evaluator_;
    NodalVector displacements_{Core::LinAlg::Initialization::zero};
    NodalVector velocities_{Core::LinAlg::Initialization::zero};
    const Core::LinAlg::SymmetricTensor<double, num_dim, num_dim> Y_ =
        Core::LinAlg::assume_symmetry(Core::LinAlg::Tensor<double, num_dim, num_dim>{
            {{1.1, 0.4, -0.3}, {0.4, -0.7, 0.2}, {-0.3, 0.2, 0.9}}});
    const double timefac_d_ = 7.0;
    const double perturbation_ = 1.0e-6;
    const double tolerance_ = 1.0e-7;
  };

  using Formulations = ::testing::Types<Discret::Elements::DisplacementBasedFormulation<celltype>,
      Discret::Elements::DisplacementBasedLinearKinematicsFormulation<celltype>,
      Discret::Elements::FBarFormulation<celltype>>;
  TYPED_TEST_SUITE(SolidEleStrainRateTest, Formulations);

  TYPED_TEST(SolidEleStrainRateTest, StrainRateMatchesFiniteDifferences)
  {
    auto& evaluator = this->evaluator_;
    NodalVector d_plus(this->displacements_);
    NodalVector d_minus(this->displacements_);
    d_plus.update(this->perturbation_, this->velocities_, 1.0);
    d_minus.update(-this->perturbation_, this->velocities_, 1.0);
    const auto fd_strain_rate = (1.0 / (2.0 * this->perturbation_)) *
                                (evaluator.strain(d_plus) - evaluator.strain(d_minus));

    const auto strain_rate = evaluator.strain_rate(this->displacements_, this->velocities_);
    for (int i = 0; i < num_dim; ++i)
      for (int j = 0; j < num_dim; ++j)
        EXPECT_NEAR(strain_rate(i, j), fd_strain_rate(i, j), this->tolerance_);
  }

  TYPED_TEST(SolidEleStrainRateTest, LinearizationsMatchFiniteDifferences)
  {
    auto& evaluator = this->evaluator_;
    const auto& Y = this->Y_;
    const double timefac_d = this->timefac_d_;

    NodalVector d_strain_d_displacements(Core::LinAlg::Initialization::zero);
    NodalVector d_strain_rate_d_displacements(Core::LinAlg::Initialization::zero);
    evaluator.evaluate(this->displacements_, this->velocities_,
        [&](const auto&, const auto&, const auto& add_contractions)
        {
          add_contractions(
              Y, timefac_d, 1.0, d_strain_d_displacements, d_strain_rate_d_displacements);
        });

    for (int k = 0; k < num_dof; ++k)
    {
      const auto strain_contraction = [&](const double delta)
      {
        NodalVector d(this->displacements_);
        d(k) += delta;
        return Core::LinAlg::ddot(evaluator.strain(d), Y);
      };
      // the velocities depend on the displacements with dv/dd = timefac_d
      const auto strain_rate_contraction = [&](const double delta)
      {
        NodalVector d(this->displacements_);
        NodalVector v(this->velocities_);
        d(k) += delta;
        v(k) += timefac_d * delta;
        return Core::LinAlg::ddot(evaluator.strain_rate(d, v), Y);
      };

      const double h = this->perturbation_;
      EXPECT_NEAR(d_strain_d_displacements(k),
          (strain_contraction(h) - strain_contraction(-h)) / (2.0 * h), this->tolerance_)
          << "dof " << k;
      EXPECT_NEAR(d_strain_rate_d_displacements(k),
          (strain_rate_contraction(h) - strain_rate_contraction(-h)) / (2.0 * h),
          timefac_d * this->tolerance_)
          << "dof " << k;
    }
  }
}  // namespace
