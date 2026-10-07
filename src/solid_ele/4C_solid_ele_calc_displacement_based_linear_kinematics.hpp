// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_SOLID_ELE_CALC_DISPLACEMENT_BASED_LINEAR_KINEMATICS_HPP
#define FOUR_C_SOLID_ELE_CALC_DISPLACEMENT_BASED_LINEAR_KINEMATICS_HPP

#include "4C_config.hpp"

#include "4C_fem_general_cell_type.hpp"
#include "4C_fem_general_cell_type_traits.hpp"
#include "4C_fem_general_element.hpp"
#include "4C_solid_ele_calc.hpp"
#include "4C_solid_ele_calc_lib.hpp"
#include "4C_solid_ele_calc_lib_strain_rate.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Discret::Elements
{
  /*!
   * @brief A container holding the linearization of the displacement based linear kinematics solid
   * element formulation
   */
  template <Core::FE::CellType celltype>
  struct DisplacementBasedLinearKinematicsLinearizationContainer
  {
    // nothing special needed
  };

  /*!
   * @brief A displacement based solid element formulation with linear kinematics (i.e. small
   * displacements)
   *
   * @tparam celltype
   */
  template <Core::FE::CellType celltype>
  struct DisplacementBasedLinearKinematicsFormulation
  {
    static constexpr bool has_gauss_point_history = false;
    static constexpr bool has_global_history = false;
    static constexpr bool has_preparation_data = false;
    static constexpr bool has_condensed_contribution = false;


    using LinearizationContainer =
        DisplacementBasedLinearKinematicsLinearizationContainer<celltype>;
    using StrainRateContainer = Discret::Elements::StrainRateContainer<celltype>;

    template <typename Evaluator>
    static inline auto evaluate(const Core::Elements::Element& ele,
        const ElementNodes<celltype>& nodal_coordinates,
        const Core::LinAlg::Tensor<double, Core::FE::dim<celltype>>& xi,
        const ShapeFunctionsAndDerivatives<celltype>& shape_functions,
        const JacobianMapping<celltype>& jacobian_mapping, Evaluator evaluator)
    {
      const DisplacementBasedLinearKinematicsLinearizationContainer<celltype> linearization{};

      Core::LinAlg::SymmetricTensor<double, Core::FE::dim<celltype>, Core::FE::dim<celltype>>
          gl_strain = evaluate_linear_gl_strain(nodal_coordinates, jacobian_mapping);

      return evaluator(Core::LinAlg::get_full(Core::LinAlg::TensorGenerators::identity<double,
                           Core::FE::dim<celltype>, Core::FE::dim<celltype>>),
          gl_strain, linearization);
    }

    static inline Core::LinAlg::Matrix<9, Core::FE::num_nodes(celltype) * Core::FE::dim<celltype>>
    evaluate_d_deformation_gradient_d_displacements(const Core::Elements::Element& ele,
        const ElementNodes<celltype>& element_nodes,
        const Core::LinAlg::Tensor<double, Core::FE::dim<celltype>>& xi,
        const ShapeFunctionsAndDerivatives<celltype>& shape_functions,
        const JacobianMapping<celltype>& jacobian_mapping,
        const Core::LinAlg::Tensor<double, Internal::num_dim<celltype>,
            Internal::num_dim<celltype>>& deformation_gradient)
    {
      // linearization is zero for small displacements
      return Core::LinAlg::Matrix<9, Core::FE::num_nodes(celltype) * Core::FE::dim<celltype>>(
          Core::LinAlg::Initialization::zero);
    }

    static inline Core::LinAlg::Matrix<9, Core::FE::dim<celltype>>
    evaluate_d_deformation_gradient_d_xi(const Core::Elements::Element& ele,
        const ElementNodes<celltype>& element_nodes,
        const Core::LinAlg::Tensor<double, Core::FE::dim<celltype>>& xi,
        const ShapeFunctionsAndDerivatives<celltype>& shape_functions,
        const JacobianMapping<celltype>& jacobian_mapping,
        const Core::LinAlg::Tensor<double, Internal::num_dim<celltype>,
            Internal::num_dim<celltype>>& deformation_gradient)
    {
      // linearization is zero for small displacements
      return Core::LinAlg::Matrix<9, Core::FE::dim<celltype>>(Core::LinAlg::Initialization::zero);
    }

    static inline Core::LinAlg::Matrix<9,
        Core::FE::num_nodes(celltype) * Core::FE::dim<celltype> * Core::FE::dim<celltype>>
    evaluate_d_deformation_gradient_d_displacements_d_xi(const Core::Elements::Element& ele,
        const ElementNodes<celltype>& element_nodes,
        const Core::LinAlg::Tensor<double, Core::FE::dim<celltype>>& xi,
        const ShapeFunctionsAndDerivatives<celltype>& shape_functions,
        const JacobianMapping<celltype>& jacobian_mapping,
        const Core::LinAlg::Tensor<double, Internal::num_dim<celltype>,
            Internal::num_dim<celltype>>& deformation_gradient)
    {
      // linearization is zero for small displacements
      return Core::LinAlg::Matrix<9,
          Core::FE::num_nodes(celltype) * Core::FE::dim<celltype> * Core::FE::dim<celltype>>(
          Core::LinAlg::Initialization::zero);
    }

    static void add_internal_force_vector(const JacobianMapping<celltype>& jacobian_mapping,
        const Core::LinAlg::Tensor<double, Core::FE::dim<celltype>, Core::FE::dim<celltype>>& F,
        const DisplacementBasedLinearKinematicsLinearizationContainer<celltype>& linearization,
        const Stress<celltype>& stress, const double integration_factor,
        Core::LinAlg::Matrix<Core::FE::num_nodes(celltype) * Core::FE::dim<celltype>, 1>&
            force_vector)
    {
      Discret::Elements::add_internal_force_vector(
          jacobian_mapping, F, stress.pk2_, integration_factor, force_vector);
    }

    static void add_stiffness_matrix(const JacobianMapping<celltype>& jacobian_mapping,
        const Core::LinAlg::Tensor<double, Core::FE::dim<celltype>, Core::FE::dim<celltype>>& F,
        const Core::LinAlg::Tensor<double, Core::FE::dim<celltype>>& xi,
        const ShapeFunctionsAndDerivatives<celltype>& shape_functions,
        const DisplacementBasedLinearKinematicsLinearizationContainer<celltype>& linearization,
        const Stress<celltype>& stress, const double integration_factor,
        Core::LinAlg::Matrix<Core::FE::num_nodes(celltype) * Core::FE::dim<celltype>,
            Core::FE::num_nodes(celltype) * Core::FE::dim<celltype>>& stiffness_matrix)
    {
      for (std::size_t i = 0; i < Core::FE::num_nodes(celltype); ++i)
      {
        const auto dNdXi_C = (integration_factor * jacobian_mapping.N_XYZ[i]) * stress.cmat_;

        for (std::size_t j = 0; j < Core::FE::num_nodes(celltype); ++j)
        {
          add_nodal_contribution<celltype>(
              i, j, dNdXi_C * jacobian_mapping.N_XYZ[j], stiffness_matrix);
        }
      }
    }

    static StrainRateContainer evaluate_strain_rate(
        const JacobianMapping<celltype>& jacobian_mapping,
        const Core::LinAlg::Tensor<double, Core::FE::dim<celltype>, Core::FE::dim<celltype>>& F,
        const DisplacementBasedLinearKinematicsLinearizationContainer<celltype>& linearization,
        const Core::LinAlg::Matrix<Core::FE::num_nodes(celltype) * Core::FE::dim<celltype>, 1>&
            nodal_velocities)
    {
      StrainRateContainer strain_rate{};
      strain_rate.deformation_gradient_rate =
          evaluate_deformation_gradient_rate<celltype>(jacobian_mapping.N_XYZ, nodal_velocities);
      // F = I, i.e., dE/dt = sym(dF/dt)
      strain_rate.strain_rate =
          evaluate_green_lagrange_strain_rate<celltype>(F, strain_rate.deformation_gradient_rate);
      return strain_rate;
    }

    static void add_d_strain_d_displacements_contraction(
        const JacobianMapping<celltype>& jacobian_mapping,
        const Core::LinAlg::Tensor<double, Core::FE::dim<celltype>, Core::FE::dim<celltype>>& F,
        const DisplacementBasedLinearKinematicsLinearizationContainer<celltype>& linearization,
        const Core::LinAlg::SymmetricTensor<double, Core::FE::dim<celltype>,
            Core::FE::dim<celltype>>& Y,
        const double factor,
        Core::LinAlg::Matrix<Core::FE::num_nodes(celltype) * Core::FE::dim<celltype>, 1>& vector)
    {
      // F = I, i.e., dE = sym(dF)
      Discret::Elements::add_internal_force_vector(jacobian_mapping, F, Y, factor, vector);
    }

    static void add_d_strain_rate_d_displacements_at_fixed_velocities_contraction(
        const JacobianMapping<celltype>& jacobian_mapping,
        const Core::LinAlg::Tensor<double, Core::FE::dim<celltype>, Core::FE::dim<celltype>>& F,
        const DisplacementBasedLinearKinematicsLinearizationContainer<celltype>& linearization,
        const StrainRateContainer& strain_rate,
        const Core::LinAlg::SymmetricTensor<double, Core::FE::dim<celltype>,
            Core::FE::dim<celltype>>& Y,
        const double factor,
        Core::LinAlg::Matrix<Core::FE::num_nodes(celltype) * Core::FE::dim<celltype>, 1>& vector)
    {
      // the strain rate is linear in the velocities and independent of the displacements
    }
  };

  template <Core::FE::CellType celltype>
  using DisplacementBasedLinearKinematicsSolidIntegrator =
      SolidEleCalc<celltype, DisplacementBasedLinearKinematicsFormulation<celltype>>;


}  // namespace Discret::Elements

FOUR_C_NAMESPACE_CLOSE
#endif