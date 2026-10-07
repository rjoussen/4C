// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_SOLID_ELE_CALC_FBAR_HPP
#define FOUR_C_SOLID_ELE_CALC_FBAR_HPP

#include "4C_config.hpp"

#include "4C_fem_general_cell_type_traits.hpp"
#include "4C_fem_general_element.hpp"
#include "4C_solid_ele_calc.hpp"
#include "4C_solid_ele_calc_lib.hpp"
#include "4C_solid_ele_calc_lib_fbar.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Discret::Elements
{
  template <Core::FE::CellType celltype>
  struct FBarPreparationData
  {
    /// jacobian mapping evaluated at element centroid
    JacobianMapping<celltype> jacobian_mapping_centroid;

    /// deformation gradient at element centroid
    SpatialMaterialMapping<celltype> spatial_material_mapping_centroid;
  };

  struct FBarHistoryData
  {
    // no history data needed
  };

  /*!
   * @brief A displacement based solid element formulation with FBAR element technology
   *
   * @tparam celltype
   */
  template <Core::FE::CellType celltype>
  struct FBarFormulation
  {
    static constexpr bool has_gauss_point_history = false;
    static constexpr bool has_global_history = false;
    static constexpr bool has_preparation_data = true;
    static constexpr bool has_condensed_contribution = false;

    using LinearizationContainer = FBarLinearizationContainer<celltype>;
    using StrainRateContainer = FBarStrainRateContainer<celltype>;
    using PreparationData = FBarPreparationData<celltype>;


    static FBarPreparationData<celltype> prepare(
        const Core::Elements::Element& ele, const ElementNodes<celltype>& nodal_coordinates)
    {
      const JacobianMapping<celltype> jacobian_mapping_centroid =
          evaluate_jacobian_mapping_centroid(nodal_coordinates);

      return {jacobian_mapping_centroid,
          evaluate_spatial_material_mapping(jacobian_mapping_centroid, nodal_coordinates)};
    }

    template <typename Evaluator>
    static auto evaluate(const Core::Elements::Element& ele,
        const ElementNodes<celltype>& nodal_coordinates,
        const Core::LinAlg::Tensor<double, Core::FE::dim<celltype>>& xi,
        const ShapeFunctionsAndDerivatives<celltype>& shape_functions,
        const JacobianMapping<celltype>& jacobian_mapping,
        const FBarPreparationData<celltype>& preparation_data, Evaluator evaluator)
    {
      const SpatialMaterialMapping<celltype> spatial_material_mapping =
          evaluate_spatial_material_mapping(jacobian_mapping, nodal_coordinates);

      // factor (detF0/detF)^1/3
      const double fbar_factor = evaluate_fbar_factor(
          preparation_data.spatial_material_mapping_centroid.determinant_deformation_gradient_,
          spatial_material_mapping.determinant_deformation_gradient_);

      const FBarLinearizationContainer<celltype> linearization = std::invoke(
          [&]()
          {
            FBarLinearizationContainer<celltype> linearization{};
            linearization.Bop =
                evaluate_strain_gradient(jacobian_mapping, spatial_material_mapping);

            linearization.Hop = evaluate_fbar_h_operator(jacobian_mapping.N_XYZ,
                preparation_data.jacobian_mapping_centroid.N_XYZ, spatial_material_mapping,
                preparation_data.spatial_material_mapping_centroid);

            linearization.fbar_factor = fbar_factor;

            linearization.cauchygreen = evaluate_cauchy_green(spatial_material_mapping);

            return linearization;
          });

      // deformation gradient F_bar and resulting strains: F_bar = (detF_0/detF)^1/3 F
      const SpatialMaterialMapping<celltype> spatial_material_mapping_bar =
          evaluate_spatial_material_mapping(jacobian_mapping, nodal_coordinates, fbar_factor);

      const Core::LinAlg::SymmetricTensor<double, Core::FE::dim<celltype>, Core::FE::dim<celltype>>
          cauchygreen_bar = evaluate_cauchy_green(spatial_material_mapping_bar);

      Core::LinAlg::SymmetricTensor<double, Core::FE::dim<celltype>, Core::FE::dim<celltype>>
          gl_strain_bar = evaluate_green_lagrange_strain(cauchygreen_bar);

      return evaluator(
          spatial_material_mapping_bar.deformation_gradient_, gl_strain_bar, linearization);
    }

    static inline Core::LinAlg::Matrix<9, Core::FE::num_nodes(celltype) * Core::FE::dim<celltype>>
    evaluate_d_deformation_gradient_d_displacements(const Core::Elements::Element& ele,
        const ElementNodes<celltype>& element_nodes,
        const Core::LinAlg::Tensor<double, Core::FE::dim<celltype>>& xi,
        const ShapeFunctionsAndDerivatives<celltype>& shape_functions,
        const JacobianMapping<celltype>& jacobian_mapping,
        const Core::LinAlg::Tensor<double, Internal::num_dim<celltype>,
            Internal::num_dim<celltype>>& deformation_gradient,
        const FBarPreparationData<celltype>& preparation_data)
    {
      FOUR_C_THROW(
          "This derivative of the deformation gradient w.r.t. the displacements is not "
          "implemented");
    }

    static inline Core::LinAlg::Matrix<9, Core::FE::dim<celltype>>
    evaluate_d_deformation_gradient_d_xi(const Core::Elements::Element& ele,
        const ElementNodes<celltype>& element_nodes,
        const Core::LinAlg::Tensor<double, Core::FE::dim<celltype>>& xi,
        const ShapeFunctionsAndDerivatives<celltype>& shape_functions,
        const JacobianMapping<celltype>& jacobian_mapping,
        const Core::LinAlg::Tensor<double, Internal::num_dim<celltype>,
            Internal::num_dim<celltype>>& deformation_gradient,
        const FBarPreparationData<celltype>& preparation_data)
    {
      FOUR_C_THROW("This derivative of the deformation gradient w.r.t. xi is not implemented");
    }

    static inline Core::LinAlg::Matrix<9,
        Core::FE::num_nodes(celltype) * Core::FE::dim<celltype> * Core::FE::dim<celltype>>
    evaluate_d_deformation_gradient_d_displacements_d_xi(const Core::Elements::Element& ele,
        const ElementNodes<celltype>& element_nodes,
        const Core::LinAlg::Tensor<double, Core::FE::dim<celltype>>& xi,
        const ShapeFunctionsAndDerivatives<celltype>& shape_functions,
        const JacobianMapping<celltype>& jacobian_mapping,
        const Core::LinAlg::Tensor<double, Internal::num_dim<celltype>,
            Internal::num_dim<celltype>>& deformation_gradient,
        const FBarPreparationData<celltype>& preparation_data)
    {
      FOUR_C_THROW(
          "This second derivative of the deformation gradient w.r.t. the displacements and xi is "
          "not implemented");
    }

    static void add_internal_force_vector(const JacobianMapping<celltype>& jacobian_mapping,
        const Core::LinAlg::Tensor<double, Core::FE::dim<celltype>, Core::FE::dim<celltype>>& F,
        const FBarLinearizationContainer<celltype>& linearization, const Stress<celltype>& stress,
        const double integration_factor, const FBarPreparationData<celltype>& preparation_data,
        Core::LinAlg::Matrix<Core::FE::num_nodes(celltype) * Core::FE::dim<celltype>, 1>&
            force_vector)
    {
      Discret::Elements::add_internal_force_vector(
          linearization.Bop, stress, integration_factor / linearization.fbar_factor, force_vector);
    }

    static void add_stiffness_matrix(const JacobianMapping<celltype>& jacobian_mapping,
        const Core::LinAlg::Tensor<double, Core::FE::dim<celltype>, Core::FE::dim<celltype>>& F,
        const Core::LinAlg::Tensor<double, Core::FE::dim<celltype>>& xi,
        const ShapeFunctionsAndDerivatives<celltype>& shape_functions,
        const FBarLinearizationContainer<celltype>& linearization, const Stress<celltype>& stress,
        const double integration_factor, const FBarPreparationData<celltype>& preparation_data,
        Core::LinAlg::Matrix<Core::FE::num_nodes(celltype) * Core::FE::dim<celltype>,
            Core::FE::num_nodes(celltype) * Core::FE::dim<celltype>>& stiffness_matrix)
    {
      Discret::Elements::add_elastic_stiffness_matrix(linearization.Bop, stress,
          integration_factor * linearization.fbar_factor, stiffness_matrix);
      Discret::Elements::add_geometric_stiffness_matrix(jacobian_mapping, stress.pk2_,
          integration_factor / linearization.fbar_factor, stiffness_matrix);

      // additional stiffness matrix needed for fbar method
      add_fbar_stiffness_matrix(linearization.Bop, linearization.Hop, linearization.fbar_factor,
          integration_factor, linearization.cauchygreen, stress, stiffness_matrix);
    }

    // The modified strain is E_bar = 1/2 (phi^2 C - I) with the compatible right Cauchy-Green
    // tensor C and d(phi) = phi/3 . Hop . dd.

    static FBarStrainRateContainer<celltype> evaluate_strain_rate(
        const JacobianMapping<celltype>& jacobian_mapping,
        const Core::LinAlg::Tensor<double, Core::FE::dim<celltype>, Core::FE::dim<celltype>>& F_bar,
        const FBarLinearizationContainer<celltype>& linearization,
        const Core::LinAlg::Matrix<Core::FE::num_nodes(celltype) * Core::FE::dim<celltype>, 1>&
            nodal_velocities,
        const FBarPreparationData<celltype>& preparation_data)
    {
      FBarStrainRateContainer<celltype> strain_rate{};
      const double phi = linearization.fbar_factor;
      const auto F = (1.0 / phi) * F_bar;

      strain_rate.deformation_gradient_rate =
          evaluate_deformation_gradient_rate<celltype>(jacobian_mapping.N_XYZ, nodal_velocities);
      strain_rate.fbar_h_operator_rate = linearization.Hop.dot(nodal_velocities);

      // dE_bar/dt = phi^2 (sym(F^T . dF/dt) + 1/3 (Hop . v) C)
      strain_rate.strain_rate =
          phi * phi *
          (evaluate_green_lagrange_strain_rate<celltype>(F, strain_rate.deformation_gradient_rate) +
              strain_rate.fbar_h_operator_rate / 3.0 * linearization.cauchygreen);

      // Hop . v = tr(F_0^-1 . dF_0/dt) - tr(F^-1 . dF/dt)
      const auto& jacobian_mapping_centroid = preparation_data.jacobian_mapping_centroid;
      add_d_trace_velocity_gradient_d_displacements<celltype>(
          preparation_data.spatial_material_mapping_centroid.deformation_gradient_,
          evaluate_deformation_gradient_rate<celltype>(
              jacobian_mapping_centroid.N_XYZ, nodal_velocities),
          jacobian_mapping_centroid.N_XYZ, 1.0, strain_rate.d_fbar_h_operator_rate_d_displacements);
      add_d_trace_velocity_gradient_d_displacements<celltype>(F,
          strain_rate.deformation_gradient_rate, jacobian_mapping.N_XYZ, -1.0,
          strain_rate.d_fbar_h_operator_rate_d_displacements);

      return strain_rate;
    }

    static void add_d_strain_d_displacements_contraction(
        const JacobianMapping<celltype>& jacobian_mapping,
        const Core::LinAlg::Tensor<double, Core::FE::dim<celltype>, Core::FE::dim<celltype>>& F_bar,
        const FBarLinearizationContainer<celltype>& linearization,
        const Core::LinAlg::SymmetricTensor<double, Core::FE::dim<celltype>,
            Core::FE::dim<celltype>>& Y,
        const double factor, const FBarPreparationData<celltype>& preparation_data,
        Core::LinAlg::Matrix<Core::FE::num_nodes(celltype) * Core::FE::dim<celltype>, 1>& vector)
    {
      // dE_bar = phi^2 . (sym(F^T . dF) + 1/3 . C . Hop . dd)
      const double phi = linearization.fbar_factor;
      Discret::Elements::add_internal_force_vector(
          jacobian_mapping, (1.0 / phi) * F_bar, Y, factor * phi * phi, vector);
      vector.update(factor * phi * phi / 3.0 * Core::LinAlg::ddot(linearization.cauchygreen, Y),
          linearization.Hop, 1.0);
    }

    static void add_d_strain_rate_d_displacements_at_fixed_velocities_contraction(
        const JacobianMapping<celltype>& jacobian_mapping,
        const Core::LinAlg::Tensor<double, Core::FE::dim<celltype>, Core::FE::dim<celltype>>& F_bar,
        const FBarLinearizationContainer<celltype>& linearization,
        const FBarStrainRateContainer<celltype>& strain_rate,
        const Core::LinAlg::SymmetricTensor<double, Core::FE::dim<celltype>,
            Core::FE::dim<celltype>>& Y,
        const double factor, const FBarPreparationData<celltype>& preparation_data,
        Core::LinAlg::Matrix<Core::FE::num_nodes(celltype) * Core::FE::dim<celltype>, 1>& vector)
    {
      // dE_bar/dt = phi^2 (sym(F^T . dF/dt) + 1/3 (Hop . v) C) with d(phi^2) = 2/3 phi^2 Hop . dd
      const double phi = linearization.fbar_factor;
      const double phi_squared = phi * phi;
      const auto F = (1.0 / phi) * F_bar;

      // derivative of phi^2
      vector.update(factor * 2.0 / 3.0 * Core::LinAlg::ddot(strain_rate.strain_rate, Y),
          linearization.Hop, 1.0);
      // derivative of sym(F^T . dF/dt)
      Discret::Elements::add_internal_force_vector(
          jacobian_mapping, strain_rate.deformation_gradient_rate, Y, factor * phi_squared, vector);
      // derivative of C
      Discret::Elements::add_internal_force_vector(jacobian_mapping, F, Y,
          factor * 2.0 / 3.0 * phi_squared * strain_rate.fbar_h_operator_rate, vector);
      // derivative of Hop . v
      vector.update(factor * phi_squared / 3.0 * Core::LinAlg::ddot(linearization.cauchygreen, Y),
          strain_rate.d_fbar_h_operator_rate_d_displacements, 1.0);
    }
  };

  template <Core::FE::CellType celltype>
  using FBarSolidIntegrator = SolidEleCalc<celltype, FBarFormulation<celltype>>;


}  // namespace Discret::Elements

FOUR_C_NAMESPACE_CLOSE
#endif
