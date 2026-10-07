// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_SOLID_ELE_CALC_LIB_STRAIN_RATE_HPP
#define FOUR_C_SOLID_ELE_CALC_LIB_STRAIN_RATE_HPP

#include "4C_config.hpp"

#include "4C_fem_general_cell_type_traits.hpp"
#include "4C_linalg_fixedsizematrix.hpp"
#include "4C_linalg_symmetric_tensor.hpp"
#include "4C_linalg_tensor.hpp"
#include "4C_solid_ele_calc_lib.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Discret::Elements
{
  /*!
   * @brief Rate of the strain that a solid formulation passes to the material, together with the
   * quantities needed for its linearization
   */
  template <Core::FE::CellType celltype>
  struct StrainRateContainer
  {
    //! rate of the strain passed to the material
    Core::LinAlg::SymmetricTensor<double, Core::FE::dim<celltype>, Core::FE::dim<celltype>>
        strain_rate{};

    //! rate of the (compatible) deformation gradient
    Core::LinAlg::Tensor<double, Core::FE::dim<celltype>, Core::FE::dim<celltype>>
        deformation_gradient_rate{};
  };

  /*!
   * @brief Evaluate the rate of the compatible deformation gradient
   *
   * dF/dt = sum_n v_n (x) dN_n/dX
   *
   * @param N_XYZ derivatives of the shape functions w.r.t. the reference coordinates
   * @param nodal_velocities nodal velocities (ordered node by node)
   */
  template <Core::FE::CellType celltype>
  Core::LinAlg::Tensor<double, Core::FE::dim<celltype>, Core::FE::dim<celltype>>
  evaluate_deformation_gradient_rate(
      const std::array<Core::LinAlg::Tensor<double, Core::FE::dim<celltype>>,
          Core::FE::num_nodes(celltype)>& N_XYZ,
      const Core::LinAlg::Matrix<Core::FE::num_nodes(celltype) * Core::FE::dim<celltype>, 1>&
          nodal_velocities)
  {
    constexpr int dim = Core::FE::dim<celltype>;
    Core::LinAlg::Tensor<double, dim, dim> deformation_gradient_rate{};
    for (int node = 0; node < Core::FE::num_nodes(celltype); ++node)
      for (int i = 0; i < dim; ++i)
        for (int j = 0; j < dim; ++j)
          deformation_gradient_rate(i, j) += nodal_velocities(node * dim + i) * N_XYZ[node](j);
    return deformation_gradient_rate;
  }

  /*!
   * @brief Evaluate the symmetric part of F^T . dF/dt, i.e., the rate of the Green-Lagrange strain
   */
  template <Core::FE::CellType celltype>
  Core::LinAlg::SymmetricTensor<double, Core::FE::dim<celltype>, Core::FE::dim<celltype>>
  evaluate_green_lagrange_strain_rate(
      const Core::LinAlg::Tensor<double, Core::FE::dim<celltype>, Core::FE::dim<celltype>>& F,
      const Core::LinAlg::Tensor<double, Core::FE::dim<celltype>, Core::FE::dim<celltype>>& F_rate)
  {
    return Core::LinAlg::assume_symmetry(
        0.5 * (Core::LinAlg::transpose(F) * F_rate + Core::LinAlg::transpose(F_rate) * F));
  }
}  // namespace Discret::Elements

FOUR_C_NAMESPACE_CLOSE

#endif
