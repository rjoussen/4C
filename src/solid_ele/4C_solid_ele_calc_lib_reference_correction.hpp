// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_SOLID_ELE_CALC_LIB_REFERENCE_CORRECTION_HPP
#define FOUR_C_SOLID_ELE_CALC_LIB_REFERENCE_CORRECTION_HPP

#include "4C_config.hpp"

#include "4C_linalg_symmetric_tensor.hpp"
#include "4C_linalg_tensor.hpp"
#include "4C_linalg_tensor_einstein.hpp"
#include "4C_linalg_tensor_generators.hpp"
#include "4C_solid_ele_calc_lib.hpp"

#include <vector>

FOUR_C_NAMESPACE_OPEN

/*!
 * @name Correction of the reference configuration of an element
 *
 * If the nodes of an element are moved without following the material (mesh regularization), the
 * deformation gradient F_e of the element with respect to its reference configuration no longer
 * describes the deformation of the material. A correction F_0 at each Gauss point maps the
 * reference configuration of the material to the one of the element, such that the deformation
 * gradient of the material is F = F_e F_0. The material is evaluated with F, while the element
 * integrates over its own reference configuration.
 */
//@{
namespace Discret::Elements
{
  //! State transferred to a Gauss point of a moved element
  struct TransferredGaussPointState
  {
    //! deformation gradient of the material
    Core::LinAlg::Tensor<double, 3, 3> deformation_gradient;
    //! history of the material (see Mat::So3Material::history_transfer_quantities())
    std::vector<double> material_history;
  };

  //! Deformation gradient of the material F = F_e F_0
  inline Core::LinAlg::Tensor<double, 3, 3> corrected_deformation_gradient(
      const Core::LinAlg::Tensor<double, 3, 3>& element_deformation_gradient,
      const Core::LinAlg::Tensor<double, 3, 3>& correction)
  {
    return element_deformation_gradient * correction;
  }

  //! Green-Lagrange strain of the material E = F_0^T E_e F_0 + (F_0^T F_0 - I) / 2
  inline Core::LinAlg::SymmetricTensor<double, 3, 3> corrected_green_lagrange_strain(
      const Core::LinAlg::SymmetricTensor<double, 3, 3>& element_strain,
      const Core::LinAlg::Tensor<double, 3, 3>& correction)
  {
    const Core::LinAlg::Tensor<double, 3, 3> correction_transpose =
        Core::LinAlg::transpose(correction);
    return Core::LinAlg::assume_symmetry(
        correction_transpose * Core::LinAlg::get_full(element_strain) * correction +
        0.5 * (correction_transpose * correction -
                  Core::LinAlg::get_full(Core::LinAlg::TensorGenerators::identity<double, 3, 3>)));
  }

  //! Ratio of the volume of the reference configuration of the element to the one of the material
  inline double reference_volume_ratio(const Core::LinAlg::Tensor<double, 3, 3>& correction)
  {
    return Core::LinAlg::det(correction);
  }

  /*!
   * @brief Transform the stress response of the material to the reference configuration of the
   * element
   *
   * The second Piola-Kirchhoff stress S and its linearization C with respect to the Green-Lagrange
   * strain of the material are pushed forward with the correction and scaled to the reference
   * volume of the element:
   *
   * S_e = F_0 S F_0^T / det F_0,  (C_e)_ijkl = (F_0)_ia (F_0)_jb (F_0)_kc (F_0)_ld C_abcd / det
   * F_0.
   *
   * S_e is the derivative of the strain energy per reference volume of the element with respect to
   * the strain of the element, and C_e its linearization.
   */
  template <Core::FE::CellType celltype>
    requires(Core::FE::dim<celltype> == 3)
  void transform_stress_to_element_reference(
      Stress<celltype>& stress, const Core::LinAlg::Tensor<double, 3, 3>& correction)
  {
    const double scale = 1.0 / reference_volume_ratio(correction);
    stress.pk2_ =
        Core::LinAlg::assume_symmetry(scale * correction * Core::LinAlg::get_full(stress.pk2_) *
                                      Core::LinAlg::transpose(correction));
    stress.cmat_ = Core::LinAlg::assume_symmetry(
        scale * Core::LinAlg::einsum<"ia", "jb", "abcd", "kc", "ld">(correction, correction,
                    Core::LinAlg::get_full(stress.cmat_), correction, correction));
  }
}  // namespace Discret::Elements
//@}

FOUR_C_NAMESPACE_CLOSE

#endif
