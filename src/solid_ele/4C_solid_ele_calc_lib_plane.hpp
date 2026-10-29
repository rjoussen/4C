// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_SOLID_ELE_CALC_LIB_PLANE_HPP
#define FOUR_C_SOLID_ELE_CALC_LIB_PLANE_HPP

#include "4C_config.hpp"

#include "4C_fem_general_cell_type.hpp"
#include "4C_solid_ele_calc_lib.hpp"

#include <Teuchos_ParameterList.hpp>

FOUR_C_NAMESPACE_OPEN

namespace Discret::Elements
{
  //! 3D tensor with the in-plane components of @p tensor_2d and the out-of-plane diagonal
  //! component @p diagonal_value
  Core::LinAlg::Tensor<double, 3, 3> make_3d_tensor(
      const Core::LinAlg::Tensor<double, 2, 2>& tensor_2d, double diagonal_value = 0.0);

  //! 3D tensor with the in-plane components of @p tensor_2d and the out-of-plane diagonal
  //! component @p diagonal_value
  Core::LinAlg::SymmetricTensor<double, 3, 3> make_3d_tensor(
      const Core::LinAlg::SymmetricTensor<double, 2, 2>& tensor_2d, double diagonal_value = 0.0);

  //! in-plane components of @p tensor_3d
  Core::LinAlg::SymmetricTensor<double, 2, 2> extract_2d_part(
      const Core::LinAlg::SymmetricTensor<double, 3, 3>& tensor_3d);

  //! in-plane components of @p tensor_3d
  Core::LinAlg::SymmetricTensor<double, 2, 2, 2, 2> extract_2d_part(
      const Core::LinAlg::SymmetricTensor<double, 3, 3, 3, 3>& tensor_3d);

  template <Core::FE::CellType celltype>
    requires(Core::FE::dim<celltype> == 2)
  void transform_to_3d(Mat::So3Material& material,
      const ElementProperties<celltype>& element_properties,
      const Core::LinAlg::Tensor<double, 2, 2>& defgrd,
      const Core::LinAlg::SymmetricTensor<double, 2, 2>& gl_strain, Teuchos::ParameterList& params,
      const Mat::EvaluationContext<2>& context, int gp, int eleGID,
      const std::function<void(const Core::LinAlg::Tensor<double, 3, 3>&,
          const Core::LinAlg::SymmetricTensor<double, 3, 3>&, const Mat::EvaluationContext<3>&)>&
          funct);

  template <Core::FE::CellType celltype>
    requires(Core::FE::dim<celltype> == 2)
  Stress<celltype> evaluate_material_stress(Mat::So3Material& material,
      const ElementProperties<celltype>& element_properties,
      const Core::LinAlg::Tensor<double, 2, 2>& defgrd,
      const Core::LinAlg::SymmetricTensor<double, 2, 2>& gl_strain, Teuchos::ParameterList& params,
      const Mat::EvaluationContext<2>& context, const int gp, const int eleGID);

  template <Core::FE::CellType celltype>
    requires(Core::FE::dim<celltype> == 2)
  void update_material(Mat::So3Material& material,
      const ElementProperties<celltype>& element_properties,
      const Core::LinAlg::Tensor<double, 2, 2>& defgrd, Teuchos::ParameterList& params,
      const Mat::EvaluationContext<2>& context, const int gp, const int eleGID);

  template <Core::FE::CellType celltype>
    requires(Core::FE::dim<celltype> == 2)
  [[nodiscard]] double evaluate_material_strain_energy(Mat::So3Material& material,
      const ElementProperties<celltype>& element_properties,
      const Core::LinAlg::SymmetricTensor<double, 2, 2>& gl_strain, Teuchos::ParameterList& params,
      const Mat::EvaluationContext<2>& context, const int gp, const int eleGID);

}  // namespace Discret::Elements

FOUR_C_NAMESPACE_CLOSE

#endif