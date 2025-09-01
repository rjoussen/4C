// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_IO_ELEMENT_VTK_CELL_TYPE_REGISTER_HPP
#define FOUR_C_IO_ELEMENT_VTK_CELL_TYPE_REGISTER_HPP

#include "4C_config.hpp"

#include "4C_fem_general_element.hpp"
#include "4C_utils_exceptions.hpp"

#include <vector>

FOUR_C_NAMESPACE_OPEN

namespace Core::IO
{
  std::pair<uint8_t, std::vector<int>> inline get_vtk_cell_type_from_element_cell_type(
      Core::FE::CellType four_c_ele_shape_type)
  {
    static_assert(29 == static_cast<int>(Core::FE::CellType::max_distype),
        "The number of element types defined by Core::FE::CellType does not match "
        "the number of vtk cell types supported by the vtu writer.");

    // the VTK element types are from the documentation of vtkCellType,
    // e.g. at http://www.vtk.org/doc/nightly/html/vtkCellType_8h.html
    // this list must be kept in sync with the element types since we use this for index translation
    switch (four_c_ele_shape_type)
    {
      case Core::FE::CellType::dis_none:
        return {0, {}};
      case Core::FE::CellType::quad4:
        return {9, {0, 1, 2, 3}};
      case Core::FE::CellType::quad6:
        return {30, {0, 1, 4, 3, 2, 5}};
      case Core::FE::CellType::quad8:
        return {23, {0, 1, 2, 3, 4, 5, 6, 7}};
      case Core::FE::CellType::quad9:
        return {28, {0, 1, 2, 3, 4, 5, 6, 7, 8}};
      case Core::FE::CellType::tri3:
        return {5, {0, 1, 2}};
      case Core::FE::CellType::tri6:
        return {22, {0, 1, 2, 3, 4, 5}};
      case Core::FE::CellType::hex8:
        return {12, {0, 1, 2, 3, 4, 5, 6, 7}};
      case Core::FE::CellType::hex16:
        return {12, {0, 1, 2, 3, 8, 9, 10, 11}};
      case Core::FE::CellType::hex18:
        return {12, {0, 1, 2, 3, 9, 10, 11, 12}};
      case Core::FE::CellType::hex20:
        return {25, {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 16, 17, 18, 19, 12, 13, 14, 15}};
      case Core::FE::CellType::hex27:
        return {29, {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 16, 17, 18, 19, 12, 13, 14, 15, 24, 22,
                        21, 23, 20, 25, 26}};
      case Core::FE::CellType::tet4:
        return {10, {0, 1, 2, 3}};
      case Core::FE::CellType::tet10:
        return {24, {0, 1, 2, 3, 4, 5, 6, 7, 8, 9}};
      case Core::FE::CellType::wedge6:
        return {13, {0, 2, 1, 3, 5, 4}};
      case Core::FE::CellType::wedge15:
        return {26, {0, 2, 1, 3, 5, 4, 8, 7, 6, 14, 13, 12, 9, 11, 10}};
      case Core::FE::CellType::pyramid5:
        return {14, {0, 1, 2, 3, 4}};
      case Core::FE::CellType::line2:
        return {3, {0, 1}};
      case Core::FE::CellType::line3:
        return {21, {0, 1, 2}};
      case Core::FE::CellType::line4:
      case Core::FE::CellType::line5:  // line5 -> mapped onto line4
      case Core::FE::CellType::line6:  // line6 -> mapped onto line4
        return {35, {0, 1, 2, 3}};
      case Core::FE::CellType::point1:
        return {1, {0}};
      default:
        FOUR_C_THROW("VTK cell type not implemented for element: {}",
            Core::FE::cell_type_to_string(four_c_ele_shape_type));
    }
  }
}  // namespace Core::IO

FOUR_C_NAMESPACE_CLOSE

#endif
