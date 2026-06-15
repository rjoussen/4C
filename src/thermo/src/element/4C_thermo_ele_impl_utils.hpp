// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_THERMO_ELE_IMPL_UTILS_HPP
#define FOUR_C_THERMO_ELE_IMPL_UTILS_HPP

/*----------------------------------------------------------------------*
 | headers                                                   dano 08/09 |
 *----------------------------------------------------------------------*/
#include "4C_config.hpp"

#include "4C_fem_general_element.hpp"
#include "4C_fem_general_utils_integration.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 |                                                           dano 08/09 |
 *----------------------------------------------------------------------*/
namespace Thermo
{
  //! Template Meta Programming version of switch over discretization type
  template <Core::FE::CellType distype>
  struct DisTypeToOptGaussRule
  {
  };
  template <>
  struct DisTypeToOptGaussRule<Core::FE::CellType::hex8>
  {
    static constexpr Core::FE::GaussRule3D rule = Core::FE::GaussRule3D::hex_8point;
  };
  template <>
  struct DisTypeToOptGaussRule<Core::FE::CellType::hex18>
  {
    static constexpr Core::FE::GaussRule3D rule = Core::FE::GaussRule3D::hex_18point;
  };
  template <>
  struct DisTypeToOptGaussRule<Core::FE::CellType::hex20>
  {
    static constexpr Core::FE::GaussRule3D rule = Core::FE::GaussRule3D::hex_27point;
  };
  template <>
  struct DisTypeToOptGaussRule<Core::FE::CellType::hex27>
  {
    static constexpr Core::FE::GaussRule3D rule = Core::FE::GaussRule3D::hex_27point;
  };
  template <>
  struct DisTypeToOptGaussRule<Core::FE::CellType::nurbs27>
  {
    static constexpr Core::FE::GaussRule3D rule = Core::FE::GaussRule3D::hex_27point;
  };
  template <>
  struct DisTypeToOptGaussRule<Core::FE::CellType::tet4>
  {
    static constexpr Core::FE::GaussRule3D rule = Core::FE::GaussRule3D::tet_4point;
  };
  template <>
  struct DisTypeToOptGaussRule<Core::FE::CellType::tet10>
  {
    static constexpr Core::FE::GaussRule3D rule = Core::FE::GaussRule3D::tet_4point;
  };
  template <>
  struct DisTypeToOptGaussRule<Core::FE::CellType::wedge6>
  {
    static constexpr Core::FE::GaussRule3D rule = Core::FE::GaussRule3D::wedge_6point;
  };
  template <>
  struct DisTypeToOptGaussRule<Core::FE::CellType::pyramid5>
  {
    static constexpr Core::FE::GaussRule3D rule = Core::FE::GaussRule3D::pyramid_8point;
  };
  template <>
  struct DisTypeToOptGaussRule<Core::FE::CellType::quad4>
  {
    static constexpr Core::FE::GaussRule2D rule = Core::FE::GaussRule2D::quad_4point;
  };
  template <>
  struct DisTypeToOptGaussRule<Core::FE::CellType::quad8>
  {
    static constexpr Core::FE::GaussRule2D rule = Core::FE::GaussRule2D::quad_9point;
  };
  template <>
  struct DisTypeToOptGaussRule<Core::FE::CellType::quad9>
  {
    static constexpr Core::FE::GaussRule2D rule = Core::FE::GaussRule2D::quad_9point;
  };
  template <>
  struct DisTypeToOptGaussRule<Core::FE::CellType::nurbs9>
  {
    static constexpr Core::FE::GaussRule2D rule = Core::FE::GaussRule2D::quad_9point;
  };
  template <>
  struct DisTypeToOptGaussRule<Core::FE::CellType::tri3>
  {
    static constexpr Core::FE::GaussRule2D rule = Core::FE::GaussRule2D::tri_3point;
  };
  template <>
  struct DisTypeToOptGaussRule<Core::FE::CellType::tri6>
  {
    static constexpr Core::FE::GaussRule2D rule = Core::FE::GaussRule2D::tri_6point;
  };
  template <>
  struct DisTypeToOptGaussRule<Core::FE::CellType::line2>
  {
    static constexpr Core::FE::GaussRule1D rule = Core::FE::GaussRule1D::line_2point;
  };
  template <>
  struct DisTypeToOptGaussRule<Core::FE::CellType::line3>
  {
    static constexpr Core::FE::GaussRule1D rule = Core::FE::GaussRule1D::line_3point;
  };

  //! Template Meta Programming version of switch over discretization type
  template <Core::FE::CellType distype>
  struct DisTypeToGaussRuleForExactSol
  {
  };
  template <>
  struct DisTypeToGaussRuleForExactSol<Core::FE::CellType::hex8>
  {
    static constexpr Core::FE::GaussRule3D rule = Core::FE::GaussRule3D::hex_27point;
  };
  template <>
  struct DisTypeToGaussRuleForExactSol<Core::FE::CellType::hex20>
  {
    static constexpr Core::FE::GaussRule3D rule = Core::FE::GaussRule3D::hex_27point;
  };
  template <>
  struct DisTypeToGaussRuleForExactSol<Core::FE::CellType::hex27>
  {
    static constexpr Core::FE::GaussRule3D rule = Core::FE::GaussRule3D::hex_27point;
  };
  template <>
  struct DisTypeToGaussRuleForExactSol<Core::FE::CellType::nurbs27>
  {
    static constexpr Core::FE::GaussRule3D rule = Core::FE::GaussRule3D::hex_27point;
  };
  template <>
  struct DisTypeToGaussRuleForExactSol<Core::FE::CellType::tet4>
  {
    static constexpr Core::FE::GaussRule3D rule = Core::FE::GaussRule3D::undefined;
  };
  template <>
  struct DisTypeToGaussRuleForExactSol<Core::FE::CellType::tet10>
  {
    static constexpr Core::FE::GaussRule3D rule = Core::FE::GaussRule3D::undefined;
  };
  template <>
  struct DisTypeToGaussRuleForExactSol<Core::FE::CellType::wedge6>
  {
    static constexpr Core::FE::GaussRule3D rule = Core::FE::GaussRule3D::undefined;
  };
  template <>
  struct DisTypeToGaussRuleForExactSol<Core::FE::CellType::pyramid5>
  {
    static constexpr Core::FE::GaussRule3D rule = Core::FE::GaussRule3D::undefined;
  };
  template <>
  struct DisTypeToGaussRuleForExactSol<Core::FE::CellType::quad4>
  {
    static constexpr Core::FE::GaussRule2D rule = Core::FE::GaussRule2D::quad_9point;
  };
  template <>
  struct DisTypeToGaussRuleForExactSol<Core::FE::CellType::quad8>
  {
    static constexpr Core::FE::GaussRule2D rule = Core::FE::GaussRule2D::quad_9point;
  };
  template <>
  struct DisTypeToGaussRuleForExactSol<Core::FE::CellType::quad9>
  {
    static constexpr Core::FE::GaussRule2D rule = Core::FE::GaussRule2D::quad_9point;
  };
  template <>
  struct DisTypeToGaussRuleForExactSol<Core::FE::CellType::tri3>
  {
    static constexpr Core::FE::GaussRule2D rule = Core::FE::GaussRule2D::undefined;
  };
  template <>
  struct DisTypeToGaussRuleForExactSol<Core::FE::CellType::tri6>
  {
    static constexpr Core::FE::GaussRule2D rule = Core::FE::GaussRule2D::undefined;
  };
  template <>
  struct DisTypeToGaussRuleForExactSol<Core::FE::CellType::line2>
  {
    static constexpr Core::FE::GaussRule1D rule = Core::FE::GaussRule1D::undefined;
  };
  template <>
  struct DisTypeToGaussRuleForExactSol<Core::FE::CellType::line3>
  {
    static constexpr Core::FE::GaussRule1D rule = Core::FE::GaussRule1D::undefined;
  };

  //! Template Meta Programming version of switch over discretization type
  template <Core::FE::CellType distype>
  struct DisTypeToNumGaussPoints
  {
  };
  template <>
  struct DisTypeToNumGaussPoints<Core::FE::CellType::hex8>
  {
    static constexpr int nquad = 8;
  };
  template <>
  struct DisTypeToNumGaussPoints<Core::FE::CellType::hex20>
  {
    static constexpr int nquad = 27;
  };
  template <>
  struct DisTypeToNumGaussPoints<Core::FE::CellType::hex27>
  {
    static constexpr int nquad = 27;
  };
  template <>
  struct DisTypeToNumGaussPoints<Core::FE::CellType::nurbs27>
  {
    static constexpr int nquad = 27;
  };
  template <>
  struct DisTypeToNumGaussPoints<Core::FE::CellType::tet4>
  {
    static constexpr int nquad = 4;
  };
  template <>
  struct DisTypeToNumGaussPoints<Core::FE::CellType::tet10>
  {
    static constexpr int nquad = 4;
  };
  template <>
  struct DisTypeToNumGaussPoints<Core::FE::CellType::wedge6>
  {
    static constexpr int nquad = 6;
  };
  template <>
  struct DisTypeToNumGaussPoints<Core::FE::CellType::pyramid5>
  {
    static constexpr int nquad = 8;
  };
  template <>
  struct DisTypeToNumGaussPoints<Core::FE::CellType::quad4>
  {
    static constexpr int nquad = 4;
  };
  template <>
  struct DisTypeToNumGaussPoints<Core::FE::CellType::quad8>
  {
    static constexpr int nquad = 9;
  };
  template <>
  struct DisTypeToNumGaussPoints<Core::FE::CellType::quad9>
  {
    static constexpr int nquad = 9;
  };
  template <>
  struct DisTypeToNumGaussPoints<Core::FE::CellType::nurbs9>
  {
    static constexpr int nquad = 9;
  };
  template <>
  struct DisTypeToNumGaussPoints<Core::FE::CellType::tri3>
  {
    static constexpr int nquad = 3;
  };
  template <>
  struct DisTypeToNumGaussPoints<Core::FE::CellType::tri6>
  {
    static constexpr int nquad = 6;
  };
  template <>
  struct DisTypeToNumGaussPoints<Core::FE::CellType::line2>
  {
    static constexpr int nquad = 2;
  };
  template <>
  struct DisTypeToNumGaussPoints<Core::FE::CellType::line3>
  {
    static constexpr int nquad = 3;
  };

  //! Template Meta Programming version of switch over discretization type
  template <Core::FE::CellType distype>
  struct DisTypeToSTRNumGaussPoints
  {
  };
  template <>
  struct DisTypeToSTRNumGaussPoints<Core::FE::CellType::hex8>
  {
    static constexpr int nquad = 8;
  };
  template <>
  struct DisTypeToSTRNumGaussPoints<Core::FE::CellType::tet4>
  {
    static constexpr int nquad = 5;
  };
  template <>
  struct DisTypeToSTRNumGaussPoints<Core::FE::CellType::tet10>
  {
    static constexpr int nquad = 11;
  };
  template <>
  struct DisTypeToSTRNumGaussPoints<Core::FE::CellType::hex27>
  {
    static constexpr int nquad = 27;
  };
  template <>
  struct DisTypeToSTRNumGaussPoints<Core::FE::CellType::hex20>
  {
    static constexpr int nquad = 27;
  };
  template <>
  struct DisTypeToSTRNumGaussPoints<Core::FE::CellType::hex18>
  {
    static constexpr int nquad = 18;
  };
  template <>
  struct DisTypeToSTRNumGaussPoints<Core::FE::CellType::nurbs27>
  {
    static constexpr int nquad = 27;
  };

}  // namespace Thermo

/*----------------------------------------------------------------------*/
FOUR_C_NAMESPACE_CLOSE

#endif
