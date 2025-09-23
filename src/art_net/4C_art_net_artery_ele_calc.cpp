// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_art_net_artery_ele_calc.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
Discret::Elements::ArteryEleCalc<distype>::ArteryEleCalc(
    const int numdofpernode, const std::string& disname)
    : funct_(), deriv_(), tderiv_(), xjm_(), xji_(), derxy_()
{
}
/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
double Discret::Elements::ArteryEleCalc<distype>::calculate_ele_length(Artery* ele)
{
  // get node coordinates and number of elements per node
  Core::Nodes::Node** nodes = ele->nodes();
  Core::LinAlg::Matrix<3, iel_> xyze;
  // TODO: does this work for line3?
  for (int inode = 0; inode < iel_; inode++)
  {
    const auto& x = nodes[inode]->x();
    xyze(0, inode) = x[0];
    xyze(1, inode) = x[1];
    if (x.size() > 2) xyze(2, inode) = x[2];
  }

  // Calculate the length of artery element
  const double L = sqrt(pow(xyze(0, 0) - xyze(0, 1), 2) + pow(xyze(1, 0) - xyze(1, 1), 2) +
                        pow(xyze(2, 0) - xyze(2, 1), 2));

  return L;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
// template classes

// 1D elements
template class Discret::Elements::ArteryEleCalc<Core::FE::CellType::line2>;

FOUR_C_NAMESPACE_CLOSE
