// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_fem_general_element_center.hpp"

#include "4C_fem_discretization.hpp"
#include "4C_fem_general_element.hpp"
#include "4C_fem_general_node.hpp"

FOUR_C_NAMESPACE_OPEN

std::vector<double> Core::FE::element_center_refe_coords(const Core::Elements::Element& ele)
{
  // get nodes of element
  const int numnodes = ele.num_node();
  const double invnumnodes = 1.0 / numnodes;

  // calculate mean of node coordinates
  std::vector<double> centercoords(3, 0.0);
  for (auto node : ele.node_range())
  {
    const auto x = node.x();
    for (size_t i = 0; i < x.size(); ++i)
    {
      centercoords[i] += x[i];
    }
  }
  for (double& centercoord : centercoords)
  {
    centercoord *= invnumnodes;
  }

  return centercoords;
}

FOUR_C_NAMESPACE_CLOSE
