// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_structure_new_mesh_regularization_ghosting.hpp"

#include "4C_comm_exporter.hpp"
#include "4C_fem_discretization.hpp"
#include "4C_fem_general_element.hpp"
#include "4C_fem_general_node.hpp"
#include "4C_linalg_map.hpp"

#include <map>
#include <set>
#include <vector>

FOUR_C_NAMESPACE_OPEN

void Solid::MeshRegularization::extend_ghosting(
    Core::FE::Discretization& discretization, const int num_layers)
{
  for (int layer = 0; layer < num_layers; ++layer)
  {
    // the owner of a node knows all elements adjacent to it
    std::map<int, std::vector<int>> node_elements;
    for (int lid = 0; lid < discretization.num_my_row_nodes(); ++lid)
    {
      const Core::Nodes::Node* node = discretization.l_row_node(lid);
      std::vector<int>& elements = node_elements[node->id()];
      for (auto element : node->adjacent_elements()) elements.push_back(element.global_id());
    }

    Core::Communication::Exporter exporter(
        *discretization.node_row_map(), *discretization.node_col_map(), discretization.get_comm());
    exporter.do_export(node_elements);

    std::set<int> column_elements;
    for (int lid = 0; lid < discretization.num_my_col_elements(); ++lid)
      column_elements.insert(discretization.l_col_element(lid)->id());
    for (const auto& [node, elements] : node_elements)
      column_elements.insert(elements.begin(), elements.end());

    const std::vector<int> ids(column_elements.begin(), column_elements.end());
    const Core::LinAlg::Map element_column_map(
        -1, static_cast<int>(ids.size()), ids.data(), 0, discretization.get_comm());
    discretization.extended_ghosting(element_column_map);
  }
}

FOUR_C_NAMESPACE_CLOSE
