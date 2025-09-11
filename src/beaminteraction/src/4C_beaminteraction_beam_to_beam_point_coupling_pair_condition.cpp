// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_beaminteraction_beam_to_beam_point_coupling_pair_condition.hpp"

#include "4C_beaminteraction_beam_to_beam_point_coupling_pair.hpp"
#include "4C_fem_condition.hpp"
#include "4C_fem_discretization.hpp"
#include "4C_geometry_pair_element.hpp"

FOUR_C_NAMESPACE_OPEN


/**
 *
 */
bool BeamInteraction::BeamToBeamPointCouplingCondition::ids_in_condition(
    const int id_line, const int id_other) const
{
  if (line_ids_.find(id_line) != line_ids_.end() and line_ids_.find(id_other) != line_ids_.end())
  {
    return true;
  }
  return false;
}

/**
 *
 */
void BeamInteraction::BeamToBeamPointCouplingCondition::clear() {}

/**
 *
 */
std::shared_ptr<BeamInteraction::BeamContactPair>
BeamInteraction::BeamToBeamPointCouplingCondition::create_contact_pair(
    const std::vector<Core::Elements::Element const*>& ele_ptrs)
{
  {
    // Check if the given elements are in this condition.
    if (!ids_in_condition(ele_ptrs[0]->id(), ele_ptrs[1]->id())) return nullptr;

    // Create the beam contact pair.
    std::shared_ptr<BeamInteraction::BeamContactPair> contact_pair =
        std::make_shared<BeamToBeamPointCouplingPair<GeometryPair::t_hermite>>(
            rotational_penalty_parameter_, positional_penalty_parameter_,
            local_parameter_coordinates_);
    // Return the newly created pair.
    return contact_pair;
  }
}

/**
 *
 */
void BeamInteraction::BeamToBeamPointCouplingCondition::build_id_sets(
    const std::shared_ptr<const Core::FE::Discretization>& discretization)
{
  // Set the IDs of the nodes to be coupled
  const std::vector<int> node_ids = *(condition_line_->get_nodes());

  if (node_ids.size() != 2)
    FOUR_C_THROW(
        "The Penalty Point Coupling Condition can only handle 2 nodes per condition! If you want "
        "to couple multiple nodes, please split them into multiple conditions, each coupling two "
        "of the beam nodes.");

  std::vector<int> element_ids(node_ids.size());
  std::vector<double> position_in_parameter_space(node_ids.size());

  int i = 0;
  for (auto node_id : node_ids)
  {
    i++;

    // This means that the node is not in the column map of this proc and the element pair will thus
    // be created on a different processor
    if (!discretization->have_global_node(node_id)) return;

    Core::Nodes::Node* node = discretization->g_node(node_id);

    Core::Elements::Element* element = node->adjacent_elements()[0].user_element();
    element_ids[i - 1] = element->id();
    if (element->node_ids()[0] == node_id)
      position_in_parameter_space[i - 1] = -1;
    else
      position_in_parameter_space[i - 1] = 1;
  }

  line_ids_ = std::set<int>(element_ids.begin(), element_ids.end());
  local_parameter_coordinates_ = {position_in_parameter_space[0], position_in_parameter_space[1]};
}

FOUR_C_NAMESPACE_CLOSE
