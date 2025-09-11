// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_beaminteraction_beam_to_beam_contact_tangentsmoothing.hpp"

#include "4C_fem_discretization.hpp"
#include "4C_fem_general_element.hpp"
#include "4C_fem_general_node.hpp"
#include "4C_fem_general_utils_fem_shapefunctions.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
BeamInteraction::B3CNeighbor::B3CNeighbor(const Core::Elements::Element* left_neighbor,
    const Core::Elements::Element* right_neighbor, int connecting_node_left,
    int connecting_node_right)
    : left_neighbor_(left_neighbor),
      right_neighbor_(right_neighbor),
      connecting_node_left_(connecting_node_left),
      connecting_node_right_(connecting_node_right)
{
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
std::shared_ptr<BeamInteraction::B3CNeighbor>
BeamInteraction::Beam3TangentSmoothing::determine_neighbors(const Core::Elements::Element* element1)
{
  const Core::Elements::Element* left_neighbor = nullptr;
  const Core::Elements::Element* right_neighbor = nullptr;
  int connecting_node_left = 0;
  int connecting_node_right = 0;

  // number of nodes element1
  int numnode = element1->num_node();

  // n_right is the local node-ID of the elements right node (at xi = 1) whereas the elements left
  // node (at xi = -1) always has the local ID 1
  int n_right;
  if (numnode == 2)
    n_right = 1;
  else
    n_right = numnode - 2;

  int globalneighborId = 0;
  int globalnodeId = 0;

  // local node 1 of element1 --> xi(element1)=-1 --> left neighbor
  globalnodeId = *element1->node_ids();

  // only one neighbor element on each side of the considered element is allowed
  if ((**(element1->nodes())).num_element() > 2)
  {
    FOUR_C_THROW(
        "ERROR: The implemented smoothing routine does not work for more than 2 adjacent Elements "
        "per node");
  }

  // Loop over all elements adjacent to the left node of element1
  for (auto neighbor_element : element1->nodes()[0]->adjacent_elements())
  {
    if (neighbor_element.global_id() != element1->id())
    {
      left_neighbor = neighbor_element.user_element();

      // Determine the connecting node
      connecting_node_left =
          (neighbor_element.user_element()->node_ids()[0] == globalnodeId) ? 0 : n_right;
    }
  }

  // ocal node n_right of element1 --> xi(element1)=1 --> right neighbor
  globalnodeId = *(element1->node_ids() + n_right);

  // only one neighbor element on each side of the considered element is allowed
  if ((**(element1->nodes() + n_right)).num_element() > 2)
  {
    FOUR_C_THROW(
        "ERROR: The implemented smoothing routine does not work for more than 2 adjacent Elements "
        "per node");
  }

  // loop over all elements adjacent to the right node of element1
  const auto rightNode = element1->nodes()[n_right];
  for (auto neighbor : rightNode->adjacent_elements())
  {
    globalneighborId = neighbor.global_id();

    if (globalneighborId != element1->id())
    {
      right_neighbor = neighbor.user_element();

      // if node 0 of the right neighbor is the connecting node:
      if (neighbor.user_element()->node_ids()[0] == globalnodeId)
      {
        connecting_node_right = 0;
      }
      else
      {
        connecting_node_right = n_right;
      }
    }
  }

  return std::make_shared<BeamInteraction::B3CNeighbor>(
      left_neighbor, right_neighbor, connecting_node_left, connecting_node_right);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int BeamInteraction::Beam3TangentSmoothing::get_boundary_node(const int nnode)
{
  if (nnode == 2)
    return 1;
  else
    return nnode - 2;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
double BeamInteraction::Beam3TangentSmoothing::get_ele_length(
    const Core::LinAlg::SerialDenseMatrix& elepos, const int nright)
{
  double length = 0.0;
  for (int i = 0; i < 3; i++)
    length += (elepos(i, nright) - elepos(i, 0)) * (elepos(i, nright) - elepos(i, 0));

  return sqrt(length);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Core::LinAlg::SerialDenseMatrix BeamInteraction::Beam3TangentSmoothing::get_nodal_derivatives(
    const int node, const int nnode, const double length, const Core::FE::CellType distype)
{
  Core::LinAlg::SerialDenseMatrix deriv1(1, nnode);

  if (node == nnode)
    Core::FE::shape_function_1d_deriv1(deriv1, -1.0 + 2.0 / (nnode - 1), distype);
  else
  {
    if (node == 1)
      Core::FE::shape_function_1d_deriv1(deriv1, -1.0, distype);
    else
      Core::FE::shape_function_1d_deriv1(deriv1, -1.0 + node * 2.0 / (nnode - 1), distype);
  }

  for (int i = 0; i < nnode; i++) deriv1(0, i) = 2.0 * deriv1(0, i) / length;

  return deriv1;
}

FOUR_C_NAMESPACE_CLOSE
