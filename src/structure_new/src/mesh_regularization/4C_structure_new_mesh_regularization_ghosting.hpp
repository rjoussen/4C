// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_STRUCTURE_NEW_MESH_REGULARIZATION_GHOSTING_HPP
#define FOUR_C_STRUCTURE_NEW_MESH_REGULARIZATION_GHOSTING_HPP

#include "4C_config.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Core::FE
{
  class Discretization;
}

namespace Solid::MeshRegularization
{
  /*!
   * @brief Number of layers of ghosted elements beyond the elements of the row nodes
   *
   * The owner of an old element evaluates the transferred state at a new Gauss point in it from
   * the patch of old elements sharing a node with it, which is known with one more layer. The
   * movement of the nodes is not limited by the ghosting: the faces and old elements of other
   * processors are searched globally and evaluated by their owners.
   */
  constexpr int num_ghost_layers = 1;

  /*!
   * @brief Extend the ghosting of a discretization by layers of elements
   *
   * Each layer adds the elements that share a node with the column elements.
   *
   * @note Has to be called before the state vectors are set up.
   */
  void extend_ghosting(Core::FE::Discretization& discretization, int num_layers);

}  // namespace Solid::MeshRegularization

FOUR_C_NAMESPACE_CLOSE

#endif
