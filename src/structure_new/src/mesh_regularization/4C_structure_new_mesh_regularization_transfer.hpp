// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_STRUCTURE_NEW_MESH_REGULARIZATION_TRANSFER_HPP
#define FOUR_C_STRUCTURE_NEW_MESH_REGULARIZATION_TRANSFER_HPP

#include "4C_config.hpp"

#include "4C_linalg_vector.hpp"

#include <vector>

FOUR_C_NAMESPACE_OPEN

namespace Core::FE
{
  class Discretization;
}

namespace Solid::MeshRegularization
{
  /*!
   * @brief Parameters of the state transfer
   */
  struct TransferParameters
  {
    //! decay c of the weights exp(-c d^2) of the history interpolation, with distances d in units
    //! of the element size
    double history_weight_decay = 10.0;
  };

  /*!
   * @brief Move the nodes of a solid discretization in the current configuration and transfer the
   * state to the moved mesh
   *
   * The current position of each node is moved by @p node_displacement, while its reference
   * position stays the initial one: the displacement is the difference of the new current position
   * and the initial position. The reference configuration therefore stays valid however large the
   * deformation is, but the elements no longer follow the material. Each moved element corrects its
   * reference configuration at its Gauss points, such that the material sees the transferred
   * deformation gradient (see Discret::Elements::corrected_deformation_gradient()).
   *
   * Further nodal fields (e.g., velocities) are interpolated at the new node positions with the
   * shape functions of the old mesh.
   *
   * The deformation gradient of the material and its history at a Gauss point of a moved element
   * are approximated with moving least squares from the Gauss points of the patch of old elements
   * around the new Gauss point: the old element containing it and the old elements with the same
   * material sharing a node with it. Elements without moved nodes keep their state exactly.
   *
   * @param discretization Solid discretization
   * @param node_displacement Movement of the nodes in the current configuration (dof row map)
   * @param displacement Displacement (dof row map), which is replaced by the displacement on the
   * moved mesh
   * @param nodal_fields Further nodal fields (dof row map) that are transferred
   * @param parameters Parameters of the transfer
   */
  void move_nodes_and_transfer_state(Core::FE::Discretization& discretization,
      const Core::LinAlg::Vector<double>& node_displacement,
      Core::LinAlg::Vector<double>& displacement,
      const std::vector<Core::LinAlg::Vector<double>*>& nodal_fields,
      const TransferParameters& parameters);

}  // namespace Solid::MeshRegularization

FOUR_C_NAMESPACE_CLOSE

#endif
