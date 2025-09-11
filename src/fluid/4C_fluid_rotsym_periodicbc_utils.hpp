// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_FLUID_ROTSYM_PERIODICBC_UTILS_HPP
#define FOUR_C_FLUID_ROTSYM_PERIODICBC_UTILS_HPP


#include "4C_config.hpp"

#include "4C_fem_condition.hpp"
#include "4C_linalg_vector.hpp"

#include <memory>

FOUR_C_NAMESPACE_OPEN

// forward declarations
namespace Core::FE
{
  class Discretization;
}  // namespace Core::FE

namespace Core::Nodes
{
  class Node;
}

namespace FLD
{
  //! return specific component of rotated global vector result
  double get_component_of_rotated_vector_field(const int idf,  ///< the component id 0,1 (and 2)
      const Core::LinAlg::Vector<double>& proc0data,           ///< the data vector
      const int lid,         ///< local dof id (w.r.t proc0data) of component idf
      const double rotangle  ///< angle of rotation (RAD)
  );

  //! Is given node a slave node of rotationally symmetric periodic boundary conditions?
  bool is_slave_node_of_rot_sym_pbc(
      const Core::FE::Discretization& discretization, const Core::Nodes::Node* node,  ///< the node
      double& rotangle  ///< the angle of slave plane rotation (RAD)
  );

  //! Access angle of rotation and convert it to RAD
  inline double get_rot_angle_from_condition(
      const Core::Conditions::Condition* cond  ///< pointer to desired periodic boundary condition
  );

  //! Get all relevant slave nodes of rotationally symmetric periodic bc's
  void get_relevant_slave_nodes_of_rot_sym_pbc(
      std::map<int, double>&
          pbcslavenodemap,  ///< map to be filled with node gids and rotation angles
      std::shared_ptr<Core::FE::Discretization> dis);  ///< discretization

}  // namespace FLD

FOUR_C_NAMESPACE_CLOSE

#endif
