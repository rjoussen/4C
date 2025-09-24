// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_fem_discretization_utils.hpp"

#include "4C_fem_discretization.hpp"
#include "4C_fem_general_element.hpp"
#include "4C_fem_general_node.hpp"
#include "4C_linalg_map.hpp"
#include "4C_utils_function.hpp"
#include "4C_utils_function_manager.hpp"

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Core::FE::evaluate_initial_field(const Core::Utils::FunctionManager& function_manager,
    const Core::FE::Discretization& discret, const std::string& fieldstring,
    Core::LinAlg::Vector<double>& fieldvector, const std::vector<int>& locids)
{
  // get initial field conditions
  std::vector<const Core::Conditions::Condition*> initfieldconditions;
  discret.get_condition("Initfield", initfieldconditions);

  //--------------------------------------------------------
  // loop through Initfield conditions and evaluate them
  //--------------------------------------------------------
  // Note that this method does not sum up but 'sets' values in fieldvector.
  // For this reason, Initfield BCs are evaluated hierarchical meaning
  // in this order (just like Dirichlet BCs):
  //                VolumeInitfield
  //                SurfaceInitfield
  //                LineInitfield
  //                PointInitfield
  // This way, lower entities override higher ones.
  const std::vector<Core::Conditions::ConditionType> evaluation_type_order = {
      Core::Conditions::VolumeInitfield, Core::Conditions::SurfaceInitfield,
      Core::Conditions::LineInitfield, Core::Conditions::PointInitfield};

  for (const auto& type : evaluation_type_order)
  {
    for (const auto& initfieldcondition : initfieldconditions)
    {
      if (initfieldcondition->type() != type) continue;
      const std::string condstring = initfieldcondition->parameters().get<std::string>("FIELD");
      if (condstring != fieldstring) continue;
      do_initial_field(function_manager, discret, *initfieldcondition, fieldvector, locids);
    }
  }
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Core::FE::do_initial_field(const Core::Utils::FunctionManager& function_manager,
    const Core::FE::Discretization& discret, const Core::Conditions::Condition& cond,
    Core::LinAlg::Vector<double>& fieldvector, const std::vector<int>& locids)
{
  const std::vector<int> cond_nodeids = *cond.get_nodes();
  if (cond_nodeids.empty()) FOUR_C_THROW("Initfield condition does not have nodal cloud.");

  // loop nodes to identify and evaluate spatial distributions
  // of Initfield boundary conditions
  const auto funct_num = cond.parameters().get<int>("FUNCT");

  for (const int cond_nodeid : cond_nodeids)
  {
    // do only nodes in my row map
    int cond_node_lid = discret.node_row_map()->lid(cond_nodeid);
    if (cond_node_lid < 0) continue;
    Core::Nodes::Node* node = discret.l_row_node(cond_node_lid);

    // call explicitly the main dofset, i.e. the first column
    std::vector<int> node_dofs = discret.dof(0, node);
    const int total_numdof = static_cast<int>(node_dofs.size());

    // Get native number of dofs at this node. There might be multiple dofsets
    // (in xfem cases), thus the size of the dofs vector might be a multiple
    // of this.
    auto elements = node->adjacent_elements();
    auto ele_with_max_dof = std::ranges::max_element(elements,
        [&](ElementRef a, ElementRef b)
        {
          return a.user_element()->num_dof_per_node(*node) <
                 b.user_element()->num_dof_per_node(*node);
        });
    const int numdof = ele_with_max_dof->user_element()->num_dof_per_node(*node);

    if ((total_numdof % numdof) != 0) FOUR_C_THROW("illegal dof set number");

    // now loop over all relevant DOFs
    for (int j = 0; j < total_numdof; ++j)
    {
      int localdof = j % numdof;

      // evaluate function if local DOF id exists
      // in the given locids vector
      for (const int locid : locids)
      {
        if (localdof == locid)
        {
          const double time = 0.0;  // dummy time here

          const double functfac =
              funct_num > 0
                  ? function_manager.function_by_id<Core::Utils::FunctionOfSpaceTime>(funct_num)
                        .evaluate(node->x(), time, localdof)
                  : 0.0;

          // assign value
          const int gid = node_dofs[j];
          const int lid = fieldvector.get_map().lid(gid);
          if (lid < 0) FOUR_C_THROW("Global id {} not on this proc in system vector", gid);
          fieldvector.get_values()[lid] = functfac;
        }
      }
    }
  }
}

FOUR_C_NAMESPACE_CLOSE
