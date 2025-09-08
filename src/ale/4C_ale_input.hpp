// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_ALE_INPUT_HPP
#define FOUR_C_ALE_INPUT_HPP

#include "4C_config.hpp"

#include "4C_io_input_spec.hpp"

#include <vector>

FOUR_C_NAMESPACE_OPEN

// forward declaration
namespace Core::Conditions
{
  class ConditionDefinition;
}

/*----------------------------------------------------------------------*/

namespace ALE
{
  /// possible types of moving boundary simulation
  enum AleDynamic
  {
    none,              ///< default
    laplace_material,  ///< Laplacian smoothing based on material configuration
    laplace_spatial,   ///< Laplacian smoothing based on spatial configuration
    springs_material,  ///< use a spring analogy based on material configuration
    springs_spatial,   ///< use a spring analogy based on spatial configuration
    solid,             ///< nonlinear pseudo-structure approach
    solid_linear       ///< linear pseudo-structure approach
  };
  /// Handling of non-converged nonlinear solver
  enum DivContAct
  {
    divcont_stop,     ///< abort simulation
    divcont_continue  ///< continue nevertheless
  };

  /// mesh tying and mesh sliding algorithm
  enum MeshTying
  {
    no_meshtying,
    meshtying,
    meshsliding
  };

  /// initial displacement field
  enum InitialDisp
  {
    initdisp_zero_disp,        ///< zero initial displacement
    initdisp_disp_by_function  ///< initial displacement from function
  };

  /// Defines all valid parameters for ale problem
  Core::IO::InputSpec set_valid_parameters();

  /// Defines ale specific conditions
  void set_valid_conditions(std::vector<Core::Conditions::ConditionDefinition>& condlist);

}  // namespace ALE

/*----------------------------------------------------------------------*/
FOUR_C_NAMESPACE_CLOSE

#endif
