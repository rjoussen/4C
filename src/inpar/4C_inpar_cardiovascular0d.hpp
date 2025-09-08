// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_INPAR_CARDIOVASCULAR0D_HPP
#define FOUR_C_INPAR_CARDIOVASCULAR0D_HPP

#include "4C_config.hpp"

#include "4C_io_input_spec.hpp"

#include <vector>

FOUR_C_NAMESPACE_OPEN

// forward declaration
namespace Core::Conditions
{
  class ConditionDefinition;
}
namespace Inpar
{
  namespace Cardiovascular0D
  {
    /// possible 0D cardiovascular-structural solvers
    enum Cardvasc0DSolveAlgo
    {
      cardvasc0dsolve_direct,  ///< build monolithic 0D cardiovascular-structural system
      cardvasc0dsolve_block,   ///< use block preconditioner for iterative solve
    };

    enum Cardvasc0DAtriumModel
    {
      atr_prescribed,
      atr_elastance_0d,
      atr_structure_3d
    };

    enum Cardvasc0DVentricleModel
    {
      ventr_prescribed,
      ventr_elastance_0d,
      ventr_structure_3d
    };

    enum Cardvasc0DRespiratoryModel
    {
      resp_none,
      resp_standard
    };

    /// set the 0Dcardiovascular parameters
    std::vector<Core::IO::InputSpec> set_valid_parameters();

    /// set specific 0Dcardiovascular conditions
    void set_valid_conditions(std::vector<Core::Conditions::ConditionDefinition>& condlist);

  }  // namespace Cardiovascular0D
}  // namespace Inpar
/*----------------------------------------------------------------------*/
FOUR_C_NAMESPACE_CLOSE

#endif
