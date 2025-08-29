// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_STI_INPUT_HPP
#define FOUR_C_STI_INPUT_HPP

#include "4C_config.hpp"

#include "4C_io_input_spec.hpp"

#include <map>
#include <vector>

FOUR_C_NAMESPACE_OPEN

// forward declaration

namespace Core::Conditions
{
  class ConditionDefinition;
}


namespace STI
{
  //! type of coupling between scatra and thermo fields
  enum class CouplingType
  {
    undefined,
    monolithic,
    oneway_scatratothermo,
    oneway_thermotoscatra,
    twoway_scatratothermo,
    twoway_scatratothermo_aitken,
    twoway_scatratothermo_aitken_dofsplit,
    twoway_thermotoscatra,
    twoway_thermotoscatra_aitken
  };

  //! type of scalar transport time integration
  enum class ScaTraTimIntType
  {
    standard,
    elch
  };

  //! set valid parameters for scatra-thermo interaction
  void set_valid_parameters(std::map<std::string, Core::IO::InputSpec>& list);

  //! set valid conditions for scatra-thermo interaction
  void set_valid_conditions(std::vector<Core::Conditions::ConditionDefinition>& condlist);
}  // namespace STI

FOUR_C_NAMESPACE_CLOSE

#endif
