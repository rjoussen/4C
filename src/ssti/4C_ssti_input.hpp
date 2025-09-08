// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_SSTI_INPUT_HPP
#define FOUR_C_SSTI_INPUT_HPP

#include "4C_config.hpp"

#include "4C_io_input_spec.hpp"

#include <vector>

FOUR_C_NAMESPACE_OPEN

namespace Core::Conditions
{
  class ConditionDefinition;
}


namespace SSTI
{
  /// Type of coupling strategy for SSI problems
  enum class SolutionScheme
  {
    monolithic
  };

  //! type of scalar transport time integration
  enum class ScaTraTimIntType
  {
    elch
  };

  /// set the ssti parameters
  std::vector<Core::IO::InputSpec> set_valid_parameters();

  /// set specific ssti conditions
  void set_valid_conditions(std::vector<Core::Conditions::ConditionDefinition>& condlist);

}  // namespace SSTI

FOUR_C_NAMESPACE_CLOSE

#endif
