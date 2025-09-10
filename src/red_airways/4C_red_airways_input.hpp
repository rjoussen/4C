// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_RED_AIRWAYS_INPUT_HPP
#define FOUR_C_RED_AIRWAYS_INPUT_HPP

#include "4C_config.hpp"

#include "4C_io_input_spec.hpp"

#include <vector>

FOUR_C_NAMESPACE_OPEN

// forward declaration
namespace Core::Conditions
{
  class ConditionDefinition;
}  // namespace Core::Conditions


namespace Airway
{
  /// types of integration schemes
  enum class RedAirwaysDyntype
  {
    OneStepTheta,
  };

  /// types of solver
  enum class RedAirwaysSolvertype
  {
    Linear,
    Nonlinear,
  };

  /// reduced airways parameters
  Core::IO::InputSpec valid_parameters();

  /// set specific reduced airways conditions
  void set_valid_conditions(std::vector<Core::Conditions::ConditionDefinition>& condlist);
}  // namespace Airway

FOUR_C_NAMESPACE_CLOSE

#endif