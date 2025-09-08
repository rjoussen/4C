// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_POROFLUID_PRESSURE_BASED_ELAST_INPUT_HPP
#define FOUR_C_POROFLUID_PRESSURE_BASED_ELAST_INPUT_HPP


#include "4C_config.hpp"

#include "4C_io_input_spec.hpp"

#include <map>
#include <vector>

FOUR_C_NAMESPACE_OPEN

namespace Core::Conditions
{
  class ConditionDefinition;
}  // namespace Core::Conditions

namespace PoroPressureBased
{
  /// type of coupling strategy for porofluid-elasticity problems
  enum class SolutionSchemePorofluidElast
  {
    undefined,
    twoway_partitioned,
    twoway_monolithic
  };

  //! relaxation methods for partitioned coupling
  enum class RelaxationMethods
  {
    none,
    constant,
    aitken
  };

  /// set the valid parameters for porofluid-elasticity problems
  std::vector<Core::IO::InputSpec> set_valid_parameters_porofluid_elast();

}  // namespace PoroPressureBased


FOUR_C_NAMESPACE_CLOSE

#endif
