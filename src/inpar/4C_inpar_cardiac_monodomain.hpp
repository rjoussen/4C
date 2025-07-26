// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_INPAR_CARDIAC_MONODOMAIN_HPP
#define FOUR_C_INPAR_CARDIAC_MONODOMAIN_HPP


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
namespace Inpar
{
  namespace ElectroPhysiology
  {
    /// set the elch parameters
    void set_valid_parameters(std::map<std::string, Core::IO::InputSpec>& list);

    /// set specific elch conditions
    void set_valid_conditions(std::vector<Core::Conditions::ConditionDefinition>& condlist);
  }  // namespace ElectroPhysiology
}  // namespace Inpar
FOUR_C_NAMESPACE_CLOSE

#endif
