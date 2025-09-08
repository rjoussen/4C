// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_SSI_INPUT_HPP
#define FOUR_C_SSI_INPUT_HPP

#include "4C_config.hpp"

#include "4C_io_input_spec.hpp"

#include <vector>

FOUR_C_NAMESPACE_OPEN

namespace Core::Conditions
{
  class ConditionDefinition;
}

namespace SSI
{
  /// Type of coupling strategy for SSI problems
  enum class SolutionSchemeOverFields
  {
    ssi_OneWay_ScatraToSolid,
    ssi_OneWay_SolidToScatra,
    // ssi_SequStagg_ScatraToSolid,
    // ssi_SequStagg_SolidToScatra,
    ssi_IterStagg,
    ssi_IterStaggFixedRel_ScatraToSolid,
    ssi_IterStaggFixedRel_SolidToScatra,
    ssi_IterStaggAitken_ScatraToSolid,
    ssi_IterStaggAitken_SolidToScatra,
    // IterStaggAitkenIrons,
    ssi_Monolithic
  };

  /// Type of coupling strategy between the two fields of the SSI problems
  enum class FieldCoupling
  {
    volume_match,
    volume_nonmatch,
    boundary_nonmatch,
    volumeboundary_match
  };

  //! type of scalar transport time integration
  enum class ScaTraTimIntType
  {
    standard,
    cardiac_monodomain,
    elch
  };

  /// set the ssi parameters
  std::vector<Core::IO::InputSpec> set_valid_parameters();

  /// set specific ssi conditions
  void set_valid_conditions(std::vector<Core::Conditions::ConditionDefinition>& condlist);

}  // namespace SSI


FOUR_C_NAMESPACE_CLOSE

#endif
