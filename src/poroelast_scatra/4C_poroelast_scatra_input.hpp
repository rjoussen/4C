// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_POROELAST_SCATRA_INPUT_HPP
#define FOUR_C_POROELAST_SCATRA_INPUT_HPP

#include "4C_config.hpp"

#include "4C_io_input_spec.hpp"


FOUR_C_NAMESPACE_OPEN


namespace PoroElastScaTra
{
  /// Type of coupling strategy for poro scatra problems
  enum SolutionSchemeOverFields
  {
    Monolithic,
    Part_ScatraToPoro,
    Part_PoroToScatra,
    Part_TwoWay
  };

  /// poroscatra parameters
  Core::IO::InputSpec valid_parameters();

}  // namespace PoroElastScaTra


FOUR_C_NAMESPACE_CLOSE

#endif
