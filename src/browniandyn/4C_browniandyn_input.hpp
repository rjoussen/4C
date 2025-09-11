// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_BROWNIANDYN_INPUT_HPP
#define FOUR_C_BROWNIANDYN_INPUT_HPP

#include "4C_config.hpp"

#include "4C_io_input_spec.hpp"
#include "4C_utils_parameter_list.fwd.hpp"


FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
namespace BrownianDynamics
{
  /// the way how damping coefficient values for beams are specified
  enum BeamDampingCoefficientSpecificationType
  {
    cylinder_geometry_approx,
    input_file,
    vague
  };

  /// brownian dynamic parameters
  Core::IO::InputSpec valid_parameters();

}  // namespace BrownianDynamics

/*----------------------------------------------------------------------*/
FOUR_C_NAMESPACE_CLOSE

#endif
