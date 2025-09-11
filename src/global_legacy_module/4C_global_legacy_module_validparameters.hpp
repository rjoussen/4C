// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_GLOBAL_LEGACY_MODULE_VALIDPARAMETERS_HPP
#define FOUR_C_GLOBAL_LEGACY_MODULE_VALIDPARAMETERS_HPP

#include "4C_config.hpp"

#include "4C_io_input_spec.hpp"

#include <iostream>
#include <memory>
#include <string>
#include <vector>

FOUR_C_NAMESPACE_OPEN

namespace Global
{
  /**
   * The valid parameters as a vector of InputSpec.
   */
  std::vector<Core::IO::InputSpec> valid_parameters();

}  // namespace Global

FOUR_C_NAMESPACE_CLOSE

#endif
