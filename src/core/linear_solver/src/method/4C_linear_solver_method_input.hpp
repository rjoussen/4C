// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_LINEAR_SOLVER_METHOD_INPUT_HPP
#define FOUR_C_LINEAR_SOLVER_METHOD_INPUT_HPP

#include "4C_config.hpp"

#include "4C_io_input_spec.hpp"

#include <map>
#include <vector>

FOUR_C_NAMESPACE_OPEN

namespace Core::LinearSolver
{
  //! the parameters for the linear solver
  std::vector<Core::IO::InputSpec> valid_parameters();

}  // namespace Core::LinearSolver

FOUR_C_NAMESPACE_CLOSE

#endif
