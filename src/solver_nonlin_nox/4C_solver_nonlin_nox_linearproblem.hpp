// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_SOLVER_NONLIN_NOX_LINEARPROBLEM_HPP
#define FOUR_C_SOLVER_NONLIN_NOX_LINEARPROBLEM_HPP

#include "4C_config.hpp"

#include "4C_linalg_sparseoperator.hpp"
#include "4C_linalg_vector.hpp"

FOUR_C_NAMESPACE_OPEN

namespace NOX::Nln
{
  class LinearProblem
  {
   public:
    std::shared_ptr<Core::LinAlg::SparseOperator> jac;
    std::shared_ptr<Core::LinAlg::Vector<double>> lhs;
    std::shared_ptr<Core::LinAlg::Vector<double>> rhs;
  };
}  // namespace NOX::Nln

FOUR_C_NAMESPACE_CLOSE

#endif
