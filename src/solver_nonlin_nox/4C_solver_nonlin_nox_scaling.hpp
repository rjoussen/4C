// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_SOLVER_NONLIN_NOX_SCALING_HPP
#define FOUR_C_SOLVER_NONLIN_NOX_SCALING_HPP

#include "4C_config.hpp"

FOUR_C_NAMESPACE_OPEN

namespace NOX::Nln
{
  // Forward declaration
  class LinearProblem;

  class Scaling
  {
   public:
    //! Virtual destructor
    virtual ~Scaling() = default;

    /**
     * @brief Compute the scaling for the linear system.
     *
     * This method is mainly to provide compatibility with the current pattern how scaling is
     * applied: if flag manualScaling is set to false, then this method is called.
     *
     * @param problem The linear problem to be scaled.
     */
    virtual void compute_scaling(const NOX::Nln::LinearProblem& problem) = 0;

    /**
     * @brief Scales the linear system.
     *
     * @param problem The linear problem to be scaled.
     */
    virtual void scale_linear_system(NOX::Nln::LinearProblem& problem) = 0;

    /**
     * @brief Remove the scaling from the linear system.
     *
     * @param problem The linear problem to be unscaled.
     */
    virtual void unscale_linear_system(NOX::Nln::LinearProblem& problem) = 0;
  };
}  // namespace NOX::Nln



FOUR_C_NAMESPACE_CLOSE

#endif
