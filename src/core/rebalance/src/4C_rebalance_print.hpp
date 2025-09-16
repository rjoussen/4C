// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_REBALANCE_PRINT_HPP
#define FOUR_C_REBALANCE_PRINT_HPP

#include "4C_config.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Core::FE
{
  class Discretization;
}  // namespace Core::FE

namespace Core::Rebalance
{
  /*! \brief Prints details on the distribution/partitioning of the distribution
   */
  void print_parallel_distribution(const Core::FE::Discretization& dis);

}  // namespace Core::Rebalance

FOUR_C_NAMESPACE_CLOSE

#endif
