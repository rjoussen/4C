// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_config.hpp"

#include "4C_utils_exceptions.hpp"
#include "4C_utils_function.hpp"

#ifndef FOUR_C_FLUID_XFLUID_FUNCTIONS_COMBUST_HPP
#define FOUR_C_FLUID_XFLUID_FUNCTIONS_COMBUST_HPP

FOUR_C_NAMESPACE_OPEN

namespace Core::Utils
{
  class FunctionManager;
}

namespace Core::FE
{
  class Discretization;
}  // namespace Core::FE

namespace Discret
{
  // namespace Utils

  namespace Utils
  {
    /// add valid combustion-specific function lines
    void add_valid_combust_functions(Core::Utils::FunctionManager& function_manager);

    /// special implementation for a level set test function "Zalesak's disk"
    class ZalesaksDiskFunction : public Core::Utils::FunctionOfSpaceTime
    {
     public:
      double evaluate(std::span<const double> x, double t, std::size_t component) const override;

      [[nodiscard]] std::size_t number_components() const override
      {
        FOUR_C_THROW("Number of components not defined for ZalesaksDiskFunction.");
      };
    };

    /// special implementation two-phase flow test case
    class CollapsingWaterColumnFunction : public Core::Utils::FunctionOfSpaceTime
    {
     public:
      double evaluate(std::span<const double> x, double t, std::size_t component) const override;

      [[nodiscard]] std::size_t number_components() const override
      {
        FOUR_C_THROW("Number of components not defined for CollapsingWaterColumnFunction.");
      };
    };
  }  // namespace Utils
}  // namespace Discret

FOUR_C_NAMESPACE_CLOSE

#endif
