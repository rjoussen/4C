// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_UTILS_SYMBOLIC_EXPRESSION_FWD_HPP
#define FOUR_C_UTILS_SYMBOLIC_EXPRESSION_FWD_HPP

#include "4C_config.hpp"

#include "4C_utils_compile_time_string.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Core::Utils
{
  /**
   * Get index of CompileTimeString @p s in the list of CompileTimeStrings @p strings. Returns -1 if
   * @p s is not found in @p strings.
   */
  template <CompileTimeString s, CompileTimeString... strings>
  consteval int index_of()
  {
    int index = 0;

    bool found = ((s == strings ? true : (++index, false)) || ...);
    if (!found) return -1;
    return index;
  }


  /**
   * @brief A variable in a symbolic expression.
   *
   * Represents a compile-time variable @p name and a runtime @p value.
   */
  template <CompileTimeString name, typename Number>
  struct VarWrapper
  {
    Number value;
  };

  /**
   * Create a variable for passing it to SymbolicExpression.
   */
  template <CompileTimeString name, typename Number>
  constexpr VarWrapper<name, Number> var(Number value)
  {
    return VarWrapper<name, Number>{value};
  }


  template <typename Number, CompileTimeString... variables>
  class SymbolicExpression;
}  // namespace Core::Utils

FOUR_C_NAMESPACE_CLOSE

#endif
