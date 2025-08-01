// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_UTILS_COMPILE_TIME_STRING_HPP
#define FOUR_C_UTILS_COMPILE_TIME_STRING_HPP

#include "4C_config.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <string_view>

FOUR_C_NAMESPACE_OPEN

namespace Core::Utils
{
  /**
   * A compile-time string class that can be used to represent strings as template parameters.
   */
  template <size_t size_with_null_terminator>
  struct CompileTimeString
  {
    static constexpr size_t size = size_with_null_terminator - 1;
    std::array<char, size> value{};

    constexpr CompileTimeString(const char (&str)[size_with_null_terminator])
    {
      std::copy_n(str, size, value.begin());
    }

    constexpr operator std::string_view() const { return {value.data(), size}; }

    template <std::size_t other_size>
    constexpr bool operator==(const CompileTimeString<other_size>& other) const
    {
      if constexpr (size_with_null_terminator != other_size)
        return false;
      else
      {
        return value == other.value;
      }
    }
  };

  // Deduction guide.
  template <size_t n>
  CompileTimeString(const char (&)[n]) -> CompileTimeString<n>;
}  // namespace Core::Utils

FOUR_C_NAMESPACE_CLOSE

#endif
