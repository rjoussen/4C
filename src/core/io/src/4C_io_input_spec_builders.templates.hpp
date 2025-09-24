// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_IO_INPUT_SPEC_BUILDERS_TEMPLATES_HPP
#define FOUR_C_IO_INPUT_SPEC_BUILDERS_TEMPLATES_HPP

#include "4C_config.hpp"

#include "4C_io_input_spec_builders.hpp"
#include "4C_utils_symbolic_expression.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Core::IO
{
  namespace InputSpecBuilders
  {
    template <typename Number, Utils::CompileTimeString... variables>
    struct SymbolicExpressionData
    {
      /**
       * An optional description of the symbolic expression.
       */
      std::string description{};

      /**
       * An optional function to store the parsed symbolic expression. By default, the result is
       * stored in an InputParameterContainer. See the in_struct() function for more details on how
       * to store the symbolic expression in a struct.
       */
      StoreFunction<Utils::SymbolicExpression<Number, variables...>> store{nullptr};
    };
  }  // namespace InputSpecBuilders
}  // namespace Core::IO


// --- template definitions --- //

template <typename Number, Core::Utils::CompileTimeString... variables>
Core::IO::InputSpec Core::IO::InputSpecBuilders::symbolic_expression(
    std::string name, SymbolicExpressionData<Number, variables...> data)
{
  auto store =
      data.store ? data.store : in_container<Utils::SymbolicExpression<Number, variables...>>(name);
  const std::type_info& stores_to = store.stores_to();

  StoreFunction<std::string> convert_to_symbolic_expression{[store](Storage& out, std::string&& in)
      {
        try
        {
          Utils::SymbolicExpression<Number, variables...> expr(std::move(in));
          return store(out, std::move(expr));
        }
        catch (const Core::Exception& e)
        {
          // Make a string from all variables.
          std::stringstream message;
          message << "could not be parsed as symbolic expression with variables: ";
          ((message << std::quoted(std::string_view(variables)) << " "), ...);

          return StoreStatus::fail(message.str());
        }
      },
      stores_to};

  return parameter<std::string>(name, {
                                          .description = data.description,
                                          .store = convert_to_symbolic_expression,
                                      });
}



FOUR_C_NAMESPACE_CLOSE

#endif
