// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_INPAR_PROBLEMTYPE_HPP
#define FOUR_C_INPAR_PROBLEMTYPE_HPP

#include "4C_config.hpp"

#include "4C_io_input_spec.hpp"
#include "4C_legacy_enum_definitions_problem_type.hpp"

#include <map>
#include <string>

FOUR_C_NAMESPACE_OPEN

namespace Inpar
{
  namespace PROBLEMTYPE
  {
    /*! \brief Define valid parameters
     *
     * @param[in/out] list Parameter list to be filled with valid parameters and their defaults
     */
    Core::IO::InputSpec set_valid_parameters();

    /// create map of problem name and problem type enum
    std::map<std::string, Core::ProblemType> string_to_problem_type_map();

    /// return problem type enum for a given problem name
    Core::ProblemType string_to_problem_type(std::string name);


  }  // namespace PROBLEMTYPE
}  // namespace Inpar

/*----------------------------------------------------------------------*/
FOUR_C_NAMESPACE_CLOSE

#endif
