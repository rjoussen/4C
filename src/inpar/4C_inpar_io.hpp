// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_INPAR_IO_HPP
#define FOUR_C_INPAR_IO_HPP

#include "4C_config.hpp"

#include "4C_io_input_spec.hpp"

#include <vector>


FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
namespace Inpar
{
  namespace IO
  {
    /*! \brief Define valid parameter for global IO control
     *
     * @param[in/out] list Parameter list to be filled with valid parameters and their defaults
     */
    std::vector<Core::IO::InputSpec> set_valid_parameters();

  }  // namespace IO
}  // namespace Inpar

/*----------------------------------------------------------------------*/
FOUR_C_NAMESPACE_CLOSE

#endif
