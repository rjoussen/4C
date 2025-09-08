// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_INPAR_IO_RUNTIME_OUTPUT_FLUID_HPP
#define FOUR_C_INPAR_IO_RUNTIME_OUTPUT_FLUID_HPP


/*----------------------------------------------------------------------*/
/* headers */
#include "4C_config.hpp"

#include "4C_io_input_spec.hpp"

#include <memory>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
namespace Inpar
{
  namespace IORuntimeOutput
  {
    namespace FLUID
    {
      /// set the valid parameters related to writing of output at runtime
      Core::IO::InputSpec set_valid_parameters();

    }  // namespace FLUID
  }  // namespace IORuntimeOutput
}  // namespace Inpar

FOUR_C_NAMESPACE_CLOSE

#endif
