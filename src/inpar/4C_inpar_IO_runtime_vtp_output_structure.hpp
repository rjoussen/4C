// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_INPAR_IO_RUNTIME_VTP_OUTPUT_STRUCTURE_HPP
#define FOUR_C_INPAR_IO_RUNTIME_VTP_OUTPUT_STRUCTURE_HPP


/*----------------------------------------------------------------------*/
/* headers */
#include "4C_config.hpp"

#include "4C_io_input_spec.hpp"

#include <memory>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
namespace Inpar
{
  namespace IORuntimeVTPStructure
  {
    /// data format for written numeric data
    enum OutputDataFormat
    {
      binary,
      ascii,
      vague
    };

    /// set the valid parameters related to writing of VTP output at runtime
    Core::IO::InputSpec set_valid_parameters();

  }  // namespace IORuntimeVTPStructure
}  // namespace Inpar

FOUR_C_NAMESPACE_CLOSE

#endif
