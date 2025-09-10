// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_inpar_IO_runtime_output_fluid.hpp"

#include "4C_io_input_spec_builders.hpp"


FOUR_C_NAMESPACE_OPEN

namespace Inpar
{
  namespace IORuntimeOutput
  {
    namespace FLUID
    {
      /*----------------------------------------------------------------------*
       *----------------------------------------------------------------------*/
      Core::IO::InputSpec valid_parameters()
      {
        using namespace Core::IO::InputSpecBuilders;

        Core::IO::InputSpec spec = group("IO/RUNTIME VTK OUTPUT/FLUID",
            {
                // whether to write output for fluid
                parameter<bool>(
                    "OUTPUT_FLUID", {.description = "write fluid output", .default_value = false}),

                // whether to write velocity state
                parameter<bool>(
                    "VELOCITY", {.description = "write velocity output", .default_value = false}),

                // whether to write pressure state
                parameter<bool>(
                    "PRESSURE", {.description = "write pressure output", .default_value = false}),

                // whether to write acceleration state
                parameter<bool>("ACCELERATION",
                    {.description = "write acceleration output", .default_value = false}),

                // whether to write displacement state
                parameter<bool>("DISPLACEMENT",
                    {.description = "write displacement output", .default_value = false}),

                // whether to write displacement state
                parameter<bool>("GRIDVELOCITY",
                    {.description = "write grid velocity output", .default_value = false}),

                // whether to write element owner
                parameter<bool>("ELEMENT_OWNER",
                    {.description = "write element owner", .default_value = false}),

                // whether to write element GIDs
                parameter<bool>("ELEMENT_GID",
                    {.description = "write 4C internal element GIDs", .default_value = false}),

                // whether to write node GIDs
                parameter<bool>("NODE_GID",
                    {.description = "write 4C internal node GIDs", .default_value = false}),
            },
            {.required = false});
        return spec;
      }
    }  // namespace FLUID
  }  // namespace IORuntimeOutput
}  // namespace Inpar

FOUR_C_NAMESPACE_CLOSE
