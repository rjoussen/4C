// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_inpar_IO_runtime_vtk_output_structure.hpp"

#include "4C_inpar_structure.hpp"
#include "4C_io_input_spec_builders.hpp"


FOUR_C_NAMESPACE_OPEN

namespace Inpar
{
  namespace IORuntimeOutput
  {
    namespace Solid
    {
      /*----------------------------------------------------------------------*
       *----------------------------------------------------------------------*/
      Core::IO::InputSpec set_valid_parameters()
      {
        using namespace Core::IO::InputSpecBuilders;

        // related sublist
        Core::IO::InputSpec spec = group("IO/RUNTIME VTK OUTPUT/STRUCTURE",
            {

                // whether to write output for structure
                parameter<bool>("OUTPUT_STRUCTURE",
                    {.description = "write structure output", .default_value = false}),

                // whether to write displacement state
                parameter<bool>("DISPLACEMENT",
                    {.description = "write displacement output", .default_value = false}),

                // whether to write velocity state
                parameter<bool>(
                    "VELOCITY", {.description = "write velocity output", .default_value = false}),

                parameter<bool>("ACCELERATION",
                    {.description = "write acceleration output", .default_value = false}),

                // whether to write element owner
                parameter<bool>("ELEMENT_OWNER",
                    {.description = "write element owner", .default_value = false}),

                // whether to write element GIDs
                parameter<bool>("ELEMENT_GID",
                    {.description = "write 4C internal element GIDs", .default_value = false}),

                // write element ghosting information
                parameter<bool>(
                    "ELEMENT_GHOSTING", {.description = "write which processors ghost the elements",
                                            .default_value = false}),

                // whether to write node GIDs
                parameter<bool>("NODE_GID",
                    {.description = "write 4C internal node GIDs", .default_value = false}),

                // write element material IDs
                parameter<bool>(
                    "ELEMENT_MAT_ID", {.description = "Output of the material id of each element",
                                          .default_value = false}),

                // whether to write stress and / or strain data
                parameter<bool>("STRESS_STRAIN",
                    {.description =
                            "Write element stress and / or strain  data. The type of stress / "
                            "strain has to be selected in the --IO input section",
                        .default_value = false}),

                // mode to write gauss point data
                parameter<Inpar::Solid::GaussPointDataOutputType>("GAUSS_POINT_DATA_OUTPUT_TYPE",
                    {.description = "Where to write gauss point data. (none, projected to nodes, "
                                    "projected to element center, raw at gauss points)",
                        .default_value = Inpar::Solid::GaussPointDataOutputType::none}),


                deprecated_selection<Inpar::Solid::OptQuantityType>("OPTIONAL_QUANTITY",
                    {
                        {"No", Inpar::Solid::optquantity_none},
                        {"membranethickness", Inpar::Solid::optquantity_membranethickness},
                    },
                    {.description = "Output of an optional quantity",
                        .default_value = Inpar::Solid::optquantity_none})},
            {.required = false});
        return spec;
      }


    }  // namespace Solid
  }  // namespace IORuntimeOutput
}  // namespace Inpar

FOUR_C_NAMESPACE_CLOSE