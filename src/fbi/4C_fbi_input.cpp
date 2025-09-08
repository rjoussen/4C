// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_fbi_input.hpp"

#include "4C_fem_condition_definition.hpp"
#include "4C_geometry_pair_input.hpp"
#include "4C_io_input_spec_builders.hpp"

FOUR_C_NAMESPACE_OPEN

std::vector<Core::IO::InputSpec> FBI::set_valid_parameters()
{
  using namespace Core::IO::InputSpecBuilders;

  /*----------------------------------------------------------------------*/
  /* parameters for beam to fluid meshtying */

  std::vector<Core::IO::InputSpec> specs;
  specs.push_back(group("FLUID BEAM INTERACTION",
      {
          deprecated_selection<BeamToFluidCoupling>("COUPLING",
              {
                  {"two-way", BeamToFluidCoupling::twoway},
                  {"fluid", BeamToFluidCoupling::fluid},
                  {"solid", BeamToFluidCoupling::solid},
              },
              {.description = "Type of FBI coupling",
                  .default_value = BeamToFluidCoupling::twoway}),

          parameter<int>("STARTSTEP",
              {.description = "Time Step at which to begin the fluid beam coupling. Usually "
                              "this will be the first step.",
                  .default_value = 0}),

          parameter<BeamToFluidPreSortStrategy>(
              "PRESORT_STRATEGY", {.description = "Presort strategy for the beam elements",
                                      .default_value = BeamToFluidPreSortStrategy::bruteforce}),

      },
      {.required = false}));

  /*----------------------------------------------------------------------*/
  std::vector<Core::IO::InputSpec> beam_to_fluid_meshtying = {

      parameter<BeamToFluidDiscretization>(
          "MESHTYING_DISCRETIZATION", {.description = "Type of employed meshtying discretization",
                                          .default_value = BeamToFluidDiscretization::none}),

      parameter<BeamToFluidConstraintEnforcement>(
          "CONSTRAINT_STRATEGY", {.description = "Type of employed constraint enforcement strategy",
                                     .default_value = BeamToFluidConstraintEnforcement::none}),

      parameter<double>("PENALTY_PARAMETER",
          {.description = "Penalty parameter for beam-to-Fluid volume meshtying",
              .default_value = 0.0}),

      parameter<double>("SEARCH_RADIUS",
          {.description =
                  "Absolute Search radius for beam-to-fluid volume meshtying. Choose carefully "
                  "to not blow up memory demand but to still find all interaction pairs!",
              .default_value = 1000.0}),


      parameter<FBI::BeamToFluidMeshtingMortarShapefunctions>("MORTAR_SHAPE_FUNCTION",
          {.description = "Shape function for the mortar Lagrange-multipliers",
              .default_value = FBI::BeamToFluidMeshtingMortarShapefunctions::none}),
  };
  // Add the geometry pair input parameters.
  GeometryPair::set_valid_parameters_line_to3_d(beam_to_fluid_meshtying);
  specs.push_back(group("FLUID BEAM INTERACTION/BEAM TO FLUID MESHTYING", beam_to_fluid_meshtying,
      {.required = false}));



  /*----------------------------------------------------------------------*/

  // Create subsection for runtime output.
  specs.push_back(group("FLUID BEAM INTERACTION/BEAM TO FLUID MESHTYING/RUNTIME VTK OUTPUT",
      {
          // Whether to write visualization output at all for beam to fluid meshtying.
          parameter<bool>(
              "WRITE_OUTPUT", {.description = "Enable / disable beam-to-fluid mesh tying output.",
                                  .default_value = false}),

          parameter<bool>("NODAL_FORCES",
              {.description = "Enable / disable output of the resulting nodal forces due "
                              "to beam to Fluid interaction.",
                  .default_value = false}),

          parameter<bool>(
              "SEGMENTATION", {.description = "Enable / disable output of segmentation points.",
                                  .default_value = false}),

          parameter<bool>("INTEGRATION_POINTS",
              {.description =
                      "Enable / disable output of used integration points. If the meshtying method "
                      "has 'forces' at the integration point, they will also be output.",
                  .default_value = false}),

          parameter<bool>("CONSTRAINT_VIOLATION",
              {.description = "Enable / disable output of the constraint violation "
                              "into a output_name.penalty csv file.",
                  .default_value = false}),

          parameter<bool>("MORTAR_LAMBDA_DISCRET",
              {.description =
                      "Enable / disable output of the discrete Lagrange multipliers at the node of "
                      "the Lagrange multiplier shape functions.",
                  .default_value = false}),

          parameter<bool>("MORTAR_LAMBDA_CONTINUOUS",
              {.description = "Enable / disable output of the continuous "
                              "Lagrange multipliers function along the beam.",
                  .default_value = false}),


          parameter<int>("MORTAR_LAMBDA_CONTINUOUS_SEGMENTS",
              {.description = "Number of segments for continuous mortar output",
                  .default_value = 5}),

      },
      {.required = false}));
  return specs;
}

void FBI::set_valid_conditions(std::vector<Core::Conditions::ConditionDefinition>& condlist)
{ /*-------------------------------------------------------------------*/ }

FOUR_C_NAMESPACE_CLOSE
