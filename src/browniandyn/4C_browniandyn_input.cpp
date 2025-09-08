// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_browniandyn_input.hpp"

#include "4C_io_input_spec_builders.hpp"

FOUR_C_NAMESPACE_OPEN


Core::IO::InputSpec BrownianDynamics::set_valid_parameters()
{
  using namespace Core::IO::InputSpecBuilders;
  Core::IO::InputSpec spec = group("BROWNIAN DYNAMICS",
      {

          parameter<bool>("BROWNDYNPROB",
              {.description = "switch Brownian dynamics on/off", .default_value = false}),

          // Reading double parameter for viscosity of background fluid
          parameter<double>("VISCOSITY", {.description = "viscosity", .default_value = 0.0}),

          // Reading double parameter for thermal energy in background fluid (temperature *
          // Boltzmann constant)
          parameter<double>("KT", {.description = "thermal energy", .default_value = 0.0}),

          // cutoff for random forces, which determines the maximal value
          parameter<double>("MAXRANDFORCE",
              {.description = "Any random force beyond MAXRANDFORCE*(standard dev.) will "
                              "be omitted and redrawn. -1.0 means no bounds.'",
                  .default_value = -1.0}),

          // time interval in which random numbers are constant
          parameter<double>("TIMESTEP",
              {.description = "Within this time interval the random numbers remain constant. -1.0 ",
                  .default_value = -1.0}),

          // the way how damping coefficient values for beams are specified
          deprecated_selection<BeamDampingCoefficientSpecificationType>(
              "BEAMS_DAMPING_COEFF_SPECIFIED_VIA",
              {
                  {"cylinder_geometry_approx", BrownianDynamics::cylinder_geometry_approx},
                  {"Cylinder_geometry_approx", BrownianDynamics::cylinder_geometry_approx},
                  {"input_file", BrownianDynamics::input_file},
                  {"Input_file", BrownianDynamics::input_file},
              },
              {.description = "In which way are damping coefficient values for beams specified?",
                  .default_value = BrownianDynamics::cylinder_geometry_approx}),

          // values for damping coefficients of beams if they are specified via input file
          // (per unit length, NOT yet multiplied by fluid viscosity)
          parameter<std::string>("BEAMS_DAMPING_COEFF_PER_UNITLENGTH",
              {.description =
                      "values for beam damping coefficients (per unit length and NOT yet "
                      "multiplied by fluid viscosity): translational perpendicular/parallel to "
                      "beam axis, rotational around axis",
                  .default_value = "0.0 0.0 0.0"})},
      {.required = false});
  return spec;
}

FOUR_C_NAMESPACE_CLOSE