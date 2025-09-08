// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_lubrication_input.hpp"

#include "4C_io_input_spec_builders.hpp"
#include "4C_utils_parameter_list.hpp"

FOUR_C_NAMESPACE_OPEN

Core::IO::InputSpec Lubrication::set_valid_parameters()
{
  using namespace Core::IO::InputSpecBuilders;
  Core::IO::InputSpec spec = group("LUBRICATION DYNAMIC",
      {

          parameter<double>(
              "MAXTIME", {.description = "Total simulation time", .default_value = 1000.0}),
          parameter<int>(
              "NUMSTEP", {.description = "Total number of time steps", .default_value = 20}),

          parameter<double>("TIMESTEP", {.description = "Time increment dt", .default_value = 0.1}),
          parameter<int>("RESULTSEVERY",
              {.description = "Increment for writing solution", .default_value = 1}),
          parameter<int>(
              "RESTARTEVERY", {.description = "Increment for writing restart", .default_value = 1}),


          deprecated_selection<Lubrication::CalcError>("CALCERROR",
              {
                  {"No", calcerror_no},
                  {"error_by_function", calcerror_byfunction},
              },
              {.description = "compute error compared to analytical solution",
                  .default_value = calcerror_no}),

          parameter<int>(
              "CALCERRORNO", {.description = "function number for lubrication error computation",
                                 .default_value = -1}),


          deprecated_selection<Lubrication::VelocityField>("VELOCITYFIELD",
              {
                  {"zero", velocity_zero},
                  {"function", velocity_function},
                  {"EHL", velocity_EHL},
              },
              {.description = "type of velocity field used for lubrication problems",
                  .default_value = velocity_zero}),

          parameter<int>(
              "VELFUNCNO", {.description = "function number for lubrication velocity field",
                               .default_value = -1}),

          deprecated_selection<Lubrication::HeightField>("HEIGHTFEILD",
              {
                  {"zero", height_zero},
                  {"function", height_function},
                  {"EHL", height_EHL},
              },
              {.description = "type of height field used for lubrication problems",
                  .default_value = height_zero}),

          parameter<int>("HFUNCNO",
              {.description = "function number for lubrication height field", .default_value = -1}),

          parameter<bool>(
              "OUTMEAN", {.description = "Output of mean values for scalars and density",
                             .default_value = false}),

          parameter<bool>(
              "OUTPUT_GMSH", {.description = "Do you want to write Gmsh postprocessing files?",
                                 .default_value = false}),

          parameter<bool>("MATLAB_STATE_OUTPUT",
              {.description = "Do you want to write the state solution to Matlab file?",
                  .default_value = false}),

          /// linear solver id used for lubrication problems
          parameter<int>("LINEAR_SOLVER",
              {.description = "number of linear solver used for the Lubrication problem",
                  .default_value = -1}),

          parameter<int>(
              "ITEMAX", {.description = "max. number of nonlin. iterations", .default_value = 10}),
          parameter<double>(
              "ABSTOLRES", {.description = "Absolute tolerance for deciding if residual "
                                           "of nonlinear problem is already zero",
                               .default_value = 1e-14}),
          parameter<double>("CONVTOL",
              {.description = "Tolerance for convergence check", .default_value = 1e-13}),

          // convergence criteria adaptivity
          parameter<bool>("ADAPTCONV", {.description = "Switch on adaptive control of linear "
                                                       "solver tolerance for nonlinear solution",
                                           .default_value = false}),
          parameter<double>("ADAPTCONV_BETTER",
              {.description =
                      "The linear solver shall be this much better than the current nonlinear "
                      "residual in the nonlinear convergence limit",
                  .default_value = 0.1}),

          deprecated_selection<ConvNorm>("NORM_PRE",
              {
                  {"Abs", convnorm_abs},
                  {"Rel", convnorm_rel},
                  {"Mix", convnorm_mix},
              },
              {.description = "type of norm for temperature convergence check",
                  .default_value = convnorm_abs}),

          deprecated_selection<ConvNorm>("NORM_RESF",
              {
                  {"Abs", convnorm_abs},
                  {"Rel", convnorm_rel},
                  {"Mix", convnorm_mix},
              },
              {.description = "type of norm for residual convergence check",
                  .default_value = convnorm_abs}),

          deprecated_selection<VectorNorm>("ITERNORM",
              {
                  {"L1", norm_l1},
                  {"L2", norm_l2},
                  {"Rms", norm_rms},
                  {"Inf", norm_inf},
              },
              {.description = "type of norm to be applied to residuals", .default_value = norm_l2}),

          /// Iterationparameters
          parameter<double>(
              "TOLPRE", {.description = "tolerance in the temperature norm of the Newton iteration",
                            .default_value = 1.0E-06}),

          parameter<double>(
              "TOLRES", {.description = "tolerance in the residual norm for the Newton iteration",
                            .default_value = 1.0E-06}),

          parameter<double>("PENALTY_CAVITATION",
              {.description = "penalty parameter for regularized cavitation", .default_value = 0.}),

          parameter<double>("GAP_OFFSET",
              {.description = "Additional offset to the fluid gap", .default_value = 0.}),

          parameter<double>("ROUGHNESS_STD_DEVIATION",
              {.description = "standard deviation of surface roughness", .default_value = 0.}),

          /// use modified reynolds equ.
          parameter<bool>("MODIFIED_REYNOLDS_EQU",
              {.description =
                      "the lubrication problem will use the modified reynolds equ. in order to "
                      "consider surface roughness",
                  .default_value = false}),

          /// Flag for considering the Squeeze term in Reynolds Equation
          parameter<bool>("ADD_SQUEEZE_TERM",
              {.description = "the squeeze term will also be considered in the Reynolds Equation",
                  .default_value = false}),

          /// Flag for considering the pure Reynolds Equation
          parameter<bool>("PURE_LUB",
              {.description = "the problem is pure lubrication", .default_value = false})},
      {.required = false});
  return spec;
}

FOUR_C_NAMESPACE_CLOSE