// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_inpar_poroscatra.hpp"

#include "4C_inpar_scatra.hpp"
#include "4C_io_input_spec_builders.hpp"
#include "4C_poroelast_input.hpp"
FOUR_C_NAMESPACE_OPEN



void Inpar::PoroScaTra::set_valid_parameters(std::map<std::string, Core::IO::InputSpec>& list)
{
  using namespace Core::IO::InputSpecBuilders;

  list["POROSCATRA CONTROL"] = group("POROSCATRA CONTROL",
      {

          // Output type
          parameter<int>(
              "RESTARTEVERY", {.description = "write restart possibility every RESTARTEVERY steps",
                                  .default_value = 1}),
          // Time loop control
          parameter<int>(
              "NUMSTEP", {.description = "maximum number of Timesteps", .default_value = 200}),
          parameter<double>(
              "MAXTIME", {.description = "total simulation time", .default_value = 1000.0}),

          parameter<double>(
              "TIMESTEP", {.description = "time step size dt", .default_value = 0.05}),
          parameter<int>("RESULTSEVERY",
              {.description = "increment for writing solution", .default_value = 1}),
          parameter<int>("ITEMAX",
              {.description = "maximum number of iterations over fields", .default_value = 10}),
          parameter<int>("ITEMIN",
              {.description = "minimal number of iterations over fields", .default_value = 1}),

          // Iterationparameters
          parameter<double>("TOLRES_GLOBAL",
              {.description = "tolerance in the residual norm for the Newton iteration",
                  .default_value = 1e-8}),
          parameter<double>("TOLINC_GLOBAL",
              {.description = "tolerance in the increment norm for the Newton iteration",
                  .default_value = 1e-8}),
          parameter<double>("TOLRES_DISP",
              {.description = "tolerance in the residual norm for the Newton iteration",
                  .default_value = 1e-8}),
          parameter<double>("TOLINC_DISP",
              {.description = "tolerance in the increment norm for the Newton iteration",
                  .default_value = 1e-8}),
          parameter<double>("TOLRES_VEL",
              {.description = "tolerance in the residual norm for the Newton iteration",
                  .default_value = 1e-8}),
          parameter<double>("TOLINC_VEL",
              {.description = "tolerance in the increment norm for the Newton iteration",
                  .default_value = 1e-8}),
          parameter<double>("TOLRES_PRES",
              {.description = "tolerance in the residual norm for the Newton iteration",
                  .default_value = 1e-8}),
          parameter<double>("TOLINC_PRES",
              {.description = "tolerance in the increment norm for the Newton iteration",
                  .default_value = 1e-8}),
          parameter<double>("TOLRES_SCALAR",
              {.description = "tolerance in the residual norm for the Newton iteration",
                  .default_value = 1e-8}),
          parameter<double>("TOLINC_SCALAR",
              {.description = "tolerance in the increment norm for the Newton iteration",
                  .default_value = 1e-8}),

          deprecated_selection<PoroElast::ConvNorm>("NORM_INC",
              {
                  {"AbsGlobal", PoroElast::convnorm_abs_global},
                  {"AbsSingleFields", PoroElast::convnorm_abs_singlefields},
              },
              {.description = "type of norm for primary variables convergence check",
                  .default_value = PoroElast::convnorm_abs_singlefields}),

          deprecated_selection<PoroElast::ConvNorm>("NORM_RESF",
              {
                  {"AbsGlobal", PoroElast::convnorm_abs_global},
                  {"AbsSingleFields", PoroElast::convnorm_abs_singlefields},
              },
              {.description = "type of norm for residual convergence check",
                  .default_value = PoroElast::convnorm_abs_singlefields}),


          deprecated_selection<PoroElast::BinaryOp>("NORMCOMBI_RESFINC",
              {
                  {"And", PoroElast::bop_and},
                  {"Or", PoroElast::bop_or},
              },
              {.description =
                      "binary operator to combine primary variables and residual force values",
                  .default_value = PoroElast::bop_and}),


          deprecated_selection<PoroElast::VectorNorm>("VECTORNORM_RESF",
              {
                  {"L1", PoroElast::norm_l1},
                  {"L1_Scaled", PoroElast::norm_l1_scaled},
                  {"L2", PoroElast::norm_l2},
                  {"Rms", PoroElast::norm_rms},
                  {"Inf", PoroElast::norm_inf},
              },
              {.description = "type of norm to be applied to residuals",
                  .default_value = PoroElast::norm_l2}),


          deprecated_selection<PoroElast::VectorNorm>("VECTORNORM_INC",
              {
                  {"L1", PoroElast::norm_l1},
                  {"L1_Scaled", PoroElast::norm_l1_scaled},
                  {"L2", PoroElast::norm_l2},
                  {"Rms", PoroElast::norm_rms},
                  {"Inf", PoroElast::norm_inf},
              },
              {.description = "type of norm to be applied to residuals",
                  .default_value = PoroElast::norm_l2}),

          // number of linear solver used for poroelasticity
          parameter<int>("LINEAR_SOLVER",
              {.description = "number of linear solver used for monolithic poroscatra problems",
                  .default_value = -1}),

          // Coupling strategy for poroscatra solvers
          deprecated_selection<SolutionSchemeOverFields>("COUPALGO",
              {
                  {"monolithic", Monolithic},
                  {"scatra_to_solid", Part_ScatraToPoro},
                  {"solid_to_scatra", Part_PoroToScatra},
                  {"two_way", Part_TwoWay},
              },
              {.description = "Coupling strategies for poroscatra solvers",
                  .default_value = Part_PoroToScatra}),


          parameter<bool>(
              "MATCHINGGRID", {.description = "is matching grid", .default_value = true})},
      {.required = false});
}

FOUR_C_NAMESPACE_CLOSE