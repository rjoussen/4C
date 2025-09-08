// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_poroelast_input.hpp"

#include "4C_inpar_fluid.hpp"
#include "4C_io_input_spec_builders.hpp"
#include "4C_linalg_equilibrate.hpp"
FOUR_C_NAMESPACE_OPEN


Core::IO::InputSpec PoroElast::set_valid_parameters()
{
  using namespace Core::IO::InputSpecBuilders;
  Core::IO::InputSpec spec = group("POROELASTICITY DYNAMIC",
      {

          // Coupling strategy for (monolithic) porous media solvers
          parameter<PoroElast::SolutionSchemeOverFields>(
              "COUPALGO", {.description = "Coupling strategies for poroelasticity solvers",
                              .default_value = PoroElast::SolutionSchemeOverFields::Monolithic}),

          // physical type of poro fluid flow (incompressible, varying density, loma, Boussinesq
          // approximation)
          deprecated_selection<Inpar::FLUID::PhysicalType>("PHYSICAL_TYPE",
              {
                  {"Poro", Inpar::FLUID::poro},
                  {"Poro_P1", Inpar::FLUID::poro_p1},
              },
              {.description = "Physical Type of Porofluid", .default_value = Inpar::FLUID::poro}),

          // physical type of poro fluid flow (incompressible, varying density, loma, Boussinesq
          // approximation)
          deprecated_selection<PoroElast::TransientEquationsOfPoroFluid>("TRANSIENT_TERMS",
              {
                  {"none", transient_none},
                  {"momentum", transient_momentum_only},
                  {"continuity", transient_continuity_only},
                  {"all", transient_all},
              },
              {.description = "which equation includes transient terms",
                  .default_value = transient_all}),

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
          parameter<int>("ITEMAX",
              {.description = "maximum number of iterations over fields", .default_value = 10}),
          parameter<int>("ITEMIN",
              {.description = "minimal number of iterations over fields", .default_value = 1}),
          parameter<int>("RESULTSEVERY",
              {.description = "increment for writing solution", .default_value = 1}),

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
          parameter<double>("TOLRES_PORO",
              {.description = "tolerance in the residual norm for the Newton iteration",
                  .default_value = 1e-8}),
          parameter<double>("TOLINC_PORO",
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
          parameter<double>("TOLRES_NCOUP",
              {.description = "tolerance in the residual norm for the Newton iteration",
                  .default_value = 1e-8}),

          deprecated_selection<PoroElast::ConvNorm>("NORM_INC",
              {
                  {"AbsGlobal", convnorm_abs_global},
                  {"AbsSingleFields", convnorm_abs_singlefields},
              },
              {.description = "type of norm for primary variables convergence check",
                  .default_value = convnorm_abs_singlefields}),

          deprecated_selection<PoroElast::ConvNorm>("NORM_RESF",
              {
                  {"AbsGlobal", convnorm_abs_global},
                  {"AbsSingleFields", convnorm_abs_singlefields},
              },
              {.description = "type of norm for residual convergence check",
                  .default_value = convnorm_abs_singlefields}),


          deprecated_selection<PoroElast::BinaryOp>("NORMCOMBI_RESFINC",
              {
                  {"And", bop_and},
                  {"Or", bop_or},
              },
              {.description =
                      "binary operator to combine primary variables and residual force values",
                  .default_value = bop_and}),


          deprecated_selection<PoroElast::VectorNorm>("VECTORNORM_RESF",
              {
                  {"L1", norm_l1},
                  {"L1_Scaled", norm_l1_scaled},
                  {"L2", norm_l2},
                  {"Rms", norm_rms},
                  {"Inf", norm_inf},
              },
              {.description = "type of norm to be applied to residuals", .default_value = norm_l2}),


          deprecated_selection<PoroElast::VectorNorm>("VECTORNORM_INC",
              {
                  {"L1", norm_l1},
                  {"L1_Scaled", norm_l1_scaled},
                  {"L2", norm_l2},
                  {"Rms", norm_rms},
                  {"Inf", norm_inf},
              },
              {.description = "type of norm to be applied to residuals", .default_value = norm_l2}),

          parameter<bool>("SECONDORDER",
              {.description = "Second order coupling at the interface.", .default_value = true}),

          parameter<bool>("CONTIPARTINT",
              {.description = "Partial integration of porosity gradient in continuity equation",
                  .default_value = false}),

          parameter<bool>("CONTACT_NO_PENETRATION",
              {.description = "No-Penetration Condition on active contact "
                              "surface in case of poro contact problem!",
                  .default_value = false}),


          parameter<bool>(
              "MATCHINGGRID", {.description = "is matching grid", .default_value = true}),

          parameter<bool>(
              "CONVECTIVE_TERM", {.description = "convective term ", .default_value = false}),

          // number of linear solver used for poroelasticity
          parameter<int>("LINEAR_SOLVER",
              {.description = "number of linear solver used for poroelasticity problems",
                  .default_value = -1}),

          // flag for equilibration of global system of equations
          parameter<Core::LinAlg::EquilibrationMethod>("EQUILIBRATION",
              {.description = "flag for equilibration of global system of equations",
                  .default_value = Core::LinAlg::EquilibrationMethod::none})},
      {.required = false});
  return spec;
}

FOUR_C_NAMESPACE_CLOSE