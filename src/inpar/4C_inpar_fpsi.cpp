// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_inpar_fpsi.hpp"

#include "4C_fem_condition_definition.hpp"
#include "4C_io_input_spec_builders.hpp"
FOUR_C_NAMESPACE_OPEN



Core::IO::InputSpec Inpar::FPSI::set_valid_parameters()
{
  using namespace Core::IO::InputSpecBuilders;
  Core::IO::InputSpec spec = group("FPSI DYNAMIC",
      {

          parameter<FpsiCouplingType>(
              "COUPALGO", {.description = "Iteration Scheme over the fields",
                              .default_value = fpsi_monolithic_plain}),

          parameter<bool>("SHAPEDERIVATIVES",
              {.description =
                      "Include linearization with respect to mesh movement in Navier Stokes "
                      "equation.\nSupported in monolithic FPSI for now.",
                  .default_value = false}),

          parameter<bool>("USESHAPEDERIVATIVES",
              {.description =
                      "Add linearization with respect to mesh movement in Navier Stokes equation "
                      "to stiffness matrix.\nSupported in monolithic FPSI for now.",
                  .default_value = false}),

          parameter<Inpar::FPSI::PartitionedCouplingMethod>(
              "PARTITIONED", {.description = "Coupling strategies for partitioned FPSI solvers.",
                                 .default_value = RobinNeumann}),

          parameter<bool>("SECONDORDER",
              {.description = "Second order coupling at the interface.", .default_value = false}),

          // Iterationparameters
          parameter<std::string>("RESTOL",
              {.description =
                      "Tolerances for single fields in the residual norm for the Newton iteration. "
                      "\nFor NORM_RESF != Abs_sys_split only the first value is used for all "
                      "fields. \nOrder of fields: porofluidvelocity, porofluidpressure, "
                      "porostructure, fluidvelocity, fluidpressure, ale",
                  .default_value = "1e-8 1e-8 1e-8 1e-8 1e-8 1e-8"}),

          parameter<std::string>("INCTOL",
              {.description = "Tolerance in the increment norm for the Newton iteration. \nFor "
                              "NORM_INC != \\*_split only the first value is used for all "
                              "fields. \nOrder of fields: porofluidvelocity, porofluidpressure, "
                              "porostructure, fluidvelocity, fluidpressure, ale",
                  .default_value = "1e-8 1e-8 1e-8 1e-8 1e-8 1e-8"}),

          deprecated_selection<Inpar::FPSI::ConvergenceNorm>("NORM_INC",
              {
                  {"Abs", absoluteconvergencenorm},
                  {"Abs_sys_split", absoluteconvergencenorm_sys_split},
                  {"Rel_sys", relativconvergencenorm_sys},
              },
              {.description = "Type of norm for primary variables convergence check.  \nAbs: "
                              "absolute values, "
                              "Abs_sys_split: absolute values with correction of systemsize for "
                              "every field "
                              "separate, Rel_sys: relative values with correction of systemsize.",
                  .default_value = absoluteconvergencenorm}),

          deprecated_selection<Inpar::FPSI::ConvergenceNorm>("NORM_RESF",
              {
                  {"Abs", absoluteconvergencenorm},
                  {"Abs_sys_split", absoluteconvergencenorm_sys_split},
                  {"Rel_sys", relativconvergencenorm_sys},
              },
              {.description = "Type of norm for primary variables convergence check. \nAbs: "
                              "absolute values, "
                              "Abs_sys_split: absolute values with correction of systemsize for "
                              "every field "
                              "separate, Rel_sys: relative values with correction of systemsize.",
                  .default_value = absoluteconvergencenorm}),

          deprecated_selection<Inpar::FPSI::BinaryOp>("NORMCOMBI_RESFINC",
              {
                  {"And", bop_and},
                  {"Or", bop_or},
              },
              {.description =
                      "binary operator to combine primary variables and residual force values",
                  .default_value = bop_and}),


          parameter<bool>(
              "LineSearch", {.description = "adapt increment in case of non-monotonic "
                                            "residual convergence or residual oscillations",
                                .default_value = false}),

          parameter<bool>(
              "FDCheck", {.description = "perform FPSIFDCheck() finite difference check",
                             .default_value = false}),

          // number of linear solver used for poroelasticity
          parameter<int>(
              "LINEAR_SOLVER", {.description = "number of linear solver used for FPSI problems",
                                   .default_value = -1}),
          parameter<int>("ITEMIN",
              {.description = "minimal number of iterations over fields", .default_value = 1}),
          parameter<int>(
              "NUMSTEP", {.description = "Total number of Timesteps", .default_value = 200}),
          parameter<int>("ITEMAX",
              {.description = "Maximum number of iterations over fields", .default_value = 100}),
          parameter<int>("RESULTSEVERY",
              {.description = "Increment for writing solution", .default_value = 1}),
          parameter<int>(
              "RESTARTEVERY", {.description = "Increment for writing restart", .default_value = 1}),
          parameter<int>("FDCheck_row",
              {.description = "print row value during fd_check", .default_value = 0}),
          parameter<int>("FDCheck_column",
              {.description = "print column value during fd_check", .default_value = 0}),


          parameter<double>("TIMESTEP", {.description = "Time increment dt", .default_value = 0.1}),
          parameter<double>(
              "MAXTIME", {.description = "Total simulation time", .default_value = 1000.0}),
          parameter<double>("CONVTOL",
              {.description = "Tolerance for iteration over fields", .default_value = 1e-6}),
          parameter<double>("ALPHABJ",
              {.description = "Beavers-Joseph-Coefficient for Slip-Boundary-Condition at "
                              "Fluid-Porous-Interface (0.1-4)",
                  .default_value = 1.0})},
      {.required = false});
  return spec;
}



void Inpar::FPSI::set_valid_conditions(std::vector<Core::Conditions::ConditionDefinition>& condlist)
{
  using namespace Core::IO::InputSpecBuilders;

  /*--------------------------------------------------------------------*/
  // FPSI

  Core::Conditions::ConditionDefinition linefpsi("DESIGN FPSI COUPLING LINE CONDITIONS",
      "fpsi_coupling", "FPSI Coupling", Core::Conditions::fpsi_coupling, true,
      Core::Conditions::geometry_type_line);
  Core::Conditions::ConditionDefinition surffpsi("DESIGN FPSI COUPLING SURF CONDITIONS",
      "fpsi_coupling", "FPSI Coupling", Core::Conditions::fpsi_coupling, true,
      Core::Conditions::geometry_type_surface);

  linefpsi.add_component(parameter<int>("coupling_id"));
  surffpsi.add_component(parameter<int>("coupling_id"));

  condlist.push_back(linefpsi);
  condlist.push_back(surffpsi);


  /*--------------------------------------------------------------------*/
  // condition for evaluation of boundary terms in fpsi problems
  // necessary where neumann term needs to be integrated in interface
  // elements which share a node with the fpsi interface. Tangential
  // Beaver-Joseph-Condition must not be overwritten by prescribed value!

  Core::Conditions::ConditionDefinition neumannintegration_surf(
      "DESIGN SURFACE NEUMANN INTEGRATION", "NeumannIntegration", "Neumann Integration",
      Core::Conditions::NeumannIntegration, true, Core::Conditions::geometry_type_surface);

  condlist.push_back(neumannintegration_surf);

  /*--------------------------------------------------------------------*/
  // condition for evaluation of boundary terms in fpsi problems

  Core::Conditions::ConditionDefinition neumannintegration_line("DESIGN LINE NEUMANN INTEGRATION",
      "NeumannIntegration", "Neumann Integration", Core::Conditions::NeumannIntegration, true,
      Core::Conditions::geometry_type_line);

  condlist.push_back(neumannintegration_line);
}

FOUR_C_NAMESPACE_CLOSE