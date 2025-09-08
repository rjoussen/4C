// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_ehl_input.hpp"

#include "4C_fem_condition_definition.hpp"
#include "4C_io_input_spec_builders.hpp"
#include "4C_utils_parameter_list.hpp"

FOUR_C_NAMESPACE_OPEN



std::vector<Core::IO::InputSpec> EHL::set_valid_parameters()
{
  using namespace Core::IO::InputSpecBuilders;

  std::vector<Core::IO::InputSpec> specs;
  specs.push_back(group("ELASTO HYDRO DYNAMIC",
      {

          // Output type
          parameter<double>("RESTARTEVERYTIME",
              {.description = "write restart possibility every RESTARTEVERY steps",
                  .default_value = 0.0}),
          parameter<int>(
              "RESTARTEVERY", {.description = "write restart possibility every RESTARTEVERY steps",
                                  .default_value = 1}),
          // Time loop control
          parameter<int>(
              "NUMSTEP", {.description = "maximum number of Timesteps", .default_value = 200}),
          parameter<double>(
              "MAXTIME", {.description = "total simulation time", .default_value = 1000.0}),

          parameter<double>(
              "TIMESTEP", {.description = "time step size dt", .default_value = -1.0}),
          parameter<bool>("DIFFTIMESTEPSIZE",
              {.description = "use different step size for lubrication and solid",
                  .default_value = false}),
          parameter<double>("RESULTSEVERYTIME",
              {.description = "increment for writing solution", .default_value = 0.0}),
          parameter<int>("RESULTSEVERY",
              {.description = "increment for writing solution", .default_value = 1}),
          parameter<int>("ITEMAX",
              {.description = "maximum number of iterations over fields", .default_value = 10}),
          parameter<int>("ITEMIN",
              {.description = "minimal number of iterations over fields", .default_value = 1}),

          // Type of coupling strategy between the two fields
          deprecated_selection<FieldCoupling>("FIELDCOUPLING",
              {
                  {"none", coupling_none},
                  {"matching", coupling_matching},
              },
              {.description = "Type of coupling strategy between fields",
                  .default_value = coupling_none}),

          // Coupling strategy for EHL solvers
          parameter<SolutionSchemeOverFields>(
              "COUPALGO", {.description = "Coupling strategies for EHL solvers",
                              .default_value = ehl_Monolithic}),

          // set unprojectable nodes to zero pressure via Dirichlet condition
          parameter<bool>("UNPROJ_ZERO_DBC",
              {.description = "set unprojectable nodes to zero pressure via Dirichlet condition",
                  .default_value = false}),

          // use dry contact model
          parameter<bool>("DRY_CONTACT_MODEL",
              {.description = "set unprojectable nodes to zero pressure via Dirichlet condition",
                  .default_value = false})},
      {.required = false}));
  /* parameters for monolithic EHL */
  specs.push_back(group("ELASTO HYDRO DYNAMIC/MONOLITHIC",
      {

          // convergence tolerance of EHL residual
          parameter<double>("CONVTOL",
              {.description = "tolerance for convergence check of EHL", .default_value = 1e-6}),
          // Iterationparameters
          parameter<double>("TOLINC",
              {.description = "tolerance for convergence check of EHL-increment in monolithic EHL",
                  .default_value = 1.0e-6}),

          deprecated_selection<ConvNorm>("NORM_RESF",
              {
                  {"Abs", convnorm_abs},
                  {"Rel", convnorm_rel},
                  {"Mix", convnorm_mix},
              },
              {.description = "type of norm for residual convergence check",
                  .default_value = convnorm_abs}),

          deprecated_selection<ConvNorm>("NORM_INC",
              {
                  {"Abs", convnorm_abs},
                  {"Rel", convnorm_rel},
                  {"Mix", convnorm_mix},
              },
              {.description = "type of norm for convergence check of primary variables in EHL",
                  .default_value = convnorm_abs}),


          deprecated_selection<BinaryOp>("NORMCOMBI_RESFINC",
              {
                  {"And", bop_and},
                  {"Or", bop_or},
                  {"Coupl_Or_Single", bop_coupl_or_single},
                  {"Coupl_And_Single", bop_coupl_and_single},
                  {"And_Single", bop_and_single},
                  {"Or_Single", bop_or_single},
              },
              {.description =
                      "binary operator to combine primary variables and residual force values",
                  .default_value = bop_coupl_and_single}),

          deprecated_selection<VectorNorm>("ITERNORM",
              {
                  {"L1", norm_l1},
                  {"L1_Scaled", norm_l1_scaled},
                  {"L2", norm_l2},
                  {"Rms", norm_rms},
                  {"Inf", norm_inf},
              },
              {.description = "type of norm to be applied to residuals",
                  .default_value = norm_rms}),


          parameter<double>(
              "PTCDT", {.description = "pseudo time step for pseudo-transient "
                                       "continuation (PTC) stabilised Newton procedure",
                           .default_value = 0.1}),

          // number of linear solver used for monolithic EHL
          parameter<int>("LINEAR_SOLVER",
              {.description = "number of linear solver used for monolithic EHL problems",
                  .default_value = -1}),

          // convergence criteria adaptivity of monolithic EHL solver
          parameter<bool>("ADAPTCONV", {.description = "Switch on adaptive control of linear "
                                                       "solver tolerance for nonlinear solution",
                                           .default_value = false}),
          parameter<double>("ADAPTCONV_BETTER",
              {.description =
                      "The linear solver shall be this much better than the current nonlinear "
                      "residual in the nonlinear convergence limit",
                  .default_value = 0.1}),

          parameter<bool>("INFNORMSCALING",
              {.description = "Scale blocks of matrix with row infnorm?", .default_value = true})},
      {.required = false}));

  /*----------------------------------------------------------------------*/
  /* parameters for partitioned EHL */
  /*----------------------------------------------------------------------*/
  specs.push_back(group("ELASTO HYDRO DYNAMIC/PARTITIONED",
      {

          // Solver parameter for relaxation of iterative staggered partitioned EHL
          parameter<double>(
              "MAXOMEGA", {.description = "largest omega allowed for Aitken relaxation",
                              .default_value = 10.0}),
          parameter<double>(
              "MINOMEGA", {.description = "smallest omega allowed for Aitken relaxation",
                              .default_value = 0.1}),
          parameter<double>(
              "STARTOMEGA", {.description = "fixed relaxation parameter", .default_value = 1.0}),

          // convergence tolerance of outer iteration loop
          parameter<double>("CONVTOL",
              {.description =
                      "tolerance for convergence check of outer iteration within partitioned EHL",
                  .default_value = 1e-6})},
      {.required = false}));
  return specs;
}


void EHL::set_valid_conditions(std::vector<Core::Conditions::ConditionDefinition>& condlist)
{
  using namespace Core::IO::InputSpecBuilders;
  /*--------------------------------------------------------------------*/
  // ehl mortar coupling

  Core::Conditions::ConditionDefinition lineehl("DESIGN LINE EHL MORTAR COUPLING CONDITIONS 2D",
      "EHLCoupling", "Line EHL (elasto-hydro-dynamic) Coupling", Core::Conditions::EHLCoupling,
      true, Core::Conditions::geometry_type_line);
  Core::Conditions::ConditionDefinition surfehl("DESIGN SURF EHL MORTAR COUPLING CONDITIONS 3D",
      "EHLCoupling", "Surface EHL (elasto-hydro-dynamic) Coupling", Core::Conditions::EHLCoupling,
      true, Core::Conditions::geometry_type_surface);

  const auto make_ehl_cond = [&condlist](Core::Conditions::ConditionDefinition& cond)
  {
    cond.add_component(parameter<int>("InterfaceID"));
    cond.add_component(deprecated_selection<std::string>(
        "Side", {"Master", "Slave"}, {.description = "interface side"}));
    cond.add_component(deprecated_selection<std::string>("Initialization", {"Inactive", "Active"},
        {.description = "initialization", .default_value = "Active"}));
    cond.add_component(
        parameter<double>("FrCoeffOrBound", {.description = "", .default_value = 0.0}));

    condlist.push_back(cond);
  };

  make_ehl_cond(lineehl);
  make_ehl_cond(surfehl);
}

FOUR_C_NAMESPACE_CLOSE