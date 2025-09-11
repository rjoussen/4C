// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_inpar_fs3i.hpp"

#include "4C_inpar_scatra.hpp"
#include "4C_io_input_spec_builders.hpp"
FOUR_C_NAMESPACE_OPEN



std::vector<Core::IO::InputSpec> Inpar::FS3I::valid_parameters()
{
  using namespace Core::IO::InputSpecBuilders;

  std::vector<Core::IO::InputSpec> specs;
  specs.push_back(group("FS3I DYNAMIC",
      {

          parameter<double>("TIMESTEP", {.description = "Time increment dt", .default_value = 0.1}),
          parameter<int>(
              "NUMSTEP", {.description = "Total number of time steps", .default_value = 20}),
          parameter<double>(
              "MAXTIME", {.description = "Total simulation time", .default_value = 1000.0}),
          parameter<int>("RESULTSEVERY",
              {.description = "Increment for writing solution", .default_value = 1}),
          parameter<int>(
              "RESTARTEVERY", {.description = "Increment for writing restart", .default_value = 1}),
          deprecated_selection<Inpar::ScaTra::SolverType>("SCATRA_SOLVERTYPE",
              {
                  {"linear", Inpar::ScaTra::solvertype_linear_incremental},
                  {"nonlinear", Inpar::ScaTra::solvertype_nonlinear},
              },
              {.description = "type of scalar transport solver",
                  .default_value = Inpar::ScaTra::solvertype_nonlinear}),
          parameter<bool>(
              "INF_PERM", {.description = "Flag for infinite permeability", .default_value = true}),
          deprecated_selection<std::string>("CONSTHERMPRESS", {"No_energy", "No_mass", "Yes"},
              {.description = "treatment of thermodynamic pressure in time",
                  .default_value = "Yes"}),

          // number of linear solver used for fs3i problems
          parameter<int>("COUPLED_LINEAR_SOLVER",
              {.description = "number of linear solver used for fs3i problem",
                  .default_value = -1}),
          parameter<int>(
              "LINEAR_SOLVER1", {.description = "number of linear solver used for fluid problem",
                                    .default_value = -1}),
          parameter<int>("LINEAR_SOLVER2",
              {.description = "number of linear solver used for structural problem",
                  .default_value = -1}),

          deprecated_selection<Inpar::ScaTra::ConvForm>("STRUCTSCAL_CONVFORM",
              {
                  {"convective", Inpar::ScaTra::convform_convective},
                  {"conservative", Inpar::ScaTra::convform_conservative},
              },
              {.description = "form of convective term of structure scalar",
                  .default_value = Inpar::ScaTra::convform_conservative}),


          deprecated_selection<Inpar::ScaTra::InitialField>("STRUCTSCAL_INITIALFIELD",
              {
                  {"zero_field", Inpar::ScaTra::initfield_zero_field},
                  {"field_by_function", Inpar::ScaTra::initfield_field_by_function},
              },
              {.description = "Initial Field for structure scalar transport problem",
                  .default_value = Inpar::ScaTra::initfield_zero_field}),

          parameter<int>("STRUCTSCAL_INITFUNCNO",
              {.description = "function number for structure scalar transport initial field",
                  .default_value = -1}),

          // Type of coupling strategy between structure and structure-scalar field
          deprecated_selection<VolumeCoupling>("STRUCTSCAL_FIELDCOUPLING",
              {
                  {"volume_matching", coupling_match},
                  {"volume_nonmatching", coupling_nonmatch},
              },
              {.description =
                      "Type of coupling strategy between structure and structure-scalar field",
                  .default_value = coupling_match}),

          // Type of coupling strategy between fluid and fluid-scalar field
          deprecated_selection<VolumeCoupling>("FLUIDSCAL_FIELDCOUPLING",
              {
                  {"volume_matching", coupling_match},
                  {"volume_nonmatching", coupling_nonmatch},
              },
              {.description = "Type of coupling strategy between fluid and fluid-scalar field",
                  .default_value = coupling_match}),

          // type of scalar transport
          deprecated_selection<Inpar::ScaTra::ImplType>("FLUIDSCAL_SCATRATYPE",
              {
                  {"Undefined", Inpar::ScaTra::impltype_undefined},
                  {"ConvectionDiffusion", Inpar::ScaTra::impltype_std},
                  {"Loma", Inpar::ScaTra::impltype_loma},
                  {"Advanced_Reaction", Inpar::ScaTra::impltype_advreac},
                  {"Chemotaxis", Inpar::ScaTra::impltype_chemo},
                  {"Chemo_Reac", Inpar::ScaTra::impltype_chemoreac},
              },
              {.description = "Type of scalar transport problem",
                  .default_value = Inpar::ScaTra::impltype_std}),

          // Restart from FSI instead of FS3I
          parameter<bool>("RESTART_FROM_PART_FSI",
              {.description = "restart from partitioned fsi problem (e.g. from "
                              "prestress calculations) instead of fs3i",
                  .default_value = false}),

      },
      {.required = false}));

  /*----------------------------------------------------------------------*/
  /* parameters for partitioned FS3I */
  /*----------------------------------------------------------------------*/
  specs.push_back(group("FS3I DYNAMIC/PARTITIONED",
      {

          // Coupling strategy for partitioned FS3I
          parameter<SolutionSchemeOverFields>(
              "COUPALGO", {.description = "Coupling strategies for FS3I solvers",
                              .default_value = fs3i_IterStagg}),

          // convergence tolerance of outer iteration loop
          parameter<double>("CONVTOL",
              {.description =
                      "tolerance for convergence check of outer iteration within partitioned FS3I",
                  .default_value = 1e-6}),

          parameter<int>("ITEMAX",
              {.description = "Maximum number of outer iterations", .default_value = 10})},
      {.required = false}));

  /*----------------------------------------------------------------------  */
  /* parameters for stabilization of the structure-scalar field             */
  /*----------------------------------------------------------------------  */

  /// HACK!
  /// reuse the parameters from scatra
  specs.push_back(group("FS3I DYNAMIC/STRUCTURE SCALAR STABILIZATION",
      {Inpar::ScaTra::all_specs_for_scatra_stabilization()}, {.required = false}));
  return specs;
}

FOUR_C_NAMESPACE_CLOSE