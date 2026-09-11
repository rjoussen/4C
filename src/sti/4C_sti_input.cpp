// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_sti_input.hpp"

#include "4C_fem_condition_definition.hpp"
#include "4C_io_input_spec_builders.hpp"
#include "4C_linalg_sparseoperator.hpp"
#include "4C_scatra_input.hpp"
FOUR_C_NAMESPACE_OPEN


std::vector<Core::IO::InputSpec> STI::valid_parameters()
{
  using namespace Core::IO::InputSpecBuilders;
  using namespace Core::IO::InputSpecBuilders::Validators;

  std::vector<Core::IO::InputSpec> specs;
  specs.push_back(group("STI DYNAMIC",
      {

          // type of scalar transport time integration
          deprecated_selection<ScaTraTimIntType>("SCATRATIMINTTYPE",
              {
                  {"Standard", ScaTraTimIntType::standard},
                  {"Elch", ScaTraTimIntType::elch},
              },
              {.description =
                      "scalar transport time integration type is needed to instantiate correct "
                      "scalar "
                      "transport time integration scheme for scatra-thermo interaction problems",
                  .default_value = ScaTraTimIntType::standard}),

          // type of coupling between scatra and thermo fields
          deprecated_selection<CouplingType>("COUPLINGTYPE",
              {
                  {"Undefined", CouplingType::undefined},
                  {"Monolithic", CouplingType::monolithic},
                  {"OneWay_ScatraToThermo", CouplingType::oneway_scatratothermo},
                  {"OneWay_ThermoToScatra", CouplingType::oneway_thermotoscatra},
                  {"TwoWay_ScatraToThermo", CouplingType::twoway_scatratothermo},
                  {"TwoWay_ScatraToThermo_Aitken", CouplingType::twoway_scatratothermo_aitken},
                  {"TwoWay_ScatraToThermo_Aitken_Dofsplit",
                      CouplingType::twoway_scatratothermo_aitken_dofsplit},
                  {"TwoWay_ThermoToScatra", CouplingType::twoway_thermotoscatra},
                  {"TwoWay_ThermoToScatra_Aitken", CouplingType::twoway_thermotoscatra_aitken},
              },
              {.description = "type of coupling between scatra and thermo fields",
                  .default_value = CouplingType::undefined}),

          // specification of initial temperature field
          parameter<ScaTra::InitialField>("THERMO_INITIALFIELD",
              {.description = "Initial Field for scatra-thermo interaction problems",
                  .default_value = ScaTra::InitialField::zero_field,
                  .validator = in_set<ScaTra::InitialField>(
                      {ScaTra::InitialField::zero_field, ScaTra::InitialField::field_by_function,
                          ScaTra::InitialField::field_by_condition})}),

          // function number for initial temperature field
          parameter<int>("THERMO_INITFUNCNO",
              {.description = "function number for initial temperature field for "
                              "scatra-thermo interaction problems",
                  .default_value = -1}),

          // ID of linear solver for temperature field
          parameter<int>("THERMO_LINEAR_SOLVER",
              {.description = "ID of linear solver for temperature field", .default_value = -1}),

          // flag for double condensation of linear equations associated with temperature field
          parameter<bool>("THERMO_CONDENSATION",
              {.description = "flag for double condensation of linear equations associated with "
                              "temperature field",
                  .default_value = false})},
      {.required = false}));
  /*----------------------------------------------------------------------*/
  // valid parameters for monolithic scatra-thermo interaction
  specs.push_back(group("STI DYNAMIC/MONOLITHIC",
      {

          // ID of linear solver for global system of equations
          parameter<int>(
              "LINEAR_SOLVER", {.description = "ID of linear solver for global system of equations",
                                   .default_value = -1}),

          // type of global system matrix in global system of equations
          deprecated_selection<Core::LinAlg::MatrixType>("MATRIXTYPE",
              {
                  {"block", Core::LinAlg::MatrixType::block_condition},
                  {"sparse", Core::LinAlg::MatrixType::sparse},
              },
              {.description = "type of global system matrix in global system of equations",
                  .default_value = Core::LinAlg::MatrixType::block_condition})},
      {.required = false}));

  /*----------------------------------------------------------------------*/
  // valid parameters for partitioned scatra-thermo interaction
  specs.push_back(group("STI DYNAMIC/PARTITIONED",
      {

          // relaxation parameter
          parameter<double>("OMEGA", {.description = "relaxation parameter", .default_value = 1.}),

          // maximum value of Aitken relaxation parameter
          parameter<double>("OMEGAMAX",
              {.description = "maximum value of Aitken relaxation parameter (0.0 = no constraint)",
                  .default_value = 0.})},
      {.required = false}));
  return specs;
}


FOUR_C_NAMESPACE_CLOSE