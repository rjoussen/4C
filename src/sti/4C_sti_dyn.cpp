// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_sti_dyn.hpp"

#include "4C_adapter_scatra_base_algorithm.hpp"
#include "4C_fem_dofset_predefineddofnumber.hpp"
#include "4C_fem_general_utils_createdis.hpp"
#include "4C_global_data.hpp"
#include "4C_scatra_resulttest_elch.hpp"
#include "4C_scatra_timint_elch.hpp"
#include "4C_sti_clonestrategy.hpp"
#include "4C_sti_monolithic.hpp"
#include "4C_sti_partitioned.hpp"
#include "4C_sti_resulttest.hpp"

#include <Teuchos_StandardParameterEntryValidators.hpp>
#include <Teuchos_TimeMonitor.hpp>

FOUR_C_NAMESPACE_OPEN

/*--------------------------------------------------------------------------------*
 | entry point for simulations of scalar-thermo interaction problems   fang 04/15 |
 *--------------------------------------------------------------------------------*/
void sti_dyn(const int& restartstep  //! time step for restart
)
{
  // access global problem
  Global::Problem* problem = Global::Problem::instance();

  // access communicator
  MPI_Comm comm = problem->get_dis("scatra")->get_comm();

  // access scatra discretization
  std::shared_ptr<Core::FE::Discretization> scatradis = problem->get_dis("scatra");

  // add dofset for velocity-related quantities to scatra discretization
  std::shared_ptr<Core::DOFSets::DofSetInterface> dofsetaux =
      std::make_shared<Core::DOFSets::DofSetPredefinedDoFNumber>(problem->n_dim() + 1, 0, 0, true);
  if (scatradis->add_dof_set(dofsetaux) != 1)
    FOUR_C_THROW("Scatra discretization has illegal number of dofsets!");

  // finalize scatra discretization
  scatradis->fill_complete();

  // safety check
  if (scatradis->num_global_nodes() == 0)
    FOUR_C_THROW(
        "The scatra discretization must not be empty, since the thermo discretization needs to be "
        "cloned from it!");

  // access thermo discretization
  std::shared_ptr<Core::FE::Discretization> thermodis = problem->get_dis("thermo");

  // add dofset for velocity-related quantities to thermo discretization
  dofsetaux =
      std::make_shared<Core::DOFSets::DofSetPredefinedDoFNumber>(problem->n_dim() + 1, 0, 0, true);
  if (thermodis->add_dof_set(dofsetaux) != 1)
    FOUR_C_THROW("Thermo discretization has illegal number of dofsets!");

  // equip thermo discretization with noderowmap for subsequent safety check
  // final fill_complete() is called at the end of discretization cloning
  thermodis->fill_complete(false, false, false);

  // safety check
  if (thermodis->num_global_nodes() != 0)
    FOUR_C_THROW(
        "The thermo discretization must be empty, since it is cloned from the scatra "
        "discretization!");

  // clone thermo discretization from scatra discretization, using clone strategy for scatra-thermo
  // interaction
  Core::FE::clone_discretization<STI::ScatraThermoCloneStrategy>(
      *scatradis, *thermodis, Global::Problem::instance()->cloning_material_map());
  thermodis->fill_complete(false, true, true);

  // add proxy of scalar transport degrees of freedom to thermo discretization and vice versa
  if (thermodis->add_dof_set(scatradis->get_dof_set_proxy()) != 2)
    FOUR_C_THROW("Thermo discretization has illegal number of dofsets!");
  if (scatradis->add_dof_set(thermodis->get_dof_set_proxy()) != 2)
    FOUR_C_THROW("Scatra discretization has illegal number of dofsets!");
  thermodis->fill_complete(true, false, false);
  scatradis->fill_complete(true, false, false);

  // add material of scatra elements to thermo elements and vice versa
  for (int i = 0; i < scatradis->num_my_col_elements(); ++i)
  {
    Core::Elements::Element* scatraele = scatradis->l_col_element(i);
    Core::Elements::Element* thermoele = thermodis->g_element(scatraele->id());

    thermoele->add_material(scatraele->material());
    scatraele->add_material(thermoele->material());
  }

  // access parameter lists for scatra-thermo interaction and scalar transport field
  const Teuchos::ParameterList& stidyn = problem->sti_dynamic_params();
  const Teuchos::ParameterList& scatradyn = problem->scalar_transport_dynamic_params();

  // extract and check ID of linear solver for scatra field
  const int solver_id_scatra = scatradyn.get<int>("LINEAR_SOLVER");
  if (solver_id_scatra == -1)
    FOUR_C_THROW(
        "No linear solver for scalar transport field was specified in input file section 'SCALAR "
        "TRANSPORT DYNAMIC'!");

  // extract and check ID of linear solver for thermo field
  const int solver_id_thermo = stidyn.get<int>("THERMO_LINEAR_SOLVER");
  if (solver_id_thermo == -1)
    FOUR_C_THROW(
        "No linear solver for temperature field was specified in input file section 'STI "
        "DYNAMIC'!");

  // instantiate coupling algorithm for scatra-thermo interaction
  std::shared_ptr<STI::Algorithm> sti_algorithm(nullptr);

  switch (Teuchos::getIntegralValue<STI::CouplingType>(stidyn, "COUPLINGTYPE"))
  {
    // monolithic algorithm
    case STI::CouplingType::monolithic:
    {
      // extract and check ID of monolithic linear solver
      const int solver_id = stidyn.sublist("MONOLITHIC").get<int>("LINEAR_SOLVER");
      if (solver_id == -1)
        FOUR_C_THROW(
            "No global linear solver was specified in input file section 'STI "
            "DYNAMIC/MONOLITHIC'!");

      sti_algorithm = std::make_shared<STI::Monolithic>(comm, stidyn, scatradyn,
          Global::Problem::instance()->solver_params(solver_id),
          Global::Problem::instance()->solver_params(solver_id_scatra),
          Global::Problem::instance()->solver_params(solver_id_thermo));

      break;
    }

    // partitioned algorithm
    case STI::CouplingType::oneway_scatratothermo:
    case STI::CouplingType::oneway_thermotoscatra:
    case STI::CouplingType::twoway_scatratothermo:
    case STI::CouplingType::twoway_scatratothermo_aitken:
    case STI::CouplingType::twoway_scatratothermo_aitken_dofsplit:
    case STI::CouplingType::twoway_thermotoscatra:
    case STI::CouplingType::twoway_thermotoscatra_aitken:
    {
      sti_algorithm = std::make_shared<STI::Partitioned>(comm, stidyn, scatradyn,
          Global::Problem::instance()->solver_params(solver_id_scatra),
          Global::Problem::instance()->solver_params(solver_id_thermo));

      break;
    }

    // unknown algorithm
    default:
    {
      FOUR_C_THROW("Unknown coupling algorithm for scatra-thermo interaction!");
    }
  }

  // read restart data if necessary
  if (restartstep) sti_algorithm->read_restart(restartstep);

  // provide scatra and thermo fields with velocities
  sti_algorithm->scatra_field()->set_velocity_field_from_function();
  sti_algorithm->thermo_field()->set_velocity_field_from_function();

  // enter time loop and solve scatra-thermo interaction problem
  sti_algorithm->time_loop();

  // summarize performance measurements
  Teuchos::TimeMonitor::summarize();

  // perform result tests
  problem->add_field_test(
      std::shared_ptr<Core::Utils::ResultTest>(new STI::STIResultTest(sti_algorithm)));
  if (Teuchos::getIntegralValue<STI::ScaTraTimIntType>(
          problem->sti_dynamic_params(), "SCATRATIMINTTYPE") == STI::ScaTraTimIntType::elch)
    problem->add_field_test(std::shared_ptr<Core::Utils::ResultTest>(new ScaTra::ElchResultTest(
        std::dynamic_pointer_cast<ScaTra::ScaTraTimIntElch>(sti_algorithm->scatra_field()))));
  else
    FOUR_C_THROW(
        "Scatra-thermo interaction is currently only available for thermodynamic electrochemistry, "
        "but not for other kinds of thermodynamic scalar transport!");
  problem->add_field_test(std::shared_ptr<Core::Utils::ResultTest>(
      new ScaTra::ScaTraResultTest(sti_algorithm->thermo_field())));
  problem->test_all(comm);

  return;
}  // sti_dyn()

FOUR_C_NAMESPACE_CLOSE
