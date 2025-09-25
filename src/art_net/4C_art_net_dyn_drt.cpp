// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_art_net_dyn_drt.hpp"

#include "4C_adapter_art_net.hpp"
#include "4C_art_net_artery_resulttest.hpp"
#include "4C_art_net_utils.hpp"
#include "4C_fem_general_utils_createdis.hpp"
#include "4C_global_data.hpp"
#include "4C_io_control.hpp"
#include "4C_utils_result_test.hpp"

#include <Teuchos_StandardParameterEntryValidators.hpp>
#include <Teuchos_TimeMonitor.hpp>

#include <cstdlib>
#include <ctime>
#include <iostream>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 * Main control routine for arterial network including various solvers:
 *
 *        o
 *
 *----------------------------------------------------------------------*/
void dyn_art_net_drt() { dyn_art_net_drt(false); }

std::shared_ptr<Adapter::ArtNet> dyn_art_net_drt(bool CoupledTo3D)
{
  if (Global::Problem::instance()->does_exist_dis("artery") == false)
  {
    return nullptr;
  }

  // define the discretization names
  const std::string artery_disname = "artery";
  const std::string scatra_disname = "artery_scatra";

  // access the problem
  Global::Problem* problem = Global::Problem::instance();

  // -------------------------------------------------------------------
  // access the discretization
  // -------------------------------------------------------------------
  std::shared_ptr<Core::FE::Discretization> actdis = nullptr;

  actdis = problem->get_dis(artery_disname);

  // -------------------------------------------------------------------
  // set degrees of freedom in the discretization
  // -------------------------------------------------------------------
  if (!actdis->filled()) actdis->fill_complete();

  // -------------------------------------------------------------------
  // If discretization is empty, then return empty time integration
  // -------------------------------------------------------------------
  if (actdis->num_global_elements() < 1)
  {
    return nullptr;
  }

  // -------------------------------------------------------------------
  // context for output and restart
  // -------------------------------------------------------------------
  std::shared_ptr<Core::IO::DiscretizationWriter> output = actdis->writer();
  output->write_mesh(0, 0.0);

  // -------------------------------------------------------------------
  // set some pointers and variables
  // -------------------------------------------------------------------
  const Teuchos::ParameterList& artdyn = problem->arterial_dynamic_params();

  // -------------------------------------------------------------------
  // create a solver
  // -------------------------------------------------------------------
  // get the solver number
  const int linsolvernumber = artdyn.get<int>("LINEAR_SOLVER");
  // check if the solver has a valid solver number
  if (linsolvernumber == (-1))
    FOUR_C_THROW(
        "no linear solver defined. Please set LINEAR_SOLVER in ARTERIAL DYNAMIC to a valid "
        "number!");

  // solution output
  if (artdyn.get<bool>("SOLVESCATRA"))
  {
    if (Core::Communication::my_mpi_rank(actdis->get_comm()) == 0)
    {
      std::cout << "<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>" << std::endl;
      std::cout << "<  ARTERY:  ScaTra coupling present  >" << std::endl;
      std::cout << "<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>" << std::endl;
    }
    std::shared_ptr<Core::FE::Discretization> scatradis = problem->get_dis(scatra_disname);
    // fill scatra discretization by cloning artery discretization
    Core::FE::clone_discretization<Arteries::ArteryScatraCloneStrategy>(
        *actdis, *scatradis, Global::Problem::instance()->cloning_material_map());
    scatradis->fill_complete();

    // the problem is one way coupled, scatra needs only artery

    // build a proxy of the structure discretization for the scatra field
    std::shared_ptr<Core::DOFSets::DofSetInterface> arterydofset = actdis->get_dof_set_proxy();

    // check if ScatraField has 2 discretizations, so that coupling is possible
    if (scatradis->add_dof_set(arterydofset) != 1)
      FOUR_C_THROW("unexpected dof sets in scatra field");

    scatradis->fill_complete({
        .assign_degrees_of_freedom = true,
        .init_elements = false,
        .do_boundary_conditions = false,
    });
  }
  else
  {
    if (Core::Communication::my_mpi_rank(actdis->get_comm()) == 0)
    {
      std::cout << "<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>" << std::endl;
      std::cout << "<  ARTERY: no ScaTra coupling present  >" << std::endl;
      std::cout << "<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>" << std::endl;
    }
  }

  // -------------------------------------------------------------------
  // algorithm construction depending on
  // time-integration (or stationary) scheme
  // -------------------------------------------------------------------
  auto timintscheme =
      Teuchos::getIntegralValue<ArtDyn::TimeIntegrationScheme>(artdyn, "DYNAMICTYPE");

  // build art net time integrator
  std::shared_ptr<Adapter::ArtNet> artnettimint = Arteries::Utils::create_algorithm(
      timintscheme, actdis, linsolvernumber, artdyn, artdyn, *output);

  // initialize
  artnettimint->init(artdyn, artdyn, scatra_disname);

  // Initialize state save vectors
  if (CoupledTo3D)
  {
    artnettimint->init_save_state();
  }

  // initial field from restart or calculated by given function
  const int restart = Global::Problem::instance()->restart();
  if (restart && !CoupledTo3D)
  {
    // read the restart information, set vectors and variables
    artnettimint->read_restart(restart);
  }
  else
  {
    // artnetexplicit.SetInitialData(init,startfuncno);
  }

  // assign materials
  // note: to be done after potential restart, as in read_restart()
  //       the secondary material is destroyed
  if (artdyn.get<bool>("SOLVESCATRA"))
    Arteries::Utils::assign_material_pointers(artery_disname, scatra_disname);

  if (!CoupledTo3D)
  {
    // call time-integration (or stationary) scheme
    std::shared_ptr<Teuchos::ParameterList> param_temp;
    artnettimint->integrate(CoupledTo3D, param_temp);

    // result test
    artnettimint->test_results();

    return artnettimint;
    //    return  nullptr;
  }
  else
  {
    return artnettimint;
  }

}  // end of dyn_art_net_drt()

FOUR_C_NAMESPACE_CLOSE
