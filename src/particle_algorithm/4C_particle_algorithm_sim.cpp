// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_particle_algorithm_sim.hpp"

#include "4C_comm_utils.hpp"
#include "4C_global_data.hpp"
#include "4C_particle_algorithm.hpp"

#include <Teuchos_RCPStdSharedPtrConversions.hpp>

#include <memory>

FOUR_C_NAMESPACE_OPEN

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
void particle_drt()
{
  // get instance of global problem
  Global::Problem* problem = Global::Problem::instance();

  // get local communicator
  MPI_Comm comm = problem->get_communicators().local_comm();

  // get parameter lists
  const Teuchos::ParameterList& params = problem->particle_params();

  // reference to vector of initial particles
  std::vector<PARTICLEENGINE::ParticleObjShrdPtr>& initialparticles = problem->particles();

  // create and init particle algorithm
  auto particlealgorithm = std::unique_ptr<PARTICLEALGORITHM::ParticleAlgorithm>(
      new PARTICLEALGORITHM::ParticleAlgorithm(comm, params));
  particlealgorithm->init(initialparticles);

  // read restart information
  const int restart = problem->restart();
  if (restart) particlealgorithm->read_restart(restart);

  // setup particle algorithm
  particlealgorithm->setup();

  // solve particle problem
  particlealgorithm->timeloop();

  // perform result tests
  {
    // create particle field specific result test objects
    std::vector<std::shared_ptr<Core::Utils::ResultTest>> allresulttests =
        particlealgorithm->create_result_tests();

    // add particle field specific result test objects
    for (auto& resulttest : allresulttests)
      if (resulttest) problem->add_field_test(resulttest);

    // perform all tests
    problem->test_all(comm);
  }

  // print summary statistics for all timers
  std::shared_ptr<const Teuchos::Comm<int>> TeuchosComm =
      Core::Communication::to_teuchos_comm<int>(comm);
  Teuchos::TimeMonitor::summarize(Teuchos::Ptr(TeuchosComm.get()), std::cout, false, true, false);
}

FOUR_C_NAMESPACE_CLOSE
