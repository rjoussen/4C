// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_lubrication_dyn.hpp"

#include "4C_fem_dofset_predefineddofnumber.hpp"
#include "4C_fem_general_utils_createdis.hpp"
#include "4C_global_data.hpp"
#include "4C_lubrication_adapter.hpp"
#include "4C_lubrication_timint_implicit.hpp"

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*
 | Main control routine for Lubrication problems            wirtz 11/15 |
 *----------------------------------------------------------------------*/
void lubrication_dyn(int restart)
{
  // access the communicator
  MPI_Comm comm = Global::Problem::instance()->get_dis("lubrication")->get_comm();

  // print problem type
  if (Core::Communication::my_mpi_rank(comm) == 0)
  {
    std::cout << "###################################################" << std::endl;
    std::cout << "# YOUR PROBLEM TYPE: " << Global::Problem::instance()->problem_name()
              << std::endl;
    std::cout << "###################################################" << std::endl;
  }

  // access the problem-specific parameter list
  const Teuchos::ParameterList& lubricationdyn =
      Global::Problem::instance()->lubrication_dynamic_params();

  // access the lubrication discretization
  std::shared_ptr<Core::FE::Discretization> lubricationdis =
      Global::Problem::instance()->get_dis("lubrication");

  lubricationdis->fill_complete();

  // we directly use the elements from the Lubrication elements section
  if (lubricationdis->num_global_nodes() == 0)
    FOUR_C_THROW("No elements in the ---LUBRICATION ELEMENTS section");

  // add proxy of velocity related degrees of freedom to lubrication discretization
  std::shared_ptr<Core::DOFSets::DofSetInterface> dofsetaux =
      std::make_shared<Core::DOFSets::DofSetPredefinedDoFNumber>(
          Global::Problem::instance()->n_dim(), 0, 0, true);
  if (lubricationdis->add_dof_set(dofsetaux) != 1)
    FOUR_C_THROW("lub discretization has illegal number of dofsets!");

  // finalize discretization
  lubricationdis->fill_complete({
      .assign_degrees_of_freedom = true,
      .init_elements = false,
      .do_boundary_conditions = false,
  });

  // get linear solver id from LUBRICATION DYNAMIC
  const int linsolvernumber = lubricationdyn.get<int>("LINEAR_SOLVER");
  if (linsolvernumber == (-1))
    FOUR_C_THROW(
        "no linear solver defined for LUBRICATION problem. Please set LINEAR_SOLVER in LUBRICATION "
        "DYNAMIC to a valid number!");

  // create instance of Lubrication basis algorithm
  std::shared_ptr<Lubrication::LubricationBaseAlgorithm> lubricationonly =
      std::make_shared<Lubrication::LubricationBaseAlgorithm>();

  // setup Lubrication basis algorithm
  lubricationonly->setup(
      lubricationdyn, lubricationdyn, Global::Problem::instance()->solver_params(linsolvernumber));

  // read the restart information, set vectors and variables
  if (restart) lubricationonly->lubrication_field()->read_restart(restart);

  // enter time loop to solve problem
  (lubricationonly->lubrication_field())->time_loop();

  // perform the result test if required
  Global::Problem::instance()->add_field_test(lubricationonly->create_lubrication_field_test());
  Global::Problem::instance()->test_all(comm);

  return;

}  // end of lubrication_dyn()

FOUR_C_NAMESPACE_CLOSE
