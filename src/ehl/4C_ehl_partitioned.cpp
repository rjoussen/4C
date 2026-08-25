// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_ehl_partitioned.hpp"

#include "4C_adapter_str_structure.hpp"
#include "4C_global_data.hpp"
#include "4C_linalg_utils_sparse_algebra_create.hpp"
#include "4C_lubrication_adapter.hpp"
#include "4C_lubrication_timint_implicit.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 | constructor                                              wirtz 12/15 |
 *----------------------------------------------------------------------*/
EHL::Partitioned::Partitioned(MPI_Comm comm, const Teuchos::ParameterList& globaltimeparams,
    const Teuchos::ParameterList& lubricationparams, const Teuchos::ParameterList& structparams,
    const std::string struct_disname, const std::string lubrication_disname)
    : Base(comm, globaltimeparams, lubricationparams, structparams, struct_disname,
          lubrication_disname),
      preincnp_(std::make_shared<Core::LinAlg::Vector<double>>(
          *lubrication_->lubrication_field()->discretization()->dof_row_map(0), true)),
      dispincnp_(std::make_shared<Core::LinAlg::Vector<double>>(*structure_->dof_row_map(0), true))
{
  // call the EHL parameter lists
  const Teuchos::ParameterList& ehlparams =
      Global::Problem::instance()->elasto_hydro_dynamic_params();
  const Teuchos::ParameterList& ehlparamspart =
      Global::Problem::instance()->elasto_hydro_dynamic_params().sublist("PARTITIONED");

  if (ehlparams.get<bool>("DIFFTIMESTEPSIZE"))
  {
    FOUR_C_THROW("Different time stepping for two way coupling not implemented yet.");
  }

  // Get the parameters for the convergence_check
  itmax_ = ehlparams.get<int>("ITEMAX");          // default: =10
  ittol_ = ehlparamspart.get<double>("CONVTOL");  // default: =1e-6

  // no dry contact in partitioned ehl
  if (dry_contact_) FOUR_C_THROW("no dry contact model in partitioned ehl");
}

/*----------------------------------------------------------------------*
 | Timeloop for EHL problems                                wirtz 12/15 |
 *----------------------------------------------------------------------*/
void EHL::Partitioned::timeloop()
{
  while (not_finished())
  {
    prepare_time_step();

    outer_loop();

    update_and_output();
  }
}


/*----------------------------------------------------------------------*
 | prepare time step                                        wirtz 12/15 |
 *----------------------------------------------------------------------*/
void EHL::Partitioned::prepare_time_step()
{
  increment_time_and_step();
  print_header();

  set_struct_solution(structure_->dispn());
  structure_->prepare_time_step();
  //  set_lubrication_solution(lubrication_->LubricationField()->Quantity()); // todo: what quantity
  lubrication_->lubrication_field()->prepare_time_step();
}


/*----------------------------------------------------------------------*
 | outer Timeloop for EHL without relaxation                wirtz 12/15 |
 *----------------------------------------------------------------------*/
void EHL::Partitioned::outer_loop()
{
  int itnum = 0;
  bool stopnonliniter = false;

  if (Core::Communication::my_mpi_rank(get_comm()) == 0)
  {
    std::cout << "\n****************************************\n          OUTER ITERATION "
                 "LOOP\n****************************************\n";
  }

  while (stopnonliniter == false)
  {
    itnum++;

    // store pressure from first solution for convergence check (like in
    // elch_algorithm: use current values)
    preincnp_->update(1.0, *lubrication_->lubrication_field()->prenp(), 0.0);
    dispincnp_->update(1.0, *structure_->dispnp(), 0.0);

    // set the external fluid force on the structure, which result from the fluid pressure
    set_lubrication_solution(lubrication_->lubrication_field()->prenp());
    if (itnum != 1) structure_->prepare_partition_step();
    // solve structural system
    do_struct_step();

    // set mesh displacement, velocity fields and film thickness
    set_struct_solution(structure_->dispnp());

    // solve lubrication equation and calculate the resulting traction, which will be applied on the
    // solids
    do_lubrication_step();
    // LubricationEvaluateSolveIterUpdate();

    // check convergence for all fields and stop iteration loop if
    // convergence is achieved overall
    stopnonliniter = convergence_check(itnum);
  }

  return;
}


/*----------------------------------------------------------------------*
 | constructor                                              wirtz 12/15 |
 *----------------------------------------------------------------------*/
void EHL::Partitioned::update_and_output()
{
  constexpr bool force_prepare = false;
  structure_->prepare_output(force_prepare);

  structure_->update();
  lubrication_->lubrication_field()->update();

  lubrication_->lubrication_field()->evaluate_error_compared_to_analytical_sol();

  structure_->output();
  lubrication_->lubrication_field()->output();
}


/*----------------------------------------------------------------------*
 | solve structure filed                                    wirtz 12/15 |
 *----------------------------------------------------------------------*/
void EHL::Partitioned::do_struct_step()
{
  if (Core::Communication::my_mpi_rank(get_comm()) == 0)
  {
    std::cout << "\n***********************\n STRUCTURE SOLVER \n***********************\n";
  }

  // Newton-Raphson iteration
  const auto structure_status = structure_->solve();
  FOUR_C_ASSERT_ALWAYS(
      structure_status == Solid::StepStatus::no_errors, "Structural solve failed.");
}


/*----------------------------------------------------------------------*
 | solve Lubrication field                                  wirtz 12/15 |
 *----------------------------------------------------------------------*/
void EHL::Partitioned::do_lubrication_step()
{
  if (Core::Communication::my_mpi_rank(get_comm()) == 0)
  {
    std::cout << "\n***********************\n  LUBRICATION SOLVER \n***********************\n";
  }

  // -------------------------------------------------------------------
  //                           solve nonlinear
  // -------------------------------------------------------------------
  lubrication_->lubrication_field()->solve();
}


/*----------------------------------------------------------------------*
 | convergence check for both fields (lubrication & structure)          |
 |                                                          wirtz 12/15 |
 *----------------------------------------------------------------------*/
bool EHL::Partitioned::convergence_check(int itnum)
{
  // convergence check based on the pressure increment
  bool stopnonliniter = false;

  //    | pressure increment |_2
  //  -------------------------------- < Tolerance
  //     | pressure+1 |_2
  //
  // AND
  //
  //    | pressure increment |_2
  //  -------------------------------- < Tolerance
  //             dt * n

  // variables to save different L2 - Norms
  // define L2-norm of incremental pressure and pressure
  double preincnorm_L2(0.0);
  double prenorm_L2(0.0);
  double dispincnorm_L2(0.0);
  double dispnorm_L2(0.0);

  // build the current pressure increment Inc T^{i+1}
  // \f Delta T^{k+1} = Inc T^{k+1} = T^{k+1} - T^{k}  \f
  preincnp_->update(1.0, *(lubrication_->lubrication_field()->prenp()), -1.0);
  dispincnp_->update(1.0, *(structure_->dispnp()), -1.0);

  // build the L2-norm of the pressure increment and the pressure
  preincnp_->norm_2(&preincnorm_L2);
  lubrication_->lubrication_field()->prenp()->norm_2(&prenorm_L2);
  dispincnp_->norm_2(&dispincnorm_L2);
  structure_->dispnp()->norm_2(&dispnorm_L2);

  // care for the case that there is (almost) zero pressure
  if (prenorm_L2 < 1e-6) prenorm_L2 = 1.0;
  if (dispnorm_L2 < 1e-6) dispnorm_L2 = 1.0;

  // print the incremental based convergence check to the screen
  if (Core::Communication::my_mpi_rank(get_comm()) == 0)
  {
    std::cout << "\n";
    std::cout
        << "***********************************************************************************\n";
    std::cout << "    OUTER ITERATION STEP    \n";
    std::cout
        << "***********************************************************************************\n";
    printf(
        "+--------------+---------------------+------------------+-----------------+---------------"
        "-------+------------------+\n");
    printf(
        "|-  step/max  -|-  tol      [norm]  -|-  pressure-inc  -|  disp-inc      -|-  "
        "pressure-rel-inc  -|-  disp-rel-inc  -|\n");
    printf(
        "|   %3d/%3d    |  %10.3E[L_2 ]   |  %10.3E      |  %10.3E     |  %10.3E          |  "
        "%10.3E      |",
        itnum, itmax_, ittol_, preincnorm_L2 / dt() / sqrt(preincnp_->global_length()),
        dispincnorm_L2 / dt() / sqrt(dispincnp_->global_length()), preincnorm_L2 / prenorm_L2,
        dispincnorm_L2 / dispnorm_L2);
    printf("\n");
    printf(
        "+--------------+---------------------+------------------+-----------------+---------------"
        "-------+------------------+\n");
  }

  // converged
  if (((preincnorm_L2 / prenorm_L2) <= ittol_) and ((dispincnorm_L2 / dispnorm_L2) <= ittol_) and
      ((dispincnorm_L2 / dt() / sqrt(dispincnp_->global_length())) <= ittol_) and
      ((preincnorm_L2 / dt() / sqrt(preincnp_->global_length())) <= ittol_))
  {
    stopnonliniter = true;
    if (Core::Communication::my_mpi_rank(get_comm()) == 0)
    {
      printf("\n");
      printf(
          "|  Outer Iteration loop converged after iteration %3d/%3d !                             "
          "                            |\n",
          itnum, itmax_);
      printf(
          "+--------------+---------------------+------------------+-----------------+-------------"
          "---------+------------------+\n");
    }
  }

  // stop if itemax is reached without convergence
  // timestep
  if ((itnum == itmax_) and
      (((preincnorm_L2 / prenorm_L2) > ittol_) or ((dispincnorm_L2 / dispnorm_L2) > ittol_) or
          ((dispincnorm_L2 / dt() / sqrt(dispincnp_->global_length())) > ittol_) or
          (preincnorm_L2 / dt() / sqrt(preincnp_->global_length())) > ittol_))
  {
    stopnonliniter = true;
    if ((Core::Communication::my_mpi_rank(get_comm()) == 0))
    {
      printf(
          "|     >>>>>> not converged in itemax steps!                                             "
          "                            |\n");
      printf(
          "+--------------+---------------------+------------------+-----------------+-------------"
          "---------+------------------+\n");
      printf("\n");
      printf("\n");
    }
    FOUR_C_THROW("The partitioned EHL solver did not converge in ITEMAX steps!");
  }

  return stopnonliniter;
}

FOUR_C_NAMESPACE_CLOSE
