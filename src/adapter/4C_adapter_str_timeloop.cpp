// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_adapter_str_timeloop.hpp"

#include "4C_comm_mpi_utils.hpp"
#include "4C_fem_discretization.hpp"
#include "4C_global_data.hpp"
#include "4C_inpar_structure.hpp"
#include "4C_io.hpp"
#include "4C_io_pstream.hpp"
#include "4C_material_time_step_request.hpp"

#include <Teuchos_StandardParameterEntryValidators.hpp>

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Adapter::StructureTimeLoop::StructureTimeLoop(
    Global::Problem& problem, std::shared_ptr<Structure> structure)
    : StructureWrapper(structure), problem_(problem)
{
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
int Adapter::StructureTimeLoop::integrate()
{
  // error checking variables
  Inpar::Solid::ConvergenceStatus convergencestatus = Inpar::Solid::conv_success;

  // target time #timen_ and step #stepn_ already set
  // time loop
  while (not_finished() and (convergencestatus == Inpar::Solid::conv_success or
                                convergencestatus == Inpar::Solid::conv_fail_repeat))
  {
    Core::Mat::TimeStepReduction::AlgorithmSupportScope support_scope(
        supports_material_time_step_reduction());
    apply_scheduled_time_step_increase();

    // call the predictor
    pre_predict();
    const StepControlResult prepare_result = prepare_time_step_control_result();
    if (handle_step_control_result(prepare_result))
    {
      convergencestatus = prepare_result.convergence_status;
      continue;
    }

    // integrate time step, i.e. do corrector steps
    // after this step we hold disn_, etc
    pre_solve();
    const StepControlResult solve_result = post_solve_control_result(solve());
    convergencestatus = solve_result.convergence_status;
    if (handle_step_control_result(solve_result))
    {
      continue;
    }

    // if everything is fine
    if (convergencestatus == Inpar::Solid::conv_success)
    {
      // calculate stresses, strains and energies
      // note: this has to be done before the update since otherwise a potential
      // material history is overwritten
      constexpr bool force_prepare = false;
      prepare_output(force_prepare);

      // update displacements, velocities, accelerations
      // after this call we will have disn_==dis_, etc
      // update time and step
      // update everything on the element level
      pre_update();
      update();
      post_update();

      // write output
      output();
      post_output();

      // print info about finished time step
      print_step();
    }
  }

  post_time_loop();

  // that's it say what went wrong
  return convergencestatus;
}

bool Adapter::StructureTimeLoop::handle_step_control_result(const StepControlResult& result)
{
  switch (result.action)
  {
    case Adapter::StepControlAction::proceed:
      return false;
    case Adapter::StepControlAction::repeat_step:
      switch (result.retry_reason)
      {
        case Adapter::RetryStepReason::none:
          return true;
        case Adapter::RetryStepReason::material_time_step_reduction:
        case Adapter::RetryStepReason::nonlinear_solver_time_step_reduction:
          apply_retry_time_step(result);
          return true;
        default:
          FOUR_C_THROW("Unknown retry-step reason.");
      }
    case Adapter::StepControlAction::stop:
      return true;
    default:
      FOUR_C_THROW("Unknown step-control action.");
  }
}

void Adapter::StructureTimeLoop::apply_retry_time_step(const StepControlResult& result)
{
  if (not result.proposed_time_step.has_value())
  {
    FOUR_C_THROW(
        "A step retry with reduced time step was requested at time t = {}, but no replacement "
        "time-step size was provided by the structure algorithm.",
        time());
  }

  const double old_dt = dt();
  const double new_dt = *result.proposed_time_step;

  if (Core::Communication::my_mpi_rank(discretization()->get_comm()) == 0)
  {
    Core::IO::cout << std::string(72, '*') << "\n"
                   << "Material evaluation requested a reduced time step at time t= " << time()
                   << "\n"
                   << "Old time step: " << old_dt << "\n"
                   << "New time step: " << new_dt << "\n"
                   << std::string(72, '*') << "\n"
                   << Core::IO::endl;
  }

  reset_step_for_time_step_retry(new_dt);
}

FOUR_C_NAMESPACE_CLOSE
