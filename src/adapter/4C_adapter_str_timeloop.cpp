// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_adapter_str_timeloop.hpp"

#include "4C_global_data.hpp"
#include "4C_structure_new_input.hpp"
#include "4C_structure_new_timint_base.hpp"

#include <Teuchos_StandardParameterEntryValidators.hpp>

#include <utility>

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Adapter::StructureTimeLoop::StructureTimeLoop(
    Global::Problem& problem, std::shared_ptr<Structure> structure)
    : StructureWrapper(std::move(structure)), problem_(problem)
{
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Adapter::StructureTimeLoop::integrate()
{
  // error checking variables
  Solid::StepStatus solve_status = Solid::StepStatus::no_errors;

  // target time #timen_ and step #stepn_ already set
  // time loop
  while (not_finished() and (solve_status == Solid::StepStatus::no_errors or
                                solve_status == Solid::StepStatus::fail_repeat))
  {
    // call the predictor
    prepare_time_step();

    // integrate time step, i.e. do corrector steps
    // after this step we hold disn_, etc
    solve_status = solve();

    // if everything is fine
    if (solve_status == Solid::StepStatus::no_errors)
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
      // dynamic rebalance is only implemented for the new time integration
      if (auto* timint = dynamic_cast<Solid::TimeInt::Base*>(structure_.get()))
        if (timint->should_perform_dynamic_rebalance()) timint->perform_dynamic_rebalance();

      // print info about finished time step
      print_step();
    }
    // todo: remove this as soon as old structure time integration is gone
    else if (Teuchos::getIntegralValue<Solid::IntegrationStrategy>(
                 problem_.structural_dynamic_params(), "INT_STRATEGY") == Solid::int_old)
    {
      solve_status = perform_error_action(solve_status);  // something went wrong update error code
                                                          // according to chosen divcont action
    }
  }

  post_time_loop();
}

FOUR_C_NAMESPACE_CLOSE
