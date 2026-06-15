// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_structure_new_timint_implicit.hpp"

#include "4C_fem_discretization.hpp"
#include "4C_inpar_structure.hpp"
#include "4C_io.hpp"
#include "4C_io_control.hpp"
#include "4C_io_pstream.hpp"
#include "4C_linalg_blocksparsematrix.hpp"
#include "4C_linalg_utils_sparse_algebra_io.hpp"
#include "4C_material_time_step_request.hpp"
#include "4C_solver_nonlin_nox_group.hpp"
#include "4C_solver_nonlin_nox_linearsystem.hpp"
#include "4C_solver_nonlin_nox_vector.hpp"
#include "4C_structure_new_impl_generic.hpp"
#include "4C_structure_new_nln_solver_factory.hpp"
#include "4C_structure_new_nln_solver_generic.hpp"
#include "4C_structure_new_predict_factory.hpp"
#include "4C_structure_new_predict_generic.hpp"
#include "4C_structure_new_timint_noxinterface.hpp"
#include "4C_structure_new_utils.hpp"

#include <NOX_Abstract_Group.H>

#include <algorithm>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Solid::TimeInt::Implicit::setup()
{
  // safety check
  check_init();

  Solid::TimeInt::Base::setup();
  validate_divergence_action_configuration();

  // ---------------------------------------------------------------------------
  // cast the base class integrator
  // ---------------------------------------------------------------------------
  implint_ptr_ = std::dynamic_pointer_cast<Solid::IMPLICIT::Generic>(integrator_ptr());

  // ---------------------------------------------------------------------------
  // build NOX interface
  // ---------------------------------------------------------------------------
  std::shared_ptr<Solid::TimeInt::NoxInterface> noxinterface_ptr =
      std::make_shared<Solid::TimeInt::NoxInterface>();
  noxinterface_ptr->init(
      data_global_state_ptr(), implint_ptr_, dbc_ptr(), Core::Utils::shared_ptr_from_ref(*this));
  noxinterface_ptr->setup();

  // ---------------------------------------------------------------------------
  // build predictor
  // ---------------------------------------------------------------------------
  const Inpar::Solid::PredEnum predtype = data_sdyn().get_predictor_type();
  predictor_ptr_ = Solid::Predict::build_predictor(predtype);
  predictor_ptr_->init(predtype, implint_ptr_, dbc_ptr(), data_global_state_ptr(), data_io_ptr(),
      data_sdyn().get_nox_params_ptr());
  predictor_ptr_->setup();

  // ---------------------------------------------------------------------------
  // build non-linear solver
  // ---------------------------------------------------------------------------
  const Inpar::Solid::NonlinSolTech nlnSolverType = data_sdyn().get_nln_solver_type();
  if (nlnSolverType == Inpar::Solid::soltech_singlestep)
    std::cout << "WARNING!!! You are trying to solve implicitly using the \"singlestep\" nonlinear "
                 "solver. This is not encouraged, since it only works for linear statics analysis. "
                 "Please use NLNSOL as \"fullnewton\" for reliable result."
              << std::endl;
  nlnsolver_ptr_ = Solid::Nln::SOLVER::build_nln_solver(nlnSolverType, data_global_state_ptr(),
      data_s_dyn_ptr(), noxinterface_ptr, implint_ptr_, Core::Utils::shared_ptr_from_ref(*this));

  // set setup flag
  issetup_ = true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Solid::TimeInt::Implicit::validate_divergence_action_configuration() const
{
  check_init_setup();

  switch (get_divergence_action())
  {
    case Inpar::Solid::divcont_stop:
    case Inpar::Solid::divcont_continue:
    case Inpar::Solid::divcont_halve_step:
    case Inpar::Solid::divcont_adapt_step:
    case Inpar::Solid::divcont_adapt_penaltycontact:
    case Inpar::Solid::divcont_repeat_simulation:
      return;
    case Inpar::Solid::divcont_repeat_step:
      FOUR_C_THROW(
          "DIVERCONT = repeat_step is not supported for structure_new implicit time integration. "
          "Use halve_step or adapt_step instead.");
    case Inpar::Solid::divcont_rand_adapt_step:
    case Inpar::Solid::divcont_rand_adapt_step_ele_err:
      FOUR_C_THROW(
          "DIVERCONT = rand_adapt_step and rand_adapt_step_ele_err are not supported for "
          "structure_new implicit time integration. Use halve_step or adapt_step instead.");
    default:
      FOUR_C_THROW("Unknown DIVERCONT case.");
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Solid::TimeInt::Implicit::set_state(const std::shared_ptr<Core::LinAlg::Vector<double>>& x)
{
  integrator_ptr()->set_state(*x);
  NOX::Nln::Vector x_nox(x, NOX::Nln::Vector::MemoryType::View);
  nln_solver().get_solution_group().setX(x_nox);
  set_state_in_sync_with_nox_group(true);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Solid::TimeInt::Implicit::prepare_partition_step()
{
  check_init_setup();
  FOUR_C_THROW("Not yet implemented!");
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Solid::TimeInt::Implicit::prepare_time_step()
{
  check_init_setup();
  // update end time \f$t_{n+1}\f$ of this time step to cope with time step size adaptivity
  /* ToDo Check if this is still necessary. I moved this part to the Update(const double endtime)
   * routine, such it becomes consistent with non-adaptive update routine! See the
   * update_step_time() routine for more information.                             hiermeier 12/2015
   *
  double& time_np = data_global_state().get_time_np();
  time_np = data_global_state().get_time_n() + (*data_global_state().get_delta_time())[0]; */

  const auto& status = prepare_time_step_with_status();
  if (status != Adapter::PrepareTimeStepStatus::successful)
  {
    FOUR_C_THROW("prepare_time_step() was not successful.");
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Adapter::PrepareTimeStepStatus Solid::TimeInt::Implicit::prepare_time_step_with_status()
{
  check_init_setup();
  ::NOX::Abstract::Group& grp = nln_solver().get_solution_group();
  retry_step_reason_ = Adapter::RetryStepReason::none;

  try
  {
    // Predictor internals may call into NOX deeply (e.g. compute/post_predict hooks). Those paths
    // must propagate MaterialTimeStepReductionRequestedFromNox unchanged so this seam can remap it
    // into status-based retry control flow.
    predictor().predict(grp);
    return Adapter::PrepareTimeStepStatus::successful;
  }
  catch (const Internal::MaterialTimeStepReductionRequestedFromNox&)
  {
    if (Core::Mat::TimeStepReduction::consume_mpi_reduced_request())
    {
      retry_step_reason_ = Adapter::RetryStepReason::material_time_step_reduction;
      return Adapter::PrepareTimeStepStatus::repeat_step;
    }
    return Adapter::PrepareTimeStepStatus::failed;
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Adapter::StepControlResult Solid::TimeInt::Implicit::prepare_time_step_control_result()
{
  Adapter::StepControlResult result;
  switch (prepare_time_step_with_status())
  {
    case Adapter::PrepareTimeStepStatus::successful:
      result.action = Adapter::StepControlAction::proceed;
      result.convergence_status = Inpar::Solid::conv_success;
      break;
    case Adapter::PrepareTimeStepStatus::repeat_step:
      result.action = Adapter::StepControlAction::repeat_step;
      result.convergence_status = Inpar::Solid::conv_fail_repeat;
      result.retry_reason = consume_retry_step_reason();
      if (result.retry_reason == Adapter::RetryStepReason::material_time_step_reduction)
      {
        result.proposed_time_step = compute_material_reduced_time_step();
      }
      break;
    case Adapter::PrepareTimeStepStatus::failed:
    default:
      result.action = Adapter::StepControlAction::stop;
      result.convergence_status = Inpar::Solid::conv_nonlin_fail;
      break;
  }
  return result;
}
/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
int Solid::TimeInt::Implicit::integrate()
{
  check_init_setup();
  FOUR_C_THROW(
      "The function is unused since the Adapter::StructureTimeLoop "
      "wrapper gives you all the flexibility you need.");
  return 0;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
int Solid::TimeInt::Implicit::integrate_step()
{
  check_init_setup();
  // do the predictor step
  ::NOX::Abstract::Group& grp = nln_solver().get_solution_group();
  predictor().predict(grp);
  return solve();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Inpar::Solid::ConvergenceStatus Solid::TimeInt::Implicit::solve()
{
  check_init_setup();
  throw_if_state_not_in_sync_with_nox_group();
  // reset the non-linear solver
  nln_solver().reset();
  // solve the non-linear problem
  return nln_solver().solve();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Solid::TimeInt::Implicit::update_state_incrementally(
    std::shared_ptr<const Core::LinAlg::Vector<double>> disiterinc)
{
  if (disiterinc == nullptr) return;

  check_init_setup();
  throw_if_state_not_in_sync_with_nox_group();
  ::NOX::Abstract::Group& grp = nln_solver().get_solution_group();

  auto* grp_ptr = dynamic_cast<NOX::Nln::Group*>(&grp);
  FOUR_C_ASSERT(grp_ptr != nullptr, "Dynamic cast failed!");

  // cast away const-qualifier for building the Nox Vector
  std::shared_ptr<Core::LinAlg::Vector<double>> mutable_disiterinc =
      Core::Utils::shared_ptr_from_ref(
          *const_cast<Core::LinAlg::Vector<double>*>(disiterinc.get()));

  // wrap the displacement vector in a NOX::Nln::Vector
  const NOX::Nln::Vector nox_disiterinc_ptr(mutable_disiterinc, NOX::Nln::Vector::MemoryType::View);

  // updated the state vector in the nox group
  grp_ptr->computeX(*grp_ptr, nox_disiterinc_ptr, 1.0);

  // Reset the state variables
  const auto& x_nox = dynamic_cast<const NOX::Nln::Vector&>(grp_ptr->getX());
  // set the consistent state in the models (e.g. structure and contact models)
  impl_int().reset_model_states(Core::LinAlg::Vector<double>(x_nox.get_linalg_vector()));
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Solid::TimeInt::Implicit::determine_stress_strain() { impl_int().determine_stress_strain(); }

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Solid::TimeInt::Implicit::evaluate(
    std::shared_ptr<const Core::LinAlg::Vector<double>> disiterinc)
{
  update_state_incrementally(disiterinc);

  evaluate();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Solid::TimeInt::Implicit::evaluate()
{
  check_init_setup();
  throw_if_state_not_in_sync_with_nox_group();
  ::NOX::Abstract::Group& grp = nln_solver().get_solution_group();

  auto* grp_ptr = dynamic_cast<NOX::Nln::Group*>(&grp);
  FOUR_C_ASSERT(grp_ptr != nullptr, "Dynamic cast failed!");

  // you definitely have to evaluate here. You might be called from a coupled
  // problem and the group might not be aware, that a different state than
  // the internally stored displacements may have changed.
  // This is a hack to get NOX to set IsValid to false.
  grp_ptr->setX(grp_ptr->getX());

  // compute the rhs vector and the stiffness matrix
  grp_ptr->compute_f_and_jacobian();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const ::NOX::Abstract::Group& Solid::TimeInt::Implicit::get_solution_group() const
{
  check_init_setup();
  return nln_solver().get_solution_group();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Inpar::Solid::ConvergenceStatus Solid::TimeInt::Implicit::perform_error_action(
    Inpar::Solid::ConvergenceStatus nonlinsoldiv)
{
  check_init_setup();
  retry_step_reason_ = Adapter::RetryStepReason::none;

  if (Core::Mat::TimeStepReduction::consume_mpi_reduced_request())
  {
    retry_step_reason_ = Adapter::RetryStepReason::material_time_step_reduction;
    return Inpar::Solid::conv_fail_repeat;
  }

  if (nonlinsoldiv == Inpar::Solid::conv_success)
  {
    check_for_material_time_step_increase();
    // Only relevant, if the input parameter DIVERCONT is used and set to divcontype_ == adapt_step:
    // In this case, the time step size is halved as consequence of a non-converging nonlinear
    // solver. After a prescribed number of converged time steps, the time step is doubled again.
    // The following methods checks, if the time step size can be increased again.
    check_for_time_step_increase(nonlinsoldiv);
    return Inpar::Solid::conv_success;
  }
  // get ID of actual processor in parallel
  const int& myrank = data_global_state().get_my_rank();

  // what to do when nonlinear solver does not converge
  switch (get_divergence_action())
  {
    case Inpar::Solid::divcont_stop:
    {
      // write restart output of last converged step before stopping
      output(true);

      // we should not get here, FOUR_C_THROW for safety
      FOUR_C_THROW("Nonlinear solver did not converge! ");
      return Inpar::Solid::conv_nonlin_fail;
      break;
    }
    case Inpar::Solid::divcont_continue:
    {
      if (myrank == 0)
      {
        Core::IO::cout
            << "\n WARNING: You are continuing your simulation although the nonlinear solver\n"
               " did not converge in the current time step.\n"
            << Core::IO::endl;
      }
      return Inpar::Solid::conv_success;
      break;
    }
    case Inpar::Solid::divcont_repeat_step:
    {
      FOUR_C_THROW(
          "DIVERCONT = repeat_step is no longer supported for structure_new implicit time "
          "integration. Use halve_step or adapt_step instead.");
    }
    case Inpar::Solid::divcont_halve_step:
    {
      if (myrank == 0)
      {
        Core::IO::cout << "Nonlinear solver failed to converge at time t= " << get_time_np()
                       << ". Divide timestep in half. " << "Old time step: " << get_delta_time()
                       << Core::IO::endl
                       << "New time step: " << 0.5 * get_delta_time() << Core::IO::endl
                       << Core::IO::endl;
      }
      retry_step_reason_ = Adapter::RetryStepReason::nonlinear_solver_time_step_reduction;
      return Inpar::Solid::conv_fail_repeat;
    }
    case Inpar::Solid::divcont_adapt_step:
    {
      if (myrank == 0)
      {
        Core::IO::cout << "Nonlinear solver failed to converge at time t= " << get_time_np()
                       << ". Divide timestep in half. " << "Old time step: " << get_delta_time()
                       << Core::IO::endl
                       << "New time step: " << 0.5 * get_delta_time() << Core::IO::endl
                       << Core::IO::endl;
      }

      set_div_con_refine_level(get_div_con_refine_level() + 1);
      set_div_con_num_fine_step(0);

      if (get_div_con_refine_level() == get_max_div_con_refine_level())
        FOUR_C_THROW(
            "Maximal divercont refinement level reached. Adapt your time basic time step size!");
      retry_step_reason_ = Adapter::RetryStepReason::nonlinear_solver_time_step_reduction;
      return Inpar::Solid::conv_fail_repeat;
    }
    case Inpar::Solid::divcont_rand_adapt_step:
    case Inpar::Solid::divcont_rand_adapt_step_ele_err:
    {
      FOUR_C_THROW(
          "DIVERCONT = rand_adapt_step and rand_adapt_step_ele_err are no longer supported for "
          "structure_new implicit time integration. Use halve_step or adapt_step instead.");
    }
    case Inpar::Solid::divcont_adapt_penaltycontact:
    {
      // adapt penalty and search parameter
      FOUR_C_THROW("Not yet implemented for new structure time integration");
      break;
    }
    case Inpar::Solid::divcont_repeat_simulation:
    {
      if (nonlinsoldiv == Inpar::Solid::conv_nonlin_fail and myrank == 0)
      {
        Core::IO::cout << "Nonlinear solver failed to converge and DIVERCONT = "
                          "repeat_simulation, hence leaving structural time integration "
                       << Core::IO::endl;
      }
      else if (nonlinsoldiv == Inpar::Solid::conv_lin_fail and myrank == 0)
      {
        Core::IO::cout << "Linear solver failed to converge and DIVERCONT = "
                          "repeat_simulation, hence leaving structural time integration "
                       << Core::IO::endl;
      }
      else if (nonlinsoldiv == Inpar::Solid::conv_ele_fail and myrank == 0)
      {
        Core::IO::cout
            << "Element failure in form of a negative Jacobian determinant and DIVERCONT = "
               "repeat_simulation, hence leaving structural time integration "
            << Core::IO::endl;
      }
      return nonlinsoldiv;  // so that time loop will be aborted
      break;
    }
    default:
      FOUR_C_THROW("Unknown DIVER_CONT case");
      return Inpar::Solid::conv_nonlin_fail;
      break;
  }
  return Inpar::Solid::conv_success;  // make compiler happy
}  // PerformErrorAction()

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Adapter::StepControlResult Solid::TimeInt::Implicit::post_solve_control_result(
    Inpar::Solid::ConvergenceStatus convergencestatus)
{
  Adapter::StepControlResult result;
  result.convergence_status = perform_error_action(convergencestatus);
  result.retry_reason = consume_retry_step_reason();

  if (result.convergence_status == Inpar::Solid::conv_success)
  {
    result.action = Adapter::StepControlAction::proceed;
  }
  else if (result.convergence_status == Inpar::Solid::conv_fail_repeat)
  {
    result.action = Adapter::StepControlAction::repeat_step;
    if (result.retry_reason == Adapter::RetryStepReason::material_time_step_reduction)
    {
      result.proposed_time_step = compute_material_reduced_time_step();
    }
    else if (result.retry_reason == Adapter::RetryStepReason::nonlinear_solver_time_step_reduction)
    {
      result.proposed_time_step = compute_nonlinear_retry_time_step(0.5);
    }
  }
  else
  {
    result.action = Adapter::StepControlAction::stop;
  }

  return result;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Adapter::RetryStepReason Solid::TimeInt::Implicit::consume_retry_step_reason()
{
  check_init_setup();

  const Adapter::RetryStepReason reason = retry_step_reason_;
  retry_step_reason_ = Adapter::RetryStepReason::none;
  return reason;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
std::optional<Inpar::Solid::MaterialTimeStepReductionSettings>
Solid::TimeInt::Implicit::material_time_step_reduction_settings() const
{
  check_init_setup();

  if (not get_data_sdyn().material_time_step_reduction_enabled()) return std::nullopt;
  return get_data_sdyn().material_time_step_reduction_settings();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Solid::TimeInt::Implicit::reset_step_for_time_step_retry(const double dtnew)
{
  check_init_setup();

  material_time_step_is_reduced_ = true;
  material_time_step_reduction_num_fine_step_ = 0;
  apply_time_step_change(dtnew, true);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Solid::TimeInt::Implicit::apply_time_step_for_next_step(const double dtnew)
{
  check_init_setup();
  apply_time_step_change(dtnew, false);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool Solid::TimeInt::Implicit::apply_scheduled_time_step_increase()
{
  check_init_setup();

  if (not pending_time_step_increase_dt_.has_value()) return false;

  const int myrank = data_global_state().get_my_rank();
  const double old_dt = get_delta_time();
  double new_dt = *pending_time_step_increase_dt_;
  const Adapter::TimeStepIncreaseReason reason = pending_time_step_increase_reason_;

  pending_time_step_increase_dt_.reset();
  pending_time_step_increase_reason_ = Adapter::TimeStepIncreaseReason::none;

  const double remaining_time = get_time_end() - get_time_n();
  FOUR_C_ASSERT_ALWAYS(remaining_time >= 0.0,
      "Remaining time to end time became negative (remaining = {}).", remaining_time);
  new_dt = std::min(new_dt, remaining_time);

  if (not(new_dt > 0.0 and old_dt > 0.0 and new_dt != old_dt)) return false;

  if (myrank == 0)
  {
    switch (reason)
    {
      case Adapter::TimeStepIncreaseReason::material_recovery:
      {
        const double applied_factor = new_dt / old_dt;
        Core::IO::cout << std::string(72, '*') << "\n"
                       << "Material time step reduction recovery: increase timestep size by factor "
                       << applied_factor << "." << Core::IO::endl
                       << "Old time step: " << old_dt << Core::IO::endl
                       << "New time step: " << new_dt << Core::IO::endl
                       << std::string(72, '*') << "\n";
        break;
      }
      case Adapter::TimeStepIncreaseReason::divergence_control_recovery:
      {
        Core::IO::cout << "Nonlinear solver successful. Double timestep size!" << Core::IO::endl;
        break;
      }
      case Adapter::TimeStepIncreaseReason::none:
      default:
        break;
    }
  }

  apply_time_step_change(new_dt, false);

  if (reason == Adapter::TimeStepIncreaseReason::material_recovery and
      new_dt >= data_sdyn().material_time_step_reduction_settings().max_timestep.value())
  {
    material_time_step_is_reduced_ = false;
    material_time_step_reduction_num_fine_step_ = 0;
  }

  return true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
std::optional<double> Solid::TimeInt::Implicit::compute_material_reduced_time_step() const
{
  check_init_setup();

  if (not get_data_sdyn().material_time_step_reduction_enabled()) return std::nullopt;

  const auto& settings = get_data_sdyn().material_time_step_reduction_settings();
  const double new_dt = get_delta_time() * settings.factor;
  if (new_dt < settings.min_timestep)
  {
    FOUR_C_THROW(
        "A material evaluation requested a reduced time step at time t = {}, but reducing dt "
        "would violate MATERIAL TIMESTEP REDUCTION/MIN_TIMESTEP. Current dt = {}, requested dt = "
        "{}, MIN_TIMESTEP = {}.",
        get_time_np(), get_delta_time(), new_dt, settings.min_timestep);
  }

  return new_dt;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Solid::TimeInt::Implicit::schedule_time_step_increase(
    const double new_dt, const Adapter::TimeStepIncreaseReason reason)
{
  check_init_setup();

  FOUR_C_ASSERT_ALWAYS(new_dt > 0.0, "Scheduled time-step increase must be positive.");

  pending_time_step_increase_dt_ = new_dt;
  pending_time_step_increase_reason_ = reason;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Solid::TimeInt::Implicit::check_for_material_time_step_increase()
{
  check_init_setup();

  if (not material_time_step_is_reduced_ or not data_sdyn().material_time_step_reduction_enabled())
    return;

  const auto& material_ts_reduction = data_sdyn().material_time_step_reduction_settings();
  const double max_timestep = material_ts_reduction.max_timestep.value();
  const double old_dt = get_delta_time();

  if (old_dt >= max_timestep)
  {
    material_time_step_is_reduced_ = false;
    material_time_step_reduction_num_fine_step_ = 0;
    return;
  }

  ++material_time_step_reduction_num_fine_step_;
  if (material_time_step_reduction_num_fine_step_ < material_ts_reduction.steps_to_increase) return;

  const double increase_factor = material_ts_reduction.increase_factor.value();
  const double candidate_dt = old_dt * increase_factor;
  const double new_dt = std::min(candidate_dt, max_timestep);

  if (new_dt > old_dt)
  {
    schedule_time_step_increase(new_dt, Adapter::TimeStepIncreaseReason::material_recovery);
  }
  material_time_step_reduction_num_fine_step_ = 0;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Solid::TimeInt::Implicit::apply_time_step_change(
    const double new_dt, const bool reset_step_for_retry)
{
  check_init_setup();

  FOUR_C_ASSERT_ALWAYS(new_dt > 0.0, "Time-step size must be positive.");

  set_delta_time(new_dt);
  set_time_np(get_time_n() + get_delta_time());

  if (reset_step_for_retry) reset_step();
  integrator().update_constant_state_contributions();
}

/*-----------------------------------------------------------------------------*
 * check, if according to divercont flag                             meier 01/15
 * time step size can be increased
 *-----------------------------------------------------------------------------*/
void Solid::TimeInt::Implicit::check_for_time_step_increase(Inpar::Solid::ConvergenceStatus status)
{
  check_init_setup();

  const int maxnumfinestep = 4;

  if (get_divergence_action() != Inpar::Solid::divcont_adapt_step)
    return;
  else if (status == Inpar::Solid::conv_success and get_div_con_refine_level() != 0)
  {
    set_div_con_num_fine_step(get_div_con_num_fine_step() + 1);

    if (get_div_con_num_fine_step() == maxnumfinestep)
    {
      set_div_con_refine_level(get_div_con_refine_level() - 1);
      set_div_con_num_fine_step(0);
      schedule_time_step_increase(
          get_delta_time() * 2.0, Adapter::TimeStepIncreaseReason::divergence_control_recovery);
    }
    return;
  }
}  // check_for_time_step_increase()

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double Solid::TimeInt::Implicit::compute_nonlinear_retry_time_step(
    const double scaling_factor) const
{
  check_init_setup();

  FOUR_C_ASSERT_ALWAYS(scaling_factor > 0.0, "Time-step scaling factor must be positive.");
  return get_delta_time() * scaling_factor;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Solid::TimeInt::Implicit::print_jacobian_in_matlab_format(
    const NOX::Nln::Group& curr_grp) const
{
  using LinalgBlockSparseMatrix =
      Core::LinAlg::BlockSparseMatrix<Core::LinAlg::DefaultBlockMatrixStrategy>;

  if (not get_data_io().is_write_jacobian_to_matlab()) return;

  // create file name
  std::stringstream filebase;

  filebase << "str_jacobian" << "_step-" << get_data_global_state().get_step_np() << "_nlniter-"
           << nlnsolver_ptr_->get_num_nln_iterations();

  std::stringstream filename;
  filename << get_data_io().get_output_ptr()->output().file_name() << "_" << filebase.str()
           << ".mtl";

  if (get_data_global_state().get_my_rank() == 0)
    std::cout << "Writing structural jacobian to \"" << filename.str() << "\"\n";

  Teuchos::RCP<const NOX::Nln::LinearSystemBase> linear_system = curr_grp.get_linear_system();

  Teuchos::RCP<const NOX::Nln::LinearSystem> nln_lin_system =
      Teuchos::rcp_dynamic_cast<const NOX::Nln::LinearSystem>(linear_system, true);

  const NOX::Nln::LinSystem::OperatorType jac_type = nln_lin_system->get_jacobian_operator_type();

  auto jac_ptr = nln_lin_system->get_jacobian_operator();

  switch (jac_type)
  {
    case NOX::Nln::LinSystem::LinalgSparseMatrix:
    {
      auto sparse_matrix = std::dynamic_pointer_cast<const Core::LinAlg::SparseMatrix>(jac_ptr);
      Core::LinAlg::print_matrix_in_matlab_format(filename.str().c_str(), *sparse_matrix);

      break;
    }
    case NOX::Nln::LinSystem::LinalgBlockSparseMatrix:
    {
      auto block_matrix = std::dynamic_pointer_cast<const LinalgBlockSparseMatrix>(jac_ptr);
      Core::LinAlg::print_block_matrix_in_matlab_format(filename.str(), *block_matrix);

      break;
    }
    default:
    {
      FOUR_C_THROW("Unsupported NOX::Nln::LinSystem::OperatorType: \"{}\"",
          NOX::Nln::LinSystem::operator_type_to_string(jac_type));
    }
  }
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Inpar::Solid::DynamicType Solid::TimeInt::Implicit::method_name() const
{
  return implint_ptr_->method_name();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
int Solid::TimeInt::Implicit::method_steps() const { return implint_ptr_->method_steps(); }

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
int Solid::TimeInt::Implicit::method_order_of_accuracy_dis() const
{
  return implint_ptr_->method_order_of_accuracy_dis();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
int Solid::TimeInt::Implicit::method_order_of_accuracy_vel() const
{
  return implint_ptr_->method_order_of_accuracy_vel();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double Solid::TimeInt::Implicit::method_lin_err_coeff_dis() const
{
  return implint_ptr_->method_lin_err_coeff_dis();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double Solid::TimeInt::Implicit::method_lin_err_coeff_vel() const
{
  return implint_ptr_->method_lin_err_coeff_vel();
}

FOUR_C_NAMESPACE_CLOSE
