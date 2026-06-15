// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_STRUCTURE_NEW_TIMINT_IMPLICIT_HPP
#define FOUR_C_STRUCTURE_NEW_TIMINT_IMPLICIT_HPP

#include "4C_config.hpp"

#include "4C_structure_new_timint_implicitbase.hpp"  // base class

#include <optional>

FOUR_C_NAMESPACE_OPEN

// forward declarations ...
namespace Solid
{
  namespace IMPLICIT
  {
    class Generic;
  }  // namespace IMPLICIT
  namespace Predict
  {
    class Generic;
  }  // namespace Predict
  namespace Nln::SOLVER
  {
    class Generic;
    namespace INTERFACE
    {
      class Required;
    }  // namespace INTERFACE
  }  // namespace Nln::SOLVER

  namespace TimeInt
  {
    /** \brief Implicit time integration strategy
     *
     * */
    class Implicit : public ImplicitBase
    {
     public:
      void setup() override;

      int integrate() override;

      int integrate_step() override;

      /// set the state of the nox group and the global state data container
      /// see class \ref Adapter::StructureNew for a detailed documentation.
      void set_state(const std::shared_ptr<Core::LinAlg::Vector<double>>& x) override;

      /*! \brief nonlinear solve
       *
       *  Do the nonlinear solve, i.e. (multiple) corrector,
       *  for the time step. All boundary conditions have
       *  been set. */
      Inpar::Solid::ConvergenceStatus solve() override;

      /** \brief Identify residual
       *
       *  This method does not predict the target solution but
       *  evaluates the residual and the stiffness matrix.
       *  In partitioned solution schemes, it is better to keep the current
       *  solution instead of evaluating the initial guess (as the predictor)
       *  does. */
      void prepare_partition_step() override;

      //! Prepare time step
      void prepare_time_step() override;

      //! Prepare the next step and report whether the outer loop should proceed or repeat.
      Adapter::PrepareTimeStepStatus prepare_time_step_with_status() override;

      //! Return the adapter-level prepare-step decision, including retry metadata when needed.
      Adapter::StepControlResult prepare_time_step_control_result() override;

      //! @name Accessors
      //! @{
      //! return the predictor
      [[nodiscard]] const Solid::Predict::Generic& predictor() const
      {
        check_init_setup();
        return *predictor_ptr_;
      }

      //! @}
      [[nodiscard]] std::shared_ptr<const Solid::Nln::SOLVER::Generic> get_nln_solver_ptr() const
      {
        return nlnsolver_ptr_;
      };

      //! do something in case nonlinear solution does not converge for some reason
      Inpar::Solid::ConvergenceStatus perform_error_action(
          Inpar::Solid::ConvergenceStatus nonlinsoldiv) override;

      //! Return the adapter-level post-solve decision, including retries and deferred dt changes.
      Adapter::StepControlResult post_solve_control_result(
          Inpar::Solid::ConvergenceStatus convergencestatus) override;

      //! The standalone structure time loop owns material timestep-reduction retry policy.
      [[nodiscard]] bool supports_material_time_step_reduction() const override { return true; }

      //! Return and clear the pending retry reason raised during prepare/solve handling.
      Adapter::RetryStepReason consume_retry_step_reason() override;

      //! Expose configured material timestep-reduction settings to adapter-level loop code.
      [[nodiscard]] std::optional<Inpar::Solid::MaterialTimeStepReductionSettings>
      material_time_step_reduction_settings() const override;

      //! Reduce dt and reset state so the current step can be retried consistently.
      void reset_step_for_time_step_retry(double dtnew) override;

      //! Apply a deferred dt change for the next step without retrying the current one.
      void apply_time_step_for_next_step(double dtnew) override;

      //! Consume and apply a previously scheduled next-step dt increase, if any.
      bool apply_scheduled_time_step_increase() override;

      //! check, if according to divercont flag time step size can be increased
      void check_for_time_step_increase(Inpar::Solid::ConvergenceStatus status);

      //! returns pointer to generic implicit object
      std::shared_ptr<Solid::IMPLICIT::Generic> impl_int_ptr()
      {
        check_init_setup();
        return implint_ptr_;
      };

      /// Update State Incrementally for coupled problems with monolithic approach
      void update_state_incrementally(
          std::shared_ptr<const Core::LinAlg::Vector<double>> disiterinc) override;

      void determine_stress_strain() override;

      ///  Evaluate routine for coupled problems with monolithic approach
      void evaluate(std::shared_ptr<const Core::LinAlg::Vector<double>> disiterinc) override;
      void evaluate() override;

      /** \brief Print structural jacobian matrix into a text file for later use
       *  in MATLAB
       *
       *  This routine can be activated via the input parameter
       *  %STRUCT_JACOBIAN_MATLAB. See the corresponding inpar section for more
       *  details.
       *
       *  The text file can be found in the user-provided output folder using the
       *  following file name extension:
       *
       *  [OUTPUT-FOLDER]/[OUTPUT FILE NAME]_str_jacobian_step-[STEP]_nlniter-[NEWTON-ITERATION].mtl
       *
       *  */
      void print_jacobian_in_matlab_format(const NOX::Nln::Group& curr_grp) const;

     protected:
      //! returns the current solution group
      [[nodiscard]] const ::NOX::Abstract::Group& get_solution_group() const override;

      Solid::IMPLICIT::Generic& impl_int()
      {
        check_init_setup();
        return *implint_ptr_;
      };

      Solid::Predict::Generic& predictor()
      {
        check_init_setup();
        return *predictor_ptr_;
      };

      std::shared_ptr<Solid::Predict::Generic> predictor_ptr()
      {
        check_init_setup();
        return predictor_ptr_;
      };

      [[nodiscard]] const Solid::Nln::SOLVER::Generic& nln_solver() const
      {
        check_init_setup();
        return *nlnsolver_ptr_;
      };

      Solid::Nln::SOLVER::Generic& nln_solver()
      {
        check_init_setup();
        return *nlnsolver_ptr_;
      };

      std::shared_ptr<Solid::Nln::SOLVER::Generic> nln_solver_ptr()
      {
        check_init_setup();
        return nlnsolver_ptr_;
      };

      //! @name Attribute access functions
      //@{

      //! Provide Name
      Inpar::Solid::DynamicType method_name() const override;

      //! Provide number of steps, e.g. a single-step method returns 1,
      //! a \f$m\f$-multistep method returns \f$m\f$
      int method_steps() const override;

      //! Give local order of accuracy of displacement part
      int method_order_of_accuracy_dis() const override;

      //! Give local order of accuracy of velocity part
      int method_order_of_accuracy_vel() const override;

      //! Return linear error coefficient of displacements
      double method_lin_err_coeff_dis() const override;

      //! Return linear error coefficient of velocities
      double method_lin_err_coeff_vel() const override;

      ///@}

     private:
      //! Reject divergence-control modes that are not supported by structure_new implicit.
      void validate_divergence_action_configuration() const;

      //! Check whether repeated successful fine steps allow recovery of a previously reduced dt.
      void check_for_material_time_step_increase();

      //! Compute the reduced retry dt requested by material evaluation, including safety checks.
      [[nodiscard]] std::optional<double> compute_material_reduced_time_step() const;

      //! Compute the reduced retry dt requested by nonlinear divergence control.
      [[nodiscard]] double compute_nonlinear_retry_time_step(double scaling_factor) const;

      //! Store a dt increase to be applied at the next loop boundary after a successful step.
      void schedule_time_step_increase(double new_dt, Adapter::TimeStepIncreaseReason reason);

      //! apply a new time-step size locally and optionally reset the current step for retry
      void apply_time_step_change(double new_dt, bool reset_step_for_retry);

      //! ptr to the implicit time integrator object
      std::shared_ptr<Solid::IMPLICIT::Generic> implint_ptr_;

      //! ptr to the predictor object
      std::shared_ptr<Solid::Predict::Generic> predictor_ptr_;

      //! ptr to the non-linear solver object
      std::shared_ptr<Solid::Nln::SOLVER::Generic> nlnsolver_ptr_;

      //! ptr to the nox group object
      std::shared_ptr<::NOX::Abstract::Group> grp_ptr_;

      Adapter::RetryStepReason retry_step_reason_ = Adapter::RetryStepReason::none;
      int material_time_step_reduction_num_fine_step_ = 0;
      std::optional<double> pending_time_step_increase_dt_;
      Adapter::TimeStepIncreaseReason pending_time_step_increase_reason_ =
          Adapter::TimeStepIncreaseReason::none;
      bool material_time_step_is_reduced_ = false;
    };
  }  // namespace TimeInt
}  // namespace Solid

FOUR_C_NAMESPACE_CLOSE

#endif
