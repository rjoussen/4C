// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_MATERIAL_TIME_STEP_REQUEST_HPP
#define FOUR_C_MATERIAL_TIME_STEP_REQUEST_HPP

#include "4C_config.hpp"

#include <mpi.h>

#include <string>

FOUR_C_NAMESPACE_OPEN

/*!
 * \namespace Core::Mat::TimeStepReduction
 * \brief Material-triggered time-step-reduction request channel.
 *
 * Material time-step-reduction requests move in two steps. A material first records a rank-local
 * request with request_time_step_reduction(). The evaluator then calls collect_local_requests() to
 * reduce these rank-local requests over MPI and publish whether any rank requested a retry. The
 * time-step algorithm finally consumes the reduced request and reduces dt before retrying.
 *
 * Requests are only valid while both a MaterialEvaluationScope and a TimeStepRetryScope are active.
 * This keeps material requests limited to synchronized evaluations whose owning algorithm can
 * actually reduce the time-step.
 */
namespace Core::Mat::TimeStepReduction
{
  /*!
   * \brief Request a retry of the current time-step with a smaller time-step size.
   *
   * Materials call this when the current state cannot be evaluated reliably, but a smaller
   * time-step may make the next attempt valid. After making the request, the material should return
   * through the normal evaluation path instead of throwing its own exception.
   *
   * A request is only made if the current evaluation and algorithm support material-triggered
   * time-step retry. Otherwise, an exception is thrown immediately and includes \p reason in the
   * error message.
   *
   * @param reason Short diagnostic text explaining why the material requests a smaller time-step.
   */
  void request(const std::string& reason);

  /*!
   * \brief Collect rank-local material requests and publish the MPI-reduced result.
   *
   * Call this exactly once after a synchronized material/element evaluation while the matching
   * MaterialEvaluationScope is still active. The function MPI-reduces the rank-local requests,
   * clears the rank-local pending request, and stores whether any rank
   * requested time-step reduction.
   *
   * @return true if at least one rank requested time-step reduction.
   */
  bool collect_local_requests();

  /*!
   * \brief Return whether an MPI-reduced time-step-reduction request is pending.
   *
   * This does not consume the request. It is useful for code that needs to forward or inspect the
   * reduced decision before the time-step policy handles it.
   */
  bool peek_mpi_reduced_request();

  /*!
   * \brief Consume the pending MPI-reduced time-step-reduction request.
   *
   * Call this exactly once inside a TimeStepRetryScope at the time-step retry policy site to mark a
   * pending MPI-reduced request as consumed. Calling this without a pending request is allowed and
   * returns false.
   *
   * @return true if a pending request was consumed, false if no request was pending.
   */
  bool consume_mpi_reduced_request();

  /*!
   * \brief Mark an algorithm as owning a time-step retry mechanism.
   *
   * Install this scope once per time-step around the algorithm region that can reduce
   * the time-step size, restore state, and repeat the time-step consistently on all MPI ranks.
   * Usually, this is inside the time-step loop in the time integrator.
   *
   * Passing \p enabled=false makes the scope inert. This is useful when a time-loop is
   * shared by multiple algorithms of which only some have a time-step retry mechanism.
   *
   * Before the scope ends, any pending MPI-reduced requests must be consumed by
   * consume_mpi_reduced_request(). This ensures that the algorithm actually handles the request.
   */
  class TimeStepRetryScope
  {
   public:
    explicit TimeStepRetryScope(bool enabled = true);
    TimeStepRetryScope(const TimeStepRetryScope&) = delete;
    TimeStepRetryScope& operator=(const TimeStepRetryScope&) = delete;
    ~TimeStepRetryScope() noexcept(false);

   private:
    // whether this scope is enabled or not.
    const bool enabled_;
    // the number of uncaught exceptions when the scope was entered. Used to detect whether the
    // scope is being exited due to an exception thrown inside the scope.
    const int uncaught_exceptions_on_entry_ = 0;
  };

  /*!
   * \brief Mark one synchronized material evaluation as request-aware.
   *
   * Install this scope around one material/element evaluation in which all ranks will later call
   * collect_local_requests(). The communicator passed to the constructor is used for the MPI
   * reduction.
   *
   * If a material request is made in this scope, collect_local_requests() must be called before the
   * scope ends. This ensures that the request is collected and published to all ranks before the
   * scope ends.
   */
  class MaterialEvaluationScope
  {
   public:
    explicit MaterialEvaluationScope(MPI_Comm comm);
    MaterialEvaluationScope(const MaterialEvaluationScope&) = delete;
    MaterialEvaluationScope& operator=(const MaterialEvaluationScope&) = delete;
    ~MaterialEvaluationScope() noexcept(false);

   private:
    // the number of uncaught exceptions when the scope was entered. Used to detect whether the
    // scope is being exited due to an exception thrown inside the scope.
    const int uncaught_exceptions_on_entry_ = 0;
  };

}  // namespace Core::Mat::TimeStepReduction

FOUR_C_NAMESPACE_CLOSE

#endif
