// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_material_time_step_request.hpp"

#include "4C_comm_mpi_utils.hpp"
#include "4C_utils_exceptions.hpp"

#include <exception>
#include <iostream>
#include <sstream>

FOUR_C_NAMESPACE_OPEN

namespace
{
  enum class Request
  {
    none,
    reduce_step,
    /// A real request was consumed. This is distinct from none so scope destructors can tell that
    /// an outstanding request was handled without treating "no request happened" as an error.
    consumed
  };

  /**
   * @brief Thread-local state for the material timestep-reduction request channel.
   *
   * The state object is intentionally private to this translation unit. The public API exposes the
   * operations that evaluator, NOX, and time-loop infrastructure need, but does not expose storage
   * mechanics. Material code stays on the small Core::Mat::request_time_step_reduction(reason) API.
   */
  struct RequestState
  {
    /// Did the current rank request time step reduction?
    Request local_request = Request::none;
    /// Did any rank request timestep reduction in the last MPI reduction?
    Request mpi_reduced_request = Request::none;
    /// The current algorithm supports material-triggered timestep reduction because it has a retry
    /// policy for it.
    bool support_scope_active = false;
    /// The current evaluation supports material-triggered timestep reduction because it will
    /// consume and MPI-reduce requests before returning.
    bool evaluation_scope_active = false;

    /// Clear both local and MPI-reduced request state without affecting active scope contracts.
    void reset_all_requests()
    {
      local_request = Request::none;
      mpi_reduced_request = Request::none;
    }

    /// Clear the rank-local request state for a new evaluator-owned material evaluation.
    void reset_local_request() { local_request = Request::none; }

    bool consume_local_request()
    {
      switch (local_request)
      {
        case Request::reduce_step:
          local_request = Request::consumed;
          return true;
        case Request::none:
          local_request = Request::consumed;
          return false;
        case Request::consumed:
          FOUR_C_THROW(
              "Attempted to consume the local material timestep-reduction request more than once. "
              "The request channel is one-shot and each request may be consumed exactly once.");
        default:
          FOUR_C_THROW("Unknown local request state while consuming.");
      }
    }
  };
  // One RequestState per thread; all request-channel operations on this thread share it.
  thread_local RequestState material_time_step_request_state;

}  // namespace

void Core::Mat::request_time_step_reduction(const std::string& reason)
{
  // A single evaluator scope can contain many material points. Once one material point has
  // requested timestep reduction, the rank-local evaluation outcome is already fixed; repeated
  // requests before the evaluator consumes the flag would only spam the same diagnostic during
  // retries or across several failing points. The first request still prints immediately at the
  // material call site, including unsupported-scope cases that throw below.
  const bool is_first_pending_request =
      (material_time_step_request_state.local_request != Request::reduce_step);
  if (is_first_pending_request)
  {
    std::ostringstream message;
    message << std::string(72, '*') << "\n"
            << "Material evaluation failed with a time step reduction request.\n"
            << "Reason: " << reason << "\n"
            << std::string(72, '*') << "\n";
    std::cout << message.str() << std::flush;
  }

  // Ensure that both algorithm and evaluation support timestep reduction. Only if they both support
  // it we don't throw.

  if (not material_time_step_request_state.evaluation_scope_active)
  {
    FOUR_C_THROW(
        "A material requested timestep reduction in an evaluation path that does not "
        "support it. Install Core::Mat::TimeStepReduction::EvaluationScope around the synchronized "
        "element/material evaluation and collect requests before scope exit.");
  }

  if (not material_time_step_request_state.support_scope_active)
  {
    FOUR_C_THROW(
        "A material requested timestep reduction, but the calling algorithm has not declared "
        "support for material-triggered timestep retry. Install "
        "Core::Mat::TimeStepReduction::AlgorithmSupportScope only around algorithms that reduce "
        "dt, restore state, and repeat the global timestep.");
  }

  // both scopes are active, set the local request flag for this rank and return silently.
  material_time_step_request_state.local_request = Request::reduce_step;
}

bool Core::Mat::TimeStepReduction::collect_local_requests(const MPI_Comm comm)
{
  if (not material_time_step_request_state.evaluation_scope_active)
  {
    FOUR_C_THROW(
        "Attempted to collect material time step reduction requests, but no evaluation scope is "
        "active. Install Core::Mat::TimeStepReduction::EvaluationScope around the synchronized "
        "element/material evaluation and collect requests before scope exit.");
  }

  if (material_time_step_request_state.mpi_reduced_request == Request::reduce_step)
  {
    FOUR_C_THROW(
        "Attempted to publish a new MPI-reduced material timestep-reduction request while a "
        "previous reduced request is still pending. The timestep policy must consume the pending "
        "request before the next reduced request is published.");
  }

  const int local_request = material_time_step_request_state.consume_local_request() ? 1 : 0;
  const bool has_global_request = (Core::Communication::max_all(local_request, comm) != 0);
  if (has_global_request)
    material_time_step_request_state.mpi_reduced_request = Request::reduce_step;
  else
    material_time_step_request_state.mpi_reduced_request = Request::none;
  return has_global_request;
}

bool Core::Mat::TimeStepReduction::has_mpi_reduced_request()
{
  switch (material_time_step_request_state.mpi_reduced_request)
  {
    case Request::none:
      return false;
    case Request::reduce_step:
      return true;
    case Request::consumed:
      FOUR_C_THROW(
          "Attempted to query an MPI-reduced material timestep-reduction request after it was "
          "already consumed. The request channel is one-shot; publish a new request before "
          "querying "
          "again.");
  }

  FOUR_C_THROW("Unknown MPI-reduced request state while querying.");
}

bool Core::Mat::TimeStepReduction::consume_mpi_reduced_request()
{
  switch (material_time_step_request_state.mpi_reduced_request)
  {
    case Request::none:
      material_time_step_request_state.mpi_reduced_request = Request::consumed;
      return false;
    case Request::reduce_step:
      material_time_step_request_state.mpi_reduced_request = Request::consumed;
      return true;
    case Request::consumed:
      FOUR_C_THROW(
          "Attempted to consume an MPI-reduced material timestep-reduction request more than once. "
          "The request channel is one-shot and each request may be consumed exactly once.");
  }

  FOUR_C_THROW("Unknown MPI-reduced request state while consuming.");
}

Core::Mat::TimeStepReduction::AlgorithmSupportScope::AlgorithmSupportScope(const bool enabled)
    : enabled_(enabled), uncaught_exceptions_on_entry_(std::uncaught_exceptions())
{
  if (enabled_)
  {
    // no nested scopes!
    if (material_time_step_request_state.support_scope_active)
    {
      FOUR_C_THROW(
          "Nested material timestep-reduction support scopes are not supported. The algorithm that "
          "owns timestep retry must install exactly one support scope.");
    }

    material_time_step_request_state.reset_all_requests();
    material_time_step_request_state.support_scope_active = true;
  }
}

Core::Mat::TimeStepReduction::AlgorithmSupportScope::~AlgorithmSupportScope() noexcept(false)
{
  if (not enabled_) return;  // inert scope, nothing to do

  if (not material_time_step_request_state.support_scope_active)
  {
    FOUR_C_THROW(
        "Material timestep-reduction support scope ended although no support scope was active.");
  }

  material_time_step_request_state.support_scope_active = false;

  if (std::uncaught_exceptions() > uncaught_exceptions_on_entry_)
  {
    // Preserve the original exception from the owning algorithm. A pending MPI-reduced request
    // cannot be resolved safely while another exception is already leaving the supported region, so
    // clear request state and let the original failure continue upward.
    material_time_step_request_state.reset_all_requests();
    return;
  }

  if (material_time_step_request_state.mpi_reduced_request == Request::reduce_step)
  {
    FOUR_C_THROW(
        "You opened a TimeStepReduction::AlgorithmSupportScope but failed to reduce the time step "
        "reduction request. If your algorithm supports material-triggered timestep reduction, you "
        "must check for requests after each evaluation and reduce dt accordingly before the next "
        "evaluation.");
  }
}

Core::Mat::TimeStepReduction::EvaluationScope::EvaluationScope()
    : uncaught_exceptions_on_entry_(std::uncaught_exceptions())
{
  if (material_time_step_request_state.evaluation_scope_active)
  {
    FOUR_C_THROW(
        "Nested material timestep-reduction evaluation scopes are not supported. The evaluator "
        "that owns the synchronized error check must install exactly one evaluation scope.");
  }

  material_time_step_request_state.reset_local_request();
  material_time_step_request_state.evaluation_scope_active = true;
}

Core::Mat::TimeStepReduction::EvaluationScope::~EvaluationScope() noexcept(false)
{
  if (not material_time_step_request_state.evaluation_scope_active)
  {
    FOUR_C_THROW(
        "Material timestep-reduction evaluation scope ended although no evaluation scope was "
        "active.");
  }

  material_time_step_request_state.evaluation_scope_active = false;

  if (std::uncaught_exceptions() > uncaught_exceptions_on_entry_)
  {
    // Preserve the original exception from the evaluator/material stack. Throwing a second
    // exception from this destructor while unwinding would either mask the real error or terminate.
    // The local request cannot be consumed and MPI-reduced on this path, so clear it before the
    // original exception continues upward.
    material_time_step_request_state.reset_local_request();
    return;
  }

  if (material_time_step_request_state.local_request != Request::consumed)
  {
    FOUR_C_THROW(
        "You opened a TimeStepReduction::EvaluationScope but failed to call "
        "TimeStepReduction::collect_local_requests().");
  }
}

FOUR_C_NAMESPACE_CLOSE
