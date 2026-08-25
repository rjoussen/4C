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
  /// Private storage container for the state of the material time-step-reduction request channel.
  struct RequestState
  {
    /// flag indicating that the current rank requested time-step reduction which has not yet been
    /// collected
    bool rank_request_pending = false;
    /// global (MPI-reduced) flag indicating that at least one rank requested time-step reduction
    /// which has not yet been consumed by the algorithm
    bool global_request_pending = false;
    /// flag indicating that the MPI-reduced request has been consumed
    bool global_request_consumed = false;
    /// flag indicating that the time-step retry scope is active
    bool time_step_retry_scope_active = false;
    /// flag indicating that the evaluation scope is active
    bool evaluation_scope_active = false;
    /// communicator used for the evaluation scope
    MPI_Comm evaluation_comm = MPI_COMM_NULL;
  };

  /// Every rank owns a single instance of the request state, private to this file.
  RequestState request_state;

}  // namespace

void Core::Mat::TimeStepReduction::request(const std::string& reason)
{
  FOUR_C_ASSERT_ALWAYS(request_state.evaluation_scope_active,
      "A material requested time-step reduction in an evaluation path that does not support it.\n"
      "Reason: {}\n"
      "Install Core::Mat::TimeStepReduction::MaterialEvaluationScope around the "
      "synchronized element/material evaluation and collect requests before scope exit.",
      reason);

  FOUR_C_ASSERT_ALWAYS(request_state.time_step_retry_scope_active,
      "A material requested time-step reduction but the calling algorithm does not support it.\n"
      "Reason: {}\n"
      "Install Core::Mat::TimeStepReduction::TimeStepRetryScope in algorithms that have a "
      "time-step retry mechanism and consume global time-step-reduction requests before scope "
      "exit.",
      reason);

  // only print the first request per rank
  if (not request_state.rank_request_pending)
  {
    auto message = std::ostringstream();
    message << std::string(72, '*') << "\n"
            << "Material evaluation failed with a time step reduction request on rank "
            << Core::Communication::my_mpi_rank(request_state.evaluation_comm) << ".\n"
            << "Reason: " << reason << "\n"
            << std::string(72, '*') << "\n";
    std::cout << message.str() << std::flush;
  }

  // register the rank-local request and leave silently
  request_state.rank_request_pending = true;
}

bool Core::Mat::TimeStepReduction::collect_local_requests()
{
  FOUR_C_ASSERT_ALWAYS(request_state.evaluation_scope_active,
      "Attempted to collect material time step reduction requests outside of an active evaluation "
      "scope.");

  FOUR_C_ASSERT_ALWAYS(not request_state.global_request_pending,
      "Attempted to MPI-reduce material time-step-reduction requests while a "
      "previous reduced request is still pending.");

  const int local_request = static_cast<int>(request_state.rank_request_pending);
  request_state.global_request_pending =
      static_cast<bool>(Core::Communication::max_all(local_request, request_state.evaluation_comm));

  // rank-local requests have been collected, they are no longer pending.
  request_state.rank_request_pending = false;

  return request_state.global_request_pending;
}

bool Core::Mat::TimeStepReduction::peek_mpi_reduced_request()
{
  return request_state.global_request_pending;
}

bool Core::Mat::TimeStepReduction::consume_mpi_reduced_request()
{
  FOUR_C_ASSERT_ALWAYS(request_state.time_step_retry_scope_active,
      "Attempted to consume an MPI-reduced material time-step-reduction request outside of an "
      "active time-step retry scope.");

  FOUR_C_ASSERT_ALWAYS(not request_state.global_request_consumed,
      "Attempted to consume the MPI-reduced material time-step-reduction request more than once in "
      "the same time-step retry scope.");

  const bool consumed_request = request_state.global_request_pending;
  request_state.global_request_consumed = true;
  request_state.global_request_pending = false;
  return consumed_request;
}

Core::Mat::TimeStepReduction::TimeStepRetryScope::TimeStepRetryScope(const bool enabled)
    : enabled_(enabled), uncaught_exceptions_on_entry_(std::uncaught_exceptions())
{
  if (not enabled_) return;

  FOUR_C_ASSERT_ALWAYS(not request_state.time_step_retry_scope_active,
      "Nested material time-step-reduction retry scopes are not supported. The algorithm that "
      "owns time-step retry must install exactly one retry scope.");

  request_state.time_step_retry_scope_active = true;
  request_state.global_request_consumed = false;
}

Core::Mat::TimeStepReduction::TimeStepRetryScope::~TimeStepRetryScope() noexcept(false)
{
  if (not enabled_) return;

  if (std::uncaught_exceptions() > uncaught_exceptions_on_entry_)
  {
    request_state.global_request_pending = false;
    request_state.global_request_consumed = false;
    request_state.time_step_retry_scope_active = false;
    return;
  }

  if (request_state.global_request_pending)
  {
    request_state.global_request_pending = false;
    request_state.global_request_consumed = false;
    request_state.time_step_retry_scope_active = false;
    FOUR_C_THROW(
        "You opened a TimeStepReduction::TimeStepRetryScope but failed to consume a pending "
        "MPI-reduced request.");
  }

  request_state.global_request_consumed = false;
  request_state.time_step_retry_scope_active = false;
}

Core::Mat::TimeStepReduction::MaterialEvaluationScope::MaterialEvaluationScope(MPI_Comm comm)
    : uncaught_exceptions_on_entry_(std::uncaught_exceptions())
{
  FOUR_C_ASSERT_ALWAYS(not request_state.evaluation_scope_active,
      "Nested material time-step-reduction evaluation scopes are not supported. The evaluator "
      "that owns the synchronized error check must install exactly one evaluation scope.");

  request_state.evaluation_comm = comm;
  request_state.evaluation_scope_active = true;
}

Core::Mat::TimeStepReduction::MaterialEvaluationScope::~MaterialEvaluationScope() noexcept(false)
{
  if (std::uncaught_exceptions() > uncaught_exceptions_on_entry_)
  {
    request_state.rank_request_pending = false;
    request_state.evaluation_scope_active = false;
    request_state.evaluation_comm = MPI_COMM_NULL;
    return;
  }

  FOUR_C_ASSERT_ALWAYS(request_state.evaluation_scope_active,
      "Material time-step-reduction evaluation scope ended although no evaluation scope was "
      "active.");

  if (request_state.rank_request_pending)
  {
    request_state.rank_request_pending = false;
    request_state.evaluation_scope_active = false;
    request_state.evaluation_comm = MPI_COMM_NULL;
    FOUR_C_THROW(
        "You opened a TimeStepReduction::MaterialEvaluationScope but failed to collect "
        "a rank-local request.");
  }

  request_state.evaluation_scope_active = false;
  request_state.evaluation_comm = MPI_COMM_NULL;
}

FOUR_C_NAMESPACE_CLOSE
