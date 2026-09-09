// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_timestepping_step_reduction_request.hpp"

#include "4C_comm_mpi_utils.hpp"
#include "4C_io_input_spec_builders.hpp"
#include "4C_utils_exceptions.hpp"

#include <exception>
#include <iostream>
#include <sstream>
#include <string>
#include <utility>

FOUR_C_NAMESPACE_OPEN

namespace
{
  /**
  @brief Internal typed exception to signal rank-local requests to \c run_and_synchronize_request().

  This signal can only be thrown inside a \c run_and_synchronize_request() call that is nested in an
  enabled \c run_and_detect_synchronized_request() call.
  */
  class RankLocalRequestSignal : public std::exception
  {
   public:
    RankLocalRequestSignal(std::string request_reason) : reason_(std::move(request_reason)) {}

    [[nodiscard]] const char* what() const noexcept override
    {
      return "rank-local material time-step reduction request";
    }

    /// Return the reason for the request, as provided through the public API.
    [[nodiscard]] const std::string& reason() const noexcept { return reason_; }

   private:
    std::string reason_;
  };

  /**
  @brief Internal typed exception to signal a synchronized request to
  \c run_and_detect_synchronized_request().

  This signal is only thrown by \c run_and_synchronize_request() as a result of at least one rank
  throwing \c RankLocalRequestSignal.
  */
  class SynchronizedRequestSignal : public std::exception
  {
   public:
    [[nodiscard]] const char* what() const noexcept override
    {
      return "synchronized material time-step reduction request";
    }
  };

  /// True while \c RankLocalRequestSignal can be caught and synchronized, i.e. inside a \c
  /// run_and_synchronize_request() call.
  bool request_synchronization_active = false;

  /// Internal RAII guard managing the \c request_synchronization_active flag.
  struct RequestSynchronizationGuard
  {
    explicit RequestSynchronizationGuard()
    {
      FOUR_C_ASSERT_ALWAYS(not request_synchronization_active,
          "Nested Core::Mat::TimeStepReduction::run_and_synchronize_request() calls are not "
          "supported.");

      request_synchronization_active = true;
    }
    RequestSynchronizationGuard(const RequestSynchronizationGuard&) = delete;
    RequestSynchronizationGuard& operator=(const RequestSynchronizationGuard&) = delete;
    ~RequestSynchronizationGuard() { request_synchronization_active = false; }
  };

  /// True while \c SynchronizedRequestSignal can be caught and detected, i.e. inside a \c
  /// run_and_detect_synchronized_request() call.
  bool synchronized_request_detection_active = false;

  /// True if material time-step reduction requests are allowed through the corresponding input
  /// option.
  bool time_step_reduction_request_enabled = false;

  /// Internal RAII guard managing the \c synchronized_request_detection_active flag.
  struct SynchronizedRequestDetectionGuard
  {
    explicit SynchronizedRequestDetectionGuard(const bool enabled)
    {
      FOUR_C_ASSERT_ALWAYS(not synchronized_request_detection_active,
          "Nested Core::Mat::TimeStepReduction::run_and_detect_synchronized_request() calls are "
          "not supported.");

      synchronized_request_detection_active = true;
      time_step_reduction_request_enabled = enabled;
    }

    SynchronizedRequestDetectionGuard(const SynchronizedRequestDetectionGuard&) = delete;
    SynchronizedRequestDetectionGuard& operator=(const SynchronizedRequestDetectionGuard&) = delete;

    ~SynchronizedRequestDetectionGuard()
    {
      time_step_reduction_request_enabled = false;
      synchronized_request_detection_active = false;
    }
  };

}  // namespace

[[nodiscard]] Core::IO::InputSpec Core::Mat::TimeStepReduction::allow_requests_input_spec()
{
  using namespace Core::IO::InputSpecBuilders;

  return parameter<bool>("ALLOW_MATERIAL_TIME_STEP_REDUCTION",
      {.description = "Reduce the time-step size according to the time_step_control settings when "
                      "a material requests time-step reduction during evaluation",
          .default_value = false});
}

void Core::Mat::TimeStepReduction::request(const std::string& reason)
{
  FOUR_C_ASSERT_ALWAYS(request_synchronization_active,
      "A material requested time-step reduction in an evaluation path that does not support it.\n"
      "Reason: {}\n"
      "Run the synchronized element/material evaluation through "
      "Core::Mat::TimeStepReduction::run_and_synchronize_request().",
      reason);

  FOUR_C_ASSERT_ALWAYS(synchronized_request_detection_active,
      "A material requested time-step reduction but the calling algorithm does not support it.\n"
      "Reason: {}\n"
      "Run the operation through "
      "Core::Mat::TimeStepReduction::run_and_detect_synchronized_request().",
      reason);

  FOUR_C_ASSERT_ALWAYS(time_step_reduction_request_enabled,
      "A material requested time-step reduction, but ALLOW_MATERIAL_TIME_STEP_REDUCTION is "
      "disabled.\n"
      "Reason: {}",
      reason);

  // Throw a rank-local exception to leave the evaluation and return to synchronization.
  throw RankLocalRequestSignal(reason);
}

void Core::Mat::TimeStepReduction::run_and_synchronize_request(
    MPI_Comm comm, const std::function<void()>& evaluation)
{
  const RequestSynchronizationGuard request_synchronization_guard;

  // Early return if not called from an enabled detection context, skipping the collective.
  if (not synchronized_request_detection_active or not time_step_reduction_request_enabled)
  {
    evaluation();
    return;
  }

  // Perform the evaluation and collect any rank-local requests.
  int local_request = 0;
  try
  {
    evaluation();
  }
  catch (const RankLocalRequestSignal& request)
  {
    local_request = 1;

    // Print a message to the console indicating the rank-local request reason.
    std::ostringstream message;
    message << std::string(72, '*') << "\n"
            << "A material requested a time-step reduction on rank "
            << Core::Communication::my_mpi_rank(comm) << ".\n"
            << "Reason: " << request.reason() << "\n"
            << std::string(72, '*') << "\n";
    std::cout << message.str() << std::flush;
  }

  // Synchronize rank-local requests to a single global request.
  const bool global_request = static_cast<bool>(Core::Communication::max_all(local_request, comm));

  // Throw collectively if any rank requested time-step reduction.
  if (global_request)
  {
    throw SynchronizedRequestSignal();
  }
}

[[nodiscard]] bool Core::Mat::TimeStepReduction::run_and_detect_synchronized_request(
    const bool enabled, const std::function<void()>& operation)
{
  const SynchronizedRequestDetectionGuard synchronized_request_detection_guard(enabled);

  try
  {
    operation();
  }
  catch (const SynchronizedRequestSignal&)
  {
    return true;
  }

  return false;
}

FOUR_C_NAMESPACE_CLOSE
