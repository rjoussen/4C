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

namespace Core::Mat
{
  /*!
   * \brief Ask the active evaluation owner to retry the global timestep with a smaller dt.
   *
   * This is the material-facing API. The request is deliberately local to the current thread and
   * MPI rank at this point; the structural model evaluator later consumes local requests, reduces
   * them across all ranks, and converts the reduced decision into a structural convergence status.
   *
   * Materials should call this only when they are about to return through their ordinary failed
   * evaluation path. They must not throw their own control-flow exception to escape element
   * evaluation, because doing so can leave other MPI ranks in different collective calls.
   *
   * The function throws a diagnostic error at the request site if the current call stack does not
   * provide both required contracts:
   * - TimeStepReduction::AlgorithmSupportScope from an algorithm that owns timestep retry; and
   * - TimeStepReduction::EvaluationScope around an evaluation whose error check consumes and
   *   MPI-reduces material requests before returning.
   *
   * @param reason Short diagnostic text explaining why a smaller timestep may make this material
   * evaluation safe. The reason is printed immediately at the request site. The MPI reduction later
   * combines only the boolean decision that at least one rank requested reduction.
   */
  void request_time_step_reduction(const std::string& reason);

}  // namespace Core::Mat

namespace Core::Mat::TimeStepReduction
{
  /*!
   * \brief Collect local material timestep-reduction requests and publish reduced request state.
   *
   * This is the evaluator-facing convenience API: it consumes the rank-local request in the active
   * EvaluationScope, performs an MPI max reduction over all ranks, publishes the reduced request
   * state, and returns whether any rank requested timestep reduction.
   *
   * @param comm communicator used for the synchronized reduction
   * @return true if at least one rank requested timestep reduction
   */
  bool collect_local_requests(MPI_Comm comm);

  /*!
   * \brief Return whether an MPI-reduced material timestep-reduction request is pending.
   *
   * This is used by transport layers, such as the new-structure NOX bridge, that must preserve the
   * information without owning the timestep-reduction policy.
   *
   * @throws if the request has already been consumed
   */
  bool has_mpi_reduced_request();

  /*!
   * \brief Consume the pending MPI-reduced material timestep-reduction decision.
   *
   * This must be called by the policy site that actually handles the request, for example by
   * reducing dt and repeating the timestep or by aborting with an explanatory error.
   *
   * Consuming without a pending request is permitted and returns false, but still marks the
   * MPI-reduced request channel as consumed to preserve one-shot semantics.
   *
   * @return true if a pending request was consumed, false if no request is pending
   * @throws if the request has already been consumed
   */
  bool consume_mpi_reduced_request();

  /*!
   * \brief RAII scope that marks an algorithm as owning material timestep-reduction retry policy.
   *
   * This long-lived scope is installed by the time loop or coupled algorithm that can actually act
   * on a pending MPI-reduced material request. It answers the policy question: if a material
   * requests a smaller global timestep, is there code above this evaluation that will reduce dt,
   * restore state, and repeat the step consistently on all MPI ranks?
   *
   * Passing \p enabled=false makes the scope inert. Unsupported algorithms can therefore use the
   * same call pattern without declaring support accidentally; materials still fail at the request
   * site with the unsupported-algorithm diagnostic.
   *
   * The support scope intentionally does not consume material requests. Consumption must happen in
   * the shorter EvaluationScope at the synchronized evaluation boundary. Keeping these two
   * lifetimes separate lets unsupported algorithms fail immediately, while missed evaluator
   * consumption is caught before the individual residual/Jacobian evaluation returns.
   */
  class AlgorithmSupportScope
  {
   public:
    explicit AlgorithmSupportScope(bool enabled = true);
    AlgorithmSupportScope(const AlgorithmSupportScope&) = delete;
    AlgorithmSupportScope& operator=(const AlgorithmSupportScope&) = delete;
    ~AlgorithmSupportScope() noexcept(false);

   private:
    const bool enabled_;
    int uncaught_exceptions_on_entry_ = 0;
  };

  /*!
   * \brief RAII scope that owns request consumption for one synchronized material evaluation.
   *
   * This short-lived scope is installed by evaluator infrastructure around one residual/Jacobian
   * fill. It answers the lifecycle question: if a material requests timestep reduction during
   * this evaluation, will the evaluator call collect_local_requests(), MPI-reduce the request, and
   * convert the failed fill into a structural status before returning to NOX?
   *
   * The scope is independent of AlgorithmSupportScope, i.e. it can be opened also if the algorithm
   * does not support material-triggered timestep reduction. In that case, material requests still
   * fail because they require both scopes to be active.
   */
  class EvaluationScope
  {
   public:
    EvaluationScope();
    EvaluationScope(const EvaluationScope&) = delete;
    EvaluationScope& operator=(const EvaluationScope&) = delete;
    ~EvaluationScope() noexcept(false);

   private:
    int uncaught_exceptions_on_entry_ = 0;
  };

}  // namespace Core::Mat::TimeStepReduction

FOUR_C_NAMESPACE_CLOSE

#endif
