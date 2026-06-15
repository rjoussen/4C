// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <gtest/gtest.h>

#include "4C_material_time_step_request.hpp"

#include "4C_comm_mpi_utils.hpp"
#include "4C_unittest_utils_assertions_test.hpp"
#include "4C_utils_exceptions.hpp"

#include <mpi.h>

#include <string>

namespace
{
  using namespace FourC;

  void reset_request_state() { Core::Mat::TimeStepReduction::AlgorithmSupportScope reset_scope; }

  void publish_mpi_reduced_request(const bool requested)
  {
    if (requested)
    {
      FOUR_C_THROW(
          "publish_mpi_reduced_request(true) is not supported in this test helper. "
          "Use an explicit AlgorithmSupportScope + request + collect flow in the test.");
    }
    Core::Mat::TimeStepReduction::EvaluationScope evaluation_scope;
    EXPECT_EQ(Core::Mat::TimeStepReduction::collect_local_requests(MPI_COMM_SELF), requested);
  }

  TEST(MaterialTimeStepRequest, InitialStateHasNoMpiReducedRequest)
  {
    reset_request_state();

    EXPECT_FALSE(Core::Mat::TimeStepReduction::has_mpi_reduced_request());
    EXPECT_FALSE(Core::Mat::TimeStepReduction::consume_mpi_reduced_request());
  }

  TEST(MaterialTimeStepRequest, RequestCanBeCollectedOnceInsideSupportedEvaluation)
  {
    reset_request_state();

    Core::Mat::TimeStepReduction::AlgorithmSupportScope support_scope;
    Core::Mat::TimeStepReduction::EvaluationScope evaluation_scope;

    // before the request
    EXPECT_FALSE(Core::Mat::TimeStepReduction::collect_local_requests(MPI_COMM_SELF));

    Core::Mat::request_time_step_reduction("material failed locally");

    // first collection after request returns true, second collection is a programming error
    EXPECT_TRUE(Core::Mat::TimeStepReduction::collect_local_requests(MPI_COMM_SELF));
    FOUR_C_EXPECT_THROW_WITH_MESSAGE(
        Core::Mat::TimeStepReduction::collect_local_requests(MPI_COMM_SELF), Core::Exception,
        "Attempted to publish a new MPI-reduced material timestep-reduction request while a "
        "previous reduced request is still pending.");
    EXPECT_TRUE(Core::Mat::TimeStepReduction::consume_mpi_reduced_request());
  }

  TEST(MaterialTimeStepRequest, RepeatedRequestsBeforeConsumptionCollapseToSingleLocalRequest)
  {
    reset_request_state();

    Core::Mat::TimeStepReduction::AlgorithmSupportScope support_scope;
    Core::Mat::TimeStepReduction::EvaluationScope evaluation_scope;

    Core::Mat::request_time_step_reduction("first failure");
    Core::Mat::request_time_step_reduction("second failure");

    EXPECT_TRUE(Core::Mat::TimeStepReduction::collect_local_requests(MPI_COMM_SELF));
    FOUR_C_EXPECT_THROW_WITH_MESSAGE(
        Core::Mat::TimeStepReduction::collect_local_requests(MPI_COMM_SELF), Core::Exception,
        "Attempted to publish a new MPI-reduced material timestep-reduction request while a "
        "previous reduced request is still pending.");
    EXPECT_TRUE(Core::Mat::TimeStepReduction::consume_mpi_reduced_request());
  }

  TEST(MaterialTimeStepRequest, CollectOutsideEvaluationScopeThrows)
  {
    reset_request_state();

    // open an algorithm support scope only.
    Core::Mat::TimeStepReduction::AlgorithmSupportScope support_scope;

    FOUR_C_EXPECT_THROW_WITH_MESSAGE(
        Core::Mat::TimeStepReduction::collect_local_requests(MPI_COMM_SELF), Core::Exception,
        "Attempted to collect material time step reduction requests, but no evaluation scope is "
        "active. "
        "Install Core::Mat::TimeStepReduction::EvaluationScope around the synchronized "
        "element/material "
        "evaluation and collect requests before scope exit.");
  }

  TEST(MaterialTimeStepRequest, RequestOutsideSupportScopeThrows)
  {
    reset_request_state();

    // open an evaluation scope only.
    Core::Mat::TimeStepReduction::EvaluationScope evaluation_scope;

    FOUR_C_EXPECT_THROW_WITH_MESSAGE(Core::Mat::request_time_step_reduction("no support scope"),
        Core::Exception,
        "A material requested timestep reduction, but the calling algorithm has not declared "
        "support for material-triggered timestep retry. Install "
        "Core::Mat::TimeStepReduction::AlgorithmSupportScope only around algorithms that reduce "
        "dt, "
        "restore state, and repeat the global timestep.");
    EXPECT_FALSE(Core::Mat::TimeStepReduction::collect_local_requests(MPI_COMM_SELF));
  }

  TEST(MaterialTimeStepRequest, RequestInsideDisabledSupportScopeThrows)
  {
    reset_request_state();

    Core::Mat::TimeStepReduction::AlgorithmSupportScope disabled_support_scope(false);
    Core::Mat::TimeStepReduction::EvaluationScope evaluation_scope;

    FOUR_C_EXPECT_THROW_WITH_MESSAGE(
        Core::Mat::request_time_step_reduction("disabled support scope"), Core::Exception,
        "A material requested timestep reduction, but the calling algorithm has not declared "
        "support for material-triggered timestep retry. Install "
        "Core::Mat::TimeStepReduction::AlgorithmSupportScope only around algorithms that reduce "
        "dt, "
        "restore state, and repeat the global timestep.");
    EXPECT_FALSE(Core::Mat::TimeStepReduction::collect_local_requests(MPI_COMM_SELF));
  }

  TEST(MaterialTimeStepRequest, RequestOutsideEvaluationScopeThrows)
  {
    reset_request_state();

    // open a support scope only.
    Core::Mat::TimeStepReduction::AlgorithmSupportScope support_scope;

    FOUR_C_EXPECT_THROW_WITH_MESSAGE(Core::Mat::request_time_step_reduction("no evaluation scope"),
        Core::Exception,
        "A material requested timestep reduction in an evaluation path that does not "
        "support it. "
        "Install Core::Mat::TimeStepReduction::EvaluationScope around the synchronized "
        "element/material "
        "evaluation and collect requests before scope exit.");
  }

  TEST(MaterialTimeStepRequest, UnconsumedRequestThrowsWhenEvaluationScopeEnds)
  {
    reset_request_state();

    FOUR_C_EXPECT_THROW_WITH_MESSAGE(
        {
          Core::Mat::TimeStepReduction::AlgorithmSupportScope support_scope;
          Core::Mat::TimeStepReduction::EvaluationScope evaluation_scope;

          Core::Mat::request_time_step_reduction("unconsumed request");
        },
        Core::Exception,
        "You opened a TimeStepReduction::EvaluationScope but failed to call "
        "TimeStepReduction::collect_local_requests().");
  }

  TEST(MaterialTimeStepRequest, EvaluationScopeDoesNotMaskOriginalExceptionDuringUnwinding)
  {
    reset_request_state();

    FOUR_C_EXPECT_THROW_WITH_MESSAGE(
        {
          Core::Mat::TimeStepReduction::AlgorithmSupportScope support_scope;
          Core::Mat::TimeStepReduction::EvaluationScope evaluation_scope;

          Core::Mat::request_time_step_reduction("request before original exception");
          FOUR_C_THROW("original evaluation exception");
        },
        Core::Exception, "original evaluation exception");
  }

  TEST(MaterialTimeStepRequest, ConsumedRequestDoesNotThrowWhenEvaluationScopeEnds)
  {
    reset_request_state();

    EXPECT_NO_THROW({
      Core::Mat::TimeStepReduction::AlgorithmSupportScope support_scope;
      Core::Mat::TimeStepReduction::EvaluationScope evaluation_scope;

      Core::Mat::request_time_step_reduction("consumed request");
      EXPECT_TRUE(Core::Mat::TimeStepReduction::collect_local_requests(MPI_COMM_SELF));
      EXPECT_TRUE(Core::Mat::TimeStepReduction::consume_mpi_reduced_request());
    });
  }

  TEST(MaterialTimeStepRequest, EvaluationScopeWithoutRequestDoesNotThrowWhenItEnds)
  {
    reset_request_state();

    EXPECT_NO_THROW({
      Core::Mat::TimeStepReduction::AlgorithmSupportScope support_scope;
      Core::Mat::TimeStepReduction::EvaluationScope evaluation_scope;
      EXPECT_FALSE(Core::Mat::TimeStepReduction::collect_local_requests(MPI_COMM_SELF));
    });
  }

  TEST(MaterialTimeStepRequest, EvaluationScopeWithNoRequestCanBeConsumedDefensively)
  {
    reset_request_state();

    EXPECT_NO_THROW({
      Core::Mat::TimeStepReduction::AlgorithmSupportScope support_scope;
      Core::Mat::TimeStepReduction::EvaluationScope evaluation_scope;

      EXPECT_FALSE(Core::Mat::TimeStepReduction::collect_local_requests(MPI_COMM_SELF));
    });
  }

  TEST(MaterialTimeStepRequest, NestedSupportScopeThrows)
  {
    reset_request_state();

    Core::Mat::TimeStepReduction::AlgorithmSupportScope support_scope;

    FOUR_C_EXPECT_THROW_WITH_MESSAGE(
        Core::Mat::TimeStepReduction::AlgorithmSupportScope nested_support_scope, Core::Exception,
        "Nested material timestep-reduction support scopes are not supported. The algorithm that "
        "owns timestep retry must install exactly one support scope.");
  }

  TEST(MaterialTimeStepRequest, NestedEvaluationScopeThrows)
  {
    reset_request_state();

    Core::Mat::TimeStepReduction::AlgorithmSupportScope support_scope;
    Core::Mat::TimeStepReduction::EvaluationScope evaluation_scope;

    FOUR_C_EXPECT_THROW_WITH_MESSAGE(
        Core::Mat::TimeStepReduction::EvaluationScope nested_evaluation_scope, Core::Exception,
        "Nested material timestep-reduction evaluation scopes are not supported. The evaluator "
        "that owns the synchronized error check must install exactly one evaluation scope.");
    EXPECT_FALSE(Core::Mat::TimeStepReduction::collect_local_requests(MPI_COMM_SELF));
  }

  TEST(MaterialTimeStepRequest, EvaluationScopeEntryClearsPreviousLocalRequest)
  {
    reset_request_state();

    Core::Mat::TimeStepReduction::AlgorithmSupportScope support_scope;

    {
      Core::Mat::TimeStepReduction::EvaluationScope evaluation_scope;
      Core::Mat::request_time_step_reduction("first evaluation");
      EXPECT_TRUE(Core::Mat::TimeStepReduction::collect_local_requests(MPI_COMM_SELF));
    }
    EXPECT_TRUE(Core::Mat::TimeStepReduction::consume_mpi_reduced_request());

    {
      Core::Mat::TimeStepReduction::EvaluationScope evaluation_scope;
      EXPECT_FALSE(Core::Mat::TimeStepReduction::collect_local_requests(MPI_COMM_SELF));
    }
  }

  TEST(MaterialTimeStepRequest, SupportScopeEntryClearsMpiReducedRequest)
  {
    reset_request_state();

    {
      Core::Mat::TimeStepReduction::AlgorithmSupportScope support_scope;
      Core::Mat::TimeStepReduction::EvaluationScope evaluation_scope;
      Core::Mat::request_time_step_reduction("publish reduced request");
      EXPECT_TRUE(Core::Mat::TimeStepReduction::collect_local_requests(MPI_COMM_SELF));
      EXPECT_TRUE(Core::Mat::TimeStepReduction::has_mpi_reduced_request());
      EXPECT_TRUE(Core::Mat::TimeStepReduction::consume_mpi_reduced_request());
    }

    {
      Core::Mat::TimeStepReduction::AlgorithmSupportScope support_scope;
      EXPECT_FALSE(Core::Mat::TimeStepReduction::has_mpi_reduced_request());
      EXPECT_FALSE(Core::Mat::TimeStepReduction::consume_mpi_reduced_request());
    }
  }

  TEST(MaterialTimeStepRequest, PeekAndConsumeMpiReducedRequest)
  {
    reset_request_state();

    publish_mpi_reduced_request(false);
    EXPECT_FALSE(Core::Mat::TimeStepReduction::has_mpi_reduced_request());
    EXPECT_FALSE(Core::Mat::TimeStepReduction::consume_mpi_reduced_request());

    {
      Core::Mat::TimeStepReduction::AlgorithmSupportScope support_scope;
      Core::Mat::TimeStepReduction::EvaluationScope evaluation_scope;
      Core::Mat::request_time_step_reduction("publish reduced request");
      EXPECT_TRUE(Core::Mat::TimeStepReduction::collect_local_requests(MPI_COMM_SELF));
      EXPECT_TRUE(Core::Mat::TimeStepReduction::has_mpi_reduced_request());
      EXPECT_TRUE(Core::Mat::TimeStepReduction::has_mpi_reduced_request());
      EXPECT_TRUE(Core::Mat::TimeStepReduction::consume_mpi_reduced_request());
      FOUR_C_EXPECT_THROW_WITH_MESSAGE(Core::Mat::TimeStepReduction::has_mpi_reduced_request(),
          Core::Exception,
          "Attempted to query an MPI-reduced material timestep-reduction request after it was "
          "already consumed. The request channel is one-shot; publish a new request before "
          "querying "
          "again.");
    }

    publish_mpi_reduced_request(false);
    EXPECT_FALSE(Core::Mat::TimeStepReduction::has_mpi_reduced_request());
    EXPECT_FALSE(Core::Mat::TimeStepReduction::consume_mpi_reduced_request());
  }

  TEST(MaterialTimeStepRequest, MpiReducedRequestCanOnlyBeConsumedOnce)
  {
    reset_request_state();

    {
      Core::Mat::TimeStepReduction::AlgorithmSupportScope support_scope;
      Core::Mat::TimeStepReduction::EvaluationScope evaluation_scope;
      Core::Mat::request_time_step_reduction("publish reduced request");
      EXPECT_TRUE(Core::Mat::TimeStepReduction::collect_local_requests(MPI_COMM_SELF));
      EXPECT_TRUE(Core::Mat::TimeStepReduction::consume_mpi_reduced_request());
      FOUR_C_EXPECT_THROW_WITH_MESSAGE(Core::Mat::TimeStepReduction::has_mpi_reduced_request(),
          Core::Exception,
          "Attempted to query an MPI-reduced material timestep-reduction request after it was "
          "already consumed. The request channel is one-shot; publish a new request before "
          "querying "
          "again.");
      FOUR_C_EXPECT_THROW_WITH_MESSAGE(Core::Mat::TimeStepReduction::consume_mpi_reduced_request(),
          Core::Exception,
          "Attempted to consume an MPI-reduced material timestep-reduction request more than once. "
          "The request channel is one-shot and each request may be consumed exactly once.");
    }
  }

  TEST(MaterialTimeStepRequest, CollectLocalRequestWhilePendingThrows)
  {
    reset_request_state();

    FOUR_C_EXPECT_THROW_WITH_MESSAGE(
        {
          Core::Mat::TimeStepReduction::AlgorithmSupportScope support_scope;
          {
            Core::Mat::TimeStepReduction::EvaluationScope evaluation_scope;
            Core::Mat::request_time_step_reduction("pending request");
            EXPECT_TRUE(Core::Mat::TimeStepReduction::collect_local_requests(MPI_COMM_SELF));
          }
          {
            Core::Mat::TimeStepReduction::EvaluationScope evaluation_scope;
            Core::Mat::TimeStepReduction::collect_local_requests(MPI_COMM_SELF);
          }
        },
        Core::Exception,
        "Attempted to publish a new MPI-reduced material timestep-reduction request while a "
        "previous reduced request is still pending. The timestep policy must consume the pending "
        "request before the next reduced request is published.");
  }

  TEST(MaterialTimeStepRequest, QueryMpiReducedRequestDoesNotConsumeSupportScopeObligation)
  {
    reset_request_state();

    FOUR_C_EXPECT_THROW_WITH_MESSAGE(
        {
          Core::Mat::TimeStepReduction::AlgorithmSupportScope support_scope;
          Core::Mat::TimeStepReduction::EvaluationScope evaluation_scope;
          Core::Mat::request_time_step_reduction("pending request");
          Core::Mat::TimeStepReduction::collect_local_requests(MPI_COMM_SELF);
          EXPECT_TRUE(Core::Mat::TimeStepReduction::has_mpi_reduced_request());
        },
        Core::Exception,
        "You opened a TimeStepReduction::AlgorithmSupportScope but failed to reduce the time step "
        "reduction request.");
  }

  TEST(MaterialTimeStepRequest, UnconsumedMpiReducedRequestThrowsWhenSupportScopeEnds)
  {
    reset_request_state();

    FOUR_C_EXPECT_THROW_WITH_MESSAGE(
        {
          Core::Mat::TimeStepReduction::AlgorithmSupportScope support_scope;
          Core::Mat::TimeStepReduction::EvaluationScope evaluation_scope;
          Core::Mat::request_time_step_reduction("pending request");
          Core::Mat::TimeStepReduction::collect_local_requests(MPI_COMM_SELF);
        },
        Core::Exception,
        "You opened a TimeStepReduction::AlgorithmSupportScope but failed to reduce the time step "
        "reduction request.");
  }

  TEST(MaterialTimeStepRequest, SupportScopeDoesNotMaskOriginalExceptionDuringUnwinding)
  {
    reset_request_state();

    FOUR_C_EXPECT_THROW_WITH_MESSAGE(
        {
          Core::Mat::TimeStepReduction::AlgorithmSupportScope support_scope;
          Core::Mat::TimeStepReduction::EvaluationScope evaluation_scope;
          Core::Mat::request_time_step_reduction("pending request");
          Core::Mat::TimeStepReduction::collect_local_requests(MPI_COMM_SELF);
          FOUR_C_THROW("original algorithm exception");
        },
        Core::Exception, "original algorithm exception");
  }

  TEST(MaterialTimeStepRequest, ConsumedMpiReducedRequestDoesNotThrowWhenSupportScopeEnds)
  {
    reset_request_state();

    EXPECT_NO_THROW({
      Core::Mat::TimeStepReduction::AlgorithmSupportScope support_scope;
      Core::Mat::TimeStepReduction::EvaluationScope evaluation_scope;
      Core::Mat::request_time_step_reduction("pending request");
      Core::Mat::TimeStepReduction::collect_local_requests(MPI_COMM_SELF);
      EXPECT_TRUE(Core::Mat::TimeStepReduction::consume_mpi_reduced_request());
    });
  }

  TEST(MaterialTimeStepRequest, SupportScopeWithoutMpiReducedRequestDoesNotThrowWhenItEnds)
  {
    reset_request_state();

    EXPECT_NO_THROW({ Core::Mat::TimeStepReduction::AlgorithmSupportScope support_scope; });
  }

  TEST(MaterialTimeStepRequest, GlobalReductionRequestIsMpiReduced)
  {
    reset_request_state();

    const int my_rank = Core::Communication::my_mpi_rank(MPI_COMM_WORLD);

    Core::Mat::TimeStepReduction::AlgorithmSupportScope support_scope;
    Core::Mat::TimeStepReduction::EvaluationScope evaluation_scope;

    if (my_rank == 1) Core::Mat::request_time_step_reduction("rank one failed");

    EXPECT_TRUE(Core::Mat::TimeStepReduction::collect_local_requests(MPI_COMM_WORLD));

    EXPECT_TRUE(Core::Mat::TimeStepReduction::has_mpi_reduced_request());
    EXPECT_TRUE(Core::Mat::TimeStepReduction::consume_mpi_reduced_request());
  }
}  // namespace
