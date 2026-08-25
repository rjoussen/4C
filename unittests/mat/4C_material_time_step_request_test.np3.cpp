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
  using namespace FourC::Core::Mat::TimeStepReduction;

  // The request channel has process-local state.
  // Tests either consume their requests or cover the scope-exit cleanup branch.
  constexpr int requesting_rank = 1;

  void request_time_step_reduction_on_rank_1(const std::string& reason)
  {
    if (Core::Communication::my_mpi_rank(MPI_COMM_WORLD) == requesting_rank)
    {
      request(reason);
    }
  }

  TEST(MaterialTimeStepRequest, CollectPeekAndConsumeWithoutRequestReturnFalse)
  {
    TimeStepRetryScope retry_scope;
    MaterialEvaluationScope evaluation_scope(MPI_COMM_WORLD);

    EXPECT_FALSE(collect_local_requests());
    EXPECT_FALSE(peek_mpi_reduced_request());
    EXPECT_FALSE(consume_mpi_reduced_request());
  }

  TEST(MaterialTimeStepRequest, ConsumeOutsideTimeStepRetryScopeThrows)
  {
    FOUR_C_EXPECT_THROW_WITH_MESSAGE(consume_mpi_reduced_request(), Core::Exception,
        "outside of an active time-step retry scope");
  }

  TEST(MaterialTimeStepRequest, CollectAndConsumeRequest)
  {
    TimeStepRetryScope retry_scope;
    MaterialEvaluationScope evaluation_scope(MPI_COMM_WORLD);

    request_time_step_reduction_on_rank_1("material failed locally");

    EXPECT_TRUE(collect_local_requests());
    EXPECT_TRUE(consume_mpi_reduced_request());
  }

  TEST(MaterialTimeStepRequest, RepeatedRequestsCollapse)
  {
    TimeStepRetryScope retry_scope;
    MaterialEvaluationScope evaluation_scope(MPI_COMM_WORLD);

    request_time_step_reduction_on_rank_1("first failure");
    request_time_step_reduction_on_rank_1("second failure");

    EXPECT_TRUE(collect_local_requests());
    EXPECT_TRUE(consume_mpi_reduced_request());
  }

  TEST(MaterialTimeStepRequest, CollectOutsideMaterialEvaluationScopeThrows)
  {
    TimeStepRetryScope retry_scope;

    FOUR_C_EXPECT_THROW_WITH_MESSAGE(
        collect_local_requests(), Core::Exception, "outside of an active evaluation scope");
  }

  TEST(MaterialTimeStepRequest, RequestOutsideTimeStepRetryScopeThrows)
  {
    MaterialEvaluationScope evaluation_scope(MPI_COMM_WORLD);

    FOUR_C_EXPECT_THROW_WITH_MESSAGE(request("no retry scope"), Core::Exception,
        "A material requested time-step reduction but the calling algorithm does not support it.\n"
        "Reason: no retry scope");
  }

  TEST(MaterialTimeStepRequest, DisabledTimeStepRetryScopeIsInert)
  {
    TimeStepRetryScope disabled_retry_scope(false);
    MaterialEvaluationScope evaluation_scope(MPI_COMM_WORLD);

    FOUR_C_EXPECT_THROW_WITH_MESSAGE(request("disabled retry scope"), Core::Exception,
        "A material requested time-step reduction but the calling algorithm does not support it.\n"
        "Reason: disabled retry scope");
  }

  TEST(MaterialTimeStepRequest, RequestOutsideMaterialEvaluationScopeThrows)
  {
    TimeStepRetryScope retry_scope;

    FOUR_C_EXPECT_THROW_WITH_MESSAGE(request("no evaluation scope"), Core::Exception,
        "A material requested time-step reduction in an evaluation path that does not support it.\n"
        "Reason: no evaluation scope");
  }

  TEST(MaterialTimeStepRequest, EvaluationScopeEnforcesCollection)
  {
    FOUR_C_EXPECT_THROW_WITH_MESSAGE(
        {
          TimeStepRetryScope retry_scope;
          MaterialEvaluationScope evaluation_scope(MPI_COMM_WORLD);

          request("uncollected request");
        },
        Core::Exception, "failed to collect a rank-local request");
  }

  TEST(MaterialTimeStepRequest, MaterialEvaluationScopeUnwindsSafely)
  {
    FOUR_C_EXPECT_THROW_WITH_MESSAGE(
        {
          TimeStepRetryScope retry_scope;
          MaterialEvaluationScope evaluation_scope(MPI_COMM_WORLD);

          request_time_step_reduction_on_rank_1("request before original exception");

          // this exception will trigger the evaluation scope destructor which would normally also
          // throw due to uncollected requests.
          FOUR_C_THROW("an unrelated evaluation exception");
        },
        Core::Exception, "an unrelated evaluation exception");
  }

  TEST(MaterialTimeStepRequest, NestedTimeStepRetryScopeThrows)
  {
    TimeStepRetryScope retry_scope;

    FOUR_C_EXPECT_THROW_WITH_MESSAGE(TimeStepRetryScope nested_retry_scope, Core::Exception,
        "Nested material time-step-reduction retry scopes are not supported");
  }

  TEST(MaterialTimeStepRequest, NestedMaterialEvaluationScopeThrows)
  {
    MaterialEvaluationScope evaluation_scope(MPI_COMM_WORLD);

    FOUR_C_EXPECT_THROW_WITH_MESSAGE(
        MaterialEvaluationScope nested_evaluation_scope(MPI_COMM_WORLD), Core::Exception,
        "Nested material time-step-reduction evaluation scopes are not supported");
  }

  TEST(MaterialTimeStepRequest, PeekAndConsumeBehaviour)
  {
    TimeStepRetryScope retry_scope;
    MaterialEvaluationScope evaluation_scope(MPI_COMM_WORLD);
    request_time_step_reduction_on_rank_1("publish reduced request");
    EXPECT_TRUE(collect_local_requests());

    // peek can be called multiple times with the same outcome
    EXPECT_TRUE(peek_mpi_reduced_request());
    EXPECT_TRUE(peek_mpi_reduced_request());

    // the first consume clears the request
    EXPECT_TRUE(consume_mpi_reduced_request());
    EXPECT_FALSE(peek_mpi_reduced_request());
  }

  TEST(MaterialTimeStepRequest, CollectWhileGlobalRequestPendingThrows)
  {
    TimeStepRetryScope retry_scope;
    {
      MaterialEvaluationScope evaluation_scope(MPI_COMM_WORLD);
      request_time_step_reduction_on_rank_1("pending request");
      EXPECT_TRUE(collect_local_requests());
    }

    // the global request is not yet consumed, hence a second collect throws:
    FOUR_C_EXPECT_THROW_WITH_MESSAGE(
        {
          MaterialEvaluationScope evaluation_scope(MPI_COMM_WORLD);
          collect_local_requests();
        },
        Core::Exception, "previous reduced request is still pending");

    // consume now to allow for safe exit of the retry scope
    EXPECT_TRUE(consume_mpi_reduced_request());
  }

  TEST(MaterialTimeStepRequest, ConsumingTwiceInSameTimeStepRetryScopeThrows)
  {
    TimeStepRetryScope retry_scope;
    MaterialEvaluationScope evaluation_scope(MPI_COMM_WORLD);
    request_time_step_reduction_on_rank_1("first request");
    EXPECT_TRUE(collect_local_requests());

    // first consume
    EXPECT_TRUE(consume_mpi_reduced_request());

    // second consume in the same retry scope throws
    FOUR_C_EXPECT_THROW_WITH_MESSAGE(consume_mpi_reduced_request(), Core::Exception,
        "more than once in the same time-step retry scope");
  }

  TEST(MaterialTimeStepRequest, TimeStepRetryScopeEnforcesConsumption)
  {
    FOUR_C_EXPECT_THROW_WITH_MESSAGE(
        {
          TimeStepRetryScope retry_scope;
          MaterialEvaluationScope evaluation_scope(MPI_COMM_WORLD);
          request_time_step_reduction_on_rank_1("pending request");
          collect_local_requests();  // marks the global request as pending, but it is not consumed
        },
        Core::Exception, "failed to consume a pending MPI-reduced request");
  }

  TEST(MaterialTimeStepRequest, TimeStepRetryScopeUnwindsSafely)
  {
    FOUR_C_EXPECT_THROW_WITH_MESSAGE(
        {
          TimeStepRetryScope retry_scope;
          MaterialEvaluationScope evaluation_scope(MPI_COMM_WORLD);
          request_time_step_reduction_on_rank_1("pending request");
          collect_local_requests();
          // this exception will trigger the retry scope destructor which would normally also throw
          // due to unconsumed requests.
          FOUR_C_THROW("original algorithm exception");
        },
        Core::Exception, "original algorithm exception");
  }

}  // namespace
