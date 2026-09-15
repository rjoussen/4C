// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <gtest/gtest.h>

#include "4C_material_step_reduction_request.hpp"

#include "4C_comm_mpi_utils.hpp"
#include "4C_unittest_utils_assertions_test.hpp"
#include "4C_utils_exceptions.hpp"

#include <mpi.h>

#include <string>
#include <tuple>

namespace
{
  using namespace FourC;
  using namespace FourC::Core::Mat::StepReduction;

  void request_reduction_on_rank_1(const std::string& reason)
  {
    if (Core::Communication::my_mpi_rank(MPI_COMM_WORLD) == 1) request(reason);
  }

  void evaluate_without_request()
  {
    run_and_synchronize_request(MPI_COMM_WORLD, [] {});
  }

  TEST(MaterialStepReductionRequest, EvaluationWithoutRequestReturnsFalse)
  {
    const bool request_occurred =
        run_and_detect_synchronized_request(true, evaluate_without_request);

    EXPECT_FALSE(request_occurred);
  }

  TEST(MaterialStepReductionRequest, RequestOnOneRankIsDetectedOnEveryRank)
  {
    const bool request_occurred = run_and_detect_synchronized_request(true,
        []
        {
          run_and_synchronize_request(
              MPI_COMM_WORLD, [] { request_reduction_on_rank_1("rank 1 requested."); });
        });

    EXPECT_TRUE(request_occurred);
  }

  TEST(MaterialStepReductionRequest, SequentialSynchronizationRegionsAreAllowed)
  {
    const bool request_occurred = run_and_detect_synchronized_request(true,
        []
        {
          run_and_synchronize_request(MPI_COMM_WORLD, [] {});
          run_and_synchronize_request(MPI_COMM_WORLD, [] {});
        });

    EXPECT_FALSE(request_occurred);
  }

  TEST(MaterialStepReductionRequest, SynchronizationWithoutRequestDoesNotRequireDetectionRegion)
  {
    run_and_synchronize_request(MPI_COMM_WORLD, [] {});
  }

  TEST(MaterialStepReductionRequest, RequestOutsideSynchronizationRegionThrows)
  {
    FOUR_C_EXPECT_THROW_WITH_MESSAGE(std::ignore = run_and_detect_synchronized_request(
                                         true, [] { request("all ranks requested."); }),
        Core::Exception,
        "A material requested step reduction in an evaluation path that does not support "
        "it.\nReason: all ranks requested.\n");
  }

  TEST(MaterialStepReductionRequest, RequestInsideSynchronizationOutsideDetectionRegionThrows)
  {
    FOUR_C_EXPECT_THROW_WITH_MESSAGE(
        { run_and_synchronize_request(MPI_COMM_WORLD, [] { request("all ranks requested."); }); },
        Core::Exception,
        "A material requested step reduction but the calling algorithm does not support "
        "it.\nReason: all ranks requested.\n");

    EXPECT_NO_THROW(run_and_synchronize_request(MPI_COMM_WORLD, [] {}));
  }

  TEST(MaterialStepReductionRequest, NestedRequestDetectionRegionsThrow)
  {
    FOUR_C_EXPECT_THROW_WITH_MESSAGE(
        std::ignore = run_and_detect_synchronized_request(
            true, [] { std::ignore = run_and_detect_synchronized_request(true, [] {}); }),
        Core::Exception,
        "Nested Core::Mat::StepReduction::run_and_detect_synchronized_request() calls are not "
        "supported");
  }

  TEST(MaterialStepReductionRequest, NestedRequestSynchronizationRegionsThrow)
  {
    FOUR_C_EXPECT_THROW_WITH_MESSAGE(
        std::ignore = run_and_detect_synchronized_request(true,
            []
            {
              run_and_synchronize_request(
                  MPI_COMM_WORLD, [] { run_and_synchronize_request(MPI_COMM_WORLD, [] {}); });
            }),
        Core::Exception,
        "Nested Core::Mat::StepReduction::run_and_synchronize_request() calls are not "
        "supported");
  }

  TEST(MaterialStepReductionRequest, UnrelatedExceptionPropagatesAndResetsEvaluationState)
  {
    FOUR_C_EXPECT_THROW_WITH_MESSAGE(std::ignore = run_and_detect_synchronized_request(true,
                                         []
                                         {
                                           run_and_synchronize_request(MPI_COMM_WORLD, []
                                               { FOUR_C_THROW("unrelated evaluation exception"); });
                                         }),
        Core::Exception, "unrelated evaluation exception");

    // A following operation can still enter a request-detection region.
    EXPECT_FALSE(run_and_detect_synchronized_request(true, evaluate_without_request));
  }

  TEST(MaterialStepReductionRequest, DisabledRequestDetectionRegionThrowsDiagnostic)
  {
    FOUR_C_EXPECT_THROW_WITH_MESSAGE(
        std::ignore = run_and_detect_synchronized_request(false,
            []
            {
              run_and_synchronize_request(MPI_COMM_WORLD, [] { request("all ranks requested."); });
            }),
        Core::Exception,
        "A material requested step reduction, but the calling algorithm disabled handling of these "
        "requests.\nReason: all ranks requested.");

    EXPECT_FALSE(run_and_detect_synchronized_request(true, evaluate_without_request));
  }
}  // namespace
