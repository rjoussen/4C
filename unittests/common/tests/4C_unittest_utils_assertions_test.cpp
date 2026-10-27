// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <gtest/gtest.h>

#include "4C_unittest_utils_assertions_test.hpp"

#include <array>
#include <limits>
#include <vector>

namespace
{
  using namespace FourC;
  using TESTING::INTERNAL::assert_near;
  using TESTING::INTERNAL::assert_near_relative;

  TEST(UnittestUtilsAssertionsTest, RelativeComparisonUsesExpectedScale)
  {
    constexpr double rel_tol = 0.5;

    const std::vector<double> smaller_vector = {0.0, 1.0};
    const std::vector<double> larger_vector = {0.0, 1.75};
    FOUR_C_EXPECT_NEAR_RELATIVE(smaller_vector, larger_vector, rel_tol);
    EXPECT_FALSE(assert_near_relative(
        "larger_vector", "smaller_vector", "rel_tol", larger_vector, smaller_vector, rel_tol));

    Core::LinAlg::Matrix<1, 2> smaller_matrix(Core::LinAlg::Initialization::zero);
    Core::LinAlg::Matrix<1, 2> larger_matrix(Core::LinAlg::Initialization::zero);
    smaller_matrix(0, 1) = 1.0;
    larger_matrix(0, 1) = 1.75;
    FOUR_C_EXPECT_NEAR_RELATIVE(smaller_matrix, larger_matrix, rel_tol);
    EXPECT_FALSE(assert_near_relative(
        "larger_matrix", "smaller_matrix", "rel_tol", larger_matrix, smaller_matrix, rel_tol));

    const Core::LinAlg::Tensor<double, 2> smaller_tensor = {{0.0, 1.0}};
    const Core::LinAlg::Tensor<double, 2> larger_tensor = {{0.0, 1.75}};
    FOUR_C_EXPECT_NEAR_RELATIVE(smaller_tensor, larger_tensor, rel_tol);
    EXPECT_FALSE(assert_near_relative(
        "larger_tensor", "smaller_tensor", "rel_tol", larger_tensor, smaller_tensor, rel_tol));
  }

  TEST(UnittestUtilsAssertionsTest, RelativeComparisonRejectsAllZeroExpectedValue)
  {
    constexpr double rel_tol = 1.0e-6;

    const std::vector<double> zero_vector = {0.0, 0.0};
    const auto vector_result = assert_near_relative(
        "zero_vector", "zero_vector", "rel_tol", zero_vector, zero_vector, rel_tol);
    EXPECT_FALSE(vector_result);
    EXPECT_THAT(vector_result.message(), ::testing::HasSubstr("all-zero entries"));

    const Core::LinAlg::Matrix<1, 2> zero_matrix(Core::LinAlg::Initialization::zero);
    EXPECT_FALSE(assert_near_relative(
        "zero_matrix", "zero_matrix", "rel_tol", zero_matrix, zero_matrix, rel_tol));

    const Core::LinAlg::Tensor<double, 2> zero_tensor = {{0.0, 0.0}};
    EXPECT_FALSE(assert_near_relative(
        "zero_tensor", "zero_tensor", "rel_tol", zero_tensor, zero_tensor, rel_tol));
  }

  TEST(UnittestUtilsAssertionsTest, AbsoluteComparisonRespectsTolerance)
  {
    constexpr double loose_tolerance = 0.5;
    constexpr double tight_tolerance = 0.125;

    const std::vector<double> actual_vector = {0.0, 0.25};
    const std::vector<double> expected_vector = {0.0, 0.0};
    FOUR_C_EXPECT_NEAR(actual_vector, expected_vector, loose_tolerance);
    EXPECT_FALSE(assert_near("actual_vector", "expected_vector", "tight_tolerance", actual_vector,
        expected_vector, tight_tolerance));

    const std::array<double, 2> actual_array = {0.0, 0.25};
    const std::array<double, 2> expected_array = {0.0, 0.0};
    FOUR_C_EXPECT_ITERABLE_NEAR(
        actual_array.data(), expected_array.data(), actual_array.size(), loose_tolerance);
    EXPECT_FALSE(assert_near("actual_array", "expected_array", "length", "tight_tolerance",
        actual_array.data(), expected_array.data(), actual_array.size(), tight_tolerance));

    Core::LinAlg::Matrix<1, 2> actual_matrix(Core::LinAlg::Initialization::zero);
    const Core::LinAlg::Matrix<1, 2> expected_matrix(Core::LinAlg::Initialization::zero);
    actual_matrix(0, 1) = 0.25;
    FOUR_C_EXPECT_NEAR(actual_matrix, expected_matrix, loose_tolerance);
    FOUR_C_EXPECT_NEAR(actual_matrix, expected_array, loose_tolerance);
    EXPECT_FALSE(assert_near("actual_matrix", "expected_matrix", "tight_tolerance", actual_matrix,
        expected_matrix, tight_tolerance));
    EXPECT_FALSE(assert_near("actual_matrix", "expected_array", "tight_tolerance", actual_matrix,
        expected_array, tight_tolerance));

    Core::LinAlg::SerialDenseMatrix actual_dense_matrix(1, 2, true);
    const Core::LinAlg::SerialDenseMatrix expected_dense_matrix(1, 2, true);
    actual_dense_matrix(0, 1) = 0.25;
    FOUR_C_EXPECT_NEAR(actual_dense_matrix, expected_dense_matrix, loose_tolerance);
    EXPECT_FALSE(assert_near("actual_dense_matrix", "expected_dense_matrix", "tight_tolerance",
        actual_dense_matrix, expected_dense_matrix, tight_tolerance));

    Core::LinAlg::SerialDenseVector actual_dense_vector(2, true);
    const Core::LinAlg::SerialDenseVector expected_dense_vector(2, true);
    actual_dense_vector(1) = 0.25;
    FOUR_C_EXPECT_NEAR(actual_dense_vector, expected_dense_vector, loose_tolerance);
    EXPECT_FALSE(assert_near("actual_dense_vector", "expected_dense_vector", "tight_tolerance",
        actual_dense_vector, expected_dense_vector, tight_tolerance));

    const Core::LinAlg::Tensor<double, 2> actual_tensor = {{0.0, 0.25}};
    const Core::LinAlg::Tensor<double, 2> expected_tensor = {{0.0, 0.0}};
    FOUR_C_EXPECT_NEAR(actual_tensor, expected_tensor, loose_tolerance);
    EXPECT_FALSE(assert_near("actual_tensor", "expected_tensor", "tight_tolerance", actual_tensor,
        expected_tensor, tight_tolerance));
  }

  TEST(UnittestUtilsAssertionsTest, ComparisonRejectsNonfiniteEntries)
  {
    constexpr double nan = std::numeric_limits<double>::quiet_NaN();
    constexpr double inf = std::numeric_limits<double>::infinity();
    constexpr double rel_tol = 1.0e-6;
    constexpr double abs_tol = 1.0e-6;

    const std::vector<double> finite_vector = {1.0, 0.0};
    const std::vector<double> nan_vector = {1.0, nan};
    const std::vector<double> inf_vector = {1.0, inf};
    const auto vector_result = assert_near_relative(
        "finite_vector", "nan_vector", "rel_tol", finite_vector, nan_vector, rel_tol);
    EXPECT_FALSE(vector_result);
    EXPECT_THAT(vector_result.message(), ::testing::HasSubstr("Expected finite entries"));
    EXPECT_FALSE(assert_near_relative(
        "nan_vector", "finite_vector", "rel_tol", nan_vector, finite_vector, rel_tol));
    EXPECT_FALSE(assert_near_relative(
        "inf_vector", "finite_vector", "rel_tol", inf_vector, finite_vector, rel_tol));
    EXPECT_FALSE(
        assert_near("nan_vector", "finite_vector", "abs_tol", nan_vector, finite_vector, abs_tol));
    EXPECT_FALSE(
        assert_near("inf_vector", "inf_vector", "abs_tol", inf_vector, inf_vector, abs_tol));

    const std::array<double, 2> finite_array = {1.0, 0.0};
    const std::array<double, 2> nan_array = {1.0, nan};
    EXPECT_FALSE(assert_near("nan_array", "finite_array", "length", "abs_tol", nan_array.data(),
        finite_array.data(), nan_array.size(), abs_tol));

    Core::LinAlg::Matrix<1, 2> finite_matrix(Core::LinAlg::Initialization::zero);
    Core::LinAlg::Matrix<1, 2> nan_matrix(Core::LinAlg::Initialization::zero);
    finite_matrix(0, 0) = 1.0;
    nan_matrix(0, 0) = 1.0;
    nan_matrix(0, 1) = nan;
    EXPECT_FALSE(assert_near_relative(
        "finite_matrix", "nan_matrix", "rel_tol", finite_matrix, nan_matrix, rel_tol));
    EXPECT_FALSE(assert_near_relative(
        "nan_matrix", "finite_matrix", "rel_tol", nan_matrix, finite_matrix, rel_tol));
    EXPECT_FALSE(
        assert_near("nan_matrix", "finite_matrix", "abs_tol", nan_matrix, finite_matrix, abs_tol));
    EXPECT_FALSE(
        assert_near("nan_matrix", "finite_array", "abs_tol", nan_matrix, finite_array, abs_tol));

    Core::LinAlg::SerialDenseMatrix finite_dense_matrix(1, 2, true);
    Core::LinAlg::SerialDenseMatrix nan_dense_matrix(1, 2, true);
    finite_dense_matrix(0, 0) = 1.0;
    nan_dense_matrix(0, 0) = 1.0;
    nan_dense_matrix(0, 1) = nan;
    EXPECT_FALSE(assert_near("nan_dense_matrix", "finite_dense_matrix", "abs_tol", nan_dense_matrix,
        finite_dense_matrix, abs_tol));

    Core::LinAlg::SerialDenseVector finite_dense_vector(2, true);
    Core::LinAlg::SerialDenseVector nan_dense_vector(2, true);
    finite_dense_vector(0) = 1.0;
    nan_dense_vector(0) = 1.0;
    nan_dense_vector(1) = nan;
    EXPECT_FALSE(assert_near("nan_dense_vector", "finite_dense_vector", "abs_tol", nan_dense_vector,
        finite_dense_vector, abs_tol));

    const Core::LinAlg::Tensor<double, 2> finite_tensor = {{1.0, 0.0}};
    const Core::LinAlg::Tensor<double, 2> nan_tensor = {{1.0, nan}};
    EXPECT_FALSE(assert_near_relative(
        "finite_tensor", "nan_tensor", "rel_tol", finite_tensor, nan_tensor, rel_tol));
    EXPECT_FALSE(assert_near_relative(
        "nan_tensor", "finite_tensor", "rel_tol", nan_tensor, finite_tensor, rel_tol));
    EXPECT_FALSE(
        assert_near("nan_tensor", "finite_tensor", "abs_tol", nan_tensor, finite_tensor, abs_tol));
  }
}  // namespace
