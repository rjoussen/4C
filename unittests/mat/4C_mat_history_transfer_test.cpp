// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <gtest/gtest.h>

#include "4C_mat_history_transfer.hpp"

#include "4C_linalg_fixedsizematrix.hpp"
#include "4C_linalg_utils_tensor_interpolation.hpp"

#include <cmath>

namespace
{
  using namespace FourC;

  constexpr double weight_decay = 10.0;

  Core::LinAlg::Matrix<3, 1> make_vector(double x, double y, double z)
  {
    Core::LinAlg::Matrix<3, 1> v;
    v(0) = x;
    v(1) = y;
    v(2) = z;
    return v;
  }

  //! points of an element patch, e.g., Gauss points of 3x3x3 elements in units of the element size
  std::vector<Core::LinAlg::Matrix<3, 1>> patch_points()
  {
    std::vector<Core::LinAlg::Matrix<3, 1>> points;
    for (int k = -2; k <= 2; ++k)
      for (int j = -2; j <= 2; ++j)
        for (int i = -2; i <= 2; ++i) points.push_back(make_vector(0.6 * i, 0.6 * j, 0.6 * k));
    return points;
  }

  //! an isochoric tensor F = R * U depending on the location
  Core::LinAlg::Matrix<3, 3> isochoric_tensor(const Core::LinAlg::Matrix<3, 1>& x)
  {
    Core::LinAlg::Matrix<3, 1> rotation_vector = make_vector(0.3 * x(0), -0.2 * x(1), 0.4 * x(2));
    const Core::LinAlg::Matrix<3, 3> R =
        Core::LinAlg::calc_rot_matrix_from_rot_vect(rotation_vector);

    const double l1 = 1.2 + 0.1 * x(0);
    const double l2 = 0.9 + 0.05 * x(1);
    Core::LinAlg::Matrix<3, 3> U(Core::LinAlg::Initialization::zero);
    U(0, 0) = l1;
    U(1, 1) = l2;
    U(2, 2) = 1.0 / (l1 * l2);

    Core::LinAlg::Matrix<3, 3> F;
    F.multiply(R, U);
    return F;
  }

  std::vector<double> flatten(const Core::LinAlg::Matrix<3, 3>& tensor)
  {
    std::vector<double> data;
    for (int r = 0; r < 3; ++r)
      for (int c = 0; c < 3; ++c) data.push_back(tensor(r, c));
    return data;
  }

  Core::LinAlg::Matrix<3, 3> unflatten(const std::vector<double>& data, int offset)
  {
    Core::LinAlg::Matrix<3, 3> tensor;
    for (int r = 0; r < 3; ++r)
      for (int c = 0; c < 3; ++c) tensor(r, c) = data[offset + 3 * r + c];
    return tensor;
  }

  TEST(HistoryTransferTest, Sizes)
  {
    const std::vector<Mat::HistoryTransferQuantity> quantities{
        {"F", Mat::HistoryTransferKind::invertible_tensor}, {"a", Mat::HistoryTransferKind::scalar},
        {"b", Mat::HistoryTransferKind::nonnegative_scalar},
        {"c", Mat::HistoryTransferKind::nearest}};
    EXPECT_EQ(Mat::history_transfer_size(quantities), 12);
  }

  TEST(HistoryTransferTest, ReturnsDataAtCoincidentPoint)
  {
    const std::vector<Mat::HistoryTransferQuantity> quantities{
        {"F", Mat::HistoryTransferKind::invertible_tensor},
        {"a", Mat::HistoryTransferKind::scalar}};
    const auto points = patch_points();
    std::vector<std::vector<double>> data;
    for (const auto& point : points)
    {
      auto values = flatten(isochoric_tensor(point));
      values.push_back(point(0) + 2.0 * point(1));
      data.push_back(values);
    }

    const auto result =
        Mat::interpolate_history(quantities, data, points, points[17], weight_decay);
    ASSERT_EQ(result.size(), data[17].size());
    for (std::size_t i = 0; i < result.size(); ++i) EXPECT_EQ(result[i], data[17][i]);
  }

  TEST(HistoryTransferTest, ReproducesTrilinearScalarField)
  {
    const std::vector<Mat::HistoryTransferQuantity> quantities{
        {"a", Mat::HistoryTransferKind::scalar}};
    const auto field = [](const Core::LinAlg::Matrix<3, 1>& x)
    { return 1.0 + 2.0 * x(0) - x(1) + 0.5 * x(2) + 0.3 * x(0) * x(1) * x(2); };

    const auto points = patch_points();
    std::vector<std::vector<double>> data;
    for (const auto& point : points) data.push_back({field(point)});

    const auto target = make_vector(0.13, -0.27, 0.41);
    const auto result = Mat::interpolate_history(quantities, data, points, target, weight_decay);
    EXPECT_NEAR(result[0], field(target), 1.0e-10);
  }

  TEST(HistoryTransferTest, ClipsNonnegativeScalars)
  {
    const std::vector<Mat::HistoryTransferQuantity> quantities{
        {"a", Mat::HistoryTransferKind::nonnegative_scalar}};
    const auto points = patch_points();
    std::vector<std::vector<double>> data;
    for (const auto& point : points) data.push_back({std::max(0.0, 2.0 * point(0))});

    // the linear extrapolation of the data would be negative far outside of the patch
    const auto result =
        Mat::interpolate_history(quantities, data, points, make_vector(-3, 0, 0), weight_decay);
    EXPECT_GE(result[0], 0.0);
  }

  TEST(HistoryTransferTest, PreservesDeterminantOfIsochoricTensors)
  {
    const std::vector<Mat::HistoryTransferQuantity> quantities{
        {"F", Mat::HistoryTransferKind::invertible_tensor}};
    const auto points = patch_points();
    std::vector<std::vector<double>> data;
    for (const auto& point : points) data.push_back(flatten(isochoric_tensor(point)));

    for (const auto& target : {make_vector(0.1, 0.2, -0.1), make_vector(0.9, -0.4, 0.2)})
    {
      const Core::LinAlg::Matrix<3, 3> F =
          unflatten(Mat::interpolate_history(quantities, data, points, target, weight_decay), 0);
      EXPECT_NEAR(F.determinant(), 1.0, 1.0e-12);

      // the weighted average of the eigenvalues smooths the linearly varying stretch
      Core::LinAlg::Matrix<3, 3> error(isochoric_tensor(target));
      error.update(-1.0, F, 1.0);
      EXPECT_LT(error.norm2(), 0.06);
    }
  }

  //! a plane strain tensor F = R * U with a rotation about z and an in-plane stretch, whose
  //! eigenvalues pass the out-of-plane eigenvalue 1 across the patch
  Core::LinAlg::Matrix<3, 3> plane_strain_tensor(const Core::LinAlg::Matrix<3, 1>& x)
  {
    Core::LinAlg::Matrix<3, 1> rotation_vector = make_vector(0.0, 0.0, 0.2 * x(0));
    const Core::LinAlg::Matrix<3, 3> R =
        Core::LinAlg::calc_rot_matrix_from_rot_vect(rotation_vector);

    const double angle = 0.3 * x(1);
    const double l1 = 1.0 + 0.05 * x(0);
    const double l2 = 1.02 - 0.04 * x(1);
    const double c = std::cos(angle);
    const double s = std::sin(angle);
    Core::LinAlg::Matrix<3, 3> U(Core::LinAlg::Initialization::zero);
    U(0, 0) = c * c * l1 + s * s * l2;
    U(1, 1) = s * s * l1 + c * c * l2;
    U(0, 1) = U(1, 0) = c * s * (l1 - l2);
    U(2, 2) = 1.0;

    Core::LinAlg::Matrix<3, 3> F;
    F.multiply(R, U);
    return F;
  }

  TEST(HistoryTransferTest, KeepsPlaneStrainTensorsPlane)
  {
    const std::vector<Mat::HistoryTransferQuantity> quantities{
        {"F", Mat::HistoryTransferKind::invertible_tensor}};
    const auto points = patch_points();
    std::vector<std::vector<double>> data;
    for (const auto& point : points) data.push_back(flatten(plane_strain_tensor(point)));

    for (const auto& target : {make_vector(0.1, 0.2, -0.1), make_vector(-0.7, 0.5, 0.3),
             make_vector(0.4, -0.9, 0.0)})
    {
      const Core::LinAlg::Matrix<3, 3> F =
          unflatten(Mat::interpolate_history(quantities, data, points, target, weight_decay), 0);
      EXPECT_NEAR(F(0, 2), 0.0, 1.0e-12);
      EXPECT_NEAR(F(1, 2), 0.0, 1.0e-12);
      EXPECT_NEAR(F(2, 0), 0.0, 1.0e-12);
      EXPECT_NEAR(F(2, 1), 0.0, 1.0e-12);
      EXPECT_NEAR(F(2, 2), 1.0, 1.0e-12);
    }
  }

  TEST(HistoryTransferTest, InterpolatesIdentityTensors)
  {
    const std::vector<Mat::HistoryTransferQuantity> quantities{
        {"F", Mat::HistoryTransferKind::invertible_tensor}};
    const auto points = patch_points();

    Core::LinAlg::Matrix<3, 3> identity(Core::LinAlg::Initialization::zero);
    for (int i = 0; i < 3; ++i) identity(i, i) = 1.0;
    const std::vector<std::vector<double>> data(points.size(), flatten(identity));

    const auto result = Mat::interpolate_history(
        quantities, data, points, make_vector(0.1, 0.2, 0.3), weight_decay);
    for (int i = 0; i < 9; ++i) EXPECT_NEAR(result[i], flatten(identity)[i], 1.0e-12);
  }

  TEST(HistoryTransferTest, SmoothsOscillations)
  {
    const std::vector<Mat::HistoryTransferQuantity> quantities{
        {"a", Mat::HistoryTransferKind::nonnegative_scalar}};
    const auto points = patch_points();

    // checkerboard oscillating around 1
    std::vector<std::vector<double>> data;
    for (std::size_t i = 0; i < points.size(); ++i) data.push_back({i % 2 == 0 ? 1.5 : 0.5});

    const auto result =
        Mat::interpolate_history(quantities, data, points, make_vector(0.1, 0.2, -0.1), 1.0);
    EXPECT_NEAR(result[0], 1.0, 0.1);
  }

  TEST(HistoryTransferTest, TakesDiscreteDataFromClosestPoint)
  {
    const std::vector<Mat::HistoryTransferQuantity> quantities{
        {"flag", Mat::HistoryTransferKind::nearest}};
    const auto points = patch_points();
    std::vector<std::vector<double>> data;
    for (std::size_t i = 0; i < points.size(); ++i) data.push_back({static_cast<double>(i)});

    // closest point is (0.6, 0.6, -0.6), i.e., i = 1, j = 1, k = -1
    const auto result = Mat::interpolate_history(
        quantities, data, points, make_vector(0.55, 0.65, -0.55), weight_decay);
    EXPECT_EQ(result[0], 1 * 25 + 3 * 5 + 3);
  }
}  // namespace
