// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_mat_history_transfer.hpp"

#include "4C_fem_general_largerotations.hpp"
#include "4C_linalg_utils_densematrix_eigen.hpp"
#include "4C_linalg_utils_scalar_interpolation.hpp"
#include "4C_utils_exceptions.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>

FOUR_C_NAMESPACE_OPEN

namespace
{
  //! distance below which the target is considered to coincide with a data point
  constexpr double coincidence_tolerance = 1.0e-8;

  //! moving least squares approximation of scalars
  std::vector<double> approximate_scalars(const std::vector<std::vector<double>>& data,
      const std::vector<Core::LinAlg::Matrix<3, 1>>& locations,
      const Core::LinAlg::Matrix<3, 1>& target, const double weight_decay)
  {
    Core::LinAlg::ScalarInterpolationParams params;
    params.exponential_decay_c = weight_decay;
    if (locations.size() < 8)
    {
      const std::vector<double> weights = Core::LinAlg::calculate_normalized_weights(locations,
          target, Core::LinAlg::ScalarInterpolationWeightingFunction::exponential, params);
      std::vector<double> result(data.front().size(), 0.0);
      for (std::size_t i = 0; i < data.size(); ++i)
        for (std::size_t j = 0; j < result.size(); ++j) result[j] += weights[i] * data[i][j];
      return result;
    }
    const Core::LinAlg::ScalarInterpolator<3, 1, 8> interpolator(
        Core::LinAlg::ScalarInterpolationType::moving_least_squares,
        Core::LinAlg::ScalarInterpolationWeightingFunction::exponential, params);
    return interpolator.get_interpolated_scalar(data, locations, target);
  }

  //! polar decomposition T = R U of a tensor with positive determinant: the rotation R and the
  //! logarithm of the stretch U
  struct PolarLogarithm
  {
    Core::LinAlg::Matrix<3, 3> rotation;
    Core::LinAlg::Matrix<3, 3> log_stretch;
  };

  PolarLogarithm polar_logarithm(const Core::LinAlg::Matrix<3, 3>& tensor)
  {
    FOUR_C_ASSERT_ALWAYS(tensor.determinant() > 0.0,
        "The history transfer of invertible tensors requires a positive determinant.");

    // U^2 = T^T T = V diag(lambda) V^T
    Core::LinAlg::Matrix<3, 3> right_cauchy_green;
    right_cauchy_green.multiply_tn(tensor, tensor);
    Core::LinAlg::Matrix<3, 3> eigenvalues;
    Core::LinAlg::Matrix<3, 3> eigenvectors;
    Core::LinAlg::syev(right_cauchy_green, eigenvalues, eigenvectors);

    PolarLogarithm result{Core::LinAlg::Matrix<3, 3>(Core::LinAlg::Initialization::zero),
        Core::LinAlg::Matrix<3, 3>(Core::LinAlg::Initialization::zero)};
    Core::LinAlg::Matrix<3, 3> inverse_stretch(Core::LinAlg::Initialization::zero);
    for (int k = 0; k < 3; ++k)
    {
      const double log_stretch = 0.5 * std::log(eigenvalues(k, k));
      const double inverse = std::exp(-log_stretch);
      for (int r = 0; r < 3; ++r)
      {
        for (int c = 0; c < 3; ++c)
        {
          const double projection = eigenvectors(r, k) * eigenvectors(c, k);
          result.log_stretch(r, c) += log_stretch * projection;
          inverse_stretch(r, c) += inverse * projection;
        }
      }
    }
    result.rotation.multiply(tensor, inverse_stretch);
    return result;
  }

  //! exponential of a symmetric tensor
  Core::LinAlg::Matrix<3, 3> symmetric_exponential(const Core::LinAlg::Matrix<3, 3>& tensor)
  {
    Core::LinAlg::Matrix<3, 3> eigenvalues;
    Core::LinAlg::Matrix<3, 3> eigenvectors;
    Core::LinAlg::syev(tensor, eigenvalues, eigenvectors);

    Core::LinAlg::Matrix<3, 3> result(Core::LinAlg::Initialization::zero);
    for (int k = 0; k < 3; ++k)
      for (int r = 0; r < 3; ++r)
        for (int c = 0; c < 3; ++c)
          result(r, c) += std::exp(eigenvalues(k, k)) * eigenvectors(r, k) * eigenvectors(c, k);
    return result;
  }

  //! rotation vector of a rotation tensor
  Core::LinAlg::Matrix<3, 1> rotation_vector(const Core::LinAlg::Matrix<3, 3>& rotation)
  {
    Core::LinAlg::Matrix<4, 1> quaternion;
    Core::LargeRotations::triadtoquaternion(rotation, quaternion);
    Core::LinAlg::Matrix<3, 1> vector;
    Core::LargeRotations::quaterniontoangle(quaternion, vector);
    return vector;
  }

  /*!
   * @brief Interpolate invertible tensors T = R U
   *
   * The logarithms of the stretches are averaged (log-Euclidean mean) and the rotations relative to
   * the rotation of the closest point as rotation vectors, both with the weights exp(-c d^2). The
   * stretch stays positive definite, its determinant is the weighted geometric mean of the
   * determinants, and components that vanish in all data (e.g., out of the plane in plane strain)
   * stay zero: unlike an interpolation of eigenvalues and eigenvectors, no eigenvectors of
   * different points are matched, which turns them out of their plane where the order of the
   * eigenvalues changes between the points.
   */
  Core::LinAlg::Matrix<3, 3> interpolate_invertible_tensors(
      const std::vector<Core::LinAlg::Matrix<3, 3>>& tensors,
      const std::vector<Core::LinAlg::Matrix<3, 1>>& locations,
      const Core::LinAlg::Matrix<3, 1>& target, const std::size_t closest, const double weight_decay)
  {
    std::vector<double> weights(tensors.size());
    double sum_of_weights = 0.0;
    for (std::size_t i = 0; i < tensors.size(); ++i)
    {
      Core::LinAlg::Matrix<3, 1> difference(locations[i]);
      difference.update(-1.0, target, 1.0);
      weights[i] = std::exp(-weight_decay * difference.dot(difference));
      sum_of_weights += weights[i];
    }

    const Core::LinAlg::Matrix<3, 3> base_rotation = polar_logarithm(tensors[closest]).rotation;
    Core::LinAlg::Matrix<3, 3> log_stretch(Core::LinAlg::Initialization::zero);
    Core::LinAlg::Matrix<3, 1> relative_rotation(Core::LinAlg::Initialization::zero);
    for (std::size_t i = 0; i < tensors.size(); ++i)
    {
      const PolarLogarithm polar = polar_logarithm(tensors[i]);
      const double weight = weights[i] / sum_of_weights;
      log_stretch.update(weight, polar.log_stretch, 1.0);

      Core::LinAlg::Matrix<3, 3> relative;
      relative.multiply_tn(base_rotation, polar.rotation);
      relative_rotation.update(weight, rotation_vector(relative), 1.0);
    }

    Core::LinAlg::Matrix<3, 3> relative;
    Core::LargeRotations::angletotriad(relative_rotation, relative);
    Core::LinAlg::Matrix<3, 3> rotation;
    rotation.multiply(base_rotation, relative);

    Core::LinAlg::Matrix<3, 3> result;
    result.multiply(rotation, symmetric_exponential(log_stretch));
    return result;
  }
}  // namespace

int Mat::history_transfer_size(const HistoryTransferKind kind)
{
  switch (kind)
  {
    case HistoryTransferKind::scalar:
    case HistoryTransferKind::nonnegative_scalar:
    case HistoryTransferKind::nearest:
      return 1;
    case HistoryTransferKind::invertible_tensor:
      return 9;
  }
  FOUR_C_THROW("Unknown history transfer kind.");
}

int Mat::history_transfer_size(const std::vector<HistoryTransferQuantity>& quantities)
{
  return std::accumulate(quantities.begin(), quantities.end(), 0,
      [](int sum, const HistoryTransferQuantity& quantity)
      { return sum + history_transfer_size(quantity.kind); });
}

std::vector<double> Mat::interpolate_history(const std::vector<HistoryTransferQuantity>& quantities,
    const std::vector<std::vector<double>>& data,
    const std::vector<Core::LinAlg::Matrix<3, 1>>& locations,
    const Core::LinAlg::Matrix<3, 1>& target, const double weight_decay)
{
  FOUR_C_ASSERT_ALWAYS(!data.empty() && data.size() == locations.size(),
      "History transfer needs the same (nonzero) number of data sets and locations.");

  const int size = history_transfer_size(quantities);
  for (const auto& point_data : data)
    FOUR_C_ASSERT_ALWAYS(static_cast<int>(point_data.size()) == size,
        "History data has {} values, but the quantities require {}.", point_data.size(), size);

  // find the closest data point
  std::size_t closest = 0;
  double closest_distance = std::numeric_limits<double>::max();
  for (std::size_t i = 0; i < locations.size(); ++i)
  {
    Core::LinAlg::Matrix<3, 1> difference(locations[i]);
    difference.update(-1.0, target, 1.0);
    if (difference.norm2() < closest_distance)
    {
      closest_distance = difference.norm2();
      closest = i;
    }
  }

  if (closest_distance < coincidence_tolerance) return data[closest];

  std::vector<double> result(size);
  int offset = 0;
  for (const HistoryTransferQuantity& quantity : quantities)
  {
    switch (quantity.kind)
    {
      case HistoryTransferKind::scalar:
      case HistoryTransferKind::nonnegative_scalar:
      {
        std::vector<std::vector<double>> scalar_data(data.size(), std::vector<double>(1));
        for (std::size_t i = 0; i < data.size(); ++i) scalar_data[i][0] = data[i][offset];

        double value = approximate_scalars(scalar_data, locations, target, weight_decay)[0];
        if (quantity.kind == HistoryTransferKind::nonnegative_scalar) value = std::max(value, 0.0);
        result[offset] = value;
        break;
      }
      case HistoryTransferKind::invertible_tensor:
      {
        std::vector<Core::LinAlg::Matrix<3, 3>> tensors(data.size());
        for (std::size_t i = 0; i < data.size(); ++i)
          for (int r = 0; r < 3; ++r)
            for (int c = 0; c < 3; ++c) tensors[i](r, c) = data[i][offset + 3 * r + c];

        const Core::LinAlg::Matrix<3, 3> interpolated =
            interpolate_invertible_tensors(tensors, locations, target, closest, weight_decay);

        for (int r = 0; r < 3; ++r)
          for (int c = 0; c < 3; ++c) result[offset + 3 * r + c] = interpolated(r, c);
        break;
      }
      case HistoryTransferKind::nearest:
      {
        result[offset] = data[closest][offset];
        break;
      }
    }
    offset += history_transfer_size(quantity.kind);
  }

  return result;
}

FOUR_C_NAMESPACE_CLOSE
