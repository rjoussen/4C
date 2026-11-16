// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_ale_mesh_quality.hpp"

#include "4C_fem_general_utils_fem_shapefunctions.hpp"
#include "4C_utils_exceptions.hpp"

#include <Sacado.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <optional>
#include <vector>

FOUR_C_NAMESPACE_OPEN

namespace
{
  template <typename T>
  using Matrix3 = std::array<std::array<T, 3>, 3>;

  template <typename T>
  T determinant(const Matrix3<T>& A)
  {
    return A[0][0] * (A[1][1] * A[2][2] - A[1][2] * A[2][1]) -
           A[0][1] * (A[1][0] * A[2][2] - A[1][2] * A[2][0]) +
           A[0][2] * (A[1][0] * A[2][1] - A[1][1] * A[2][0]);
  }

  Matrix3<double> inverse(const Matrix3<double>& A)
  {
    const double det = determinant(A);
    Matrix3<double> inv;
    for (int i = 0; i < 3; ++i)
    {
      for (int j = 0; j < 3; ++j)
      {
        // cofactor of the entry (j, i)
        const int r0 = (j + 1) % 3, r1 = (j + 2) % 3, c0 = (i + 1) % 3, c1 = (i + 2) % 3;
        inv[i][j] = (A[r0][c0] * A[r1][c1] - A[r0][c1] * A[r1][c0]) / det;
      }
    }
    return inv;
  }

  /*!
   * @brief Corner data of a cell type: derivatives of the shape functions at the corners and the
   * inverse Jacobian of the ideal element
   */
  struct CornerData
  {
    int num_nodes;
    //! derivatives of the shape functions at each corner (corner, direction, node)
    std::vector<Core::LinAlg::SerialDenseMatrix> derivatives;
    Matrix3<double> ideal_jacobian_inverse;
    //! Jacobian of the ideal element with unit edge length
    Matrix3<double> unit_ideal_jacobian;
    //! volume of the parameter space
    double parameter_space_volume;
  };

  CornerData make_corner_data(const Core::FE::CellType celltype)
  {
    CornerData data;
    std::vector<std::array<double, 3>> corners;
    switch (celltype)
    {
      case Core::FE::CellType::hex8:
      {
        data.num_nodes = 8;
        corners = {{-1, -1, -1}, {1, -1, -1}, {1, 1, -1}, {-1, 1, -1}, {-1, -1, 1}, {1, -1, 1},
            {1, 1, 1}, {-1, 1, 1}};
        // the ideal element is a cube, whose Jacobian is a multiple of the identity
        data.ideal_jacobian_inverse = {{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}};
        // the parameter space cube has edge length 2
        data.unit_ideal_jacobian = {{{0.5, 0, 0}, {0, 0.5, 0}, {0, 0, 0.5}}};
        data.parameter_space_volume = 8.0;
        break;
      }
      case Core::FE::CellType::tet4:
      {
        data.num_nodes = 4;
        corners = {{0, 0, 0}};
        // Jacobian of the regular tetrahedron with unit edges: columns (1,0,0), (1/2,sqrt(3)/2,0),
        // (1/2,sqrt(3)/6,sqrt(2/3))
        const double s3 = std::sqrt(3.0);
        const double s23 = std::sqrt(2.0 / 3.0);
        const Matrix3<double> W = {{{1.0, 0.5, 0.5}, {0.0, s3 / 2.0, s3 / 6.0}, {0.0, 0.0, s23}}};
        // inverse of the upper triangular matrix W
        Matrix3<double> W_inv{};
        W_inv[0][0] = 1.0 / W[0][0];
        W_inv[1][1] = 1.0 / W[1][1];
        W_inv[2][2] = 1.0 / W[2][2];
        W_inv[0][1] = -W[0][1] / (W[0][0] * W[1][1]);
        W_inv[1][2] = -W[1][2] / (W[1][1] * W[2][2]);
        W_inv[0][2] = (W[0][1] * W[1][2] - W[0][2] * W[1][1]) / (W[0][0] * W[1][1] * W[2][2]);
        data.ideal_jacobian_inverse = W_inv;
        data.unit_ideal_jacobian = W;
        data.parameter_space_volume = 1.0 / 6.0;
        break;
      }
      default:
        FOUR_C_THROW("The mesh quality is not implemented for cell type {}.",
            Core::FE::cell_type_to_string(celltype));
    }

    for (const auto& corner : corners)
    {
      Core::LinAlg::SerialDenseMatrix derivatives(3, data.num_nodes);
      Core::FE::shape_function_3d_deriv1(derivatives, corner[0], corner[1], corner[2], celltype);
      data.derivatives.push_back(derivatives);
    }
    return data;
  }

  const CornerData& corner_data(const Core::FE::CellType celltype)
  {
    static const CornerData hex8 = make_corner_data(Core::FE::CellType::hex8);
    static const CornerData tet4 = make_corner_data(Core::FE::CellType::tet4);
    switch (celltype)
    {
      case Core::FE::CellType::hex8:
        return hex8;
      case Core::FE::CellType::tet4:
        return tet4;
      default:
        FOUR_C_THROW("The mesh quality is not implemented for cell type {}.",
            Core::FE::cell_type_to_string(celltype));
    }
  }

  //! Jacobian of the mapping from parameter space at a corner
  template <typename T>
  Matrix3<T> corner_jacobian(const Core::LinAlg::SerialDenseMatrix& derivatives,
      std::span<const double> coordinates, const std::vector<T>& displacements)
  {
    Matrix3<T> J{};
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j) J[i][j] = 0.0;

    for (int node = 0; node < derivatives.numCols(); ++node)
      for (int i = 0; i < 3; ++i)
      {
        const T x = coordinates[3 * node + i] + displacements[3 * node + i];
        for (int j = 0; j < 3; ++j) J[i][j] += x * derivatives(j, node);
      }
    return J;
  }

  template <typename T>
  Matrix3<T> multiply(const Matrix3<T>& A, const Matrix3<double>& B)
  {
    Matrix3<T> C{};
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
      {
        C[i][j] = 0.0;
        for (int k = 0; k < 3; ++k) C[i][j] += A[i][k] * B[k][j];
      }
    return C;
  }

  /*!
   * @brief Mean ratio distortion |T|^2 / (3 det(T)^(2/3)) >= 1 of the Jacobian relative to the
   * ideal element
   *
   * @return std::nullopt for inverted corners
   */
  template <typename T>
  std::optional<T> mean_ratio_distortion(const Matrix3<T>& J, const Matrix3<double>& W_inv)
  {
    const Matrix3<T> relative = multiply(J, W_inv);
    const T det = determinant(relative);
    if (!(Sacado::ScalarValue<T>::eval(det) > 0.0)) return std::nullopt;

    T frobenius_squared = 0.0;
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j) frobenius_squared += relative[i][j] * relative[i][j];

    return frobenius_squared / (3.0 * std::pow(det, 2.0 / 3.0));
  }

  template <typename T>
  T potential(const CornerData& data, std::span<const double> reference_coordinates,
      const std::vector<T>& displacements, std::span<const double> target_coordinates,
      const ALE::MeshQualityParameters& parameters)
  {
    const std::vector<double> zero(displacements.size(), 0.0);
    const double weight = 1.0 / static_cast<double>(data.derivatives.size());

    // the target shape and size of each corner: those of the target configuration if given,
    // otherwise the ideal element shape and the size in the reference configuration
    std::vector<Matrix3<double>> target_jacobian_inverse;
    std::vector<double> target_det;
    // maps the displacement at each corner to units of the element size, U J^-1 with the Jacobian
    // U of the ideal element with unit edges and the Jacobian J of the target if given, otherwise
    // of the reference configuration: |U J^-1 u| = |u| / h for any rotated ideal element of size h
    std::vector<Matrix3<double>> motion_scaling;
    for (const auto& derivatives : data.derivatives)
    {
      if (target_coordinates.empty())
      {
        const Matrix3<double> reference_jacobian =
            corner_jacobian(derivatives, reference_coordinates, zero);
        target_jacobian_inverse.push_back(data.ideal_jacobian_inverse);
        target_det.push_back(determinant(reference_jacobian));
        if (target_det.back() > 0.0)
          motion_scaling.push_back(multiply(data.unit_ideal_jacobian, inverse(reference_jacobian)));
      }
      else
      {
        const Matrix3<double> target_jacobian =
            corner_jacobian(derivatives, target_coordinates, zero);
        target_det.push_back(determinant(target_jacobian));
        if (target_det.back() > 0.0)
        {
          target_jacobian_inverse.push_back(inverse(target_jacobian));
          motion_scaling.push_back(
              multiply(data.unit_ideal_jacobian, target_jacobian_inverse.back()));
        }
      }
    }

    // the potential is not defined for elements whose target is inverted
    if (std::ranges::any_of(target_det, [](double det) { return det <= 0.0; }))
      return std::numeric_limits<double>::quiet_NaN();

    const bool corner_per_node = static_cast<int>(data.derivatives.size()) == data.num_nodes;

    T energy = 0.0;
    for (std::size_t c = 0; c < data.derivatives.size(); ++c)
    {
      const Matrix3<T> J =
          corner_jacobian(data.derivatives[c], reference_coordinates, displacements);
      const auto distortion = mean_ratio_distortion(J, target_jacobian_inverse[c]);
      if (!distortion) return std::numeric_limits<double>::quiet_NaN();

      const T volume_ratio = determinant(J) / target_det[c];
      // for growth only, (v - 1)^3 / v^2 for v > 1 has a continuous Hessian at v = 1
      T size = volume_ratio + 1.0 / volume_ratio - 2.0;
      if (parameters.size_growth_only)
      {
        if (Sacado::ScalarValue<T>::eval(volume_ratio) > 1.0)
        {
          const T excess = volume_ratio - 1.0;
          size = excess * excess * excess / (volume_ratio * volume_ratio);
        }
        else
          size = 0.0;
      }

      // (mu^p - 1) / p equals mu - 1 near ideal shapes and grows faster for strong distortions
      T shape = *distortion - 1.0;
      if (parameters.shape_exponent != 1.0)
      {
        using std::pow;
        shape = (pow(*distortion, parameters.shape_exponent) - 1.0) / parameters.shape_exponent;
      }

      energy += weight * (shape + parameters.size_weight * size);
    }

    if (parameters.motion_weight != 0.0)
    {
      // the displacement of each node relative to the element size at its corner (a single
      // corner for elements with a constant Jacobian)
      for (int node = 0; node < data.num_nodes; ++node)
      {
        const Matrix3<double>& scaling = motion_scaling[corner_per_node ? node : 0];
        T squared_motion = 0.0;
        for (int i = 0; i < 3; ++i)
        {
          T component = 0.0;
          for (int j = 0; j < 3; ++j) component += scaling[i][j] * displacements[3 * node + j];
          squared_motion += component * component;
        }
        energy += parameters.motion_weight / data.num_nodes * squared_motion;
      }
    }
    return energy;
  }
}  // namespace

double ALE::mean_ratio_quality(
    const Core::FE::CellType celltype, std::span<const double> coordinates)
{
  const CornerData& data = corner_data(celltype);
  const std::vector<double> zero(coordinates.size(), 0.0);

  double quality = 1.0;
  for (const auto& derivatives : data.derivatives)
  {
    const Matrix3<double> J = corner_jacobian(derivatives, coordinates, zero);
    const auto distortion = mean_ratio_distortion(J, data.ideal_jacobian_inverse);
    if (!distortion) return -1.0;
    quality = std::min(quality, 1.0 / *distortion);
  }
  return quality;
}

double ALE::element_size(const Core::FE::CellType celltype, std::span<const double> coordinates)
{
  const CornerData& data = corner_data(celltype);
  const std::vector<double> zero(coordinates.size(), 0.0);
  double sum = 0.0;
  for (const auto& derivatives : data.derivatives)
    sum += determinant(corner_jacobian(derivatives, coordinates, zero));
  return sum / static_cast<double>(data.derivatives.size()) * data.parameter_space_volume;
}

double ALE::evaluate_mesh_quality_potential(const Core::FE::CellType celltype,
    std::span<const double> reference_coordinates, std::span<const double> displacements,
    const MeshQualityParameters& parameters, Core::LinAlg::SerialDenseVector* residual,
    Core::LinAlg::SerialDenseMatrix* stiffness, std::span<const double> target_coordinates)
{
  const CornerData& data = corner_data(celltype);
  const int num_dofs = 3 * data.num_nodes;
  FOUR_C_ASSERT_ALWAYS(
      static_cast<int>(reference_coordinates.size()) == num_dofs &&
          static_cast<int>(displacements.size()) == num_dofs &&
          (target_coordinates.empty() || static_cast<int>(target_coordinates.size()) == num_dofs),
      "Expected {} nodal coordinates and displacements.", num_dofs);

  // the value alone needs no automatic differentiation
  if (!residual && !stiffness)
  {
    const std::vector<double> u(displacements.begin(), displacements.end());
    return potential(data, reference_coordinates, u, target_coordinates, parameters);
  }

  using FAD = Sacado::Fad::DFad<double>;
  using FAD2 = Sacado::Fad::DFad<FAD>;

  std::vector<FAD2> u(num_dofs);
  for (int i = 0; i < num_dofs; ++i) u[i] = FAD2(num_dofs, i, FAD(num_dofs, i, displacements[i]));

  const FAD2 energy = potential(data, reference_coordinates, u, target_coordinates, parameters);
  const double value = energy.val().val();

  if (residual)
  {
    residual->size(num_dofs);
    for (int i = 0; i < num_dofs; ++i)
      (*residual)(i) =
          std::isnan(value) ? std::numeric_limits<double>::quiet_NaN() : energy.dx(i).val();
  }
  if (stiffness)
  {
    stiffness->shape(num_dofs, num_dofs);
    if (!std::isnan(value))
      for (int i = 0; i < num_dofs; ++i)
        for (int j = 0; j < num_dofs; ++j) (*stiffness)(i, j) = energy.dx(i).dx(j);
  }

  return value;
}

FOUR_C_NAMESPACE_CLOSE
