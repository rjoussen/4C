// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_MAT_HISTORY_TRANSFER_HPP
#define FOUR_C_MAT_HISTORY_TRANSFER_HPP

#include "4C_config.hpp"

#include "4C_linalg_fixedsizematrix.hpp"

#include <string>
#include <vector>

FOUR_C_NAMESPACE_OPEN

namespace Mat
{
  /*!
   * @brief Kind of a history quantity, which determines how the quantity is interpolated when
   * history is transferred between meshes (e.g., after moving the nodes of a mesh).
   *
   * Each kind defines the properties that are preserved by the interpolation.
   */
  enum class HistoryTransferKind
  {
    //! scalar without constraints (1 value)
    scalar,
    //! non-negative scalar, e.g., an accumulated plastic strain (1 value)
    nonnegative_scalar,
    //! invertible second-order tensor T = R U with positive determinant (9 values, row-major).
    //! The rotations and the logarithms of the stretches U are interpolated separately, such that
    //! U stays positive definite, the determinant is preserved if all data has the same
    //! determinant (e.g., det = 1 for the inverse plastic deformation gradient in isochoric
    //! plasticity), and components that vanish in all data stay zero (e.g., in plane strain).
    invertible_tensor,
    //! discrete data, which is taken from the closest data point (1 value)
    nearest,
  };

  /*!
   * @brief A named history quantity of a given kind
   */
  struct HistoryTransferQuantity
  {
    std::string name;
    HistoryTransferKind kind;
  };

  /*!
   * @brief Number of values of a history quantity of the given kind
   */
  int history_transfer_size(HistoryTransferKind kind);

  /*!
   * @brief Total number of values of the given history quantities
   */
  int history_transfer_size(const std::vector<HistoryTransferQuantity>& quantities);

  /*!
   * @brief Interpolate history data given at a set of points to a target point
   *
   * The data is approximated with moving least squares and exponential weights exp(-c d^2) of the
   * distances d to the target (Satheesh et al., A novel mesh regularization approach based on
   * finite element distortion potentials, Comput Methods Appl Mech Eng, 2024). Scalars are
   * approximated with a trilinear basis if at least 8 points are given (otherwise with a weighted
   * average), tensors by weighted averages of the rotation vectors and of the logarithms of the
   * stretches (log-Euclidean mean). Unlike the R-LOG scheme of Satheesh et al., which interpolates
   * the eigenvalues and eigenvectors of the stretches, this does not need to match eigenvectors of
   * different points, which turns them out of their plane where the order of the eigenvalues
   * changes between the points. If the target coincides with one of the points, the data of this point is
   * returned unchanged.
   *
   * @param quantities Description of the history data. The data of each point contains the values
   * of all quantities in this order.
   * @param data History data at the points
   * @param locations Locations of the points (e.g., Gauss points of an element patch)
   * @param target Location of the target point
   * @param weight_decay Decay c of the weights in the units of the locations
   * @return Interpolated history data at the target point
   */
  std::vector<double> interpolate_history(const std::vector<HistoryTransferQuantity>& quantities,
      const std::vector<std::vector<double>>& data,
      const std::vector<Core::LinAlg::Matrix<3, 1>>& locations,
      const Core::LinAlg::Matrix<3, 1>& target, double weight_decay);

}  // namespace Mat

FOUR_C_NAMESPACE_CLOSE

#endif
