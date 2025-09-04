// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_utils_cubic_spline_interpolation.hpp"

#include "4C_utils_exceptions.hpp"

#include <Teuchos_SerialDenseSolver.hpp>

#include <utility>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
Core::Utils::CubicSplineInterpolation::CubicSplineInterpolation(
    std::vector<double> x, std::vector<double> y)
    : a_(std::move(y)), x_(std::move(x))
// zeroth-order coefficients a_ equal the function values y
{
  // safety checks
  if (x_.size() != a_.size())
    FOUR_C_THROW("Length of vectors provided to cubic spline interpolation must match!");
  if (not std::is_sorted(x_.begin(), x_.end()))
    FOUR_C_THROW("Data points must be sorted in ascending order!");

  // temp variables for solution of the linear system
  const int N = static_cast<int>(x_.size());
  Core::LinAlg::SerialDenseMatrix A(N, N);
  Core::LinAlg::SerialDenseVector c(N), b(N);

  build_matrix_and_rhs(N, A, b);
  solve_linear_system(A, c, b);
  setup_internal_vectors(c);
}

/*----------------------------------------------------------------------*/
void Core::Utils::CubicSplineInterpolation::build_matrix_and_rhs(
    const int N, Core::LinAlg::SerialDenseMatrix& A, Core::LinAlg::SerialDenseVector& b) const
{
  // fill everything except the boundary condition lines
  for (int i = 1; i < N - 1; ++i)
  {
    // a_ = y, see initialization list in constructor
    const double delta_y_p = a_[i + 1] - a_[i];
    const double delta_y_m = a_[i] - a_[i - 1];
    const double delta_x_p = x_[i + 1] - x_[i];
    const double delta_x_m = x_[i] - x_[i - 1];
    const double delta_x_pm = x_[i + 1] - x_[i - 1];

    const double alpha_i = delta_x_m / delta_x_pm;
    const double beta_i = delta_x_p / delta_x_pm;
    const double gamma_i = 3.0 / delta_x_pm * ((delta_y_p / delta_x_p) - (delta_y_m / delta_x_m));

    A(i, i - 1) = alpha_i;
    A(i, i) = 2.0;
    A(i, i + 1) = beta_i;
    b(i) = gamma_i;
  }

  // set values due to natural boundary conditions
  A(0, 0) = A(N - 1, N - 1) = 2.0;
  b(0) = b(N - 1) = 0.0;
}

/*----------------------------------------------------------------------*/
double Core::Utils::CubicSplineInterpolation::evaluate(const double x) const
{
  // safety check
  if (x < x_.front() or x > x_.back())
    FOUR_C_THROW("Sampling point x = {} lies outside sampling point range!", x);

  auto greater_equal_x = [x](const double val) { return val >= x; };
  auto right_position = std::find_if(x_.begin(), x_.end(), greater_equal_x);
  FOUR_C_ASSERT(right_position != x_.begin(), "Internal error.");
  // the left side of the sought interval is found by deleting 1
  const auto left_position = std::distance(x_.begin(), right_position) - 1;

  const double delta_x = x - x_[left_position];
  return a_[left_position] + b_[left_position] * delta_x + c_[left_position] * delta_x * delta_x +
         d_[left_position] * delta_x * delta_x * delta_x;
}

/*----------------------------------------------------------------------*/
double Core::Utils::CubicSplineInterpolation::evaluate_derivative(
    const double x, const int deriv_order) const
{
  // safety check
  if (x < x_.front() or x > x_.back())
    FOUR_C_THROW("Sampling point x = {} lies outside sampling point range!", x);

  auto greater_equal_x = [x](const double val) { return val >= x; };
  auto right_position = std::find_if(x_.begin(), x_.end(), greater_equal_x);
  FOUR_C_ASSERT(right_position != x_.begin(), "Internal error.");
  // the left side of the sought interval is found by deleting 1
  const auto left_position = std::distance(x_.begin(), right_position) - 1;

  switch (deriv_order)
  {
    case 1:
    {
      const double delta_x = x - x_[left_position];
      return b_[left_position] + 2.0 * c_[left_position] * delta_x +
             3.0 * d_[left_position] * delta_x * delta_x;
    }
    case 2:
    {
      const double delta_x = x - x_[left_position];
      return 2.0 * c_[left_position] + 6.0 * d_[left_position] * delta_x;
    }
  }

  FOUR_C_THROW(
      "Evaluation of {} derivative is not implemented for CubicSplineInterpolation", deriv_order);
}

/*----------------------------------------------------------------------*/
void Core::Utils::CubicSplineInterpolation::solve_linear_system(Core::LinAlg::SerialDenseMatrix& A,
    Core::LinAlg::SerialDenseVector& c, Core::LinAlg::SerialDenseVector& b) const
{
  // solve for third-order coefficients for cubic spline interpolation
  using ordinalType = Core::LinAlg::SerialDenseMatrix::ordinalType;
  using scalarType = Core::LinAlg::SerialDenseMatrix::scalarType;
  Teuchos::SerialDenseSolver<ordinalType, scalarType> solver;
  solver.setMatrix(Teuchos::rcpFromRef(A.base()));
  solver.setVectors(Teuchos::rcpFromRef(c), Teuchos::rcpFromRef(b));
  solver.factorWithEquilibration(true);
  solver.solveToRefinedSolution(true);
  if (solver.factor() or solver.solve())
    FOUR_C_THROW("Solution of linear system of equations failed!");
}

/*----------------------------------------------------------------------*/
void Core::Utils::CubicSplineInterpolation::setup_internal_vectors(
    const Core::LinAlg::SerialDenseVector& c)
{
  const std::size_t system_size = x_.size();

  // resize of the internal member vectors
  c_.resize(system_size, 0.0);
  b_.resize(system_size - 1, 0.0), d_.resize(system_size - 1, 0.0);

  for (auto i = 0; i < static_cast<int>(system_size); ++i) c_[i] = c(i);

  for (std::size_t i = 0; i < system_size - 1; ++i)
  {
    // a_ = y, see initialization list in constructor
    const double delta_y_p = a_[i + 1] - a_[i];
    const double delta_x_p = x_[i + 1] - x_[i];

    b_[i] = delta_y_p / delta_x_p - 1.0 / 3.0 * delta_x_p * (2.0 * c_[i] + c_[i + 1]);
    d_[i] = (c_[i + 1] - c_[i]) / (3.0 * delta_x_p);
  }
}
FOUR_C_NAMESPACE_CLOSE
