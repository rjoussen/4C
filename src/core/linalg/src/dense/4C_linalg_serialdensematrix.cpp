// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_linalg_serialdensematrix.hpp"

#include "4C_utils_exceptions.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Core::LinAlg
{
  SerialDenseMatrix::SerialDenseMatrix(int rows, int cols) : mat_(rows, cols) {}
  SerialDenseMatrix::SerialDenseMatrix(int rows, int cols, bool zeroOut) : mat_(rows, cols, zeroOut)
  {
  }
  SerialDenseMatrix::SerialDenseMatrix(const Base& src, Teuchos::ETransp trans) : mat_(src, trans)
  {
  }
  SerialDenseMatrix::SerialDenseMatrix(
      Teuchos::DataAccess cv, double* values, int stride, int rows, int cols)
      : mat_(cv, values, stride, rows, cols)
  {
  }

  int SerialDenseMatrix::num_rows() const { return mat_.numRows(); }
  int SerialDenseMatrix::num_cols() const { return mat_.numCols(); }
  int SerialDenseMatrix::numRows() const { return mat_.numRows(); }
  int SerialDenseMatrix::numCols() const { return mat_.numCols(); }
  int SerialDenseMatrix::stride() const { return mat_.stride(); }
  double SerialDenseMatrix::normInf() const { return mat_.normInf(); }
  double SerialDenseMatrix::normOne() const { return mat_.normOne(); }

  double* SerialDenseMatrix::values() const { return mat_.values(); }

  void SerialDenseMatrix::scale(double alpha)
  {
    [[maybe_unused]] int err = mat_.scale(alpha);
    FOUR_C_ASSERT(err == 0, "Internal error: Scalar multiplication did not work.");
  }
  void SerialDenseMatrix::put_scalar(double val)
  {
    [[maybe_unused]] int err = mat_.putScalar(val);
    FOUR_C_ASSERT(err == 0, "Internal error: Matrix update with scalar did not work.");
  }
  void SerialDenseMatrix::putScalar(double val)
  {
    [[maybe_unused]] int err = mat_.putScalar(val);
    FOUR_C_ASSERT(err == 0, "Internal error: Matrix update with scalar did not work.");
  }
  void SerialDenseMatrix::shape(int rows, int cols)
  {
    [[maybe_unused]] int err = mat_.shape(rows, cols);
    FOUR_C_ASSERT(
        err == 0, "Internal error: Matrix shape update and reinitialization to zero did not work.");
  }
  void SerialDenseMatrix::reshape(int rows, int cols)
  {
    [[maybe_unused]] int err = mat_.reshape(rows, cols);
    FOUR_C_ASSERT(err == 0, "Internal error: Matrix reshaping did not work.");
  }
  void SerialDenseMatrix::assign(const SerialDenseMatrix& source) { mat_.assign(source.base()); }

  SerialDenseMatrix& SerialDenseMatrix::operator+=(const SerialDenseMatrix& other)
  {
    mat_ += other.mat_;
    return *this;
  }

  SerialDenseMatrix::Base& SerialDenseMatrix::base() { return mat_; }
  const SerialDenseMatrix::Base& SerialDenseMatrix::base() const { return mat_; }


  void update(double alpha, const Core::LinAlg::SerialDenseMatrix& A, double beta,
      Core::LinAlg::SerialDenseMatrix& B)
  {
    B.scale(beta);
    Core::LinAlg::SerialDenseMatrix Acopy(A.base());
    Acopy.scale(alpha);
    B += Acopy;
  }

  long double det_long(const SerialDenseMatrix& matrix)
  {
    int n = matrix.num_cols();

    if (n == 1)
    {
      return matrix(0, 0);
    }
    else if (n == 2)
    {
      return (long double)matrix(0, 0) * (long double)matrix(1, 1) -
             (long double)matrix(0, 1) * (long double)matrix(1, 0);
    }
    else if (n > 2)
    {
      long double out_det = 0;
      int sign = 1;

      for (int i_col = 0; i_col < n; i_col++)
      {
        SerialDenseMatrix temp(n - 1, n - 1);

        for (int c_col = 0; c_col < i_col; c_col++)
          for (int row = 1; row < n; row++) temp(row - 1, c_col) = matrix(row, c_col);

        for (int c_col = i_col + 1; c_col < n; c_col++)
          for (int row = 1; row < n; row++) temp(row - 1, c_col - 1) = matrix(row, c_col);

        out_det += (long double)sign * (long double)matrix(0, i_col) * det_long(temp);
        sign *= -1;
      }
      return out_det;
    }
    else
      return 0;
  }

  void zero(SerialDenseMatrix& mat, int Length)
  {
    int cnt = 0;
    for (int j = 0; j < mat.num_cols(); ++j)
    {
      for (int i = 0; i < mat.num_rows(); ++i)
      {
        if (cnt++ < Length)
          mat(i, j) = 0.0;
        else
          return;
      }
    }
  }

  void copy(const double* vec, SerialDenseMatrix& mat)
  {
    int cnt = 0;
    for (int j = 0; j < mat.num_cols(); ++j)
      for (int i = 0; i < mat.num_rows(); ++i) mat(i, j) = vec[cnt++];
  }

}  // namespace Core::LinAlg

FOUR_C_NAMESPACE_CLOSE
