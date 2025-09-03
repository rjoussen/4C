// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_LINALG_UTILS_DENSEMATRIX_SVD_HPP
#define FOUR_C_LINALG_UTILS_DENSEMATRIX_SVD_HPP

#include "4C_config.hpp"

#include "4C_linalg_fixedsizematrix.hpp"
#include "4C_linalg_serialdensematrix.hpp"

#include <Teuchos_LAPACK.hpp>

FOUR_C_NAMESPACE_OPEN

namespace Core::LinAlg
{
  /*!
   \brief Compute singular value decomposition (SVD) of a real M-by-N matrix A
   A = U * SIGMA * transpose(V)

   \param A (in/out):    Matrix to be decomposed
   \param U (in/out):    M-by-M orthogonal matrix
   \param SIGMA (in/out):M-by-N matrix which is zero except for its min(m,n) diagonal elements
   \param Vt (in/out):   V is a N-by-N orthogonal matrix, actually returned is V^T
   */
  void svd(const Core::LinAlg::SerialDenseMatrix& A, Core::LinAlg::SerialDenseMatrix& Q,
      Core::LinAlg::SerialDenseMatrix& SIGMA, Core::LinAlg::SerialDenseMatrix& Vt);

  /*!
   \brief Singular value decomposition (SVD) of a real M-by-N matrix in fixed
   size format

   A = Q * S * VT

   \tparam row Number of rows
   \tparam col Number of columns

   \param A (in):        M-by-N matrix to be decomposed
   \param Q (out):       M-by-M orthogonal matrix
   \param S (out):       M-by-N matrix which is zero except for its min(m,n) diagonal elements
   \param VT (out):      N-by-N orthogonal matrix (transpose of V)
   */
  template <unsigned int rows, unsigned int cols>
  void svd(const Core::LinAlg::Matrix<rows, cols>& A, Core::LinAlg::Matrix<rows, rows>& Q,
      Core::LinAlg::Matrix<rows, cols>& S, Core::LinAlg::Matrix<cols, cols>& VT)
  {
    Matrix<rows, cols> tmp(A.data(), false);  // copy, because content of matrix is destroyed
    const char jobu = 'A';                    // compute and return all M columns of U
    const char jobvt = 'A';                   // compute and return all N rows of V^T
    std::vector<double> s(std::min(rows, cols));
    int info;
    int lwork = std::max(3 * std::min(rows, cols) + std::max(rows, cols), 5 * std::min(rows, cols));
    std::vector<double> work(lwork);
    double rwork;

    Teuchos::LAPACK<int, double> lapack;
    lapack.GESVD(jobu, jobvt, rows, cols, tmp.data(), tmp.m(), s.data(), Q.data(), Q.m(), VT.data(),
        VT.m(), work.data(), lwork, &rwork, &info);

    if (info) FOUR_C_THROW("Lapack's dgesvd returned {}", info);

    for (unsigned int i = 0; i < std::min(rows, cols); ++i)
    {
      for (unsigned int j = 0; j < std::min(rows, cols); ++j)
      {
        S(i, j) = (i == j) * s[i];  // 0 for off-diagonal, otherwise s
      }
    }
    return;
  }

}  // namespace Core::LinAlg

FOUR_C_NAMESPACE_CLOSE

#endif
