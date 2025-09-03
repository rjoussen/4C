// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_linalg_utils_densematrix_svd.hpp"

#include "4C_utils_exceptions.hpp"

#include <Teuchos_ScalarTraits.hpp>

FOUR_C_NAMESPACE_OPEN

void Core::LinAlg::svd(const Core::LinAlg::SerialDenseMatrix& A, Core::LinAlg::SerialDenseMatrix& Q,
    Core::LinAlg::SerialDenseMatrix& S, Core::LinAlg::SerialDenseMatrix& VT)
{
  Core::LinAlg::SerialDenseMatrix tmp(A);  // copy, because content of A is destroyed
  const char jobu = 'A';                   // compute and return all M columns of U
  const char jobvt = 'A';                  // compute and return all N rows of V^T
  const int n = tmp.numCols();
  const int m = tmp.numRows();
  std::vector<double> s(std::min(n, m));
  int info;
  int lwork = std::max(3 * std::min(m, n) + std::max(m, n), 5 * std::min(m, n));
  std::vector<double> work(lwork);
  double rwork;

  Teuchos::LAPACK<int, double> lapack;
  lapack.GESVD(jobu, jobvt, m, n, tmp.values(), tmp.stride(), s.data(), Q.values(), Q.stride(),
      VT.values(), VT.stride(), work.data(), lwork, &rwork, &info);

  if (info) FOUR_C_THROW("Lapack's dgesvd returned {}", info);

  // 0 for off-diagonal, otherwise s
  S.put_scalar(0.0);
  for (int i = 0; i < std::min(n, m); ++i)
  {
    S(i, i) = s[i];
  }
}

FOUR_C_NAMESPACE_CLOSE
