// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_linalg_fixedsizematrix_tensor_products.hpp"

#include "4C_linalg_symmetric_tensor.hpp"
#include "4C_linalg_tensor.hpp"
#include "4C_linalg_tensor_symmetric_einstein.hpp"

#include <Sacado.hpp>

using FAD = Sacado::Fad::DFad<double>;


FOUR_C_NAMESPACE_OPEN

void Core::LinAlg::FourTensorOperations::add_elasticity_tensor_product(
    Core::LinAlg::Matrix<6, 6>& C, const double scalar_AB, const Core::LinAlg::Matrix<3, 3>& A,
    const Core::LinAlg::Matrix<3, 3>& B, const double scalar_this)
{
  // everything in Voigt-Notation
  Core::LinAlg::Matrix<6, 1> A_voigt;
  Core::LinAlg::Matrix<6, 1> B_voigt;

  A_voigt(0, 0) = A(0, 0);
  A_voigt(1, 0) = A(1, 1);
  A_voigt(2, 0) = A(2, 2);
  /* Voigts vector notation on strain entities usually implies 2 times ()12 ()23 ()13
   * however, this is not the case here to arrive at the consistent elasticity */
  A_voigt(3, 0) = A(1, 0);
  A_voigt(4, 0) = A(2, 1);
  A_voigt(5, 0) = A(2, 0);

  B_voigt(0, 0) = B(0, 0);
  B_voigt(1, 0) = B(1, 1);
  B_voigt(2, 0) = B(2, 2);
  B_voigt(3, 0) = B(1, 0);
  B_voigt(4, 0) = B(2, 1);
  B_voigt(5, 0) = B(2, 0);

  C.multiply_nt(scalar_AB, A_voigt, B_voigt, scalar_this);
}

void Core::LinAlg::FourTensorOperations::add_symmetric_elasticity_tensor_product(
    Core::LinAlg::Matrix<6, 6>& C, const double scalar_AB, const Core::LinAlg::Matrix<3, 3>& A,
    const Core::LinAlg::Matrix<3, 3>& B, const double scalar_this)
{
#ifdef FOUR_C_ENABLE_ASSERTIONS
  // check sizes
  if (A.num_rows() != A.num_cols() || B.num_rows() != B.num_cols() || A.num_rows() != 3 ||
      B.num_rows() != 3)
  {
    FOUR_C_THROW("2nd order tensors must be 3 by 3");
  }
  if (C.num_rows() != C.num_cols() || C.num_rows() != 6)
    FOUR_C_THROW("4th order tensor must be 6 by 6");
#endif

  // everything in Voigt-Notation
  Core::LinAlg::Matrix<6, 1> A_voigt;
  Core::LinAlg::Matrix<6, 1> B_voigt;

  A_voigt(0, 0) = A(0, 0);
  A_voigt(1, 0) = A(1, 1);
  A_voigt(2, 0) = A(2, 2);
  /* Voigts vector notation on strain entities usually implies 2 times ()12 ()23 ()13
   * however, this is not the case here to arrive at the consistent elasticity */
  A_voigt(3, 0) = A(1, 0);
  A_voigt(4, 0) = A(2, 1);
  A_voigt(5, 0) = A(2, 0);

  B_voigt(0, 0) = B(0, 0);
  B_voigt(1, 0) = B(1, 1);
  B_voigt(2, 0) = B(2, 2);
  B_voigt(3, 0) = B(1, 0);
  B_voigt(4, 0) = B(2, 1);
  B_voigt(5, 0) = B(2, 0);

  C.multiply_nt(scalar_AB, A_voigt, B_voigt, scalar_this);
  C.multiply_nt(scalar_AB, B_voigt, A_voigt, 1.0);
}

void Core::LinAlg::FourTensorOperations::add_kronecker_tensor_product(Core::LinAlg::Matrix<6, 6>& C,
    const double scalar_AB, const Core::LinAlg::Matrix<3, 3>& A,
    const Core::LinAlg::Matrix<3, 3>& B, const double scalar_this)
{
  const double scalar_AB_half = scalar_AB * 0.5;
  C(0, 0) = scalar_this * C(0, 0) + scalar_AB * (A(0, 0) * B(0, 0));  // C1111
  C(0, 1) = scalar_this * C(0, 1) + scalar_AB * (A(0, 1) * B(0, 1));  // C1122
  C(0, 2) = scalar_this * C(0, 2) + scalar_AB * (A(0, 2) * B(0, 2));  // C1133
  C(0, 3) =
      scalar_this * C(0, 3) + scalar_AB_half * (A(0, 0) * B(0, 1) + A(0, 1) * B(0, 0));  // C1112
  C(0, 4) =
      scalar_this * C(0, 4) + scalar_AB_half * (A(0, 1) * B(0, 2) + A(0, 2) * B(0, 1));  // C1123
  C(0, 5) =
      scalar_this * C(0, 5) + scalar_AB_half * (A(0, 0) * B(0, 2) + A(0, 2) * B(0, 0));  // C1113

  C(1, 0) = scalar_this * C(1, 0) + scalar_AB * (A(1, 0) * B(1, 0));  // C2211
  C(1, 1) = scalar_this * C(1, 1) + scalar_AB * (A(1, 1) * B(1, 1));  // C2222
  C(1, 2) = scalar_this * C(1, 2) + scalar_AB * (A(1, 2) * B(1, 2));  // C2233
  C(1, 3) =
      scalar_this * C(1, 3) + scalar_AB_half * (A(1, 0) * B(1, 1) + A(1, 1) * B(1, 0));  // C2212
  C(1, 4) =
      scalar_this * C(1, 4) + scalar_AB_half * (A(1, 1) * B(1, 2) + A(1, 2) * B(1, 1));  // C2223
  C(1, 5) =
      scalar_this * C(1, 5) + scalar_AB_half * (A(1, 0) * B(1, 2) + A(1, 2) * B(1, 0));  // C2213

  C(2, 0) = scalar_this * C(2, 0) + scalar_AB * (A(2, 0) * B(2, 0));  // C3311
  C(2, 1) = scalar_this * C(2, 1) + scalar_AB * (A(2, 1) * B(2, 1));  // C3322
  C(2, 2) = scalar_this * C(2, 2) + scalar_AB * (A(2, 2) * B(2, 2));  // C3333
  C(2, 3) =
      scalar_this * C(2, 3) + scalar_AB_half * (A(2, 0) * B(2, 1) + A(2, 1) * B(2, 0));  // C3312
  C(2, 4) =
      scalar_this * C(2, 4) + scalar_AB_half * (A(2, 1) * B(2, 2) + A(2, 2) * B(2, 1));  // C3323
  C(2, 5) =
      scalar_this * C(2, 5) + scalar_AB_half * (A(2, 0) * B(2, 2) + A(2, 2) * B(2, 0));  // C3313

  C(3, 0) = scalar_this * C(3, 0) + scalar_AB * (A(0, 0) * B(1, 0));  // C1211
  C(3, 1) = scalar_this * C(3, 1) + scalar_AB * (A(0, 1) * B(1, 1));  // C1222
  C(3, 2) = scalar_this * C(3, 2) + scalar_AB * (A(0, 2) * B(1, 2));  // C1233
  C(3, 3) =
      scalar_this * C(3, 3) + scalar_AB_half * (A(0, 0) * B(1, 1) + A(0, 1) * B(1, 0));  // C1212
  C(3, 4) =
      scalar_this * C(3, 4) + scalar_AB_half * (A(0, 1) * B(1, 2) + A(0, 2) * B(1, 1));  // C1223
  C(3, 5) =
      scalar_this * C(3, 5) + scalar_AB_half * (A(0, 0) * B(1, 2) + A(0, 2) * B(1, 0));  // C1213

  C(4, 0) = scalar_this * C(4, 0) + scalar_AB * (A(1, 0) * B(2, 0));  // C2311
  C(4, 1) = scalar_this * C(4, 1) + scalar_AB * (A(1, 1) * B(2, 1));  // C2322
  C(4, 2) = scalar_this * C(4, 2) + scalar_AB * (A(1, 2) * B(2, 2));  // C2333
  C(4, 3) =
      scalar_this * C(4, 3) + scalar_AB_half * (A(1, 0) * B(2, 1) + A(1, 1) * B(2, 0));  // C2312
  C(4, 4) =
      scalar_this * C(4, 4) + scalar_AB_half * (A(1, 1) * B(2, 2) + A(1, 2) * B(2, 1));  // C2323
  C(4, 5) =
      scalar_this * C(4, 5) + scalar_AB_half * (A(1, 0) * B(2, 2) + A(1, 2) * B(2, 0));  // C2313

  C(5, 0) = scalar_this * C(5, 0) + scalar_AB * (A(0, 0) * B(2, 0));  // C1311
  C(5, 1) = scalar_this * C(5, 1) + scalar_AB * (A(0, 1) * B(2, 1));  // C1322
  C(5, 2) = scalar_this * C(5, 2) + scalar_AB * (A(0, 2) * B(2, 2));  // C1333
  C(5, 3) =
      scalar_this * C(5, 3) + scalar_AB_half * (A(0, 0) * B(2, 1) + A(0, 1) * B(2, 0));  // C1312
  C(5, 4) =
      scalar_this * C(5, 4) + scalar_AB_half * (A(0, 1) * B(2, 2) + A(0, 2) * B(2, 1));  // C1323
  C(5, 5) =
      scalar_this * C(5, 5) + scalar_AB_half * (A(0, 0) * B(2, 2) + A(0, 2) * B(2, 0));  // C1313
}


void Core::LinAlg::FourTensorOperations::add_kronecker_tensor_product(Core::LinAlg::Matrix<6, 9>& C,
    const double scalar_AB, const Core::LinAlg::Matrix<3, 3>& A,
    const Core::LinAlg::Matrix<3, 3>& B, const double scalar_this)
{
  const double scalar_AB_half = scalar_AB * 0.5;
  C(0, 0) = scalar_this * C(0, 0) + scalar_AB * (A(0, 0) * B(0, 0));  // C1111
  C(0, 1) = scalar_this * C(0, 1) + scalar_AB * (A(0, 1) * B(0, 1));  // C1122
  C(0, 2) = scalar_this * C(0, 2) + scalar_AB * (A(0, 2) * B(0, 2));  // C1133
  C(0, 3) =
      scalar_this * C(0, 3) + scalar_AB_half * (A(0, 0) * B(0, 1) + A(0, 1) * B(0, 0));  // C1112
  C(0, 4) =
      scalar_this * C(0, 4) + scalar_AB_half * (A(0, 1) * B(0, 2) + A(0, 2) * B(0, 1));  // C1123
  C(0, 5) =
      scalar_this * C(0, 5) + scalar_AB_half * (A(0, 0) * B(0, 2) + A(0, 2) * B(0, 0));  // C1113
  C(0, 6) =
      scalar_this * C(0, 6) + scalar_AB_half * (A(0, 1) * B(0, 0) + A(0, 0) * B(0, 1));  // C1121
  C(0, 7) =
      scalar_this * C(0, 7) + scalar_AB_half * (A(0, 2) * B(0, 1) + A(0, 1) * B(0, 2));  // C1132
  C(0, 8) =
      scalar_this * C(0, 8) + scalar_AB_half * (A(0, 2) * B(0, 0) + A(0, 0) * B(0, 2));  // C1131


  C(1, 0) = scalar_this * C(1, 0) + scalar_AB * (A(1, 0) * B(1, 0));  // C2211
  C(1, 1) = scalar_this * C(1, 1) + scalar_AB * (A(1, 1) * B(1, 1));  // C2222
  C(1, 2) = scalar_this * C(1, 2) + scalar_AB * (A(1, 2) * B(1, 2));  // C2233
  C(1, 3) =
      scalar_this * C(1, 3) + scalar_AB_half * (A(1, 0) * B(1, 1) + A(1, 1) * B(1, 0));  // C2212
  C(1, 4) =
      scalar_this * C(1, 4) + scalar_AB_half * (A(1, 1) * B(1, 2) + A(1, 2) * B(1, 1));  // C2223
  C(1, 5) =
      scalar_this * C(1, 5) + scalar_AB_half * (A(1, 0) * B(1, 2) + A(1, 2) * B(1, 0));  // C2213
  C(1, 6) =
      scalar_this * C(1, 6) + scalar_AB_half * (A(1, 1) * B(1, 0) + A(1, 0) * B(1, 1));  // C2221
  C(1, 7) =
      scalar_this * C(1, 7) + scalar_AB_half * (A(1, 2) * B(1, 1) + A(1, 1) * B(1, 2));  // C2232
  C(1, 8) =
      scalar_this * C(1, 8) + scalar_AB_half * (A(1, 2) * B(1, 0) + A(1, 0) * B(1, 2));  // C2231



  C(2, 0) = scalar_this * C(2, 0) + scalar_AB * (A(2, 0) * B(2, 0));  // C3311
  C(2, 1) = scalar_this * C(2, 1) + scalar_AB * (A(2, 1) * B(2, 1));  // C3322
  C(2, 2) = scalar_this * C(2, 2) + scalar_AB * (A(2, 2) * B(2, 2));  // C3333
  C(2, 3) =
      scalar_this * C(2, 3) + scalar_AB_half * (A(2, 0) * B(2, 1) + A(2, 1) * B(2, 0));  // C3312
  C(2, 4) =
      scalar_this * C(2, 4) + scalar_AB_half * (A(2, 1) * B(2, 2) + A(2, 2) * B(2, 1));  // C3323
  C(2, 5) =
      scalar_this * C(2, 5) + scalar_AB_half * (A(2, 0) * B(2, 2) + A(2, 2) * B(2, 0));  // C3313
  C(2, 6) =
      scalar_this * C(2, 6) + scalar_AB_half * (A(2, 0) * B(2, 1) + A(2, 1) * B(2, 0));  // C3321
  C(2, 7) =
      scalar_this * C(2, 7) + scalar_AB_half * (A(2, 2) * B(2, 1) + A(2, 1) * B(2, 2));  // C3332
  C(2, 8) =
      scalar_this * C(2, 8) + scalar_AB_half * (A(2, 2) * B(2, 0) + A(2, 0) * B(2, 2));  // C3331



  C(3, 0) = scalar_this * C(3, 0) + scalar_AB * (A(0, 0) * B(1, 0));  // C1211
  C(3, 1) = scalar_this * C(3, 1) + scalar_AB * (A(0, 1) * B(1, 1));  // C1222
  C(3, 2) = scalar_this * C(3, 2) + scalar_AB * (A(0, 2) * B(1, 2));  // C1233
  C(3, 3) =
      scalar_this * C(3, 3) + scalar_AB_half * (A(0, 0) * B(1, 1) + A(0, 1) * B(1, 0));  // C1212
  C(3, 4) =
      scalar_this * C(3, 4) + scalar_AB_half * (A(0, 1) * B(1, 2) + A(0, 2) * B(1, 1));  // C1223
  C(3, 5) =
      scalar_this * C(3, 5) + scalar_AB_half * (A(0, 0) * B(1, 2) + A(0, 2) * B(1, 0));  // C1213
  C(3, 6) =
      scalar_this * C(3, 6) + scalar_AB_half * (A(0, 1) * B(1, 0) + A(0, 0) * B(1, 1));  // C1221
  C(3, 7) =
      scalar_this * C(3, 7) + scalar_AB_half * (A(0, 2) * B(1, 1) + A(0, 1) * B(1, 2));  // C1232
  C(3, 8) =
      scalar_this * C(3, 8) + scalar_AB_half * (A(0, 2) * B(1, 0) + A(0, 0) * B(1, 2));  // C1231



  C(4, 0) = scalar_this * C(4, 0) + scalar_AB * (A(1, 0) * B(2, 0));  // C2311
  C(4, 1) = scalar_this * C(4, 1) + scalar_AB * (A(1, 1) * B(2, 1));  // C2322
  C(4, 2) = scalar_this * C(4, 2) + scalar_AB * (A(1, 2) * B(2, 2));  // C2333
  C(4, 3) =
      scalar_this * C(4, 3) + scalar_AB_half * (A(1, 0) * B(2, 1) + A(1, 1) * B(2, 0));  // C2312
  C(4, 4) =
      scalar_this * C(4, 4) + scalar_AB_half * (A(1, 1) * B(2, 2) + A(1, 2) * B(2, 1));  // C2323
  C(4, 5) =
      scalar_this * C(4, 5) + scalar_AB_half * (A(1, 0) * B(2, 2) + A(1, 2) * B(2, 0));  // C2313
  C(4, 6) =
      scalar_this * C(4, 6) + scalar_AB_half * (A(1, 1) * B(2, 0) + A(1, 0) * B(2, 1));  // C2321
  C(4, 7) =
      scalar_this * C(4, 7) + scalar_AB_half * (A(1, 2) * B(2, 1) + A(1, 1) * B(2, 2));  // C2332
  C(4, 8) =
      scalar_this * C(4, 8) + scalar_AB_half * (A(1, 2) * B(2, 0) + A(1, 0) * B(2, 2));  // C2331



  C(5, 0) = scalar_this * C(5, 0) + scalar_AB * (A(0, 0) * B(2, 0));  // C1311
  C(5, 1) = scalar_this * C(5, 1) + scalar_AB * (A(0, 1) * B(2, 1));  // C1322
  C(5, 2) = scalar_this * C(5, 2) + scalar_AB * (A(0, 2) * B(2, 2));  // C1333
  C(5, 3) =
      scalar_this * C(5, 3) + scalar_AB_half * (A(0, 0) * B(2, 1) + A(0, 1) * B(2, 0));  // C1312
  C(5, 4) =
      scalar_this * C(5, 4) + scalar_AB_half * (A(0, 1) * B(2, 2) + A(0, 2) * B(2, 1));  // C1323
  C(5, 5) =
      scalar_this * C(5, 5) + scalar_AB_half * (A(0, 0) * B(2, 2) + A(0, 2) * B(2, 0));  // C1313
  C(5, 5) =
      scalar_this * C(5, 5) + scalar_AB_half * (A(0, 0) * B(2, 2) + A(0, 2) * B(2, 0));  // C1313
  C(5, 6) =
      scalar_this * C(5, 6) + scalar_AB_half * (A(0, 1) * B(2, 0) + A(0, 0) * B(2, 1));  // C1321
  C(5, 7) =
      scalar_this * C(5, 7) + scalar_AB_half * (A(0, 2) * B(2, 1) + A(0, 1) * B(2, 2));  // C1332
  C(5, 8) =
      scalar_this * C(5, 7) + scalar_AB_half * (A(0, 2) * B(2, 0) + A(0, 0) * B(2, 2));  // C1331
}

template <typename T>
void Core::LinAlg::FourTensorOperations::add_holzapfel_product(
    Core::LinAlg::Matrix<6, 6, T>& cmat, const Core::LinAlg::Matrix<6, 1, T>& invc, const T scalar)
{
#ifdef FOUR_C_ENABLE_ASSERTIONS
  if (cmat.num_rows() != 6 or cmat.num_cols() != 6 or invc.num_rows() != 6)
    FOUR_C_THROW("Wrong dimensions in function add_holzapfel_product");
#endif

  // and the 'boeppel-product' for the expression d(invc)/dc (see Holzapfel p. 254)
  cmat(0, 0) += scalar * invc(0) * invc(0);
  cmat(0, 1) += scalar * invc(3) * invc(3);
  cmat(0, 2) += scalar * invc(5) * invc(5);
  cmat(0, 3) += scalar * invc(0) * invc(3);
  cmat(0, 4) += scalar * invc(3) * invc(5);
  cmat(0, 5) += scalar * invc(0) * invc(5);

  cmat(1, 0) += scalar * invc(3) * invc(3);
  cmat(1, 1) += scalar * invc(1) * invc(1);
  cmat(1, 2) += scalar * invc(4) * invc(4);
  cmat(1, 3) += scalar * invc(3) * invc(1);
  cmat(1, 4) += scalar * invc(1) * invc(4);
  cmat(1, 5) += scalar * invc(3) * invc(4);

  cmat(2, 0) += scalar * invc(5) * invc(5);
  cmat(2, 1) += scalar * invc(4) * invc(4);
  cmat(2, 2) += scalar * invc(2) * invc(2);
  cmat(2, 3) += scalar * invc(5) * invc(4);
  cmat(2, 4) += scalar * invc(4) * invc(2);
  cmat(2, 5) += scalar * invc(5) * invc(2);

  cmat(3, 0) += scalar * invc(0) * invc(3);
  cmat(3, 1) += scalar * invc(3) * invc(1);
  cmat(3, 2) += scalar * invc(5) * invc(4);
  cmat(3, 3) += scalar * 0.5 * (invc(0) * invc(1) + invc(3) * invc(3));
  cmat(3, 4) += scalar * 0.5 * (invc(3) * invc(4) + invc(5) * invc(1));
  cmat(3, 5) += scalar * 0.5 * (invc(0) * invc(4) + invc(5) * invc(3));

  cmat(4, 0) += scalar * invc(3) * invc(5);
  cmat(4, 1) += scalar * invc(1) * invc(4);
  cmat(4, 2) += scalar * invc(4) * invc(2);
  cmat(4, 3) += scalar * 0.5 * (invc(3) * invc(4) + invc(5) * invc(1));
  cmat(4, 4) += scalar * 0.5 * (invc(1) * invc(2) + invc(4) * invc(4));
  cmat(4, 5) += scalar * 0.5 * (invc(3) * invc(2) + invc(4) * invc(5));

  cmat(5, 0) += scalar * invc(0) * invc(5);
  cmat(5, 1) += scalar * invc(3) * invc(4);
  cmat(5, 2) += scalar * invc(5) * invc(2);
  cmat(5, 3) += scalar * 0.5 * (invc(0) * invc(4) + invc(5) * invc(3));
  cmat(5, 4) += scalar * 0.5 * (invc(3) * invc(2) + invc(4) * invc(5));
  cmat(5, 5) += scalar * 0.5 * (invc(0) * invc(2) + invc(5) * invc(5));
}

template <typename T>
Core::LinAlg::SymmetricTensor<T, 3, 3, 3, 3> Core::LinAlg::FourTensorOperations::holzapfel_product(
    const Core::LinAlg::SymmetricTensor<T, 3, 3>& invc)
{
  auto c_inv_ac_bd = einsum_sym<"ik", "jl">(invc, invc);
  auto c_inv_ad_bc = einsum_sym<"il", "jk">(invc, invc);

  return (0.5) * (c_inv_ac_bd + c_inv_ad_bc);
}

void Core::LinAlg::FourTensorOperations::add_symmetric_holzapfel_product(
    Core::LinAlg::Matrix<6, 6>& X, const Core::LinAlg::Matrix<3, 3>& A,
    const Core::LinAlg::Matrix<3, 3>& B, const double fac)
{
  X(0, 0) += 4 * fac * A(0, 0) * B(0, 0);
  X(0, 3) += fac * (2 * A(0, 0) * B(1, 0) + 2 * A(1, 0) * B(0, 0));
  X(0, 5) += fac * (2 * A(0, 0) * B(2, 0) + 2 * A(2, 0) * B(0, 0));
  X(0, 1) += 4 * fac * A(1, 0) * B(1, 0);
  X(0, 4) += fac * (2 * A(1, 0) * B(2, 0) + 2 * A(2, 0) * B(1, 0));
  X(0, 2) += 4 * fac * A(2, 0) * B(2, 0);

  X(3, 0) += fac * (2 * A(0, 0) * B(0, 1) + 2 * A(0, 1) * B(0, 0));
  X(3, 3) += fac * (A(0, 0) * B(1, 1) + A(1, 0) * B(0, 1) + A(1, 1) * B(0, 0) + A(0, 1) * B(1, 0));
  X(3, 5) += fac * (A(0, 0) * B(2, 1) + A(2, 0) * B(0, 1) + A(2, 1) * B(0, 0) + A(0, 1) * B(2, 0));
  X(3, 1) += fac * (2 * A(1, 0) * B(1, 1) + 2 * A(1, 1) * B(1, 0));
  X(3, 4) += fac * (A(1, 0) * B(2, 1) + A(2, 0) * B(1, 1) + A(2, 1) * B(1, 0) + A(1, 1) * B(2, 0));
  X(3, 2) += fac * (2 * A(2, 0) * B(2, 1) + 2 * A(2, 1) * B(2, 0));

  X(5, 0) += fac * (2 * A(0, 0) * B(0, 2) + 2 * A(0, 2) * B(0, 0));
  X(5, 3) += fac * (A(0, 0) * B(1, 2) + A(1, 0) * B(0, 2) + A(1, 2) * B(0, 0) + A(0, 2) * B(1, 0));
  X(5, 5) += fac * (A(0, 0) * B(2, 2) + A(2, 0) * B(0, 2) + A(2, 2) * B(0, 0) + A(0, 2) * B(2, 0));
  X(5, 1) += fac * (2 * A(1, 0) * B(1, 2) + 2 * A(1, 2) * B(1, 0));
  X(5, 4) += fac * (A(1, 0) * B(2, 2) + A(2, 0) * B(1, 2) + A(2, 2) * B(1, 0) + A(1, 2) * B(2, 0));
  X(5, 2) += fac * (2 * A(2, 0) * B(2, 2) + 2 * A(2, 2) * B(2, 0));

  X(1, 0) += 4 * fac * A(0, 1) * B(0, 1);
  X(1, 3) += fac * (2 * A(0, 1) * B(1, 1) + 2 * A(1, 1) * B(0, 1));
  X(1, 5) += fac * (2 * A(0, 1) * B(2, 1) + 2 * A(2, 1) * B(0, 1));
  X(1, 1) += 4 * fac * A(1, 1) * B(1, 1);
  X(1, 4) += fac * (2 * A(1, 1) * B(2, 1) + 2 * A(2, 1) * B(1, 1));
  X(1, 2) += 4 * fac * A(2, 1) * B(2, 1);

  X(4, 0) += fac * (2 * A(0, 1) * B(0, 2) + 2 * A(0, 2) * B(0, 1));
  X(4, 3) += fac * (A(0, 1) * B(1, 2) + A(1, 1) * B(0, 2) + A(1, 2) * B(0, 1) + A(0, 2) * B(1, 1));
  X(4, 5) += fac * (A(0, 1) * B(2, 2) + A(2, 1) * B(0, 2) + A(2, 2) * B(0, 1) + A(0, 2) * B(2, 1));
  X(4, 1) += fac * (2 * A(1, 1) * B(1, 2) + 2 * A(1, 2) * B(1, 1));
  X(4, 4) += fac * (A(1, 1) * B(2, 2) + A(2, 1) * B(1, 2) + A(2, 2) * B(1, 1) + A(1, 2) * B(2, 1));
  X(4, 2) += fac * (2 * A(2, 1) * B(2, 2) + 2 * A(2, 2) * B(2, 1));

  X(2, 0) += 4 * fac * A(0, 2) * B(0, 2);
  X(2, 3) += fac * (2 * A(0, 2) * B(1, 2) + 2 * A(1, 2) * B(0, 2));
  X(2, 5) += fac * (2 * A(0, 2) * B(2, 2) + 2 * A(2, 2) * B(0, 2));
  X(2, 1) += 4 * fac * A(1, 2) * B(1, 2);
  X(2, 4) += fac * (2 * A(1, 2) * B(2, 2) + 2 * A(2, 2) * B(1, 2));
  X(2, 2) += 4 * fac * A(2, 2) * B(2, 2);
}

template <typename T>
void Core::LinAlg::FourTensorOperations::add_right_non_symmetric_holzapfel_product(
    Core::LinAlg::Matrix<6, 9, T>& out, Core::LinAlg::Matrix<3, 3, T> const& A,
    Core::LinAlg::Matrix<3, 3, T> const& B, T const fac)
{
  out(0, 0) += 2 * fac * A(0, 0) * B(0, 0);
  out(0, 3) += 2 * fac * A(0, 0) * B(0, 1);
  out(0, 5) += 2 * fac * A(0, 0) * B(0, 2);
  out(0, 6) += 2 * fac * A(0, 1) * B(0, 0);
  out(0, 1) += 2 * fac * A(0, 1) * B(0, 1);
  out(0, 4) += 2 * fac * A(0, 1) * B(0, 2);
  out(0, 8) += 2 * fac * A(0, 2) * B(0, 0);
  out(0, 7) += 2 * fac * A(0, 2) * B(0, 1);
  out(0, 2) += 2 * fac * A(0, 2) * B(0, 2);

  out(1, 0) += 2 * fac * A(1, 0) * B(1, 0);
  out(1, 3) += 2 * fac * A(1, 0) * B(1, 1);
  out(1, 5) += 2 * fac * A(1, 0) * B(1, 2);
  out(1, 6) += 2 * fac * A(1, 1) * B(1, 0);
  out(1, 1) += 2 * fac * A(1, 1) * B(1, 1);
  out(1, 4) += 2 * fac * A(1, 1) * B(1, 2);
  out(1, 8) += 2 * fac * A(1, 2) * B(1, 0);
  out(1, 7) += 2 * fac * A(1, 2) * B(1, 1);
  out(1, 2) += 2 * fac * A(1, 2) * B(1, 2);

  out(2, 0) += 2 * fac * A(2, 0) * B(2, 0);
  out(2, 3) += 2 * fac * A(2, 0) * B(2, 1);
  out(2, 5) += 2 * fac * A(2, 0) * B(2, 2);
  out(2, 6) += 2 * fac * A(2, 1) * B(2, 0);
  out(2, 1) += 2 * fac * A(2, 1) * B(2, 1);
  out(2, 4) += 2 * fac * A(2, 1) * B(2, 2);
  out(2, 8) += 2 * fac * A(2, 2) * B(2, 0);
  out(2, 7) += 2 * fac * A(2, 2) * B(2, 1);
  out(2, 2) += 2 * fac * A(2, 2) * B(2, 2);

  out(3, 0) += fac * (A(0, 0) * B(1, 0) + A(1, 0) * B(0, 0));
  out(3, 3) += fac * (A(0, 0) * B(1, 1) + A(1, 0) * B(0, 1));
  out(3, 5) += fac * (A(0, 0) * B(1, 2) + A(1, 0) * B(0, 2));
  out(3, 6) += fac * (A(0, 1) * B(1, 0) + A(1, 1) * B(0, 0));
  out(3, 1) += fac * (A(0, 1) * B(1, 1) + A(1, 1) * B(0, 1));
  out(3, 4) += fac * (A(0, 1) * B(1, 2) + A(1, 1) * B(0, 2));
  out(3, 8) += fac * (A(0, 2) * B(1, 0) + A(1, 2) * B(0, 0));
  out(3, 7) += fac * (A(0, 2) * B(1, 1) + A(1, 2) * B(0, 1));
  out(3, 2) += fac * (A(0, 2) * B(1, 2) + A(1, 2) * B(0, 2));

  out(4, 0) += fac * (A(1, 0) * B(2, 0) + A(2, 0) * B(1, 0));
  out(4, 3) += fac * (A(1, 0) * B(2, 1) + A(2, 0) * B(1, 1));
  out(4, 5) += fac * (A(1, 0) * B(2, 2) + A(2, 0) * B(1, 2));
  out(4, 6) += fac * (A(1, 1) * B(2, 0) + A(2, 1) * B(1, 0));
  out(4, 1) += fac * (A(1, 1) * B(2, 1) + A(2, 1) * B(1, 1));
  out(4, 4) += fac * (A(1, 1) * B(2, 2) + A(2, 1) * B(1, 2));
  out(4, 8) += fac * (A(1, 2) * B(2, 0) + A(2, 2) * B(1, 0));
  out(4, 7) += fac * (A(1, 2) * B(2, 1) + A(2, 2) * B(1, 1));
  out(4, 2) += fac * (A(1, 2) * B(2, 2) + A(2, 2) * B(1, 2));

  out(5, 0) += fac * (A(0, 0) * B(2, 0) + A(2, 0) * B(0, 0));
  out(5, 3) += fac * (A(0, 0) * B(2, 1) + A(2, 0) * B(0, 1));
  out(5, 5) += fac * (A(0, 0) * B(2, 2) + A(2, 0) * B(0, 2));
  out(5, 6) += fac * (A(0, 1) * B(2, 0) + A(2, 1) * B(0, 0));
  out(5, 1) += fac * (A(0, 1) * B(2, 1) + A(2, 1) * B(0, 1));
  out(5, 4) += fac * (A(0, 1) * B(2, 2) + A(2, 1) * B(0, 2));
  out(5, 8) += fac * (A(0, 2) * B(2, 0) + A(2, 2) * B(0, 0));
  out(5, 7) += fac * (A(0, 2) * B(2, 1) + A(2, 2) * B(0, 1));
  out(5, 2) += fac * (A(0, 2) * B(2, 2) + A(2, 2) * B(0, 2));
}

template <typename T>
void Core::LinAlg::FourTensorOperations::add_right_non_symmetric_holzapfel_product_strain_like(
    Core::LinAlg::Matrix<6, 9, T>& out, Core::LinAlg::Matrix<3, 3, T> const& A,
    Core::LinAlg::Matrix<3, 3, T> const& B, T const fac)
{
  out(0, 0) += 2 * fac * A(0, 0) * B(0, 0);
  out(0, 3) += 2 * fac * A(0, 0) * B(0, 1);
  out(0, 5) += 2 * fac * A(0, 0) * B(0, 2);
  out(0, 6) += 2 * fac * A(0, 1) * B(0, 0);
  out(0, 1) += 2 * fac * A(0, 1) * B(0, 1);
  out(0, 4) += 2 * fac * A(0, 1) * B(0, 2);
  out(0, 8) += 2 * fac * A(0, 2) * B(0, 0);
  out(0, 7) += 2 * fac * A(0, 2) * B(0, 1);
  out(0, 2) += 2 * fac * A(0, 2) * B(0, 2);

  out(1, 0) += 2 * fac * A(1, 0) * B(1, 0);
  out(1, 3) += 2 * fac * A(1, 0) * B(1, 1);
  out(1, 5) += 2 * fac * A(1, 0) * B(1, 2);
  out(1, 6) += 2 * fac * A(1, 1) * B(1, 0);
  out(1, 1) += 2 * fac * A(1, 1) * B(1, 1);
  out(1, 4) += 2 * fac * A(1, 1) * B(1, 2);
  out(1, 8) += 2 * fac * A(1, 2) * B(1, 0);
  out(1, 7) += 2 * fac * A(1, 2) * B(1, 1);
  out(1, 2) += 2 * fac * A(1, 2) * B(1, 2);

  out(2, 0) += 2 * fac * A(2, 0) * B(2, 0);
  out(2, 3) += 2 * fac * A(2, 0) * B(2, 1);
  out(2, 5) += 2 * fac * A(2, 0) * B(2, 2);
  out(2, 6) += 2 * fac * A(2, 1) * B(2, 0);
  out(2, 1) += 2 * fac * A(2, 1) * B(2, 1);
  out(2, 4) += 2 * fac * A(2, 1) * B(2, 2);
  out(2, 8) += 2 * fac * A(2, 2) * B(2, 0);
  out(2, 7) += 2 * fac * A(2, 2) * B(2, 1);
  out(2, 2) += 2 * fac * A(2, 2) * B(2, 2);

  out(3, 0) += 2 * fac * (A(0, 0) * B(1, 0) + A(1, 0) * B(0, 0));
  out(3, 3) += 2 * fac * (A(0, 0) * B(1, 1) + A(1, 0) * B(0, 1));
  out(3, 5) += 2 * fac * (A(0, 0) * B(1, 2) + A(1, 0) * B(0, 2));
  out(3, 6) += 2 * fac * (A(0, 1) * B(1, 0) + A(1, 1) * B(0, 0));
  out(3, 1) += 2 * fac * (A(0, 1) * B(1, 1) + A(1, 1) * B(0, 1));
  out(3, 4) += 2 * fac * (A(0, 1) * B(1, 2) + A(1, 1) * B(0, 2));
  out(3, 8) += 2 * fac * (A(0, 2) * B(1, 0) + A(1, 2) * B(0, 0));
  out(3, 7) += 2 * fac * (A(0, 2) * B(1, 1) + A(1, 2) * B(0, 1));
  out(3, 2) += 2 * fac * (A(0, 2) * B(1, 2) + A(1, 2) * B(0, 2));

  out(4, 0) += 2 * fac * (A(1, 0) * B(2, 0) + A(2, 0) * B(1, 0));
  out(4, 3) += 2 * fac * (A(1, 0) * B(2, 1) + A(2, 0) * B(1, 1));
  out(4, 5) += 2 * fac * (A(1, 0) * B(2, 2) + A(2, 0) * B(1, 2));
  out(4, 6) += 2 * fac * (A(1, 1) * B(2, 0) + A(2, 1) * B(1, 0));
  out(4, 1) += 2 * fac * (A(1, 1) * B(2, 1) + A(2, 1) * B(1, 1));
  out(4, 4) += 2 * fac * (A(1, 1) * B(2, 2) + A(2, 1) * B(1, 2));
  out(4, 8) += 2 * fac * (A(1, 2) * B(2, 0) + A(2, 2) * B(1, 0));
  out(4, 7) += 2 * fac * (A(1, 2) * B(2, 1) + A(2, 2) * B(1, 1));
  out(4, 2) += 2 * fac * (A(1, 2) * B(2, 2) + A(2, 2) * B(1, 2));

  out(5, 0) += 2 * fac * (A(0, 0) * B(2, 0) + A(2, 0) * B(0, 0));
  out(5, 3) += 2 * fac * (A(0, 0) * B(2, 1) + A(2, 0) * B(0, 1));
  out(5, 5) += 2 * fac * (A(0, 0) * B(2, 2) + A(2, 0) * B(0, 2));
  out(5, 6) += 2 * fac * (A(0, 1) * B(2, 0) + A(2, 1) * B(0, 0));
  out(5, 1) += 2 * fac * (A(0, 1) * B(2, 1) + A(2, 1) * B(0, 1));
  out(5, 4) += 2 * fac * (A(0, 1) * B(2, 2) + A(2, 1) * B(0, 2));
  out(5, 8) += 2 * fac * (A(0, 2) * B(2, 0) + A(2, 2) * B(0, 0));
  out(5, 7) += 2 * fac * (A(0, 2) * B(2, 1) + A(2, 2) * B(0, 1));
  out(5, 2) += 2 * fac * (A(0, 2) * B(2, 2) + A(2, 2) * B(0, 2));
}

void Core::LinAlg::FourTensorOperations::add_left_non_symmetric_holzapfel_product(
    Core::LinAlg::Matrix<9, 6>& out, Core::LinAlg::Matrix<3, 3> const& A,
    Core::LinAlg::Matrix<3, 3> const& B, double const fac)
{
  out(0, 0) += fac * 2 * A(0, 0) * B(0, 0);
  out(0, 3) += fac * (A(0, 0) * B(0, 1) + A(0, 1) * B(0, 0));
  out(0, 5) += fac * (A(0, 0) * B(0, 2) + A(0, 2) * B(0, 0));
  out(0, 1) += fac * 2 * A(0, 1) * B(0, 1);
  out(0, 4) += fac * (A(0, 1) * B(0, 2) + A(0, 2) * B(0, 1));
  out(0, 2) += fac * 2 * A(0, 2) * B(0, 2);

  out(3, 0) += fac * 2 * A(0, 0) * B(1, 0);
  out(3, 3) += fac * (A(0, 0) * B(1, 1) + A(0, 1) * B(1, 0));
  out(3, 5) += fac * (A(0, 0) * B(1, 2) + A(0, 2) * B(1, 0));
  out(3, 1) += fac * 2 * A(0, 1) * B(1, 1);
  out(3, 4) += fac * (A(0, 1) * B(1, 2) + A(0, 2) * B(1, 1));
  out(3, 2) += fac * 2 * A(0, 2) * B(1, 2);

  out(5, 0) += fac * 2 * A(0, 0) * B(2, 0);
  out(5, 3) += fac * (A(0, 0) * B(2, 1) + A(0, 1) * B(2, 0));
  out(5, 5) += fac * (A(0, 0) * B(2, 2) + A(0, 2) * B(2, 0));
  out(5, 1) += fac * 2 * A(0, 1) * B(2, 1);
  out(5, 4) += fac * (A(0, 1) * B(2, 2) + A(0, 2) * B(2, 1));
  out(5, 2) += fac * 2 * A(0, 2) * B(2, 2);

  out(6, 0) += fac * 2 * A(1, 0) * B(0, 0);
  out(6, 3) += fac * (A(1, 0) * B(0, 1) + A(1, 1) * B(0, 0));
  out(6, 5) += fac * (A(1, 0) * B(0, 2) + A(1, 2) * B(0, 0));
  out(6, 1) += fac * 2 * A(1, 1) * B(0, 1);
  out(6, 4) += fac * (A(1, 1) * B(0, 2) + A(1, 2) * B(0, 1));
  out(6, 2) += fac * 2 * A(1, 2) * B(0, 2);

  out(1, 0) += fac * 2 * A(1, 0) * B(1, 0);
  out(1, 3) += fac * (A(1, 0) * B(1, 1) + A(1, 1) * B(1, 0));
  out(1, 5) += fac * (A(1, 0) * B(1, 2) + A(1, 2) * B(1, 0));
  out(1, 1) += fac * 2 * A(1, 1) * B(1, 1);
  out(1, 4) += fac * (A(1, 1) * B(1, 2) + A(1, 2) * B(1, 1));
  out(1, 2) += fac * 2 * A(1, 2) * B(1, 2);

  out(4, 0) += fac * 2 * A(1, 0) * B(2, 0);
  out(4, 3) += fac * (A(1, 0) * B(2, 1) + A(1, 1) * B(2, 0));
  out(4, 5) += fac * (A(1, 0) * B(2, 2) + A(1, 2) * B(2, 0));
  out(4, 1) += fac * 2 * A(1, 1) * B(2, 1);
  out(4, 4) += fac * (A(1, 1) * B(2, 2) + A(1, 2) * B(2, 1));
  out(4, 2) += fac * 2 * A(1, 2) * B(2, 2);

  out(8, 0) += fac * 2 * A(2, 0) * B(0, 0);
  out(8, 3) += fac * (A(2, 0) * B(0, 1) + A(2, 1) * B(0, 0));
  out(8, 5) += fac * (A(2, 0) * B(0, 2) + A(2, 2) * B(0, 0));
  out(8, 1) += fac * 2 * A(2, 1) * B(0, 1);
  out(8, 4) += fac * (A(2, 1) * B(0, 2) + A(2, 2) * B(0, 1));
  out(8, 2) += fac * 2 * A(2, 2) * B(0, 2);

  out(7, 0) += fac * 2 * A(2, 0) * B(1, 0);
  out(7, 3) += fac * (A(2, 0) * B(1, 1) + A(2, 1) * B(1, 0));
  out(7, 5) += fac * (A(2, 0) * B(1, 2) + A(2, 2) * B(1, 0));
  out(7, 1) += fac * 2 * A(2, 1) * B(1, 1);
  out(7, 4) += fac * (A(2, 1) * B(1, 2) + A(2, 2) * B(1, 1));
  out(7, 2) += fac * 2 * A(2, 2) * B(1, 2);

  out(2, 0) += fac * 2 * A(2, 0) * B(2, 0);
  out(2, 3) += fac * (A(2, 0) * B(2, 1) + A(2, 1) * B(2, 0));
  out(2, 5) += fac * (A(2, 0) * B(2, 2) + A(2, 2) * B(2, 0));
  out(2, 1) += fac * 2 * A(2, 1) * B(2, 1);
  out(2, 4) += fac * (A(2, 1) * B(2, 2) + A(2, 2) * B(2, 1));
  out(2, 2) += fac * 2 * A(2, 2) * B(2, 2);
}

void Core::LinAlg::FourTensorOperations::add_adbc_tensor_product(const double fac,
    const Core::LinAlg::Matrix<3, 3>& A, const Core::LinAlg::Matrix<3, 3>& B,
    Core::LinAlg::Matrix<9, 9>& out)
{
  out(0, 0) += fac * A(0, 0) * B(0, 0);
  out(0, 3) += fac * A(0, 1) * B(0, 0);
  out(0, 5) += fac * A(0, 2) * B(0, 0);
  out(0, 6) += fac * A(0, 0) * B(0, 1);
  out(0, 1) += fac * A(0, 1) * B(0, 1);
  out(0, 4) += fac * A(0, 2) * B(0, 1);
  out(0, 8) += fac * A(0, 0) * B(0, 2);
  out(0, 7) += fac * A(0, 1) * B(0, 2);
  out(0, 2) += fac * A(0, 2) * B(0, 2);

  out(3, 0) += fac * A(0, 0) * B(1, 0);
  out(3, 3) += fac * A(0, 1) * B(1, 0);
  out(3, 5) += fac * A(0, 2) * B(1, 0);
  out(3, 6) += fac * A(0, 0) * B(1, 1);
  out(3, 1) += fac * A(0, 1) * B(1, 1);
  out(3, 4) += fac * A(0, 2) * B(1, 1);
  out(3, 8) += fac * A(0, 0) * B(1, 2);
  out(3, 7) += fac * A(0, 1) * B(1, 2);
  out(3, 2) += fac * A(0, 2) * B(1, 2);

  out(5, 0) += fac * A(0, 0) * B(2, 0);
  out(5, 3) += fac * A(0, 1) * B(2, 0);
  out(5, 5) += fac * A(0, 2) * B(2, 0);
  out(5, 6) += fac * A(0, 0) * B(2, 1);
  out(5, 1) += fac * A(0, 1) * B(2, 1);
  out(5, 4) += fac * A(0, 2) * B(2, 1);
  out(5, 8) += fac * A(0, 0) * B(2, 2);
  out(5, 7) += fac * A(0, 1) * B(2, 2);
  out(5, 2) += fac * A(0, 2) * B(2, 2);

  out(6, 0) += fac * A(1, 0) * B(0, 0);
  out(6, 3) += fac * A(1, 1) * B(0, 0);
  out(6, 5) += fac * A(1, 2) * B(0, 0);
  out(6, 6) += fac * A(1, 0) * B(0, 1);
  out(6, 1) += fac * A(1, 1) * B(0, 1);
  out(6, 4) += fac * A(1, 2) * B(0, 1);
  out(6, 8) += fac * A(1, 0) * B(0, 2);
  out(6, 7) += fac * A(1, 1) * B(0, 2);
  out(6, 2) += fac * A(1, 2) * B(0, 2);

  out(1, 0) += fac * A(1, 0) * B(1, 0);
  out(1, 3) += fac * A(1, 1) * B(1, 0);
  out(1, 5) += fac * A(1, 2) * B(1, 0);
  out(1, 6) += fac * A(1, 0) * B(1, 1);
  out(1, 1) += fac * A(1, 1) * B(1, 1);
  out(1, 4) += fac * A(1, 2) * B(1, 1);
  out(1, 8) += fac * A(1, 0) * B(1, 2);
  out(1, 7) += fac * A(1, 1) * B(1, 2);
  out(1, 2) += fac * A(1, 2) * B(1, 2);

  out(4, 0) += fac * A(1, 0) * B(2, 0);
  out(4, 3) += fac * A(1, 1) * B(2, 0);
  out(4, 5) += fac * A(1, 2) * B(2, 0);
  out(4, 6) += fac * A(1, 0) * B(2, 1);
  out(4, 1) += fac * A(1, 1) * B(2, 1);
  out(4, 4) += fac * A(1, 2) * B(2, 1);
  out(4, 8) += fac * A(1, 0) * B(2, 2);
  out(4, 7) += fac * A(1, 1) * B(2, 2);
  out(4, 2) += fac * A(1, 2) * B(2, 2);

  out(8, 0) += fac * A(2, 0) * B(0, 0);
  out(8, 3) += fac * A(2, 1) * B(0, 0);
  out(8, 5) += fac * A(2, 2) * B(0, 0);
  out(8, 6) += fac * A(2, 0) * B(0, 1);
  out(8, 1) += fac * A(2, 1) * B(0, 1);
  out(8, 4) += fac * A(2, 2) * B(0, 1);
  out(8, 8) += fac * A(2, 0) * B(0, 2);
  out(8, 7) += fac * A(2, 1) * B(0, 2);
  out(8, 2) += fac * A(2, 2) * B(0, 2);

  out(7, 0) += fac * A(2, 0) * B(1, 0);
  out(7, 3) += fac * A(2, 1) * B(1, 0);
  out(7, 5) += fac * A(2, 2) * B(1, 0);
  out(7, 6) += fac * A(2, 0) * B(1, 1);
  out(7, 1) += fac * A(2, 1) * B(1, 1);
  out(7, 4) += fac * A(2, 2) * B(1, 1);
  out(7, 8) += fac * A(2, 0) * B(1, 2);
  out(7, 7) += fac * A(2, 1) * B(1, 2);
  out(7, 2) += fac * A(2, 2) * B(1, 2);

  out(2, 0) += fac * A(2, 0) * B(2, 0);
  out(2, 3) += fac * A(2, 1) * B(2, 0);
  out(2, 5) += fac * A(2, 2) * B(2, 0);
  out(2, 6) += fac * A(2, 0) * B(2, 1);
  out(2, 1) += fac * A(2, 1) * B(2, 1);
  out(2, 4) += fac * A(2, 2) * B(2, 1);
  out(2, 8) += fac * A(2, 0) * B(2, 2);
  out(2, 7) += fac * A(2, 1) * B(2, 2);
  out(2, 2) += fac * A(2, 2) * B(2, 2);
}

void Core::LinAlg::FourTensorOperations::add_non_symmetric_product(double const& fac,
    Core::LinAlg::Matrix<3, 3> const& A, Core::LinAlg::Matrix<3, 3> const& B,
    Core::LinAlg::Matrix<9, 9>& out)
{
  out(0, 0) += fac * A(0, 0) * B(0, 0);
  out(0, 3) += fac * A(0, 0) * B(1, 0);
  out(0, 5) += fac * A(0, 0) * B(2, 0);
  out(0, 6) += fac * A(0, 1) * B(0, 0);
  out(0, 1) += fac * A(0, 1) * B(1, 0);
  out(0, 4) += fac * A(0, 1) * B(2, 0);
  out(0, 8) += fac * A(0, 2) * B(0, 0);
  out(0, 7) += fac * A(0, 2) * B(1, 0);
  out(0, 2) += fac * A(0, 2) * B(2, 0);

  out(3, 0) += fac * A(0, 0) * B(0, 1);
  out(3, 3) += fac * A(0, 0) * B(1, 1);
  out(3, 5) += fac * A(0, 0) * B(2, 1);
  out(3, 6) += fac * A(0, 1) * B(0, 1);
  out(3, 1) += fac * A(0, 1) * B(1, 1);
  out(3, 4) += fac * A(0, 1) * B(2, 1);
  out(3, 8) += fac * A(0, 2) * B(0, 1);
  out(3, 7) += fac * A(0, 2) * B(1, 1);
  out(3, 2) += fac * A(0, 2) * B(2, 1);

  out(5, 0) += fac * A(0, 0) * B(0, 2);
  out(5, 3) += fac * A(0, 0) * B(1, 2);
  out(5, 5) += fac * A(0, 0) * B(2, 2);
  out(5, 6) += fac * A(0, 1) * B(0, 2);
  out(5, 1) += fac * A(0, 1) * B(1, 2);
  out(5, 4) += fac * A(0, 1) * B(2, 2);
  out(5, 8) += fac * A(0, 2) * B(0, 2);
  out(5, 7) += fac * A(0, 2) * B(1, 2);
  out(5, 2) += fac * A(0, 2) * B(2, 2);

  out(6, 0) += fac * A(1, 0) * B(0, 0);
  out(6, 3) += fac * A(1, 0) * B(1, 0);
  out(6, 5) += fac * A(1, 0) * B(2, 0);
  out(6, 6) += fac * A(1, 1) * B(0, 0);
  out(6, 1) += fac * A(1, 1) * B(1, 0);
  out(6, 4) += fac * A(1, 1) * B(2, 0);
  out(6, 8) += fac * A(1, 2) * B(0, 0);
  out(6, 7) += fac * A(1, 2) * B(1, 0);
  out(6, 2) += fac * A(1, 2) * B(2, 0);

  out(1, 0) += fac * A(1, 0) * B(0, 1);
  out(1, 3) += fac * A(1, 0) * B(1, 1);
  out(1, 5) += fac * A(1, 0) * B(2, 1);
  out(1, 6) += fac * A(1, 1) * B(0, 1);
  out(1, 1) += fac * A(1, 1) * B(1, 1);
  out(1, 4) += fac * A(1, 1) * B(2, 1);
  out(1, 8) += fac * A(1, 2) * B(0, 1);
  out(1, 7) += fac * A(1, 2) * B(1, 1);
  out(1, 2) += fac * A(1, 2) * B(2, 1);

  out(4, 0) += fac * A(1, 0) * B(0, 2);
  out(4, 3) += fac * A(1, 0) * B(1, 2);
  out(4, 5) += fac * A(1, 0) * B(2, 2);
  out(4, 6) += fac * A(1, 1) * B(0, 2);
  out(4, 1) += fac * A(1, 1) * B(1, 2);
  out(4, 4) += fac * A(1, 1) * B(2, 2);
  out(4, 8) += fac * A(1, 2) * B(0, 2);
  out(4, 7) += fac * A(1, 2) * B(1, 2);
  out(4, 2) += fac * A(1, 2) * B(2, 2);

  out(8, 0) += fac * A(2, 0) * B(0, 0);
  out(8, 3) += fac * A(2, 0) * B(1, 0);
  out(8, 5) += fac * A(2, 0) * B(2, 0);
  out(8, 6) += fac * A(2, 1) * B(0, 0);
  out(8, 1) += fac * A(2, 1) * B(1, 0);
  out(8, 4) += fac * A(2, 1) * B(2, 0);
  out(8, 8) += fac * A(2, 2) * B(0, 0);
  out(8, 7) += fac * A(2, 2) * B(1, 0);
  out(8, 2) += fac * A(2, 2) * B(2, 0);

  out(7, 0) += fac * A(2, 0) * B(0, 1);
  out(7, 3) += fac * A(2, 0) * B(1, 1);
  out(7, 5) += fac * A(2, 0) * B(2, 1);
  out(7, 6) += fac * A(2, 1) * B(0, 1);
  out(7, 1) += fac * A(2, 1) * B(1, 1);
  out(7, 4) += fac * A(2, 1) * B(2, 1);
  out(7, 8) += fac * A(2, 2) * B(0, 1);
  out(7, 7) += fac * A(2, 2) * B(1, 1);
  out(7, 2) += fac * A(2, 2) * B(2, 1);

  out(2, 0) += fac * A(2, 0) * B(0, 2);
  out(2, 3) += fac * A(2, 0) * B(1, 2);
  out(2, 5) += fac * A(2, 0) * B(2, 2);
  out(2, 6) += fac * A(2, 1) * B(0, 2);
  out(2, 1) += fac * A(2, 1) * B(1, 2);
  out(2, 4) += fac * A(2, 1) * B(2, 2);
  out(2, 8) += fac * A(2, 2) * B(0, 2);
  out(2, 7) += fac * A(2, 2) * B(1, 2);
  out(2, 2) += fac * A(2, 2) * B(2, 2);
}

template <int dim>
void Core::LinAlg::FourTensorOperations::multiply_four_tensor_matrix(
    Core::LinAlg::FourTensor<dim>& four_tensor_result,
    const Core::LinAlg::FourTensor<dim>& four_tensor, const Core::LinAlg::Matrix<dim, dim>& matrix,
    const bool clear_result_tensor)
{
  if (clear_result_tensor) four_tensor_result.clear();
  for (int i = 0; i < dim; ++i)
  {
    for (int j = 0; j < dim; ++j)
    {
      for (int k = 0; k < dim; ++k)
      {
        for (int l = 0; l < dim; ++l)
        {
          for (int m = 0; m < dim; ++m)
          {  // C^ijkl = A^ijkm * B_m^l
            four_tensor_result(i, j, k, l) += four_tensor(i, j, k, m) * matrix(m, l);
          }
        }
      }
    }
  }
}

template <int dim>
void Core::LinAlg::FourTensorOperations::multiply_matrix_four_tensor(
    Core::LinAlg::FourTensor<dim>& four_tensor_result, const Core::LinAlg::Matrix<dim, dim>& matrix,
    const Core::LinAlg::FourTensor<dim>& four_tensor, const bool clear_result_tensor)
{
  if (clear_result_tensor) four_tensor_result.clear();
  for (int i = 0; i < dim; ++i)
  {
    for (int j = 0; j < dim; ++j)
    {
      for (int k = 0; k < dim; ++k)
      {
        for (int l = 0; l < dim; ++l)
        {
          for (int m = 0; m < dim; ++m)
          {
            // C^ijkl = B^i_m * A^mjkl
            four_tensor_result(i, j, k, l) += matrix(i, m) * four_tensor(m, j, k, l);
          }
        }
      }
    }
  }
}

template <int dim>
void Core::LinAlg::FourTensorOperations::multiply_matrix_four_tensor_by_second_index(
    Core::LinAlg::FourTensor<dim>& four_tensor_result, const Core::LinAlg::Matrix<dim, dim>& matrix,
    const Core::LinAlg::FourTensor<dim>& four_tensor, const bool clear_result_tensor)
{
  if (clear_result_tensor) four_tensor_result.clear();
  for (int i = 0; i < dim; ++i)
  {
    for (int j = 0; j < dim; ++j)
    {
      for (int k = 0; k < dim; ++k)
      {
        for (int l = 0; l < dim; ++l)
        {
          for (int m = 0; m < dim; ++m)
          {
            // C^ijkl = B_m^j * A^imkl
            four_tensor_result(i, j, k, l) += matrix(m, j) * four_tensor(i, m, k, l);
          }
        }
      }
    }
  }
}

template <int dim>
void Core::LinAlg::FourTensorOperations::multiply_four_tensor_four_tensor(
    Core::LinAlg::FourTensor<dim>& four_tensor_result,
    const Core::LinAlg::FourTensor<dim>& four_tensor_1,
    const Core::LinAlg::FourTensor<dim>& four_tensor_2, const bool clear_result_tensor)
{
  if (clear_result_tensor) four_tensor_result.clear();
  for (int i = 0; i < dim; ++i)
  {
    for (int j = 0; j < dim; ++j)
    {
      for (int k = 0; k < dim; ++k)
      {
        for (int l = 0; l < dim; ++l)
        {
          // C^ijkl = A^ij_ab * B^abkl
          for (int a = 0; a < dim; ++a)
            for (int b = 0; b < dim; ++b)
              four_tensor_result(i, j, k, l) +=
                  four_tensor_1(i, j, a, b) * four_tensor_2(a, b, k, l);
        }
      }
    }
  }
}

void Core::LinAlg::FourTensorOperations::add_dyadic_product_matrix_matrix(
    Core::LinAlg::FourTensor<3>& four_tensor_result, const Core::LinAlg::Matrix<3, 3>& matrix_A,
    const Core::LinAlg::Matrix<3, 3>& matrix_B)
{
  for (unsigned i = 0; i < 3; ++i)
    for (unsigned j = 0; j < 3; ++j)
      for (unsigned k = 0; k < 3; ++k)
        for (unsigned l = 0; l < 3; ++l)
          four_tensor_result(i, j, k, l) += matrix_A(i, j) * matrix_B(k, l);
}

void Core::LinAlg::FourTensorOperations::add_dyadic_product_matrix_matrix(
    Core::LinAlg::FourTensor<3>& four_tensor_result, const double scale,
    const Core::LinAlg::Matrix<3, 3>& matrix_A, const Core::LinAlg::Matrix<3, 3>& matrix_B)
{
  for (unsigned i = 0; i < 3; ++i)
    for (unsigned j = 0; j < 3; ++j)
      for (unsigned k = 0; k < 3; ++k)
        for (unsigned l = 0; l < 3; ++l)
          four_tensor_result(i, j, k, l) += scale * matrix_A(i, j) * matrix_B(k, l);
}

void Core::LinAlg::FourTensorOperations::add_contraction_matrix_four_tensor(
    Core::LinAlg::Matrix<3, 3>& matrix_result, const Core::LinAlg::Matrix<3, 3>& matrix,
    const Core::LinAlg::FourTensor<3>& four_tensor)
{
  for (unsigned k = 0; k < 3; ++k)
    for (unsigned l = 0; l < 3; ++l)
      for (unsigned i = 0; i < 3; ++i)
        for (unsigned j = 0; j < 3; ++j)
          matrix_result(k, l) += matrix(i, j) * four_tensor(i, j, k, l);
}

void Core::LinAlg::FourTensorOperations::add_contraction_matrix_four_tensor(
    Core::LinAlg::Matrix<3, 3>& matrix_result, const double scale,
    const Core::LinAlg::FourTensor<3>& four_tensor, const Core::LinAlg::Matrix<3, 3>& matrix)
{
  for (unsigned i = 0; i < 3; ++i)
    for (unsigned j = 0; j < 3; ++j)
      for (unsigned k = 0; k < 3; ++k)
        for (unsigned l = 0; l < 3; ++l)
          matrix_result(i, j) += scale * four_tensor(i, j, k, l) * matrix(k, l);
}

double Core::LinAlg::FourTensorOperations::contract_matrix_matrix(
    const Core::LinAlg::Matrix<3, 3>& matrix_A, const Core::LinAlg::Matrix<3, 3>& matrix_B)
{
  double scalarContraction = 0.0;
  for (unsigned i = 0; i < 3; ++i)
    for (unsigned j = 0; j < 3; ++j) scalarContraction += matrix_A(i, j) * matrix_B(i, j);

  return scalarContraction;
}

// explicit instantiation of template functions
template void Core::LinAlg::FourTensorOperations::add_holzapfel_product<double>(
    Core::LinAlg::Matrix<6, 6, double>&, const Core::LinAlg::Matrix<6, 1, double>&,
    const double scalar);
template void Core::LinAlg::FourTensorOperations::add_holzapfel_product<FAD>(
    Core::LinAlg::Matrix<6, 6, FAD>&, const Core::LinAlg::Matrix<6, 1, FAD>&, const FAD scalar);
template Core::LinAlg::SymmetricTensor<double, 3, 3, 3, 3>
Core::LinAlg::FourTensorOperations::holzapfel_product<double>(
    const Core::LinAlg::SymmetricTensor<double, 3, 3>&);
template Core::LinAlg::SymmetricTensor<FAD, 3, 3, 3, 3>
Core::LinAlg::FourTensorOperations::holzapfel_product<FAD>(
    const Core::LinAlg::SymmetricTensor<FAD, 3, 3>&);
template void Core::LinAlg::FourTensorOperations::add_right_non_symmetric_holzapfel_product<double>(
    Core::LinAlg::Matrix<6, 9, double>&, Core::LinAlg::Matrix<3, 3, double> const&,
    Core::LinAlg::Matrix<3, 3, double> const&, double const);
template void Core::LinAlg::FourTensorOperations::add_right_non_symmetric_holzapfel_product<FAD>(
    Core::LinAlg::Matrix<6, 9, FAD>&, Core::LinAlg::Matrix<3, 3, FAD> const&,
    Core::LinAlg::Matrix<3, 3, FAD> const&, FAD const);
template void
Core::LinAlg::FourTensorOperations::add_right_non_symmetric_holzapfel_product_strain_like<double>(
    Core::LinAlg::Matrix<6, 9, double>& out, Core::LinAlg::Matrix<3, 3, double> const& A,
    Core::LinAlg::Matrix<3, 3, double> const& B, double const fac);
template void
Core::LinAlg::FourTensorOperations::add_right_non_symmetric_holzapfel_product_strain_like<FAD>(
    Core::LinAlg::Matrix<6, 9, FAD>& out, Core::LinAlg::Matrix<3, 3, FAD> const& A,
    Core::LinAlg::Matrix<3, 3, FAD> const& B, FAD const fac);
template void Core::LinAlg::FourTensorOperations::multiply_four_tensor_matrix<3>(
    Core::LinAlg::FourTensor<3>& four_tensor_result, const Core::LinAlg::FourTensor<3>& four_tensor,
    const Core::LinAlg::Matrix<3, 3>& matrix, const bool clear_result_tensor);
template void Core::LinAlg::FourTensorOperations::multiply_matrix_four_tensor<3>(
    Core::LinAlg::FourTensor<3>& four_tensor_result, const Core::LinAlg::Matrix<3, 3>& matrix,
    const Core::LinAlg::FourTensor<3>& four_tensor, const bool clear_result_tensor);
template void Core::LinAlg::FourTensorOperations::multiply_matrix_four_tensor_by_second_index<3>(
    Core::LinAlg::FourTensor<3>& four_tensor_result, const Core::LinAlg::Matrix<3, 3>& matrix,
    const Core::LinAlg::FourTensor<3>& four_tensor, const bool clear_result_tensor);
template void Core::LinAlg::FourTensorOperations::multiply_four_tensor_four_tensor<3>(
    Core::LinAlg::FourTensor<3>& four_tensor_result,
    const Core::LinAlg::FourTensor<3>& four_tensor_1,
    const Core::LinAlg::FourTensor<3>& four_tensor_2, const bool clear_result_tensor);



FOUR_C_NAMESPACE_CLOSE
