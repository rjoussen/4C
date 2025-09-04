// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_LINALG_UTILS_DENSEMATRIX_MULTIPLY_HPP
#define FOUR_C_LINALG_UTILS_DENSEMATRIX_MULTIPLY_HPP

#include "4C_config.hpp"

#include "4C_linalg_serialdensematrix.hpp"
#include "4C_linalg_serialdensevector.hpp"
#include "4C_utils_exceptions.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Core::LinAlg::Internal
{
  /*!
   \brief Utility function to get type string of a matrix/vector object
   */
  template <typename VectorOrMatrix>
  inline std::string get_matrix_or_vector_string();

  template <>
  inline std::string get_matrix_or_vector_string<Core::LinAlg::SerialDenseMatrix>()
  {
    return "Matrix";
  }

  template <>
  inline std::string get_matrix_or_vector_string<Core::LinAlg::SerialDenseVector::Base>()
  {
    return "Vector";
  }

  /*!
   \brief Utility function to get the correct case string of a matrix/vector object
   */
  template <typename VectorOrMatrix>
  inline std::string get_matrix_or_vector_case(char ch);

  template <>
  inline std::string get_matrix_or_vector_case<Core::LinAlg::SerialDenseMatrix>(char ch)
  {
    char c = std::toupper(ch);
    std::string s;
    s = c;
    return s;
  }

  template <>
  inline std::string get_matrix_or_vector_case<Core::LinAlg::SerialDenseVector::Base>(char ch)
  {
    char c = std::tolower(ch);
    std::string s;
    s = c;
    return s;
  }

  /*!
   \brief Utility function to get type string of transposition operation
   */
  template <bool transpose>
  inline std::string get_transpose_string();

  template <>
  inline std::string get_transpose_string<false>()
  {
    return "";
  }

  template <>
  inline std::string get_transpose_string<true>()
  {
    return "^T";
  }

  /*!
   \brief Utility function to check for error code of LINALG multiplication
   */
  template <bool transpose1, bool transpose2, typename VectorOrMatrix1, typename VectorOrMatrix2>
  inline void check_error_code_in_debug(
      int errorCode, const VectorOrMatrix1& a, const VectorOrMatrix2& b, const VectorOrMatrix2& c)
  {
    FOUR_C_ASSERT(errorCode == 0,
        "Error code ({}) is returned. Something went wrong with the {}-{} multiplication {} = {}{} "
        "* {}{}. "
        "Dimensions of {}, {}, and {} are ({}x{}), ({}x{}), and ({}x{}) respectively.",
        errorCode, get_matrix_or_vector_string<VectorOrMatrix1>(),
        get_matrix_or_vector_string<VectorOrMatrix2>(),
        get_matrix_or_vector_case<VectorOrMatrix2>('c'),
        get_matrix_or_vector_case<VectorOrMatrix1>('a'), get_transpose_string<transpose1>(),
        get_matrix_or_vector_case<VectorOrMatrix2>('b'), get_transpose_string<transpose2>(),
        get_matrix_or_vector_case<VectorOrMatrix1>('a'),
        get_matrix_or_vector_case<VectorOrMatrix2>('b'),
        get_matrix_or_vector_case<VectorOrMatrix2>('c'), a.numRows(), a.numCols(), b.numRows(),
        b.numCols(), c.numRows(), c.numCols());
  }
}  // namespace Core::LinAlg::Internal

namespace Core::LinAlg
{
  /*!
   \brief Matrix-vector multiplication c = A*b

   \param A (in):        Matrix A
   \param b (in):        Vector b
   \param c (out):       Vector c
   */
  inline int multiply(
      SerialDenseVector::Base& c, const SerialDenseMatrix& A, const SerialDenseVector::Base& b)
  {
    const int err = c.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, 1.0, A.base(), b, 0.0);
    Internal::check_error_code_in_debug<false, false>(err, A, b, c);
    return err;
  }

  /*!
   \brief Matrix-vector multiplication c = alpha*A*b + beta*c

   \param A (in):        Matrix A
   \param b (in):        Vector b
   \param c (out):       Vector c
   */
  inline int multiply(double beta, SerialDenseVector::Base& c, double alpha,
      const SerialDenseMatrix& A, const SerialDenseVector::Base& b)
  {
    const int err = c.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, alpha, A.base(), b, beta);
    Internal::check_error_code_in_debug<false, false>(err, A, b, c);
    return err;
  }

  /*!
   \brief Matrix-vector multiplication c = A^T*b

   \param A (in):        Matrix A
   \param b (in):        Vector b
   \param c (out):       Vector c
   */
  inline int multiply_tn(
      SerialDenseVector::Base& c, const SerialDenseMatrix& A, const SerialDenseVector::Base& b)
  {
    const int err = c.multiply(Teuchos::TRANS, Teuchos::NO_TRANS, 1.0, A.base(), b, 0.0);
    Internal::check_error_code_in_debug<true, false>(err, A, b, c);
    return err;
  }

  /*!
   \brief Matrix-vector multiplication c = alpha*A^T*b + beta*c

   \param A (in):        Matrix A
   \param b (in):        Vector b
   \param c (out):       Vector c
   */
  inline int multiply_tn(double beta, SerialDenseVector::Base& c, double alpha,
      const SerialDenseMatrix& A, const SerialDenseVector::Base& b)
  {
    const int err = c.multiply(Teuchos::TRANS, Teuchos::NO_TRANS, alpha, A.base(), b, beta);
    Internal::check_error_code_in_debug<true, false>(err, A, b, c);
    return err;
  }

  /*!
   \brief Matrix-matrix multiplication C = A*B

   \param A (in):        Matrix A
   \param b (in):        Matrix B
   \param c (out):       Matrix C
   */
  inline int multiply(SerialDenseMatrix& C, const SerialDenseMatrix& A, const SerialDenseMatrix& B)
  {
    const int err =
        C.base().multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, 1.0, A.base(), B.base(), 0.0);
    Internal::check_error_code_in_debug<false, false>(err, A, B, C);
    return err;
  }

  /*!
   \brief Matrix-matrix multiplication C = alpha*A*B + beta*C

   \param A (in):        Matrix A
   \param b (in):        Matrix B
   \param c (out):       Matrix C
   */
  inline int multiply(double beta, SerialDenseMatrix& C, double alpha, const SerialDenseMatrix& A,
      const SerialDenseMatrix& B)
  {
    const int err =
        C.base().multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, alpha, A.base(), B.base(), beta);
    Internal::check_error_code_in_debug<false, false>(err, A, B, C);
    return err;
  }

  /*!
   \brief Matrix-matrix multiplication C = A^T*B

   \param A (in):        Matrix A
   \param b (in):        Matrix B
   \param c (out):       Matrix C
   */
  inline int multiply_tn(
      SerialDenseMatrix& C, const SerialDenseMatrix& A, const SerialDenseMatrix& B)
  {
    const int err =
        C.base().multiply(Teuchos::TRANS, Teuchos::NO_TRANS, 1.0, A.base(), B.base(), 0.0);
    Internal::check_error_code_in_debug<true, false>(err, A, B, C);
    return err;
  }

  /*!
   \brief Matrix-matrix multiplication C = alpha*A^T*B + beta*C

   \param A (in):        Matrix A
   \param b (in):        Matrix B
   \param c (out):       Matrix C
   */
  inline int multiply_tn(double beta, SerialDenseMatrix& C, double alpha,
      const SerialDenseMatrix& A, const SerialDenseMatrix& B)
  {
    const int err =
        C.base().multiply(Teuchos::TRANS, Teuchos::NO_TRANS, alpha, A.base(), B.base(), beta);
    Internal::check_error_code_in_debug<true, false>(err, A, B, C);
    return err;
  }

  /*!
   \brief Matrix-matrix multiplication C = A*B^T

   \param A (in):        Matrix A
   \param b (in):        Matrix B
   \param c (out):       Matrix C
   */
  inline int multiply_nt(
      SerialDenseMatrix& C, const SerialDenseMatrix& A, const SerialDenseMatrix& B)
  {
    const int err =
        C.base().multiply(Teuchos::NO_TRANS, Teuchos::TRANS, 1.0, A.base(), B.base(), 0.0);
    Internal::check_error_code_in_debug<false, true>(err, A, B, C);
    return err;
  }

  /*!
   \brief Matrix-matrix multiplication C = alpha*A*B^T + beta*C

   \param A (in):        Matrix A
   \param b (in):        Matrix B
   \param c (out):       Matrix C
   */
  inline int multiply_nt(double beta, SerialDenseMatrix& C, double alpha,
      const SerialDenseMatrix& A, const SerialDenseMatrix& B)
  {
    const int err =
        C.base().multiply(Teuchos::NO_TRANS, Teuchos::TRANS, alpha, A.base(), B.base(), beta);
    Internal::check_error_code_in_debug<false, true>(err, A, B, C);
    return err;
  }

  /*!
   \brief Matrix-matrix multiplication C = A^T*B^T

   \param A (in):        Matrix A
   \param b (in):        Matrix B
   \param c (out):       Matrix C
   */
  inline int multiply_tt(
      SerialDenseMatrix& C, const SerialDenseMatrix& A, const SerialDenseMatrix& B)
  {
    const int err = C.base().multiply(Teuchos::TRANS, Teuchos::TRANS, 1.0, A.base(), B.base(), 0.0);
    Internal::check_error_code_in_debug<true, true>(err, A, B, C);
    return err;
  }

  /*!
   \brief Matrix-matrix multiplication C = alpha*A^T*B^T + beta*C

   \param A (in):        Matrix A
   \param b (in):        Matrix B
   \param c (out):       Matrix C
   */
  inline int multiply_tt(double beta, SerialDenseMatrix& C, double alpha,
      const SerialDenseMatrix& A, const SerialDenseMatrix& B)
  {
    const int err =
        C.base().multiply(Teuchos::TRANS, Teuchos::TRANS, alpha, A.base(), B.base(), beta);
    Internal::check_error_code_in_debug<true, true>(err, A, B, C);
    return err;
  }
}  // namespace Core::LinAlg

FOUR_C_NAMESPACE_CLOSE

#endif
