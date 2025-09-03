// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_LINALG_SERIALDENSEMATRIX_HPP
#define FOUR_C_LINALG_SERIALDENSEMATRIX_HPP

#include "4C_config.hpp"

#include <Teuchos_DataAccess.hpp>
#include <Teuchos_SerialDenseMatrix.hpp>

FOUR_C_NAMESPACE_OPEN

namespace Core::LinAlg
{
  /*!
   * \brief A wrapper around Teuchos::SerialDenseMatrix
   */
  class SerialDenseMatrix
  {
   public:
    using ordinalType = int;
    using scalarType = double;
    using Base = Teuchos::SerialDenseMatrix<ordinalType, scalarType>;

    // --- constructors ---
    SerialDenseMatrix() = default;
    SerialDenseMatrix(int rows, int cols);
    SerialDenseMatrix(int rows, int cols, bool zeroOut);
    SerialDenseMatrix(const Base& src, Teuchos::ETransp trans = Teuchos::NO_TRANS);
    SerialDenseMatrix(Teuchos::DataAccess cv, double* values, int stride, int rows, int cols);

    // NOLINTBEGIN(readability-identifier-naming)

    // --- size queries ---
    int num_rows() const;
    int num_cols() const;
    int numRows() const;
    int numCols() const;
    int stride() const;
    double normInf() const;
    double normOne() const;

    // --- element access ---
    double& operator()(int i, int j) { return mat_(i, j); };
    const double& operator()(int i, int j) const { return mat_(i, j); };

    // --- data access ---
    //! Returns a pointer to the raw data in column-major order.
    double* values() const;

    // --- modifiers ---
    void scale(double alpha);
    void put_scalar(double val);
    void putScalar(double val);
    void shape(int rows, int cols);
    void reshape(int rows, int cols);
    void assign(const SerialDenseMatrix& source);

    // NOLINTEND(readability-identifier-naming)

    // --- algebraic operations ---
    SerialDenseMatrix& operator+=(const SerialDenseMatrix& other);

    // --- access to underlying Trilinos object ---
    Base& base();
    const Base& base() const;

   private:
    Base mat_;
  };

  // ===============================================================
  //  Free functions operating on SerialDenseMatrix
  // ===============================================================

  /*!
   * \brief Update matrix components with scaled values of A,
   *        B = alpha*A + beta*B
   */
  void update(double alpha, const SerialDenseMatrix& A, double beta, SerialDenseMatrix& B);

  /*!
   * \brief Zero out first n elements of the matrix (column-major order).
   */
  void zero(SerialDenseMatrix& mat, int length);

  /*!
   * \brief Determinant computation using recursive expansion (Sarrus rule).
   * \note Only safe for n <= 4 due to complexity.
   */
  long double det_long(const SerialDenseMatrix& matrix);

  /*!
   * \brief Utility function to copy raw vector data into a SerialDenseMatrix (column-major).
   */
  void copy(const double* vec, SerialDenseMatrix& mat);

  //! Output stream operator
  inline std::ostream& operator<<(std::ostream& out, const SerialDenseMatrix& mat)
  {
    mat.base().print(out);
    return out;
  }

}  // namespace Core::LinAlg

FOUR_C_NAMESPACE_CLOSE

#endif