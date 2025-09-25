// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_LINALG_SERIALDENSEVECTOR_HPP
#define FOUR_C_LINALG_SERIALDENSEVECTOR_HPP


#include "4C_config.hpp"

#include <Teuchos_SerialDenseVector.hpp>

#include <span>

FOUR_C_NAMESPACE_OPEN

namespace Core::LinAlg
{
  /*!
 \brief A class that wraps Teuchos::SerialDenseVector

      This is done in favor of typedef to allow forward declaration
 */
  class SerialDenseVector : public Teuchos::SerialDenseVector<int, double>
  {
   public:
    /// Base type definition
    using Base = Teuchos::SerialDenseVector<int, double>;

    /// Using the base class constructor
    using Base::SerialDenseVector;

    //! Return number of rows. Use our case style to better facilitate generic code.
    [[nodiscard]] int num_rows() const { return this->numRows(); }

    //! Return number of columns. Use our case style to better facilitate generic code.
    //!
    //! @note This function exists because of a design decision in Trilinos where a vector is
    //! implemented as a matrix with one column.
    [[nodiscard]] int num_cols() const { return this->numCols(); }

    /**
     * Get the vector as a std::span.
     */
    [[nodiscard]] std::span<const double> as_span() const
    {
      return std::span(this->values(), this->length());
    }

    /**
     * Get the vector as a std::span.
     */
    [[nodiscard]] std::span<double> as_span() { return std::span(this->values(), this->length()); }
  };

  // type definition for serial integer vector
  using IntSerialDenseVector = Teuchos::SerialDenseVector<int, int>;

  /*!
    \brief Update vector components with scaled values of a,
           b = alpha*a + beta*b
    */
  void update(double alpha, const SerialDenseVector& a, double beta, SerialDenseVector& b);

  // wrapper function to compute Norm of vector
  double norm2(const SerialDenseVector& v);

  // output stream operator
  inline std::ostream& operator<<(std::ostream& out, const SerialDenseVector& vec)
  {
    vec.print(out);
    return out;
  }

  // output stream operator
  inline std::ostream& operator<<(std::ostream& out, const IntSerialDenseVector& vec)
  {
    vec.print(out);
    return out;
  }
}  // namespace Core::LinAlg


FOUR_C_NAMESPACE_CLOSE

#endif
