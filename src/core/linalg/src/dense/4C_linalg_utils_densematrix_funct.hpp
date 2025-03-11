// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/*----------------------------------------------------------------------*/
/*! \file

\brief A collection of methods to evaluate matrix functions, such as the square root, the
exponential and the logarithm, along with specific derivatives in namespace Core::LinAlg

\level 0
*/
/*----------------------------------------------------------------------*/

#ifndef FOUR_C_LINALG_UTILS_DENSEMATRIX_FUNCT_HPP
#define FOUR_C_LINALG_UTILS_DENSEMATRIX_FUNCT_HPP

#include "4C_config.hpp"

#include "4C_linalg_fixedsizematrix.hpp"
#include "4C_linalg_utils_densematrix_eigen.hpp"

#include <cstddef>
#include <optional>

FOUR_C_NAMESPACE_OPEN

namespace Core::LinAlg
{
  /// enum class: error types for the calculation of matrix functions
  enum class MatrixFunctErrorType
  {
    NoErrors,          // evaluation without errors
    UnsuitableMethod,  // unsuitable computation method (e.g., the matrix norm is too large to use a
                       // Taylor series description)
    FailedComputation,  // although the computation method is suitable, the evaluation fails
                        // nonetheless (e.g., the maximum number of
                        // iterations is exceeded)
  };

  /// matrix function error type to string
  std::string matrix_funct_err_to_string(MatrixFunctErrorType err_type);


  /// enum class: computation method used for the calculation of the matrix square root
  enum class MatrixSqrtCalcMethod
  {
    DBIterScaledProductForm,  // Product form of Denman and Beavers iteration (scaled), as described
                              // in Higham, Functions of Matrices, Chapter 6: Matrix Square Root,
                              // (6.29)
  };

  /// conversion function for matrix square root calculation method: from string to enum class
  MatrixSqrtCalcMethod matrix_sqrt_calc_string_to_method(std::string matrix_sqrt_calc_string);

  /*!
   * @brief Computes the matrix square root of a given real matrix with a specified computation
   * method.
   *
   *
   * @param[in] input input matrix
   * @param[in,out] err_status Error status in the evaluation
   * @param[out] num_of_iters number of iterations required to compute the matrix square root
   * @param[in] calc_method utilized computation method
   * @returns matrix square root of the input matrix
   */
  template <unsigned int dim>
  Matrix<dim, dim> matrix_sqrt(const Matrix<dim, dim>& input, MatrixFunctErrorType& err_status,
      unsigned int* num_of_iters = nullptr,
      const MatrixSqrtCalcMethod calc_method = MatrixSqrtCalcMethod::DBIterScaledProductForm);

  /// enum class: computation method used for the calculation of the matrix exponential
  enum class MatrixExpCalcMethod
  {
    Default,       // default computation, employing either of the other methods below depending on
                   // different characteristics, such as the matrix norm
    TaylorSeries,  // computation using the Taylor series,
    SpectralDecomp,  // computation using the spectral decomposition,
  };

  /// conversion function for matrix exponential calculation method: from string to enum class
  MatrixExpCalcMethod matrix_exp_calc_string_to_method(std::string matrix_exp_calc_string);

  /*!
   * @brief Computes the matrix exponential of a given real matrix, using the Taylor series or a
   * spectral decomposition
   *
   * @note Computation method depends on the norm of the input matrix.
   *
   * @param[in] input input matrix
   * @param[in,out] err_status Error status in the evaluation
   * @param[in] calc_method utilized computation method
   * @returns matrix exponential of the input matrix
   */
  template <unsigned int dim>
  Matrix<dim, dim> matrix_exp(const Matrix<dim, dim>& input, MatrixFunctErrorType& err_status,
      const MatrixExpCalcMethod calc_method = MatrixExpCalcMethod::Default);

  /// enum class: computation method used for the calculation of the matrix logarithm
  enum class MatrixLogCalcMethod
  {
    DefaultSeries,   // default series computation, employing one of the series descriptions (or the
                     // spectral decomposition) below, depending on different characteristics, such
                     // as the matrix norm
    TaylorSeries,    // computation using the Taylor series,
    GregorySeries,   // computation using the Gregory series,
    SpectralDecomp,  // computation using the spectral decomposition,
    InvScalSquare,   // computation using the inverse scaling and squaring algorithm presented in
                     // Higham: Functions of Matrices, Chapter 11: Matrix Logarithm, Algorithm 11.10
    PadePartFract,   // Pade approximation using a partial fraction expansion, as presented in
                     // Higham: Functions of Matrices, Chapter 11: Matrix Logarithm, Eq. 11.18
  };

  /// conversion function for matrix logarithm calculation method: from string to enum class
  MatrixLogCalcMethod matrix_log_calc_string_to_method(std::string matrix_log_calc_string);

  /*!
   * @brief Computes the (principal) matrix logarithm of a given real matrix, using either the
   * Taylor series or the Gregory series or a spectral decomposition
   *
   *
   * For further information, refer to:
   *    -# Higham, Functions of Matrices: Theory and Computation, Society for Industrial and
   * Applied Mathematics, 2008
   *
   * @note The current implementation only takes input matrices into account, where all eigenvalues
   * possess positive real parts.
   *
   * @param[in] input input matrix
   * @param[in,out] err_status Error status in the evaluation
   * @param[in] calc_method utilized computation method
   * @param[in, out] pade_order Pade approximation order (only for
   * algorithms employing this). This is an input or output parameter
   * depending on the employed method (e.g.: input for Pade
   * approximation, output for inverse scaling and squaring method)
   * @returns principal matrix logarithm of the input matrix
   */
  template <unsigned int dim>
  Matrix<dim, dim> matrix_log(const Matrix<dim, dim>& input, MatrixFunctErrorType& err_status,
      const MatrixLogCalcMethod calc_method = MatrixLogCalcMethod::DefaultSeries,
      unsigned int* pade_order = nullptr);

  /// enum class: computation method used for the calculation of the first derivative of the matrix
  /// exponential for a general, not necessarily symmetric matrix
  enum class GenMatrixExpFirstDerivCalcMethod
  {
    Default,       // default computation, employing either of the other methods below depending on
                   // different characteristics, such as the matrix norm
    TaylorSeries,  // computation using the Taylor series,
  };

  /// conversion function for the calculation method of the first derivative of a general
  /// matrix' exponential: from string to enum class
  GenMatrixExpFirstDerivCalcMethod genmatrix_exp_1st_deriv_calc_string_to_method(
      std::string genmatrix_exp_1st_deriv_calc_string);

  /*!
   * @brief Computes the first derivative of the matrix exponential (general, not necessarily
   * symmetric 3x3 matrix) with respect to its argument
   *
   * For further information, refer to:
   *    -# deSouza, Computational Methods for Plasticity: Theory and Applications, Wiley & Sons,
   * 2008, Section B.2
   *
   * @param[in] input input 3x3 matrix
   * @param[in,out] err_status Error status in the evaluation
   * @param[in] calc_method utilized computation method
   * @return first derivative of input matrix exponential w.r.t. input matrix, specified in Voigt
   * notation
   */
  Matrix<9, 9> matrix_3x3_exp_1st_deriv(const Matrix<3, 3>& input, MatrixFunctErrorType& err_status,
      GenMatrixExpFirstDerivCalcMethod calc_method = GenMatrixExpFirstDerivCalcMethod::Default);

  /// enum class: computation method used for the calculation of the first derivative of the matrix
  /// logarithm for a general, not necessarily symmetric matrix
  enum class GenMatrixLogFirstDerivCalcMethod
  {
    DefaultSeries,  // default series computation, employing either of the series below depending on
                    // different characteristics, such as the matrix norm
    TaylorSeries,   // computation using the Taylor series,
    GregorySeries,  // computation using the Gregory series,
    PadePartFract,  // Pade approximation using a partial fraction expansion, as presented in
                    // Higham: Functions of Matrices, Chapter 11: Matrix Logarithm, Eq. 11.18
  };

  /// conversion function for the calculation method of the first derivative of a general
  /// matrix' logarithm: from string to enum class
  GenMatrixLogFirstDerivCalcMethod genmatrix_log_1st_deriv_calc_string_to_method(
      std::string genmatrix_log_1st_deriv_calc_string);

  /*!
   * @brief Computes the derivative of the matrix logarithm (general, not necessarily symmetric
   * matrix) with respect to its argument, using either the Taylor series or
   * the Gregory series
   *
   * For further information, refer to:
   *    -# Higham, Functions of Matrices: Theory and Computation, Society for Industrial and Applied
   * Mathematics, 2008
   *
   * @param[in] input input 3x3 matrix
   * @param[in,out] err_status Error status in the evaluation
   * @param[in] calc_method utilized computation method
   * @param[in] pade_order Pade approximation order (only for
   * algorithms employing this)
   * @return derivative of input matrix logarithm w.r.t. input
   * matrix, specified in Voigt notation
   */
  Matrix<9, 9> matrix_3x3_log_1st_deriv(const Matrix<3, 3>& input, MatrixFunctErrorType& err_status,
      GenMatrixLogFirstDerivCalcMethod calc_method =
          GenMatrixLogFirstDerivCalcMethod::DefaultSeries,
      const unsigned int* pade_order = nullptr);

  /// enum class: computation method used for the calculation of the first derivative of the matrix
  /// exponential for a symmetric matrix
  enum class SymMatrixExpFirstDerivCalcMethod
  {
    Default,       // default computation, employing either of the other methods below depending on
                   // different characteristics, such as the matrix norm
    TaylorSeries,  // computation using the Taylor series,
    EigenprojBased,  // computation based on eigenprojections, as shown in deSouza, Computational
                     // Methods for Plasticity: Theory and Applications, Wiley & Sons,
                     // 2008, Section A.5
  };

  /// conversion function for the calculation method of the first derivative of a symmetric
  /// matrix' exponential: from string to enum class
  SymMatrixExpFirstDerivCalcMethod symmatrix_exp_1st_deriv_calc_string_to_method(
      std::string symmatrix_exp_1st_deriv_calc_string);

  /*!
   * @brief Computes the derivative of the matrix exponential (symmetric 3x3 matrix) with respect to
   * its argument
   *
   * For further information, refer to:
   *    -# deSouza, Computational Methods for Plasticity: Theory and Applications, Wiley & Sons,
   * 2008, Section A.5 / Box B.2
   *
   * @param[in] input input 3x3 matrix
   * @param[in,out] err_status Error status in the evaluation
   * @param[in] calc_method utilized computation method
   * @return derivative of input matrix exponential w.r.t. input matrix, specified in Voigt
   * stress-stress form
   */
  Matrix<6, 6> sym_matrix_3x3_exp_1st_deriv(const Matrix<3, 3>& input,
      MatrixFunctErrorType& err_status,
      SymMatrixExpFirstDerivCalcMethod calc_method = SymMatrixExpFirstDerivCalcMethod::Default);

  /*!
   * @brief Computes the exponential of a symmetric matrix along with the first and second
   * derivatives with respect to the argument, in Voigt notation
   *
   *  For further information, refer to:
   *    -# Ortiz et al., The computation of the exponential and logarithmic mappings and their first
   * and second linearizations, Int. J. Number. Meth. Engng 52, 2001
   *
   * @param[in] input input 3x3 matrix
   * @param[out] exp exponential of the input matrix
   * @param[out] dexp_mat first derivative of exponential w.r.t. matrix
   * @param[in,out] err_status Error status in the evaluation
   * @param[out] ddexp_mat second derivative of exponential w.r.t. matrix
   * notation
   */
  void sym_matrix_3x3_exp_2nd_deriv_voigt(const Core::LinAlg::Matrix<3, 3>& input,
      Core::LinAlg::Matrix<3, 3>& exp, Core::LinAlg::Matrix<6, 6>& dexp_mat,
      Core::LinAlg::Matrix<6, 6>* ddexp_mat, MatrixFunctErrorType& err_status);

}  // namespace Core::LinAlg

FOUR_C_NAMESPACE_CLOSE

#endif
