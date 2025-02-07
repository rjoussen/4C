// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/*----------------------------------------------------------------------*/
/*! \file

\brief A collection of methods to evaluate the matrix exponential and logarithm, along with specific
derivatives in namespace Core::LinAlg

\level 0
*/
/*----------------------------------------------------------------------*/

#ifndef FOUR_C_LINALG_UTILS_DENSEMATRIX_EXP_LOG_HPP
#define FOUR_C_LINALG_UTILS_DENSEMATRIX_EXP_LOG_HPP

#include "4C_config.hpp"

#include "4C_linalg_fixedsizematrix.hpp"
#include "4C_linalg_utils_densematrix_eigen.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Core::LinAlg
{

  /*!
   * @brief Computes the matrix exponential of a given real matrix, using the Taylor series or a
   * spectral decomposition
   *
   * @note Computation method depends on the norm of the input matrix.
   *
   * @param[in] input input matrix
   * @param[in,out] err_status Error status: 0 if the evaluation was
   * successful, -1 otherwise
   * @returns matrix exponential of the input matrix
   */
  template <unsigned int dim>
  Matrix<dim, dim> matrix_exp(const Matrix<dim, dim>& input, int& err_status);

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
   * @param[in,out] err_status Error status: 0 if the evaluation was
   * successful, -1 otherwise
   * @returns principal matrix logarithm of the input matrix
   */
  template <unsigned int dim>
  Matrix<dim, dim> matrix_log(const Matrix<dim, dim>& input, int& err_status);

  /*!
   * @brief Computes the first derivative of the matrix exponential (general, not necessarily
   * symmetric 3x3 matrix) with respect to its argument
   *
   * For further information, refer to:
   *    -# deSouza, Computational Methods for Plasticity: Theory and Applications, Wiley & Sons,
   * 2008, Section B.2
   *
   * @param[in] input input 3x3 matrix
   * @param[in,out] err_status Error status: 0 if the evaluation was
   * successful, -1 otherwise
   * @return first derivative of input matrix exponential w.r.t. input matrix, specified in Voigt
   * notation
   */
  Matrix<9, 9> matrix_3x3_exp_1st_deriv(const Matrix<3, 3>& input, int& err_status);

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
   * @param[in,out] err_status Error status: 0 if the evaluation was
   * successful, -1 otherwise
   * @return derivative of input matrix logarithm w.r.t. input matrix, specified in Voigt
   * notation
   */
  Matrix<9, 9> matrix_3x3_log_1st_deriv(const Matrix<3, 3>& input, int& err_status);


  /*!
   * @brief Computes the derivative of the matrix exponential (symmetric 3x3 matrix) with respect to
   * its argument
   *
   * For further information, refer to:
   *    -# deSouza, Computational Methods for Plasticity: Theory and Applications, Wiley & Sons,
   * 2008, Section A.5
   *
   * @param[in] input input 3x3 matrix
   * @param[in,out] err_status Error status: 0 if the evaluation was
   * successful, -1 otherwise
   * @return derivative of input matrix exponential w.r.t. input matrix, specified in Voigt
   * stress-stress form
   */
  Matrix<6, 6> sym_matrix_3x3_exp_1st_deriv(const Matrix<3, 3>& input, int& err_status);

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
   * @param[in,out] err_status Error status: 0 if the evaluation was
   * successful, -1 otherwise
   * @param[out] ddexp_mat second derivative of exponential w.r.t. matrix
   * notation
   */
  void sym_matrix_3x3_exp_2nd_deriv_voigt(const Core::LinAlg::Matrix<3, 3>& input,
      Core::LinAlg::Matrix<3, 3>& exp, Core::LinAlg::Matrix<6, 6>& dexp_mat,
      Core::LinAlg::Matrix<6, 6>* ddexp_mat, int& err_status);

}  // namespace Core::LinAlg

FOUR_C_NAMESPACE_CLOSE

#endif
