// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

// This file is part of 4C multiphysics licensed under the // GNU Lesser General Public License v3.0
// or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_LINALG_FIXEDSIZEMATRIX_TENSOR_PRODUCTS_HPP
#define FOUR_C_LINALG_FIXEDSIZEMATRIX_TENSOR_PRODUCTS_HPP

#include "4C_config.hpp"

#include "4C_linalg_fixedsizematrix.hpp"
#include "4C_linalg_four_tensor.hpp"
#include "4C_linalg_symmetric_tensor.hpp"
#include "4C_linalg_tensor.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Core::LinAlg::FourTensorOperations
{
  /*!
   * @brief Multiply two 2nd order tensors A x B and add the result to a 4th order symmetric
   * material tensor in matrix notation
   *
   * In tensor index notation this method does
   * \f[
   * C_{IJKL} := \text{scalar_this} \cdot C_{IJKL}
   *          + \frac{1}{2} \cdot \text{scalar_AB} \cdot \left( A_{IJ} \cdot B_{KL} \right)
   * \f]
   *
   * Compute the "elasticity tensor product" A x B of two 2nd order tensors (in matrix notation) and
   * add the result to a 4th order tensor (in Voigt matrix notation!) using the symmetry-conditions
   * inherent to elasticity tensors.
   *
   * @param[in,out] C       Material tangent matrix to be modified
   * @param[in] scalar_AB    Scalar to multiply with A x B
   * @param[in] A           Dense matrix (3 x 3) as 2nd order tensor A
   * @param[in] B           Dense matrix (3 x 3) as 2nd order tensor B
   * @param[in] scalar_this  Scalar to multiply with C before adding A x B
   */
  void add_elasticity_tensor_product(Core::LinAlg::Matrix<6, 6>& C, const double scalar_AB,
      const Core::LinAlg::Matrix<3, 3>& A, const Core::LinAlg::Matrix<3, 3>& B,
      const double scalar_this);

  /*!
   * @brief Multiply two 2nd order tensors (A x B + B x A) and add the result to a 4th order
   * symmetric material tensor in matrix notation
   *
   * In tensor index notation this method does
   * \f[
   * C_{IJKL} := \text{scalar_this} \cdot C_{IJKL}
   *          + \text{scalar_AB} \cdot \left( A_{IJ} \cdot B_{KL} + B_{IJ} \cdot A_{KL} \right)
   * \f]
   *
   * Compute the "elasticity tensor product" (A x B + B x A) of two 2nd order tensors (in matrix
   * notation) and add the result to a 4th order tensor (in Voigt matrix notation!) using the
   * symmetry-conditions inherent to elasticity tensors
   *
   * @param[in,out] C       Material tangent matrix to be modified
   * @param[in] scalar_AB    Scalar to multiply with (A x B + B x A)
   * @param[in] A           Dense matrix (3 x 3) as 2nd order tensor A
   * @param[in] B           Dense matrix (3 x 3) as 2nd order tensor B
   * @param[in] scalar_this  Scalar to multiply with C before adding (A x B + B x A)
   */
  void add_symmetric_elasticity_tensor_product(Core::LinAlg::Matrix<6, 6>& C,
      const double scalar_AB, const Core::LinAlg::Matrix<3, 3>& A,
      const Core::LinAlg::Matrix<3, 3>& B, const double scalar_this);

  /*!
   * @brief Multiply two 2nd order tensors A o B and add the result to a 4th order symmetric
   * material tensor in matrix notation.
   *
   * In tensor index notation this method does
   * \f[
   * C_{IJKL} := \text{scalar_this} \cdot C_{IJKL} + \frac{1}{2} \cdot \text{scalar_AB} \cdot \left(
   *             A_{IK} \cdot B_{JL} + A_{IL} \cdot B_{JK} \right) \f]
   *
   * Compute the "material tensor product" A o B (also known as kronecker-tensor-product) of two 2nd
   * order tensors (in matrix notation) and add the result to a 4th order tensor (also in matrix
   * notation) using the symmetry-conditions inherent to material tensors, or tangent matrices,
   * respectively AND the Voigt notation of E,S, and C with the famous factor 2! The implementation
   * uses the fixed size matrix.
   *
   * @param[in,out] C       Material tangent matrix to be modified
   * @param[in] scalar_AB    Scalar to multiply with A o B
   * @param[in] A           Dense matrix (3 x 3) as 2nd order tensor A
   * @param[in] B           Dense matrix (3 x 3) as 2nd order tensor B
   * @param[in] scalar_this  Scalar to multiply with C before adding A o B
   */
  void add_kronecker_tensor_product(Core::LinAlg::Matrix<6, 6>& C, const double scalar_AB,
      const Core::LinAlg::Matrix<3, 3>& A, const Core::LinAlg::Matrix<3, 3>& B,
      const double scalar_this);

  /*!
   * @brief Multiply two 2nd order tensors A o B and add the result to a 4th order material tensor
   * in matrix notation, possessing left minor symmetry.
   *
   * In tensor index notation this method does
   * \f[
   * C_{IJKL} := \text{scalar_this} \cdot C_{IJKL} + \frac{1}{2} \cdot \text{scalar_AB} \cdot \left(
   *             A_{IK} \cdot B_{JL} + A_{IL} \cdot B_{JK} \right) \f]
   *
   *
   * @param[in,out] C       Material tangent matrix to be modified
   * @param[in] scalar_AB    Scalar to multiply with A o B
   * @param[in] A           Dense matrix (3 x 3) as 2nd order tensor A
   * @param[in] B           Dense matrix (3 x 3) as 2nd order tensor B
   * @param[in] scalar_this  Scalar to multiply with C before adding A o B
   */
  void add_kronecker_tensor_product(Core::LinAlg::Matrix<6, 9>& C, const double scalar_AB,
      const Core::LinAlg::Matrix<3, 3>& A, const Core::LinAlg::Matrix<3, 3>& B,
      const double scalar_this);


  /*!
   * @brief Add 'Holzapfel product' contribution to constitutive tensor using Voigt notation
   *
   * This function adds the following contribution to the given constitutive
   * matrix \f$\mathbb{C}\f$ (#cmat) based on the inverse of the right Cauchy-Green vector
   * \f$\boldsymbol{C}^{-1}\f$ (#invc):
   *
   * In compact tensor notation this method performs
   * \f[
   *    \mathbb{C} := \mathbb{C} + \text{scalar} \cdot \boldsymbol{C}^{-1} \odot
   *                           \boldsymbol{C}^{-1} \f]
   *
   * which is related to the derivative of the inverse right Cauchy--Green tensor
   * with respect to its non-inverted sibling, i.e.
   * \f[
   *    \frac{\partial \boldsymbol{C}^{-1}}{\partial \boldsymbol{C}}
   *    = - \boldsymbol{C}^{-1} \odot \boldsymbol{C}^{-1}
   * \f]
   * as written by Holzapfel [1], p. 254.
   *
   * In tensor index notation this \f$\odot\f$-product means
   * \f[
   * \mathbb{C}^{ABCD} := \mathbb{C}^{ABCD} + \text{scalar} \left( {C^{-1}}^{AC} {C^{-1}}^{BD} +
   * {C^{-1}}^{AD} {C^{-1}}^{BC} \right) \f]
   *
   * References:
   * [1] G.A. Holzapfel, "Nonlinear solid mechanics", Wiley, 2000.
   *
   * @param[in,out] cmat Material tangent matrix(6x6) to be modified, its rows are denoted
   *                     stress-like 6-Voigt, i.e. \f$[ ()_{11}, \; ()_{22}, \; ()_{33}, \; ()_{12},
   *                     \; ()_{23}, \; ()_{31} ]\f$ and its columns strain-like 6-Voigt, i.e. \f$[
   *                     ()_{11}, \; ()_{22}, \; ()_{33}, \; 2()_{12}, \; 2()_{23}, \; 2()_{31} ]\f$
   *
   * @param[in] invC     Inverse right Cauchy-Green vector in stress-like 6-Voigt notation
   *                     \f$[ C^{-1}_{11} \; C^{-1}_{22} \; C^{-1}_{33} \; C^{-1}_{12} \;
   *                     C^{-1}_{23} \; C^{-1}_{31}]\f$
   *
   * @param[in] scalar   Scalar, corresponding to \f$\delta_7\f$ in Holzapfel [1], p.261
   */
  template <typename T>
  void add_holzapfel_product(Core::LinAlg::Matrix<6, 6, T>& cmat,
      const Core::LinAlg::Matrix<6, 1, T>& invc, const T scalar);

  /*!
   * @brief Calculate the 'Holzapfel product' for a given symmetric 2nd-order tensor
   *
   * This function computes the fourth-order tensor product:
   * \f[
   *    \mathbb{H} = \boldsymbol{C}^{-1} \odot \boldsymbol{C}^{-1}
   * \f]
   * where \f$\odot\f$ denotes the symmetric tensor product defined in index notation as:
   * \f[
   *    \mathbb{H}^{ABCD} = \frac{1}{2} \left(
   *      {C^{-1}}^{AC} {C^{-1}}^{BD} + {C^{-1}}^{AD} {C^{-1}}^{BC}
   *    \right)
   * \f]
   *
   * This is often used to compute the derivative of the inverse right Cauchy-Green tensor
   * with respect to itself:
   * \f[
   *    \mathbb{H} = -\frac{\partial \boldsymbol{C}^{-1}}{\partial \boldsymbol{C}}
   * \f]
   * as derived by Holzapfel [1], p. 254.
   *
   * References:
   * [1] G.A. Holzapfel, "Nonlinear solid mechanics", Wiley, 2000.
   *
   * @tparam T            Scalar type (e.g., double, float)
   * @param[in] invC      Any symmetric 2nd-order Tensor. Mostly the
   *                      inverse right Cauchy-Green tensor in 3D space.
   * @return              Fourth-order symmetric tensor representing
   *                      the Holzapfel product in 3D space.
   */
  template <typename T>
  Core::LinAlg::SymmetricTensor<T, 3, 3, 3, 3> holzapfel_product(
      const Core::LinAlg::SymmetricTensor<T, 3, 3>& invc);

  /*!
   * @brief Add symmetric Holzapfel product to a 4th order tensor in matrix notation
   *
   * In tensor index notation this does:
   * \f[
   *    X_{abcd} += \text{fac} \cdot \left( A_{ca} \cdot B_{db} + A_{da} \cdot B_{cb} + A_{db}
   * \cdot B_{ca} + A_{cb} \cdot B_{da} \right)
   * \f]
   *
   * The result has minor symmetries, but no major symmetry, i.e. symmetric w.r.t. A and B
   *
   * @param[in,out] X 4th order tensor \f[X_{abcd}\f]
   * @param[in] A     2nd order tensor \f[A_{ef}\f]
   * @param[in] B     2nd order tensor \f[B_{gh}\f]
   * @param[in] fac   Scalar factor
   */
  void add_symmetric_holzapfel_product(Core::LinAlg::Matrix<6, 6>& X,
      const Core::LinAlg::Matrix<3, 3>& A, const Core::LinAlg::Matrix<3, 3>& B, const double fac);


  /*!
   * @brief Add right non-symmetric Holzapfel product to a 4th order tensor in matrix notation
   *
   * In tensor index notation this does:
   * \f[
   *    \text{out}_{IJKL} += \text{fac} \cdot \left( A_{IK} \cdot B_{JL}
   *                       + A_{JK} \cdot B_{IL} \right)
   * \f]
   *
   * Note that the result is symmetric within the first two indices. Thus the first two indices
   * are stored in stress-like Voigt notation. The 2nd index pair has no symmetries and is
   * therefore stored as 9-vector. Be aware that this is NOT symmetric w.r.t. A and B !
   *
   * @param[in,out] out     4th order tensor \f[\text{out}_{ABCD}\f]
   * @param[in] A           2nd order tensor \f[A_{EF}\f]
   * @param[in] B           2nd order tensor \f[B_{GH}\f]
   * @param[in] fac         Scalar factor
   */
  template <typename T>
  void add_right_non_symmetric_holzapfel_product(Core::LinAlg::Matrix<6, 9, T>& out,
      Core::LinAlg::Matrix<3, 3, T> const& A, Core::LinAlg::Matrix<3, 3, T> const& B, T const fac);


  /*!
   * @brief Add right non-symmetric Holzapfel product, where the symmetric part is stored
   * strain-like to a 4th order tensor in matrix notation
   *
   * In tensor index notation this does:
   * \f[
   *    \text{out}_{IJKL} += \text{fac} \cdot \left( A_{IK} \cdot B_{JL}
   *                      + A_{JK} \cdot B_{IL} \right)
   * \f]
   *
   * Note that the result is symmetric within the first two indices. The first two indices are
   * stored in strain-like Voigt notation. The 2nd index pair has no symmetries and is therefore
   * stored as 9-vector. Be aware that this is NOT symmetric w.r.t. A and B !
   *
   * @param[in,out] out     4th order tensor \f[\text{out}_{ABCD}\f]
   * @param[in] A           2nd order tensor \f[A_{EF}\f]
   * @param[in] B           2nd order tensor \f[B_{GH}\f]
   * @param[in] fac         Scalar factor
   */
  template <typename T>
  void add_right_non_symmetric_holzapfel_product_strain_like(Core::LinAlg::Matrix<6, 9, T>& out,
      Core::LinAlg::Matrix<3, 3, T> const& A, Core::LinAlg::Matrix<3, 3, T> const& B, T const fac);


  /*!
   * @brief Add left non-symmetric Holzapfel product to a 4th order tensor in matrix notation
   *
   * In tensor index notation this does:
   * \f[
   *    \text{out}_{IJKL} += \text{fac} \cdot \left( A_{IK} \cdot B_{JL}
   *                       + A_{IL} \cdot B_{JK} \right)
   * \f]
   *
   * Note that the result is symmetric within the second two indices. Thus the second two indices
   * are stored in stress-like Voigt notation. The 1st index pair has no symmetries and is
   * therefore stored as 9-vector. Be aware that this is NOT symmetric w.r.t. A and B !
   *
   * @param[in,out] out     4th order tensor \f[\text{out}_{ABCD}\f]
   * @param[in] A           2nd order tensor \f[A_{EF}\f]
   * @param[in] B           2nd order tensor \f[B_{GH}\f]
   * @param[in] fac         Scalar factor
   */
  void add_left_non_symmetric_holzapfel_product(Core::LinAlg::Matrix<9, 6>& out,
      Core::LinAlg::Matrix<3, 3> const& A, Core::LinAlg::Matrix<3, 3> const& B, double const fac);

  /*!
   * @brief Add non-symmetric product to a 4th order tensor in matrix notation
   *
   * In tensor index notation this does:
   * \f[
   *    \text{out}_{ABCD} = \text{fac} \cdot \left( A_{AD} \cdot B_{BC} \right)
   * \f]
   *
   * There are no symmetries and the matrix is therefore stored in a 9x9 matrix. A and B are not
   * required to be symmetric!!!
   *
   * @param[in] fac         Scalar factor
   * @param[in] A           2nd order tensor \f[A_{AB}\f]
   * @param[in] B           2nd order tensor \f[B_{CD}\f]
   * @param[in,out] out     4th order tensor \f[\text{out}_{EFGH}\f]
   */
  void add_adbc_tensor_product(const double fac, const Core::LinAlg::Matrix<3, 3>& A,
      const Core::LinAlg::Matrix<3, 3>& B, Core::LinAlg::Matrix<9, 9>& out);

  /*!
   * @brief Add non-symmetric product to a 4th order tensor in matrix notation
   *
   * In tensor index notation this does:
   * \f[
   *    \text{out}_{ABCD} = \text{fac} \cdot \left( A_{AC} \cdot B_{DB} \right)
   * \f]
   *
   * There are no symmetries and the matrix is therefore stored in a 9x9 matrix. A and B are not
   * required to be symmetric!!!
   *
   * @param[in] fac         Scalar factor
   * @param[in] A           2nd order tensor \f[A_{AB}\f]
   * @param[in] B           2nd order tensor \f[B_{CD}\f]
   * @param[in,out] out     4th order tensor \f[\text{out}_{EFGH}\f]
   */
  void add_non_symmetric_product(double const& fac, Core::LinAlg::Matrix<3, 3> const& A,
      Core::LinAlg::Matrix<3, 3> const& B, Core::LinAlg::Matrix<9, 9>& out);

  /*!
   * @brief Multiply 4th order tensor with 2nd order matrix
   *
   * @param[in,out] four_tensor_result   Result of 4th order tensor - matrix multiplication
   *                                   \f$C^{ijkl} += A^{ijkm} B_m^l\f$
   * @param[in] four_tensor             4th order tensor \f$A^{ijkm}\f$
   * @param[in] matrix                 2-tensor \f$B_m^l\f$
   * @param[in] clear_result_tensor      if true: \f$C^{ijkl}\f$ is cleared when method is called
   */
  template <int dim>
  void multiply_four_tensor_matrix(Core::LinAlg::FourTensor<dim>& four_tensor_result,
      const Core::LinAlg::FourTensor<dim>& four_tensor,
      const Core::LinAlg::Matrix<dim, dim>& matrix, const bool clear_result_tensor = true);

  /*!
   * @brief Multiply 2nd order matrix with 4th order tensor
   *
   * @param[in,out] four_tensor_result   Result of matrix - 4th order tensor multiplication
   *                                   \f$C^{ijkl} += B^i_m A^{mjkl}\f$
   * @param[in] matrix                 2-tensor \f$B^i_m\f$
   * @param[in] four_tensor             4th order tensor \f$A^{mjkl}\f$
   * @param[in] clear_result_tensor      if true: \f$C^{ijkl}\f$ is cleared when method is called
   */
  template <int dim>
  void multiply_matrix_four_tensor(Core::LinAlg::FourTensor<dim>& four_tensor_result,
      const Core::LinAlg::Matrix<dim, dim>& matrix,
      const Core::LinAlg::FourTensor<dim>& four_tensor, const bool clear_result_tensor = true);

  /*!
   * @brief Multiply 2nd order matrix with 4th order tensor by the second index
   *
   * @param[in,out] four_tensor_result   Result of matrix - 4th order tensor multiplication
   *                                   \f$C^{ijkl} += B_m^j A^{imkl}\f$
   * @param[in] matrix                 2nd order tensor \f$B_m^j\f$
   * @param[in] four_tensor             4th order tensor \f$A^{imkl}\f$
   * @param[in] clear_result_tensor      if true: C_ijkl is cleared when method is called
   */
  template <int dim>
  void multiply_matrix_four_tensor_by_second_index(
      Core::LinAlg::FourTensor<dim>& four_tensor_result,
      const Core::LinAlg::Matrix<dim, dim>& matrix,
      const Core::LinAlg::FourTensor<dim>& four_tensor, const bool clear_result_tensor = true);

  /*!
   * @brief Multiply 4th order tensor with 4th order tensor by the first two indices
   *
   * @param[in,out] four_tensor_result   Result of 4th order tensor - 4th order tensor
   * multiplication \f$C^{ijkl} += A^{ij}_{ab} B^{abkl}\f$
   * @param[in] four_tensor_1            4th order tensor \f$A^{ij}_{ab}\f$
   * @param[in] four_tensor_2            4th order tensor \f$B^{abkl}\f$
   * @param[in] clear_result_tensor      if true: C_ijkl is cleared when method is called
   */
  template <int dim>
  void multiply_four_tensor_four_tensor(Core::LinAlg::FourTensor<dim>& four_tensor_result,
      const Core::LinAlg::FourTensor<dim>& four_tensor_1,
      const Core::LinAlg::FourTensor<dim>& four_tensor_2, const bool clear_result_tensor = true);

  /*!
   * @brief Adds the dyadic product of two 2nd order tensors to a 4th order tensor
   *
   * @param[in,out] four_tensor_result    4th order tensor \f$C^{ijkl} = A^{ij} B^{kl}\f$
   * @param[in]     matrix_A             2nd order tensor \f$A^{ij}\f$
   * @param[in]     matrix_B             2nd order tensor \f$B^{kl}\f$
   */
  void add_dyadic_product_matrix_matrix(Core::LinAlg::FourTensor<3>& four_tensor_result,
      const Core::LinAlg::Matrix<3, 3>& matrix_A, const Core::LinAlg::Matrix<3, 3>& matrix_B);

  /*!
   * @brief Adds the dyadic product of two 2nd order tensors to a 4th order tensor
   *
   * @param[in,out] four_tensor_result    4th order tensor \f$C^{ijkl} = A^{ij} B^{kl}\f$
   * @param[in]     matrix_A             2nd order tensor \f$A^{ij}\f$
   * @param[in]     matrix_B             2nd order tensor \f$B^{kl}\f$
   * @param[in]     scale               scaling factor
   */
  void add_dyadic_product_matrix_matrix(Core::LinAlg::FourTensor<3>& four_tensor_result,
      const double scale, const Core::LinAlg::Matrix<3, 3>& matrix_A,
      const Core::LinAlg::Matrix<3, 3>& matrix_B);

  /*!
   * @brief Adds the contraction of a 2nd order tensor and a 4th order tensor
   * to a 2nd order tensor
   *
   * @param[in,out] matrix_result    2nd order tensor \f$C^{kl} = A_{ij} B^{ijkl}\f$
   * @param[in]     matrix          2nd order tensor \f$A_{ij}\f$
   * @param[in]     four_tensor      4th order tensor \f$B^{ijkl}\f$
   */
  void add_contraction_matrix_four_tensor(Core::LinAlg::Matrix<3, 3>& matrix_result,
      const Core::LinAlg::Matrix<3, 3>& matrix, const Core::LinAlg::FourTensor<3>& four_tensor);

  /*!
   * @brief Adds the contraction of a 4th order tensor and a 2nd order tensor
   * to a 2nd order tensor
   *
   * @param[in,out] matrix_result    2nd order tensor \f$C^{ij} = A_{ijkl} B^{kl}\f$
   * @param[in]     scale           the scaling factor
   * @param[in]     four_tensor      4th order tensor \f$B^{kl}\f$
   * @param[in]     matrix          2nd order tensor \f$A_{ijkl}\f$
   */
  void add_contraction_matrix_four_tensor(Core::LinAlg::Matrix<3, 3>& matrix_result,
      const double scale, const Core::LinAlg::FourTensor<3>& four_tensor,
      const Core::LinAlg::Matrix<3, 3>& matrix);

  /*!
   * @brief Returns the double contraction of two 2nd order tensors
   *
   * @param[out]    scalarContraction   0th order tensor \f$s = A_{ij} B^{ij}\f$
   * @param[in]     matrix_A             2nd order tensor \f$A_{ij}\f$
   * @param[in]     matrix_B             2nd order tensor \f$B^{ij}\f$
   */
  double contract_matrix_matrix(
      const Core::LinAlg::Matrix<3, 3>& matrix_A, const Core::LinAlg::Matrix<3, 3>& matrix_B);

}  // namespace Core::LinAlg::FourTensorOperations


FOUR_C_NAMESPACE_CLOSE

#endif