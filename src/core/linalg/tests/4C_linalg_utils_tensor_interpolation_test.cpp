// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <gtest/gtest.h>

#include "4C_linalg_utils_tensor_interpolation.hpp"

#include "4C_linalg_fixedsizematrix.hpp"
#include "4C_unittest_utils_assertions_test.hpp"

#include <Teuchos_ParameterList.hpp>

#include <array>
#include <cmath>
#include <iomanip>
#include <ios>
#include <utility>
#include <vector>

FOUR_C_NAMESPACE_OPEN


namespace
{
  TEST(LinalgTensorInterpolationTest, 2Sym3x3MatrixInterp2D)
  {
    // see Satheesh et al., 2024, 10.1002/nme.7373, Section 3.1.1, R-LOG interpolation

    Core::LinAlg::Matrix<3, 3> temp3x3(true);

    // create interpolation parameter list and append exponential decay
    // factor
    Teuchos::ParameterList param_list;
    param_list.set("weighting_expdecay_factor", 10.0);

    // create tensor interpolator
    Core::LinAlg::TensorInterpolation::SecondOrderTensorInterpolator<1> interp{1,
        Core::LinAlg::TensorInterpolation::RotInterpType::RInterp,
        Core::LinAlg::TensorInterpolation::EigenvalInterpType::LOG, param_list};

    std::vector<Core::LinAlg::Matrix<3, 3>> ref_matrices;
    std::vector<Core::LinAlg::Matrix<1, 1>> ref_locs;

    Core::LinAlg::Matrix<3, 3> lambda_T_1(true);
    lambda_T_1(0, 0) = 10.0;
    lambda_T_1(1, 1) = 1.0;
    lambda_T_1(2, 2) = 1.0;
    Core::LinAlg::Matrix<3, 3> Q_T_1(true);
    double angle_T_1 = M_PI / 4.0;
    Q_T_1(0, 0) = std::cos(angle_T_1);
    Q_T_1(0, 1) = -std::sin(angle_T_1);
    Q_T_1(1, 0) = std::sin(angle_T_1);
    Q_T_1(1, 1) = std::cos(angle_T_1);
    Q_T_1(2, 2) = 1.0;
    Core::LinAlg::Matrix<3, 3> T_1(true);
    temp3x3.multiply_tn(1.0, Q_T_1, lambda_T_1, 0.0);
    T_1.multiply_nn(1.0, temp3x3, Q_T_1, 0.0);
    Core::LinAlg::Matrix<1, 1> loc_T_1(true);
    loc_T_1(0, 0) = -5.0;

    Core::LinAlg::Matrix<3, 3> lambda_T_2(true);
    lambda_T_2(0, 0) = 20.0;
    lambda_T_2(1, 1) = 4.0;
    lambda_T_2(2, 2) = 1.0;
    Core::LinAlg::Matrix<3, 3> Q_T_2(true);
    double angle_T_2 = -0.99 * M_PI / 4.0;
    Q_T_2(0, 0) = std::cos(angle_T_2);
    Q_T_2(0, 1) = -std::sin(angle_T_2);
    Q_T_2(1, 0) = std::sin(angle_T_2);
    Q_T_2(1, 1) = std::cos(angle_T_2);
    Q_T_2(2, 2) = 1.0;
    Core::LinAlg::Matrix<3, 3> T_2(true);
    temp3x3.multiply_tn(1.0, Q_T_2, lambda_T_2, 0.0);
    T_2.multiply_nn(1.0, temp3x3, Q_T_2, 0.0);
    Core::LinAlg::Matrix<1, 1> loc_T_2(true);
    loc_T_2(0, 0) = 5.0;


    ref_matrices.push_back(T_1);
    ref_matrices.push_back(T_2);
    ref_locs.push_back(loc_T_1);
    ref_locs.push_back(loc_T_2);

    Core::LinAlg::Matrix<1, 1> loc_interp(true);

    Core::LinAlg::Matrix<3, 3> T_interp =
        interp.get_interpolated_matrix(ref_matrices, ref_locs, loc_interp);

    // reference result
    Core::LinAlg::Matrix<3, 3> lambda_T_ref(true);
    lambda_T_ref(0, 0) =
        std::pow(lambda_T_1(0, 0), 1.0 / 2.0) * std::pow(lambda_T_2(0, 0), 1.0 / 2.0);
    lambda_T_ref(1, 1) =
        std::pow(lambda_T_1(1, 1), 1.0 / 2.0) * std::pow(lambda_T_2(1, 1), 1.0 / 2.0);
    lambda_T_ref(2, 2) =
        std::pow(lambda_T_1(2, 2), 1.0 / 2.0) * std::pow(lambda_T_2(2, 2), 1.0 / 2.0);
    Core::LinAlg::Matrix<3, 3> Q_T_ref(true);
    double angle_T_ref = (1 - 0.99) / 2.0 * M_PI / 4.0;
    Q_T_ref(0, 0) = std::cos(angle_T_ref);
    Q_T_ref(0, 1) = -std::sin(angle_T_ref);
    Q_T_ref(1, 0) = std::sin(angle_T_ref);
    Q_T_ref(1, 1) = std::cos(angle_T_ref);
    Q_T_ref(2, 2) = 1.0;
    Core::LinAlg::Matrix<3, 3> T_ref(true);
    temp3x3.multiply_tn(1.0, Q_T_ref, lambda_T_ref, 0.0);
    T_ref.multiply_nn(1.0, temp3x3, Q_T_ref, 0.0);


    FOUR_C_EXPECT_NEAR(T_interp, T_ref, 1.0e-9);
  }


  TEST(LinalgTensorInterpolationTest, 2NonSym3x3MatrixInterp2D)
  {
    // see Satheesh et al., 2024, 10.1002/nme.7373, Section 3.2.1, R-LOG interpolation

    Core::LinAlg::Matrix<3, 3> temp3x3(true);

    // create interpolation parameter list and append exponential decay
    // factor
    Teuchos::ParameterList param_list;
    param_list.set("weighting_expdecay_factor", 10.0);

    // create tensor interpolator
    Core::LinAlg::TensorInterpolation::SecondOrderTensorInterpolator<1> interp{1,
        Core::LinAlg::TensorInterpolation::RotInterpType::RInterp,
        Core::LinAlg::TensorInterpolation::EigenvalInterpType::LOG, param_list};

    std::vector<Core::LinAlg::Matrix<3, 3>> ref_matrices;
    std::vector<Core::LinAlg::Matrix<1, 1>> ref_locs;

    Core::LinAlg::Matrix<3, 3> lambda_T_1(true);
    lambda_T_1(0, 0) = 10.0;
    lambda_T_1(1, 1) = 1.0;
    lambda_T_1(2, 2) = 1.0;
    Core::LinAlg::Matrix<3, 3> Q_T_1(true);
    double angle_T_1_U = M_PI / 4.0;
    Q_T_1(0, 0) = std::cos(angle_T_1_U);
    Q_T_1(0, 1) = -std::sin(angle_T_1_U);
    Q_T_1(1, 0) = std::sin(angle_T_1_U);
    Q_T_1(1, 1) = std::cos(angle_T_1_U);
    Q_T_1(2, 2) = 1.0;
    Core::LinAlg::Matrix<3, 3> R_T_1(true);
    double angle_T_1_R = M_PI / 2.0;
    R_T_1(0, 0) = std::cos(angle_T_1_R);
    R_T_1(0, 1) = -std::sin(angle_T_1_R);
    R_T_1(1, 0) = std::sin(angle_T_1_R);
    R_T_1(1, 1) = std::cos(angle_T_1_R);
    R_T_1(2, 2) = 1.0;
    Core::LinAlg::Matrix<3, 3> T_1(true);
    T_1.multiply_tn(1.0, Q_T_1, lambda_T_1, 0.0);
    temp3x3.multiply_nn(1.0, T_1, Q_T_1, 0.0);
    T_1.multiply_nn(1.0, R_T_1, temp3x3, 0.0);
    Core::LinAlg::Matrix<1, 1> loc_T_1(true);
    loc_T_1(0, 0) = -5.0;

    Core::LinAlg::Matrix<3, 3> lambda_T_2(true);
    lambda_T_2(0, 0) = 20.0;
    lambda_T_2(1, 1) = 4.0;
    lambda_T_2(2, 2) = 1.0;
    Core::LinAlg::Matrix<3, 3> Q_T_2(true);
    double angle_T_2 = -0.99 * M_PI / 4.0;
    Q_T_2(0, 0) = std::cos(angle_T_2);
    Q_T_2(0, 1) = -std::sin(angle_T_2);
    Q_T_2(1, 0) = std::sin(angle_T_2);
    Q_T_2(1, 1) = std::cos(angle_T_2);
    Q_T_2(2, 2) = 1.0;
    Core::LinAlg::Matrix<3, 3> T_2(true);
    temp3x3.multiply_tn(1.0, Q_T_2, lambda_T_2, 0.0);
    T_2.multiply_nn(1.0, temp3x3, Q_T_2, 0.0);
    Core::LinAlg::Matrix<1, 1> loc_T_2(true);
    loc_T_2(0, 0) = 5.0;


    ref_matrices.push_back(T_1);
    ref_matrices.push_back(T_2);
    ref_locs.push_back(loc_T_1);
    ref_locs.push_back(loc_T_2);

    Core::LinAlg::Matrix<1, 1> loc_interp(true);

    Core::LinAlg::Matrix<3, 3> T_interp =
        interp.get_interpolated_matrix(ref_matrices, ref_locs, loc_interp);

    // reference result
    Core::LinAlg::Matrix<3, 3> lambda_T_ref(true);
    lambda_T_ref(0, 0) =
        std::pow(lambda_T_1(0, 0), 1.0 / 2.0) * std::pow(lambda_T_2(0, 0), 1.0 / 2.0);
    lambda_T_ref(1, 1) =
        std::pow(lambda_T_1(1, 1), 1.0 / 2.0) * std::pow(lambda_T_2(1, 1), 1.0 / 2.0);
    lambda_T_ref(2, 2) =
        std::pow(lambda_T_1(2, 2), 1.0 / 2.0) * std::pow(lambda_T_2(2, 2), 1.0 / 2.0);
    Core::LinAlg::Matrix<3, 3> Q_T_ref(true);
    double angle_T_ref_U = (1 - 0.99) / 2.0 * M_PI / 4.0;
    Q_T_ref(0, 0) = std::cos(angle_T_ref_U);
    Q_T_ref(0, 1) = -std::sin(angle_T_ref_U);
    Q_T_ref(1, 0) = std::sin(angle_T_ref_U);
    Q_T_ref(1, 1) = std::cos(angle_T_ref_U);
    Q_T_ref(2, 2) = 1.0;
    Core::LinAlg::Matrix<3, 3> R_T_ref(true);
    double angle_T_ref_R = M_PI / 4.0;
    R_T_ref(0, 0) = std::cos(angle_T_ref_R);
    R_T_ref(0, 1) = -std::sin(angle_T_ref_R);
    R_T_ref(1, 0) = std::sin(angle_T_ref_R);
    R_T_ref(1, 1) = std::cos(angle_T_ref_R);
    R_T_ref(2, 2) = 1.0;
    Core::LinAlg::Matrix<3, 3> T_ref(true);
    T_ref.multiply_tn(1.0, Q_T_ref, lambda_T_ref, 0.0);
    temp3x3.multiply_nn(1.0, T_ref, Q_T_ref, 0.0);
    T_ref.multiply_nn(1.0, R_T_ref, temp3x3, 0.0);

    FOUR_C_EXPECT_NEAR(T_interp, T_ref, 1.0e-9);
  }

  TEST(LinalgTensorInterpolationTest, General1DInterpTest)
  {
    // set left matrix (loc = 0)
    Core::LinAlg::Matrix<3, 3> left_matrix{true};
    left_matrix(0, 0) = 1.802440e-01;
    left_matrix(0, 1) = 7.461039e-01;
    left_matrix(0, 2) = -3.148670e-10;
    left_matrix(1, 0) = -5.701622e-01;
    left_matrix(1, 1) = 3.263754e+00;
    left_matrix(1, 2) = -8.927429e-10;
    left_matrix(2, 0) = 1.019112e-10;
    left_matrix(2, 1) = 1.152076e-10;
    left_matrix(2, 2) = 9.865120e-01;
    // perform polar decomposition of left matrix
    Core::LinAlg::Matrix<3, 3> left_matrix_R{true};
    Core::LinAlg::Matrix<3, 3> left_matrix_U{true};
    Core::LinAlg::Matrix<3, 3> left_matrix_lambda{true};
    std::array<std::pair<double, Core::LinAlg::Matrix<3, 1>>, 3> left_matrix_spectral;
    Core::LinAlg::matrix_3x3_polar_decomposition(
        left_matrix, left_matrix_R, left_matrix_U, left_matrix_lambda, left_matrix_spectral);
    // build Q left matrix
    Core::LinAlg::Matrix<3, 3> left_matrix_Q{true};
    left_matrix_Q(0, 0) = left_matrix_spectral[0].second(0);
    left_matrix_Q(0, 1) = left_matrix_spectral[0].second(1);
    left_matrix_Q(0, 2) = left_matrix_spectral[0].second(2);
    left_matrix_Q(1, 0) = left_matrix_spectral[1].second(0);
    left_matrix_Q(1, 1) = left_matrix_spectral[1].second(1);
    left_matrix_Q(1, 2) = left_matrix_spectral[1].second(2);
    left_matrix_Q(2, 0) = left_matrix_spectral[2].second(0);
    left_matrix_Q(2, 1) = left_matrix_spectral[2].second(1);
    left_matrix_Q(2, 2) = left_matrix_spectral[2].second(2);
    // get the rotation vectors
    Core::LinAlg::Matrix<3, 1> left_matrix_R_vect =
        Core::LinAlg::calc_rot_vect_from_rot_matrix(left_matrix_R);
    Core::LinAlg::Matrix<3, 1> left_matrix_Q_vect =
        Core::LinAlg::calc_rot_vect_from_rot_matrix(left_matrix_Q);

    // set right matrix (loc = 1)
    Core::LinAlg::Matrix<3, 3> right_matrix{true};
    right_matrix(0, 0) = 4.337834e-01;
    right_matrix(0, 1) = 6.330147e-01;
    right_matrix(0, 2) = -9.099803e-11;
    right_matrix(1, 0) = 6.330147e-01;
    right_matrix(1, 1) = 3.260640e+00;
    right_matrix(1, 2) = -1.181847e-10;
    right_matrix(2, 0) = -9.099803e-11;
    right_matrix(2, 1) = -1.181847e-10;
    right_matrix(2, 2) = 9.864815e-01;
    // perform polar decomposition of right matrix
    Core::LinAlg::Matrix<3, 3> right_matrix_R{true};
    Core::LinAlg::Matrix<3, 3> right_matrix_U{true};
    Core::LinAlg::Matrix<3, 3> right_matrix_lambda{true};
    std::array<std::pair<double, Core::LinAlg::Matrix<3, 1>>, 3> right_matrix_spectral;
    Core::LinAlg::matrix_3x3_polar_decomposition(
        right_matrix, right_matrix_R, right_matrix_U, right_matrix_lambda, right_matrix_spectral);
    // build Q matrix
    Core::LinAlg::Matrix<3, 3> right_matrix_Q{true};
    right_matrix_Q(0, 0) = right_matrix_spectral[0].second(0);
    right_matrix_Q(0, 1) = right_matrix_spectral[0].second(1);
    right_matrix_Q(0, 2) = right_matrix_spectral[0].second(2);
    right_matrix_Q(1, 0) = right_matrix_spectral[1].second(0);
    right_matrix_Q(1, 1) = right_matrix_spectral[1].second(1);
    right_matrix_Q(1, 2) = right_matrix_spectral[1].second(2);
    right_matrix_Q(2, 0) = right_matrix_spectral[2].second(0);
    right_matrix_Q(2, 1) = right_matrix_spectral[2].second(1);
    right_matrix_Q(2, 2) = right_matrix_spectral[2].second(2);
    // get the rotation vectors
    Core::LinAlg::Matrix<3, 1> right_matrix_R_vect =
        Core::LinAlg::calc_rot_vect_from_rot_matrix(right_matrix_R);
    Core::LinAlg::Matrix<3, 1> right_matrix_Q_vect =
        Core::LinAlg::calc_rot_vect_from_rot_matrix(right_matrix_Q);

    // create interpolation parameter list and append exponential decay
    // factor
    Teuchos::ParameterList param_list;
    param_list.set("weighting_expdecay_factor", 10.0);

    // create tensor interpolator
    Core::LinAlg::TensorInterpolation::SecondOrderTensorInterpolator<1> interp{1,
        Core::LinAlg::TensorInterpolation::RotInterpType::RInterp,
        Core::LinAlg::TensorInterpolation::EigenvalInterpType::LOG, param_list};

    // get interp matrix (loc = specified)
    double loc = 5.695328e-01;
    std::vector<Core::LinAlg::Matrix<3, 3>> ref_matrices{left_matrix, right_matrix};
    std::vector<double> ref_locs{0.0, 1.0};
    Core::LinAlg::Matrix<3, 3> interp_matrix =
        interp.get_interpolated_matrix(ref_matrices, ref_locs, loc);

    // perform polar decomposition of interp matrix
    Core::LinAlg::Matrix<3, 3> interp_matrix_R{true};
    Core::LinAlg::Matrix<3, 3> interp_matrix_U{true};
    Core::LinAlg::Matrix<3, 3> interp_matrix_lambda{true};
    std::array<std::pair<double, Core::LinAlg::Matrix<3, 1>>, 3> interp_matrix_spectral;
    Core::LinAlg::matrix_3x3_polar_decomposition(interp_matrix, interp_matrix_R, interp_matrix_U,
        interp_matrix_lambda, interp_matrix_spectral);
    // build Q matrix
    Core::LinAlg::Matrix<3, 3> interp_matrix_Q{true};
    interp_matrix_Q(0, 0) = interp_matrix_spectral[0].second(0);
    interp_matrix_Q(0, 1) = interp_matrix_spectral[0].second(1);
    interp_matrix_Q(0, 2) = interp_matrix_spectral[0].second(2);
    interp_matrix_Q(1, 0) = interp_matrix_spectral[1].second(0);
    interp_matrix_Q(1, 1) = interp_matrix_spectral[1].second(1);
    interp_matrix_Q(1, 2) = interp_matrix_spectral[1].second(2);
    interp_matrix_Q(2, 0) = interp_matrix_spectral[2].second(0);
    interp_matrix_Q(2, 1) = interp_matrix_spectral[2].second(1);
    interp_matrix_Q(2, 2) = interp_matrix_spectral[2].second(2);
    // get the rotation vectors
    Core::LinAlg::Matrix<3, 1> interp_matrix_R_vect =
        Core::LinAlg::calc_rot_vect_from_rot_matrix(interp_matrix_R);
    Core::LinAlg::Matrix<3, 1> interp_matrix_Q_vect =
        Core::LinAlg::calc_rot_vect_from_rot_matrix(interp_matrix_Q);

    // define corresponding weighting during interpolation
    double c = 10.0;  // THE SAME VALUE HAS TO BE SPECIFIED IN THE INTERPOLATION ROUTINE AS WELL
    double left_unnorm = std::exp(-c * loc * loc);
    double right_unnorm = std::exp(-c * (1.0 - loc) * (1.0 - loc));
    double left_weight = left_unnorm / (left_unnorm + right_unnorm);
    double right_weight = right_unnorm / (left_unnorm + right_unnorm);

    // verify eigenvalues
    Core::LinAlg::Matrix<3, 3> lambda_ref{true};
    lambda_ref(2, 2) = std::exp(left_weight * std::log(left_matrix_spectral[0].first) +
                                right_weight * std::log(right_matrix_spectral[0].first));
    lambda_ref(1, 1) = std::exp(left_weight * std::log(left_matrix_spectral[1].first) +
                                right_weight * std::log(right_matrix_spectral[1].first));
    lambda_ref(0, 0) = std::exp(left_weight * std::log(left_matrix_spectral[2].first) +
                                right_weight * std::log(right_matrix_spectral[2].first));

    FOUR_C_EXPECT_NEAR(interp_matrix_lambda, lambda_ref, 1.0e-9);

    // compute relative rotation matrices
    Core::LinAlg::Matrix<3, 3>& Q_master = right_matrix_Q;
    Core::LinAlg::Matrix<3, 3>& R_master = right_matrix_R;
    if (loc < 0.5)
    {
      Q_master = left_matrix_Q;
      R_master = left_matrix_R;
    }
    Core::LinAlg::Matrix<3, 3> left_matrix_Q_rel{true};
    left_matrix_Q_rel.multiply_tn(1.0, Q_master, left_matrix_Q, 0.0);
    Core::LinAlg::Matrix<3, 3> right_matrix_Q_rel{true};
    right_matrix_Q_rel.multiply_tn(1.0, Q_master, right_matrix_Q, 0.0);
    Core::LinAlg::Matrix<3, 3> interp_matrix_Q_rel{true};
    interp_matrix_Q_rel.multiply_tn(1.0, Q_master, interp_matrix_Q, 0.0);
    Core::LinAlg::Matrix<3, 3> left_matrix_R_rel{true};
    left_matrix_R_rel.multiply_tn(1.0, R_master, left_matrix_R, 0.0);
    Core::LinAlg::Matrix<3, 3> right_matrix_R_rel{true};
    right_matrix_R_rel.multiply_tn(1.0, R_master, right_matrix_R, 0.0);
    Core::LinAlg::Matrix<3, 3> interp_matrix_R_rel{true};
    interp_matrix_R_rel.multiply_tn(1.0, R_master, interp_matrix_R, 0.0);
    // compute relative rotation vectors
    Core::LinAlg::Matrix<3, 1> left_matrix_Q_rel_vect =
        Core::LinAlg::calc_rot_vect_from_rot_matrix(left_matrix_Q_rel);
    Core::LinAlg::Matrix<3, 1> right_matrix_Q_rel_vect =
        Core::LinAlg::calc_rot_vect_from_rot_matrix(right_matrix_Q_rel);
    Core::LinAlg::Matrix<3, 1> interp_matrix_Q_rel_vect =
        Core::LinAlg::calc_rot_vect_from_rot_matrix(interp_matrix_Q_rel);
    Core::LinAlg::Matrix<3, 1> left_matrix_R_rel_vect =
        Core::LinAlg::calc_rot_vect_from_rot_matrix(left_matrix_R_rel);
    Core::LinAlg::Matrix<3, 1> right_matrix_R_rel_vect =
        Core::LinAlg::calc_rot_vect_from_rot_matrix(right_matrix_R_rel);
    Core::LinAlg::Matrix<3, 1> interp_matrix_R_rel_vect =
        Core::LinAlg::calc_rot_vect_from_rot_matrix(interp_matrix_R_rel);


    // verify rotation vectors
    Core::LinAlg::Matrix<3, 1> Q_rel_vect_ref{true};
    Core::LinAlg::Matrix<3, 1> R_rel_vect_ref{true};
    // build P-matrix
    Core::LinAlg::Matrix<2, 2> P{true};
    P(0, 0) = left_weight + right_weight;
    P(0, 1) = right_weight;
    P(1, 0) = right_weight;
    P(1, 1) = right_weight;
    // get inverse of the P-matrix
    Core::LinAlg::Matrix<2, 2> iP{true};
    iP.invert(P);
    // declare coefficient <a_i> vector and right-hand side <b_i> vector
    Core::LinAlg::Matrix<2, 1> a_i{true};
    Core::LinAlg::Matrix<2, 1> b_i{true};
    // build Q vector (reference)
    for (int i = 0; i < 3; ++i)
    {
      // build b_i vector
      b_i(0) = left_weight * left_matrix_Q_rel_vect(i) + right_weight * right_matrix_Q_rel_vect(i);
      b_i(1) = right_weight * right_matrix_Q_rel_vect(i);

      // get corresponding coefficients
      a_i.multiply(1.0, iP, b_i, 0.0);

      // compute corresponding rotation vector contribution (interpolation)
      Q_rel_vect_ref(i) = 1.0 * a_i(0) + loc * a_i(1);
    }
    // build R vector (reference)
    for (int i = 0; i < 3; ++i)
    {
      // build b_i vector
      b_i(0) = left_weight * left_matrix_R_rel_vect(i) + right_weight * right_matrix_R_rel_vect(i);
      b_i(1) = right_weight * right_matrix_R_rel_vect(i);

      // get corresponding coefficients
      a_i.multiply(1.0, iP, b_i, 0.0);

      // compute corresponding rotation vector contribution (interpolation)
      R_rel_vect_ref(i) = 1.0 * a_i(0) + loc * a_i(1);
    }

    FOUR_C_EXPECT_NEAR(interp_matrix_Q_rel_vect, Q_rel_vect_ref, 1.0e-9);
    FOUR_C_EXPECT_NEAR(interp_matrix_R_rel_vect, R_rel_vect_ref, 1.0e-9);
  }



}  // namespace

FOUR_C_NAMESPACE_CLOSE
