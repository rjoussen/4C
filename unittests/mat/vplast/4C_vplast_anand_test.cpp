// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <gtest/gtest.h>

#include "4C_global_data.hpp"
#include "4C_linalg_fixedsizematrix.hpp"
#include "4C_mat_material_factory.hpp"
#include "4C_mat_vplast_anand.hpp"
#include "4C_mat_vplast_law.hpp"
#include "4C_unittest_utils_assertions_test.hpp"
#include "4C_utils_exceptions.hpp"
#include "4C_utils_singleton_owner.hpp"

#include <memory>
namespace
{
  using namespace FourC;

  class AnandTest : public ::testing::Test
  {
   protected:
    void SetUp() override
    {
      // set up equivalent stress and plastic strain
      equiv_stress_ = 1.0;
      equiv_plastic_strain_ = 1.0;

      // manually create viscoplastic law (Anand, parameters from Anand et al. 2019)
      Core::IO::InputParameterContainer vplast_law_Anand_data;
      vplast_law_Anand_data.add("STRAIN_RATE_PREFAC", 0.013888191563123096);
      vplast_law_Anand_data.add("STRAIN_RATE_SENS", 0.015);
      vplast_law_Anand_data.add("INIT_FLOW_RES", 0.95);
      vplast_law_Anand_data.add("HARDEN_RATE_PREFAC", 10.0);
      vplast_law_Anand_data.add("HARDEN_RATE_SENS", 2.0);
      vplast_law_Anand_data.add("FLOW_RES_SAT_FAC", 2.0);
      vplast_law_Anand_data.add("FLOW_RES_SAT_EXP", 0.05);
      params_vplast_law_Anand_ = std::dynamic_pointer_cast<Mat::Viscoplastic::PAR::Anand>(
          std::shared_ptr(Mat::make_parameter(
              1, Core::Materials::MaterialType::mvl_Anand, vplast_law_Anand_data)));
      vplast_law_Anand_ =
          std::make_shared<Mat::Viscoplastic::Anand>(params_vplast_law_Anand_.get());

      // define setup parameter for InelasticDefGradTransvIsotropElastViscoplast
      Core::IO::InputParameterContainer setup_vplast_law_Anand;  // can stay empty

      // call setup method for Anand
      int numgp = 8;  // HEX8 element, although not really relevant for the tested methods
      vplast_law_Anand_->setup(numgp, setup_vplast_law_Anand);

      // call pre_evaluate
      vplast_law_Anand_->pre_evaluate(0);
    }

    // equivalent stress
    double equiv_stress_;
    // equivalent stress
    double equiv_plastic_strain_;
    // reference solution for the stress ratio (Anand)
    double stress_ratio_Anand_solution_;
    // reference solution for the plastic strain rate (Anand)
    double plastic_strain_rate_Anand_solution_;
    // reference solution for the plastic strain rate derivatives, w.r.t. equivalent stress and
    // plastic strain (Anand)
    Core::LinAlg::Matrix<2, 1> deriv_plastic_strain_rate_Anand_solution_;
    // pointer to Anand
    std::shared_ptr<Mat::Viscoplastic::Anand> vplast_law_Anand_;
    // pointer to parameters of Anand
    std::shared_ptr<Mat::Viscoplastic::PAR::Anand> params_vplast_law_Anand_;

    Core::Utils::SingletonOwnerRegistry::ScopeGuard guard;
  };

  TEST_F(AnandTest, TestEvaluateStressRatio)
  {
    // set reference solution
    stress_ratio_Anand_solution_ = 0.8813854742856903;

    // compute solution from the viscoplasticity law
    double stress_ratio_Anand =
        vplast_law_Anand_->evaluate_stress_ratio(equiv_stress_, equiv_plastic_strain_);

    // compare solutions
    EXPECT_NEAR(stress_ratio_Anand_solution_, stress_ratio_Anand, 1.0e-10);
  }

  TEST_F(AnandTest, TestEvaluatePlasticStrainRate)
  {
    // set reference solution
    plastic_strain_rate_Anand_solution_ = 3.069292557971517e-06;

    // declare error status and overflow check boolean
    Mat::ViscoplastErrorType err_status = Mat::ViscoplastErrorType::NoErrors;

    // compute solution from the viscoplasticity law
    double plastic_strain_rate_Anand = vplast_law_Anand_->evaluate_plastic_strain_rate(
        equiv_stress_, equiv_plastic_strain_, 1.0, 1.0e30, err_status, true);

    if (err_status != Mat::ViscoplastErrorType::NoErrors)
      FOUR_C_THROW("Error encountered during testing of TestEvaluatePlasticStrainRate");


    // compare solutions
    EXPECT_NEAR(plastic_strain_rate_Anand_solution_, plastic_strain_rate_Anand, 1.0e-10);
  }

  TEST_F(AnandTest, TestEvaluatePlasticStrainRateDerivatives)
  {
    // set reference solution
    deriv_plastic_strain_rate_Anand_solution_(0, 0) = 6.301050522594067e-05;
    deriv_plastic_strain_rate_Anand_solution_(1, 0) = -3.33954141454838e-06;

    // declare error status and overflow check boolean
    Mat::ViscoplastErrorType err_status = Mat::ViscoplastErrorType::NoErrors;

    // call method for plastic strain rate evaluation in order to update the history variables, and
    // make the material ready for the derivative evaluation
    vplast_law_Anand_->evaluate_plastic_strain_rate(
        equiv_stress_, equiv_plastic_strain_, 1.0, 1.0e30, err_status, true);
    if (err_status != Mat::ViscoplastErrorType::NoErrors)
      FOUR_C_THROW("Error encountered during testing of TestEvaluatePlasticStrainRate");



    // compute solution from the viscoplasticity law
    Core::LinAlg::Matrix<2, 1> deriv_plastic_strain_rate_Anand =
        vplast_law_Anand_->evaluate_derivatives_of_plastic_strain_rate(
            equiv_stress_, equiv_plastic_strain_, 1.0, 1.0e30, err_status, false);

    if (err_status != Mat::ViscoplastErrorType::NoErrors)
      FOUR_C_THROW("Error encountered during testing of TestEvaluatePlasticStrainRateDerivatives");

    // compare solutions
    FOUR_C_EXPECT_NEAR(
        deriv_plastic_strain_rate_Anand_solution_, deriv_plastic_strain_rate_Anand, 1.0e-10);
  }

}  // namespace
