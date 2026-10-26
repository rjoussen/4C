// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <gtest/gtest.h>

#include "4C_inelastic_defgrad_factors_test_utils.hpp"
#include "4C_io_input_parameter_container.templates.hpp"
#include "4C_linalg_serialdensematrix.hpp"
#include "4C_linalg_tensor_conversion.hpp"
#include "4C_linalg_utils_densematrix_funct.hpp"
#include "4C_mat_inelastic_defgrad_factors_service.hpp"
#include "4C_mat_multiplicative_split_defgrad_elasthyper.hpp"
#include "4C_multiplicative_split_defgrad_elasthyper_test_utils.hpp"
#include "4C_structure_new_enum_lists.hpp"
#include "4C_unittest_utils_assertions_test.hpp"
#include "4C_utils_enum.hpp"
#include "4C_utils_singleton_owner.hpp"

#include <cmath>
#include <format>
#include <memory>
#include <optional>
#include <string>
#include <tuple>
#include <vector>

namespace
{
  using namespace FourC;
  namespace MultSplitTestUtils = MultiplicativeSplitDefgradElastHyperTestUtils;
  namespace ViscoplastUtils = Mat::InelasticDefgradTransvIsotropElastViscoplastUtils;

  /// state of the Gauss point
  enum class PlasticState
  {
    elastic,  //< stress below the yield stress, no inelastic deformation
    plastic   //< viscoplastic flow
  };

  /// coupling of the material to the temperature
  enum class TemperatureCoupling
  {
    /// no temperature in the parameter list and no thermal expansion (like solid-only simulations)
    without_temperature,
    /// temperature different from the reference temperature and thermal expansion
    with_temperature
  };

  constexpr int matid_fiber_reader = MultSplitTestUtils::first_free_matid;
  constexpr int matid_viscoplastic_law = MultSplitTestUtils::first_free_matid + 1;
  constexpr int matid_viscoplastic_factor = MultSplitTestUtils::first_free_matid + 2;
  constexpr int matid_multiplicative_split = MultSplitTestUtils::first_free_matid + 3;

  constexpr double ref_temperature = 293.0;
  constexpr double temperature = 313.0;
  constexpr double time_step_size = 1.0e-6;

  /// set up the multiplicative split with an ElastHyper material with the given elastic summands,
  /// the viscoplastic factor with the given material behavior and the thermal expansion of the
  /// given temperature coupling
  void set_up_materials(const MultSplitTestUtils::ElasticSummands elastic_summands,
      const ViscoplastUtils::MatBehavior mat_behavior,
      const TemperatureCoupling temperature_coupling)
  {
    MultSplitTestUtils::set_up_elast_hyper(elastic_summands);

    // only provides the fiber direction of the viscoplastic factor
    Global::Problem::instance()->materials()->insert(
        matid_fiber_reader, Mat::make_parameter(matid_fiber_reader,
                                Core::Materials::MaterialType::mes_couptransverselyisotropic,
                                MultSplitTestUtils::transv_iso_data(1.0, 1.0, 1.0,
                                    MultSplitTestUtils::matid_structural_tensor_strategy)));

    Core::IO::InputParameterContainer viscoplastic_law_data;
    viscoplastic_law_data.add("STRAIN_RATE_PREFAC", 1.0);
    viscoplastic_law_data.add("STRAIN_RATE_EXP_FAC", 0.014);
    viscoplastic_law_data.add("INIT_YIELD_STRENGTH", 792.0);
    viscoplastic_law_data.add("ISOTROP_HARDEN_PREFAC", 510.0);
    viscoplastic_law_data.add("ISOTROP_HARDEN_EXP", 0.26);
    viscoplastic_law_data.add("REF_TEMPERATURE", ref_temperature);
    viscoplastic_law_data.add("MELT_TEMPERATURE", 1793.0);
    viscoplastic_law_data.add("TEMPERATURE_SENS", 1.03);
    Global::Problem::instance()->materials()->insert(matid_viscoplastic_law,
        Mat::make_parameter(matid_viscoplastic_law,
            Core::Materials::MaterialType::mvl_reformulated_Johnson_Cook, viscoplastic_law_data));

    Core::IO::InputParameterContainer viscoplastic_factor_data;
    viscoplastic_factor_data.add("FIBER_READER_ID", matid_fiber_reader);
    viscoplastic_factor_data.add("TAYLOR_QUINNEY_COEFFICIENT", 0.85);
    viscoplastic_factor_data.add("LINEARIZATION", ViscoplastUtils::LinearizationType::analytic);
    viscoplastic_factor_data.add(
        "MATRIX_EXP_CALC_METHOD", Core::LinAlg::MatrixExpCalcMethod::automatic);
    viscoplastic_factor_data.add(
        "MATRIX_EXP_DERIV_CALC_METHOD", Core::LinAlg::GenMatrixExpFirstDerivCalcMethod::automatic);
    viscoplastic_factor_data.add(
        "MATRIX_LOG_CALC_METHOD", Core::LinAlg::MatrixLogCalcMethod::inv_scal_square);
    viscoplastic_factor_data.add("MATRIX_LOG_DERIV_CALC_METHOD",
        Core::LinAlg::GenMatrixLogFirstDerivCalcMethod::pade_part_fract);
    viscoplastic_factor_data.add("MAT_BEHAVIOR", mat_behavior);
    viscoplastic_factor_data.add(
        "TIME_INTEGRATION_HIST_VARS", ViscoplastUtils::TimIntType::logarithmic);
    viscoplastic_factor_data.add("VISCOPLAST_LAW_ID", matid_viscoplastic_law);
    viscoplastic_factor_data.add<std::optional<double>>("YIELD_COND_A", 1.0);
    viscoplastic_factor_data.add<std::optional<double>>("YIELD_COND_B", 2.0);
    viscoplastic_factor_data.add<std::optional<double>>("YIELD_COND_F", 2.5);
    viscoplastic_factor_data.group("LOCAL_SUBSTEPPING").add("USE_SUBSTEPPING", false);
    viscoplastic_factor_data.group("LOCAL_SUBSTEPPING").add("MAX_SUBSTEPPING_HALVE_NUM", 0);
    viscoplastic_factor_data.add("LOCAL_NEWTON",
        ViscoplastUtils::LocalNewtonParams{
            .res_tol = 1.0e-14,
            .incr_tol = 1.0e-14,
            .conv_check = ViscoplastUtils::LocalNewtonConvCheck::residual_and_increment_ratio,
            .diver_cont = ViscoplastUtils::LocalNewtonDiverCont::stop,
            .max_iter = 100,
            .max_exceedance_fact_res_tol = 1.0e1,
            .max_exceedance_fact_incr_tol = 1.0e1,
        });
    viscoplastic_factor_data.add("ERROR_REGISTRATION_SETTINGS",
        ViscoplastUtils::ErrorRegistrationSettings{.register_plastic_strain_incr_overflow = true,
            .max_plastic_strain_incr = std::exp(30.0),
            .register_plastic_strain_deriv_incr_overflow = false,
            .max_plastic_strain_deriv_incr = std::exp(30.0)});
    viscoplastic_factor_data.add(
        "ADAPTIVE_ESTIMATE_INTERPOLATION", InelasticDefgradFactorsTestUtils::set_up_aei_params(
                                               {.use_adaptive_estimate_interpolation = false}));
    Global::Problem::instance()->materials()->insert(matid_viscoplastic_factor,
        Mat::make_parameter(matid_viscoplastic_factor,
            Core::Materials::MaterialType::mfi_transv_isotrop_elast_viscoplast,
            viscoplastic_factor_data));

    Core::IO::InputParameterContainer multiplicative_split_data;
    multiplicative_split_data.add("MATID_ELAST", MultSplitTestUtils::matid_elast_hyper);
    multiplicative_split_data.add("NUMFACINEL", 1);
    multiplicative_split_data.add("INELDEFGRADFACIDS", std::vector<int>{matid_viscoplastic_factor});
    multiplicative_split_data.add("REF_TEMPERATURE", ref_temperature);
    multiplicative_split_data.add("THERMAL_EXPANSION_COEFFICIENT",
        temperature_coupling == TemperatureCoupling::with_temperature ? 1.0e-5 : 0.0);
    Global::Problem::instance()->materials()->insert(matid_multiplicative_split,
        Mat::make_parameter(matid_multiplicative_split,
            Core::Materials::MaterialType::m_multiplicative_split_defgrad_elasthyper,
            multiplicative_split_data));

    MultSplitTestUtils::set_up_generic_parameters();
  }

  /// Create a parameter list that carries the temperature, if given, and the correct
  /// differentiation type
  Teuchos::ParameterList create_parameter_list_with_temperature(
      const std::optional<double> temperature)
  {
    Teuchos::ParameterList params;
    if (temperature.has_value()) params.set<double>("temperature", *temperature);
    // only read by evaluate_d_stress_d_scalar()
    params.set<int>("differentiationtype", static_cast<int>(Solid::DifferentiationType::temp));
    return params;
  }

  /// deformation gradient of the given state and material behavior
  Core::LinAlg::Tensor<double, 3, 3> get_deformation_gradient(
      const ViscoplastUtils::MatBehavior mat_behavior, const PlasticState state)
  {
    switch (state)
    {
      case PlasticState::elastic:
        return {{{1.0005, 0.00005, 0.00015}, {0.0002, 1.001, 0.0001}, {0.0003, 0.00025, 0.9985}}};
      case PlasticState::plastic:
      {
        // shear entry such that the equivalent stress is moderately above the yield stress
        const double shear_entry =
            mat_behavior == ViscoplastUtils::MatBehavior::isotropic ? 0.0065 : 0.004;
        return {{{1.0, shear_entry, 0.0005}, {0.0002, 1.0003, 0.0001}, {0.0003, 0.0002, 0.9998}}};
      }
      default:
        FOUR_C_THROW("Unexpected plastic state.");
    }
  }

  /*!
   * Checks the multiplicative split material with the transversely isotropic elastic viscoplastic
   * inelastic factor through its public interface only:
   *
   * - the local Newton iteration only runs in the plastic state
   * - the stiffness tensor equals the finite difference of the stress w.r.t. F
   *
   * The history is not updated, i.e., all evaluations start from the initial state.
   */
  class MultiplicativeSplitDefgradElastHyperViscoplastTest
      : public ::testing::TestWithParam<std::tuple<MultSplitTestUtils::ElasticSummands,
            ViscoplastUtils::MatBehavior, PlasticState, TemperatureCoupling>>
  {
   protected:
    Core::Utils::SingletonOwnerRegistry::ScopeGuard guard_;

    void SetUp() override
    {
      const auto [elastic_summands, mat_behavior, plastic_state, temperature_coupling] = GetParam();
      set_up_materials(elastic_summands, mat_behavior, temperature_coupling);
      multiplicative_split_ =
          MultSplitTestUtils::create_material<Mat::MultiplicativeSplitDefgradElastHyper>(
              matid_multiplicative_split);
      elast_hyper_ = MultSplitTestUtils::create_material<Mat::So3Material>(
          MultSplitTestUtils::matid_elast_hyper);

      plastic_state_ = plastic_state;
      F_ = get_deformation_gradient(mat_behavior, plastic_state);
      if (temperature_coupling == TemperatureCoupling::with_temperature) temperature_ = temperature;
    }

    /// evaluate the stress of the multiplicative split material at F and the given temperature
    [[nodiscard]] Core::LinAlg::SymmetricTensor<double, 3, 3> evaluate_stress(
        const Core::LinAlg::Tensor<double, 3, 3>& F, const std::optional<double> temperature) const
    {
      return MultSplitTestUtils::evaluate(
          *multiplicative_split_, F, create_parameter_list_with_temperature(temperature), context_)
          .stress;
    }

    const Mat::EvaluationContext<3> context_{.total_time = &time_step_size,
        .time_step_size = &time_step_size,
        .xi = nullptr,
        .ref_coords = nullptr};
    PlasticState plastic_state_{};
    Core::LinAlg::Tensor<double, 3, 3> F_;
    /// temperature in the parameter list, none without temperature coupling
    std::optional<double> temperature_;

    std::shared_ptr<Mat::MultiplicativeSplitDefgradElastHyper> multiplicative_split_;
    std::shared_ptr<Mat::So3Material> elast_hyper_;
  };

  TEST_P(MultiplicativeSplitDefgradElastHyperViscoplastTest, LocalNewtonRunsOnlyInPlasticState)
  {
    [[maybe_unused]] const auto S = evaluate_stress(F_, temperature_);

    Core::LinAlg::SerialDenseMatrix local_newton_iterations(1, 1);
    multiplicative_split_->evaluate_output_data("local_newton_iters", local_newton_iterations);
    if (plastic_state_ == PlasticState::elastic)
      EXPECT_EQ(local_newton_iterations(0, 0), 0.0);
    else
      EXPECT_GT(local_newton_iterations(0, 0), 0.0);
  }

  TEST_P(MultiplicativeSplitDefgradElastHyperViscoplastTest, CmatEqualsFiniteDifference)
  {
    const Teuchos::ParameterList params = create_parameter_list_with_temperature(temperature_);
    const Core::LinAlg::SymmetricTensor<double, 3, 3, 3, 3> cmat =
        MultSplitTestUtils::evaluate(*multiplicative_split_, F_, params, context_).cmat;

    // analytic
    const Core::LinAlg::Tensor<double, 3, 3, 3, 3> dSdF =
        MultSplitTestUtils::compute_stress_deriv_wrt_defgrad(cmat, F_);

    // finite difference
    const auto evaluate_stress_at_fixed_temperature =
        [&](const Core::LinAlg::Tensor<double, 3, 3>& F)
    { return evaluate_stress(F, temperature_); };
    const Core::LinAlg::Tensor<double, 3, 3, 3, 3> dSdF_fd =
        MultSplitTestUtils::finite_difference_stress_deriv_wrt_defgrad(
            evaluate_stress_at_fixed_temperature, F_);

    FOUR_C_EXPECT_NEAR_RELATIVE(dSdF, dSdF_fd, 1.0e-6);
  }

  INSTANTIATE_TEST_SUITE_P(ElasticSummandsMatBehaviorsAndStates,
      MultiplicativeSplitDefgradElastHyperViscoplastTest,
      ::testing::Combine(
          ::testing::ValuesIn(EnumTools::enum_values<MultSplitTestUtils::ElasticSummands>()),
          ::testing::ValuesIn(EnumTools::enum_values<ViscoplastUtils::MatBehavior>()),
          ::testing::ValuesIn(EnumTools::enum_values<PlasticState>()),
          ::testing::ValuesIn(EnumTools::enum_values<TemperatureCoupling>())),
      [](const auto& info)
      {
        return std::format("{}_{}_{}_{}", EnumTools::enum_name(std::get<0>(info.param)),
            EnumTools::enum_name(std::get<1>(info.param)),
            EnumTools::enum_name(std::get<2>(info.param)),
            EnumTools::enum_name(std::get<3>(info.param)));
      });

  /*!
   * The same setup, instantiated with temperature coupling only. Checks additionally:
   *
   * - the derivative of the stress w.r.t. the temperature equals its finite difference
   * - the derivatives of the Taylor-Quinney heat source w.r.t. C and T equal their finite
   *   differences
   */
  class MultiplicativeSplitDefgradElastHyperViscoplastWithTemperatureTest
      : public MultiplicativeSplitDefgradElastHyperViscoplastTest
  {
  };

  TEST_P(MultiplicativeSplitDefgradElastHyperViscoplastWithTemperatureTest,
      DStressDTemperatureEqualsFiniteDifference)
  {
    // analytic
    const Core::LinAlg::SymmetricTensor<double, 3, 3> dSdT =
        multiplicative_split_->evaluate_d_stress_d_scalar(F_,
            MultSplitTestUtils::green_lagrange(F_),
            create_parameter_list_with_temperature(temperature_), context_, MultSplitTestUtils::gp,
            MultSplitTestUtils::ele_gid);

    // finite difference
    constexpr double h = 1.0e-3;
    const Core::LinAlg::SymmetricTensor<double, 3, 3> dSdT_fd =
        (evaluate_stress(F_, *temperature_ + h) - evaluate_stress(F_, *temperature_ - h)) /
        (2.0 * h);

    FOUR_C_EXPECT_NEAR_RELATIVE(dSdT, dSdT_fd, 1.0e-6);
  }

  TEST_P(MultiplicativeSplitDefgradElastHyperViscoplastWithTemperatureTest,
      HeatSourceDerivativesEqualFiniteDifference)
  {
    const auto heat_source = [&](const Core::LinAlg::Tensor<double, 3, 3>& F, const double T)
    {
      const Core::LinAlg::Matrix<3, 3> FM = Core::LinAlg::make_matrix(F);
      return multiplicative_split_->evaluate_additional_heat_source(
          context_, MultSplitTestUtils::gp, MultSplitTestUtils::ele_gid, &FM, T);
    };
    const Mat::HeatSource Q = heat_source(F_, *temperature_);

    // analytic: dQ/dF_kl = dQ/dC_mn dC_mn/dF_kl = 2 F_km dQ/dC_ml
    Core::LinAlg::Matrix<6, 1> dQdC_voigt;
    dQdC_voigt.update_t(Q.derivative_wrt_cauchy_green);
    const Core::LinAlg::Tensor<double, 3, 3> dQdF =
        2.0 * Core::LinAlg::einsum<"km", "ml">(F_,
                  Core::LinAlg::make_symmetric_tensor_from_stress_like_voigt_matrix(dQdC_voigt));

    // finite difference, the heat source is strongly nonlinear, a small step keeps the truncation
    // error low
    constexpr double h = 1.0e-7;
    Core::LinAlg::Tensor<double, 3, 3> dQdF_fd{};
    for (std::size_t k = 0; k < 3; ++k)
    {
      for (std::size_t l = 0; l < 3; ++l)
      {
        Core::LinAlg::Tensor<double, 3, 3> F_plus = F_;
        Core::LinAlg::Tensor<double, 3, 3> F_minus = F_;
        F_plus(k, l) += h;
        F_minus(k, l) -= h;
        dQdF_fd(k, l) =
            (heat_source(F_plus, *temperature_).value - heat_source(F_minus, *temperature_).value) /
            (2.0 * h);
      }
    }

    constexpr double h_T = 1.0e-3;
    const double dQdT_fd =
        (heat_source(F_, *temperature_ + h_T).value - heat_source(F_, *temperature_ - h_T).value) /
        (2.0 * h_T);
    const double dQdT_analytic = Q.derivative_wrt_temperature;

    // without plastic flow, the heat source and its derivatives vanish
    if (plastic_state_ == PlasticState::elastic)
    {
      FOUR_C_EXPECT_NEAR(dQdF, dQdF_fd, 1.0e-12);
      EXPECT_NEAR(dQdT_analytic, dQdT_fd, 1.0e-12);
    }
    else
    {
      FOUR_C_EXPECT_NEAR_RELATIVE(dQdF, dQdF_fd, 1.0e-6);
      EXPECT_NEAR(dQdT_analytic, dQdT_fd, 1.0e-6 * std::abs(dQdT_fd));
    }
  }

  INSTANTIATE_TEST_SUITE_P(ElasticSummandsMatBehaviorsAndStates,
      MultiplicativeSplitDefgradElastHyperViscoplastWithTemperatureTest,
      ::testing::Combine(
          ::testing::ValuesIn(EnumTools::enum_values<MultSplitTestUtils::ElasticSummands>()),
          ::testing::ValuesIn(EnumTools::enum_values<ViscoplastUtils::MatBehavior>()),
          ::testing::ValuesIn(EnumTools::enum_values<PlasticState>()),
          ::testing::Values(TemperatureCoupling::with_temperature)),
      [](const auto& info)
      {
        return std::format("{}_{}_{}_{}", EnumTools::enum_name(std::get<0>(info.param)),
            EnumTools::enum_name(std::get<1>(info.param)),
            EnumTools::enum_name(std::get<2>(info.param)),
            EnumTools::enum_name(std::get<3>(info.param)));
      });

  /// the same setup, instantiated for the elastic state only
  class MultiplicativeSplitDefgradElastHyperViscoplastElasticStateTest
      : public MultiplicativeSplitDefgradElastHyperViscoplastTest
  {
  };

  /// In the elastic state without temperature, i.e., without inelastic and thermal deformation, the
  /// stress and the stiffness tensor equal those of an ElastHyper material with the same summands.
  TEST_P(
      MultiplicativeSplitDefgradElastHyperViscoplastElasticStateTest, StressAndCmatEqualElastHyper)
  {
    const Teuchos::ParameterList params = create_parameter_list_with_temperature(std::nullopt);
    const MultSplitTestUtils::StressResponse response =
        MultSplitTestUtils::evaluate(*multiplicative_split_, F_, params, context_);
    const MultSplitTestUtils::StressResponse response_elast_hyper =
        MultSplitTestUtils::evaluate(*elast_hyper_, F_, params, context_);

    FOUR_C_EXPECT_NEAR_RELATIVE(response.stress, response_elast_hyper.stress, 1.0e-12);
    FOUR_C_EXPECT_NEAR_RELATIVE(response.cmat, response_elast_hyper.cmat, 1.0e-12);
  }

  INSTANTIATE_TEST_SUITE_P(ElasticSummandsAndMatBehaviors,
      MultiplicativeSplitDefgradElastHyperViscoplastElasticStateTest,
      ::testing::Combine(
          ::testing::ValuesIn(EnumTools::enum_values<MultSplitTestUtils::ElasticSummands>()),
          ::testing::ValuesIn(EnumTools::enum_values<ViscoplastUtils::MatBehavior>()),
          ::testing::Values(PlasticState::elastic),
          ::testing::Values(TemperatureCoupling::without_temperature)),
      [](const auto& info)
      {
        return std::format("{}_{}_{}_{}", EnumTools::enum_name(std::get<0>(info.param)),
            EnumTools::enum_name(std::get<1>(info.param)),
            EnumTools::enum_name(std::get<2>(info.param)),
            EnumTools::enum_name(std::get<3>(info.param)));
      });
}  // namespace
