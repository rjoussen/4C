// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_config.hpp"

#include "4C_global_data.hpp"
#include "4C_io_input_parameter_container.templates.hpp"
#include "4C_linalg_serialdensematrix.hpp"
#include "4C_linalg_symmetric_tensor.hpp"
#include "4C_linalg_tensor.hpp"
#include "4C_linalg_tensor_generators.hpp"
#include "4C_linalg_utils_densematrix_funct.hpp"
#include "4C_mat_inelastic_defgrad_factors_service.hpp"
#include "4C_mat_material_factory.hpp"
#include "4C_mat_multiplicative_split_defgrad_elasthyper.hpp"
#include "4C_mat_par_bundle.hpp"
#include "4C_solid_ele_fibers.hpp"
#include "4C_ssi_input.hpp"
#include "4C_structure_new_enum_lists.hpp"
#include "4C_utils_exceptions.hpp"
#include "4C_utils_singleton_owner.hpp"

#include <benchmark/benchmark.h>
#include <Teuchos_ParameterList.hpp>

#include <cmath>
#include <memory>
#include <optional>
#include <string>
#include <vector>

FOUR_C_NAMESPACE_OPEN

namespace
{
  namespace ViscoplastUtils = Mat::InelasticDefgradTransvIsotropElastViscoplastUtils;
  namespace AEI = ViscoplastUtils::AdaptiveEstimateInterpolation;

  /// state of the Gauss point
  enum class ViscoplasticState
  {
    /// stress below the yield stress, no inelastic deformation
    elastic,
    /// viscoplastic flow
    plastic
  };

  /// coupling of the material to the temperature
  enum class TemperatureCoupling
  {
    /// no temperature in the parameter list and no thermal expansion (structure only)
    off,
    /// temperature equal to the reference temperature and thermal expansion, i.e., the same
    /// mechanical state as without coupling plus the cost of the temperature coupling
    on
  };

  /// adaptive estimate interpolation of the starting point of the local Newton
  enum class AdaptiveEstimateInterpolation
  {
    off,
    on
  };

  /// evaluated material functions
  enum class Evaluation
  {
    /// stress and elasticity tensor
    stress,
    /// stress and elasticity tensor followed by the derivative of the stress w.r.t. the
    /// temperature at the same state, as in a monolithic TSI iteration
    stress_and_temperature_derivative
  };

  constexpr int matid_coup_neo_hooke = 1;
  constexpr int matid_fiber_reader = 2;
  constexpr int matid_structural_tensor_strategy = 3;
  constexpr int matid_viscoplastic_law = 4;
  constexpr int matid_viscoplastic_factor = 5;
  constexpr int matid_multiplicative_split = 6;

  constexpr double ref_temperature = 293.0;
  constexpr double time_step_size = 1.0e-6;

  constexpr int gp = 0;
  constexpr int ele_gid = 0;

  void set_up_materials_in_global_problem(const ViscoplastUtils::MatBehavior mat_behavior,
      const AdaptiveEstimateInterpolation adaptive_estimate_interpolation,
      const double thermal_expansion_coefficient)
  {
    Global::Problem& problem = *Global::Problem::instance();
    problem.materials()->set_read_from_problem(0);

    Core::IO::InputParameterContainer coup_neo_hooke_data;
    coup_neo_hooke_data.add("YOUNG", 2.0e5);
    coup_neo_hooke_data.add("NUE", 0.29);
    problem.materials()->insert(matid_coup_neo_hooke,
        Mat::make_parameter(matid_coup_neo_hooke, Core::Materials::MaterialType::mes_coupneohooke,
            coup_neo_hooke_data));

    Core::IO::InputParameterContainer fiber_reader_data;
    fiber_reader_data.add("ALPHA", 1.0);
    fiber_reader_data.add("BETA", 1.0);
    fiber_reader_data.add("GAMMA", 1.0);
    fiber_reader_data.add("ANGLE", 0.0);
    fiber_reader_data.add("STR_TENS_ID", matid_structural_tensor_strategy);
    fiber_reader_data.add<std::string>("STRATEGY", "Standard");
    fiber_reader_data.add<std::string>("DISTR", "none");
    fiber_reader_data.add("C1", 1.0);
    fiber_reader_data.add("C2", 0.0);
    fiber_reader_data.add("C3", 0.0);
    fiber_reader_data.add("C4", 1.0e16);
    fiber_reader_data.add("FIBER", 1);
    fiber_reader_data.add("INIT", 1);
    problem.materials()->insert(matid_structural_tensor_strategy,
        Mat::make_parameter(matid_structural_tensor_strategy,
            Core::Materials::MaterialType::mes_structuraltensorstratgy, fiber_reader_data));
    problem.materials()->insert(matid_fiber_reader,
        Mat::make_parameter(matid_fiber_reader,
            Core::Materials::MaterialType::mes_couptransverselyisotropic, fiber_reader_data));

    Core::IO::InputParameterContainer viscoplastic_law_data;
    viscoplastic_law_data.add("STRAIN_RATE_PREFAC", 1.0);
    viscoplastic_law_data.add("STRAIN_RATE_EXP_FAC", 0.014);
    viscoplastic_law_data.add("INIT_YIELD_STRENGTH", 792.0);
    viscoplastic_law_data.add("ISOTROP_HARDEN_PREFAC", 510.0);
    viscoplastic_law_data.add("ISOTROP_HARDEN_EXP", 0.26);
    viscoplastic_law_data.add("REF_TEMPERATURE", ref_temperature);
    viscoplastic_law_data.add("MELT_TEMPERATURE", 1793.0);
    viscoplastic_law_data.add("TEMPERATURE_SENS", 1.03);
    problem.materials()->insert(matid_viscoplastic_law,
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
    viscoplastic_factor_data.add("MAX_PLASTIC_STRAIN_DERIV_INCR", std::exp(30.0));
    viscoplastic_factor_data.add("MAX_PLASTIC_STRAIN_INCR", std::exp(30.0));
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
    viscoplastic_factor_data.add("ADAPTIVE_ESTIMATE_INTERPOLATION",
        AEI::AEIParams{
            .use_adaptive_estimate_interpolation =
                adaptive_estimate_interpolation == AdaptiveEstimateInterpolation::on,
            .elastic_predictor_zero_component_threshold = 1.0e-13,
            .plastic_predictor_construction =
                {
                    .elastic_stretch_eigenval_type =
                        AEI::PrelimPlasticPredictor::ElasticStretchEigenvalType::scale_unit,
                    .elastic_stretch_eigenvect_type = AEI::PrelimPlasticPredictor::
                        ElasticStretchEigenvectType::from_elastic_predictor,
                    .elastic_rotation_type =
                        AEI::PrelimPlasticPredictor::ElasticRotationType::from_elastic_predictor,
                    .max_iter = 50,
                    .relative_understress_tol = 1.0e-6,
                    .interval_scanning_param = 0.5,
                },
            .estimate_interpolation =
                {
                    .starting_point_type = AEI::StartingPointType::equiv_stress_history,
                    .user_set_starting_point = 0.5,
                    .max_iter = 50,
                    .interval_scanning_param = 0.5,
                },
            .hardening =
                {
                    .method = AEI::HardeningManagementMethod::integrate_via_evolution_equations,
                    .max_iter_integration = 50,
                    .tol_integration = 1.0e-8,
                },
            .reestimation =
                {
                    .max_num_reestimations = 10,
                    .interval_scanning_param = 0.5,
                },
        });
    problem.materials()->insert(matid_viscoplastic_factor,
        Mat::make_parameter(matid_viscoplastic_factor,
            Core::Materials::MaterialType::mfi_transv_isotrop_elast_viscoplast,
            viscoplastic_factor_data));

    // the viscoplastic factor finds its parent material in the material bundle
    Core::IO::InputParameterContainer multiplicative_split_data;
    multiplicative_split_data.add("NUMMATEL", 1);
    multiplicative_split_data.add("MATIDSEL", std::vector<int>{matid_coup_neo_hooke});
    multiplicative_split_data.add("NUMFACINEL", 1);
    multiplicative_split_data.add("INELDEFGRADFACIDS", std::vector<int>{matid_viscoplastic_factor});
    multiplicative_split_data.add("DENS", 1.0);
    multiplicative_split_data.add("REF_TEMPERATURE", ref_temperature);
    multiplicative_split_data.add("THERMAL_EXPANSION_COEFFICIENT", thermal_expansion_coefficient);
    problem.materials()->insert(matid_multiplicative_split,
        Mat::make_parameter(matid_multiplicative_split,
            Core::Materials::MaterialType::m_multiplicative_split_defgrad_elasthyper,
            multiplicative_split_data));

    auto parameter_list = std::make_shared<Teuchos::ParameterList>();
    parameter_list->sublist("STRUCTURAL DYNAMIC", false).set("MASSLIN", Solid::MassLin::ml_none);
    parameter_list->sublist("SSI CONTROL")
        .set("COUPALGO", SSI::SolutionSchemeOverFields::ssi_IterStagg);
    problem.set_parameter_list(parameter_list);
  }

  /// create @p num_materials independent instances, each with a single Gauss point
  std::vector<std::shared_ptr<Mat::MultiplicativeSplitDefgradElastHyper>> create_materials(
      const std::size_t num_materials)
  {
    Discret::Elements::Fibers fibers;
    fibers.element_fibers.emplace_back(
        Core::LinAlg::Tensor<double, 3>{{1.0 / 3.0, 2.0 / 3.0, 2.0 / 3.0}});

    std::vector<std::shared_ptr<Mat::MultiplicativeSplitDefgradElastHyper>> materials;
    materials.reserve(num_materials);
    for (std::size_t i = 0; i < num_materials; ++i)
    {
      auto material = std::dynamic_pointer_cast<Mat::MultiplicativeSplitDefgradElastHyper>(
          Mat::factory(matid_multiplicative_split));
      material->setup(1, fibers, std::nullopt);
      materials.emplace_back(std::move(material));
    }
    return materials;
  }

  Core::LinAlg::Tensor<double, 3, 3> deformation_gradient(
      const ViscoplasticState state, const ViscoplastUtils::MatBehavior mat_behavior)
  {
    switch (state)
    {
      case ViscoplasticState::elastic:
        return {{{1.0005, 0.00005, 0.00015}, {0.0002, 1.001, 0.0001}, {0.0003, 0.00025, 0.9985}}};
      case ViscoplasticState::plastic:
      {
        // simple shear such that the equivalent stress is moderately above the yield stress
        const double shear =
            mat_behavior == ViscoplastUtils::MatBehavior::isotropic ? 0.0065 : 0.004;
        return {{{1.0, shear, 0.0005}, {0.0002, 1.0003, 0.0001}, {0.0003, 0.0002, 0.9998}}};
      }
    }
    return {};
  }

  /// local Newton iterations of the benchmarked state; a change alters the benchmarked workload
  int expected_local_newton_iterations(const ViscoplasticState state,
      const ViscoplastUtils::MatBehavior mat_behavior,
      const AdaptiveEstimateInterpolation adaptive_estimate_interpolation)
  {
    if (state == ViscoplasticState::elastic) return 0;
    if (adaptive_estimate_interpolation == AdaptiveEstimateInterpolation::on) return 7;
    return mat_behavior == ViscoplastUtils::MatBehavior::isotropic ? 17 : 14;
  }

  Core::LinAlg::SymmetricTensor<double, 3, 3> green_lagrange_strain(
      const Core::LinAlg::Tensor<double, 3, 3>& F)
  {
    return 0.5 * (Core::LinAlg::assume_symmetry(Core::LinAlg::dot(Core::LinAlg::transpose(F), F)) -
                     Core::LinAlg::TensorGenerators::identity<double, 3, 3>);
  }

  /*!
   * Evaluates the multiplicative split with the viscoplastic factor at a fixed state.
   *
   * The factor caches its inverse inelastic deformation gradient and the thermo-mechanical
   * coupling quantities, i.e., a repeated evaluation at the same state would skip the local
   * Newton. Hence, each iteration evaluates a fresh material instance once, as in a new Newton
   * iteration of the structure. The instances are taken from a small pool, which is refilled
   * outside of the timing, to keep the memory independent of the number of iterations.
   */
  template <ViscoplasticState viscoplastic_state, TemperatureCoupling temperature_coupling,
      Evaluation evaluation,
      ViscoplastUtils::MatBehavior mat_behavior = ViscoplastUtils::MatBehavior::isotropic,
      AdaptiveEstimateInterpolation adaptive_estimate_interpolation =
          AdaptiveEstimateInterpolation::off>
  void evaluate_viscoplastic_multiplicative_split(benchmark::State& state)
  {
    Core::Utils::SingletonOwnerRegistry::ScopeGuard guard;
    set_up_materials_in_global_problem(mat_behavior, adaptive_estimate_interpolation,
        temperature_coupling == TemperatureCoupling::on ? 1.0e-5 : 0.0);

    double total_time = time_step_size;
    double dt = time_step_size;
    const Mat::EvaluationContext<3> context{
        .total_time = &total_time, .time_step_size = &dt, .xi = nullptr, .ref_coords = nullptr};

    Teuchos::ParameterList params;
    if (temperature_coupling == TemperatureCoupling::on)
      params.set<double>("temperature", ref_temperature);
    Teuchos::ParameterList params_temperature_derivative = params;
    params_temperature_derivative.set<int>(
        "differentiationtype", static_cast<int>(Solid::DifferentiationType::temp));

    const Core::LinAlg::Tensor<double, 3, 3> F =
        deformation_gradient(viscoplastic_state, mat_behavior);
    const Core::LinAlg::SymmetricTensor<double, 3, 3> E = green_lagrange_strain(F);
    Core::LinAlg::SymmetricTensor<double, 3, 3> stress{};
    Core::LinAlg::SymmetricTensor<double, 3, 3, 3, 3> cmat{};

    // check the state on a separate instance via the number of local Newton iterations
    {
      const auto material = create_materials(1).front();
      material->evaluate(&F, E, params, context, stress, cmat, gp, ele_gid);
      Core::LinAlg::SerialDenseMatrix local_newton_iterations(1, 1);
      material->evaluate_output_data("local_newton_iters", local_newton_iterations);
      const int iterations = static_cast<int>(local_newton_iterations(0, 0));
      const int expected_iterations = expected_local_newton_iterations(
          viscoplastic_state, mat_behavior, adaptive_estimate_interpolation);
      // fail hard (also in the dry run of the regular tests) instead of reporting a wrong workload
      FOUR_C_ASSERT_ALWAYS(iterations == expected_iterations,
          "Expected {} local Newton iterations, got {}.", expected_iterations, iterations);
    }

    constexpr std::size_t pool_size = 64;
    auto materials = create_materials(pool_size);
    auto material = materials.begin();
    for (auto _ : state)
    {
      if (material == materials.end())
      {
        state.PauseTiming();
        materials = create_materials(pool_size);
        material = materials.begin();
        state.ResumeTiming();
      }

      (*material)->evaluate(&F, E, params, context, stress, cmat, gp, ele_gid);
      benchmark::DoNotOptimize(stress);
      benchmark::DoNotOptimize(cmat);
      if constexpr (evaluation == Evaluation::stress_and_temperature_derivative)
      {
        Core::LinAlg::SymmetricTensor<double, 3, 3> dstress_dT =
            (*material)->evaluate_d_stress_d_scalar(
                F, E, params_temperature_derivative, context, gp, ele_gid);
        benchmark::DoNotOptimize(dstress_dT);
      }
      ++material;
    }
  }
}  // namespace

BENCHMARK(evaluate_viscoplastic_multiplicative_split<ViscoplasticState::elastic,
    TemperatureCoupling::off, Evaluation::stress>);
BENCHMARK(evaluate_viscoplastic_multiplicative_split<ViscoplasticState::elastic,
    TemperatureCoupling::on, Evaluation::stress>);
BENCHMARK(evaluate_viscoplastic_multiplicative_split<ViscoplasticState::plastic,
    TemperatureCoupling::off, Evaluation::stress>);
BENCHMARK(evaluate_viscoplastic_multiplicative_split<ViscoplasticState::plastic,
    TemperatureCoupling::on, Evaluation::stress>);
BENCHMARK(evaluate_viscoplastic_multiplicative_split<ViscoplasticState::plastic,
    TemperatureCoupling::on, Evaluation::stress_and_temperature_derivative>);
BENCHMARK(evaluate_viscoplastic_multiplicative_split<ViscoplasticState::elastic,
    TemperatureCoupling::on, Evaluation::stress_and_temperature_derivative>);
BENCHMARK(evaluate_viscoplastic_multiplicative_split<ViscoplasticState::plastic,
    TemperatureCoupling::off, Evaluation::stress, ViscoplastUtils::MatBehavior::transv_isotropic>);
BENCHMARK(evaluate_viscoplastic_multiplicative_split<ViscoplasticState::plastic,
    TemperatureCoupling::on, Evaluation::stress, ViscoplastUtils::MatBehavior::transv_isotropic>);
BENCHMARK(evaluate_viscoplastic_multiplicative_split<ViscoplasticState::plastic,
    TemperatureCoupling::off, Evaluation::stress, ViscoplastUtils::MatBehavior::isotropic,
    AdaptiveEstimateInterpolation::on>);

FOUR_C_NAMESPACE_CLOSE
