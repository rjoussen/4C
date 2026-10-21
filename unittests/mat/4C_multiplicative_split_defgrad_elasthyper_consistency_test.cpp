// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <gtest/gtest.h>

#include "4C_global_data.hpp"
#include "4C_inelastic_defgrad_factors_test_utils.hpp"
#include "4C_io_input_field.hpp"
#include "4C_io_input_parameter_container.templates.hpp"
#include "4C_linalg_symmetric_tensor.hpp"
#include "4C_linalg_tensor.hpp"
#include "4C_linalg_tensor_conversion.hpp"
#include "4C_linalg_tensor_einstein.hpp"
#include "4C_linalg_tensor_generators.hpp"
#include "4C_linalg_utils_densematrix_funct.hpp"
#include "4C_mat_inelastic_defgrad_factors_service.hpp"
#include "4C_mat_material_factory.hpp"
#include "4C_mat_multiplicative_split_defgrad_elasthyper.hpp"
#include "4C_mat_par_bundle.hpp"
#include "4C_mat_so3_material.hpp"
#include "4C_material_parameter_base.hpp"
#include "4C_solid_ele_fibers.hpp"
#include "4C_ssi_input.hpp"
#include "4C_structure_new_enum_lists.hpp"
#include "4C_unittest_utils_assertions_test.hpp"
#include "4C_utils_singleton_owner.hpp"

#include <algorithm>
#include <cmath>
#include <functional>
#include <memory>
#include <optional>
#include <string>
#include <tuple>
#include <vector>

namespace
{
  using namespace FourC;

  /// elastic summands used for the elastic part of the multiplicative split
  enum class ElasticSummands
  {
    coup_neo_hooke,
    iso_neo_hooke_vol_sussman_bathe,
    coup_neo_hooke_transv_iso
  };

  /// concentration-driven inelastic deformation gradient factor
  enum class InelasticFactor
  {
    /// F_in = lambda I with lambda^3 = 1 + g
    lin_scalar_iso,
    /// F_in = I + g n x n
    lin_scalar_aniso,
    /// F_in = F_in,no_growth F_in,lin_scalar_iso with F_in,no_growth = I
    no_growth_and_lin_scalar_iso
  };

  /// state of the concentration-driven inelastic deformation
  enum class GrowthState
  {
    /// det(F_in) != 1
    swelling,
    /// det(F_in) = 1, but dF_in/dC != 0
    reference_concentration
  };

  constexpr int matid_coup_neo_hooke = 1;
  constexpr int matid_iso_neo_hooke = 2;
  constexpr int matid_vol_sussman_bathe = 3;
  constexpr int matid_transv_iso = 4;
  constexpr int matid_structural_tensor_strategy = 5;
  constexpr int matid_lin_scalar_iso = 10;
  constexpr int matid_lin_scalar_aniso = 11;
  constexpr int matid_no_growth = 12;
  constexpr int matid_elast_hyper = 20;

  constexpr double growth_fac = 1.0e-4;
  constexpr double ref_conc = 1.0e3;
  const std::vector<double> growth_direction = {1.0, 2.0, 2.0};

  constexpr int gp = 0;
  constexpr int ele_gid = 0;

  Core::LinAlg::SymmetricTensor<double, 3, 3> green_lagrange(
      const Core::LinAlg::Tensor<double, 3, 3>& F)
  {
    return 0.5 * (Core::LinAlg::assume_symmetry(Core::LinAlg::dot(Core::LinAlg::transpose(F), F)) -
                     Core::LinAlg::TensorGenerators::identity<double, 3, 3>);
  }

  double max_abs(const auto& tensor)
  {
    return std::abs(
        std::ranges::max(tensor.container(), {}, [](const double v) { return std::abs(v); }));
  }

  double tolerance(const auto& value, const auto& reference)
  {
    return 1.0e-6 * std::max(max_abs(value), max_abs(reference)) + 1.0e-12;
  }

  /// dS_ij/dF_pq = 1/2 cmat_ijkl dC_kl/dF_pq = cmat_ijkq F_pk
  Core::LinAlg::Tensor<double, 3, 3, 3, 3> dstress_ddefgrad(
      const Core::LinAlg::SymmetricTensor<double, 3, 3, 3, 3>& cmat,
      const Core::LinAlg::Tensor<double, 3, 3>& F)
  {
    return Core::LinAlg::einsum<"ijkq", "pk">(cmat, F);
  }

  /// central finite difference of the stress w.r.t. the deformation gradient
  Core::LinAlg::Tensor<double, 3, 3, 3, 3> finite_difference_dstress_ddefgrad(
      const std::function<Core::LinAlg::SymmetricTensor<double, 3, 3>(
          const Core::LinAlg::Tensor<double, 3, 3>&)>& stress,
      const Core::LinAlg::Tensor<double, 3, 3>& F)
  {
    constexpr double h = 1.0e-6;
    Core::LinAlg::Tensor<double, 3, 3, 3, 3> dSdF_fd{};
    for (std::size_t p = 0; p < 3; ++p)
    {
      for (std::size_t q = 0; q < 3; ++q)
      {
        Core::LinAlg::Tensor<double, 3, 3> F_plus = F;
        Core::LinAlg::Tensor<double, 3, 3> F_minus = F;
        F_plus(p, q) += h;
        F_minus(p, q) -= h;

        const Core::LinAlg::SymmetricTensor<double, 3, 3> dS =
            (stress(F_plus) - stress(F_minus)) / (2.0 * h);
        for (std::size_t i = 0; i < 3; ++i)
          for (std::size_t j = 0; j < 3; ++j) dSdF_fd(i, j, p, q) = dS(i, j);
      }
    }
    return dSdF_fd;
  }

  /*!
   * Checks the multiplicative split material through its public interface only:
   *
   * - the stress equals the pull-back of the stress of an ElastHyper material with the same
   *   summands, evaluated at the elastic deformation gradient F_e = F F_in^{-1}:
   *   S = det(F_in) F_in^{-1} S_e(F_e) F_in^{-T}
   * - the elasticity tensor equals the finite difference of the stress w.r.t. F
   * - the derivative of the stress w.r.t. the concentration equals its finite difference
   *
   * The inelastic deformation gradient is driven by the growth g = growth_fac (c det(F) -
   * ref_conc), i.e. it depends on C through det(F).
   */
  class MultiplicativeSplitDefgradElastHyperConsistencyTest
      : public ::testing::TestWithParam<std::tuple<ElasticSummands, GrowthState, InelasticFactor>>
  {
   protected:
    void SetUp() override
    {
      F_ = {{{1.1, 0.01, 0.03}, {0.04, 1.2, 0.02}, {0.06, 0.05, 1.3}}};

      Global::Problem& problem = (*Global::Problem::instance());
      problem.materials()->set_read_from_problem(0);

      Core::IO::InputParameterContainer coup_neo_hooke_data;
      coup_neo_hooke_data.add("YOUNG", 1.5e2);
      coup_neo_hooke_data.add("NUE", 0.3);
      problem.materials()->insert(matid_coup_neo_hooke,
          Mat::make_parameter(matid_coup_neo_hooke, Core::Materials::MaterialType::mes_coupneohooke,
              coup_neo_hooke_data));

      Core::IO::InputParameterContainer iso_neo_hooke_data;
      iso_neo_hooke_data.add("MUE", Core::IO::InputField<double>(5.0e1));
      problem.materials()->insert(matid_iso_neo_hooke,
          Mat::make_parameter(matid_iso_neo_hooke, Core::Materials::MaterialType::mes_isoneohooke,
              iso_neo_hooke_data));

      Core::IO::InputParameterContainer vol_sussman_bathe_data;
      vol_sussman_bathe_data.add("KAPPA", 1.0e2);
      problem.materials()->insert(matid_vol_sussman_bathe,
          Mat::make_parameter(matid_vol_sussman_bathe,
              Core::Materials::MaterialType::mes_volsussmanbathe, vol_sussman_bathe_data));

      Core::IO::InputParameterContainer transv_iso_data;
      transv_iso_data.add("ALPHA", 1.0e1);
      transv_iso_data.add("BETA", 5.0);
      transv_iso_data.add("GAMMA", 3.0e1);
      transv_iso_data.add("ANGLE", 0.0);
      transv_iso_data.add("STR_TENS_ID", matid_structural_tensor_strategy);
      transv_iso_data.add<std::string>("STRATEGY", "Standard");
      transv_iso_data.add<std::string>("DISTR", "none");
      transv_iso_data.add("C1", 1.0);
      transv_iso_data.add("C2", 0.0);
      transv_iso_data.add("C3", 0.0);
      transv_iso_data.add("C4", 1.0e16);
      transv_iso_data.add("FIBER", 1);
      transv_iso_data.add("INIT", 1);
      problem.materials()->insert(matid_structural_tensor_strategy,
          Mat::make_parameter(matid_structural_tensor_strategy,
              Core::Materials::MaterialType::mes_structuraltensorstratgy, transv_iso_data));
      problem.materials()->insert(matid_transv_iso,
          Mat::make_parameter(matid_transv_iso,
              Core::Materials::MaterialType::mes_couptransverselyisotropic, transv_iso_data));

      Core::IO::InputParameterContainer lin_scalar_iso_data;
      lin_scalar_iso_data.add("SCALAR1", 1);
      lin_scalar_iso_data.add("SCALAR1_MolarGrowthFac", growth_fac);
      lin_scalar_iso_data.add("SCALAR1_RefConc", ref_conc);
      problem.materials()->insert(matid_lin_scalar_iso,
          Mat::make_parameter(matid_lin_scalar_iso,
              Core::Materials::MaterialType::mfi_lin_scalar_iso, lin_scalar_iso_data));

      problem.materials()->insert(matid_no_growth,
          Mat::make_parameter(matid_no_growth, Core::Materials::MaterialType::mfi_no_growth,
              Core::IO::InputParameterContainer{}));

      Core::IO::InputParameterContainer lin_scalar_aniso_data;
      lin_scalar_aniso_data.add("SCALAR1", 1);
      lin_scalar_aniso_data.add("SCALAR1_MolarGrowthFac", growth_fac);
      lin_scalar_aniso_data.add("SCALAR1_RefConc", ref_conc);
      lin_scalar_aniso_data.add("GrowthDirection", growth_direction);
      problem.materials()->insert(matid_lin_scalar_aniso,
          Mat::make_parameter(matid_lin_scalar_aniso,
              Core::Materials::MaterialType::mfi_lin_scalar_aniso, lin_scalar_aniso_data));

      const std::vector<int> elastic_matids = std::invoke(
          [&]() -> std::vector<int>
          {
            switch (std::get<0>(GetParam()))
            {
              case ElasticSummands::coup_neo_hooke:
                return {matid_coup_neo_hooke};
              case ElasticSummands::iso_neo_hooke_vol_sussman_bathe:
                return {matid_iso_neo_hooke, matid_vol_sussman_bathe};
              case ElasticSummands::coup_neo_hooke_transv_iso:
                return {matid_coup_neo_hooke, matid_transv_iso};
            }
            return {};
          });

      // elastic material of the multiplicative split, also evaluated separately as reference
      Core::IO::InputParameterContainer elast_hyper_data;
      elast_hyper_data.add("NUMMAT", static_cast<int>(elastic_matids.size()));
      elast_hyper_data.add("MATIDS", elastic_matids);
      elast_hyper_data.add("DENS", 1.0);
      elast_hyper_data.add("POLYCONVEX", 0);
      problem.materials()->insert(
          matid_elast_hyper, Mat::make_parameter(matid_elast_hyper,
                                 Core::Materials::MaterialType::m_elasthyper, elast_hyper_data));

      auto parameter_list_pointer = std::make_shared<Teuchos::ParameterList>();
      parameter_list_pointer->sublist("STRUCTURAL DYNAMIC", false)
          .set("MASSLIN", Solid::MassLin::ml_none);
      parameter_list_pointer->sublist("SSI CONTROL")
          .set("COUPALGO", SSI::SolutionSchemeOverFields::ssi_IterStagg);
      problem.set_parameter_list(parameter_list_pointer);

      Core::IO::InputParameterContainer multiplicative_split_data;
      multiplicative_split_data.add("NUMMATEL", static_cast<int>(elastic_matids.size()));
      multiplicative_split_data.add("MATIDSEL", elastic_matids);
      const std::vector<int> inelastic_matids = std::invoke(
          [&]() -> std::vector<int>
          {
            switch (std::get<2>(GetParam()))
            {
              case InelasticFactor::lin_scalar_iso:
                return {matid_lin_scalar_iso};
              case InelasticFactor::lin_scalar_aniso:
                return {matid_lin_scalar_aniso};
              case InelasticFactor::no_growth_and_lin_scalar_iso:
                return {matid_no_growth, matid_lin_scalar_iso};
            }
            return {};
          });
      multiplicative_split_data.add("NUMFACINEL", static_cast<int>(inelastic_matids.size()));
      multiplicative_split_data.add("INELDEFGRADFACIDS", inelastic_matids);
      multiplicative_split_data.add("DENS", 1.0);
      multiplicative_split_data.add("REF_TEMPERATURE", 293.0);
      multiplicative_split_data.add("THERMAL_EXPANSION_COEFFICIENT", 0.0);
      multiplicative_split_params_ =
          std::make_shared<Mat::PAR::MultiplicativeSplitDefgradElastHyper>(
              Core::Mat::PAR::Parameter::Data{.parameters = multiplicative_split_data});
      multiplicative_split_ = std::make_shared<Mat::MultiplicativeSplitDefgradElastHyper>(
          multiplicative_split_params_.get());

      elast_hyper_ = std::dynamic_pointer_cast<Mat::So3Material>(Mat::factory(matid_elast_hyper));
      ASSERT_NE(elast_hyper_, nullptr);

      Discret::Elements::Fibers fibers;
      fibers.element_fibers.emplace_back(
          Core::LinAlg::Tensor<double, 3>{{1.0 / 3.0, 2.0 / 3.0, 2.0 / 3.0}});
      multiplicative_split_->setup(1, fibers, std::nullopt);
      elast_hyper_->setup(1, fibers, std::nullopt);

      switch (std::get<1>(GetParam()))
      {
        case GrowthState::swelling:
          concentration_ = 2.0 * ref_conc;
          break;
        case GrowthState::reference_concentration:
          concentration_ = ref_conc / Core::LinAlg::det(F_);
          break;
      }
    }

    [[nodiscard]] static Teuchos::ParameterList evaluation_parameters(const double concentration)
    {
      Teuchos::ParameterList params;
      params.set<std::shared_ptr<std::vector<double>>>(
          "scalars", std::make_shared<std::vector<double>>(std::vector<double>{concentration}));
      return params;
    }

    void evaluate(const Core::LinAlg::Tensor<double, 3, 3>& F, const double concentration,
        Core::LinAlg::SymmetricTensor<double, 3, 3>& stress,
        Core::LinAlg::SymmetricTensor<double, 3, 3, 3, 3>& cmat) const
    {
      const Teuchos::ParameterList params = evaluation_parameters(concentration);
      stress = {};
      cmat = {};
      multiplicative_split_->evaluate(
          &F, green_lagrange(F), params, context_, stress, cmat, gp, ele_gid);
    }

    [[nodiscard]] Core::LinAlg::SymmetricTensor<double, 3, 3> evaluate_stress(
        const Core::LinAlg::Tensor<double, 3, 3>& F, const double concentration) const
    {
      Core::LinAlg::SymmetricTensor<double, 3, 3> stress{};
      Core::LinAlg::SymmetricTensor<double, 3, 3, 3, 3> cmat{};
      evaluate(F, concentration, stress, cmat);
      return stress;
    }

    /// inverse inelastic deformation gradient
    [[nodiscard]] Core::LinAlg::Tensor<double, 3, 3> inverse_inelastic_defgrad() const
    {
      const double growth = growth_fac * (concentration_ * Core::LinAlg::det(F_) - ref_conc);
      const auto id =
          Core::LinAlg::get_full(Core::LinAlg::TensorGenerators::identity<double, 3, 3>);
      if (std::get<2>(GetParam()) != InelasticFactor::lin_scalar_aniso)
        return (1.0 / std::cbrt(1.0 + growth)) * id;

      const Core::LinAlg::Tensor<double, 3> n = {
          {growth_direction[0], growth_direction[1], growth_direction[2]}};
      return Core::LinAlg::inv(id + (growth / (n * n)) * Core::LinAlg::dyadic(n, n));
    }

    Core::LinAlg::Tensor<double, 3, 3> F_;
    double concentration_ = 0.0;
    const Mat::EvaluationContext<3> context_{};

    std::shared_ptr<Mat::PAR::MultiplicativeSplitDefgradElastHyper> multiplicative_split_params_;
    std::shared_ptr<Mat::MultiplicativeSplitDefgradElastHyper> multiplicative_split_;
    std::shared_ptr<Mat::So3Material> elast_hyper_;

    Core::Utils::SingletonOwnerRegistry::ScopeGuard guard_;
  };

  TEST_P(MultiplicativeSplitDefgradElastHyperConsistencyTest, StressEqualsPulledBackElasticStress)
  {
    const Core::LinAlg::Tensor<double, 3, 3> iFin = inverse_inelastic_defgrad();
    const double detFin = 1.0 / Core::LinAlg::det(iFin);

    // elastic stress from the ElastHyper material at F_e = F F_in^{-1}
    const Core::LinAlg::Tensor<double, 3, 3> Fe = Core::LinAlg::dot(F_, iFin);
    const Teuchos::ParameterList params = evaluation_parameters(concentration_);
    Core::LinAlg::SymmetricTensor<double, 3, 3> Se{};
    Core::LinAlg::SymmetricTensor<double, 3, 3, 3, 3> cmat_e{};
    elast_hyper_->evaluate(&Fe, green_lagrange(Fe), params, context_, Se, cmat_e, gp, ele_gid);

    const Core::LinAlg::SymmetricTensor<double, 3, 3> S_ref =
        detFin * Core::LinAlg::assume_symmetry(
                     Core::LinAlg::dot(Core::LinAlg::dot(iFin, Se), Core::LinAlg::transpose(iFin)));

    const Core::LinAlg::SymmetricTensor<double, 3, 3> S = evaluate_stress(F_, concentration_);

    FOUR_C_EXPECT_NEAR(S, S_ref, tolerance(S, S_ref));
  }

  TEST_P(MultiplicativeSplitDefgradElastHyperConsistencyTest, CmatEqualsFiniteDifference)
  {
    Core::LinAlg::SymmetricTensor<double, 3, 3> stress{};
    Core::LinAlg::SymmetricTensor<double, 3, 3, 3, 3> cmat{};
    evaluate(F_, concentration_, stress, cmat);

    const Core::LinAlg::Tensor<double, 3, 3, 3, 3> dSdF = dstress_ddefgrad(cmat, F_);
    const Core::LinAlg::Tensor<double, 3, 3, 3, 3> dSdF_fd =
        finite_difference_dstress_ddefgrad([&](const Core::LinAlg::Tensor<double, 3, 3>& F)
            { return evaluate_stress(F, concentration_); }, F_);

    FOUR_C_EXPECT_NEAR(dSdF, dSdF_fd, tolerance(dSdF, dSdF_fd));
  }

  TEST_P(MultiplicativeSplitDefgradElastHyperConsistencyTest,
      DStressDConcentrationEqualsFiniteDifference)
  {
    Teuchos::ParameterList params = evaluation_parameters(concentration_);
    params.set<int>("differentiationtype", static_cast<int>(Solid::DifferentiationType::elch));
    const Core::LinAlg::SymmetricTensor<double, 3, 3> dSdc =
        multiplicative_split_->evaluate_d_stress_d_scalar(
            F_, green_lagrange(F_), params, context_, gp, ele_gid);

    const double h = 1.0e-6 * concentration_;
    const Core::LinAlg::SymmetricTensor<double, 3, 3> dSdc_fd =
        (evaluate_stress(F_, concentration_ + h) - evaluate_stress(F_, concentration_ - h)) /
        (2.0 * h);

    FOUR_C_EXPECT_NEAR(dSdc, dSdc_fd, tolerance(dSdc, dSdc_fd));
  }

  INSTANTIATE_TEST_SUITE_P(ElasticSummandsAndGrowthStates,
      MultiplicativeSplitDefgradElastHyperConsistencyTest,
      ::testing::Combine(::testing::Values(ElasticSummands::coup_neo_hooke,
                             ElasticSummands::iso_neo_hooke_vol_sussman_bathe,
                             ElasticSummands::coup_neo_hooke_transv_iso),
          ::testing::Values(GrowthState::swelling, GrowthState::reference_concentration),
          ::testing::Values(InelasticFactor::lin_scalar_iso, InelasticFactor::lin_scalar_aniso,
              InelasticFactor::no_growth_and_lin_scalar_iso)),
      [](const auto& info)
      {
        const std::string summands = std::invoke(
            [&]() -> std::string
            {
              switch (std::get<0>(info.param))
              {
                case ElasticSummands::coup_neo_hooke:
                  return "CoupNeoHooke";
                case ElasticSummands::iso_neo_hooke_vol_sussman_bathe:
                  return "IsoNeoHookeVolSussmanBathe";
                case ElasticSummands::coup_neo_hooke_transv_iso:
                  return "CoupNeoHookeTransvIso";
              }
              return "";
            });
        const std::string growth = std::get<1>(info.param) == GrowthState::swelling
                                       ? "Swelling"
                                       : "ReferenceConcentration";
        const std::string factor = std::invoke(
            [&]() -> std::string
            {
              switch (std::get<2>(info.param))
              {
                case InelasticFactor::lin_scalar_iso:
                  return "Iso";
                case InelasticFactor::lin_scalar_aniso:
                  return "Aniso";
                case InelasticFactor::no_growth_and_lin_scalar_iso:
                  return "NoGrowthAndIso";
              }
              return "";
            });
        return summands + "_" + growth + "_" + factor;
      });
  namespace ViscoplastUtils = Mat::InelasticDefgradTransvIsotropElastViscoplastUtils;

  /// state of the viscoplastic inelastic deformation
  enum class ViscoplasticState
  {
    /// stress below the yield stress, no inelastic deformation
    elastic,
    /// viscoplastic flow
    plastic
  };

  /// elastic summands used for the elastic part of the viscoplastic material
  enum class ViscoplasticElasticSummands
  {
    coup_neo_hooke,
    coup_neo_hooke_transv_iso
  };

  /// thermal expansion of the viscoplastic material
  enum class ThermalExpansion
  {
    off,
    on
  };

  constexpr int vp_matid_coup_neo_hooke = 1;
  constexpr int vp_matid_transv_iso = 2;
  constexpr int vp_matid_fiber_reader = 3;
  constexpr int vp_matid_structural_tensor_strategy = 4;
  constexpr int vp_matid_viscoplastic_law = 5;
  constexpr int vp_matid_viscoplastic_factor = 6;
  constexpr int vp_matid_elast_hyper = 7;
  constexpr int vp_matid_multiplicative_split = 8;

  /*!
   * Checks the multiplicative split material with the transversely isotropic elastic viscoplastic
   * inelastic factor through its public interface only:
   *
   * - the inelastic deformation is zero in the elastic state and nonzero in the plastic state
   * - the elasticity tensor equals the finite difference of the stress w.r.t. F
   * - the derivative of the stress w.r.t. the temperature equals its finite difference
   * - the derivatives of the Taylor-Quinney heat source w.r.t. C and T equal their finite
   *   differences
   *
   * The history is not updated, i.e., all evaluations start from the initial state.
   */
  class MultiplicativeSplitDefgradElastHyperViscoplasticConsistencyTest
      : public ::testing::TestWithParam<std::tuple<ViscoplastUtils::MatBehavior, ViscoplasticState,
            ViscoplasticElasticSummands, ThermalExpansion>>
  {
   protected:
    // declared first to be destroyed last, i.e., after the materials referencing its parameters
    Core::Utils::SingletonOwnerRegistry::ScopeGuard guard_;

    void SetUp() override
    {
      const auto [mat_behavior, state, elastic_summands, thermal_expansion] = GetParam();
      if (thermal_expansion == ThermalExpansion::on) thermal_expansion_coefficient_ = 1.0e-5;

      Global::Problem& problem = (*Global::Problem::instance());
      problem.materials()->set_read_from_problem(0);

      Core::IO::InputParameterContainer coup_neo_hooke_data;
      coup_neo_hooke_data.add("YOUNG", 2.0e5);
      coup_neo_hooke_data.add("NUE", 0.29);
      problem.materials()->insert(vp_matid_coup_neo_hooke,
          Mat::make_parameter(vp_matid_coup_neo_hooke,
              Core::Materials::MaterialType::mes_coupneohooke, coup_neo_hooke_data));

      const auto transv_iso_data = [](const double alpha, const double beta, const double gamma)
      {
        Core::IO::InputParameterContainer data;
        data.add("ALPHA", alpha);
        data.add("BETA", beta);
        data.add("GAMMA", gamma);
        data.add("ANGLE", 0.0);
        data.add("STR_TENS_ID", vp_matid_structural_tensor_strategy);
        data.add<std::string>("STRATEGY", "Standard");
        data.add<std::string>("DISTR", "none");
        data.add("C1", 1.0);
        data.add("C2", 0.0);
        data.add("C3", 0.0);
        data.add("C4", 1.0e16);
        data.add("FIBER", 1);
        data.add("INIT", 1);
        return data;
      };
      problem.materials()->insert(vp_matid_structural_tensor_strategy,
          Mat::make_parameter(vp_matid_structural_tensor_strategy,
              Core::Materials::MaterialType::mes_structuraltensorstratgy,
              transv_iso_data(1.0, 1.0, 1.0)));
      problem.materials()->insert(
          vp_matid_transv_iso, Mat::make_parameter(vp_matid_transv_iso,
                                   Core::Materials::MaterialType::mes_couptransverselyisotropic,
                                   transv_iso_data(1.0e4, 5.0e3, 2.0e4)));
      problem.materials()->insert(
          vp_matid_fiber_reader, Mat::make_parameter(vp_matid_fiber_reader,
                                     Core::Materials::MaterialType::mes_couptransverselyisotropic,
                                     transv_iso_data(1.0, 1.0, 1.0)));

      Core::IO::InputParameterContainer viscoplastic_law_data;
      viscoplastic_law_data.add("STRAIN_RATE_PREFAC", 1.0);
      viscoplastic_law_data.add("STRAIN_RATE_EXP_FAC", 0.014);
      viscoplastic_law_data.add("INIT_YIELD_STRENGTH", 792.0);
      viscoplastic_law_data.add("ISOTROP_HARDEN_PREFAC", 510.0);
      viscoplastic_law_data.add("ISOTROP_HARDEN_EXP", 0.26);
      viscoplastic_law_data.add("REF_TEMPERATURE", ref_temperature_);
      viscoplastic_law_data.add("MELT_TEMPERATURE", 1793.0);
      viscoplastic_law_data.add("TEMPERATURE_SENS", 1.03);
      problem.materials()->insert(vp_matid_viscoplastic_law,
          Mat::make_parameter(vp_matid_viscoplastic_law,
              Core::Materials::MaterialType::mvl_reformulated_Johnson_Cook, viscoplastic_law_data));

      Core::IO::InputParameterContainer viscoplastic_factor_data;
      viscoplastic_factor_data.add("FIBER_READER_ID", vp_matid_fiber_reader);
      viscoplastic_factor_data.add("TAYLOR_QUINNEY_COEFFICIENT", 0.85);
      viscoplastic_factor_data.add("LINEARIZATION", ViscoplastUtils::LinearizationType::analytic);
      viscoplastic_factor_data.add(
          "MATRIX_EXP_CALC_METHOD", Core::LinAlg::MatrixExpCalcMethod::automatic);
      viscoplastic_factor_data.add("MATRIX_EXP_DERIV_CALC_METHOD",
          Core::LinAlg::GenMatrixExpFirstDerivCalcMethod::automatic);
      viscoplastic_factor_data.add(
          "MATRIX_LOG_CALC_METHOD", Core::LinAlg::MatrixLogCalcMethod::inv_scal_square);
      viscoplastic_factor_data.add("MATRIX_LOG_DERIV_CALC_METHOD",
          Core::LinAlg::GenMatrixLogFirstDerivCalcMethod::pade_part_fract);
      viscoplastic_factor_data.add("MAT_BEHAVIOR", mat_behavior);
      viscoplastic_factor_data.add("MAX_PLASTIC_STRAIN_DERIV_INCR", std::exp(30.0));
      viscoplastic_factor_data.add("MAX_PLASTIC_STRAIN_INCR", std::exp(30.0));
      viscoplastic_factor_data.add(
          "TIME_INTEGRATION_HIST_VARS", ViscoplastUtils::TimIntType::logarithmic);
      viscoplastic_factor_data.add("VISCOPLAST_LAW_ID", vp_matid_viscoplastic_law);
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
      problem.materials()->insert(vp_matid_viscoplastic_factor,
          Mat::make_parameter(vp_matid_viscoplastic_factor,
              Core::Materials::MaterialType::mfi_transv_isotrop_elast_viscoplast,
              viscoplastic_factor_data));

      const std::vector<int> elastic_matids =
          elastic_summands == ViscoplasticElasticSummands::coup_neo_hooke
              ? std::vector<int>{vp_matid_coup_neo_hooke}
              : std::vector<int>{vp_matid_coup_neo_hooke, vp_matid_transv_iso};
      Core::IO::InputParameterContainer elast_hyper_data;
      elast_hyper_data.add("NUMMAT", static_cast<int>(elastic_matids.size()));
      elast_hyper_data.add("MATIDS", elastic_matids);
      elast_hyper_data.add("DENS", 1.0);
      elast_hyper_data.add("POLYCONVEX", 0);
      problem.materials()->insert(
          vp_matid_elast_hyper, Mat::make_parameter(vp_matid_elast_hyper,
                                    Core::Materials::MaterialType::m_elasthyper, elast_hyper_data));

      // the viscoplastic factor finds its parent material in the material bundle
      Core::IO::InputParameterContainer multiplicative_split_data;
      multiplicative_split_data.add("NUMMATEL", static_cast<int>(elastic_matids.size()));
      multiplicative_split_data.add("MATIDSEL", elastic_matids);
      multiplicative_split_data.add("NUMFACINEL", 1);
      multiplicative_split_data.add(
          "INELDEFGRADFACIDS", std::vector<int>{vp_matid_viscoplastic_factor});
      multiplicative_split_data.add("DENS", 1.0);
      multiplicative_split_data.add("REF_TEMPERATURE", ref_temperature_);
      multiplicative_split_data.add(
          "THERMAL_EXPANSION_COEFFICIENT", thermal_expansion_coefficient_);
      problem.materials()->insert(vp_matid_multiplicative_split,
          Mat::make_parameter(vp_matid_multiplicative_split,
              Core::Materials::MaterialType::m_multiplicative_split_defgrad_elasthyper,
              multiplicative_split_data));

      auto parameter_list_pointer = std::make_shared<Teuchos::ParameterList>();
      parameter_list_pointer->sublist("STRUCTURAL DYNAMIC", false)
          .set("MASSLIN", Solid::MassLin::ml_none);
      parameter_list_pointer->sublist("SSI CONTROL")
          .set("COUPALGO", SSI::SolutionSchemeOverFields::ssi_IterStagg);
      problem.set_parameter_list(parameter_list_pointer);

      multiplicative_split_ = std::dynamic_pointer_cast<Mat::MultiplicativeSplitDefgradElastHyper>(
          Mat::factory(vp_matid_multiplicative_split));
      ASSERT_NE(multiplicative_split_, nullptr);
      elast_hyper_ =
          std::dynamic_pointer_cast<Mat::So3Material>(Mat::factory(vp_matid_elast_hyper));
      ASSERT_NE(elast_hyper_, nullptr);

      Discret::Elements::Fibers fibers;
      fibers.element_fibers.emplace_back(
          Core::LinAlg::Tensor<double, 3>{{1.0 / 3.0, 2.0 / 3.0, 2.0 / 3.0}});
      multiplicative_split_->setup(1, fibers, std::nullopt);
      elast_hyper_->setup(1, fibers, std::nullopt);

      switch (state)
      {
        case ViscoplasticState::elastic:
          F_ = {{{1.0005, 0.00005, 0.00015}, {0.0002, 1.001, 0.0001}, {0.0003, 0.00025, 0.9985}}};
          break;
        case ViscoplasticState::plastic:
        {
          // simple shear such that the equivalent stress is moderately above the yield stress
          const double shear =
              mat_behavior == ViscoplastUtils::MatBehavior::isotropic ? 0.0065 : 0.004;
          F_ = {{{1.0, shear, 0.0005}, {0.0002, 1.0003, 0.0001}, {0.0003, 0.0002, 0.9998}}};
          break;
        }
      }
    }

    [[nodiscard]] static Teuchos::ParameterList evaluation_parameters(const double temperature)
    {
      Teuchos::ParameterList params;
      params.set<double>("temperature", temperature);
      return params;
    }

    [[nodiscard]] Core::LinAlg::SymmetricTensor<double, 3, 3> evaluate_stress(
        const Core::LinAlg::Tensor<double, 3, 3>& F, const double temperature) const
    {
      Core::LinAlg::SymmetricTensor<double, 3, 3> stress{};
      Core::LinAlg::SymmetricTensor<double, 3, 3, 3, 3> cmat{};
      multiplicative_split_->evaluate(&F, green_lagrange(F), evaluation_parameters(temperature),
          context_, stress, cmat, gp, ele_gid);
      return stress;
    }

    const double ref_temperature_ = 293.0;
    const double temperature_ = 313.0;
    double thermal_expansion_coefficient_ = 0.0;
    double total_time_ = 1.0e-6;
    double time_step_size_ = 1.0e-6;
    const Mat::EvaluationContext<3> context_{.total_time = &total_time_,
        .time_step_size = &time_step_size_,
        .xi = nullptr,
        .ref_coords = nullptr};
    Core::LinAlg::Tensor<double, 3, 3> F_;

    std::shared_ptr<Mat::MultiplicativeSplitDefgradElastHyper> multiplicative_split_;
    std::shared_ptr<Mat::So3Material> elast_hyper_;
  };

  TEST_P(MultiplicativeSplitDefgradElastHyperViscoplasticConsistencyTest,
      InelasticDeformationMatchesState)
  {
    const Core::LinAlg::SymmetricTensor<double, 3, 3> S = evaluate_stress(F_, temperature_);

    // purely elastic response, i.e., F_in = I: S = S_e(C) - S_e(C_T) with C_T = (1 + 2 alpha dT) I
    const auto elastic_stress = [&](const Core::LinAlg::Tensor<double, 3, 3>& F)
    {
      Core::LinAlg::SymmetricTensor<double, 3, 3> Se{};
      Core::LinAlg::SymmetricTensor<double, 3, 3, 3, 3> cmat_e{};
      elast_hyper_->evaluate(&F, green_lagrange(F), evaluation_parameters(temperature_), context_,
          Se, cmat_e, gp, ele_gid);
      return Se;
    };
    const Core::LinAlg::Tensor<double, 3, 3> FT =
        std::sqrt(1.0 + 2.0 * thermal_expansion_coefficient_ * (temperature_ - ref_temperature_)) *
        Core::LinAlg::get_full(Core::LinAlg::TensorGenerators::identity<double, 3, 3>);
    const Core::LinAlg::SymmetricTensor<double, 3, 3> S_elastic =
        elastic_stress(F_) - elastic_stress(FT);

    const double deviation = max_abs(Core::LinAlg::SymmetricTensor<double, 3, 3>(S - S_elastic));
    if (std::get<1>(GetParam()) == ViscoplasticState::elastic)
      EXPECT_LE(deviation, 1.0e-10 * max_abs(S_elastic));
    else
      EXPECT_GE(deviation, 1.0e-2 * max_abs(S_elastic));
  }

  TEST_P(
      MultiplicativeSplitDefgradElastHyperViscoplasticConsistencyTest, CmatEqualsFiniteDifference)
  {
    Core::LinAlg::SymmetricTensor<double, 3, 3> stress{};
    Core::LinAlg::SymmetricTensor<double, 3, 3, 3, 3> cmat{};
    multiplicative_split_->evaluate(&F_, green_lagrange(F_), evaluation_parameters(temperature_),
        context_, stress, cmat, gp, ele_gid);

    const Core::LinAlg::Tensor<double, 3, 3, 3, 3> dSdF = dstress_ddefgrad(cmat, F_);
    const Core::LinAlg::Tensor<double, 3, 3, 3, 3> dSdF_fd =
        finite_difference_dstress_ddefgrad([&](const Core::LinAlg::Tensor<double, 3, 3>& F)
            { return evaluate_stress(F, temperature_); }, F_);

    FOUR_C_EXPECT_NEAR(dSdF, dSdF_fd, tolerance(dSdF, dSdF_fd));
  }

  TEST_P(MultiplicativeSplitDefgradElastHyperViscoplasticConsistencyTest,
      DStressDTemperatureEqualsFiniteDifference)
  {
    Teuchos::ParameterList params = evaluation_parameters(temperature_);
    params.set<int>("differentiationtype", static_cast<int>(Solid::DifferentiationType::temp));
    const Core::LinAlg::SymmetricTensor<double, 3, 3> dSdT =
        multiplicative_split_->evaluate_d_stress_d_scalar(
            F_, green_lagrange(F_), params, context_, gp, ele_gid);

    constexpr double h = 1.0e-3;
    const Core::LinAlg::SymmetricTensor<double, 3, 3> dSdT_fd =
        (evaluate_stress(F_, temperature_ + h) - evaluate_stress(F_, temperature_ - h)) / (2.0 * h);

    FOUR_C_EXPECT_NEAR(dSdT, dSdT_fd, tolerance(dSdT, dSdT_fd));
  }

  TEST_P(MultiplicativeSplitDefgradElastHyperViscoplasticConsistencyTest,
      HeatSourceDerivativesEqualFiniteDifference)
  {
    const auto heat_source = [&](const Core::LinAlg::Tensor<double, 3, 3>& F, const double T)
    {
      const Core::LinAlg::Matrix<3, 3> FM = Core::LinAlg::make_matrix(F);
      return multiplicative_split_->evaluate_additional_heat_source(context_, gp, ele_gid, &FM, T);
    };
    const Mat::HeatSource Q = heat_source(F_, temperature_);

    // dQ/dF_pq = dQ/dC_kl dC_kl/dF_pq = 2 F_pk dQ/dC_kq
    Core::LinAlg::Matrix<6, 1> dQdC_voigt;
    dQdC_voigt.update_t(Q.derivative_wrt_cauchy_green);
    const Core::LinAlg::Tensor<double, 3, 3> dQdF =
        2.0 * Core::LinAlg::einsum<"pk", "kq">(F_,
                  Core::LinAlg::make_symmetric_tensor_from_stress_like_voigt_matrix(dQdC_voigt));

    // the heat source is strongly nonlinear, a small step keeps the truncation error low
    constexpr double h = 1.0e-7;
    Core::LinAlg::Tensor<double, 3, 3> dQdF_fd{};
    for (std::size_t p = 0; p < 3; ++p)
    {
      for (std::size_t q = 0; q < 3; ++q)
      {
        Core::LinAlg::Tensor<double, 3, 3> F_plus = F_;
        Core::LinAlg::Tensor<double, 3, 3> F_minus = F_;
        F_plus(p, q) += h;
        F_minus(p, q) -= h;
        dQdF_fd(p, q) =
            (heat_source(F_plus, temperature_).value - heat_source(F_minus, temperature_).value) /
            (2.0 * h);
      }
    }
    FOUR_C_EXPECT_NEAR(dQdF, dQdF_fd, tolerance(dQdF, dQdF_fd));

    constexpr double h_T = 1.0e-3;
    const double dQdT_fd =
        (heat_source(F_, temperature_ + h_T).value - heat_source(F_, temperature_ - h_T).value) /
        (2.0 * h_T);
    EXPECT_NEAR(Q.derivative_wrt_temperature, dQdT_fd,
        1.0e-6 * std::max(std::abs(Q.derivative_wrt_temperature), std::abs(dQdT_fd)) + 1.0e-12);
  }

  INSTANTIATE_TEST_SUITE_P(MatBehaviorsStatesAndElasticSummands,
      MultiplicativeSplitDefgradElastHyperViscoplasticConsistencyTest,
      ::testing::Combine(::testing::Values(ViscoplastUtils::MatBehavior::isotropic,
                             ViscoplastUtils::MatBehavior::transv_isotropic),
          ::testing::Values(ViscoplasticState::elastic, ViscoplasticState::plastic),
          ::testing::Values(ViscoplasticElasticSummands::coup_neo_hooke,
              ViscoplasticElasticSummands::coup_neo_hooke_transv_iso),
          ::testing::Values(ThermalExpansion::off, ThermalExpansion::on)),
      [](const auto& info)
      {
        const std::string behavior =
            std::get<0>(info.param) == ViscoplastUtils::MatBehavior::isotropic ? "Isotropic"
                                                                               : "TransvIsotropic";
        const std::string state =
            std::get<1>(info.param) == ViscoplasticState::elastic ? "Elastic" : "Plastic";
        const std::string summands =
            std::get<2>(info.param) == ViscoplasticElasticSummands::coup_neo_hooke
                ? "CoupNeoHooke"
                : "CoupNeoHookeTransvIso";
        const std::string thermal_expansion =
            std::get<3>(info.param) == ThermalExpansion::on ? "_ThermalExpansion" : "";
        return behavior + "_" + state + "_" + summands + thermal_expansion;
      });
}  // namespace
