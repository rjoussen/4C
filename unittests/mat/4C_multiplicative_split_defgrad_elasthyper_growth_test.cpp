// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <gtest/gtest.h>

#include "4C_linalg_tensor.hpp"
#include "4C_mat_multiplicative_split_defgrad_elasthyper.hpp"
#include "4C_multiplicative_split_defgrad_elasthyper_test_utils.hpp"
#include "4C_structure_new_enum_lists.hpp"
#include "4C_unittest_utils_assertions_test.hpp"
#include "4C_utils_enum.hpp"
#include "4C_utils_singleton_owner.hpp"

#include <format>
#include <memory>
#include <string>
#include <tuple>
#include <vector>

namespace
{
  using namespace FourC;
  namespace MultSplitTestUtils = MultiplicativeSplitDefgradElastHyperTestUtils;

  /// inelastic factor configurations used in the test
  enum class InelasticFactors
  {
    lin_scalar_iso,               //< single lin_scalar_iso factor
    lin_scalar_aniso,             //< single lin_scalar_aniso factor
    no_growth_and_lin_scalar_iso  //< composition of no_growth and lin_scalar_iso factors
  };

  constexpr int matid_lin_scalar_iso = MultSplitTestUtils::first_free_matid;
  constexpr int matid_lin_scalar_aniso = MultSplitTestUtils::first_free_matid + 1;
  constexpr int matid_no_growth = MultSplitTestUtils::first_free_matid + 2;
  constexpr int matid_multiplicative_split = MultSplitTestUtils::first_free_matid + 3;

  // input configuration for the inelastic factors modeling growth
  constexpr double growth_factor = 5.27e-7;
  constexpr double reference_concentration = 46456.0;
  const std::vector<double> growth_direction = {1.0, 0.5, 0.3};

  /// material ids of the given inelastic factors, in the order of their composition
  std::vector<int> inelastic_factor_ids(const InelasticFactors inelastic_factors)
  {
    switch (inelastic_factors)
    {
      case InelasticFactors::lin_scalar_iso:
        return {matid_lin_scalar_iso};
      case InelasticFactors::lin_scalar_aniso:
        return {matid_lin_scalar_aniso};
      case InelasticFactors::no_growth_and_lin_scalar_iso:
        return {matid_no_growth, matid_lin_scalar_iso};
      default:
        FOUR_C_THROW("Unexpected inelastic factors.");
    }
  }

  /// set up the multiplicative split with an ElastHyper material with the given elastic summands
  /// and the given inelastic factors
  void set_up_materials(const MultSplitTestUtils::ElasticSummands elastic_summands,
      const InelasticFactors inelastic_factors)
  {
    MultSplitTestUtils::set_up_elast_hyper(elastic_summands);

    // insert all factors into the global problem
    Core::IO::InputParameterContainer lin_scalar_data;
    lin_scalar_data.add("SCALAR1", 1);
    lin_scalar_data.add("SCALAR1_MolarGrowthFac", growth_factor);
    lin_scalar_data.add("SCALAR1_RefConc", reference_concentration);
    Global::Problem::instance()->materials()->insert(matid_lin_scalar_iso,
        Mat::make_parameter(matid_lin_scalar_iso, Core::Materials::MaterialType::mfi_lin_scalar_iso,
            lin_scalar_data));

    // for the lin_scalar_aniso factor, the growth direction is required additionally
    lin_scalar_data.add("GrowthDirection", growth_direction);
    Global::Problem::instance()->materials()->insert(matid_lin_scalar_aniso,
        Mat::make_parameter(matid_lin_scalar_aniso,
            Core::Materials::MaterialType::mfi_lin_scalar_aniso, lin_scalar_data));

    Global::Problem::instance()->materials()->insert(matid_no_growth,
        Mat::make_parameter(matid_no_growth, Core::Materials::MaterialType::mfi_no_growth, {}));

    // only reference the factors of the current test
    const std::vector<int> factor_ids = inelastic_factor_ids(inelastic_factors);
    Core::IO::InputParameterContainer multiplicative_split_data;
    const std::vector<int> elastic_summand_ids =
        MultSplitTestUtils::elastic_summand_ids(elastic_summands);
    multiplicative_split_data.add("NUMMATEL", static_cast<int>(elastic_summand_ids.size()));
    multiplicative_split_data.add("MATIDSEL", elastic_summand_ids);
    multiplicative_split_data.add("NUMFACINEL", static_cast<int>(factor_ids.size()));
    multiplicative_split_data.add("INELDEFGRADFACIDS", factor_ids);
    multiplicative_split_data.add("DENS", 1.0);
    multiplicative_split_data.add("REF_TEMPERATURE", 293.0);
    multiplicative_split_data.add("THERMAL_EXPANSION_COEFFICIENT", 0.0);
    Global::Problem::instance()->materials()->insert(matid_multiplicative_split,
        Mat::make_parameter(matid_multiplicative_split,
            Core::Materials::MaterialType::m_multiplicative_split_defgrad_elasthyper,
            multiplicative_split_data));

    MultSplitTestUtils::set_up_generic_parameters();
  }

  /// Create a parameter list that carries the concentration and the correct differentiation type
  Teuchos::ParameterList create_parameter_list_with_concentration(const double concentration)
  {
    Teuchos::ParameterList params;
    params.set<std::shared_ptr<std::vector<double>>>(
        "scalars", std::make_shared<std::vector<double>>(std::vector<double>{concentration}));
    // only read by evaluate_d_stress_d_scalar()
    params.set<int>("differentiationtype", static_cast<int>(Solid::DifferentiationType::elch));
    return params;
  }

  /*!
   * Checks the multiplicative split material with concentration-driven growth through its public
   * interface only:
   *
   * - the stiffness tensor equals the finite difference of the stress w.r.t. F
   * - the derivative of the stress w.r.t. the concentration equals its finite difference
   */
  class MultiplicativeSplitDefgradElastHyperGrowthTest
      : public ::testing::TestWithParam<
            std::tuple<MultSplitTestUtils::ElasticSummands, InelasticFactors>>
  {
   protected:
    Core::Utils::SingletonOwnerRegistry::ScopeGuard guard_;

    void SetUp() override
    {
      const auto [elastic_summands, inelastic_factors] = GetParam();
      set_up_materials(elastic_summands, inelastic_factors);
      multiplicative_split_ =
          MultSplitTestUtils::create_material<Mat::MultiplicativeSplitDefgradElastHyper>(
              matid_multiplicative_split);
    }

    /// evaluate the stress of the multiplicative split material at F and the given concentration
    [[nodiscard]] Core::LinAlg::SymmetricTensor<double, 3, 3> evaluate_stress(
        const Core::LinAlg::Tensor<double, 3, 3>& F, const double concentration) const
    {
      return MultSplitTestUtils::evaluate(
          *multiplicative_split_, F, create_parameter_list_with_concentration(concentration), {})
          .stress;
    }

    const Core::LinAlg::Tensor<double, 3, 3> F_ = {
        {{1.1, 0.01, 0.03}, {0.04, 1.2, 0.02}, {0.06, 0.05, 1.3}}};
    /// concentration with the material concentration \f$ c \det \boldsymbol{F} = 2 c_\mathrm{ref}
    /// \f$, i.e., the growth \f$ g = k c_\mathrm{ref} > 0 \f$
    const double concentration_ = 2.0 * reference_concentration / Core::LinAlg::det(F_);

    std::shared_ptr<Mat::MultiplicativeSplitDefgradElastHyper> multiplicative_split_;
  };

  TEST_P(MultiplicativeSplitDefgradElastHyperGrowthTest, CmatEqualsFiniteDifference)
  {
    const Teuchos::ParameterList params = create_parameter_list_with_concentration(concentration_);
    const Core::LinAlg::SymmetricTensor<double, 3, 3, 3, 3> cmat =
        MultSplitTestUtils::evaluate(*multiplicative_split_, F_, params, {}).cmat;

    // analytic
    const Core::LinAlg::Tensor<double, 3, 3, 3, 3> dSdF =
        MultSplitTestUtils::compute_stress_deriv_wrt_defgrad(cmat, F_);

    // finite difference
    const auto evaluate_stress_at_fixed_concentration =
        [&](const Core::LinAlg::Tensor<double, 3, 3>& F)
    { return evaluate_stress(F, concentration_); };
    const Core::LinAlg::Tensor<double, 3, 3, 3, 3> dSdF_fd =
        MultSplitTestUtils::finite_difference_stress_deriv_wrt_defgrad(
            evaluate_stress_at_fixed_concentration, F_);

    FOUR_C_EXPECT_NEAR_RELATIVE(dSdF, dSdF_fd, 1.0e-6);
  }

  TEST_P(
      MultiplicativeSplitDefgradElastHyperGrowthTest, DStressDConcentrationEqualsFiniteDifference)
  {
    // analytic
    const Core::LinAlg::SymmetricTensor<double, 3, 3> dSdc =
        multiplicative_split_->evaluate_d_stress_d_scalar(F_,
            MultSplitTestUtils::green_lagrange(F_),
            create_parameter_list_with_concentration(concentration_), {}, MultSplitTestUtils::gp,
            MultSplitTestUtils::ele_gid);

    // finite difference
    const double h = 1.0e-6 * concentration_;
    const Core::LinAlg::SymmetricTensor<double, 3, 3> dSdc_fd =
        (evaluate_stress(F_, concentration_ + h) - evaluate_stress(F_, concentration_ - h)) /
        (2.0 * h);

    FOUR_C_EXPECT_NEAR_RELATIVE(dSdc, dSdc_fd, 1.0e-6);
  }

  INSTANTIATE_TEST_SUITE_P(ElasticSummandsAndInelasticFactors,
      MultiplicativeSplitDefgradElastHyperGrowthTest,
      ::testing::Combine(
          ::testing::ValuesIn(EnumTools::enum_values<MultSplitTestUtils::ElasticSummands>()),
          ::testing::ValuesIn(EnumTools::enum_values<InelasticFactors>())),
      [](const auto& info)
      {
        return std::format("{}_{}", EnumTools::enum_name(std::get<0>(info.param)),
            EnumTools::enum_name(std::get<1>(info.param)));
      });
}  // namespace
