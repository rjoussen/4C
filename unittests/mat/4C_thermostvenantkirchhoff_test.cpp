// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <gtest/gtest.h>

#include "4C_linalg_symmetric_tensor.hpp"
#include "4C_linalg_tensor.hpp"
#include "4C_linalg_tensor_generators.hpp"
#include "4C_mat_thermostvenantkirchhoff.hpp"
#include "4C_material_parameter_base.hpp"
#include "4C_unittest_utils_assertions_test.hpp"

#include <Teuchos_ParameterList.hpp>

#include <vector>

namespace
{
  using namespace FourC;

  class ThermoStVenantKirchhoffTest : public ::testing::Test
  {
   protected:
    void SetUp() override
    {
      Core::IO::InputParameterContainer container;
      // temperature-dependent Young's modulus E(T) = sum_i a_i T^i
      container.add("YOUNG", std::vector<double>{2.1e5, -3.0e2, 1.5, -4.0e-3});
      container.add("NUE", 0.3);
      container.add("DENS", 1.0);
      container.add("THEXPANS", 1.2e-5);
      container.add("INITTEMP", initial_temperature_);

      parameters_ = std::make_shared<Mat::PAR::ThermoStVenantKirchhoff>(
          Core::Mat::PAR::Parameter::Data{.parameters = container});
      material_ = std::make_shared<Mat::ThermoStVenantKirchhoff>(parameters_.get());
    }

    [[nodiscard]] Mat::HeatSource evaluate_heat_source(
        const double temperature, const Core::LinAlg::SymmetricTensor<double, 3, 3>& strain)
    {
      return material_->evaluate_mechanical_heat_source(temperature,
          Mat::KinematicState::from_linear_strain(strain, strain_rate_), context_, 0, 0);
    }

    const double initial_temperature_ = 290.0;
    const double temperature_ = 330.0;
    const double time_step_size_ = 0.1;
    const double total_time_ = 1.0;
    const Mat::EvaluationContext<3> context_{.total_time = &total_time_,
        .time_step_size = &time_step_size_,
        .xi = nullptr,
        .ref_coords = nullptr};

    const Core::LinAlg::SymmetricTensor<double, 3, 3> strain_ =
        Core::LinAlg::assume_symmetry(Core::LinAlg::Tensor<double, 3, 3>{
            {{1.0e-3, 4.0e-4, -2.0e-4}, {4.0e-4, -5.0e-4, 3.0e-4}, {-2.0e-4, 3.0e-4, 2.0e-3}}});
    const Core::LinAlg::SymmetricTensor<double, 3, 3> strain_rate_ =
        Core::LinAlg::assume_symmetry(Core::LinAlg::Tensor<double, 3, 3>{
            {{2.0e-2, -1.0e-2, 5.0e-3}, {-1.0e-2, 3.0e-2, 7.0e-3}, {5.0e-3, 7.0e-3, -4.0e-2}}});

    std::shared_ptr<Mat::PAR::ThermoStVenantKirchhoff> parameters_;
    std::shared_ptr<Mat::ThermoStVenantKirchhoff> material_;
  };

  TEST_F(ThermoStVenantKirchhoffTest, HeatingUsesTemperatureDerivativeOfStress)
  {
    // the thermoelastic heating is T . partial S/partial T : dE/dt. Since the material has no
    // internal variables, partial S/partial T equals the total derivative used by the structure.
    Teuchos::ParameterList params;
    params.set<double>("temperature", temperature_);
    const auto dS_dT = material_->evaluate_d_stress_d_scalar(
        Core::LinAlg::get_full(Core::LinAlg::TensorGenerators::identity<double, 3, 3>), strain_,
        params, context_, 0, 0);

    const Mat::HeatSource heat_source = evaluate_heat_source(temperature_, strain_);

    FOUR_C_EXPECT_NEAR(heat_source.derivative_wrt_strain_rate, temperature_ * dS_dT, 1.0e-8);
    EXPECT_NEAR(heat_source.value, temperature_ * Core::LinAlg::ddot(dS_dT, strain_rate_), 1.0e-8);
  }

  TEST_F(ThermoStVenantKirchhoffTest, HeatSourceDerivativesMatchFiniteDifferences)
  {
    const Mat::HeatSource heat_source = evaluate_heat_source(temperature_, strain_);

    // derivative w.r.t. temperature
    const double dT = 1.0e-4;
    const double fd_derivative_wrt_temperature =
        (evaluate_heat_source(temperature_ + dT, strain_).value -
            evaluate_heat_source(temperature_ - dT, strain_).value) /
        (2.0 * dT);
    EXPECT_NEAR(heat_source.derivative_wrt_temperature, fd_derivative_wrt_temperature,
        1.0e-6 * std::abs(fd_derivative_wrt_temperature));

    // derivative w.r.t. strain (a perturbation of the off-diagonal entry of a symmetric tensor
    // perturbs both E_ij and E_ji)
    const double dE = 1.0e-7;
    for (unsigned i = 0; i < 3; ++i)
    {
      for (unsigned j = i; j < 3; ++j)
      {
        auto strain_plus = strain_;
        auto strain_minus = strain_;
        strain_plus(i, j) += dE;
        strain_minus(i, j) -= dE;
        const double scale = i == j ? 1.0 : 2.0;
        const double fd_derivative = (evaluate_heat_source(temperature_, strain_plus).value -
                                         evaluate_heat_source(temperature_, strain_minus).value) /
                                     (2.0 * dE * scale);
        EXPECT_NEAR(heat_source.derivative_wrt_strain(i, j), fd_derivative,
            1.0e-6 * std::abs(fd_derivative) + 1.0e-6);
      }
    }
  }
}  // namespace
