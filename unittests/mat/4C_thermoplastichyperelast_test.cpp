// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <gtest/gtest.h>

#include "4C_linalg_symmetric_tensor.hpp"
#include "4C_linalg_tensor.hpp"
#include "4C_mat_thermoplastichyperelast.hpp"
#include "4C_material_parameter_base.hpp"

#include <optional>

namespace
{
  using namespace FourC;

  class ThermoPlasticHyperElastTest : public ::testing::Test
  {
   protected:
    void SetUp() override
    {
      Core::IO::InputParameterContainer container;
      container.add("YOUNG", 2.1e5);
      container.add("NUE", 0.3);
      container.add("DENS", 1.0);
      container.add("CTE", 1.2e-5);
      container.add("INITTEMP", initial_temperature_);
      container.add("YIELD", 300.0);
      container.add("ISOHARD", 100.0);
      container.add("SATHARDENING", 400.0);
      container.add("HARDEXPO", 10.0);
      container.add("YIELDSOFT", 0.002);
      container.add("HARDSOFT", 0.002);
      container.add("TOL", 1.0e-8);

      parameters_ = std::make_shared<Mat::PAR::ThermoPlasticHyperElast>(
          Core::Mat::PAR::Parameter::Data{.parameters = container});
      material_ = std::make_shared<Mat::ThermoPlasticHyperElast>(parameters_.get());
      // no plastic history, i.e., the heat source only consists of the thermoelastic heating
      material_->setup(1, {}, {});
    }

    [[nodiscard]] Mat::HeatSource evaluate_heat_source(
        const double temperature, const Core::LinAlg::SymmetricTensor<double, 3, 3>& strain)
    {
      return material_->evaluate_mechanical_heat_source(temperature,
          Mat::KinematicState{
              .strain = strain, .strain_rate = strain_rate_, .defgrad = std::nullopt},
          context_, 0, 0);
    }

    const double initial_temperature_ = 290.0;
    const double temperature_ = 330.0;
    const double time_step_size_ = 0.1;
    const double total_time_ = 1.0;
    const Mat::EvaluationContext<3> context_{.total_time = &total_time_,
        .time_step_size = &time_step_size_,
        .xi = nullptr,
        .ref_coords = nullptr};

    // finite Green-Lagrange strain, i.e., J differs noticeably from 1
    const Core::LinAlg::SymmetricTensor<double, 3, 3> strain_ =
        Core::LinAlg::assume_symmetry(Core::LinAlg::Tensor<double, 3, 3>{
            {{5.0e-2, 2.0e-2, -1.0e-2}, {2.0e-2, -3.0e-2, 1.5e-2}, {-1.0e-2, 1.5e-2, 8.0e-2}}});
    const Core::LinAlg::SymmetricTensor<double, 3, 3> strain_rate_ =
        Core::LinAlg::assume_symmetry(Core::LinAlg::Tensor<double, 3, 3>{
            {{2.0e-2, -1.0e-2, 5.0e-3}, {-1.0e-2, 3.0e-2, 7.0e-3}, {5.0e-3, 7.0e-3, -4.0e-2}}});

    std::shared_ptr<Mat::PAR::ThermoPlasticHyperElast> parameters_;
    std::shared_ptr<Mat::ThermoPlasticHyperElast> material_;
  };

  TEST_F(ThermoPlasticHyperElastTest, HeatSourceDerivativesMatchFiniteDifferences)
  {
    const Mat::HeatSource heat_source = evaluate_heat_source(temperature_, strain_);
    ASSERT_NE(heat_source.value, 0.0);

    // derivative w.r.t. temperature
    const double dT = 1.0e-4;
    const double fd_derivative_wrt_temperature =
        (evaluate_heat_source(temperature_ + dT, strain_).value -
            evaluate_heat_source(temperature_ - dT, strain_).value) /
        (2.0 * dT);
    EXPECT_NEAR(heat_source.derivative_wrt_temperature, fd_derivative_wrt_temperature,
        1.0e-6 * std::abs(fd_derivative_wrt_temperature));

    // derivative w.r.t. strain, i.e., the strain derivative of the stress-temperature modulus (a
    // perturbation of the off-diagonal entry of a symmetric tensor perturbs both E_ij and E_ji)
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
