// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <gtest/gtest.h>

#include "4C_mat_elast_coupneohooke.hpp"
#include "4C_material_parameter_base.hpp"
#include "4C_unittest_utils_assertions_test.hpp"
#include "4C_utils_exceptions.hpp"

#include <limits>

namespace
{
  using namespace FourC;

  Mat::Elastic::CoupNeoHooke create_coup_neo_hooke(
      Mat::Elastic::PAR::CoupNeoHooke& material_parameters)
  {
    return {&material_parameters};
  }

  Mat::Elastic::PAR::CoupNeoHooke create_coup_neo_hooke_parameters()
  {
    Core::IO::InputParameterContainer parameters;
    parameters.add("YOUNG", 2.0);
    parameters.add("NUE", 0.3);

    return Mat::Elastic::PAR::CoupNeoHooke(
        Core::Mat::PAR::Parameter::Data{.parameters = parameters});
  }

  void evaluate_principal_derivatives(Mat::Elastic::CoupNeoHooke& summand, const double i3)
  {
    Core::LinAlg::Matrix<3, 1> prinv(Core::LinAlg::Initialization::zero);
    prinv(0) = 3.0;
    prinv(1) = 3.0;
    prinv(2) = i3;

    Core::LinAlg::Matrix<3, 1> dPI(Core::LinAlg::Initialization::zero);
    Core::LinAlg::Matrix<6, 1> ddPII(Core::LinAlg::Initialization::zero);

    summand.add_derivatives_principal(dPI, ddPII, prinv, 0, 0);
  }

  TEST(CoupNeoHookeTest, PrincipalDerivativesThrowForNonPositiveI3)
  {
    auto material_parameters = create_coup_neo_hooke_parameters();
    auto summand = create_coup_neo_hooke(material_parameters);

    FOUR_C_EXPECT_THROW_WITH_MESSAGE(evaluate_principal_derivatives(summand, 0.0), Core::Exception,
        "Error in principal derivative computation. Derivative of strain energy wrt. I3 = -inf");
  }

  TEST(CoupNeoHookeTest, PrincipalDerivativesThrowForOverflowingI3)
  {
    auto material_parameters = create_coup_neo_hooke_parameters();
    auto summand = create_coup_neo_hooke(material_parameters);

    FOUR_C_EXPECT_THROW_WITH_MESSAGE(
        evaluate_principal_derivatives(summand, std::numeric_limits<double>::min()),
        Core::Exception,
        "Error in principal derivative computation. Derivative of strain energy wrt. I3 = -inf");
  }
}  // namespace
