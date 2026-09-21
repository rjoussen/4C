// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <gtest/gtest.h>

#include "4C_contact_tsi_interface.hpp"

#include "4C_contact_friction_node.hpp"
#include "4C_fem_discretization.hpp"
#include "4C_global_data.hpp"
#include "4C_io_control.hpp"
#include "4C_linalg_sparsematrix.hpp"

#include <array>
#include <cmath>

namespace
{
  using namespace FourC;

  /**
   * @brief Test individual TSI contact Jacobian contributions without contact search or
   * integration.
   *
   * Continuum notation follows Seitz (2019), Computational Methods for Thermo-Elasto-Plastic
   * Contact, Section 2.7. Source and target correspond to sides (1) and (2), respectively.
   * The code uses lambda = -t_c, hence lambda_n = -p_n is positive in compression.
   *
   * One source node (displacement DOFs 0, 1, 2) and one target node (DOFs 3, 4, 5) provide
   * the rows needed for assembly. We prescribe their mortar data directly:
   * \f[
   *   D=2,\quad M=3,\quad \Delta t=0.5,\quad n=(0,0,1)^T,\quad
   *   \lambda=(-2,-3,5)^T,\quad j=(0.4,0.2,0)^T .
   * \f]
   * Here D is the diagonal dual-mortar source weight, M a target coupling entry, lambda the
   * mechanical multiplier, and j the mortar-weighted relative displacement increment.
   * These are synthetic, independent assembly inputs, not an integrated contact patch or a
   * converged Coulomb state. In particular, D and M need not satisfy partition of unity here.
   *
   * With the tangential projector P, the signed dissipation term used by the implementation is
   * \f[
   *   P=I-n\otimes n,\qquad w=-\frac{\lambda^T Pj}{\Delta t D}.
   * \f]
   * The discrete frictional power H=-w approximates t_tau . v_tau in Eq. (2.129). It is
   * not the total contact dissipation D_c of Eq. (2.120), which includes heat transfer.
   * Our off-solution fixture deliberately has w=1.4 (and hence nonphysical negative H).
   * Using D different from one exposes a missing normalization. Two nonzero tangential
   * components exercise both sliding directions.
   *
   * Each test perturbs only the dependency whose Jacobian contribution is being tested.
   * Derivatives are checked with the central difference
   * \f[
   *   f'(x)\approx\frac{f(x+\varepsilon)-f(x-\varepsilon)}{2\varepsilon},
   *   \qquad \varepsilon=10^{-6}.
   * \f]
   * Thermal quantities temporarily use the first displacement DOF of each node (0 and 3).
   * The coupled strategy later maps these entries onto the actual thermal DOFs.
   */
  class TSIInterfaceTest : public ::testing::Test
  {
   protected:
    void SetUp() override
    {
      Global::Problem::instance()->set_spatial_approximation_type(
          Core::FE::ShapeFunctionType::polynomial);
      Global::Problem::instance()->set_output_control_file(
          std::make_shared<Core::IO::OutputControl>(MPI_COMM_WORLD, "tsi",
              Core::FE::ShapeFunctionType::polynomial, "", "", "tsi_interface_test", 3, 0, 1, false,
              false));

      Teuchos::ParameterList params;
      params.sublist("PARALLEL REDISTRIBUTION")
          .set("GHOSTING_STRATEGY", Mortar::ExtendGhosting::redundant_target);
      params.set("SEARCH_ALGORITHM", Mortar::search_bfele);
      params.set("SEARCH_PARAM", 1.0);
      params.set("SEARCH_USE_AUX_POS", false);
      params.set("NURBS", false);
      params.set("LM_SHAPEFCN", Mortar::shape_dual);
      params.set("NONSMOOTH_GEOMETRIES", false);
      params.set("Two_half_pass", false);
      params.set("CONSTRAINT_DIRECTIONS", CONTACT::ConstraintDirection::ntt);
      params.set("FRICTION", CONTACT::FrictionType::coulomb);
      params.set("PROBTYPE", CONTACT::Problemtype::tsi);
      params.set("TIMESTEP", dt);
      params.set("FRCOEFF", 0.3);
      params.set("FRBOUND", 0.0);
      params.set("GP_SLIP_INCR", false);
      params.set("FRLESS_FIRST", false);
      params.set("SEMI_SMOOTH_CN", 1.0);
      params.set("SEMI_SMOOTH_CT", 1.0);

      data = std::make_shared<CONTACT::InterfaceDataContainer>();
      interface =
          std::make_shared<CONTACT::TSIInterface>(data, 0, MPI_COMM_WORLD, 3, params, false);
      const std::array<double, 3> coordinates{0., 0., 0.};
      source = std::make_shared<CONTACT::FriNode>(
          0, coordinates, 0, std::vector<int>{0, 1, 2}, true, true, true);
      auto target = std::make_shared<CONTACT::FriNode>(
          1, coordinates, 0, std::vector<int>{3, 4, 5}, false, false, true);
      interface->add_mortar_node(source);
      interface->add_mortar_node(target);
      // Attach the nodes to the discretization before initializing their contact data.
      // No elements are needed: these tests supply already integrated mortar entries.
      interface->discret().fill_complete(Core::FE::OptionsFillComplete::none());
      source->initialize_data_container();
      source->initialize_tsi_data_container(0., 100.);
      source->mo_data().get_d().resize(1);
      source->mo_data().get_d()[0] = 2.;
      source->mo_data().get_m()[1] = 3.;
      source->mo_data().n()[2] = 1.;
      source->mo_data().lm()[0] = -2.;
      source->mo_data().lm()[1] = -3.;
      source->mo_data().lm()[2] = 5.;
      source->fri_data().jump()[0] = 0.4;
      source->fri_data().jump()[1] = 0.2;
      source->data().txi()[0] = 1.;
      source->data().teta()[1] = 1.;
      source->data().get_deriv_n().resize(3);
      source->data().get_deriv_txi().resize(3);
      source->data().get_deriv_teta().resize(3);
      source->fri_data().get_deriv_jump().resize(3);

      const int node_id = 0;
      data->active_nodes() = std::make_shared<Core::LinAlg::Map>(1, 1, &node_id, 0, MPI_COMM_WORLD);
      data->slip_nodes() = data->active_nodes();
      // In local n/t/t constraint coordinates, rows 1 and 2 contain the two slip equations.
      // They are constraint row labels, not the Cartesian directions of the nodal normal.
      const std::array<int, 2> tangent_dofs{1, 2};
      data->slip_t() =
          std::make_shared<Core::LinAlg::Map>(2, 2, tangent_dofs.data(), 0, MPI_COMM_WORLD);
      dofs = std::make_shared<Core::LinAlg::Map>(6, 0, MPI_COMM_WORLD);
      data->cn_values() = std::make_shared<Core::LinAlg::Vector<double>>(*data->active_nodes());
      data->ct_values() = std::make_shared<Core::LinAlg::Vector<double>>(*data->active_nodes());
      // Deliberately differ from the input values, as at contact edges and corners.
      data->cn_values()->put_scalar(2.);
      data->ct_values()->put_scalar(3.);
    }

    void TearDown() override { Global::Problem::instance()->set_output_control_file(nullptr); }

    /// Evaluate the signed term w directly, independently of its analytical Jacobian assembly.
    double dissipation() const
    {
      const Core::LinAlg::Matrix<3, 1> lm(source->mo_data().lm(), true);
      const Core::LinAlg::Matrix<3, 1> n(source->mo_data().n(), true);
      const Core::LinAlg::Matrix<3, 1> jump(source->fri_data().jump(), true);
      return -(lm.dot(jump) - lm.dot(n) * jump.dot(n)) / (dt * source->mo_data().get_d()[0]);
    }

    /// Extract column c as A e_c, including completion of the finite-element matrix assembly.
    Core::LinAlg::Vector<double> column(Core::LinAlg::SparseMatrix& matrix, int col) const
    {
      matrix.complete(*dofs, *dofs);
      Core::LinAlg::Vector<double> basis(*dofs, true), result(*dofs, true);
      basis.get_values()[col] = 1.;
      matrix.multiply(false, basis, result);
      return result;
    }

    static constexpr double dt = 0.5;
    static constexpr double eps = 1.e-6;
    std::shared_ptr<CONTACT::InterfaceDataContainer> data;
    std::shared_ptr<CONTACT::TSIInterface> interface;
    std::shared_ptr<CONTACT::FriNode> source;
    std::shared_ptr<Core::LinAlg::Map> dofs;
  };

  /**
   * @brief Check the derivative of the inverse source weight inside the dissipation term.
   *
   * Prescribe D(x)=2+0.7x, holding lambda, n and j fixed. Then
   * \f[
   *   w_{,x}=\frac{\lambda^T Pj}{\Delta t D^2}D_{,x}
   *         =-\frac{w}{D}D_{,x}.
   * \f]
   * assemble_dm_lin_diss supplies D w_{,x} and M w_{,x}, with the outer weights held fixed
   * at their base values. Their derivatives belong to the next test. The normalization
   * derivative is scalar: adding it inside a loop over three spatial components triples it.
   */
  TEST_F(TSIInterfaceTest, DissipationNormalizationDerivative)
  {
    source->data().get_deriv_d()[0][0] = 0.7;
    Core::LinAlg::SparseMatrix dlin(*dofs, 6, true, false, Core::LinAlg::SparseMatrix::FE_MATRIX);
    Core::LinAlg::SparseMatrix mlin(*dofs, 6, true, false, Core::LinAlg::SparseMatrix::FE_MATRIX);
    interface->assemble_dm_lin_diss(&dlin, &mlin, nullptr, nullptr, 1.);

    source->mo_data().get_d()[0] = 2. + eps * 0.7;
    const double plus = dissipation();
    source->mo_data().get_d()[0] = 2. - eps * 0.7;
    const double minus = dissipation();
    const double derivative = (plus - minus) / (2. * eps);
    EXPECT_NEAR(column(dlin, 0).local_values_as_span()[0], 2. * derivative, 1.e-9);
    EXPECT_NEAR(column(mlin, 0).local_values_as_span()[3], 3. * derivative, 1.e-9);
  }

  /**
   * @brief Check the outer mortar-weight derivatives in the product rule.
   *
   * This test freezes w at the base state and varies D(x)=2+0.7x and M(x)=3-0.4x:
   * \f[
   *   \left.\frac{d(Dw)}{dx}\right|_{w\ \mathrm{fixed}}=D_{,x}w,\qquad
   *   \left.\frac{d(-Mw)}{dx}\right|_{w\ \mathrm{fixed}}=-M_{,x}w .
   * \f]
   * The minus sign is the convention of assemble_lin_dm_x for target-side M terms.
   * Its caller selects the sign needed by the thermal residual. The frozen w must still
   * contain the base-state normalization 1/D; accidentally using D=1 fails this test.
   */
  TEST_F(TSIInterfaceTest, MortarDerivativeUsesNormalizedDissipation)
  {
    source->data().get_deriv_d()[0][0] = 0.7;
    source->data().get_deriv_m()[1][0] = -0.4;
    Core::LinAlg::SparseMatrix dlin(*dofs, 6, true, false, Core::LinAlg::SparseMatrix::FE_MATRIX);
    Core::LinAlg::SparseMatrix mlin(*dofs, 6, true, false, Core::LinAlg::SparseMatrix::FE_MATRIX);
    interface->assemble_lin_dm_x(
        &dlin, &mlin, 1., CONTACT::TSIInterface::LinDM_Diss, data->active_nodes());

    const double w = dissipation();
    const double derivative_d = ((2. + eps * 0.7) * w - (2. - eps * 0.7) * w) / (2. * eps);
    const double derivative_m = -((3. - eps * 0.4) * w - (3. + eps * 0.4) * w) / (2. * eps);
    EXPECT_NEAR(column(dlin, 0).local_values_as_span()[0], derivative_d, 1.e-9);
    EXPECT_NEAR(column(mlin, 0).local_values_as_span()[3], derivative_m, 1.e-9);
  }

  /**
   * @brief Differentiate the mechanical slip residual with respect to source temperature.
   *
   * With tangential vectors expressed in the two local tangent directions, define
   * \f[
   *   a=\lambda_\tau+c_t j_\tau,\qquad
   *   C_\tau=\|a\|\lambda_\tau-\mu(\vartheta_c)(\lambda_n-c_n g)a,\qquad
   *   \frac{\partial C_\tau}{\partial T^{(1)}}
   *     =-\mu'(\vartheta_c)(\lambda_n-c_n g)a .
   * \f]
   * Here lambda_n=lambda.n=-p_n, g is the weighted gap, and c_n/c_t are the nodal semismooth
   * parameters. T^(1)=20 is strictly above projected T^(2)=10, so vartheta_c=T^(1).
   * The maximum-temperature law is smooth here and only the source-temperature derivative
   * is involved. With mu_0=0.3, T_0=0 and T_d=100, Eq. (2.125) becomes
   * \f[
   *   \mu(\vartheta_c)=0.3\left(\frac{\vartheta_c-100}{100}\right)^2 .
   * \f]
   * The nodal values c_n=2 and c_t=3 deliberately differ from the input values (both one),
   * as can happen when parameters are scaled at contact edges/corners. Compare the TSI
   * tangent against finite differences of the base mechanical residual, not a duplicate
   * implementation of the tangent. Both tangential constraint rows must agree.
   */
  TEST_F(TSIInterfaceTest, SlipTemperatureDerivativeUsesNodalParameters)
  {
    source->tsi_data().temp() = 20.;
    source->tsi_data().temp_target() = 10.;
    source->data().getg() = -0.1;
    Core::LinAlg::SparseMatrix lm(*dofs, 6), disp(*dofs, 6), temp(*dofs, 6);
    Core::LinAlg::Vector<double> rhs(*dofs, true);
    interface->assemble_lin_slip(lm, disp, temp, rhs);
    auto tangent = column(temp, 0);

    auto residual = [&](double temperature)
    {
      source->tsi_data().temp() = temperature;
      Core::LinAlg::SparseMatrix lm_dummy(*dofs, 6), disp_dummy(*dofs, 6);
      Core::LinAlg::Vector<double> result(*dofs, true);
      interface->CONTACT::Interface::assemble_lin_slip(lm_dummy, disp_dummy, result);
      result.scale(-1.);  // The assembly returns the negative residual.
      return result;
    };
    auto plus = residual(20. + eps);
    auto minus = residual(20. - eps);
    for (int row : {1, 2})
      EXPECT_NEAR(tangent.local_values_as_span()[row],
          (plus.local_values_as_span()[row] - minus.local_values_as_span()[row]) / (2. * eps),
          1.e-8);
  }

  /**
   * @brief Check both slip equations while the contact frame and projected temperature vary.
   *
   * A scalar displacement x rotates n and txi about the y axis and changes the weighted
   * gap and jump. Prescribe the hotter projected target temperature as
   * \f$\vartheta_c=40+0.8x\f$, including its derivative in deriv_temp_target_disp.
   * This exercises the displacement derivative of mu as well as the mechanical slip tangent;
   * it does not test how the interpolator obtains the projected temperature.
   * Perturb all three multipliers independently.
   * The active/slip set stays fixed: this checks a smooth branch of the Newton equations.
   */
  TEST_F(TSIInterfaceTest, SlipJacobianWithRotatingFrameAndHotterTarget)
  {
    auto set_displacement = [&](double x)
    {
      source->mo_data().n()[0] = std::sin(x);
      source->mo_data().n()[2] = std::cos(x);
      source->data().txi()[0] = std::cos(x);
      source->data().txi()[2] = -std::sin(x);
      source->fri_data().jump()[0] = 0.4 + 0.7 * x;
      source->fri_data().jump()[1] = 0.2 - 0.3 * x;
      source->fri_data().jump()[2] = 0.1 + 0.2 * x;
      source->data().getg() = -0.1 + 0.6 * x;
      source->tsi_data().temp_target() = 40. + 0.8 * x;
    };
    source->tsi_data().temp() = 20.;
    set_displacement(0.);
    for (int d = 0; d < 3; ++d)
    {
      source->data().get_deriv_n()[d].resize(1);
      source->data().get_deriv_txi()[d].resize(1);
      source->data().get_deriv_teta()[d].resize(1);
      source->data().get_deriv_n()[d][0] = d == 0 ? 1. : 0.;
      source->data().get_deriv_txi()[d][0] = d == 2 ? -1. : 0.;
      source->data().get_deriv_teta()[d][0] = 0.;
      source->fri_data().get_deriv_jump()[d][0] = std::array{0.7, -0.3, 0.2}[d];
    }
    source->data().get_deriv_g()[0] = 0.6;
    source->tsi_data().deriv_temp_target_disp()[0] = 0.8;
    source->tsi_data().deriv_temp_target_temp()[3] = 1.;

    Core::LinAlg::SparseMatrix lm(*dofs, 6), disp(*dofs, 6), temp(*dofs, 6);
    Core::LinAlg::Vector<double> rhs(*dofs, true);
    interface->assemble_lin_slip(lm, disp, temp, rhs);
    auto residual = [&]()
    {
      Core::LinAlg::SparseMatrix lm_dummy(*dofs, 6), disp_dummy(*dofs, 6);
      Core::LinAlg::Vector<double> result(*dofs, true);
      interface->CONTACT::Interface::assemble_lin_slip(lm_dummy, disp_dummy, result);
      result.scale(-1.);
      return result;
    };
    auto check = [&](Core::LinAlg::SparseMatrix& matrix, int col, auto perturb)
    {
      const auto tangent = column(matrix, col);
      perturb(eps);
      const auto plus = residual();
      perturb(-eps);
      const auto minus = residual();
      perturb(0.);
      for (int row : {1, 2})
        EXPECT_NEAR(tangent.local_values_as_span()[row],
            (plus.local_values_as_span()[row] - minus.local_values_as_span()[row]) / (2. * eps),
            1.e-8)
            << "row " << row << ", column " << col;
    };
    check(disp, 0, set_displacement);
    check(temp, 3, [&](double t) { source->tsi_data().temp_target() = 40. + t; });
    for (int d = 0; d < 3; ++d)
    {
      const double base = source->mo_data().lm()[d];
      check(lm, d, [&](double increment) { source->mo_data().lm()[d] = base + increment; });
    }
  }

  /**
   * @brief Check all four Jacobian blocks of the thermal contact equation.
   *
   * With q=thermo_lm, w=dissipation(), beta_c=2/3 and delta_c=1/3, assemble_lin_conduct
   * differentiates
   * \f$C_T=Dq-\delta_c Dw-\beta_c(\lambda\cdot n)(DT^{(1)}-MT^{(2)})\f$.
   * Perturb the geometry, both temperatures, q and every component of lambda separately.
   * An inclined, rotating normal and a jump with a normal component expose errors that
   * a flat sliding patch cannot detect. The geometry perturbation also changes D and M.
   */
  TEST_F(TSIInterfaceTest, ThermalContactJacobian)
  {
    data->i_mortar().set("HEATTRANSSLAVE", 1.).set("HEATTRANSMASTER", 2.);
    data->i_mortar().set("LM_QUAD", Mortar::lagmult_undefined);
    auto* target = dynamic_cast<CONTACT::Node*>(interface->discret().g_node(1));
    target->initialize_tsi_data_container(0., 100.);
    source->tsi_data().temp() = 20.;
    target->tsi_data().temp() = 40.;
    source->tsi_data().thermo_lm() = 2.;
    auto set_displacement = [&](double x)
    {
      source->mo_data().n()[0] = std::sin(0.3 + x);
      source->mo_data().n()[2] = std::cos(0.3 + x);
      source->mo_data().get_d()[0] = 2. + 0.7 * x;
      source->mo_data().get_m()[1] = 3. - 0.4 * x;
      source->fri_data().jump()[0] = 0.4 + 0.7 * x;
      source->fri_data().jump()[1] = 0.2 - 0.3 * x;
      source->fri_data().jump()[2] = 0.1 + 0.2 * x;
    };
    set_displacement(0.);
    source->data().get_deriv_d()[0][0] = 0.7;
    source->data().get_deriv_m()[1][0] = -0.4;
    for (int d = 0; d < 3; ++d)
    {
      source->data().get_deriv_n()[d].resize(1);
      source->data().get_deriv_n()[d][0] = std::array{std::cos(0.3), 0., -std::sin(0.3)}[d];
      source->fri_data().get_deriv_jump()[d][0] = std::array{0.7, -0.3, 0.2}[d];
    }
    Core::LinAlg::SparseMatrix disp(*dofs, 6, true, false, Core::LinAlg::SparseMatrix::FE_MATRIX);
    Core::LinAlg::SparseMatrix temp(*dofs, 6, true, false, Core::LinAlg::SparseMatrix::FE_MATRIX);
    Core::LinAlg::SparseMatrix q(*dofs, 6, true, false, Core::LinAlg::SparseMatrix::FE_MATRIX);
    Core::LinAlg::SparseMatrix lm(*dofs, 6, true, false, Core::LinAlg::SparseMatrix::FE_MATRIX);
    interface->assemble_lin_conduct(disp, temp, q, lm);
    auto residual = [&]()
    {
      const double d = source->mo_data().get_d()[0];
      const double m = source->mo_data().get_m()[1];
      const Core::LinAlg::Matrix<3, 1> normal(source->mo_data().n(), true);
      const Core::LinAlg::Matrix<3, 1> multiplier(source->mo_data().lm(), true);
      return d * source->tsi_data().thermo_lm() - d * dissipation() / 3. -
             2. / 3. * multiplier.dot(normal) *
                 (d * source->tsi_data().temp() - m * target->tsi_data().temp());
    };
    auto check = [&](Core::LinAlg::SparseMatrix& matrix, int col, auto perturb)
    {
      const auto tangent = column(matrix, col);
      perturb(eps);
      const double plus = residual();
      perturb(-eps);
      const double minus = residual();
      perturb(0.);
      EXPECT_NEAR(tangent.local_values_as_span()[0], (plus - minus) / (2. * eps), 1.e-7)
          << "column " << col;
    };
    check(disp, 0, set_displacement);
    check(temp, 0, [&](double t) { source->tsi_data().temp() = 20. + t; });
    check(temp, 3, [&](double t) { target->tsi_data().temp() = 40. + t; });
    check(q, 0, [&](double increment) { source->tsi_data().thermo_lm() = 2. + increment; });
    for (int d = 0; d < 3; ++d)
    {
      const double base = source->mo_data().lm()[d];
      check(lm, d, [&](double increment) { source->mo_data().lm()[d] = base + increment; });
    }
  }

}  // namespace
