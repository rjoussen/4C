// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <gtest/gtest.h>

#include "4C_contact_tsi_interface.hpp"

#include "4C_contact_abstract_data_container.hpp"
#include "4C_contact_friction_node.hpp"
#include "4C_contact_lagrange_strategy_tsi.hpp"
#include "4C_fem_discretization.hpp"
#include "4C_global_data.hpp"
#include "4C_io_control.hpp"
#include "4C_linalg_blocksparsematrix.hpp"
#include "4C_linalg_mapextractor.hpp"
#include "4C_linalg_serialdensesolver.hpp"
#include "4C_linalg_sparsematrix.hpp"
#include "4C_linalg_utils_sparse_algebra_manipulation.hpp"

#include <array>
#include <cmath>
#include <optional>

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
    void check_coupled_strategy(
        bool finite_difference, bool active = true, bool equilibrium = false);

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

  TEST_F(TSIInterfaceTest, RecreatedStrategyStartsWithZeroContactNorms)
  {
    auto params = data->i_mortar();
    params.set("SYSTEM", CONTACT::SystemType::condensed);
    params.set("STRATEGY", CONTACT::SolvingStrategy::lagmult);
    params.sublist("PARALLEL REDISTRIBUTION")
        .set("PARALLEL_REDIST", Mortar::ParallelRedist::redist_none);
    Core::LinAlg::Map nodes(2, 0, MPI_COMM_WORLD);
    std::optional<CONTACT::LagrangeStrategyTsi> strategy;
    for (int construction = 0; construction < 2; ++construction)
    {
      strategy.emplace(std::make_shared<CONTACT::AbstractStrategyDataContainer>(), dofs.get(),
          &nodes, params, std::vector<std::shared_ptr<CONTACT::Interface>>{interface}, 3,
          MPI_COMM_WORLD, 0., 7);
      EXPECT_EQ(strategy->mech_contact_res_, 0.);
      EXPECT_EQ(strategy->mech_contact_incr_, 0.);
      EXPECT_EQ(strategy->thermo_contact_incr_, 0.);
      // Reuse storage containing nonzero norms from a previous strategy.
      strategy->mech_contact_res_ = 1.;
      strategy->mech_contact_incr_ = 2.;
      strategy->thermo_contact_incr_ = 3.;
      strategy.reset();
    }
  }

  using ContactVectors = CONTACT::TsiContactLinearization::Vectors;

  ContactVectors apply_linearization(
      const CONTACT::TsiContactLinearization& system, const ContactVectors& increment)
  {
    ContactVectors result;
    for (int row = 0; row < 4; ++row)
    {
      result[row] =
          std::make_shared<Core::LinAlg::Vector<double>>(system.residual[row]->get_map(), true);
      for (int col = 0; col < 4; ++col)
      {
        if (!system.matrix[row][col]) continue;
        const auto& block = *system.matrix[row][col];
        Core::LinAlg::Vector<double> x(block.domain_map(), true), y(block.range_map(), true);
        Core::LinAlg::export_to(*increment[col], x);
        block.multiply(false, x, y);
        CONTACT::Utils::add_vector(y, *result[row]);
      }
    }
    return result;
  }

  ContactVectors solve_linearization(const CONTACT::TsiContactLinearization& system,
      const std::array<std::shared_ptr<const Core::LinAlg::Map>, 2>& dbc)
  {
    const auto comm = system.residual[0]->get_map().get_comm();
    std::array<std::vector<int>, 4> gids;
    std::array<int, 5> offset{};
    ContactVectors basis, result;
    for (int block = 0; block < 4; ++block)
    {
      const auto& map = system.residual[block]->get_map();
      std::vector<int> local;
      for (int i = 0; i < map.num_my_elements(); ++i) local.push_back(map.gid(i));
      for (const auto& rank : Core::Communication::all_gather(local, comm))
        gids[block].insert(gids[block].end(), rank.begin(), rank.end());
      offset[block + 1] = offset[block] + gids[block].size();
      basis[block] = std::make_shared<Core::LinAlg::Vector<double>>(map, true);
      result[block] = std::make_shared<Core::LinAlg::Vector<double>>(map, true);
    }
    const int size = offset[4];
    Core::LinAlg::SerialDenseMatrix matrix(size, size, true);
    Core::LinAlg::SerialDenseVector rhs(size), solution(size);
    std::vector<int> constrained(size, 0);
    for (int row = 0; row < 4; ++row)
      for (int i = 0; i < static_cast<int>(gids[row].size()); ++i)
      {
        const int lid = system.residual[row]->get_map().lid(gids[row][i]);
        rhs[offset[row] + i] = lid < 0 ? 0. : -system.residual[row]->local_values_as_span()[lid];
        if (row < 2 && dbc[row] && dbc[row]->lid(gids[row][i]) >= 0)
          constrained[offset[row] + i] = 1;
      }
    rhs = Core::Communication::sum_all(rhs, comm);
    constrained = Core::Communication::sum_all(constrained, comm);
    for (int col = 0; col < 4; ++col)
      for (int j = 0; j < static_cast<int>(gids[col].size()); ++j)
      {
        const int lid = basis[col]->get_map().lid(gids[col][j]);
        if (lid >= 0) basis[col]->get_values()[lid] = 1.;
        const auto action = apply_linearization(system, basis);
        std::vector<double> column(size, 0.);
        for (int row = 0; row < 4; ++row)
          for (int i = 0; i < static_cast<int>(gids[row].size()); ++i)
          {
            const int local = action[row]->get_map().lid(gids[row][i]);
            if (local >= 0) column[offset[row] + i] = action[row]->local_values_as_span()[local];
          }
        column = Core::Communication::sum_all(column, comm);
        for (int i = 0; i < size; ++i) matrix(i, offset[col] + j) = column[i];
        if (lid >= 0) basis[col]->get_values()[lid] = 0.;
      }
    for (int i = 0; i < size; ++i)
      if (constrained[i])
      {
        for (int j = 0; j < size; ++j) matrix(i, j) = matrix(j, i) = 0.;
        matrix(i, i) = 1.;
        rhs[i] = 0.;
      }
    Core::LinAlg::SerialDenseSolver solver;
    solver.set_matrix(matrix);
    solver.set_vectors(solution, rhs);
    solver.factor_with_equilibration(true);
    FOUR_C_ASSERT_ALWAYS(solver.solve() == 0, "Uncondensed TSI verification solve failed");
    for (int block = 0; block < 4; ++block)
      for (int i = 0; i < static_cast<int>(gids[block].size()); ++i)
      {
        const int lid = result[block]->get_map().lid(gids[block][i]);
        if (lid >= 0) result[block]->get_values()[lid] = solution[offset[block] + i];
      }
    return result;
  }

  // Supply already integrated mortar data to the production coupled assembly.
  // The interface tests exercise geometry derivatives; integration regressions use the real mesh.
  class AssembledTSIStrategy : public CONTACT::LagrangeStrategyTsi
  {
   public:
    AssembledTSIStrategy(const Teuchos::ParameterList& params,
        const std::shared_ptr<CONTACT::TSIInterface>& interface,
        const std::shared_ptr<Core::LinAlg::Map>& dis,
        const std::shared_ptr<Core::LinAlg::Map>& nodes,
        const std::shared_ptr<Core::LinAlg::Map>& source,
        const std::shared_ptr<Core::LinAlg::Map>& target,
        const std::shared_ptr<Core::LinAlg::Map>& thermal, double alpha, double theta, bool active,
        double mortar_weight = 2.)
        : CONTACT::LagrangeStrategyTsi(std::make_shared<CONTACT::AbstractStrategyDataContainer>(),
              dis.get(), nodes.get(), params, {interface}, 3, MPI_COMM_WORLD, alpha, 7)
    {
      gdisprowmap_ = gstdofrowmap_ = dis;
      gsdofrowmap_ = source;
      gactivedofs_ = active ? source : std::make_shared<Core::LinAlg::Map>(0, 0, MPI_COMM_WORLD);
      gtdofrowmap_ = target;
      const int node = 0;
      gsnoderowmap_ = std::make_shared<Core::LinAlg::Map>(1, 1, &node, 0, MPI_COMM_WORLD);
      gactivenodes_ = gactiven_ =
          active ? gsnoderowmap_ : std::make_shared<Core::LinAlg::Map>(0, 0, MPI_COMM_WORLD);
      dmatrix_ = std::make_shared<Core::LinAlg::SparseMatrix>(*source, 1);
      mmatrix_ = std::make_shared<Core::LinAlg::SparseMatrix>(*source, 1);
      for (int d = 0; d < 3; ++d)
      {
        dmatrix_->assemble(mortar_weight, d, d);
        mmatrix_->assemble(mortar_weight, d, d + 3);
      }
      dmatrix_->complete(*source, *source);
      mmatrix_->complete(*target, *source);
      tsi_alpha_ = theta;
      fscn_ = std::make_shared<Core::LinAlg::Vector<double>>(*dis, true);
      ftcn_ = std::make_shared<Core::LinAlg::Vector<double>>(*thermal, true);
      fscn_->put_scalar(0.1);
      ftcn_->put_scalar(0.2);
    }

    void set_coupled_multipliers(const Core::LinAlg::Vector<double>& mechanical,
        const Core::LinAlg::Vector<double>& thermal, Coupling::Adapter::Coupling& coupling)
    {
      z_ = std::make_shared<Core::LinAlg::Vector<double>>(mechanical);
      z_thermo_ = std::make_shared<Core::LinAlg::Vector<double>>(thermal);
      store_nodal_quantities(Mortar::StrategyBase::lmupdate, coupling);
      store_nodal_quantities(Mortar::StrategyBase::lmThermo, coupling);
    }

    std::shared_ptr<const Core::LinAlg::Vector<double>> thermal_multiplier() const
    {
      return z_thermo_;
    }

   protected:
    void prepare_coupled_contact(const Core::LinAlg::Vector<double>&,
        const Core::LinAlg::Vector<double>&, Coupling::Adapter::Coupling&) override
    {
    }
  };

  void TSIInterfaceTest::check_coupled_strategy(
      bool finite_difference, bool active, bool equilibrium)
  {
    auto map = [](std::initializer_list<int> ids)
    {
      return std::make_shared<Core::LinAlg::Map>(
          ids.size(), ids.size(), ids.begin(), 0, MPI_COMM_WORLD);
    };
    auto source_dofs = map({0, 1, 2}), target_dofs = map({3, 4, 5});
    auto thermal = map({6, 7}), thermal_source = map({6}), nodes = map({0, 1});
    auto full = map({0, 1, 2, 3, 4, 5, 6, 7}), empty = map({});
    auto coupling = std::make_shared<Coupling::Adapter::Coupling>();
    auto thermal_on_structure = map({0, 3});
    coupling->setup_coupling(thermal, thermal, thermal_on_structure, thermal_on_structure);
    data->s_node_row_map() = data->slave_node_col_map() = data->active_nodes();
    data->slave_dof_row_map() = data->slave_dof_col_map() = source_dofs;
    data->master_dof_row_map() = data->master_dof_col_map() = target_dofs;
    data->active_n() = data->active_nodes();
    data->active_t() = data->slip_t();
    data->i_mortar().set("HEATTRANSSLAVE", 1.).set("HEATTRANSMASTER", 2.);
    data->i_mortar().set("LM_QUAD", Mortar::lagmult_undefined);
    source->data().get_deriv_d()[0];
    source->data().get_deriv_m()[1];
    source->mo_data().get_m()[1] = 2.;
    source->data().getg() = -0.1;
    source->data().get_deriv_g()[2] = 2.;
    source->data().get_deriv_g()[5] = -2.;
    for (int d = 0; d < 2; ++d)
    {
      source->fri_data().get_deriv_jump()[d][d] = 2.;
      source->fri_data().get_deriv_jump()[d][d + 3] = -2.;
    }
    source->tsi_data().temp() = 20.;
    source->tsi_data().temp_target() = 40.;
    source->tsi_data().deriv_temp_target_temp()[3] = 1.;
    auto* target = dynamic_cast<CONTACT::Node*>(interface->discret().g_node(1));
    target->initialize_tsi_data_container(0., 100.);
    target->tsi_data().temp() = 40.;
    if (!active)
    {
      data->active_nodes() = data->slip_nodes() = data->active_n() = data->active_t() =
          data->slip_t() = empty;
    }
    if (equilibrium)
    {
      // Pythagorean tangential traction and slip give an exact Coulomb equilibrium.
      // A non-binary mortar inverse exposes cancellation when absolute forces
      // are projected onto the slip equations before their residual is formed.
      source->mo_data().get_d()[0] = source->mo_data().get_m()[1] = 3.;
      source->data().getg() = 0.;
      source->fri_data().jump()[0] = 3.;
      source->fri_data().jump()[1] = 4.;
      source->tsi_data().temp() = source->tsi_data().temp_target() = 0.;
      target->tsi_data().temp() = 0.;
      data->i_mortar().set("HEATTRANSMASTER", 1.).set("FRCOEFF", 0.5);
    }
    auto params = data->i_mortar();
    params.set("CONDENSED_LM_INCREMENTS", true);
    params.set("SYSTEM", CONTACT::SystemType::condensed);
    params.set("STRATEGY", CONTACT::SolvingStrategy::lagmult);
    params.sublist("PARALLEL REDISTRIBUTION")
        .set("PARALLEL_REDIST", Mortar::ParallelRedist::redist_none);
    Core::LinAlg::MultiMapExtractor extractor(*full, {dofs, thermal});
    // Nonzero bulk coupling and nontrivial time weights exercise all elimination terms.
    for (const auto weights : {std::array{0., 1.}, std::array{0.25, 0.6}})
    {
      if (equilibrium && weights[0] != 0.) continue;
      AssembledTSIStrategy strategy(params, interface, dofs, nodes, source_dofs, target_dofs,
          thermal, weights[0], weights[1], active, equilibrium ? 3. : 2.);
      EXPECT_EQ(strategy.mech_contact_res_, 0.);
      EXPECT_EQ(strategy.mech_contact_incr_, 0.);
      EXPECT_EQ(strategy.thermo_contact_incr_, 0.);
      Core::LinAlg::Vector<double> lm(*source_dofs, true), q(*thermal_source, true);
      for (int d = 0; d < 3; ++d) lm.get_values()[d] = std::array{-2., -3., 5.}[d];
      q.get_values()[0] = 1.7;
      if (equilibrium)
      {
        for (int d = 0; d < 3; ++d) lm.get_values()[d] = std::array{30000., 40000., 100000.}[d];
        q.get_values()[0] = -250000. / 3.;
      }
      strategy.set_coupled_multipliers(lm, q, *coupling);
      auto matrix = std::make_shared<
          Core::LinAlg::BlockSparseMatrix<Core::LinAlg::DefaultBlockMatrixStrategy>>(
          extractor, extractor, 8, true, false);
      for (int i = 0; i < 6; ++i)
      {
        matrix->matrix(0, 0).assemble(10. + i, i, i);
        for (int j = 0; j < 2; ++j)
        {
          matrix->matrix(0, 1).assemble(0.03 * (i + 1) * (j + 1), i, 6 + j);
          matrix->matrix(1, 0).assemble(-0.02 * (i + 1) * (j + 1), 6 + j, i);
        }
      }
      for (int i = 0; i < 2; ++i) matrix->matrix(1, 1).assemble(4. + i, 6 + i, 6 + i);
      matrix->complete();
      auto rhs = std::make_shared<Core::LinAlg::Vector<double>>(*full, true);
      for (int i = 0; i < 8; ++i) rhs->get_values()[i] = 0.2 * (i + 1);
      if (equilibrium)
      {
        for (int d = 0; d < 3; ++d)
        {
          rhs->get_values()[d] = 3. * lm.local_values_as_span()[d];
          rhs->get_values()[d + 3] = -rhs->get_values()[d];
        }
        rhs->get_values()[6] = rhs->get_values()[7] = -250000.;
      }
      auto dis = std::make_shared<Core::LinAlg::Vector<double>>(*dofs, true);
      auto temp = std::make_shared<Core::LinAlg::Vector<double>>(*thermal, true);
      std::shared_ptr<Core::LinAlg::BlockSparseMatrixBase> bulk =
          matrix->clone(Core::LinAlg::DataAccess::Copy);
      bulk->complete();
      const auto bulk_rhs = std::make_shared<Core::LinAlg::Vector<double>>(*rhs);
      if (equilibrium)
      {
        // Balanced bulk/contact loads and zero contact constraints need no Newton step.
        strategy.evaluate(matrix, rhs, coupling, dis, temp);
        double norm = 0.;
        rhs->norm_2(&norm);
        EXPECT_LT(norm, 1.e-8);
        continue;
      }
      CONTACT::TsiContactLinearization uncondensed;
      strategy.evaluate(matrix, rhs, coupling, dis, temp, &uncondensed);
      for (int d = 0; d < 3; ++d)
      {
        EXPECT_NEAR(uncondensed.residual[0]->local_values_as_span()[d],
            -0.2 * (d + 1) + weights[0] * 0.1 +
                (active ? 1. - weights[0] : 0.) * 2. * lm.local_values_as_span()[d],
            1.e-13);
        EXPECT_NEAR(uncondensed.residual[0]->local_values_as_span()[d + 3],
            -0.2 * (d + 4) + weights[0] * 0.1 -
                (active ? 1. - weights[0] : 0.) * 2. * lm.local_values_as_span()[d],
            1.e-13);
      }
      if (finite_difference)
      {
        CONTACT::TsiContactLinearization::Vectors direction;
        for (int block = 0; block < 4; ++block)
          direction[block] = std::make_shared<Core::LinAlg::Vector<double>>(
              uncondensed.residual[block]->get_map(), true);
        // Every column, including temperature-sensitive friction, rather than a
        // single random direction. Old time forces stay fixed in every probe.
        for (int col = 0; col < 4; ++col)
          for (int local = 0; local < direction[col]->local_length(); ++local)
          {
            direction[col]->get_values()[local] = 1.;
            const auto exact = apply_linearization(uncondensed, direction);
            auto probe = [&](double h)
            {
              Core::LinAlg::Vector<double> perturbed_lm(lm), perturbed_q(q);
              perturbed_lm.update(h, *direction[2], 1.);
              perturbed_q.update(h, *direction[3], 1.);
              strategy.set_coupled_multipliers(perturbed_lm, perturbed_q, *coupling);
              const auto u = direction[0]->local_values_as_span();
              const auto t = direction[1]->local_values_as_span();
              source->data().getg() = -0.1 + h * 2. * (u[2] - u[5]);
              source->fri_data().jump()[0] = 0.4 + h * 2. * (u[0] - u[3]);
              source->fri_data().jump()[1] = 0.2 + h * 2. * (u[1] - u[4]);
              source->tsi_data().temp() = 20. + h * t[0];
              source->tsi_data().temp_target() = target->tsi_data().temp() = 40. + h * t[1];
              std::shared_ptr<Core::LinAlg::BlockSparseMatrixBase> perturbed_matrix =
                  bulk->clone(Core::LinAlg::DataAccess::Copy);
              perturbed_matrix->complete();
              auto perturbed_rhs = std::make_shared<Core::LinAlg::Vector<double>>(*bulk_rhs);
              for (int row = 0; row < 2; ++row)
              {
                Core::LinAlg::Vector<double> residual_change(
                    uncondensed.residual[row]->get_map(), true);
                for (int column = 0; column < 2; ++column)
                {
                  Core::LinAlg::Vector<double> term(residual_change.get_map(), true);
                  bulk->matrix(row, column).multiply(false, *direction[column], term);
                  residual_change.update(-h, term, 1.);
                }
                CONTACT::Utils::add_vector(residual_change, *perturbed_rhs);
              }
              CONTACT::TsiContactLinearization snapshot;
              strategy.evaluate(perturbed_matrix, perturbed_rhs, coupling, dis, temp, &snapshot);
              return snapshot;
            };
            const auto plus = probe(eps), minus = probe(-eps);
            probe(0.);
            for (int row = 0; row < 4; ++row)
            {
              Core::LinAlg::Vector<double> error(*plus.residual[row]);
              error.update(-1., *minus.residual[row], 1.);
              error.scale(0.5 / eps);
              error.update(-1., *exact[row], 1.);
              double norm = 0.;
              error.norm_2(&norm);
              EXPECT_LT(norm, 2.e-7) << "block " << row << ',' << col << ", column " << local;
            }
            direction[col]->get_values()[local] = 0.;
          }
        continue;
      }
      const std::array<std::shared_ptr<const Core::LinAlg::Map>, 2> dbc{map({3}), empty};
      const auto expected = solve_linearization(uncondensed, dbc);
      CONTACT::TsiContactLinearization condensed;
      for (int row = 0; row < 2; ++row)
      {
        condensed.residual[row] = extractor.extract_vector(*rhs, row);
        condensed.residual[row]->scale(-1.);
        for (int col = 0; col < 2; ++col)
          condensed.matrix[row][col] =
              std::make_shared<Core::LinAlg::SparseMatrix>(matrix->matrix(row, col));
      }
      condensed.residual[2] = std::make_shared<Core::LinAlg::Vector<double>>(*empty, true);
      condensed.residual[3] = std::make_shared<Core::LinAlg::Vector<double>>(*empty, true);
      auto actual = solve_linearization(condensed, dbc);
      strategy.recover_coupled(actual[0], actual[1], coupling);
      actual[2] = std::make_shared<Core::LinAlg::Vector<double>>(*strategy.lagrange_multiplier());
      actual[3] = std::make_shared<Core::LinAlg::Vector<double>>(*strategy.thermal_multiplier());
      if (!active)
      {
        double norm = 0.;
        actual[2]->norm_2(&norm);
        EXPECT_EQ(norm, 0.);
        actual[3]->norm_2(&norm);
        EXPECT_EQ(norm, 0.);
        for (int i = 0; i < 2; ++i)
          EXPECT_NEAR(uncondensed.residual[1]->local_values_as_span()[i],
              -0.2 * (7 + i) + (1. - weights[1]) * 0.2, 1.e-13);
      }
      actual[2]->update(-1., lm, 1.);
      actual[3]->update(-1., q, 1.);
      for (int block = 0; block < 4; ++block)
      {
        Core::LinAlg::Vector<double> difference(expected[block]->get_map(), true);
        Core::LinAlg::export_to(*actual[block], difference);
        difference.update(-1., *expected[block], 1.);
        double error = 0.;
        difference.norm_2(&error);
        EXPECT_LT(error, 1.e-10) << "block " << block;
      }
    }
  }

  TEST_F(TSIInterfaceTest, CondensationPreservesHighPressureEquilibrium)
  {
    check_coupled_strategy(false, true, true);
  }

  TEST_F(TSIInterfaceTest, AssembledCoupledJacobianMatchesEveryFiniteDifferenceColumn)
  {
    check_coupled_strategy(true);
  }

  TEST_F(TSIInterfaceTest, CondensationAndRecoveryMatchFullMultiplierSolve)
  {
    check_coupled_strategy(false);
  }

}  // namespace
