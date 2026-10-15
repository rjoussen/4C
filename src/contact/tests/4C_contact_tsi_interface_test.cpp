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

  constexpr double dt = 0.5;

  Teuchos::ParameterList contact_params()
  {
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
    return params;
  }

  /**
   * @brief Compare TSI contact Jacobian contributions with central finite differences.
   *
   * One source node (DOFs 0, 1, 2) and one target node (DOFs 3, 4, 5) with prescribed mortar
   * data, so no contact search or integration is needed. As in the interface assembly, the
   * source DOFs also index the three contact Lagrange multipliers, and the first DOF of each
   * node (0 and 3) stands for its temperature and thermal multiplier. D = 2 and the nodal
   * c_n = 2, c_t = 3 deliberately differ from one and from the input values, so that a missing
   * normalization or the wrong parameters show up.
   */
  class TSIInterfaceTest : public ::testing::Test
  {
   protected:
    void SetUp() override
    {
      // the interface constructor needs an output control
      Global::Problem::instance()->set_output_control_file(
          std::make_shared<Core::IO::OutputControl>(MPI_COMM_WORLD, "tsi",
              Core::FE::ShapeFunctionType::polynomial, "", "", "tsi_interface_test", 3, 0, 1, false,
              false));

      interface = std::make_shared<CONTACT::TSIInterface>(
          data, 0, MPI_COMM_WORLD, 3, contact_params(), false);
      interface->add_mortar_node(source);
      interface->add_mortar_node(target);
      // No elements are needed, the mortar entries are prescribed below. The node data
      // containers are sized from the adjacent elements, so initialize them after fill_complete.
      interface->discret().fill_complete(Core::FE::OptionsFillComplete::none());
      source->initialize_data_container();
      source->initialize_tsi_data_container(0., 100.);
      // mortar entries D_00 and M_01, normal, multipliers (lambda = -t) and weighted jump
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
      // one derivative map per spatial component, filled by the tests that need them
      source->data().get_deriv_n().resize(3);
      source->data().get_deriv_txi().resize(3);
      source->data().get_deriv_teta().resize(3);
      source->fri_data().get_deriv_jump().resize(3);

      // the source node is active and slipping
      const int node_id = 0;
      data->active_nodes() = std::make_shared<Core::LinAlg::Map>(1, 1, &node_id, 0, MPI_COMM_WORLD);
      data->slip_nodes() = data->active_nodes();
      // rows of the two tangential slip equations in n/t/t constraint coordinates
      const std::array<int, 2> tangent_dofs{1, 2};
      data->slip_t() =
          std::make_shared<Core::LinAlg::Map>(2, 2, tangent_dofs.data(), 0, MPI_COMM_WORLD);
      data->cn_values() = std::make_shared<Core::LinAlg::Vector<double>>(*data->active_nodes());
      data->ct_values() = std::make_shared<Core::LinAlg::Vector<double>>(*data->active_nodes());
      // nodal parameters, e.g. scaled at contact edges, instead of the input values 1.0
      data->cn_values()->put_scalar(2.);
      data->ct_values()->put_scalar(3.);
    }

    void TearDown() override { Global::Problem::instance()->set_output_control_file(nullptr); }

    /// Signed dissipation term -lambda . (I - n x n) . jump / (dt D) as used in the thermal
    /// residual
    double dissipation()
    {
      const Core::LinAlg::Matrix<3, 1> lm(source->mo_data().lm(), true);
      const Core::LinAlg::Matrix<3, 1> n(source->mo_data().n(), true);
      const Core::LinAlg::Matrix<3, 1> jump(source->fri_data().jump(), true);
      return -(lm.dot(jump) - lm.dot(n) * jump.dot(n)) / (dt * source->mo_data().get_d()[0]);
    }

    /// Column col of the (completed) matrix
    Core::LinAlg::Vector<double> column(Core::LinAlg::SparseMatrix& matrix, int col) const
    {
      matrix.complete(dofs, dofs);
      Core::LinAlg::Vector<double> basis(dofs, true), result(dofs, true);
      basis.get_values()[col] = 1.;
      matrix.multiply(false, basis, result);
      return result;
    }

    static constexpr double eps = 1.e-6;
    static constexpr std::array<double, 3> coordinates{0., 0., 0.};
    std::shared_ptr<CONTACT::InterfaceDataContainer> data =
        std::make_shared<CONTACT::InterfaceDataContainer>();
    std::shared_ptr<CONTACT::FriNode> source = std::make_shared<CONTACT::FriNode>(
        0, coordinates, 0, std::vector<int>{0, 1, 2}, true, true, true);
    std::shared_ptr<CONTACT::FriNode> target = std::make_shared<CONTACT::FriNode>(
        1, coordinates, 0, std::vector<int>{3, 4, 5}, false, false, true);
    std::shared_ptr<CONTACT::TSIInterface> interface;
    Core::LinAlg::Map dofs{6, 0, MPI_COMM_WORLD};
  };

  /**
   * @brief Temperature derivative of the slip condition must use the nodal c_n and c_t.
   *
   * The source is hotter than the projected target, so the friction coefficient depends on
   * the source temperature. The reference is a finite difference of the mechanical slip
   * residual.
   */
  TEST_F(TSIInterfaceTest, SlipTemperatureDerivativeUsesNodalParameters)
  {
    source->tsi_data().temp() = 20.;
    source->tsi_data().temp_target() = 10.;
    source->data().getg() = -0.1;
    Core::LinAlg::SparseMatrix lm(dofs, 6), disp(dofs, 6), temp(dofs, 6);
    Core::LinAlg::Vector<double> rhs(dofs, true);
    interface->assemble_lin_slip(lm, disp, temp, rhs);
    auto tangent = column(temp, 0);

    auto residual = [&](double temperature)
    {
      source->tsi_data().temp() = temperature;
      Core::LinAlg::SparseMatrix lm_dummy(dofs, 6), disp_dummy(dofs, 6);
      Core::LinAlg::Vector<double> result(dofs, true);
      interface->CONTACT::Interface::assemble_lin_slip(lm_dummy, disp_dummy, result);
      result.scale(-1.);  // the assembly returns the negative residual
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
   * @brief Displacement and multiplier derivatives of the thermal contact condition.
   *
   * assemble_lin_conduct linearizes
   * C_T = D q - delta D w - beta (lambda . n) (D T1 - M T2)
   * with w = dissipation(), beta = 2/3 and delta = 1/3. One displacement x changes D, M, the
   * normal and the jump at once, which covers the 1/D normalization of w and its derivative.
   */
  TEST_F(TSIInterfaceTest, ThermalContactJacobian)
  {
    data->i_mortar().set("HEATTRANSSLAVE", 1.).set("HEATTRANSMASTER", 2.);
    data->i_mortar().set("LM_QUAD", Mortar::lagmult_undefined);
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
    // derivatives of the quantities in set_displacement with respect to DOF 0
    source->data().get_deriv_d()[0][0] = 0.7;
    source->data().get_deriv_m()[1][0] = -0.4;
    for (int d = 0; d < 3; ++d)
    {
      source->data().get_deriv_n()[d].resize(1);
      source->data().get_deriv_n()[d][0] = std::array{std::cos(0.3), 0., -std::sin(0.3)}[d];
      source->fri_data().get_deriv_jump()[d][0] = std::array{0.7, -0.3, 0.2}[d];
    }
    Core::LinAlg::SparseMatrix disp(dofs, 6, true, false, Core::LinAlg::SparseMatrix::FE_MATRIX);
    Core::LinAlg::SparseMatrix temp(dofs, 6, true, false, Core::LinAlg::SparseMatrix::FE_MATRIX);
    Core::LinAlg::SparseMatrix q(dofs, 6, true, false, Core::LinAlg::SparseMatrix::FE_MATRIX);
    Core::LinAlg::SparseMatrix lm(dofs, 6, true, false, Core::LinAlg::SparseMatrix::FE_MATRIX);
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
    // compare column col with the central difference of the residual, where perturb(x) sets the
    // varied quantity to its base value plus x
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
    for (int d = 0; d < 3; ++d)
    {
      const double base = source->mo_data().lm()[d];
      check(lm, d, [&](double increment) { source->mo_data().lm()[d] = base + increment; });
    }
  }

}  // namespace
