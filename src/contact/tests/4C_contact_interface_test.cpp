// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <gtest/gtest.h>

#include "4C_contact_interface.hpp"

#include "4C_contact_friction_node.hpp"
#include "4C_fem_discretization.hpp"
#include "4C_global_data.hpp"
#include "4C_io_control.hpp"
#include "4C_linalg_sparsematrix.hpp"

#include <array>

namespace
{
  using namespace FourC;

  // Exercise the shared mechanical Coulomb assembler without thermal data or coupling.
  // Prescribe nodal contact data and call the slip branch directly, without contact search
  // or active-set selection. Nodal vectors (normal, multiplier, jump) use Cartesian x/y/z.
  // ConstraintDirection changes the assembled equation rows, not those vector components.
  class CoulombInterfaceTest : public ::testing::Test
  {
   protected:
    static constexpr int x = 0, y = 1, z = 2;
    static constexpr std::array<int, 2> tangent_components{x, y};
    static constexpr double friction_coefficient = 0.3;
    static constexpr double normal_weight = 2.;
    static constexpr double tangential_weight = 3.;

    void SetUp() override
    {
      Global::Problem::instance()->set_spatial_approximation_type(
          Core::FE::ShapeFunctionType::polynomial);
      Global::Problem::instance()->set_output_control_file(
          std::make_shared<Core::IO::OutputControl>(MPI_COMM_WORLD, "contact",
              Core::FE::ShapeFunctionType::polynomial, "", "", "coulomb_interface_test", 3, 0, 1,
              false, false));

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
      params.set("PROBTYPE", CONTACT::Problemtype::structure);
      params.set("FRCOEFF", friction_coefficient);
      params.set("FRBOUND", 0.0);
      params.set("GP_SLIP_INCR", false);
      params.set("FRLESS_FIRST", false);
      params.set("SEMI_SMOOTH_CN", normal_weight);
      params.set("SEMI_SMOOTH_CT", tangential_weight);

      data = std::make_shared<CONTACT::InterfaceDataContainer>();
      interface = std::make_shared<CONTACT::Interface>(data, 0, MPI_COMM_WORLD, 3, params, false);
      const std::array<double, 3> coordinates{0., 0., 0.};
      source = std::make_shared<CONTACT::FriNode>(
          0, coordinates, 0, std::vector<int>{0, 1, 2}, true, true, true);
      interface->add_mortar_node(source);
      // Attach the nodes to the discretization before initializing their contact data.
      // No elements are needed: these tests supply already integrated mortar entries.
      interface->discret().fill_complete(Core::FE::OptionsFillComplete::none());
      source->initialize_data_container();
      // Choose n=e_z, t1=e_x, t2=e_y. Thus lambda_z is normal, lambda_x/y tangential.
      // This deliberately differs from the local equation order (n,t1,t2), so the two
      // coordinate checks also exercise the mapping from tangents to Cartesian rows.
      source->mo_data().n()[z] = 1.;
      source->data().txi()[x] = 1.;
      source->data().teta()[y] = 1.;
      source->mo_data().lm()[x] = -2.;
      source->mo_data().lm()[y] = -3.;
      source->mo_data().lm()[z] = 5.;
      source->fri_data().jump()[x] = 0.4;
      source->fri_data().jump()[y] = 0.2;
      source->data().get_deriv_n().resize(3);
      source->data().get_deriv_txi().resize(3);
      source->data().get_deriv_teta().resize(3);
      source->fri_data().get_deriv_jump().resize(3);

      const int node_id = 0;
      data->active_nodes() = std::make_shared<Core::LinAlg::Map>(1, 1, &node_id, 0, MPI_COMM_WORLD);
      data->slip_nodes() = data->active_nodes();
      // In local n/t/t constraint coordinates, rows 1 and 2 contain the two slip equations.
      // They are constraint row labels, not the Cartesian directions of the nodal normal.
      const std::array<int, 2> local_tangent_rows{1, 2};
      data->slip_t() =
          std::make_shared<Core::LinAlg::Map>(2, 2, local_tangent_rows.data(), 0, MPI_COMM_WORLD);
      dofs = std::make_shared<Core::LinAlg::Map>(3, 0, MPI_COMM_WORLD);
      data->cn_values() = std::make_shared<Core::LinAlg::Vector<double>>(*data->active_nodes());
      data->ct_values() = std::make_shared<Core::LinAlg::Vector<double>>(*data->active_nodes());
      data->cn_values()->put_scalar(normal_weight);
      data->ct_values()->put_scalar(tangential_weight);
    }

    void TearDown() override { Global::Problem::instance()->set_output_control_file(nullptr); }

    /// Extract column c as A e_c, including completion of the finite-element matrix assembly.
    Core::LinAlg::Vector<double> column(Core::LinAlg::SparseMatrix& matrix, int col) const
    {
      matrix.complete(*dofs, *dofs);
      Core::LinAlg::Vector<double> basis(*dofs, true), result(*dofs, true);
      basis.get_values()[col] = 1.;
      matrix.multiply(false, basis, result);
      return result;
    }

    // t=0/1 identifies the first/second tangent, not a Cartesian component.
    // Local rows: (normal,t1,t2)=(0,1,2); Cartesian rows: (t1,t2,normal)=(x,y,z).
    int tangential_equation_row(int t) const
    {
      return data->constraint_direction() == CONTACT::ConstraintDirection::xyz
                 ? tangent_components[t]
                 : t + 1;
    }

    std::shared_ptr<CONTACT::InterfaceDataContainer> data;
    std::shared_ptr<CONTACT::Interface> interface;
    std::shared_ptr<CONTACT::FriNode> source;
    std::shared_ptr<Core::LinAlg::Map> dofs;
  };

  TEST_F(CoulombInterfaceTest, SlipResidualPreservesSmallTransverseMotionAtHighPressure)
  {
    // Place the tangential multiplier on the friction bound: |(30000,40000)|=50000.
    // The small prescribed slip perturbation tests evaluation of this fixed slip branch.
    constexpr double friction_bound = 50000.;
    constexpr std::array<double, 2> slip_increment{3.713e-5, -2.915e-5};
    source->data().getg() = 0.;
    source->mo_data().lm()[x] = 30000.;
    source->mo_data().lm()[y] = 40000.;
    source->mo_data().lm()[z] = friction_bound / friction_coefficient;
    source->fri_data().jump()[x] = slip_increment[0];
    source->fri_data().jump()[y] = slip_increment[1];
    // For a=lambda_t+ct*j_t and b=mu*lambda_n, RHS=-|a|*lambda_t+b*a.
    // Project onto e=(-4/5,3/5), perpendicular to lambda_t=(30000,40000).
    // Since e.lambda_t=0, the exact transverse RHS is simply b*ct*(e.j_t).
    // This supplies an independent reference without evaluating |a| or subtracting
    // O(1e9) terms. The original expanded assembly loses about 1.5e-7 here.
    constexpr std::array<double, 2> transverse_direction{-0.8, 0.6};
    const double transverse_slip =
        transverse_direction[0] * slip_increment[0] + transverse_direction[1] * slip_increment[1];
    const double expected = friction_bound * tangential_weight * transverse_slip;
    for (const auto direction :
        {CONTACT::ConstraintDirection::ntt, CONTACT::ConstraintDirection::xyz})
    {
      SCOPED_TRACE(
          direction == CONTACT::ConstraintDirection::ntt ? "Local rows" : "Cartesian rows");
      data->constraint_direction() = direction;
      Core::LinAlg::SparseMatrix multiplier_tangent(*dofs, 3), displacement_tangent(*dofs, 3);
      Core::LinAlg::Vector<double> rhs(*dofs, true);
      interface->assemble_lin_slip(multiplier_tangent, displacement_tangent, rhs);
      const double transverse_rhs =
          transverse_direction[0] * rhs.local_values_as_span()[tangential_equation_row(0)] +
          transverse_direction[1] * rhs.local_values_as_span()[tangential_equation_row(1)];
      EXPECT_NEAR(transverse_rhs, expected, 1.e-12);
    }
  }

  TEST_F(CoulombInterfaceTest, FrictionlessAndZeroTrialSlipBranches)
  {
    source->data().getg() = 0.;
    for (bool zero_trial : {false, true})
    {
      data->i_mortar().set("FRCOEFF", zero_trial ? friction_coefficient : 0.);
      // Either mu=0, or a=lambda_t+ct*j_t=0. Both branches enforce lambda_t=0.
      if (zero_trial)
        for (const int component : tangent_components)
          source->fri_data().jump()[component] =
              -source->mo_data().lm()[component] / tangential_weight;
      for (const auto direction :
          {CONTACT::ConstraintDirection::ntt, CONTACT::ConstraintDirection::xyz})
      {
        SCOPED_TRACE(
            direction == CONTACT::ConstraintDirection::ntt ? "Local rows" : "Cartesian rows");
        SCOPED_TRACE(zero_trial ? "Zero trial" : "Zero friction");
        data->constraint_direction() = direction;
        Core::LinAlg::SparseMatrix lm(*dofs, 3), disp(*dofs, 3);
        Core::LinAlg::Vector<double> rhs(*dofs, true);
        interface->assemble_lin_slip(lm, disp, rhs);
        for (int t = 0; t < 2; ++t)
        {
          EXPECT_NEAR(rhs.local_values_as_span()[tangential_equation_row(t)],
              -source->mo_data().lm()[tangent_components[t]], 1.e-14);
          EXPECT_NEAR(
              column(lm, tangent_components[t]).local_values_as_span()[tangential_equation_row(t)],
              1., 1.e-14);
        }
      }
    }
  }

}  // namespace
