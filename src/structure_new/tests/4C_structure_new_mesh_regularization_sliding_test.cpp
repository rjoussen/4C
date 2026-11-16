// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <gtest/gtest.h>

#include "4C_structure_new_mesh_regularization_sliding.hpp"

#include "4C_fem_discretization.hpp"
#include "4C_io_pstream.hpp"
#include "4C_unittest_utils_create_discretization_helper_test.hpp"

#include <array>

namespace
{
  using namespace FourC;

  class BoundarySlidingTest : public testing::Test
  {
   public:
    BoundarySlidingTest()
    {
      comm_ = MPI_COMM_WORLD;
      discretization_ = std::make_shared<Core::FE::Discretization>("dummy", comm_, 3);
      Core::IO::cout.setup(false, false, false, Core::IO::standard, comm_, 0, 0, "dummyFilePrefix");
      // unit cube with 2x2x2 elements, node id = 9 i + 3 j + k at (i, j, k) / 2
      TESTING::fill_discretization_hyper_cube(*discretization_, 2, comm_);
    }

    void TearDown() override { Core::IO::cout.close(); }

    //! displacement with the given movement of one node
    Core::LinAlg::Vector<double> move_node(int node_id, double dx, double dy, double dz) const
    {
      Core::LinAlg::Vector<double> displacement(*discretization_->dof_row_map(), true);
      const std::vector<int> dofs = discretization_->dof(0, discretization_->g_node(node_id));
      const double movement[3] = {dx, dy, dz};
      for (int d = 0; d < 3; ++d)
        displacement.get_values()[discretization_->dof_row_map()->lid(dofs[d])] = movement[d];
      return displacement;
    }

    double value(const Core::LinAlg::Vector<double>& displacement, int node_id, int d) const
    {
      const std::vector<int> dofs = discretization_->dof(0, discretization_->g_node(node_id));
      return displacement.local_values_as_span()[discretization_->dof_row_map()->lid(dofs[d])];
    }

   protected:
    std::shared_ptr<Core::FE::Discretization> discretization_;
    MPI_Comm comm_;
  };

  TEST_F(BoundarySlidingTest, ConstrainsNormalMovementOnly)
  {
    if (Core::Communication::num_mpi_ranks(comm_) > 1) GTEST_SKIP();

    Solid::MeshRegularization::BoundarySliding sliding(*discretization_, 30.0);
    sliding.setup();

    // squared length of the tangential part of a movement of a node
    const auto tangential =
        [&](int node_id, double dx, double dy, double dz, const std::array<bool, 3>& fixed = {})
    {
      Core::LinAlg::Matrix<3, 1> movement;
      movement(0) = dx;
      movement(1) = dy;
      movement(2) = dz;
      Core::LinAlg::Matrix<3, 1> projected;
      projected.multiply(sliding.tangent_projection(node_id, fixed), movement);
      return projected.dot(projected);
    };

    constexpr int interior = 13;
    constexpr int face = 4;    // center of the face x = 0
    constexpr int edge = 1;    // middle of the edge x = y = 0
    constexpr int corner = 0;  // corner x = y = z = 0

    EXPECT_DOUBLE_EQ(tangential(interior, 0.1, 0.2, 0.3), 0.14);
    EXPECT_DOUBLE_EQ(tangential(interior, 0.1, 0.2, 0.3, {false, true, false}), 0.1);

    EXPECT_DOUBLE_EQ(tangential(face, 0.0, 0.1, -0.2), 0.05);
    EXPECT_DOUBLE_EQ(tangential(face, 0.1, 0.0, 0.0), 0.0);

    EXPECT_DOUBLE_EQ(tangential(edge, 0.0, 0.0, 0.1), 0.01);
    EXPECT_DOUBLE_EQ(tangential(edge, 0.1, 0.1, 0.0), 0.0);

    EXPECT_DOUBLE_EQ(tangential(corner, 0.1, 0.2, 0.3), 0.0);
  }

  TEST_F(BoundarySlidingTest, KeepsFeaturesOfTheInitialMesh)
  {
    if (Core::Communication::num_mpi_ranks(comm_) > 1) GTEST_SKIP();

    Solid::MeshRegularization::BoundarySliding sliding(*discretization_, 30.0);

    // kink the face x = 0 by moving its center node outwards: the adjacent faces enclose 77 degrees
    discretization_->l_col_node(discretization_->node_col_map()->lid(4))->set_pos({-0.4, 0.5, 0.5});
    sliding.setup();

    // the center node still slides on its face instead of becoming a corner: by symmetry, its
    // averaged normal is the x direction
    const Core::LinAlg::Matrix<3, 3> projection = sliding.tangent_projection(4);
    EXPECT_NEAR(projection(0, 0), 0.0, 1.0e-14);
    EXPECT_NEAR(projection(1, 1), 1.0, 1.0e-14);
    EXPECT_NEAR(projection(2, 2), 1.0, 1.0e-14);
  }

  TEST_F(BoundarySlidingTest, ProjectsOntoBoundary)
  {
    if (Core::Communication::num_mpi_ranks(comm_) > 1) GTEST_SKIP();

    Solid::MeshRegularization::BoundarySliding sliding(*discretization_, 30.0);
    sliding.setup();

    const Core::LinAlg::Map no_dirichlet_dofs(0, 0, comm_);

    // a face node moved off the face is projected back onto it, keeping its tangential movement
    Core::LinAlg::Vector<double> displacement = move_node(4, 0.05, 0.1, -0.2);
    sliding.project_onto_boundary(displacement, no_dirichlet_dofs);
    EXPECT_NEAR(value(displacement, 4, 0), 0.0, 1.0e-12);
    EXPECT_NEAR(value(displacement, 4, 1), 0.1, 1.0e-12);
    EXPECT_NEAR(value(displacement, 4, 2), -0.2, 1.0e-12);

    // an edge node moved off the edge is projected back onto it
    displacement = move_node(1, -0.05, 0.03, 0.2);
    sliding.project_onto_boundary(displacement, no_dirichlet_dofs);
    EXPECT_NEAR(value(displacement, 1, 0), 0.0, 1.0e-12);
    EXPECT_NEAR(value(displacement, 1, 1), 0.0, 1.0e-12);
    EXPECT_NEAR(value(displacement, 1, 2), 0.2, 1.0e-12);

    // corners are fixed
    displacement = move_node(0, 0.05, 0.03, 0.2);
    sliding.project_onto_boundary(displacement, no_dirichlet_dofs);
    for (int d = 0; d < 3; ++d) EXPECT_EQ(value(displacement, 0, d), 0.0);

    // Dirichlet dofs keep their values
    const int dirichlet_dof = discretization_->dof(0, discretization_->g_node(4))[0];
    const Core::LinAlg::Map dirichlet_dofs(-1, 1, &dirichlet_dof, 0, comm_);
    displacement = move_node(4, 0.05, 0.1, -0.2);
    sliding.project_onto_boundary(displacement, dirichlet_dofs);
    EXPECT_EQ(value(displacement, 4, 0), 0.05);
    EXPECT_NEAR(value(displacement, 4, 1), 0.1, 1.0e-12);
  }

  TEST(BoundarySlidingWalkTest, WalksAlongTheBoundary)
  {
    const MPI_Comm comm = MPI_COMM_WORLD;
    if (Core::Communication::num_mpi_ranks(comm) > 1) GTEST_SKIP();

    Core::IO::cout.setup(false, false, false, Core::IO::standard, comm, 0, 0, "dummyFilePrefix");
    // unit cube with 4x4x4 elements, node id = 25 i + 5 j + k at (i, j, k) / 4
    Core::FE::Discretization discretization("dummy", comm, 3);
    TESTING::fill_discretization_hyper_cube(discretization, 4, comm);

    Solid::MeshRegularization::BoundarySliding sliding(discretization, 30.0);
    sliding.setup();
    const Core::LinAlg::Map no_dirichlet_dofs(0, 0, comm);

    // the node at (0.25, 0.5, 0) slides on the face z = 0
    constexpr int node = 35;
    const std::vector<int> dofs = discretization.dof(0, discretization.g_node(node));
    const auto project = [&](double dx, double dy)
    {
      Core::LinAlg::Vector<double> displacement(*discretization.dof_row_map(), true);
      displacement.get_values()[discretization.dof_row_map()->lid(dofs[0])] = dx;
      displacement.get_values()[discretization.dof_row_map()->lid(dofs[1])] = dy;
      sliding.project_onto_boundary(displacement, no_dirichlet_dofs);
      return displacement;
    };
    const auto value = [&](const Core::LinAlg::Vector<double>& displacement, int d)
    { return displacement.local_values_as_span()[discretization.dof_row_map()->lid(dofs[d])]; };

    // more than two elements away from the faces around the node
    Core::LinAlg::Vector<double> displacement = project(0.6, 0.0);
    EXPECT_NEAR(value(displacement, 0), 0.6, 1.0e-12);
    EXPECT_NEAR(value(displacement, 1), 0.0, 1.0e-12);

    // back and across the edge of the cube, a feature edge: projected onto the edge
    displacement = project(0.0, 0.7);
    EXPECT_NEAR(value(displacement, 0), 0.0, 1.0e-12);
    EXPECT_NEAR(value(displacement, 1), 0.5, 1.0e-12);

    Core::IO::cout.close();
  }
}  // namespace
