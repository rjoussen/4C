// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <gtest/gtest.h>

#include "4C_structure_new_mesh_regularization_sliding.hpp"

#include "4C_fem_discretization.hpp"
#include "4C_fem_general_node.hpp"
#include "4C_io_pstream.hpp"
#include "4C_structure_new_mesh_regularization_ghosting.hpp"
#include "4C_unittest_utils_create_discretization_helper_test.hpp"

#include <vector>

namespace
{
  using namespace FourC;

  // unit cube with 8x8x8 elements, node id = 81 i + 9 j + k at (i, j, k) / 8
  constexpr int subdivisions = 8;

  int node_id(int i, int j, int k)
  {
    return i * (subdivisions + 1) * (subdivisions + 1) + j * (subdivisions + 1) + k;
  }

  class BoundarySlidingParallelTest : public testing::Test
  {
   public:
    BoundarySlidingParallelTest()
    {
      comm_ = MPI_COMM_WORLD;
      discretization_ = std::make_shared<Core::FE::Discretization>("dummy", comm_, 3);
      Core::IO::cout.setup(false, false, false, Core::IO::standard, comm_, 0, 0, "dummyFilePrefix");
      TESTING::fill_discretization_hyper_cube(*discretization_, subdivisions, comm_);

      // two slabs: rank 0 owns the nodes with i <= 4, rank 1 the nodes with i >= 5
      const int rank = Core::Communication::my_mpi_rank(comm_);
      std::vector<int> row_nodes;
      std::vector<int> col_nodes;
      for (int i = 0; i <= subdivisions; ++i)
      {
        for (int j = 0; j <= subdivisions; ++j)
        {
          for (int k = 0; k <= subdivisions; ++k)
          {
            if ((rank == 0) == (i <= 4)) row_nodes.push_back(node_id(i, j, k));
            if (rank == 0 ? i <= 5 : i >= 4) col_nodes.push_back(node_id(i, j, k));
          }
        }
      }
      const Core::LinAlg::Map node_row_map(
          -1, static_cast<int>(row_nodes.size()), row_nodes.data(), 0, comm_);
      const Core::LinAlg::Map node_col_map(
          -1, static_cast<int>(col_nodes.size()), col_nodes.data(), 0, comm_);
      discretization_->redistribute({node_row_map, node_col_map});

      // as for the solid discretization of the mesh regularization
      Solid::MeshRegularization::extend_ghosting(
          *discretization_, Solid::MeshRegularization::num_ghost_layers);
    }

    void TearDown() override { Core::IO::cout.close(); }

   protected:
    std::shared_ptr<Core::FE::Discretization> discretization_;
    MPI_Comm comm_;
  };

  TEST_F(BoundarySlidingParallelTest, SlidesOverTheFacesOfOtherProcessors)
  {
    Solid::MeshRegularization::BoundarySliding sliding(*discretization_, 30.0);
    sliding.setup();

    // a node of rank 0 on the face z = 0 at the processor boundary, at (0.5, 0.5, 0); it slides on
    // the face z = 0 over the faces of the elements of rank 1, which are fetched from rank 1
    const int moved_node = node_id(4, 4, 0);
    const bool owned = discretization_->node_row_map()->lid(moved_node) >= 0;
    std::vector<int> lids;
    if (owned)
      for (int dof : discretization_->dof(0, discretization_->g_node(moved_node)))
        lids.push_back(discretization_->dof_row_map()->lid(dof));
    const Core::LinAlg::Map no_dirichlet_dofs(0, 0, comm_);

    const auto project = [&](double dx, double dy)
    {
      Core::LinAlg::Vector<double> displacement(*discretization_->dof_row_map(), true);
      if (owned)
      {
        displacement.get_values()[lids[0]] = dx;
        displacement.get_values()[lids[1]] = dy;
      }
      sliding.project_onto_boundary(displacement, no_dirichlet_dofs);
      std::vector<double> values;
      if (owned)
        for (int lid : lids) values.push_back(displacement.local_values_as_span()[lid]);
      return values;
    };

    // far into the faces of rank 1, beyond the elements ghosted on rank 0: the movement is kept
    const std::vector<double> far = project(0.45, -0.1);
    if (owned)
    {
      EXPECT_NEAR(far[0], 0.45, 1.0e-12);
      EXPECT_NEAR(far[1], -0.1, 1.0e-12);
      EXPECT_NEAR(far[2], 0.0, 1.0e-12);
    }

    // beyond the face z = 0: the node stops at its feature edge at x = 1 instead of sliding onto
    // the face x = 1, which belongs to another patch
    const std::vector<double> beyond = project(0.6, -0.1);
    if (owned)
    {
      EXPECT_NEAR(beyond[0], 0.5, 1.0e-12);
      EXPECT_NEAR(beyond[1], -0.1, 1.0e-12);
      EXPECT_NEAR(beyond[2], 0.0, 1.0e-12);
    }
  }
}  // namespace
