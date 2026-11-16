// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <gtest/gtest.h>

#include "4C_structure_new_mesh_regularization_sliding.hpp"

#include "4C_comm_mpi_utils.hpp"
#include "4C_fem_discretization.hpp"
#include "4C_global_data.hpp"
#include "4C_io_pstream.hpp"
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

  TEST(BoundarySlidingParallelTest, PatchesSpanProcessors)
  {
    const MPI_Comm comm = MPI_COMM_WORLD;
    Core::IO::cout.setup(false, false, false, Core::IO::standard, comm, 0, 0, "dummyFilePrefix");
    auto discretization = std::make_shared<Core::FE::Discretization>("dummy", comm, 3);
    TESTING::fill_discretization_hyper_cube(*discretization, subdivisions, comm);

    // two slabs: rank 0 owns the nodes with i <= 4, rank 1 the nodes with i >= 5
    const int rank = Core::Communication::my_mpi_rank(comm);
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
        -1, static_cast<int>(row_nodes.size()), row_nodes.data(), 0, comm);
    const Core::LinAlg::Map node_col_map(
        -1, static_cast<int>(col_nodes.size()), col_nodes.data(), 0, comm);
    discretization->redistribute({node_row_map, node_col_map});

    const Solid::MeshRegularization::BoundarySliding sliding(
        *Global::Problem::instance(), discretization, 30.0);

    // nodes of the face z = 0 on both ranks belong to the same patch
    const std::vector<int> patches = sliding.patches(node_id(rank == 0 ? 2 : 6, 4, 0));
    ASSERT_EQ(patches.size(), 1);
    EXPECT_EQ(Core::Communication::max_all(patches.front(), comm),
        -Core::Communication::max_all(-patches.front(), comm));

    Core::IO::cout.close();
  }
}  // namespace
