// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <gtest/gtest.h>

#include "4C_structure_new_mesh_regularization_sliding.hpp"

#include "4C_fem_discretization.hpp"
#include "4C_global_data.hpp"
#include "4C_io_pstream.hpp"
#include "4C_unittest_utils_create_discretization_helper_test.hpp"

namespace
{
  using namespace FourC;

  TEST(BoundarySlidingTest, GroupsTheBoundaryFacesIntoPatches)
  {
    const MPI_Comm comm = MPI_COMM_WORLD;
    if (Core::Communication::num_mpi_ranks(comm) > 1) GTEST_SKIP();

    Core::IO::cout.setup(false, false, false, Core::IO::standard, comm, 0, 0, "dummyFilePrefix");
    // unit cube with 2x2x2 elements, node id = 9 i + 3 j + k at (i, j, k) / 2
    auto discretization = std::make_shared<Core::FE::Discretization>("dummy", comm, 3);
    TESTING::fill_discretization_hyper_cube(*discretization, 2, comm);

    const Solid::MeshRegularization::BoundarySliding sliding(
        *Global::Problem::instance(), discretization, 30.0);

    constexpr int interior = 13;
    constexpr int face = 4;    // center of the face x = 0
    constexpr int edge = 1;    // middle of the edge x = y = 0
    constexpr int corner = 0;  // corner x = y = z = 0
    EXPECT_TRUE(sliding.patches(interior).empty());
    EXPECT_EQ(sliding.patches(face).size(), 1);
    EXPECT_EQ(sliding.patches(edge).size(), 2);
    EXPECT_EQ(sliding.patches(corner).size(), 3);

    // the faces x = 0 and y = 0 are the patches at the edge and the corner
    const std::vector<int> face_patch = sliding.patches(face);
    const std::vector<int> edge_patches = sliding.patches(edge);
    EXPECT_NE(std::ranges::find(edge_patches, face_patch.front()), edge_patches.end());
    for (int patch : edge_patches)
    {
      const std::vector<int> corner_patches = sliding.patches(corner);
      EXPECT_NE(std::ranges::find(corner_patches, patch), corner_patches.end());
    }

    Core::IO::cout.close();
  }
}  // namespace
