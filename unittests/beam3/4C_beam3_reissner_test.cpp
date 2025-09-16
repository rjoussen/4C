// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <gtest/gtest.h>

#include "4C_beam3_reissner.hpp"

#include "4C_fem_discretization_nullspace.hpp"
#include "4C_fem_general_element.hpp"
#include "4C_linalg_sparsematrix.hpp"
#include "4C_unittest_utils_assertions_test.hpp"

#include <array>

const double testTolerance = 1e-14;

namespace
{
  using namespace FourC;

  class Beam3r : public ::testing::Test
  {
   public:
    Beam3r()
    {
      testdis_ = std::make_shared<Core::FE::Discretization>("Beam3r", MPI_COMM_WORLD, 3);

      std::vector<std::vector<double>> xrefe{{-0.05, 0.05, 0.3}, {0.45, -0.05, 0.1}, {0.2, 0, 0.2}};
      std::vector<double> xrefe_full{-0.05, 0.05, 0.3, 0.45, -0.05, 0.1};

      testdis_->add_node(std::make_shared<Core::Nodes::Node>(0, xrefe[0], 0));
      testdis_->add_node(std::make_shared<Core::Nodes::Node>(1, xrefe[1], 0));
      testdis_->add_node(std::make_shared<Core::Nodes::Node>(2, xrefe[2], 0));

      testele_ = std::make_shared<Discret::Elements::Beam3r>(0, 0);
      std::array<int, 3> node_ids{0, 1, 2};
      testele_->set_node_ids(3, node_ids.data());

      // setup internal beam element parameters
      std::vector<double> rotrefe(9);
      rotrefe[0] = -0.03493077177287424;
      rotrefe[1] = 0.3794316793480091;
      rotrefe[2] = -0.18138091085085048;
      rotrefe[3] = -0.03493077177287424;
      rotrefe[4] = 0.3794316793480091;
      rotrefe[5] = -0.18138091085085048;
      rotrefe[6] = -0.03493077177287424;
      rotrefe[7] = 0.3794316793480091;
      rotrefe[8] = -0.18138091085085048;

      testele_->set_centerline_hermite(true);
      testele_->set_up_reference_geometry<3, 2, 2>(xrefe_full, rotrefe);

      // create 1 element discretization
      testdis_->add_element(testele_);
      testdis_->fill_complete(true, false, false);
    }

   protected:
    //! dummy discretization for holding element and node pointers
    std::shared_ptr<Core::FE::Discretization> testdis_;
    //! the beam3r element to be tested
    std::shared_ptr<Discret::Elements::Beam3r> testele_;
  };

  /**
   * Test reference length calculation of Simo-Reissner beam
   */
  TEST_F(Beam3r, RefLength)
  {
    EXPECT_NEAR(testele_->ref_length(), 0.54772255750516607, testTolerance);
  }

  /**
   * Test nodal nullspace calculation of Simo-Reissner beam. The given matrix A multiplied by the
   * nullspace B has to be zero A*B=0.
   * The given beam element has three nodes for the interpolation of the rotations and two for the
   * interpolation of the centerline displacement. This results in two nodes with 9 DOFs
   * (3 displacements + 3 rotations + 3 tangents) and one with 3 DOFs (3 rotations).
   * The matrix size thus is 21.
   */
  TEST_F(Beam3r, ComputeNullSpace)
  {
    EXPECT_EQ(testele_->num_node(), 3);
    EXPECT_EQ(testele_->num_centerline_nodes(), 2);

    Core::LinAlg::SparseMatrix A(*testdis_->dof_row_map(), 21);
    {
      const std::array<int, 21> indices = {
          0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20};
      std::array<double, 21> values;

      values = {1.97165966578126, -0.179241787798296, -0.358483575596592, 0, 0, 0,
          0.179987079111915, -0.0163624617374468, -0.0327249234748937, -1.97165966578126,
          0.179241787798296, 0.358483575596592, 0, 0, 0, 0.179987079111915, -0.0163624617374468,
          -0.0327249234748937, 0, -0.143393430238637, 0.0716967151193184};
      A.insert_global_values(0, indices.size(), values.data(), indices.data());

      values = {-0.179241787798296, 1.11129908434944, 0.0716967151193185, 0, 0, 0,
          -0.0163624617374468, 0.101447262772170, 0.00654498469497874, 0.179241787798296,
          -1.11129908434944, -0.0716967151193185, 0, 0, 0, -0.0163624617374468, 0.101447262772170,
          0.00654498469497874, 0.143393430238637, 0, 0.358483575596592};
      A.insert_global_values(1, indices.size(), values.data(), indices.data());

      values = {-0.358483575596592, 0.0716967151193185, 1.21884415702841, 0, 0, 0,
          -0.0327249234748937, 0.00654498469497874, 0.111264739814639, 0.358483575596592,
          -0.0716967151193185, -1.21884415702841, 0, 0, 0, -0.0327249234748937, 0.00654498469497874,
          0.111264739814639, -0.0716967151193184, -0.358483575596592, 0};
      A.insert_global_values(2, indices.size(), values.data(), indices.data());

      values = {0, 0, 0, 0.00806588045092352, 0.00597472625994320, 0.0119494525198864, 0,
          -0.0130899693899575, 0.00654498469497874, 0, 0, 0, 0.000298736312997262, 0, 0, 0, 0, 0,
          -0.00238989050397757, 0, 0};
      A.insert_global_values(3, indices.size(), values.data(), indices.data());

      values = {0, 0, 0, 0.00597472625994321, 0.0367445664986509, -0.00238989050397728,
          0.0130899693899575, 0, 0.0327249234748937, 0, 0, 0, 0, 0.000298736312997262, 0, 0, 0, 0,
          0, -0.00238989050397757, 0};
      A.insert_global_values(4, indices.size(), values.data(), indices.data());

      values = {0, 0, 0, 0.0119494525198864, -0.00238989050397728, 0.0331597307426850,
          -0.00654498469497874, -0.0327249234748937, 0, 0, 0, 0, 0, 0, 0.000298736312997262, 0, 0,
          0, 0, 0, -0.00238989050397757};
      A.insert_global_values(5, indices.size(), values.data(), indices.data());

      values = {0.179987079111915, -0.0163624617374468, -0.0327249234748937, 0, 0.0130899693899575,
          -0.00654498469497873, 0.0821524860742190, -0.00746840782492901, -0.0149368156498580,
          -0.179987079111915, 0.0163624617374468, 0.0327249234748937, 0, 0, 0, 0.0164304972148438,
          -0.00149368156498580, -0.00298736312997160, 0, -0.0130899693899575, 0.00654498469497873};
      A.insert_global_values(6, indices.size(), values.data(), indices.data());

      values = {-0.0163624617374468, 0.101447262772170, 0.00654498469497874, -0.0130899693899575, 0,
          -0.0327249234748937, -0.00746840782492901, 0.0463041285145598, 0.00298736312997160,
          0.0163624617374468, -0.101447262772170, -0.00654498469497874, 0, 0, 0,
          -0.00149368156498580, 0.00926082570291197, 0.000597472625994321, 0.0130899693899575, 0,
          0.0327249234748937};
      A.insert_global_values(7, indices.size(), values.data(), indices.data());

      values = {-0.0327249234748937, 0.00654498469497874, 0.111264739814639, 0.00654498469497874,
          0.0327249234748937, 0, -0.0149368156498580, 0.00298736312997160, 0.0507851732095172,
          0.0327249234748937, -0.00654498469497874, -0.111264739814639, 0, 0, 0,
          -0.00298736312997160, 0.000597472625994321, 0.0101570346419034, -0.00654498469497874,
          -0.0327249234748937, 0};
      A.insert_global_values(8, indices.size(), values.data(), indices.data());

      values = {-1.97165966578126, 0.179241787798296, 0.358483575596592, 0, 0, 0,
          -0.179987079111915, 0.0163624617374468, 0.0327249234748937, 1.97165966578126,
          -0.179241787798296, -0.358483575596592, 0, 0, 0, -0.179987079111915, 0.0163624617374468,
          0.0327249234748937, 0, 0.143393430238637, -0.0716967151193184};
      A.insert_global_values(9, indices.size(), values.data(), indices.data());

      values = {0.179241787798296, -1.11129908434944, -0.0716967151193185, 0, 0, 0,
          0.0163624617374468, -0.101447262772170, -0.00654498469497874, -0.179241787798296,
          1.11129908434944, 0.0716967151193185, 0, 0, 0, 0.0163624617374468, -0.101447262772170,
          -0.00654498469497874, -0.143393430238637, 0, -0.358483575596592};
      A.insert_global_values(10, indices.size(), values.data(), indices.data());

      values = {0.358483575596592, -0.0716967151193185, -1.21884415702841, 0, 0, 0,
          0.0327249234748937, -0.00654498469497874, -0.111264739814639, -0.358483575596592,
          0.0716967151193185, 1.21884415702841, 0, 0, 0, 0.0327249234748937, -0.00654498469497874,
          -0.111264739814639, 0.0716967151193184, 0.358483575596592, 0};
      A.insert_global_values(11, indices.size(), values.data(), indices.data());

      values = {0, 0, 0, 0.000298736312997262, 0, 0, 0, 0, 0, 0, 0, 0, 0.00806588045092352,
          0.00597472625994320, 0.0119494525198864, 0, -0.0130899693899575, 0.00654498469497874,
          -0.00238989050397757, 0, 0};
      A.insert_global_values(12, indices.size(), values.data(), indices.data());

      values = {0, 0, 0, 0, 0.000298736312997262, 0, 0, 0, 0, 0, 0, 0, 0.00597472625994321,
          0.0367445664986509, -0.00238989050397728, 0.0130899693899575, 0, 0.0327249234748937, 0,
          -0.00238989050397757, 0};
      A.insert_global_values(13, indices.size(), values.data(), indices.data());

      values = {0, 0, 0, 0, 0, 0.000298736312997262, 0, 0, 0, 0, 0, 0, 0.0119494525198864,
          -0.00238989050397728, 0.0331597307426850, -0.00654498469497874, -0.0327249234748937, 0, 0,
          0, -0.00238989050397757};
      A.insert_global_values(14, indices.size(), values.data(), indices.data());

      values = {0.179987079111915, -0.0163624617374468, -0.0327249234748937, 0, 0, 0,
          0.0164304972148438, -0.00149368156498580, -0.00298736312997160, -0.179987079111915,
          0.0163624617374468, 0.0327249234748937, 0, 0.0130899693899575, -0.00654498469497873,
          0.0821524860742190, -0.00746840782492901, -0.0149368156498580, 0, -0.0130899693899575,
          0.00654498469497873};
      A.insert_global_values(15, indices.size(), values.data(), indices.data());

      values = {-0.0163624617374468, 0.101447262772170, 0.00654498469497874, 0, 0, 0,
          -0.00149368156498580, 0.00926082570291197, 0.000597472625994321, 0.0163624617374468,
          -0.101447262772170, -0.00654498469497874, -0.0130899693899575, 0, -0.0327249234748937,
          -0.00746840782492901, 0.0463041285145598, 0.00298736312997160, 0.0130899693899575, 0,
          0.0327249234748937};
      A.insert_global_values(16, indices.size(), values.data(), indices.data());

      values = {-0.0327249234748937, 0.00654498469497874, 0.111264739814639, 0, 0, 0,
          -0.00298736312997160, 0.000597472625994321, 0.0101570346419034, 0.0327249234748937,
          -0.00654498469497874, -0.111264739814639, 0.00654498469497874, 0.0327249234748937, 0,
          -0.0149368156498580, 0.00298736312997160, 0.0507851732095172, -0.00654498469497874,
          -0.0327249234748937, 0};
      A.insert_global_values(17, indices.size(), values.data(), indices.data());

      values = {0, 0.143393430238637, -0.0716967151193184, -0.00238989050397757, 0, 0, 0,
          0.0130899693899575, -0.00654498469497874, 0, -0.143393430238637, 0.0716967151193184,
          -0.00238989050397757, 0, 0, 0, 0.0130899693899575, -0.00654498469497874,
          0.0286786860477280, 0.0238989050397728, 0.0477978100795456};
      A.insert_global_values(18, indices.size(), values.data(), indices.data());

      values = {-0.143393430238637, 0, -0.358483575596592, 0, -0.00238989050397757, 0,
          -0.0130899693899575, 0, -0.0327249234748937, 0.143393430238637, 0, 0.358483575596592, 0,
          -0.00238989050397757, 0, -0.0130899693899575, 0, -0.0327249234748937, 0.0238989050397728,
          0.143393430238638, -0.00955956201590912};
      A.insert_global_values(19, indices.size(), values.data(), indices.data());

      values = {0.0716967151193184, 0.358483575596592, 0, 0, 0, -0.00238989050397757,
          0.00654498469497873, 0.0327249234748937, 0, -0.0716967151193184, -0.358483575596592, 0, 0,
          0, -0.00238989050397757, 0.00654498469497873, 0.0327249234748937, 0, 0.0477978100795457,
          -0.00955956201590912, 0.129054087214774};
      A.insert_global_values(20, indices.size(), values.data(), indices.data());
    }
    A.complete();

    int numdof, dimnsp;
    testele_->element_type().nodal_block_information(testele_.get(), numdof, dimnsp);

    const auto B =
        Core::FE::compute_null_space(*testdis_, numdof, dimnsp, *testdis_->dof_row_map());
    Core::LinAlg::MultiVector<double> zero(*testdis_->dof_row_map(), 6);

    A.multiply(false, *B, zero);

    // Check the inf-norm of the nullspace vectors
    // norm is initialized with ones just to avoid having accidental zeros
    std::array<double, 6> inf_norm{1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
    zero.NormInf(inf_norm.data());

    EXPECT_NEAR(inf_norm[0], 0.0, 1e-12);
    EXPECT_NEAR(inf_norm[1], 0.0, 1e-12);
    EXPECT_NEAR(inf_norm[2], 0.0, 1e-12);
    EXPECT_NEAR(inf_norm[3], 0.0, 1e-12);
    EXPECT_NEAR(inf_norm[4], 0.0, 1e-12);
    EXPECT_NEAR(inf_norm[5], 0.0, 1e-12);
  }

}  // namespace
