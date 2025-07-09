// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <gtest/gtest.h>

#include "4C_linalg_fevector.hpp"

#include "4C_comm_mpi_utils.hpp"
#include "4C_linalg_map.hpp"
#include "4C_linalg_multi_vector.hpp"
#include "4C_linalg_sparsematrix.hpp"
#include "4C_linalg_vector.hpp"

#include <memory>

FOUR_C_NAMESPACE_OPEN

namespace
{
  class FEVectorTest : public testing::Test
  {
   public:
    MPI_Comm comm_;
    std::shared_ptr<Core::LinAlg::Map> map_;
    int NumGlobalElements = 10;

   protected:
    FEVectorTest()
    {
      // set up communicator
      comm_ = MPI_COMM_WORLD;

      // set up a map
      map_ = std::make_shared<Core::LinAlg::Map>(NumGlobalElements, 0, comm_);
    }
  };

  TEST_F(FEVectorTest, DeepCopying)
  {
    Core::LinAlg::FEVector<double> a(*map_, true);
    a.put_scalar(1.0);

    Core::LinAlg::FEVector<double> b(*map_, true);
    // copy assign
    b = a;
    b.put_scalar(2.0);
    double norm_a = 0.0;
    double norm_b = 0.0;
    a.norm_2(&norm_a);
    b.norm_2(&norm_b);

    EXPECT_FLOAT_EQ(norm_a, 1.0 * std::sqrt(NumGlobalElements));
    EXPECT_FLOAT_EQ(norm_b, 2.0 * std::sqrt(NumGlobalElements));

    // copy constructor
    Core::LinAlg::FEVector<double> c(a);
    c.put_scalar(3.0);
    double norm_c = 0.0;
    c.norm_2(&norm_c);
    EXPECT_FLOAT_EQ(norm_c, 3.0 * std::sqrt(NumGlobalElements));
  }

  TEST_F(FEVectorTest, PutScalar)
  {
    // initialize with false value
    double norm_of_test_vector = 0.0;

    // copy zero vector into new interface
    Core::LinAlg::FEVector<double> test_vector = Core::LinAlg::FEVector<double>(*map_, true);

    test_vector.put_scalar(2.0);

    // check result
    test_vector.norm_2(&norm_of_test_vector);
    ASSERT_FLOAT_EQ(NumGlobalElements * 2.0 * 2.0, norm_of_test_vector * norm_of_test_vector);
  }

  TEST_F(FEVectorTest, Update)
  {
    Core::LinAlg::FEVector<double> a = Core::LinAlg::FEVector<double>(*map_, true);
    a.put_scalar(1.0);

    Core::LinAlg::FEVector<double> b = Core::LinAlg::FEVector<double>(*map_, true);
    b.put_scalar(1.0);

    // update the vector
    b.update(2.0, a, 3.0);

    // initialize with false value
    double b_norm = 0.0;

    // check norm of vector
    b.norm_2(&b_norm);
    ASSERT_FLOAT_EQ(NumGlobalElements * (2.0 + 3.0) * (2.0 + 3.0), b_norm * b_norm);

    Core::LinAlg::FEVector<double> c = Core::LinAlg::FEVector<double>(*map_, true);
    c.update(1, a, -1, b, 0);

    // initialize with false value
    double c_norm = 0.0;

    // check norm of vector
    c.norm_1(&c_norm);
    ASSERT_FLOAT_EQ(4 * NumGlobalElements, c_norm);
  }


  TEST_F(FEVectorTest, View)
  {
    Epetra_FEVector a(map_->get_epetra_block_map(), true);
    a.PutScalar(1.0);
    // Scope in which a is modified by the view
    {
      Core::LinAlg::View a_view(a);

      double norm = 0.0;
      ((Core::LinAlg::FEVector<double>&)a_view).norm_2(&norm);
      EXPECT_EQ(norm, std::sqrt(NumGlobalElements));

      ((Core::LinAlg::FEVector<double>&)a_view).put_scalar(2.0);
    }
    const Epetra_FEVector& a_const = a;
    Core::LinAlg::View a_view_const(a_const);
    // Change must be reflected in a
    double norm = 0.0;
    static_cast<const Core::LinAlg::FEVector<double>&>(a_view_const).norm_2(&norm);
    EXPECT_EQ(norm, 2.0 * std::sqrt(NumGlobalElements));
  }

  std::vector<double> means_multi_vector(const Core::LinAlg::MultiVector<double>& mv)
  {
    std::vector<double> means(mv.NumVectors());
    mv.MeanValue(means.data());
    return means;
  }


  TEST_F(FEVectorTest, MultiVectorImplicitConversionView)
  {
    Core::LinAlg::FEVector<double> a(*map_, true);
    a.put_scalar(1.0);

    // This views the data that is in a. It does not copy the data.
    // This results in the same behavior as inheritance would give.
    EXPECT_EQ(means_multi_vector(a)[0], 1.0);

    // This copies the data.
    Core::LinAlg::MultiVector<double> mv = a;
    a.put_scalar(2.0);

    // mv should still be 1.0 because we only modified a.
    EXPECT_EQ(means_multi_vector(mv)[0], 1.0);
  }

  TEST_F(FEVectorTest, MultiVectorImplicitConversionCopy)
  {
    auto a = std::make_shared<Core::LinAlg::FEVector<double>>(*map_, true);
    a->put_scalar(1.0);

    // This copies the data.
    Core::LinAlg::MultiVector<double> mv = *a;
    a->put_scalar(2.0);
    // Explicitly deallocate a to make sure that mv is not a view.
    a = nullptr;

    // mv should still be 1.0 because we only modified a.
    EXPECT_EQ(means_multi_vector(mv)[0], 1.0);
  }

  TEST_F(FEVectorTest, MultiVectorImplicitConversionRef)
  {
    Core::LinAlg::FEVector<double> a(*map_, true);
    a.put_scalar(1.0);

    Core::LinAlg::MultiVector<double>& mv = a;
    mv.PutScalar(2.0);
    EXPECT_EQ(means_multi_vector(a)[0], 2.0);

    // Reassigning to a must keep mv valid: move assign
    a = Core::LinAlg::FEVector<double>(*map_, true);
    EXPECT_EQ(means_multi_vector(mv)[0], 0.0);
    a.put_scalar(3.0);
    EXPECT_EQ(means_multi_vector(mv)[0], 3.0);

    // Reassigning to a must keep mv valid: copy assign
    Core::LinAlg::FEVector<double> b(*map_, true);
    a = b;
    EXPECT_EQ(means_multi_vector(mv)[0], 0.0);
    a.put_scalar(4.0);
    EXPECT_EQ(means_multi_vector(mv)[0], 4.0);
  }

  TEST_F(FEVectorTest, AssignToRef)
  {
    Core::LinAlg::FEVector<double> a(*map_, true);
    a.put_scalar(1.0);
    EXPECT_EQ(means_multi_vector(a)[0], 1.0);
    Core::LinAlg::MultiVector<double>& mv = a;
    // Actually assign an MV to a via the ref. Note that this would throw in Trilinos if not using a
    // single column.
    mv = Core::LinAlg::MultiVector<double>(*map_, 1, true);
    EXPECT_EQ(means_multi_vector(mv)[0], 0.0);
  }


}  // namespace

FOUR_C_NAMESPACE_CLOSE
