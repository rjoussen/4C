// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <gtest/gtest.h>

#include "4C_linear_solver_method_direct.hpp"

#include "4C_linalg_map.hpp"
#include "4C_linalg_multi_vector.hpp"
#include "4C_linalg_sparsematrix.hpp"
#include "4C_linear_solver_method.hpp"

#include <mpi.h>
#include <Teuchos_ParameterList.hpp>

#include <memory>

namespace
{
  using namespace FourC;

  constexpr int num_rows = 20;

  // Tridiagonal matrix. The additional corner entry changes the sparsity pattern.
  void fill_matrix(Core::LinAlg::SparseMatrix& matrix, double diagonal, bool rebuild, bool corner)
  {
    if (rebuild)
      matrix.reset();
    else
      matrix.zero();

    for (int i = 0; i < num_rows; ++i)
    {
      matrix.set_value(diagonal, i, i);
      if (i > 0) matrix.set_value(-1.0, i, i - 1);
      if (i < num_rows - 1) matrix.set_value(-1.0, i, i + 1);
    }
    if (corner) matrix.set_value(-0.5, 0, num_rows - 1);

    if (!matrix.filled()) matrix.complete();
  }

  double relative_residual(const Core::LinAlg::SparseMatrix& matrix,
      const Core::LinAlg::MultiVector<double>& x, const Core::LinAlg::MultiVector<double>& b)
  {
    Core::LinAlg::MultiVector<double> residual(b.get_map(), 1, true);
    matrix.multiply(false, x, residual);
    residual.update(1.0, b, -1.0);

    double residual_norm = 0.0;
    double rhs_norm = 0.0;
    residual.norm_2(&residual_norm);
    b.norm_2(&rhs_norm);
    return residual_norm / rhs_norm;
  }

  class DirectSolverSymbolicReuse : public ::testing::TestWithParam<Core::LinearSolver::SolverType>
  {
   protected:
    void SetUp() override
    {
      if (GetParam() == Core::LinearSolver::SolverType::MUMPS && !Amesos2::query("MUMPS"))
        GTEST_SKIP() << "Amesos2 is built without MUMPS";
    }
    DirectSolverSymbolicReuse()
        : map_(num_rows, 0, MPI_COMM_WORLD),
          matrix_(std::make_shared<Core::LinAlg::SparseMatrix>(map_, 3)),
          b_(std::make_shared<Core::LinAlg::MultiVector<double>>(map_, 1, true)),
          x_(map_, 1, true)
    {
      for (int i = 0; i < num_rows; ++i) b_->replace_global_value(i, 0, 1.0 + i);
    }

    // Solve with the flags a nonlinear solver sets: refactor in every iteration, reset at the first
    void solve_and_check(Core::LinearSolver::DirectSolver& solver)
    {
      solver.reset();
      solver.setup(matrix_, b_, /*refactor=*/true);
      solver.solve(x_);
      EXPECT_LT(relative_residual(*matrix_, x_, *b_), 1.0e-12);
    }

    Core::LinAlg::Map map_;
    std::shared_ptr<Core::LinAlg::SparseMatrix> matrix_;
    std::shared_ptr<Core::LinAlg::MultiVector<double>> b_;
    Core::LinAlg::MultiVector<double> x_;
  };

  TEST_P(DirectSolverSymbolicReuse, ChangedValuesKeepSymbolicFactorization)
  {
    Core::LinearSolver::DirectSolver solver(GetParam(), Teuchos::ParameterList());
    ASSERT_TRUE(solver.reuses_symbolic_factorization());

    fill_matrix(*matrix_, 3.0, /*rebuild=*/true, /*corner=*/false);
    solve_and_check(solver);
    EXPECT_EQ(solver.num_symbolic_factorizations(), 1);

    // same matrix object and pattern, other values
    fill_matrix(*matrix_, 4.0, /*rebuild=*/false, /*corner=*/false);
    solve_and_check(solver);
    fill_matrix(*matrix_, 5.0, /*rebuild=*/false, /*corner=*/false);
    solve_and_check(solver);
    EXPECT_EQ(solver.num_symbolic_factorizations(), 1);
  }

  TEST_P(DirectSolverSymbolicReuse, ChangedPatternRedoesSymbolicFactorization)
  {
    Core::LinearSolver::DirectSolver solver(GetParam(), Teuchos::ParameterList());

    fill_matrix(*matrix_, 3.0, /*rebuild=*/true, /*corner=*/false);
    solve_and_check(solver);
    const int num_symbolic_before = solver.num_symbolic_factorizations();

    // The rebuilt matrix may live at the same address, but its pattern differs
    fill_matrix(*matrix_, 3.0, /*rebuild=*/true, /*corner=*/true);
    solve_and_check(solver);
    EXPECT_GT(solver.num_symbolic_factorizations(), num_symbolic_before);

    // and back to the original pattern
    fill_matrix(*matrix_, 4.0, /*rebuild=*/true, /*corner=*/false);
    solve_and_check(solver);
  }

  TEST_P(DirectSolverSymbolicReuse, OptOutRedoesSymbolicFactorization)
  {
    Teuchos::ParameterList params;
    params.set<bool>("reuse_symbolic_factorization", false);
    Core::LinearSolver::DirectSolver solver(GetParam(), params);
    ASSERT_FALSE(solver.reuses_symbolic_factorization());

    fill_matrix(*matrix_, 3.0, /*rebuild=*/true, /*corner=*/false);
    solve_and_check(solver);
    fill_matrix(*matrix_, 4.0, /*rebuild=*/false, /*corner=*/false);
    solve_and_check(solver);
    EXPECT_EQ(solver.num_symbolic_factorizations(), 2);
  }

  INSTANTIATE_TEST_SUITE_P(KLU2UMFPACKAndMUMPS, DirectSolverSymbolicReuse,
      ::testing::Values(Core::LinearSolver::SolverType::KLU2,
          Core::LinearSolver::SolverType::UMFPACK, Core::LinearSolver::SolverType::MUMPS));

  TEST(DirectSolverSymbolicReuseSolverTypes, SuperluDoesNotReuse)
  {
    Core::LinearSolver::DirectSolver solver(
        Core::LinearSolver::SolverType::Superlu, Teuchos::ParameterList());
    EXPECT_FALSE(solver.reuses_symbolic_factorization());
  }
}  // namespace
