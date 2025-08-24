// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_linear_solver_method_direct.hpp"

#include "4C_linalg_krylov_projector.hpp"
#include "4C_linalg_utils_sparse_algebra_math.hpp"

#if FOUR_C_TRILINOS_INTERNAL_VERSION_GE(2025, 3)
#else
#include <Amesos_Klu.h>
#include <Amesos_Superludist.h>
#include <Amesos_Umfpack.h>
#include <Epetra_LinearProblem.h>
#endif

FOUR_C_NAMESPACE_OPEN

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
Core::LinearSolver::DirectSolver::DirectSolver(std::string solvertype)
    : solvertype_(solvertype), factored_(false), solver_(nullptr), projector_(nullptr)
{
#if FOUR_C_TRILINOS_INTERNAL_VERSION_GE(2025, 3)
#else
  reindexer_ = nullptr;
  linear_problem_ = std::make_shared<Epetra_LinearProblem>();
#endif
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
void Core::LinearSolver::DirectSolver::setup(std::shared_ptr<Core::LinAlg::SparseOperator> matrix,
    std::shared_ptr<Core::LinAlg::MultiVector<double>> x,
    std::shared_ptr<Core::LinAlg::MultiVector<double>> b, const bool refactor, const bool reset,
    std::shared_ptr<Core::LinAlg::KrylovProjector> projector)
{
  auto crsA = std::dynamic_pointer_cast<Core::LinAlg::SparseMatrix>(matrix);

  // 1. merge the block system matrix into a standard sparse matrix if necessary
  if (!crsA)
  {
    std::shared_ptr<Core::LinAlg::BlockSparseMatrixBase> Ablock =
        std::dynamic_pointer_cast<Core::LinAlg::BlockSparseMatrixBase>(matrix);

    int matrixDim = Ablock->full_range_map().num_global_elements();
    if (matrixDim > 50000)
      std::cout << "\n WARNING: Direct linear solver is merging matrix, this is very expensive! \n";

    crsA = Ablock->merge();
  }

  // 2. project the linear system if close to being singular and set the final matrix and vectors
  projector_ = projector;
  if (projector_ != nullptr)
  {
    Core::LinAlg::SparseMatrix A_view(*crsA);
    crsA = projector_->project(A_view);

    projector_->apply_pt(*b);
  }

  x_ = x;
  b_ = b;
  a_ = crsA;

#if FOUR_C_TRILINOS_INTERNAL_VERSION_GE(2025, 3)
#else
  linear_problem_->SetRHS(&b_->get_epetra_multi_vector());
  linear_problem_->SetLHS(&x_->get_epetra_multi_vector());
  linear_problem_->SetOperator(a_->epetra_operator());

  if (reindexer_ and not(reset or refactor)) reindexer_->fwd();
#endif

  // 3. create linear solver
  if (reset or refactor or not is_factored())
  {
#if FOUR_C_TRILINOS_INTERNAL_VERSION_GE(2025, 3)
    std::string solver_type;
    Teuchos::ParameterList params("Amesos2");
#else
    reindexer_ = std::make_shared<EpetraExt::LinearProblem_Reindex2>(nullptr);
#endif

    if (solvertype_ == "umfpack")
    {
#if FOUR_C_TRILINOS_INTERNAL_VERSION_GE(2025, 3)
      solver_type = "Umfpack";
      auto umfpack_params = Teuchos::sublist(Teuchos::rcpFromRef(params), solver_type);
      umfpack_params->set("IsContiguous", false, "Are GIDs Contiguous");
#else
      solver_ = std::make_shared<Amesos_Umfpack>((*reindexer_)(*linear_problem_));
#endif
    }
    else if (solvertype_ == "superlu")
    {
#if FOUR_C_TRILINOS_INTERNAL_VERSION_GE(2025, 3)
      solver_type = "SuperLU_DIST";
      auto superludist_params = Teuchos::sublist(Teuchos::rcpFromRef(params), solver_type);
      superludist_params->set("Equil", true, "Whether to equilibrate the system before solve");
      superludist_params->set("RowPerm", "LargeDiag_MC64", "Row ordering");
      superludist_params->set("ReplaceTinyPivot", true, "Replace tiny pivot");
      superludist_params->set("IsContiguous", false, "Are GIDs Contiguous");
#else
      solver_ = std::make_shared<Amesos_Superludist>((*reindexer_)(*linear_problem_));
#endif
    }
    else
    {
#if FOUR_C_TRILINOS_INTERNAL_VERSION_GE(2025, 3)
      solver_type = "KLU2";
      auto klu_params = Teuchos::sublist(Teuchos::rcpFromRef(params), solver_type);
      klu_params->set("IsContiguous", false, "Are GIDs Contiguous");
#else
      solver_ = std::make_shared<Amesos_Klu>((*reindexer_)(*linear_problem_));
#endif
    }

#if FOUR_C_TRILINOS_INTERNAL_VERSION_GE(2025, 3)
    solver_ = Amesos2::create<Epetra_CrsMatrix, Epetra_MultiVector>(solver_type,
        Teuchos::rcp_dynamic_cast<Epetra_CrsMatrix>(Teuchos::rcpFromRef(*a_->epetra_operator())),
        Teuchos::rcpFromRef(x_->get_epetra_multi_vector()),
        Teuchos::rcpFromRef(b_->get_epetra_multi_vector()));

    solver_->setParameters(Teuchos::rcpFromRef(params));
#endif

    factored_ = false;
  }
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
int Core::LinearSolver::DirectSolver::solve()
{
  if (not is_factored())
  {
#if FOUR_C_TRILINOS_INTERNAL_VERSION_GE(2025, 3)
    solver_->symbolicFactorization();
    solver_->numericFactorization();
#else
    solver_->SymbolicFactorization();
    solver_->NumericFactorization();
#endif
    factored_ = true;
  }

#if FOUR_C_TRILINOS_INTERNAL_VERSION_GE(2025, 3)
  solver_->solve();
#else
  solver_->Solve();
#endif

  if (projector_ != nullptr) projector_->apply_p(*x_);

  return 0;
}

FOUR_C_NAMESPACE_CLOSE
