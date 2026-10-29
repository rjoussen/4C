// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_linear_solver_method_direct.hpp"

#include "4C_comm_mpi_utils.hpp"
#include "4C_linalg_blocksparsematrix.hpp"
#include "4C_linalg_multi_vector.hpp"
#include "4C_linalg_sparseoperator.hpp"
#include "4C_linalg_vector.hpp"
#include "4C_linear_solver_method_projector.hpp"

#include <magic_enum/magic_enum.hpp>

#include <algorithm>
#include <memory>
#include <ranges>

FOUR_C_NAMESPACE_OPEN

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
Core::LinearSolver::DirectSolver::DirectSolver(
    Core::LinearSolver::SolverType solvertype, const Teuchos::ParameterList& params)
    : solvertype_(solvertype),
      auto_detect_contiguous_dof_ids_(!params.isParameter("auto_detect_contiguous_dof_ids") ||
                                      params.get<bool>("auto_detect_contiguous_dof_ids")),
      reuse_symbolic_((solvertype != SolverType::Superlu) &&
                      (!params.isParameter("reuse_symbolic_factorization") ||
                          params.get<bool>("reuse_symbolic_factorization"))),
      factored_(false),
      solver_(nullptr),
      projector_(nullptr)
{
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
bool Core::LinearSolver::DirectSolver::has_contiguous_dof_ids() const
{
  // Amesos2 'ROOTED' needs the GIDs to be exactly [base, base + N) in LID order: the matrix is
  // redistributed by GID while a single vector is copied in LID order. Epetra only flags maps as
  // linear if they are not built from a list of GIDs, so the GIDs are checked, too.
  const auto is_ordered_range = [](const Core::LinAlg::Map& map)
  {
    if (map.linear_map()) return true;

    return std::ranges::all_of(std::views::iota(0, map.num_my_elements()),
        [&map](int lid) { return map.gid(lid) == map.index_base() + lid; });
  };

  return a_->filled() && Core::Communication::num_mpi_ranks(a_->get_comm()) == 1 &&
         a_->num_global_rows() == a_->num_global_cols() &&
         a_->col_map().same_as(a_->domain_map()) && a_->row_map().same_as(a_->range_map()) &&
         b_->get_map().same_as(a_->range_map()) && is_ordered_range(a_->row_map()) &&
         is_ordered_range(a_->domain_map());
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
bool Core::LinearSolver::DirectSolver::can_reuse_amesos2_solver() const
{
  // These conditions are identical on all ranks.
  if (!reuse_symbolic_ || solver_.is_null() || projector_ != nullptr || !a_->filled()) return false;

  FOUR_C_ASSERT(system_snapshot_.graph.has_value(),
      "The system the Amesos2 solver was created on is unknown.");

  // Amesos2 holds on to the matrix it was created on.
  const bool same_matrix = system_snapshot_.matrix_ == &a_->epetra_matrix();

  // Comparing maps can be collective, so all comparisons are done unconditionally on every rank.
  const bool same_graph = system_snapshot_.graph->is_identical_to(a_->get_current_graph());
  const bool same_rhs_map = system_snapshot_.rhs_map->same_as(b_->get_map());

  const bool same = same_matrix && same_graph && same_rhs_map;

  // Factorizations are collective, so all ranks have to take the same decision.
  return Core::Communication::min_all(same ? 1 : 0, a_->get_comm()) == 1;
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
void Core::LinearSolver::DirectSolver::record_system_snapshot()
{
  system_snapshot_.matrix_ = &a_->epetra_matrix();
  system_snapshot_.rhs_ = &b_->get_epetra_multi_vector();
  system_snapshot_.graph.emplace(a_->get_current_graph());
  system_snapshot_.rhs_map = b_->get_map();
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
void Core::LinearSolver::DirectSolver::reset()
{
  a_ = nullptr;
  b_ = nullptr;
  projector_ = nullptr;
  factored_ = false;

  if (reuse_symbolic_) return;

  solver_ = Teuchos::null;
  system_snapshot_ = LinearSystemSnapshot{};
  symbolic_done_ = false;
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
void Core::LinearSolver::DirectSolver::setup(std::shared_ptr<Core::LinAlg::SparseOperator> matrix,
    std::shared_ptr<Core::LinAlg::MultiVector<double>> b, const bool refactor,
    std::shared_ptr<Core::LinAlg::LinearSystemProjector> projector)
{
  // 1. project the linear system if close to being singular and set the final matrix and vectors
  projector_ = projector;
  std::shared_ptr<Core::LinAlg::MultiVector<double>> rhs = b;
  if (projector_ != nullptr)
  {
    FOUR_C_ASSERT_ALWAYS(b->num_vectors() == 1,
        "Expecting only one solution vector during projector call! Got {} vectors.",
        b->num_vectors());
    LinAlg::Vector<double> projected_b = projector_->to_reduced(b->get_vector(0));
    rhs = std::make_shared<Core::LinAlg::MultiVector<double>>(projected_b.as_multi_vector());

    matrix = projector_->to_reduced(*matrix);
  }


  // 2. merge the block system matrix into a standard sparse matrix if necessary
  auto crsA = std::dynamic_pointer_cast<Core::LinAlg::SparseMatrix>(matrix);
  if (!crsA)
  {
    std::shared_ptr<Core::LinAlg::BlockSparseMatrixBase> Ablock =
        std::dynamic_pointer_cast<Core::LinAlg::BlockSparseMatrixBase>(matrix);

    int matrixDim = Ablock->full_range_map().num_global_elements();
    if (matrixDim > 50000 and Communication::my_mpi_rank(matrix->domain_map().get_comm()) == 0)
      std::cout << "\n WARNING: Direct linear solver is merging matrix, this is very expensive! \n";

    crsA = Ablock->merge();
  }

  b_ = rhs;
  a_ = crsA;

  // 3. create linear solver
  if (refactor or not is_factored() or system_snapshot_.matrix_ != &a_->epetra_matrix() or
      system_snapshot_.rhs_ != &b_->get_epetra_multi_vector())
  {
    if (can_reuse_amesos2_solver())
    {
      // Keep the Amesos2 solver with its symbolic factorization. Only the right-hand side vector
      // may differ, and the values are read again in the numeric factorization.
      solver_->setB(Teuchos::rcpFromRef(b_->get_epetra_multi_vector()));
      system_snapshot_.rhs_ = &b_->get_epetra_multi_vector();

      factored_ = false;
      return;
    }

    std::string solver_type;
    Teuchos::ParameterList params("Amesos2");

    const bool is_contiguous = auto_detect_contiguous_dof_ids_ && has_contiguous_dof_ids();

    switch (solvertype_)
    {
      case SolverType::KLU2:
      {
        solver_type = "KLU2";
        auto& klu_params = params.sublist(solver_type);
        klu_params.set("IsContiguous", is_contiguous, "Are GIDs Contiguous");
        break;
      }
      case SolverType::MUMPS:
      {
        solver_type = "MUMPS";
        auto& mumps_params = params.sublist(solver_type);
        mumps_params.set("IsContiguous", is_contiguous, "Are GIDs Contiguous");
        break;
      }
      case SolverType::UMFPACK:
      {
        solver_type = "Umfpack";
        auto& umfpack_params = params.sublist(solver_type);
        umfpack_params.set("IsContiguous", is_contiguous, "Are GIDs Contiguous");
        break;
      }
      case SolverType::Superlu:
      {
        solver_type = "SuperLU_DIST";
        auto& superludist_params = params.sublist(solver_type);
        superludist_params.set("Equil", true, "Whether to equilibrate the system before solve");
        superludist_params.set("RowPerm", "LargeDiag_MC64", "Row ordering");
        superludist_params.set("ReplaceTinyPivot", true, "Replace tiny pivot");
        superludist_params.set("IsContiguous", is_contiguous, "Are GIDs Contiguous");
        break;
      }
      default:
        FOUR_C_THROW("Unsupported solver type {}!", magic_enum::enum_name(solvertype_));
    }

    FOUR_C_ASSERT_ALWAYS(Amesos2::query(solver_type),
        "Requested direct solver {} is not available in Amesos2! Check your Amesos2 installation "
        "and choose an available solver!",
        solver_type);

    solver_ = Amesos2::create<Epetra_CrsMatrix, Epetra_MultiVector>(
        solver_type, Teuchos::rcpFromRef(a_->epetra_matrix()));
    solver_->setB(Teuchos::rcpFromRef(b_->get_epetra_multi_vector()));

    solver_->setParameters(Teuchos::make_rcp<Teuchos::ParameterList>(std::move(params)));

    record_system_snapshot();

    symbolic_done_ = false;
    factored_ = false;
  }
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
int Core::LinearSolver::DirectSolver::solve(Core::LinAlg::MultiVector<double>& x)
{
  std::unique_ptr<LinAlg::Vector<double>> projected_x = nullptr;
  if (projector_ != nullptr)
  {
    FOUR_C_ASSERT_ALWAYS(x.num_vectors() == 1,
        "Expecting only one solution vector during projector call! Got {} vectors.",
        x.num_vectors());

    // Create empty x in reduced space
    projected_x = std::make_unique<LinAlg::Vector<double>>(a_->domain_map(), true);

    solver_->setX(Teuchos::rcpFromRef(projected_x->as_multi_vector().get_epetra_multi_vector()));
  }
  else
  {
    solver_->setX(Teuchos::rcpFromRef(x.get_epetra_multi_vector()));
  }

  if (not is_factored())
  {
    if (not symbolic_done_)
    {
      solver_->symbolicFactorization();
      symbolic_done_ = true;
      ++num_symbolic_factorizations_;
    }
    solver_->numericFactorization();

    factored_ = true;
  }

  solver_->solve();

  if (projector_ != nullptr)
  {
    x.get_vector(0) = projector_->to_full(*projected_x);
  }

  return 0;
}

FOUR_C_NAMESPACE_CLOSE