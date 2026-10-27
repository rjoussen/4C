// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_LINEAR_SOLVER_METHOD_DIRECT_HPP
#define FOUR_C_LINEAR_SOLVER_METHOD_DIRECT_HPP

#include "4C_config.hpp"

#include "4C_linalg_graph.hpp"
#include "4C_linalg_map.hpp"
#include "4C_linalg_sparseoperator.hpp"
#include "4C_linear_solver_method.hpp"
#include "4C_linear_solver_method_projector.hpp"

#include <Amesos2.hpp>

#include <optional>

FOUR_C_NAMESPACE_OPEN

namespace Core::LinearSolver
{
  /// direct linear solver (using amesos)
  class DirectSolver : public SolverTypeBase
  {
   public:
    /*! \brief Construct the solver
     *
     * @param solvertype Amesos2 solver to use
     * @param params Solver parameters. Amesos2 is told that the GIDs are contiguous if the system
     * is serial and its GIDs are exactly [0, N) in LID order. This detection is switched off if
     * params contains the bool "auto_detect_contiguous_dof_ids" set to false.
     *
     * For KLU2, UMFPACK and MUMPS the Amesos2 solver and its symbolic step are kept as long as the
     * sparsity pattern of the matrix is unchanged. For KLU2 and MUMPS this keeps the symbolic
     * factorization. Amesos2 releases the one of UMFPACK after every solve, so UMFPACK only saves
     * the reading of the matrix structure. This is switched off if params contains the bool
     * "reuse_symbolic_factorization" set to false.
     */
    explicit DirectSolver(
        Core::LinearSolver::SolverType solvertype, const Teuchos::ParameterList& params);

    /*! \brief Forget the system and the factorization
     *
     * If the symbolic factorization is reused, the Amesos2 solver is kept and setup() only keeps it
     * if the new system has the same graph and maps.
     */
    void reset() override;

    /*! \brief Setup the solver object
     *
     * @param matrix Matrix of the linear system
     * @param b Right-hand side vector of the linear system
     * @param refactor Boolean flag to enforce a refactorization of the matrix
     * @param projector Krylov projector
     */
    void setup(std::shared_ptr<Core::LinAlg::SparseOperator> matrix,
        std::shared_ptr<Core::LinAlg::MultiVector<double>> b, const bool refactor,
        std::shared_ptr<Core::LinAlg::LinearSystemProjector> projector = nullptr) override;

    int solve(Core::LinAlg::MultiVector<double>& x) override;

    [[nodiscard]] bool is_factored() const { return factored_; }

    //! Whether the symbolic factorization is kept across solves if the sparsity pattern is
    //! unchanged
    [[nodiscard]] bool reuses_symbolic_factorization() const { return reuse_symbolic_; }

    //! Number of symbolic factorizations performed so far
    [[nodiscard]] int num_symbolic_factorizations() const { return num_symbolic_factorizations_; }

   private:
    //! Whether Amesos2 can be told that the dof ids of the current system a_ and b_ are contiguous
    [[nodiscard]] bool has_contiguous_dof_ids() const;

    /*! \brief Whether the existing Amesos2 solver, including its symbolic factorization, is valid
     * for the current system a_ and b_
     *
     * This is the case if the matrix is the one the Amesos2 solver was created on and its
     * graph and maps are unchanged. The result is the same on all ranks.
     */
    [[nodiscard]] bool can_reuse_amesos2_solver() const;

    //! Store the graph and the right-hand side map of the current system a_ and b_ for later
    //! comparison
    void record_created_system();

    //! type/implementation of Amesos solver to be used
    const Core::LinearSolver::SolverType solvertype_;

    //! whether contiguous dof ids are detected automatically (opt-out)
    const bool auto_detect_contiguous_dof_ids_;

    //! whether the symbolic factorization is kept while the sparsity pattern is unchanged
    const bool reuse_symbolic_;

    //! flag indicating whether a valid factorization is stored
    bool factored_;

    //! flag indicating whether the Amesos2 solver has a symbolic factorization of the pattern
    bool symbolic_done_ = false;

    //! number of symbolic factorizations performed so far
    int num_symbolic_factorizations_ = 0;

    //! right hand side vector
    std::shared_ptr<Core::LinAlg::MultiVector<double>> b_;

    //! system of equations
    std::shared_ptr<Core::LinAlg::SparseMatrix> a_;

    //! an abstract Amesos2 solver that can be any of the concrete implementations
    Teuchos::RCP<Amesos2::Solver<Epetra_CrsMatrix, Epetra_MultiVector>> solver_;

    //! System the Amesos2 solver was created on
    struct CreatedSystem
    {
      //! identity of the Epetra matrix (non-owning)
      const Epetra_CrsMatrix* matrix_ = nullptr;

      //! identity of the Epetra right-hand side (non-owning)
      const Epetra_MultiVector* rhs_ = nullptr;

      std::optional<Core::LinAlg::Graph> graph;
      std::optional<Core::LinAlg::Map> rhs_map;
    };
    CreatedSystem created_system_;

    /*! \brief A projector applied before solving the linear systems
     *
     * Instead of solving Ax=b a projected system of the form P'APu=P'b is solved.
     * With P being the linear system projector.
     *
     * For example, for nearly singular systems, a Krylov projector can be used:
     *
     * P. Bochev and R. B. Lehoucq: On the Finite Element Solution of the Pure Neumann Problem,
     * SIAM Review, 47(1):50-66, 2005, http://dx.doi.org/10.1137/S0036144503426074
     *
     */
    std::shared_ptr<Core::LinAlg::LinearSystemProjector> projector_;
  };
}  // namespace Core::LinearSolver

FOUR_C_NAMESPACE_CLOSE

#endif
