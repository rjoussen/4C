// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_POROELAST_MONOLITHICSPLIT_NOPENETRATION_HPP
#define FOUR_C_POROELAST_MONOLITHICSPLIT_NOPENETRATION_HPP

#include "4C_config.hpp"

#include "4C_coupling_adapter_converter.hpp"
#include "4C_poroelast_monolithicsplit.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Adapter
{
  class CouplingNonLinMortar;
}

namespace PoroElast
{
  //! monolithic structure split for condensing DOFs, when using the brinkman-equation
  class MonolithicSplitNoPenetration : public MonolithicSplit
  {
   public:
    explicit MonolithicSplitNoPenetration(MPI_Comm comm, const Teuchos::ParameterList& timeparams,
        std::shared_ptr<Core::LinAlg::MapExtractor> porosity_splitter);

    //! Setup the monolithic system
    void setup_system() override;

    //! setup composed right hand side from field solvers
    void setup_rhs(bool firstcall = false) override;

    //! setup composed system matrix from field solvers
    void setup_system_matrix(Core::LinAlg::BlockSparseMatrixBase& mat) override;

    //! take current results for converged and save for next time step
    void update() override;

    //! read restart data
    void read_restart(const int step) override;

    //! contains text to print_newton_iter
    void print_newton_iter_text_stream(std::ostringstream& oss) override;

    //! contains header to print_newton_iter
    void print_newton_iter_header_stream(std::ostringstream& oss) override;


   protected:
    //! @name Apply current field state to system

    //! Evaluate mechanical-fluid system matrix
    void apply_str_coupl_matrix(
        std::shared_ptr<Core::LinAlg::SparseOperator> k_sf  //!< mechanical-fluid stiffness matrix
        ) override;

    //! Evaluate fluid-mechanical system matrix
    void apply_fluid_coupl_matrix(
        std::shared_ptr<Core::LinAlg::SparseOperator> k_fs  //!< fluid-mechanical tangent matrix
        ) override;

    //! recover Lagrange multiplier \f$\lambda_\Gamma\f$ at the interface at the end of each
    //! iteration step (i.e. condensed forces onto the structure) needed for rhs in next newton step
    void recover_lagrange_multiplier_after_newton_step(
        std::shared_ptr<const Core::LinAlg::Vector<double>> x) override;

    //! recover Lagrange multiplier \f$\lambda_\Gamma\f$ at the interface at the end of each time
    //! step (i.e. condensed forces onto the structure) needed for rhs in next time step
    void recover_lagrange_multiplier_after_time_step() override;

    //! output
    void output(bool forced_writerestart = false) override;

    //! setup of coupling object and system matrices
    void setup_coupling_and_matrices() override;

    //! start a new time step
    void prepare_time_step() override;

    //! convergence check for Newton solver
    void build_convergence_norms() override;

    //!@}
   private:
    //! build block vector from field vectors, e.g. rhs, increment vector
    void setup_vector(Core::LinAlg::Vector<double>& f,  //!< vector of length of all dofs
        std::shared_ptr<const Core::LinAlg::Vector<double>>
            sv,  //!< vector containing only structural dofs
        std::shared_ptr<const Core::LinAlg::Vector<double>>
            fv  //!< vector containing only fluid dofs
        ) override;

    //! @name Global matrices and vectors

    std::shared_ptr<Core::LinAlg::SparseOperator> k_struct_;
    std::shared_ptr<Core::LinAlg::SparseOperator> k_fluid_;

    std::shared_ptr<Core::LinAlg::SparseMatrix> k_lambda_;
    std::shared_ptr<Core::LinAlg::SparseMatrix> k_d_;
    std::shared_ptr<Core::LinAlg::SparseMatrix> k_inv_d_;

    std::shared_ptr<Core::LinAlg::SparseMatrix> k_dn_;

    std::shared_ptr<Core::LinAlg::SparseMatrix> k_lambdainv_d_;

    std::shared_ptr<Core::LinAlg::BlockSparseMatrixBase> k_porodisp_;
    std::shared_ptr<Core::LinAlg::SparseOperator> k_porofluid_;

    std::shared_ptr<Core::LinAlg::Vector<double>> nopenetration_rhs_;

    //! transform object for k_D matrix \f$D\f$
    std::shared_ptr<Coupling::Adapter::MatrixColTransform> k_d_transform_;
    //! transform object for k_D matrix \f$D\f$
    std::shared_ptr<Coupling::Adapter::MatrixRowTransform> k_inv_d_transform_;

    //! transform object for linearization of k_D matrix \f$D\f$
    std::shared_ptr<Coupling::Adapter::MatrixColTransform> k_d_lin_transform_;

    //! Lagrange multiplier \f$\lambda_\Gamma^{n+1}\f$ at the interface (ie condensed forces onto
    //! the structure) evaluated at actual iteration step \f$t_{n+1}\f$ but needed for next
    //! iteration step
    std::shared_ptr<Core::LinAlg::Vector<double>> lambdanp_;

    //!@}

    //! @name Some quantities to recover the Lagrange multiplier at the end of each iteration step

    //! block \f$F_{\Gamma I,i+1}\f$ of fluid matrix at current iteration \f$i+1\f$
    std::shared_ptr<const Core::LinAlg::SparseMatrix> fgicur_;

    //! block \f$S_{\Gamma\Gamma,i+1}\f$ of fluid matrix at current iteration \f$i+1\f$
    std::shared_ptr<const Core::LinAlg::SparseMatrix> fggcur_;

    //! block \f$Cfs_{\Gamma\Gamma,i+1}\f$ of fs-coupling matrix at current iteration \f$i+1\f$
    std::shared_ptr<const Core::LinAlg::SparseMatrix> cfsgicur_;

    //! block \f$Cfs_{\Gamma\Gamma,i+1}\f$ of fs-coupling matrix at current iteration \f$i+1\f$
    std::shared_ptr<const Core::LinAlg::SparseMatrix> cfsggcur_;

    std::shared_ptr<Core::LinAlg::Vector<double>> rhs_fgcur_;

    //!@}

    //! @name Iterative solution technique

    //!< norm of no-penetration constraint
    double normrhs_nopenetration_;
    //!@}

    std::shared_ptr<Adapter::CouplingNonLinMortar> mortar_adapter_;
  };

}  // namespace PoroElast


FOUR_C_NAMESPACE_CLOSE

#endif
