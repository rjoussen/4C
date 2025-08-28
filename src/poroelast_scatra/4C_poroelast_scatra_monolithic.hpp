// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_POROELAST_SCATRA_MONOLITHIC_HPP
#define FOUR_C_POROELAST_SCATRA_MONOLITHIC_HPP

/*----------------------------------------------------------------------*
 | header inclusions                                                     |
 *----------------------------------------------------------------------*/
#include "4C_config.hpp"

#include "4C_linalg_mapextractor.hpp"
#include "4C_poroelast_scatra_base.hpp"
#include "4C_poroelast_scatra_input.hpp"
#include "4C_utils_parameter_list.fwd.hpp"

#include <Teuchos_Time.hpp>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 |  forward declarations                                              |
 *----------------------------------------------------------------------*/
namespace Core::LinAlg
{
  //  class SparseMatrix;
  //  class SparseOperator;
  //
  //  class BlockSparseMatrixBase;
  //  class Solver;
  class MapExtractor;
  class MultiMapExtractor;
}  // namespace Core::LinAlg

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/

namespace PoroElastScaTra
{
  /// base class of all monolithic porous media - scalar transport - interaction algorithms
  class PoroScatraMono : public PoroScatraBase
  {
   public:
    explicit PoroScatraMono(MPI_Comm comm,
        const Teuchos::ParameterList& timeparams);  // Problem builder

    //! Main time loop.
    void timeloop() override;

    //! read and set fields needed for restart
    void read_restart(int restart) override;

    //! prepare time step
    void prepare_time_step(bool printheader = true) override;

    //! is convergence reached of iterative solution technique?
    bool converged();

    /*! do the setup for the monolithic system


     1.) setup coupling
     2.) get maps for all blocks in the system (and for the whole system as well)
     create combined map
     3.) create system matrix


     \note We want to do this setup after reading the restart information, not
     directly in the constructor. This is necessary since during restart (if
     read_mesh is called), the dofmaps for the blocks might get invalid.
     */
    //! Setup the monolithic Poroelasticity system
    void setup_system() override;

    //! setup composed right hand side from field solvers
    virtual void setup_rhs(bool firstcall = false);

    /// setup composed system matrix from field solvers
    virtual void setup_system_matrix();

    //! evaluate all fields at x^n+1 with x^n+1 = x_n + stepinc
    virtual void evaluate(std::shared_ptr<const Core::LinAlg::Vector<double>>
            stepinc  //!< increment between time step n and n+1
    );

    //! solve one time step
    void solve() override;

    //! take current results for converged and save for next time step
    void update() override;

    //! write output
    void output() override;

    // Setup solver for monolithic system
    bool setup_solver() override;

    //! @name Access methods

    //! composed system matrix
    std::shared_ptr<Core::LinAlg::SparseMatrix> system_matrix();

    //! right hand side vector
    std::shared_ptr<Core::LinAlg::Vector<double>> rhs() { return rhs_; };

    //! full monolithic dof row map
    std::shared_ptr<const Core::LinAlg::Map> dof_row_map() const;

    //! unique map of all dofs that should be constrained with DBC
    std::shared_ptr<const Core::LinAlg::Map> combined_dbc_map() const;

    //@}

   protected:
    //! extractor to communicate between full monolithic map and block maps
    std::shared_ptr<const Core::LinAlg::MultiMapExtractor> extractor() const
    {
      return blockrowdofmap_;
    }

    //! extractor for DBCs
    const std::shared_ptr<Core::LinAlg::MapExtractor>& dbc_extractor() const { return dbcmaps_; }

    //! set full monolithic dof row map
    /*!
     A subclass calls this method (from its constructor) and thereby
     defines the number of blocks, their maps and the block order. The block
     maps must be row maps by themselves and must not contain identical GIDs.
     */
    void set_dof_row_maps(const std::vector<std::shared_ptr<const Core::LinAlg::Map>>& maps);

    //! @name Apply current field state to system

    //! Evaluate off diagonal matrix in poro row
    void evaluate_od_block_mat_poro();

    //! Evaluate off diagonal matrix in scatra row
    void evaluate_od_block_mat_scatra();

    //@}


   private:
    //! build block vector from field vectors, e.g. rhs, increment vector
    void setup_vector(Core::LinAlg::Vector<double>& f,  //!< vector of length of all dofs
        std::shared_ptr<const Core::LinAlg::Vector<double>>
            pv,  //!< vector containing only structural dofs
        std::shared_ptr<const Core::LinAlg::Vector<double>>
            sv  //!< vector containing only fluid dofs
    );

    //! perform one time step (setup + solve + output)
    void do_time_step();

    //! calculate stress, strains, energies ...
    void prepare_output() override;

    //! @name helper methods for Newton loop

    void build_convergence_norms();

    //! solve linear system
    void linear_solve();

    //@}

    //! @name Newton Output

    //! print to screen information about residual forces and displacements
    void print_newton_iter();

    //! contains text to print_newton_iter
    void print_newton_iter_text(FILE* ofile  //!< output file handle
    );

    //! contains header to print_newton_iter
    void print_newton_iter_header(FILE* ofile  //!< output file handle
    );

    //! print statistics of converged Newton-Raphson iteration
    void print_newton_conv();

    //@}

    void fd_check();

    //! @name Printing and output
    //@{

    int printscreen_;  //!< print infos to standard out every printscreen_ steps
    bool printiter_;   //!< print intermediate iterations during solution

    //@}

    //! @name General purpose algorithm members
    //@{
    std::shared_ptr<Core::LinAlg::Solver> solver_;  //!< linear algebraic solver
    double solveradaptolbetter_;                    //!< tolerance to which is adapted ?
    bool solveradapttol_;                           //!< adapt solver tolerance
    //@}

    //! @name Iterative solution technique

    enum PoroElast::ConvNorm normtypeinc_;   //!< convergence check for increments
    enum PoroElast::ConvNorm normtypefres_;  //!< convergence check for residual forces
    enum PoroElast::BinaryOp combincfres_;  //!< binary operator to combine increments and residuals
    enum PoroElast::VectorNorm vectornormfres_;  //!< type of norm for residual
    enum PoroElast::VectorNorm vectornorminc_;   //!< type of norm for increments

    double tolinc_;   //!< tolerance residual increment
    double tolfres_;  //!< tolerance force residual

    double tolinc_struct_;   //!< tolerance residual increment for structure displacements
    double tolfres_struct_;  //!< tolerance force residual for structure displacements

    double tolinc_velocity_;   //!< tolerance residual increment for fluid velocity field
    double tolfres_velocity_;  //!< tolerance force residual for fluid velocity field

    double tolinc_pressure_;   //!< tolerance residual increment for fluid pressure field
    double tolfres_pressure_;  //!< tolerance force residual for fluid pressure field

    //    double tolinc_porosity_;     //!< tolerance residual increment for porosity field
    //    double tolfres_porosity_;    //!< tolerance force residual for porosity field

    double tolinc_scalar_;   //!< tolerance residual increment for scalar field
    double tolfres_scalar_;  //!< tolerance force residual for scalar field

    int itermax_;  //!< maximally permitted iterations
    int itermin_;  //!< minimally requested iteration

    //! norm of residual forces
    double normrhs_{std::numeric_limits<double>::infinity()};
    //! norm of residual unknowns
    double norminc_{std::numeric_limits<double>::infinity()};

    //! norm of residual forces (fluid velocity)
    double normrhsfluidvel_{std::numeric_limits<double>::infinity()};
    //! norm of residual unknowns (fluid velocity)
    double normincfluidvel_{std::numeric_limits<double>::infinity()};
    //! norm of residual forces (fluid pressure)
    double normrhsfluidpres_{std::numeric_limits<double>::infinity()};
    //! norm of residual unknowns (fluid pressure)
    double normincfluidpres_{std::numeric_limits<double>::infinity()};
    //! norm of residual forces (fluid )
    double normrhsfluid_{std::numeric_limits<double>::infinity()};
    //! norm of residual unknowns (fluid )
    double normincfluid_{std::numeric_limits<double>::infinity()};
    //!< norm of residual forces (structure)
    double normrhsstruct_{std::numeric_limits<double>::infinity()};
    //!< norm of residual unknowns (structure)
    double normincstruct_{std::numeric_limits<double>::infinity()};
    //!< norm of residual forces (scatra)
    double normrhsscalar_{std::numeric_limits<double>::infinity()};
    //!< norm of residual unknowns (scatra)
    double normincscalar_{std::numeric_limits<double>::infinity()};

    Teuchos::Time timer_;  //!< timer for solution technique

    int iter_;  //!< iteration step

    std::shared_ptr<Core::LinAlg::Vector<double>>
        iterinc_;  //!< increment between Newton steps k and k+1

    std::shared_ptr<Core::LinAlg::Vector<double>> zeros_;  //!< a zero vector of full length

    //@}

    //! @name variables of monolithic system

    //! block systemmatrix
    std::shared_ptr<Core::LinAlg::BlockSparseMatrixBase> systemmatrix_;

    //! rhs of monolithic system
    std::shared_ptr<Core::LinAlg::Vector<double>> rhs_;

    //! structure-scatra coupling matrix
    std::shared_ptr<Core::LinAlg::SparseMatrix> k_pss_;
    //! fluid-scatra coupling matrix
    std::shared_ptr<Core::LinAlg::SparseMatrix> k_pfs_;

    //! scatra-structure coupling matrix
    std::shared_ptr<Core::LinAlg::SparseMatrix> k_sps_;
    //! scatra-fluid coupling matrix
    std::shared_ptr<Core::LinAlg::SparseMatrix> k_spf_;

    //! dof row map split in (field) blocks
    std::shared_ptr<Core::LinAlg::MultiMapExtractor> blockrowdofmap_;

    //! scatra row map as map extractor (used to build coupling matrixes)
    Core::LinAlg::MultiMapExtractor scatrarowdofmap_;
    Core::LinAlg::MultiMapExtractor pororowdofmap_;

    //! dirichlet map of monolithic system
    std::shared_ptr<Core::LinAlg::MapExtractor> dbcmaps_;
    //@}

    //! flag for direct solver
    bool directsolve_;

  };  // class PoroScatraMono

}  // namespace PoroElastScaTra

FOUR_C_NAMESPACE_CLOSE

#endif
