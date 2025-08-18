// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_FSI_FLUIDFLUIDMONOLITHIC_STRUCTURESPLIT_NONOX_HPP
#define FOUR_C_FSI_FLUIDFLUIDMONOLITHIC_STRUCTURESPLIT_NONOX_HPP

#include "4C_config.hpp"

#include "4C_coupling_adapter_converter.hpp"
#include "4C_fsi_monolithic_nonox.hpp"

FOUR_C_NAMESPACE_OPEN

namespace FSI
{
  /// monolithic Fluid-Fluid FSI algorithm (structuresplit)
  /*!
    Here the structural matrix is split whereas the fluid matrix is taken as
    it is.

  */
  class FluidFluidMonolithicStructureSplitNoNOX : public MonolithicNoNOX
  {
    friend class FSI::FSIResultTest;

   public:
    /// constructor
    explicit FluidFluidMonolithicStructureSplitNoNOX(
        MPI_Comm comm, const Teuchos::ParameterList& timeparams);

    /*! do the setup for the monolithic system


    1.) setup coupling; right now, we use matching meshes at the interface
    2.) create combined map
    3.) create block system matrix


    */
    void setup_system() override;

   protected:
    //! @name Apply current field state to system

    /// setup composed right hand side from field solvers
    void setup_rhs(Core::LinAlg::Vector<double>& f, bool firstcall) override;

    /// setup composed system block matrix
    void setup_system_matrix() override;
    //@}

    /// create merged map of DOF in the final system from all fields
    void create_combined_dof_row_map() override;

    /// Extract initial guess from fields
    void initial_guess(std::shared_ptr<Core::LinAlg::Vector<double>> ig) override;

    /// apply infnorm scaling to linear block system
    virtual void scale_system(
        Core::LinAlg::BlockSparseMatrixBase& mat, Core::LinAlg::Vector<double>& b);

    /// undo infnorm scaling from scaled solution
    virtual void unscale_solution(Core::LinAlg::BlockSparseMatrixBase& mat,
        Core::LinAlg::Vector<double>& x, Core::LinAlg::Vector<double>& b);

    /// create merged map with Dirichlet-constrained DOF from all fields
    std::shared_ptr<Core::LinAlg::Map> combined_dbc_map() override;

    //! Extract the three field vectors from a given composed vector
    //!
    //! In analogy to NOX, x is step increment \f$\Delta x\f$
    //! that brings us from \f$t^{n}\f$ to \f$t^{n+1}\f$:
    //! \f$x^{n+1} = x^{n} + \Delta x\f$
    //!
    //! Iteration increments, that are needed internally in the single fields,
    //! have to be computed somewhere else.
    //!
    //! \param x  (i) composed vector that contains all field vectors
    //! \param sx (o) structural displacements
    //! \param fx (o) fluid velocities and pressure
    //! \param ax (o) ale displacements
    void extract_field_vectors(std::shared_ptr<const Core::LinAlg::Vector<double>> x,
        std::shared_ptr<const Core::LinAlg::Vector<double>>& sx,
        std::shared_ptr<const Core::LinAlg::Vector<double>>& fx,
        std::shared_ptr<const Core::LinAlg::Vector<double>>& ax) override;

    /// compute the Lagrange multiplier (FSI stresses) for the current time step
    void recover_lagrange_multiplier() override;

    /// compute the residual and incremental norms required for convergence check
    void build_convergence_norms() override;

    /// read restart data
    void read_restart(int step) override;

    /// output of fluid, structure & ALE-quantities and Lagrange multiplier
    void output() override;

    /*!
     * In case of a change in the fluid DOF row maps during the Newton loop (full Newton approach),
     * reset vectors accordingly.

     */
    void handle_fluid_dof_map_change_in_newton() override;

    /*!
     * Determine a change in fluid DOF map
     * \param (in) : DOF map of fluid increment vector
     * \return : true, in case of a mismatch between map of increment vector
     * and inner fluid DOF map after evaluation

     */
    bool has_fluid_dof_map_changed(const Core::LinAlg::Map& fluidincrementmap) override;

   private:
    /// build block vector from field vectors
    void setup_vector(Core::LinAlg::Vector<double>& f, const Core::LinAlg::Vector<double>& sv,
        const Core::LinAlg::Vector<double>& fv, const Core::LinAlg::Vector<double>& av,
        double fluidscale);

    /// access type-cast pointer to problem-specific fluid-wrapper
    const std::shared_ptr<Adapter::FluidFluidFSI>& fluid_field() { return MonolithicNoNOX::fluid_; }

    /// block system matrix
    // std::shared_ptr<OverlappingBlockMatrix> systemmatrix_;

    /// @name Matrix block transform objects
    /// Handle row and column map exchange for matrix blocks

    std::shared_ptr<Coupling::Adapter::MatrixRowColTransform> sggtransform_;
    std::shared_ptr<Coupling::Adapter::MatrixRowTransform> sgitransform_;
    std::shared_ptr<Coupling::Adapter::MatrixColTransform> sigtransform_;
    std::shared_ptr<Coupling::Adapter::MatrixColTransform> aigtransform_;

    std::shared_ptr<Coupling::Adapter::MatrixColTransform> fmiitransform_;
    std::shared_ptr<Coupling::Adapter::MatrixColTransform> fmgitransform_;

    std::shared_ptr<Coupling::Adapter::MatrixColTransform> fsaigtransform_;
    std::shared_ptr<Coupling::Adapter::MatrixColTransform> fsmgitransform_;

    ///@}

    /// @name infnorm scaling

    std::shared_ptr<Core::LinAlg::Vector<double>> srowsum_;
    std::shared_ptr<Core::LinAlg::Vector<double>> scolsum_;
    std::shared_ptr<Core::LinAlg::Vector<double>> arowsum_;
    std::shared_ptr<Core::LinAlg::Vector<double>> acolsum_;

    //@}

    /// @name Some quantities to recover the Lagrange multiplier at the end of each time step

    //! Lagrange multiplier \f$\lambda_\Gamma^n\f$ at the interface (ie condensed forces onto the
    //! structure) evaluated at old time step \f$t_n\f$ but needed for next time step \f$t_{n+1}\f$
    // lambda lives at the slave side (here at structure)
    std::shared_ptr<Core::LinAlg::Vector<double>> lambda_;

    //! interface force \f$f_{\Gamma,i+1}^{S,n+1}\f$ onto the structure at current iteration
    //! \f$i+1\f$
    // std::shared_ptr<const Core::LinAlg::Vector<double>> fgcur_;

    //! interface force \f$f_{\Gamma,i}^{S,n+1}\f$ onto the structure at previous iteration \f$i\f$
    // std::shared_ptr<const Core::LinAlg::Vector<double>> fgpre_;

    //! inner structural displacement increment \f$\Delta(\Delta d_{I,i+1}^{n+1})\f$ at current
    //! iteration \f$i+1\f$
    std::shared_ptr<Core::LinAlg::Vector<double>> ddiinc_;

    //! inner displacement solution of the structure at previous iteration
    std::shared_ptr<const Core::LinAlg::Vector<double>> solipre_;

    //! structural interface displacement increment \f$\Delta(\Delta d_{\Gamma,i+1}^{n+1})\f$ at
    //! current iteration \f$i+1\f$
    std::shared_ptr<Core::LinAlg::Vector<double>> ddginc_;

    //! interface displacement solution of the structure at previous iteration
    std::shared_ptr<const Core::LinAlg::Vector<double>> solgpre_;

    //! block \f$S_{\Gamma I,i+1}\f$ of structural matrix at current iteration \f$i+1\f$
    std::shared_ptr<const Core::LinAlg::SparseMatrix> sgicur_;

    //! block \f$S_{\Gamma\Gamma,i+1}\f$ of structural matrix at current iteration \f$i+1\f$
    std::shared_ptr<const Core::LinAlg::SparseMatrix> sggcur_;
    //@}
  };
}  // namespace FSI

FOUR_C_NAMESPACE_CLOSE

#endif
