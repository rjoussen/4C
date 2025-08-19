// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_POROELAST_MONOLITHICSTRUCTURESPLIT_HPP
#define FOUR_C_POROELAST_MONOLITHICSTRUCTURESPLIT_HPP

#include "4C_config.hpp"

#include "4C_coupling_adapter_converter.hpp"
#include "4C_poroelast_monolithicsplit.hpp"

FOUR_C_NAMESPACE_OPEN


namespace PoroElast
{
  //! monolithic structure split for condensing DOFs, when using the brinkman-equation
  class MonolithicStructureSplit : public MonolithicSplit
  {
   public:
    explicit MonolithicStructureSplit(MPI_Comm comm, const Teuchos::ParameterList& timeparams,
        std::shared_ptr<Core::LinAlg::MapExtractor> porosity_splitter);

    /*! do the setup for the monolithic system


     1.) setup coupling
     2.) get maps for all blocks in the system (and for the whole system as well)
     create combined map
     3.) create system matrix


     \note We want to do this setup after reading the restart information, not
     directly in the constructor. This is necessary since during restart (if
     read_mesh is called), the dofmaps for the blocks might get invalid.
     */
    //! Setup the monolithic system
    void setup_system() override;

    //! setup composed right hand side from field solvers
    void setup_rhs(bool firstcall = false) override;

    //! setup composed system matrix from field solvers
    void setup_system_matrix(Core::LinAlg::BlockSparseMatrixBase& mat) override;

   private:
    //! build block vector from field vectors
    void setup_vector(Core::LinAlg::Vector<double>& f, const Core::LinAlg::Vector<double>& sv,
        const Core::LinAlg::Vector<double>& fv, double fluidscale);

    //! extract the field vectors from a given composed vector
    /*!
     \param x  (i) composed vector that contains all field vectors
     \param sx (o) structural vector (e.g. displacements)
     \param fx (o) fluid vector (e.g. velocities and pressure)
     */
    void extract_field_vectors(std::shared_ptr<const Core::LinAlg::Vector<double>> x,
        std::shared_ptr<const Core::LinAlg::Vector<double>>& sx,
        std::shared_ptr<const Core::LinAlg::Vector<double>>& fx, bool firstcall = false) override;

    //! recover Lagrange multiplier \f$\lambda_\Gamma\f$ at the interface at the end of each time
    //! step (i.e. condensed forces onto the structure) needed for rhs in next time step
    void recover_lagrange_multiplier_after_time_step() override;

    //! @name matrix transformation
    //! transform object for structure interface matrix \f$S_{\Gamma \Gamma}\f$
    std::shared_ptr<Coupling::Adapter::MatrixRowColTransform> sggtransform_;
    //! transform object for structure interface matrix \f$S_{\Gamma I}\f$
    std::shared_ptr<Coupling::Adapter::MatrixRowTransform> sgitransform_;
    //! transform object for structure interface matrix \f$S_{I \Gamma}\f$
    std::shared_ptr<Coupling::Adapter::MatrixColTransform> sigtransform_;
    //! transform object for structure coupling matrix \f$C_{\Gamma \Gamma}^S\f$
    std::shared_ptr<Coupling::Adapter::MatrixRowTransform> csggtransform_;
    //! transform object for fluid coupling matrix \f$C_{\Gamma \Gamma}^G\f$
    std::shared_ptr<Coupling::Adapter::MatrixColTransform> cfggtransform_;
    //! transform object for structure coupling matrix \f$C_{\Gamma I}^S\f$
    std::shared_ptr<Coupling::Adapter::MatrixRowTransform> csgitransform_;
    //! transform object for fluid coupling matrix \f$C_{I \Gamma}^F\f$
    std::shared_ptr<Coupling::Adapter::MatrixColTransform> cfigtransform_;

    //!@}

    //! @name Some quantities to recover the Lagrange multiplier at the end of each time step

    //! block \f$S_{\Gamma I,i+1}\f$ of structural matrix at current iteration \f$i+1\f$
    std::shared_ptr<const Core::LinAlg::SparseOperator> sgicur_;

    //! block \f$S_{\Gamma\Gamma,i+1}\f$ of structural matrix at current iteration \f$i+1\f$
    std::shared_ptr<const Core::LinAlg::SparseOperator> sggcur_;

    //! block \f$S_{\Gamma\Gamma,i+1}\f$ of structural matrix at current iteration \f$i+1\f$
    std::shared_ptr<const Core::LinAlg::SparseOperator> cgicur_;

    //! block \f$S_{\Gamma\Gamma,i+1}\f$ of structural matrix at current iteration \f$i+1\f$
    std::shared_ptr<const Core::LinAlg::SparseOperator> cggcur_;

    //!@}
  };

}  // namespace PoroElast

FOUR_C_NAMESPACE_CLOSE

#endif
