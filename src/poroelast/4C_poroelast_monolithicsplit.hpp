// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_POROELAST_MONOLITHICSPLIT_HPP
#define FOUR_C_POROELAST_MONOLITHICSPLIT_HPP


#include "4C_config.hpp"

#include "4C_coupling_adapter.hpp"
#include "4C_poroelast_monolithic.hpp"

FOUR_C_NAMESPACE_OPEN

namespace PoroElast
{
  //! base class for monolithic approaches, when the system is split for some reason
  //! (most of the time this means some dofs are condensed)
  class MonolithicSplit : public Monolithic
  {
   public:
    explicit MonolithicSplit(MPI_Comm comm, const Teuchos::ParameterList& timeparams,
        std::shared_ptr<Core::LinAlg::MapExtractor> porosity_splitter);

    //! Setup the monolithic system (depends on which field is split)
    void setup_system() override = 0;

    //! setup composed right hand side from field solvers (depends on which field is split)
    void setup_rhs(bool firstcall = false) override = 0;

    //! setup composed system matrix from field solvers (depends on which field is split)
    void setup_system_matrix(Core::LinAlg::BlockSparseMatrixBase& mat) override = 0;

    //! start a new time step
    void prepare_time_step() override;

    //! recover Lagrange multiplier \f$\lambda_\Gamma\f$ at the interface at the end of each time
    //! step (i.e. condensed forces onto the structure) needed for rhs in next time step
    void recover_lagrange_multiplier_after_time_step() override = 0;

    //! inner newton iteration
    void solve() override;

   protected:
    //! @name transfer helpers

    //! field transform (interface only)
    std::shared_ptr<Core::LinAlg::Vector<double>> structure_to_fluid_at_interface(
        const Core::LinAlg::Vector<double>& iv) const;

    //! field transform (interface only)
    std::shared_ptr<Core::LinAlg::Vector<double>> fluid_to_structure_at_interface(
        const Core::LinAlg::Vector<double>& iv) const;

    //!@}

    //! combined DBC map
    //! unique map of all dofs that should be constrained with DBC
    void build_combined_dbc_map() override;

    //! map containing the dofs with Dirichlet BC and FSI Coupling Condition on structure side
    std::shared_ptr<Core::LinAlg::Map> fsidbc_map();

    //! setup of coupling object and systemmatrixes
    virtual void setup_coupling_and_matrices();

    //! coupling of fluid and structure (interface only), only needed by algorithms, who perform a
    //! split, i.e. structure or fluid split.
    std::shared_ptr<Coupling::Adapter::Coupling> icoupfs_;

    //! flag indicating whether there are no slip conditions to be evaluated at the interface
    bool evaluateinterface_;

    //! map containing DOFs with both fsi- and DBC conditions
    std::shared_ptr<Core::LinAlg::Map> fsibcmap_;

    //! map extractor DOFs with both fsi- and DBC conditions
    std::shared_ptr<Core::LinAlg::MapExtractor> fsibcextractor_;

    //! @name Some quantities to recover the Lagrange multiplier at the end of each time step

    //! Lagrange multiplier \f$\lambda_\Gamma^n\f$ at the interface (ie condensed forces onto the
    //! structure) evaluated at old time step \f$t_n\f$ but needed for next time step \f$t_{n+1}\f$
    std::shared_ptr<Core::LinAlg::Vector<double>> lambda_;

    //! interface force \f$f_{\Gamma,i+1}^{S,n+1}\f$ onto the structure at current iteration
    //! \f$i+1\f$
    std::shared_ptr<const Core::LinAlg::Vector<double>> fgcur_;

    //! inner structural displacement increment \f$\Delta(\Delta d_{I,i+1}^{n+1})\f$ at current
    //! iteration \f$i+1\f$
    std::shared_ptr<Core::LinAlg::Vector<double>> ddiinc_;

    //! inner fluid velocity increment \f$\Delta(\Delta u_{I,i+1}^{n+1})\f$ at current iteration
    //! \f$i+1\f$
    std::shared_ptr<Core::LinAlg::Vector<double>> duiinc_;

    //! inner displacement solution of the structure at previous iteration
    std::shared_ptr<const Core::LinAlg::Vector<double>> solipre_;

    //! inner velocity/pressure solution of the fluid at previous iteration
    std::shared_ptr<const Core::LinAlg::Vector<double>> solivelpre_;

    //! structural interface displacement increment \f$\Delta(\Delta d_{\Gamma,i+1}^{n+1})\f$ at
    //! current iteration \f$i+1\f$
    std::shared_ptr<Core::LinAlg::Vector<double>> ddginc_;

    //! fluid interface velocity increment \f$\Delta(\Delta d_{\Gamma,i+1}^{n+1})\f$ at current
    //! iteration \f$i+1\f$
    std::shared_ptr<Core::LinAlg::Vector<double>> duginc_;

    //! interface displacement solution of the structure at previous iteration
    std::shared_ptr<const Core::LinAlg::Vector<double>> solgpre_;

    //! interface displacement solution of the fluid at previous iteration
    std::shared_ptr<const Core::LinAlg::Vector<double>> solgvelpre_;

    //!@}

    //! interface increment for dof with dirichlet condition on fsi-interface
    std::shared_ptr<Core::LinAlg::Vector<double>> ddi_;
  };
}  // namespace PoroElast

FOUR_C_NAMESPACE_CLOSE

#endif
