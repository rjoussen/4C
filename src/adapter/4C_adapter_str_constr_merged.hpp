// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_ADAPTER_STR_CONSTR_MERGED_HPP
#define FOUR_C_ADAPTER_STR_CONSTR_MERGED_HPP

/*----------------------------------------------------------------------*/
/* headers */
#include "4C_config.hpp"

#include "4C_adapter_str_fsiwrapper.hpp"
#include "4C_linalg_map.hpp"
#include "4C_linalg_vector.hpp"
#include "4C_utils_parameter_list.fwd.hpp"

#include <memory>

FOUR_C_NAMESPACE_OPEN


// forward declarations
namespace Solid
{
  namespace Aux
  {
    class MapExtractor;
  }
}  // namespace Solid

namespace Core::LinAlg
{
  class MapExtractor;
}

/*----------------------------------------------------------------------*/
namespace Adapter
{
  /*====================================================================*/
  /*!
   * \brief Adapter to constrained structural time integration.
   * This class wraps one of the standard adapters for structural time
   * integration. The results are modified and/or merged to account for the
   * additional degrees of freedom of the lagrange multipliers.
   *

   */
  class StructureConstrMerged : public FSIStructureWrapper
  {
   public:
    /// Constructor
    StructureConstrMerged(std::shared_ptr<Structure> stru);

    /// setup this object
    void setup() override;

    /// initial guess of Newton's method
    std::shared_ptr<const Core::LinAlg::Vector<double>> initial_guess() override;

    /// right-hand-side of Newton's method
    std::shared_ptr<const Core::LinAlg::Vector<double>> rhs() override;

    /// unknown displacements at \f$t_{n+1}\f$
    std::shared_ptr<const Core::LinAlg::Vector<double>> dispnp() const override;

    /// known displacements at \f$t_{n}\f$
    std::shared_ptr<const Core::LinAlg::Vector<double>> dispn() const override;

    /*! \brief known velocity at \f$t_{n}\f$
     *
     *  Lagrange multiplier does not have a time derivative. Though we need a map
     *  including the Lagrange multiplier, thus, we include it and set it to zero.
     */
    std::shared_ptr<const Core::LinAlg::Vector<double>> veln() const override;

    /*! known acceleration at \f$t_{n}\f$
     *
     *  Lagrange multiplier does not have a time derivative. Though we need a map
     *  including the Lagrange multiplier, thus, we include it and set it to zero.
     */
    std::shared_ptr<const Core::LinAlg::Vector<double>> accn() const override;

    /// dof map of vector of unknowns
    std::shared_ptr<const Core::LinAlg::Map> dof_row_map() override;

    /// apply interface forces to structural solver
    ///
    /// This prepares a new solve of the structural field within one time
    /// step. The middle values are newly created.
    ///
    /// \note This is not yet the most efficient implementation.
    void apply_interface_forces_temporary_deprecated(
        std::shared_ptr<Core::LinAlg::Vector<double>> iforce) override;

    /// direct access to system matrix
    std::shared_ptr<Core::LinAlg::SparseMatrix> system_matrix() override;

    /// direct access to system matrix
    std::shared_ptr<Core::LinAlg::BlockSparseMatrixBase> block_system_matrix() override;

    /// update displacement and evaluate elements
    void evaluate(std::shared_ptr<const Core::LinAlg::Vector<double>>
            dispstepinc  ///< solution increment between time step n and n+1
        ) override;

    //! Return MapExtractor for Dirichlet boundary conditions
    std::shared_ptr<const Core::LinAlg::MapExtractor> get_dbc_map_extractor() override
    {
      return structure_->get_dbc_map_extractor();
    };

    /// domain map of system matrix
    const Core::LinAlg::Map& domain_map() const override;

    /// are there any algebraic constraints?
    bool have_constraint() override { return structure_->have_constraint(); };

    /// Return bool indicating if constraints are defined
    std::shared_ptr<Constraints::ConstrManager> get_constraint_manager() override
    {
      return structure_->get_constraint_manager();
    };

    Inpar::Solid::StcScale get_stc_algo() override { return structure_->get_stc_algo(); };

    std::shared_ptr<Core::LinAlg::SparseMatrix> get_stc_mat() override
    {
      FOUR_C_THROW("FSI with merged structural constraints does not work in combination with STC!");
      return structure_->get_stc_mat();
    }

    //! Update iteration
    //! Add residual increment to Lagrange multipliers stored in Constraint manager
    void update_iter_incr_constr(
        std::shared_ptr<Core::LinAlg::Vector<double>> lagrincr  ///< Lagrange multiplier increment
        ) override
    {
      structure_->update_iter_incr_constr(lagrincr);
    }

    /// @name Apply interface forces

    //@}

    /// Integrate from t1 to t2
    int integrate() override { return structure_->integrate(); }

   private:
    /// the constraint map setup for full <-> stuct+constr transition
    std::shared_ptr<Core::LinAlg::MapExtractor> conmerger_;

    /// the complete non-overlapping degree of freedom row map for structure and lagrange
    /// multipliers
    std::shared_ptr<Core::LinAlg::Map> dofrowmap_;

    /// @name local copies of input parameters
    //{@
    std::shared_ptr<Core::FE::Discretization> discret_;  ///< the discretization
    std::shared_ptr<Teuchos::ParameterList>
        sdynparams_;  ///< dynamic control flags ... used, but could/should be circumvented
    std::shared_ptr<Teuchos::ParameterList> xparams_;         ///< eXtra input parameters
    std::shared_ptr<Core::LinAlg::Solver> solver_;            ///< the linear solver
    std::shared_ptr<Core::IO::DiscretizationWriter> output_;  ///< the output writer

    //@}

    /// flag indicating if setup() was called
    bool issetup_;

  };  // class StructureConstrained

}  // namespace Adapter

/*----------------------------------------------------------------------*/
FOUR_C_NAMESPACE_CLOSE

#endif
