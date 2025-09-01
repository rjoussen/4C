// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_CONSTRAINT_MULTIPOINTCONSTRAINT_HPP
#define FOUR_C_CONSTRAINT_MULTIPOINTCONSTRAINT_HPP

#include "4C_config.hpp"

#include "4C_constraint.hpp"

FOUR_C_NAMESPACE_OPEN



namespace Constraints
{
  /*!
  \brief This pure virtual class is the common interface for multi point constraints.
  It is derived from the basic constraint class.
  */
  class MPConstraint : public Constraint
  {
   public:
    /*!
    \brief Standard Constructor
    */
    MPConstraint(
        std::shared_ptr<Core::FE::Discretization> discr,  ///< discretization constraint lives on
        const std::string& conditionname,  ///< Name of condition to create constraint from
        int& minID,                        ///< minimum constraint or monitor ID so far
        int& maxID                         ///< maximum constraint or monitor ID so far
    );

    /*!
        \brief Alternative Constructor
    */
    MPConstraint(
        std::shared_ptr<Core::FE::Discretization> discr,  ///< discretization constraint lives on
        const std::string& conditionname  ///< Name of condition to create constraint from
    );

    /*!
        \brief Destructor
    */
    virtual ~MPConstraint() = default;

    /// Set state of the underlying constraint discretization
    void set_constr_state(const std::string& state,  ///< name of state to set
        const Core::LinAlg::Vector<double>& V        ///< values to set
    );

    /// initialization routine called by the manager ctor to get correct reference base values and
    /// activating the right conditions at the beginning
    virtual void initialize(
        Teuchos::ParameterList&
            params,  ///< parameter list to communicate between elements and discretization
        std::shared_ptr<Core::LinAlg::Vector<double>>
            systemvector3  ///< distributed vector that may be filled
                           ///< by assembly of element contributions
        ) = 0;

    /// initialization routine called at restart to activate the right conditions
    virtual void initialize(const double& time  ///< current time
        ) = 0;

    //! Evaluate routine to call from outside. In here the right action is determined and the
    //! #evaluate_constraint routine is called
    virtual void evaluate(
        Teuchos::ParameterList&
            params,  ///< parameter list to communicate between elements and discretization
        std::shared_ptr<Core::LinAlg::SparseOperator>
            systemmatrix1,  ///< sparse matrix that may be filled by assembly of element
                            ///< contributions
        std::shared_ptr<Core::LinAlg::SparseOperator>
            systemmatrix2,  ///< sparse (rectangular) matrix that may be filled by assembly of
                            ///< element contributions
        std::shared_ptr<Core::LinAlg::Vector<double>>
            systemvector1,  ///< distributed vector that may be filled by
                            ///< assembly of element contributions
        std::shared_ptr<Core::LinAlg::Vector<double>>
            systemvector2,  ///< distributed vector that may be filled by
                            ///< assembly of element contributions
        std::shared_ptr<Core::LinAlg::Vector<double>>
            systemvector3  ///< distributed vector that may be filled
                           ///< by assembly of element contributions
        ) = 0;

    //! Is there a constraint defined in this class?
    bool have_constraint() { return constrtype_ != none; };

   protected:
    // cctor
    MPConstraint(const MPConstraint& old);


    //! additional discretization consisting of constraint elements
    std::map<int, std::shared_ptr<Core::FE::Discretization>> constraintdis_;

    //! Evaluate constraint discretization and assemble the results
    virtual void evaluate_constraint(
        std::shared_ptr<Core::FE::Discretization> disc,  ///< discretization to evaluate
        Teuchos::ParameterList&
            params,  ///< parameter list to communicate between elements and discretization
        std::shared_ptr<Core::LinAlg::SparseOperator>
            systemmatrix1,  ///< sparse matrix that may be filled by assembly of element
                            ///< contributions
        std::shared_ptr<Core::LinAlg::SparseOperator>
            systemmatrix2,  ///< sparse (rectangular) matrix that may be filled by assembly of
                            ///< element contributions
        std::shared_ptr<Core::LinAlg::Vector<double>>
            systemvector1,  ///< distributed vector that may be filled by
                            ///< assembly of element contributions
        std::shared_ptr<Core::LinAlg::Vector<double>>
            systemvector2,  ///< distributed vector that may be filled by
                            ///< assembly of element contributions
        std::shared_ptr<Core::LinAlg::Vector<double>>
            systemvector3  ///< distributed vector that may be filled
                           ///< by assembly of element contributions
        ) = 0;

    //! creating a new discretization based on conditions containing constraint elements
    virtual std::map<int, std::shared_ptr<Core::FE::Discretization>>
    create_discretization_from_condition(std::shared_ptr<Core::FE::Discretization> actdisc,
        std::span<const Core::Conditions::Condition*>
            constrcond,                   ///< conditions as discretization basis
        const std::string& discret_name,  ///< name of new discretization
        const std::string& element_name,  ///< name of element type to create
        int& startID) = 0;


    //    /// find col node map so that we can evaluate the constraint elements
    //    std::shared_ptr<Core::LinAlg::Map> ComputeNodeColMap(
    //             const std::shared_ptr<Core::FE::Discretization> sourcedis,  ///< standard
    //             discretization we want to redistribute const
    //             std::shared_ptr<Core::FE::Discretization> constraintdis
    //             ///< constraint discretization prescribing ghosting ) const;

  };  // class
}  // namespace Constraints
FOUR_C_NAMESPACE_CLOSE

#endif
