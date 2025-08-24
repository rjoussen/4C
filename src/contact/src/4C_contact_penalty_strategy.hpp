// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_CONTACT_PENALTY_STRATEGY_HPP
#define FOUR_C_CONTACT_PENALTY_STRATEGY_HPP

#include "4C_config.hpp"

#include "4C_contact_abstract_strategy.hpp"

FOUR_C_NAMESPACE_OPEN


namespace CONTACT
{
  // forward declaration
  // class WearInterface;
  /*!
   \brief Contact solving strategy with regularization of Lagrangian multipliers,
   also known as Penalty Method or regularization. An Augmented Lagrangian version
   based on the Uzawa algorithm is included, too.

   This is a specialization of the abstract contact algorithm as defined in AbstractStrategy.
   For a more general documentation of the involved functions refer to CONTACT::AbstractStrategy.

   Refer also to the Semesterarbeit of Bernd Budich, 2009

   */
  class PenaltyStrategy : public AbstractStrategy
  {
   public:
    /*!
    \brief Standard constructor

    \param[in] dof_row_map Dof row map of underlying problem
    \param[in] NodeRowMap Node row map of underlying problem
    \param[in] params List of contact/parameters
    \param[in] interface All contact interface objects
    \param[in] spatialDim Spatial dimension of the problem
    \param[in] comm Communicator
    \param[in] alphaf Mid-point for Generalized-alpha time integration
    \param[in] maxdof Highest DOF number in global problem
    */
    PenaltyStrategy(const Core::LinAlg::Map* dof_row_map, const Core::LinAlg::Map* NodeRowMap,
        Teuchos::ParameterList params, std::vector<std::shared_ptr<CONTACT::Interface>> interface,
        const int spatialDim, const MPI_Comm& comm, const double alphaf, const int maxdof);

    /*!
    \brief Shared data constructor

    \param[in] strategyData Data container object
    \param[in] dof_row_map Dof row map of underlying problem
    \param[in] NodeRowMap Node row map of underlying problem
    \param[in] params List of contact/parameters
    \param[in] interface All contact interface objects
    \param[in] spatialDim Spatial dimension of the problem
    \param[in] comm Communicator
    \param[in] alphaf Mid-point for Generalized-alpha time integration
    \param[in] maxdof Highest DOF number in global problem
    */
    PenaltyStrategy(const std::shared_ptr<CONTACT::AbstractStrategyDataContainer>& data_ptr,
        const Core::LinAlg::Map* dof_row_map, const Core::LinAlg::Map* NodeRowMap,
        Teuchos::ParameterList params, std::vector<std::shared_ptr<CONTACT::Interface>> interface,
        const int spatialDim, const MPI_Comm& comm, const double alphaf, const int maxdof);

    //! @name Access methods

    /*!
    \brief Return L2-norm of active constraints

    */
    double constraint_norm() const override { return constrnorm_; }

    /*!
    \brief Return L2-norm of slip constraints

    */
    double constraint_norm_tan() const { return constrnormtan_; }


    /*!
    \brief Return initial penalty parameter for non-penetration

    */
    double initial_penalty() const override { return initialpenalty_; }

    /*!
    \brief Return initial penalty parameter for tangential direction

    */
    double initial_penalty_tan() const { return initialpenaltytan_; }

    //@}

    //! @name Evaluation methods

    /*!
    \brief Save nodal kappa-coefficients

    Before starting with the time integration, we have to calculate a nodal scaling factor,
    which will compensate the different integration area for computing the nodal weighted
    gap. Omitting this scaling, nodes on edges or boundaries would have a smaller weighted
    gap, even in case of a uniform physical gap. Hence, this scaling is of crucial importance
    for a penalty strategy since the weighted gap determines the lagrangian multipliers.

    */
    void save_reference_state(std::shared_ptr<const Core::LinAlg::Vector<double>> dis) override;

    /*!
    \brief Evaluate relative movement of contact bodies in predictor

    This is a tiny control routine for evaluating the relative movement of
    contact bodies in the predictor of an implicit time integration scheme.
    This evaluation (resetting) is ONLY necessary for penalty strategy and
    Uzawa augmented lagrange strategy, thus this tiny routine here.

    */
    void predict_relative_movement() override;

    /*!
    \brief Initialize general contact variables for next Newton step

    For a penalty strategy this involves the derivative matrix for the regularized lagrange
    multipliers.

    */
    void initialize() override;

    /*!
    \brief Evaluate contact

    For a penalty strategy this includes the evaluation of regularized forces
    in normal and tangential direction and results in a simple addition of extra
    stiffness contributions to kteff and extra contact forces to feff.

    */
    void evaluate_contact(std::shared_ptr<Core::LinAlg::SparseOperator>& kteff,
        std::shared_ptr<Core::LinAlg::Vector<double>>& feff) override;

    /*!
    \brief Evaluate frictional contact

    This includes the evaluation of of the frictional contact forces.

    */
    void evaluate_friction(std::shared_ptr<Core::LinAlg::SparseOperator>& kteff,
        std::shared_ptr<Core::LinAlg::Vector<double>>& feff) override;

    /*!
    \brief Reset penalty parameter to initial value

    When applying an Uzawa Augmented Lagrangian version of the penalty approach,
    the penalty parameter is sometimes updated during the Uzawa steps in
    order to accelerate convergence of the constraint norm. This increase
    in penalty stiffness can be dealt with, because at the time it is applied
    the constraint norm is already quite low. Yet, for a new time step, we have
    to come back to the initial penalty parameter. Thus, this method is called
    at the beginning of each time step and resets the penalty parameter to its initial value.

    */
    void reset_penalty() override;

    void modify_penalty() override;

    /*!
    \brief Initialize Uzawa step


    This method is called at the beginning of the second, third, ... Uzawa
    iterarion in order to create an of an out-of-balance force again. First,
    the contact force and stiffness terms are removed from feff and kteff.
    Then the LM and derivatives are updated (Uzawa AugmentedLagrange) and the new
    contact forces and stiffness terms are created by calling initialize()
    and finally evaluate().

    */
    void initialize_uzawa(std::shared_ptr<Core::LinAlg::SparseOperator>& kteff,
        std::shared_ptr<Core::LinAlg::Vector<double>>& feff) override;

    /*!
    \brief Compute L2-norm of active constraints

    In a classical penalty approach, the constraint norm is only monitored.
    When applying an Uzawa Augmented Lagrangian version, the constraint norm is the
    relevant stopping criterion of the Uzawa iteration. In order to accelerate
    convergence, a heuristic update formula for the penalty parameter is applied
    in this method, too.

    */
    void update_constraint_norm(int uzawaiter = 0) override;

    /*!
    \brief Store Lagrange multipliers for next Uzawa step

    A method ONLY called for the Uzawa Augmented Lagrangian version of the penalty method.
    At the end of an Uzawa step, the converged Lagrange multiplier value is stored
    in the variable zuzawa_, which is then used in the next Uzawa step.

    */
    void update_uzawa_augmented_lagrange() override;

    /*! \brief Compute force terms
     *
     *  \param cparams (in): parameter interface between the contact objects and the structural time
     * integration*/
    void evaluate_force(CONTACT::ParamsInterface& cparams) override;

    /*! \brief Compute force and stiffness terms
     *
     *  \param cparams (in): parameter interface between the contact objects and the structural time
     * integration*/
    void evaluate_force_stiff(CONTACT::ParamsInterface& cparams) override;

    /*! \brief Assemble force and stiffness terms to global vector and matrix */
    void assemble();

    /*! \brief Run at the beginning of the evaluate() routine
     *         set force evaluation flag
     *
     */
    void pre_evaluate(CONTACT::ParamsInterface& cparams) override;

    /*! \brief Run in the end of the evaluate() routine to reset
     *         force evaluation flag
     *
     *
     */
    void post_evaluate(CONTACT::ParamsInterface& cparams) override;


    /*! \brief Return the desired right-hand-side block pointer (read-only)
     *
     *  \remark Please note, that a nullptr pointer is returned, if no active contact
     *  contributions are present.
     *
     *  \param bt (in): Desired vector block type, e.g. displ, constraint,*/
    std::shared_ptr<const Core::LinAlg::Vector<double>> get_rhs_block_ptr(
        const enum CONTACT::VecBlockType& bt) const override;

    /*! \brief Return the desired matrix block pointer (read-only)
     *
     *  \remark Please note, that a nullptr pointer is returned, if no active contact
     *  contributions are present.
     *
     *  \param bt (in): Desired matrix block type, e.g. displ_displ, displ_lm, ...
     *  \param cparams (in): contact parameter interface (read-only) */
    std::shared_ptr<Core::LinAlg::SparseMatrix> get_matrix_block_ptr(
        const enum CONTACT::MatBlockType& bt,
        const CONTACT::ParamsInterface* cparams = nullptr) const override;

    //@}

    //! @name Empty functions (Lagrange contact)

    // All these functions only have functionality in Lagrange contact simulations,
    // thus they are defined empty here in the case of Penalty contact.
    std::shared_ptr<const Core::LinAlg::Map> get_old_active_row_nodes() const override
    {
      return nullptr;
    };
    std::shared_ptr<const Core::LinAlg::Map> get_old_slip_row_nodes() const override
    {
      return nullptr;
    };
    bool active_set_converged() const override { return true; }
    int active_set_steps() const override { return 0; }
    void reset_active_set() override {}
    void recover(std::shared_ptr<Core::LinAlg::Vector<double>> disi) override { return; };
    void build_saddle_point_system(std::shared_ptr<Core::LinAlg::SparseOperator> kdd,
        std::shared_ptr<Core::LinAlg::Vector<double>> fd,
        std::shared_ptr<Core::LinAlg::Vector<double>> sold,
        std::shared_ptr<Core::LinAlg::MapExtractor> dbcmaps,
        std::shared_ptr<Core::LinAlg::SparseOperator>& blockMat,
        std::shared_ptr<Core::LinAlg::Vector<double>>& blocksol,
        std::shared_ptr<Core::LinAlg::Vector<double>>& blockrhs) override
    {
      FOUR_C_THROW(
          "A penalty approach does not have Lagrange multiplier DOFs. So, saddle point system "
          "makes no sense here.");
    }
    void update_displacements_and_l_mincrements(std::shared_ptr<Core::LinAlg::Vector<double>> sold,
        std::shared_ptr<const Core::LinAlg::Vector<double>> blocksol) override
    {
      FOUR_C_THROW(
          "A penalty approach does not have Lagrange multiplier DOFs. So, saddle point system "
          "makes no sense here.");
    }
    void evaluate_constr_rhs() override {}
    void update_active_set() override {}
    void update_active_set_semi_smooth(const bool firstStepPredictor = false) override {}
    bool is_penalty() const override { return true; };
    void reset_lagrange_multipliers(
        const CONTACT::ParamsInterface& cparams, const Core::LinAlg::Vector<double>& xnew) override
    {
    }
    bool is_saddle_point_system() const override { return false; }
    bool is_condensed_system() const override { return false; }
    bool is_nitsche() const override { return false; }

    /*! \brief recover the current state
     *
     *  The main task of this method is to recover the Lagrange multiplier solution.
     *  The Lagrange multiplier solution will be stored inside the corresponding strategy
     *  and is necessary for different internal evaluation methods. If the Lagrange multiplier
     *  is condensed, this method is the right place to recover it from the displacement solution.
     *  If it is not condensed (saddle-point system) use the ResetLagrangeMultiplier routine
     * instead.
     *
     *  \param cparams (in): parameter interface between the contact objects and the structural time
     * integration \param xold    (in): old solution vector of the NOX solver \param dir     (in):
     * current search direction (in general NOT the actual step, keep in mind that the step length
     * can differ from 1.0) \param xnew    (in): new solution vector of the NOX solver
     *

     *  */
    void run_post_compute_x(const CONTACT::ParamsInterface& cparams,
        const Core::LinAlg::Vector<double>& xold, const Core::LinAlg::Vector<double>& dir,
        const Core::LinAlg::Vector<double>& xnew) override
    {
    }
    std::shared_ptr<const Core::LinAlg::Vector<double>> lagrange_multiplier_n(
        const bool& redist) const override;
    std::shared_ptr<const Core::LinAlg::Vector<double>> lagrange_multiplier_np(
        const bool& redist) const override;
    std::shared_ptr<const Core::LinAlg::Vector<double>> lagrange_multiplier_old() const override;
    std::shared_ptr<const Core::LinAlg::Map> lm_dof_row_map_ptr(const bool& redist) const override;

   protected:
    //! derived
    std::vector<std::shared_ptr<CONTACT::Interface>>& interfaces() override { return interface_; }

    //! derived
    const std::vector<std::shared_ptr<CONTACT::Interface>>& interfaces() const override
    {
      return interface_;
    }

    // don't want = operator and cctor
    PenaltyStrategy operator=(const PenaltyStrategy& old) = delete;
    PenaltyStrategy(const PenaltyStrategy& old) = delete;

    std::vector<std::shared_ptr<Interface>> interface_;  // contact interfaces
    std::shared_ptr<Core::LinAlg::SparseMatrix>
        linzmatrix_;            // global matrix LinZ with derivatives of LM
    double constrnorm_;         // L2-norm of normal contact constraints
    double constrnormtan_;      // L2-norm of tangential contact constraints
    double initialpenalty_;     // initial penalty parameter
    double initialpenaltytan_;  // initial tangential penalty parameter
    bool evalForceCalled_;      //< flag for evaluate force call

    std::shared_ptr<Core::LinAlg::Vector<double>> fc_;  //< contact penalty force
    std::shared_ptr<Core::LinAlg::SparseMatrix> kc_;    //< contact penalty stiffness

  };  // class PenaltyStrategy
}  // namespace CONTACT

FOUR_C_NAMESPACE_CLOSE

#endif
