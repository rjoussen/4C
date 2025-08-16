// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_CONTACT_MONOCOUPLED_LAGRANGE_STRATEGY_HPP
#define FOUR_C_CONTACT_MONOCOUPLED_LAGRANGE_STRATEGY_HPP

#include "4C_config.hpp"

#include "4C_contact_lagrange_strategy.hpp"
#include "4C_utils_exceptions.hpp"

FOUR_C_NAMESPACE_OPEN


namespace CONTACT
{
  /*!
   \brief Contact solving strategy with (standard/dual) Lagrangian multipliers.

   This is a specialization of the abstract contact algorithm as defined in AbstractStrategy.
   For a more general documentation of the involved functions refer to CONTACT::AbstractStrategy.

   */
  class MonoCoupledLagrangeStrategy : public LagrangeStrategy
  {
   public:
    /*!
    \brief Standard Constructor

    */
    MonoCoupledLagrangeStrategy(
        const std::shared_ptr<CONTACT::AbstractStrategyDataContainer>& data_ptr,
        const Core::LinAlg::Map* dof_row_map, const Core::LinAlg::Map* NodeRowMap,
        Teuchos::ParameterList params, std::vector<std::shared_ptr<CONTACT::Interface>> interface,
        int dim, MPI_Comm comm, double alphaf, int maxdof);


    // Overload CONTACT::AbstractStrategy::apply_force_stiff_cmt as this is called in the structure
    // --> to early for monolithically coupled algorithms!
    void apply_force_stiff_cmt(std::shared_ptr<Core::LinAlg::Vector<double>> dis,
        std::shared_ptr<Core::LinAlg::SparseOperator>& kt,
        std::shared_ptr<Core::LinAlg::Vector<double>>& f, const int step, const int iter,
        bool predictor) override
    {
      if (has_to_evaluate_ && 0)
        FOUR_C_THROW(
            "MonoCoupledLagrangeStrategy::You have to call apply_force_stiff_cmt_coupled() for "
            "Contact "
            "Evaluation!");  // what to do in the predictor?
      has_to_evaluate_ = true;
      return;
    };

    // Overload CONTACT::LagrangeStrategy::recover as this is called in the structure --> no
    // enough information available for monolithically coupled algorithms!
    void recover(std::shared_ptr<Core::LinAlg::Vector<double>> disi) override
    {
      if (has_to_recover_ && 0)
        FOUR_C_THROW(
            "MonoCoupledLagrangeStrategy::You have to call RecoverCoupled() for Contact Recovery!");
      has_to_recover_ = true;
      return;
    };

    //! @name Access methods

    //@}

    //! @name Evaluation methods

    // Alternative Method to CONTACT::AbstractStrategy::apply_force_stiff_cmt for monolithically
    // coupled algorithms
    virtual void apply_force_stiff_cmt_coupled(std::shared_ptr<Core::LinAlg::Vector<double>> dis,
        std::shared_ptr<Core::LinAlg::SparseOperator>& k_ss,
        std::map<int, std::shared_ptr<Core::LinAlg::SparseOperator>*> k_sx,
        std::shared_ptr<Core::LinAlg::Vector<double>>& rhs_s, const int step, const int iter,
        bool predictor);

    // Alternative Method to CONTACT::AbstractStrategy::apply_force_stiff_cmt for monolithically
    // coupled algorithms
    virtual void apply_force_stiff_cmt_coupled(std::shared_ptr<Core::LinAlg::Vector<double>> dis,
        std::shared_ptr<Core::LinAlg::SparseOperator>& k_ss,
        std::shared_ptr<Core::LinAlg::SparseOperator>& k_sx,
        std::shared_ptr<Core::LinAlg::Vector<double>>& rhs_s, const int step, const int iter,
        bool predictor);

    // Alternative Method to CONTACT::LagrangeStrategy::recover as this is called in the structure
    // --> no enough information available for monolithically coupled algorithms!
    /*!
    \brief Recovery method

    We only recover the Lagrange multipliers here, which had been
    statically condensated during the setup of the global problem!
    Optionally satisfaction or violation of the contact boundary
    conditions can be checked, too.*/
    virtual void recover_coupled(std::shared_ptr<Core::LinAlg::Vector<double>> disi,
        std::map<int, std::shared_ptr<Core::LinAlg::Vector<double>>> inc);

    virtual void recover_coupled(std::shared_ptr<Core::LinAlg::Vector<double>> disi,
        std::shared_ptr<Core::LinAlg::Vector<double>> inc);

    void evaluate_off_diag_contact(std::shared_ptr<Core::LinAlg::SparseOperator>& kteff,
        int Column_Block_Id);  // condensation for all off diagonal matrixes k_s? in monolithically
                               // coupled problems!

   protected:
    // don't want = operator and cctor
    MonoCoupledLagrangeStrategy operator=(const MonoCoupledLagrangeStrategy& old) = delete;
    MonoCoupledLagrangeStrategy(const MonoCoupledLagrangeStrategy& old) = delete;

    void save_coupling_matrices(std::shared_ptr<Core::LinAlg::SparseMatrix> dhat,
        std::shared_ptr<Core::LinAlg::SparseMatrix> mhataam,
        std::shared_ptr<Core::LinAlg::SparseMatrix> invda) override;

    std::map<int, std::shared_ptr<Core::LinAlg::SparseOperator>>
        csx_s_;  // offdiagonal coupling stiffness blocks on slave side!

    std::shared_ptr<Core::LinAlg::SparseMatrix> dhat_;
    std::shared_ptr<Core::LinAlg::SparseMatrix> mhataam_;
    std::shared_ptr<Core::LinAlg::SparseMatrix> invda_;

    std::shared_ptr<Core::LinAlg::Vector<double>>
        lambda_;  // current vector of Lagrange multipliers(for poro no pen.) at t_n+1
    std::shared_ptr<Core::LinAlg::Vector<double>>
        lambdaold_;  // old vector of Lagrange multipliers(for poro no pen.) at t_n

    //! pure usage safety flags
    bool has_to_evaluate_;  // checks if apply_force_stiff_cmt_coupled() after every call of
                            // apply_force_stiff_cmt()
    bool has_to_recover_;   // checks if RecoverCoupled() after every call of recover()

  };  // class MonoCoupledLagrangeStrategy
}  // namespace CONTACT


FOUR_C_NAMESPACE_CLOSE

#endif
