// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_CONTACT_NOXINTERFACE_HPP
#define FOUR_C_CONTACT_NOXINTERFACE_HPP

#include "4C_config.hpp"

#include "4C_solver_nonlin_nox_constraint_interface_required.hpp"

#include <vector>

FOUR_C_NAMESPACE_OPEN

namespace CONTACT
{
  class AbstractStrategy;

  class NoxInterface : public NOX::Nln::CONSTRAINT::Interface::Required
  {
   public:
    /// constructor
    NoxInterface();

    /// initialize important member variables
    void init(const std::shared_ptr<CONTACT::AbstractStrategy>& strategy_ptr);

    /** \brief Setup important new member variables
     *
     *  Supposed to be overloaded by derived classes. */
    virtual void setup();

    /// @name Supported basic interface functions
    /// @{
    //! Returns the constraint right-hand-side norms [derived]
    double get_constraint_rhs_norms(const Core::LinAlg::Vector<double>& F,
        NOX::Nln::StatusTest::QuantityType checkQuantity, ::NOX::Abstract::Vector::NormType type,
        bool isScaled) const override;

    /// Returns the root mean square (abbr.: RMS) of the Lagrange multiplier updates [derived]
    double get_lagrange_multiplier_update_rms(const Core::LinAlg::Vector<double>& xNew,
        const Core::LinAlg::Vector<double>& xOld, double aTol, double rTol,
        NOX::Nln::StatusTest::QuantityType checkQuantity,
        bool disable_implicit_weighting) const override;

    /// Returns the increment norm of the largange multiplier DoFs
    double get_lagrange_multiplier_update_norms(const Core::LinAlg::Vector<double>& xNew,
        const Core::LinAlg::Vector<double>& xOld, NOX::Nln::StatusTest::QuantityType checkQuantity,
        ::NOX::Abstract::Vector::NormType type, bool isScaled) const override;

    /// Returns the previous solution norm of the largange multiplier DoFs
    double get_previous_lagrange_multiplier_norms(const Core::LinAlg::Vector<double>& xOld,
        NOX::Nln::StatusTest::QuantityType checkQuantity, ::NOX::Abstract::Vector::NormType type,
        bool isScaled) const override;

    /// Returns the active set info [derived]
    enum ::NOX::StatusTest::StatusType get_active_set_info(
        enum NOX::Nln::StatusTest::QuantityType checkQuantity, int& activesetsize) const override;

    /// Returns the current active set map
    Teuchos::RCP<const Core::LinAlg::Map> get_current_active_set_map(
        enum NOX::Nln::StatusTest::QuantityType checkQuantity) const override;

    /// Returns the old active set map of the previous Newton step
    Teuchos::RCP<const Core::LinAlg::Map> get_old_active_set_map(
        enum NOX::Nln::StatusTest::QuantityType checkQuantity) const override;
    /// @}

    /// @name Merit function support functions
    /// @{

    double get_model_value(NOX::Nln::MeritFunction::MeritFctName name) const override;

    double get_linearized_model_terms(const Core::LinAlg::Vector<double>& dir,
        const enum NOX::Nln::MeritFunction::MeritFctName name,
        const enum NOX::Nln::MeritFunction::LinOrder linorder,
        const enum NOX::Nln::MeritFunction::LinType lintype) const override;

    /// @}

   protected:
    /// get the init indicator state
    inline const bool& is_init() const { return isinit_; };

    /// get the setup indicator state
    inline const bool& is_setup() const { return issetup_; };

    /// Check if init() has been called
    inline void check_init() const
    {
      if (not is_init()) FOUR_C_THROW("Call init() first!");
    };

    /// Check if init() and setup() have been called, yet.
    inline void check_init_setup() const
    {
      if (not is_init() or not is_setup()) FOUR_C_THROW("Call init() and setup() first!");
    };

    /// Access the underlying strategy
    const CONTACT::AbstractStrategy& strategy() const
    {
      check_init();
      return *strategy_ptr_;
    };

   protected:
    /// flag indicating if init() has been called
    bool isinit_;

    /// flag indicating if setup() has been called
    bool issetup_;

   private:
    std::shared_ptr<CONTACT::AbstractStrategy> strategy_ptr_;

    std::vector<std::shared_ptr<Core::LinAlg::Map>> cycling_maps_;
  };
}  // namespace CONTACT


FOUR_C_NAMESPACE_CLOSE

#endif
