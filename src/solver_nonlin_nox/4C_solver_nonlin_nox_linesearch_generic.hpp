// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_SOLVER_NONLIN_NOX_LINESEARCH_GENERIC_HPP
#define FOUR_C_SOLVER_NONLIN_NOX_LINESEARCH_GENERIC_HPP

#include "4C_config.hpp"

#include "4C_solver_nonlin_nox_inner_statustest_interface_required.hpp"

#include <NOX_LineSearch_Generic.H>

FOUR_C_NAMESPACE_OPEN

namespace NOX
{
  namespace Nln
  {
    namespace LineSearch
    {
      class PrePostOperator;

      class Generic : public ::NOX::LineSearch::Generic,
                      public NOX::Nln::Inner::StatusTest::Interface::Required
      {
       public:
        //! constructor
        Generic() = default;


        //! @name NOX::Nln::LineSearch::Generic
        //! @{
        //! returns the slope in the current search direction
        virtual const ::NOX::Abstract::Vector& get_search_direction() const = 0;

        //! returns the stepSize
        virtual double get_step_length() const = 0;

        //! sets the stepSize
        virtual void set_step_length(double step) = 0;
        //! @}

        //! @name ::NOX::LineSearch::Generic
        //! @{
        bool compute(::NOX::Abstract::Group& grp, double& step, const ::NOX::Abstract::Vector& dir,
            const ::NOX::Solver::Generic& s) override = 0;
        //! @}

        //! @name NOX::Nln::Inner::StatusTest::Interface::Required
        //! @{
        //! get the number of line search iterations
        int get_num_iterations() const override = 0;

        //! get the merit function
        const ::NOX::MeritFunction::Generic& get_merit_function() const override = 0;
        //! @}

       protected:
        Teuchos::RCP<PrePostOperator> prePostOperatorPtr_ = Teuchos::null;
      };
    }  // namespace LineSearch
  }  // namespace Nln
}  // namespace NOX

FOUR_C_NAMESPACE_CLOSE

#endif
