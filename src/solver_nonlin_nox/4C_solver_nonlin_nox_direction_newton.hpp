// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_SOLVER_NONLIN_NOX_DIRECTION_NEWTON_HPP
#define FOUR_C_SOLVER_NONLIN_NOX_DIRECTION_NEWTON_HPP

#include "4C_config.hpp"

#include "4C_solver_nonlin_nox_forward_decl.hpp"

#include <NOX_Direction_Newton.H>  // base class

FOUR_C_NAMESPACE_OPEN

namespace NOX
{
  namespace Nln
  {
    namespace Direction
    {
      class Newton : public ::NOX::Direction::Newton
      {
       public:
        //! Constructor
        Newton(const Teuchos::RCP<::NOX::GlobalData>& gd, Teuchos::ParameterList& params);

        bool compute(::NOX::Abstract::Vector& dir, ::NOX::Abstract::Group& group,
            const ::NOX::Solver::Generic& solver) override;

       private:
        // throw NOX error
        void throw_error(const std::string& functionName, const std::string& errorMsg);

       private:
        //! NOX_Utils pointer
        Teuchos::RCP<::NOX::Utils> utils_;
      };
    }  // namespace Direction
  }  // namespace Nln
}  // namespace NOX

FOUR_C_NAMESPACE_CLOSE

#endif
