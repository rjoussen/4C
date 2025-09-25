// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_SOLVER_NONLIN_NOX_GROUP_BASE_HPP
#define FOUR_C_SOLVER_NONLIN_NOX_GROUP_BASE_HPP

#include "4C_config.hpp"

#include <NOX_Epetra_Group.H>

FOUR_C_NAMESPACE_OPEN

namespace NOX
{
  namespace Nln
  {
    class GroupBase : public ::NOX::Epetra::Group
    {
     public:
      GroupBase(Teuchos::ParameterList& printParams,
          const Teuchos::RCP<::NOX::Epetra::Interface::Required>& i, const ::NOX::Epetra::Vector& x,
          const Teuchos::RCP<::NOX::Epetra::LinearSystem>& linSys);

      GroupBase(const NOX::Nln::GroupBase& source, ::NOX::CopyType type);
    };  // class GroupBase
  }  // namespace Nln
}  // namespace NOX

FOUR_C_NAMESPACE_CLOSE

#endif
