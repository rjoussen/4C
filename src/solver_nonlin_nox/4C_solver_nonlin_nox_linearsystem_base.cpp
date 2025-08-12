// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_solver_nonlin_nox_linearsystem_base.hpp"

#include "4C_utils_exceptions.hpp"

#include <NOX_Epetra_Scaling.H>

FOUR_C_NAMESPACE_OPEN

Teuchos::RCP<::NOX::Epetra::Scaling> NOX::Nln::LinearSystemBase::getScaling()
{
  FOUR_C_THROW("This is a mock and if you need it, most probably you are doing something wrong.");
  return Teuchos::null;
}

void NOX::Nln::LinearSystemBase::resetScaling(const Teuchos::RCP<::NOX::Epetra::Scaling>& s)
{
  (void)s;  // avoid unused parameter warning
  FOUR_C_THROW("This is a mock and if you need it, most probably you are doing something wrong.");
}

FOUR_C_NAMESPACE_CLOSE
