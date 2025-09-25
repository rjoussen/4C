// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_solver_nonlin_nox_group_base.hpp"

FOUR_C_NAMESPACE_OPEN

NOX::Nln::GroupBase::GroupBase(Teuchos::ParameterList& printParams,
    const Teuchos::RCP<::NOX::Epetra::Interface::Required>& i, const ::NOX::Epetra::Vector& x,
    const Teuchos::RCP<::NOX::Epetra::LinearSystem>& linSys)
    : ::NOX::Epetra::Group(printParams, i, x, linSys)
{
}

NOX::Nln::GroupBase::GroupBase(const NOX::Nln::GroupBase& source, ::NOX::CopyType type)
    : ::NOX::Epetra::Group(source, type)
{
}

FOUR_C_NAMESPACE_CLOSE
