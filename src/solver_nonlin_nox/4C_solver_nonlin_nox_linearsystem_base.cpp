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

void NOX::Nln::LinearSystemBase::setPrecOperatorForSolve(
    const Teuchos::RCP<const Epetra_Operator>& solvePrecOp)
{
  (void)solvePrecOp;

  FOUR_C_THROW("setPrecOperatorForSolve() is not implemented for NOX::Nln::LinearSystemBase.");
}

bool NOX::Nln::LinearSystemBase::isPreconditionerConstructed() const { return false; }

bool NOX::Nln::LinearSystemBase::hasPreconditioner() const { return false; }

Teuchos::RCP<const Epetra_Operator> NOX::Nln::LinearSystemBase::getGeneratedPrecOperator() const
{
  return Teuchos::null;
}

Teuchos::RCP<Epetra_Operator> NOX::Nln::LinearSystemBase::getGeneratedPrecOperator()
{
  return Teuchos::null;
}

bool NOX::Nln::LinearSystemBase::createPreconditioner(
    const ::NOX::Epetra::Vector& x, Teuchos::ParameterList& p, bool recomputeGraph) const
{
  (void)x;
  (void)p;
  (void)recomputeGraph;

  return false;
}

bool NOX::Nln::LinearSystemBase::destroyPreconditioner() const { return false; }

bool NOX::Nln::LinearSystemBase::recomputePreconditioner(
    const ::NOX::Epetra::Vector& x, Teuchos::ParameterList& linearSolverParams) const
{
  (void)x;
  (void)linearSolverParams;

  return false;
}

::NOX::Epetra::LinearSystem::PreconditionerReusePolicyType
NOX::Nln::LinearSystemBase::getPreconditionerPolicy(bool advanceReuseCounter)
{
  (void)advanceReuseCounter;

  return ::NOX::Epetra::LinearSystem::PRPT_REBUILD;
}

bool NOX::Nln::LinearSystemBase::applyRightPreconditioning(bool useTranspose,
    Teuchos::ParameterList& params, const ::NOX::Epetra::Vector& input,
    ::NOX::Epetra::Vector& result) const
{
  (void)useTranspose;
  (void)params;

  if (&result != &input) result = input;

  return true;
}

FOUR_C_NAMESPACE_CLOSE
