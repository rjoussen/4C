// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_solver_nonlin_nox_direction_newton.hpp"

#include "4C_solver_nonlin_nox_group.hpp"

#include <NOX_GlobalData.H>
#include <NOX_Utils.H>

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::Nln::Direction::Newton::Newton(
    const Teuchos::RCP<::NOX::GlobalData>& gd, Teuchos::ParameterList& p)
    : ::NOX::Direction::Newton(gd, p), utils_(gd->getUtils())
{
  // empty constructor
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool NOX::Nln::Direction::Newton::compute(::NOX::Abstract::Vector& dir,
    ::NOX::Abstract::Group& group, const ::NOX::Solver::Generic& solver)
{
  ::NOX::Abstract::Group::ReturnType status;

  // dynamic cast of the nox_abstract_group
  NOX::Nln::Group* nlnSoln = dynamic_cast<NOX::Nln::Group*>(&group);

  if (nlnSoln == nullptr)
  {
    throw_error("compute", "dynamic_cast to nox_nln_group failed!");
  }

  // Compute F and Jacobian at current solution at once.
  status = nlnSoln->compute_f_and_jacobian();
  if (status != ::NOX::Abstract::Group::Ok)
    throw_error("compute", "Unable to compute F and/or Jacobian at once");

  // ------------------------------------------------
  // call base class version
  // ------------------------------------------------
  // NOTE: The base class calls first computeF() and in a second attempt
  // computeJacobian(). In many cases this is uneconomical and instead we do it
  // here at once. In this way the base class compute function calls perform a
  // direct return, because the isValid flags are already set to true!
  try
  {
    return ::NOX::Direction::Newton::compute(dir, group, solver);
  }
  catch (const char* e)
  {
    if (utils_->isPrintType(::NOX::Utils::Warning))
    {
      utils_->out() << e;
    }
  }
  return false;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::Nln::Direction::Newton::throw_error(
    const std::string& functionName, const std::string& errorMsg)
{
  if (utils_->isPrintType(::NOX::Utils::Error))
    utils_->err() << "NOX::Nln::Direction::Newton::" << functionName << " - " << errorMsg
                  << std::endl;
  throw "NOX Error";
}

FOUR_C_NAMESPACE_CLOSE
