// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_SOLVER_NONLIN_NOX_LINEARSYSTEM_BASE_HPP
#define FOUR_C_SOLVER_NONLIN_NOX_LINEARSYSTEM_BASE_HPP

#include "4C_config.hpp"

#include <NOX_Epetra_LinearSystem.H>

// forward declaration
namespace NOX
{
  namespace Epetra
  {
    class Scaling;
  }  // namespace Epetra
}  // namespace NOX

FOUR_C_NAMESPACE_OPEN

namespace NOX
{
  namespace Nln
  {
    class Scaling;

    /*
     * \brief Base class for NOX linear systems.
     *
     * This mocks the scaling related methods of the NOX::Epetra::LinearSystem class. Most probably,
     * it will be removed as soon as Epetra is removed from 4C.
     */
    class LinearSystemBase : public ::NOX::Epetra::LinearSystem
    {
      Teuchos::RCP<::NOX::Epetra::Scaling> getScaling() override;

      void resetScaling(const Teuchos::RCP<::NOX::Epetra::Scaling>& s) override;
    };
  }  // namespace Nln
}  // namespace NOX

FOUR_C_NAMESPACE_CLOSE

#endif
