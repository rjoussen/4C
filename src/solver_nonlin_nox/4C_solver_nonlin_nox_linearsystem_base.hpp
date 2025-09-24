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

    /**
     * \brief Base class for NOX linear systems.
     *
     * Currently, this class is a temporary base to proceed with a smooth step-by-step clean up of
     * the Epetra related interface methods. Most probably, it will be removed as soon as Epetra is
     * removed from 4C.
     */
    class LinearSystemBase : public ::NOX::Epetra::LinearSystem
    {
      /**
       * \brief Get native Epetra scaling object.
       *
       * This method throws an exception and is a temporary mock.
       */
      Teuchos::RCP<::NOX::Epetra::Scaling> getScaling() final;

      /**
       * \brief Get native Epetra scaling object.
       *
       * This method throws an exception and is a temporary mock.
       */
      void resetScaling(const Teuchos::RCP<::NOX::Epetra::Scaling>& s) final;

      /**
       * \brief Set preconditioner operator for solve.
       *
       * This method does nothing and is a temporary mock.
       */
      void setPrecOperatorForSolve(const Teuchos::RCP<const Epetra_Operator>& solvePrecOp) final;

      /**
       * \brief Indicates whether a preconditioner has been constructed.
       *
       * This method does nothing and is a temporary mock.
       */
      bool isPreconditionerConstructed() const final;

      /**
       * \brief Indicates whether the linear system has a preconditioner.
       *
       * This method does nothing and is a temporary mock.
       */
      bool hasPreconditioner() const final;

      /**
       * \brief Return preconditioner operator.
       *
       * This method returns null and is a temporary mock.
       */
      Teuchos::RCP<const Epetra_Operator> getGeneratedPrecOperator() const final;

      /**
       * \brief Return preconditioner operator.
       *
       * This method returns null and is a temporary mock.
       */
      Teuchos::RCP<Epetra_Operator> getGeneratedPrecOperator() final;

      /**
       * \brief Explicitly constructs a preconditioner based on the solution vector x and the
       * parameter list p.
       *
       * This method does nothing and is a temporary mock.
       */
      bool createPreconditioner(const ::NOX::Epetra::Vector& x, Teuchos::ParameterList& p,
          bool recomputeGraph) const final;

      /**
       * \brief Deletes the preconditioner.
       *
       * This method does nothing and is a temporary mock.
       */
      bool destroyPreconditioner() const final;

      /**
       * \brief Recalculates the preconditioner using an already allocated graph.
       *
       * This method does nothing and is a temporary mock.
       */
      bool recomputePreconditioner(
          const ::NOX::Epetra::Vector& x, Teuchos::ParameterList& linearSolverParams) const final;

      /**
       * \brief Evaluates the preconditioner policy at the current state.
       *
       * This method returns a value and is a temporary mock.
       */
      ::NOX::Epetra::LinearSystem::PreconditionerReusePolicyType getPreconditionerPolicy(
          bool advanceReuseCounter = true) final;

      /**
       * \brief Apply right preconditiong to the given input vector.
       *
       * This method does nothing and is a temporary mock.
       */
      bool applyRightPreconditioning(bool useTranspose, Teuchos::ParameterList& params,
          const ::NOX::Epetra::Vector& input, ::NOX::Epetra::Vector& result) const final;

      /**
       * \brief Set Jacobian operator for solve.
       *
       * This method does nothing and is a temporary mock.
       */
      void setJacobianOperatorForSolve(const Teuchos::RCP<const Epetra_Operator>& solveJacOp) final;
    };
  }  // namespace Nln
}  // namespace NOX

FOUR_C_NAMESPACE_CLOSE

#endif
