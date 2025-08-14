// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_FSI_NOX_LINEARSYSTEM_GCR_HPP
#define FOUR_C_FSI_NOX_LINEARSYSTEM_GCR_HPP

#include "4C_config.hpp"

#include "4C_solver_nonlin_nox_linearsystem_base.hpp"
#include "4C_solver_nonlin_nox_scaling.hpp"
#include "4C_utils_parameter_list.fwd.hpp"

#include <NOX_Common.H>
#include <NOX_Epetra_Interface_Jacobian.H>
#include <NOX_Epetra_Interface_Required.H>
#include <NOX_Epetra_LinearSystem.H>
#include <NOX_Epetra_Vector.H>
#include <NOX_Utils.H>
#include <Teuchos_Time.hpp>

#include <memory>
#include <vector>

FOUR_C_NAMESPACE_OPEN

namespace NOX
{
  namespace FSI
  {
    /// generalized conjugate residual linear system solver for NOX
    /*!
      No preconditioner supported.
      Reuse of internal search directions supported.
     */
    class LinearSystemGCR : public NOX::Nln::LinearSystemBase
    {
     protected:
      //! List of types of epetra objects that can be used for the Jacobian and/or Preconditioner.
      enum OperatorType
      {
        //! An Epetra_Operator derived object.
        EpetraOperator,
        //! An Epetra_RowMatrix derived object.
        EpetraRowMatrix,
        //! An Epetra_VbrMatrix object.
        EpetraVbrMatrix,
        //! An Epetra_CrsMatrix object.
        EpetraCrsMatrix
      };

     public:
      //! Constructor with a user supplied Jacobian Operator.
      LinearSystemGCR(Teuchos::ParameterList& printParams,
          Teuchos::ParameterList& linearSolverParams,
          const Teuchos::RCP<::NOX::Epetra::Interface::Required>& iReq,
          const Teuchos::RCP<::NOX::Epetra::Interface::Jacobian>& iJac,
          const Teuchos::RCP<Epetra_Operator>& J, const ::NOX::Epetra::Vector& cloneVector,
          const Teuchos::RCP<::NOX::Epetra::Scaling> scalingObject = Teuchos::null);

      //! Reset the linear solver parameters.
      virtual void reset(Teuchos::ParameterList& linearSolverParams);

      /*!
        \brief Applies Jacobian to the given input vector and puts the answer in the result.

        Computes
        \f[ v = J u, \f]
        where \f$J\f$ is the Jacobian, \f$u\f$ is the input vector,
        and \f$v\f$ is the result vector.  Returns true if successful.
      */
      bool applyJacobian(
          const ::NOX::Epetra::Vector& input, ::NOX::Epetra::Vector& result) const override;

      /*!
        \brief Applies Jacobian-Transpose to the given input vector and puts the answer in the
        result.

        Computes
        \f[ v = J^T u, \f]
        where \f$J\f$ is the Jacobian, \f$u\f$ is the input vector, and \f$v\f$ is the result
        vector.  Returns true if successful.

      */
      bool applyJacobianTranspose(
          const ::NOX::Epetra::Vector& input, ::NOX::Epetra::Vector& result) const override;

      /*!
        \brief Applies the inverse of the Jacobian matrix to the given
        input vector and puts the answer in result.

        Computes
        \f[ v = J^{-1} u, \f]
        where \f$J\f$ is the Jacobian, \f$u\f$ is the input vector,
        and \f$v\f$ is the result vector.

        The parameter list contains the linear solver options.
      */
      bool applyJacobianInverse(Teuchos::ParameterList& params, const ::NOX::Epetra::Vector& input,
          ::NOX::Epetra::Vector& result) override;

      /*!
        \brief Apply right preconditiong to the given input vector

        Let \f$M\f$ be a right preconditioner for the Jacobian \f$J\f$; in
        other words, \f$M\f$ is a matrix such that
        \f[ JM \approx I. \f]

        Compute
        \f[ u = M^{-1} v, \f]
        where \f$u\f$ is the input vector and \f$v\f$ is the result vector.

        If <em>useTranspose</em> is true, then the transpose of the
        preconditioner is applied:
        \f[ u = {M^{-1}}^T v, \f]
        The transpose preconditioner is currently only required for
        Tensor methods.

        The parameter list contains the linear solver options.
      */
      bool applyRightPreconditioning(bool useTranspose, Teuchos::ParameterList& params,
          const ::NOX::Epetra::Vector& input, ::NOX::Epetra::Vector& result) const override;

      //! Evaluates the Jacobian based on the solution vector x.
      bool computeJacobian(const ::NOX::Epetra::Vector& x) override;

      /*!
        \brief Explicitly constructs a preconditioner based on the solution vector x and the
        parameter list p.

        The user has the option of recomputing the graph when a new preconditioner is created. The
        Group controls the isValid flag for the preconditioner and will control when to call this.
      */
      bool createPreconditioner(const ::NOX::Epetra::Vector& x, Teuchos::ParameterList& p,
          bool recomputeGraph) const override;

      /*!
        \brief Deletes the preconditioner.

        The Group controls the isValid flag for the preconditioner and will control when to call
        this.
      */
      bool destroyPreconditioner() const override;

      /*! \brief Recalculates the preconditioner using an already allocated graph.

      Use this to compute a new preconditioner while using the same
      graph for the preconditioner.  This avoids deleting and
      reallocating the memory required for the preconditioner and
      results in a big speed-up for large-scale jobs.
      */
      bool recomputePreconditioner(const ::NOX::Epetra::Vector& x,
          Teuchos::ParameterList& linearSolverParams) const override;

      /*! \brief  Evaluates the preconditioner policy at the current state.

      NOTE: This can change values between nonlienar iterations.  It is
      not a static value.
      */
      PreconditionerReusePolicyType getPreconditionerPolicy(
          bool advanceReuseCounter = true) override;

      //! Indicates whether a preconditioner has been constructed
      bool isPreconditionerConstructed() const override;

      //! Indicates whether the linear system has a preconditioner
      bool hasPreconditioner() const override;

      //! Return Jacobian operator
      Teuchos::RCP<const Epetra_Operator> getJacobianOperator() const override;

      //! Return Jacobian operator
      Teuchos::RCP<Epetra_Operator> getJacobianOperator() override;

      //! Return preconditioner operator
      /*!
       * Note:  This should only be called if hasPreconditioner() returns true.
       */
      Teuchos::RCP<const Epetra_Operator> getGeneratedPrecOperator() const override;

      //! Return preconditioner operator
      Teuchos::RCP<Epetra_Operator> getGeneratedPrecOperator() override;

      //! Set Jacobian operator for solve
      void setJacobianOperatorForSolve(
          const Teuchos::RCP<const Epetra_Operator>& solveJacOp) override;

      //! Set preconditioner operator for solve
      /*!
       * Note:  This should only be called if hasPreconditioner() returns true.
       */
      void setPrecOperatorForSolve(const Teuchos::RCP<const Epetra_Operator>& solvePrecOp) override;

      //! Returns the type of operator that is passed into the group constructors.
      /*! Uses dynamic casting to identify the underlying object type. */
      virtual OperatorType get_operator_type(const Epetra_Operator& o);

     protected:
      /// generalized conjugate residual solver
      /*!
        Implemented following GMRESR without inner GMRES.

        H. A. Van der Vorst and C. Vuik, GMRESR: A family of nested
        GMRES methods, Num. Lin. Alg. Appl., 1 (1994),
        pp. 369--386. http://citeseer.ist.psu.edu/vandervorst91gmresr.html
       */
      int solve_gcr(
          const ::NOX::Epetra::Vector& b, ::NOX::Epetra::Vector& x, int& maxit, double& tol);

      /// GMRES solver
      /*!
        Implementation taken from netlib.

        Barrett, R. and Berry, M. and Chan, T. and Demmel, J. and
        Donato, J. and Dongarra J. and Eijkhout, V. and Pozo, R. and
        Romine Ch. and van der Vorst, H.: Templates for the Solution of
        Linear Systems: Building Blocks for Iterative Methods, SIAM
        (1993)
      */
      int solve_gmres(const ::NOX::Epetra::Vector& b, ::NOX::Epetra::Vector& x, int& max_iter,
          double& tol, int m);

      /// helper for GMRES
      void apply_plane_rotation(double& dx, double& dy, double& cs, double& sn);

      /// helper for GMRES
      void generate_plane_rotation(double& dx, double& dy, double& cs, double& sn);

      virtual void throw_error(const std::string& functionName, const std::string& errorMsg) const;

     protected:
      //! Printing Utilities object
      ::NOX::Utils utils;

      //! Reference to the user supplied Jacobian interface functions
      Teuchos::RCP<::NOX::Epetra::Interface::Jacobian> jacInterfacePtr;

      //! Type of operator for the Jacobian.
      OperatorType jacType;

      //! Pointer to the Jacobian operator.
      mutable Teuchos::RCP<Epetra_Operator> jacPtr;

      //! Scaling object supplied by the user
      Teuchos::RCP<::NOX::Epetra::Scaling> scaling;

      //! An extra temporary vector, only allocated if needed.
      mutable std::shared_ptr<::NOX::Epetra::Vector> tmpVectorPtr;

      mutable double conditionNumberEstimate;

      //! If set to true, solver information is printed to the "Output" sublist of the "Linear
      //! Solver" list.
      bool outputSolveDetails;

      //! Zero out the initial guess for linear solves performed through applyJacobianInverse calls
      //! (i.e. zero out the result vector before the linear solve).
      bool zeroInitialGuess;

      //! Stores the parameter "Compute Scaling Manually".
      bool manualScaling;

      //! Teuchos_Time object
      Teuchos::Time timer;

      //! Total time spent in applyJacobianInverse (sec.).
      mutable double timeApplyJacbianInverse;

      std::vector<std::shared_ptr<::NOX::Epetra::Vector>> u_;

      std::vector<std::shared_ptr<::NOX::Epetra::Vector>> c_;
    };

  }  // namespace FSI
}  // namespace NOX

FOUR_C_NAMESPACE_CLOSE

#endif
