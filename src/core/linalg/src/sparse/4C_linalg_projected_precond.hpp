// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_LINALG_PROJECTED_PRECOND_HPP
#define FOUR_C_LINALG_PROJECTED_PRECOND_HPP


#include "4C_config.hpp"

#include "4C_comm_mpi_utils.hpp"
#include "4C_linear_solver_method_projector.hpp"

#include <Epetra_Operator.h>

#include <memory>

// forward declarations

FOUR_C_NAMESPACE_OPEN

namespace Core::LinAlg
{
  class KrylovProjector;


  /*!

  A common interface for ifpack, ml and simpler preconditioners.
  This interface allows a modification of the vector returned
  by the ApplyInverse call, which is necessary to do a solution on
  a Krylov space krylovized to certain (for example rigid body)
  modes.

  The linalg preconditioner interface class holds a pointer (rcp)
  to the actual preconditioner. All methods implemented to support
  the Epetra_Operator interface just call the corresponding functions
  of the actual preconditioner.

  Only the ApplyInverse method is modified and performs the
  projection if desired.

  See linalg_projected_operator.H for related docu and code.

  */
  class LinalgPrecondOperator : public Epetra_Operator
  {
   public:
    /*!
    \brief Standard Constructor

    */
    LinalgPrecondOperator(std::shared_ptr<Epetra_Operator> precond, bool project,
        std::shared_ptr<Core::LinAlg::LinearSystemProjector> projector);



    //! @name Attribute set methods required to support the Epetra_Operator interface

    int SetUseTranspose(bool UseTranspose) override
    {
      return (precond_->SetUseTranspose(UseTranspose));
    }
    //! @}

    //! @name Mathematical functions required to support the Epetra_Operator interface (pass
    //! through)
    int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const override
    {
      return (precond_->Apply(X, Y));
    }

    double NormInf() const override { return (precond_->NormInf()); }
    //! @}

    //! @name Attribute access functions required to support the Epetra_Operator interface
    const char* Label() const override { return (precond_->Label()); }

    bool UseTranspose() const override { return (precond_->UseTranspose()); }

    bool HasNormInf() const override { return (precond_->HasNormInf()); }

    const Epetra_Comm& Comm() const override { return (precond_->Comm()); }

    const Epetra_Map& OperatorDomainMap() const override { return (precond_->OperatorDomainMap()); }

    const Epetra_Map& OperatorRangeMap() const override { return (precond_->OperatorRangeMap()); }

    //! @}

    //! @name Mathematical functions required to support the Epetra_Operator interface (modified)
    /*
      (Modified) ApplyInverse call

      This method calls ApplyInverse on the actual preconditioner and, the
      solution is krylovized against a set of weight vectors provided in a
      multivector.

      This is done using a projector P defined by

                                      T
                                     x * w
                          P  x = x - ------ c
                                      T
                                     w * c

      w is the vector of weights, c a vector of ones (in the dofs under
      consideration) corresponding to the matrix kernel.

      The system we are solving with this procedure is not Au=b for u (since A
      might be singular), but we are solving

                          / T \         T
                         | P A | P u = P b ,
                          \   /

      for the projection of the solution Pu, i.e. in the preconditioned case


                                                            -+
             / T   \     /      -1 \          T              |
            | P * A | * |  P * M    | * xi = P  * b          |
             \     /     \         /                         |
                                                  -1         |
                                         x = P * M  * xi     |
                                                            -+


      Hence, P is always associated with the apply inverse call of the
      preconditioner (the right bracket) and always called after the call
      to ApplyInverse.


      Properties of P are:

      1) c defines the kernel of P, i.e. P projects out the matrix kernel

                            T
                           c * w
                P c = c - ------- c = c - c = 0
                            T
                           w * c

      2) The space spanned by P x is krylov to the weight vector

                         /      T      \              T
       T   /   \     T  |      x * w    |    T       x * w     T       T       T
      w * | P x | = w * | x - ------- c | = w * x - ------- * w * c = w * x - w * x = 0
           \   /        |       T       |             T
                         \     w * c   /             w * c


      This modified Apply call is for singular matrices A when c is
      a vector defining A's nullspace. The preceding projection
      operator ensures
                              |           |
                             -+-         -+-T
                    A u = A u     where u    * c =0,

      even if A*c != 0 (for numerical inaccuracies during the computation
      of A)

      See the following article for further reading:

      @article{1055401,
       author = {Bochev,, Pavel and Lehoucq,, R. B.},
       title = {On the Finite Element Solution of the Pure Neumann Problem},
       journal = {SIAM Rev.},
       volume = {47},
       number = {1},
       year = {2005},
       issn = {0036-1445},
       pages = {50--66},
       doi = {http://dx.doi.org/10.1137/S0036144503426074},
       publisher = {Society for Industrial and Applied Mathematics},
       address = {Philadelphia, PA, USA},
       }

    */
    int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const override;


   private:
    //! flag whether to do a projection or just pass through
    bool project_;

    //! the actual preconditioner
    std::shared_ptr<Epetra_Operator> precond_;

    //! projector
    std::shared_ptr<Core::LinAlg::LinearSystemProjector> projector_;
  };

}  // namespace Core::LinAlg

FOUR_C_NAMESPACE_CLOSE

#endif
