// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_STI_PARTITIONED_HPP
#define FOUR_C_STI_PARTITIONED_HPP

#include "4C_config.hpp"

#include "4C_sti_algorithm.hpp"
#include "4C_sti_input.hpp"

FOUR_C_NAMESPACE_OPEN

namespace STI
{
  //! partitioned coupling algorithm for scatra-thermo interaction
  class Partitioned : public Algorithm
  {
   public:
    //! constructor
    explicit Partitioned(MPI_Comm comm,        //! communicator
        const Teuchos::ParameterList& stidyn,  //! parameter list for scatra-thermo interaction
        const Teuchos::ParameterList&
            scatradyn,  //! scalar transport parameter list for scatra and thermo fields
        const Teuchos::ParameterList&
            solverparams_scatra,  //! solver parameter list for scatra field
        const Teuchos::ParameterList&
            solverparams_thermo  //! solver parameter list for thermo field
    );

   private:
    //! convergence check for iterative staggered TSI solver
    bool exit_outer_coupling() const;

    //! evaluate time step using outer coupling iteration
    void solve() override;

    //! evaluate time step using one-way coupling iteration
    void solve_one_way() const;

    //! evaluate time step using two-way coupling iteration
    void solve_two_way();

    //! type of coupling between scatra and thermo fields
    const STI::CouplingType couplingtype_;

    //! maximum value of Aitken relaxation parameter
    const double omegamax_;
  };  // class Partitioned : public Algorithm
}  // namespace STI
FOUR_C_NAMESPACE_CLOSE

#endif
