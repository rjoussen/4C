// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_CONTACT_INTEGRATOR_FACTORY_HPP
#define FOUR_C_CONTACT_INTEGRATOR_FACTORY_HPP

#include "4C_config.hpp"

#include "4C_contact_input.hpp"
#include "4C_fem_general_element.hpp"

// forward declaration
#include "4C_utils_parameter_list.fwd.hpp"

FOUR_C_NAMESPACE_OPEN

namespace CONTACT
{
  class Integrator;
  namespace INTEGRATOR
  {
    class Factory
    {
     public:
      /*! \brief Build the desired contact integrator
       *

       * */
      std::shared_ptr<CONTACT::Integrator> build_integrator(
          const CONTACT::SolvingStrategy& sol_type, Teuchos::ParameterList& mortar_params,
          const Core::FE::CellType& slave_type, MPI_Comm comm) const;
    };  // class Factory

    // non-member function, please call this one from outside!
    std::shared_ptr<CONTACT::Integrator> build_integrator(const CONTACT::SolvingStrategy& sol_type,
        Teuchos::ParameterList& mortar_params, const Core::FE::CellType& slave_type, MPI_Comm comm);
  }  // namespace INTEGRATOR
}  // namespace CONTACT


FOUR_C_NAMESPACE_CLOSE

#endif
