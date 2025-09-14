// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_STRUCTURE_NEW_NLN_SOLVER_FACTORY_HPP
#define FOUR_C_STRUCTURE_NEW_NLN_SOLVER_FACTORY_HPP

#include "4C_config.hpp"

#include "4C_inpar_structure.hpp"

#include <Teuchos_ParameterList.hpp>

FOUR_C_NAMESPACE_OPEN

namespace Solid
{
  class Integrator;
  namespace TimeInt
  {
    class BaseDataGlobalState;
    class BaseDataSDyn;
    class Base;
    class NoxInterface;
  }  // namespace TimeInt

  namespace Nln::SOLVER
  {
    class Generic;

    /*! Non-member function, which relates to the Solid::Nln::SOLVER::Factory class
     *  Please call this method from outside! */
    std::shared_ptr<Solid::Nln::SOLVER::Generic> build_nln_solver(
        const enum Inpar::Solid::NonlinSolTech& nlnSolType,
        const std::shared_ptr<Solid::TimeInt::BaseDataGlobalState>& gstate,
        const std::shared_ptr<Solid::TimeInt::BaseDataSDyn>& sdyn,
        const std::shared_ptr<Solid::TimeInt::NoxInterface>& noxinterface,
        const std::shared_ptr<Solid::Integrator>& integrator,
        const std::shared_ptr<const Solid::TimeInt::Base>& timint);
  }  // namespace Nln::SOLVER
}  // namespace Solid

FOUR_C_NAMESPACE_CLOSE

#endif
