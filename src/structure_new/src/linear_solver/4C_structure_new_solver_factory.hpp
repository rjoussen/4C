// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_STRUCTURE_NEW_SOLVER_FACTORY_HPP
#define FOUR_C_STRUCTURE_NEW_SOLVER_FACTORY_HPP

#include "4C_config.hpp"

#include "4C_utils_parameter_list.fwd.hpp"

#include <map>
#include <memory>
#include <set>

FOUR_C_NAMESPACE_OPEN

// forward declarations...
namespace Inpar
{
  namespace Solid
  {
    enum ModelType : int;
  }
}  // namespace Inpar
namespace CONTACT
{
  enum class SolvingStrategy;
  enum class SystemType;
}  // namespace CONTACT
namespace Core::LinAlg
{
  class Solver;
}
namespace Core::FE
{
  class Discretization;
}  // namespace Core::FE
namespace Solid
{
  namespace SOLVER
  {
    /*! \brief Factory to build the desired linear solver std::map corresponding
     *  to the active model types.
     *
     *  */
    class Factory
    {
     private:
      using LinSolMap =
          std::map<enum Inpar::Solid::ModelType, std::shared_ptr<Core::LinAlg::Solver>>;

     public:
      //! constructor
      Factory() = default;

      //! destructor
      virtual ~Factory() = default;

      //! build the desired linear solvers
      std::shared_ptr<LinSolMap> build_lin_solvers(
          const std::set<enum Inpar::Solid::ModelType>& modeltypes,
          const Teuchos::ParameterList& sdyn, Core::FE::Discretization& actdis) const;

      //! create the meshtying/contact linear solver
      static std::shared_ptr<Core::LinAlg::Solver> build_meshtying_contact_lin_solver(
          Core::FE::Discretization& actdis, enum CONTACT::SolvingStrategy sol_type,
          enum CONTACT::SystemType sys_type, const int lin_solver_id);

     private:
      //! create the structural linear solver (should be called by default)
      std::shared_ptr<Core::LinAlg::Solver> build_structure_lin_solver(
          const Teuchos::ParameterList& sdyn, Core::FE::Discretization& actdis) const;

      //! create the meshtying/contact linear solver
      std::shared_ptr<Core::LinAlg::Solver> build_meshtying_contact_lin_solver(
          Core::FE::Discretization& actdis) const;

      //! create the Lagrange/penalty enforced constraint linear solver
      std::shared_ptr<Core::LinAlg::Solver> build_lag_pen_constraint_lin_solver(
          const Teuchos::ParameterList& sdyn, Core::FE::Discretization& actdis) const;

      //! create the Windkessel linear solver
      std::shared_ptr<Core::LinAlg::Solver> build_cardiovascular0_d_lin_solver(
          const Teuchos::ParameterList& sdyn, Core::FE::Discretization& actdis) const;

    };  // class Factory

    /*! Non-member function, which relates to the Solid::SOLVER::Factory class
     *  Please call this method from outside! */
    std::shared_ptr<std::map<enum Inpar::Solid::ModelType, std::shared_ptr<Core::LinAlg::Solver>>>
    build_lin_solvers(const std::set<enum Inpar::Solid::ModelType>& modeltypes,
        const Teuchos::ParameterList& sdyn, Core::FE::Discretization& actdis);
  }  // namespace SOLVER
}  // namespace Solid

FOUR_C_NAMESPACE_CLOSE

#endif
