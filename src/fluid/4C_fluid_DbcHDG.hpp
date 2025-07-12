// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_FLUID_DBCHDG_HPP
#define FOUR_C_FLUID_DBCHDG_HPP


#include "4C_config.hpp"

#include "4C_fem_discretization.hpp"
#include "4C_fem_discretization_utils.hpp"
#include "4C_linalg_utils_densematrix_communication.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Discret::Utils
{
  class Dbc;
  // class DbcInfo;
}  // namespace Discret::Utils

namespace FLD
{
  namespace Utils
  {
    /** \brief Specialized Dbc evaluation class for HDG discretizations
     *
     *  */
    class DbcHdgFluid : public Core::FE::Dbc
    {
     public:
      /// constructor
      DbcHdgFluid() {};

     protected:
      void read_dirichlet_condition(const Teuchos::ParameterList& params,
          const Core::FE::Discretization& discret, const Core::Conditions::Condition& cond,
          double time, DbcInfo& info, const std::shared_ptr<std::set<int>>* dbcgids,
          int hierarchical_order) const override;

      void read_dirichlet_condition(const Teuchos::ParameterList& params,
          const Core::FE::DiscretizationFaces& discret, const Core::Conditions::Condition& cond,
          double time, DbcInfo& info, const std::shared_ptr<std::set<int>>* dbcgids,
          int hierarchical_order) const;

      void do_dirichlet_condition(const Teuchos::ParameterList& params,
          const Core::FE::Discretization& discret, const Core::Conditions::Condition& cond,
          double time, const std::shared_ptr<Core::LinAlg::Vector<double>>* systemvectors,
          const Core::LinAlg::Vector<int>& toggle,
          const std::shared_ptr<std::set<int>>* dbcgids) const override;

      void do_dirichlet_condition(const Teuchos::ParameterList& params,
          const Core::FE::DiscretizationFaces& discret, const Core::Conditions::Condition& cond,
          double time, const std::shared_ptr<Core::LinAlg::Vector<double>>* systemvectors,
          const Core::LinAlg::Vector<int>& toggle) const;
    };  // class DbcHDG_Fluid
  }  // namespace Utils

}  // namespace FLD

FOUR_C_NAMESPACE_CLOSE

#endif
