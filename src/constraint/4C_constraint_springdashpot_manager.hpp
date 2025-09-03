// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_CONSTRAINT_SPRINGDASHPOT_MANAGER_HPP
#define FOUR_C_CONSTRAINT_SPRINGDASHPOT_MANAGER_HPP

#include "4C_config.hpp"

#include "4C_linalg_vector.hpp"
#include "4C_utils_parameter_list.fwd.hpp"

#include <memory>

FOUR_C_NAMESPACE_OPEN

namespace Core::FE
{
  class Discretization;
}  // namespace Core::FE

namespace Core::LinAlg
{
  class SparseMatrix;
}  // namespace Core::LinAlg

namespace Core::IO
{
  class DiscretizationWriter;
  class DiscretizationReader;
}  // namespace Core::IO

namespace Constraints
{
  class SpringDashpot;

  class SpringDashpotManager
  {
   public:
    /*!
      \brief Constructor
    */
    SpringDashpotManager(std::shared_ptr<Core::FE::Discretization> dis);

    /*!
     \brief Return if there are spring dashpots
    */
    bool have_spring_dashpot() const { return havespringdashpot_; };

    //! add contribution of spring dashpot BC to residual vector and stiffness matrix
    void stiffness_and_internal_forces(std::shared_ptr<Core::LinAlg::SparseMatrix> stiff,
        std::shared_ptr<Core::LinAlg::Vector<double>> fint,
        std::shared_ptr<Core::LinAlg::Vector<double>> disn, Core::LinAlg::Vector<double>& veln,
        Teuchos::ParameterList parlist);

    //! update for each new time step
    void update();

    //! output of gap, normal, and nodal stiffness
    void output(Core::IO::DiscretizationWriter& output, Core::FE::Discretization& discret,
        Core::LinAlg::Vector<double>& disp);

    //! output of prestressing offset for restart
    void output_restart(std::shared_ptr<Core::IO::DiscretizationWriter> output_restart,
        Core::FE::Discretization& discret, Core::LinAlg::Vector<double>& disp);

    /*!
     \brief Read restart information
    */
    void read_restart(Core::IO::DiscretizationReader& reader, const double& time);

    //! reset spring after having done a MULF prestressing update (mhv 12/2015)
    void reset_prestress(Core::LinAlg::Vector<double>& disold);

   private:
    std::shared_ptr<Core::FE::Discretization> actdisc_;    ///< standard discretization
    std::vector<std::shared_ptr<SpringDashpot>> springs_;  ///< all spring dashpot instances

    bool havespringdashpot_;  ///< are there any spring dashpot BCs at all?
    int n_conds_;             ///< number of spring dashpot conditions
  };  // class
}  // namespace Constraints
FOUR_C_NAMESPACE_CLOSE

#endif
