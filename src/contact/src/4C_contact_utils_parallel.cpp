// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_contact_utils_parallel.hpp"

#include "4C_contact_input.hpp"
#include "4C_global_data.hpp"
#include "4C_inpar_mortar.hpp"
#include "4C_inpar_structure.hpp"
#include "4C_legacy_enum_definitions_problem_type.hpp"

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_StandardParameterEntryValidators.hpp>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool CONTACT::Utils::use_safe_redistribute_and_ghosting(const Teuchos::ParameterList& contactParams)
{
  /* Limit the use of the new safe "redistribute & ghosting" branch to our core contact
   * capabilities. If your case of interest is missing here, feel free to migrate your scenario to
   * the new safe branch.
   */
  bool use_safe_ghosting_branch = false;
  {
    const Teuchos::ParameterList& sdyn = Global::Problem::instance()->structural_dynamic_params();
    const auto intstrat =
        Teuchos::getIntegralValue<Inpar::Solid::IntegrationStrategy>(sdyn, "INT_STRATEGY");

    if (intstrat == Inpar::Solid::int_old)
    {
      /* Enable new safe ghosting only for interface discretization type "mortar"
       *
       * There's a conflict with create_volume_ghosting(). This affects all Nitsche-type algorithms
       * and also classical Penalty with Gauss-point-to-segment (GPTS).
       *
       * In theory, penalty with GPTS should work just fine, because it should never need a volume
       * ghosting. However, penalty with GPTS is implemented in the NitscheStrategy, which
       * always requires volume ghosting.
       *
       * Other cases require volume ghosting as well and, thus, have to stick to the old code
       * branch. Everything porous media related has to stick to the old code branch as well.
       */
      if (Teuchos::getIntegralValue<Inpar::Mortar::AlgorithmType>(contactParams, "ALGORITHM") ==
              Inpar::Mortar::algorithm_mortar &&
          (Global::Problem::instance()->get_problem_type() != Core::ProblemType::poroelast &&
              Global::Problem::instance()->get_problem_type() != Core::ProblemType::poroscatra))
        use_safe_ghosting_branch = true;
    }
    else
    {
      // Use old code path, if structure uses the new time integration.
    }
  }

  return use_safe_ghosting_branch;
}

FOUR_C_NAMESPACE_CLOSE
