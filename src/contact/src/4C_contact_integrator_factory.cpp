// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_contact_integrator_factory.hpp"

#include "4C_contact_ehl_integrator.hpp"
#include "4C_contact_nitsche_integrator.hpp"
#include "4C_contact_nitsche_integrator_fpi.hpp"
#include "4C_contact_nitsche_integrator_fsi.hpp"
#include "4C_contact_nitsche_integrator_poro.hpp"
#include "4C_contact_nitsche_integrator_ssi.hpp"
#include "4C_contact_nitsche_integrator_ssi_elch.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
std::shared_ptr<CONTACT::Integrator> CONTACT::INTEGRATOR::Factory::build_integrator(
    const CONTACT::SolvingStrategy& sol_type, Teuchos::ParameterList& mortar_params,
    const Core::FE::CellType& slave_type, MPI_Comm comm) const
{
  std::shared_ptr<CONTACT::Integrator> integrator = nullptr;
  switch (sol_type)
  {
    case CONTACT::SolvingStrategy::nitsche:
    {
      if (mortar_params.get<CONTACT::Problemtype>("PROBTYPE") == CONTACT::Problemtype::ssi)
      {
        integrator =
            std::make_shared<CONTACT::IntegratorNitscheSsi>(mortar_params, slave_type, comm);
      }
      else if (mortar_params.get<CONTACT::Problemtype>("PROBTYPE") ==
               CONTACT::Problemtype::ssi_elch)
      {
        integrator =
            std::make_shared<CONTACT::IntegratorNitscheSsiElch>(mortar_params, slave_type, comm);
      }
      else if (mortar_params.get<CONTACT::Problemtype>("PROBTYPE") ==
                   CONTACT::Problemtype::poroelast ||
               mortar_params.get<CONTACT::Problemtype>("PROBTYPE") ==
                   CONTACT::Problemtype::poroscatra)
      {
        integrator =
            std::make_shared<CONTACT::IntegratorNitschePoro>(mortar_params, slave_type, comm);
      }
      else if (mortar_params.get<CONTACT::Problemtype>("PROBTYPE") == CONTACT::Problemtype::fsi)
      {
        integrator =
            std::make_shared<CONTACT::IntegratorNitscheFsi>(mortar_params, slave_type, comm);
      }
      else if (mortar_params.get<CONTACT::Problemtype>("PROBTYPE") == CONTACT::Problemtype::fpi)
      {
        integrator =
            std::make_shared<CONTACT::IntegratorNitscheFpi>(mortar_params, slave_type, comm);
      }
      else
      {
        integrator = std::make_shared<CONTACT::IntegratorNitsche>(mortar_params, slave_type, comm);
      }
      break;
    }
    case CONTACT::SolvingStrategy::penalty:
    case CONTACT::SolvingStrategy::multiscale:
    {
      if (Teuchos::getIntegralValue<Inpar::Mortar::AlgorithmType>(mortar_params, "ALGORITHM") ==
          Inpar::Mortar::algorithm_gpts)
        integrator = std::make_shared<CONTACT::IntegratorNitsche>(mortar_params, slave_type, comm);
      else
        integrator = std::make_shared<CONTACT::Integrator>(mortar_params, slave_type, comm);
      break;
    }
    case CONTACT::SolvingStrategy::lagmult:
    case CONTACT::SolvingStrategy::uzawa:
    {
      integrator = std::make_shared<CONTACT::Integrator>(mortar_params, slave_type, comm);
      break;
    }
    case CONTACT::SolvingStrategy::ehl:
    {
      integrator = std::make_shared<CONTACT::IntegratorEhl>(mortar_params, slave_type, comm);

      break;
    }
    default:
    {
      FOUR_C_THROW("Unsupported solving strategy! (stype = {})",
          CONTACT::solving_strategy_to_string(sol_type));
    }
  }  // end switch

  return integrator;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
std::shared_ptr<CONTACT::Integrator> CONTACT::INTEGRATOR::build_integrator(
    const CONTACT::SolvingStrategy& sol_type, Teuchos::ParameterList& mortar_params,
    const Core::FE::CellType& slave_type, MPI_Comm comm)
{
  Factory factory;
  return factory.build_integrator(sol_type, mortar_params, slave_type, comm);
}

FOUR_C_NAMESPACE_CLOSE
