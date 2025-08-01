// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_structure_new_nln_solver_factory.hpp"

#include "4C_structure_new_nln_solver_nox.hpp"
#include "4C_structure_new_nln_solver_utils.hpp"
#include "4C_structure_new_timint_base.hpp"
#include "4C_utils_enum.hpp"

#include <NOX_Utils.H>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Solid::Nln::SOLVER::Factory::Factory()
{
  // empty
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
std::shared_ptr<Solid::Nln::SOLVER::Generic> Solid::Nln::SOLVER::Factory::build_nln_solver(
    const enum Inpar::Solid::NonlinSolTech& nlnSolType,
    const std::shared_ptr<Solid::TimeInt::BaseDataGlobalState>& gstate,
    const std::shared_ptr<Solid::TimeInt::BaseDataSDyn>& sdyn,
    const std::shared_ptr<Solid::TimeInt::NoxInterface>& noxinterface,
    const std::shared_ptr<Solid::Integrator>& integrator,
    const std::shared_ptr<const Solid::TimeInt::Base>& timint) const
{
  Teuchos::ParameterList nox_params;

  switch (nlnSolType)
  {
    case Inpar::Solid::soltech_newtonfull:
    {
      // ---------------------------------------------------------------------------
      // Set-up the full Newton method
      // ---------------------------------------------------------------------------
      // Non-linear solver
      nox_params.set("Nonlinear Solver", "Line Search Based");

      // Direction
      Teuchos::ParameterList& pdir = nox_params.sublist("Direction");
      pdir.set("Method", "Newton");

      // Newton
      Teuchos::ParameterList& pnewton = pdir.sublist("Newton");

      // Adapt forcing term
      std::string forcing_term_method = sdyn->get_nox_params()
                                            .sublist("Direction")
                                            .sublist("Newton")
                                            .get<std::string>("Forcing Term Method");
      pnewton.set("Forcing Term Method", forcing_term_method);

      // Line Search
      Teuchos::ParameterList& plinesearch = nox_params.sublist("Line Search");
      plinesearch.set("Method", "Full Step");
    }
    break;

    case Inpar::Solid::soltech_ptc:
    {
      // ---------------------------------------------------------------------------
      // Set-up the pseudo transient method
      // ---------------------------------------------------------------------------
      // Non-linear solver
      nox_params.set("Nonlinear Solver", "Pseudo Transient");

      // Direction
      Teuchos::ParameterList& pdir = nox_params.sublist("Direction");
      pdir.set("Method", "Newton");

      /* The following parameters create a NOX::Nln::Solver::Nox
       * solver which is equivalent to the old 4C implementation.
       *
       * If you are keen on using the new features, please use the corresponding
       * input section "STRUCT NOX/Pseudo Transient" in your input file. */
      Teuchos::ParameterList& pptc = nox_params.sublist("Pseudo Transient");

      pptc.set<double>("deltaInit", sdyn->get_initial_ptc_pseudo_time_step());
      pptc.set<double>("deltaMax", std::numeric_limits<double>::max());
      pptc.set<double>("deltaMin", 0.0);
      pptc.set<int>("Maximum Number of Pseudo-Transient Iterations", (sdyn->get_iter_max() + 1));
      pptc.set<int>("Maximum Number of Pseudo-Transient Iterations", 51);
      pptc.set<std::string>("Time Step Control", "SER");
      pptc.set<double>("SER_alpha", 1.0);
      pptc.set<double>("ScalingFactor", 1.0);
      pptc.set<std::string>("Norm Type for TSC", "Max Norm");
      pptc.set<std::string>("Build scaling operator", "every timestep");
      pptc.set<std::string>("Scaling Type", "Identity");
    }
    break;

    case Inpar::Solid::soltech_singlestep:
    {
      nox_params.set("Nonlinear Solver", "Single Step");

      // Set the printing parameters in the "Printing" sublist
      Teuchos::ParameterList& printParams = nox_params.sublist("Printing");
      printParams.set("Output Precision", 3);
      printParams.set("Output Processor", 0);
      printParams.set("Output Information",
          ::NOX::Utils::OuterIteration + ::NOX::Utils::OuterIterationStatusTest +
              ::NOX::Utils::InnerIteration + ::NOX::Utils::Details + ::NOX::Utils::Warning);

      // Set NOX solving parameters: "Update Jacobian" is essential to compute the Newton direction,
      // i.e. x=A^-1 b. For details behind how it's handled, see NOX_Solver_SingleStep.C
      Teuchos::ParameterList& sssParams = nox_params.sublist("Single Step Solver");
      sssParams.set("Update Jacobian", true);
      sssParams.set("Print Norms", true);
      sssParams.set("Compute Relative Norm", true);

      // Sublist for linear solver for the single step
      Teuchos::ParameterList& lsParams = sssParams.sublist("Linear Solver");
      lsParams.set("Adaptive Control", false);
      lsParams.set("Number of Nonlinear Iterations", 1);
    }
    break;

    default:
      break;
  }

  // ---------------------------------------------------------------------------
  // STATUS TEST
  // ---------------------------------------------------------------------------
  /* This is only necessary for the special case, that you use no xml-file for
   * the definition of your convergence tests, but you use the input file instead.
   */
  if (not is_xml_status_test_file(sdyn->get_nox_params().sublist("Status Test")))
  {
    std::set<enum NOX::Nln::StatusTest::QuantityType> qtypes;
    Solid::Nln::SOLVER::create_quantity_types(qtypes, *sdyn);

    // remove the unsupported quantity of status test:
    // EAS is removed since it is an element
    // quantity and not a nodal dof
    qtypes.erase(NOX::Nln::StatusTest::quantity_eas);

    Solid::Nln::SOLVER::set_status_test_params(nox_params.sublist("Status Test"), *sdyn, qtypes);
  }

  return std::make_shared<Solid::Nln::SOLVER::Nox>(
      nox_params, gstate, sdyn, noxinterface, integrator, timint);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
std::shared_ptr<Solid::Nln::SOLVER::Generic> Solid::Nln::SOLVER::build_nln_solver(
    const enum Inpar::Solid::NonlinSolTech& nlnSolType,
    const std::shared_ptr<Solid::TimeInt::BaseDataGlobalState>& gstate,
    const std::shared_ptr<Solid::TimeInt::BaseDataSDyn>& sdyn,
    const std::shared_ptr<Solid::TimeInt::NoxInterface>& noxinterface,
    const std::shared_ptr<Solid::Integrator>& integrator,
    const std::shared_ptr<const Solid::TimeInt::Base>& timint)
{
  Factory factory;
  return factory.build_nln_solver(nlnSolType, gstate, sdyn, noxinterface, integrator, timint);
}

FOUR_C_NAMESPACE_CLOSE
