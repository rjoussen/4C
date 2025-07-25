// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_beaminteraction_beam_to_solid_surface_contact_params.hpp"

#include "4C_beaminteraction_beam_to_solid_surface_visualization_output_params.hpp"
#include "4C_geometry_pair_input.hpp"
#include "4C_global_data.hpp"

#include <Teuchos_StandardParameterEntryValidators.hpp>

FOUR_C_NAMESPACE_OPEN


/**
 *
 */
BeamInteraction::BeamToSolidSurfaceContactParams::BeamToSolidSurfaceContactParams()
    : BeamToSolidParamsBase(),
      contact_type_(Inpar::BeamToSolid::BeamToSolidSurfaceContact::none),
      penalty_law_(Inpar::BeamToSolid::BeamToSolidSurfaceContactPenaltyLaw::none),
      penalty_parameter_g0_(0.0),
      mortar_contact_configuration_(
          Inpar::BeamToSolid::BeamToSolidSurfaceContactMortarDefinedIn::none),
      output_params_ptr_(nullptr)
{
  // Empty Constructor.
}


/**
 *
 */
void BeamInteraction::BeamToSolidSurfaceContactParams::init()
{
  // Teuchos parameter list for beam contact
  const Teuchos::ParameterList& beam_to_solid_contact_params_list =
      Global::Problem::instance()->beam_interaction_params().sublist(
          "BEAM TO SOLID SURFACE CONTACT");

  // Set the common beam-to-solid parameters.
  set_base_params(beam_to_solid_contact_params_list);

  // Get parameters form input file.
  {
    contact_type_ = Teuchos::getIntegralValue<Inpar::BeamToSolid::BeamToSolidSurfaceContact>(
        beam_to_solid_contact_params_list, "CONTACT_TYPE");

    penalty_law_ =
        Teuchos::getIntegralValue<Inpar::BeamToSolid::BeamToSolidSurfaceContactPenaltyLaw>(
            beam_to_solid_contact_params_list, "PENALTY_LAW");

    penalty_parameter_g0_ = beam_to_solid_contact_params_list.get<double>("PENALTY_PARAMETER_G0");

    mortar_contact_configuration_ =
        Teuchos::getIntegralValue<Inpar::BeamToSolid::BeamToSolidSurfaceContactMortarDefinedIn>(
            beam_to_solid_contact_params_list, "MORTAR_CONTACT_DEFINED_IN");
  }

  // Setup the output parameter object.
  {
    output_params_ptr_ = std::make_shared<BeamToSolidSurfaceVisualizationOutputParams>();
    output_params_ptr_->init();
    output_params_ptr_->setup();
  }

  isinit_ = true;
}


/**
 *
 */
int BeamInteraction::BeamToSolidSurfaceContactParams::get_fad_order() const

{
  switch (get_contact_type())
  {
    case Inpar::BeamToSolid::BeamToSolidSurfaceContact::gap_variation:
      return 1;
      break;
    case Inpar::BeamToSolid::BeamToSolidSurfaceContact::potential:
      return 2;
      break;
    default:
      FOUR_C_THROW("Got unexpected contact type.");
      return 0;
  }
}

FOUR_C_NAMESPACE_CLOSE
