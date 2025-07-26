// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_beaminteraction_beam_to_solid_surface_meshtying_params.hpp"

#include "4C_beaminteraction_beam_to_solid_surface_visualization_output_params.hpp"
#include "4C_geometry_pair_input.hpp"
#include "4C_global_data.hpp"

#include <Teuchos_StandardParameterEntryValidators.hpp>

FOUR_C_NAMESPACE_OPEN


/**
 *
 */
BeamInteraction::BeamToSolidSurfaceMeshtyingParams::BeamToSolidSurfaceMeshtyingParams()
    : BeamToSolidParamsBase(),
      coupling_type_(Inpar::BeamToSolid::BeamToSolidSurfaceCoupling::none),
      output_params_ptr_(nullptr),
      rotational_coupling_penalty_parameter_(-1.0),
      rotational_coupling_triad_construction_(
          Inpar::BeamToSolid::BeamToSolidSurfaceRotationCoupling::none)
{
  // Empty Constructor.
}


/**
 *
 */
void BeamInteraction::BeamToSolidSurfaceMeshtyingParams::init()
{
  // Teuchos parameter list for beam contact
  const Teuchos::ParameterList& beam_to_solid_surface_meshtying_params_list =
      Global::Problem::instance()->beam_interaction_params().sublist(
          "BEAM TO SOLID SURFACE MESHTYING");

  // Set the common beam-to-solid parameters.
  set_base_params(beam_to_solid_surface_meshtying_params_list);

  // Get parameters form input file.
  {
    // Type of coupling evaluation to be used.
    coupling_type_ = Teuchos::getIntegralValue<Inpar::BeamToSolid::BeamToSolidSurfaceCoupling>(
        beam_to_solid_surface_meshtying_params_list, "COUPLING_TYPE");

    // Parameters for rotational coupling.
    rotational_coupling_ =
        beam_to_solid_surface_meshtying_params_list.get<bool>("ROTATIONAL_COUPLING");
    rotational_coupling_penalty_parameter_ =
        beam_to_solid_surface_meshtying_params_list.get<double>(
            "ROTATIONAL_COUPLING_PENALTY_PARAMETER");
    rotational_coupling_triad_construction_ =
        Teuchos::getIntegralValue<Inpar::BeamToSolid::BeamToSolidSurfaceRotationCoupling>(
            beam_to_solid_surface_meshtying_params_list, "ROTATIONAL_COUPLING_SURFACE_TRIAD");

    if (rotational_coupling_)
    {
      switch (coupling_type_)
      {
        case Inpar::BeamToSolid::BeamToSolidSurfaceCoupling::consistent_fad:
        case Inpar::BeamToSolid::BeamToSolidSurfaceCoupling::displacement_fad:
        case Inpar::BeamToSolid::BeamToSolidSurfaceCoupling::
            reference_configuration_forced_to_zero_fad:
          break;
        default:
          FOUR_C_THROW(
              "Beam-to-solid surface coupling with rotational coupling is only implemented in "
              "combination with the \"fad\" variants of surface coupling.");
      }
    }
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
std::shared_ptr<BeamInteraction::BeamToSolidSurfaceVisualizationOutputParams>
BeamInteraction::BeamToSolidSurfaceMeshtyingParams::get_visualization_output_params_ptr()
{
  return output_params_ptr_;
};

FOUR_C_NAMESPACE_CLOSE
