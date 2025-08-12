// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_beaminteraction_contact_runtime_visualization_output_params.hpp"

#include "4C_global_data.hpp"
#include "4C_utils_exceptions.hpp"

FOUR_C_NAMESPACE_OPEN

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
BeamInteraction::BeamContactRuntimeVisualizationOutputParams::
    BeamContactRuntimeVisualizationOutputParams(const double restart_time)
    : isinit_(false),
      issetup_(false),
      visualization_parameters_(Core::IO::visualization_parameters_factory(
          Global::Problem::instance()->io_params().sublist("RUNTIME VTK OUTPUT"),
          *Global::Problem::instance()->output_control_file(), restart_time)),
      output_interval_steps_(-1),
      output_every_iteration_(false),
      output_forces_(false),
      output_gaps_(false)
{
  // empty constructor
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void BeamInteraction::BeamContactRuntimeVisualizationOutputParams::init()
{
  issetup_ = false;
  // empty for now

  isinit_ = true;
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void BeamInteraction::BeamContactRuntimeVisualizationOutputParams::setup()
{
  throw_error_if_not_init();

  // Teuchos parameter list for beam contact
  const Teuchos::ParameterList& beam_contact_visualization_output_paramslist =
      Global::Problem::instance()->beam_contact_params().sublist("RUNTIME VTK OUTPUT");

  /****************************************************************************/
  // get and check required parameters
  /****************************************************************************/
  output_interval_steps_ = beam_contact_visualization_output_paramslist.get<int>("INTERVAL_STEPS");

  output_every_iteration_ =
      beam_contact_visualization_output_paramslist.get<bool>("EVERY_ITERATION");
  visualization_parameters_.every_iteration_ = output_every_iteration_;

  /****************************************************************************/
  output_forces_ = beam_contact_visualization_output_paramslist.get<bool>("CONTACT_FORCES");

  /****************************************************************************/
  output_gaps_ = beam_contact_visualization_output_paramslist.get<bool>("GAPS");

  output_angles_ = beam_contact_visualization_output_paramslist.get<bool>("CONTACT_ANGLE");
  output_types_ = beam_contact_visualization_output_paramslist.get<bool>("CONTACT_TYPE");


  issetup_ = true;
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void BeamInteraction::BeamContactRuntimeVisualizationOutputParams::
    throw_error_if_not_init_and_setup() const
{
  if (!is_init() or !is_setup()) FOUR_C_THROW("Call init() and setup() first!");
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void BeamInteraction::BeamContactRuntimeVisualizationOutputParams::throw_error_if_not_init() const
{
  if (!is_init()) FOUR_C_THROW("init() has not been called, yet!");
}

FOUR_C_NAMESPACE_CLOSE
