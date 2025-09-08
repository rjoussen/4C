// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_BEAMCONTACT_INPUT_HPP
#define FOUR_C_BEAMCONTACT_INPUT_HPP

#include "4C_config.hpp"

#include "4C_fem_condition.hpp"
#include "4C_io_input_spec.hpp"
#include "4C_utils_parameter_list.fwd.hpp"

#include <vector>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/


// Forward declaration.
namespace Core::Conditions
{
  class ConditionDefinition;
}

namespace BeamContact
{
  /// Type of employed solving strategy
  /// (this enum represents the input file parameter BEAMS_STRATEGY)
  enum Strategy
  {
    bstr_none,    ///< no beam contact
    bstr_penalty  ///< penalty method
  };

  /// Type of employed solving strategy
  /// (this enum represents the input file parameter BEAMS_STRATEGY)
  enum Modelevaluator
  {
    bstr_old,      ///<  use old beaminteraction model
    bstr_standard  ///<  use new beamcontact submodel
  };

  /// Application of a smoothed tangent field
  /// (this enum represents the input file parameter BEAMS_SMOOTHING)
  enum Smoothing
  {
    bsm_none,  ///< no smoothing
    bsm_cpp    ///< smoothing only for the closest point projection; element evaluation without
               ///< smoothing"
  };

  /// Application of a smoothed tangent field
  /// (this enum represents the input file parameter BEAMS_SMOOTHING)
  enum PenaltyLaw
  {
    pl_lp,     ///< linear penalty law
    pl_qp,     ///< quadratic penalty law
    pl_lnqp,   ///< linear penalty law with quadratic regularization for negative gaps
    pl_lpqp,   ///< linear penalty law with quadratic regularization for positive gaps
    pl_lpcp,   ///< linear penalty law with cubic regularization for positive gaps
    pl_lpdqp,  ///< linear penalty law with double quadratic regularization for positive gaps
    pl_lpep    ///< linear penalty law with exponential regularization for positive gaps
  };

  /// Beam Contact Octree and Bounding Box Type
  /// (this enum represents the input file parameter BEAMS_OCTREEBBOX)
  enum OctreeType
  {
    boct_none,  ///< no bounding box -> no octree
    boct_aabb,  ///< axis aligned bounding boxes
    boct_cobb,  ///< cylindrical oriented bounding boxes
    boct_spbb   ///< spherical bounding boxes
  };

  /// set the beam contact parameters
  std::vector<Core::IO::InputSpec> set_valid_parameters();

  /**
   * \brief Set beam beam-to-beam specific conditions.
   */
  void set_valid_conditions(std::vector<Core::Conditions::ConditionDefinition>& condlist);
}  // namespace BeamContact

/*----------------------------------------------------------------------*/
FOUR_C_NAMESPACE_CLOSE

#endif
