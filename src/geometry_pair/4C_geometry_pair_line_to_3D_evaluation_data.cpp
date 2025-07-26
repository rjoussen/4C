// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_geometry_pair_line_to_3D_evaluation_data.hpp"

#include "4C_utils_parameter_list.hpp"

FOUR_C_NAMESPACE_OPEN

/**
 *
 */
GeometryPair::LineTo3DEvaluationData::LineTo3DEvaluationData(
    const Teuchos::ParameterList& input_parameter_list)
    : GeometryEvaluationDataBase(input_parameter_list),
      strategy_(GeometryPair::LineTo3DStrategy::none),
      gauss_rule_(Core::FE::GaussRule1D::undefined),
      integration_points_circumference_(-1),
      gauss_point_projection_tracker_(),
      n_search_points_(0),
      not_all_gauss_points_project_valid_action_(
          GeometryPair::NotAllGaussPointsProjectValidAction::fail),
      segment_tracker_()
{
  // Get parameters from the input file.
  {
    strategy_ = Teuchos::getIntegralValue<GeometryPair::LineTo3DStrategy>(
        input_parameter_list, "GEOMETRY_PAIR_STRATEGY");

    n_search_points_ = input_parameter_list.get<int>("GEOMETRY_PAIR_SEGMENTATION_SEARCH_POINTS");
    not_all_gauss_points_project_valid_action_ =
        Teuchos::getIntegralValue<GeometryPair::NotAllGaussPointsProjectValidAction>(
            input_parameter_list,
            "GEOMETRY_PAIR_SEGMENTATION_NOT_ALL_GAUSS_POINTS_PROJECT_VALID_ACTION");

    gauss_rule_ = GeometryPair::int_to_gauss_rule1_d(input_parameter_list.get<int>("GAUSS_POINTS"));

    integration_points_circumference_ =
        input_parameter_list.get<int>("INTEGRATION_POINTS_CIRCUMFERENCE");
  }

  // Initialize evaluation data structures.
  clear();
}

/**
 *
 */
void GeometryPair::LineTo3DEvaluationData::clear()
{
  // Call reset on the base method.
  GeometryEvaluationDataBase::clear();

  // Initialize evaluation data structures.
  {
    // Tracker for gauss point projection method.
    gauss_point_projection_tracker_.clear();

    // Segment tracker for segmentation.
    segment_tracker_.clear();
  }
}

/**
 *
 */
void GeometryPair::LineTo3DEvaluationData::reset_tracker()
{
  for (auto& data : gauss_point_projection_tracker_)
    std::fill(data.second.begin(), data.second.end(), false);

  for (auto& data : segment_tracker_) data.second.clear();
}

FOUR_C_NAMESPACE_CLOSE
