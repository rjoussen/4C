// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_GEOMETRY_PAIR_INPUT_HPP
#define FOUR_C_GEOMETRY_PAIR_INPUT_HPP


#include "4C_config.hpp"

#include "4C_fem_general_utils_integration.hpp"

// Forward declaration.
#include "4C_io_input_spec.hpp"
#include "4C_utils_parameter_list.hpp"

#include <map>

FOUR_C_NAMESPACE_OPEN

namespace GeometryPair
{
  /**
   * Method to be used for line to 3D pairs.
   */
  enum class LineTo3DStrategy
  {
    none,                                                  ///< Default value.
    gauss_point_projection_without_boundary_segmentation,  ///< Project the Gauss points.
    gauss_point_projection_boundary_segmentation,  ///< Project the Gauss points, segment lines
                                                   ///< that poke out of the volume.
    gauss_point_projection_cross_section,  ///< Project the Gauss points on the beam surface.
    segmentation                           ///< Segment the line with the volume.
  };

  /**
   * Flag on what to do if not all Gauss points of a segment project
   */
  enum class NotAllGaussPointsProjectValidAction
  {
    fail,     // Raise an error
    warning,  // Print a warning message to the terminal and continue to evaluate the pair
  };

  /**
   * Method to be used for normals in surface pairs.
   */
  enum class SurfaceNormals
  {
    standard,        ///< Use normals on the actual surface.
    extended_volume  ///< Use normals from extended volume projection.
  };

  /**
   * \brief Map number of gauss points to 1D gauss rule
   */
  Core::FE::GaussRule1D int_to_gauss_rule1_d(const int n_gauss_points);

  /**
   * \brief Set valid input parameters for line to 3D geometry pairs.
   *
   * This function adds all available input parameters to a given section. Call this function with
   * the parameter list of the input section you want the line to 3D parameters to be added to.
   *
   * @param (out) Parameter list to add the line to 3D parameters to.
   */
  void valid_parameters_line_to3_d(std::vector<Core::IO::InputSpec>& list);

  /**
   * \brief Set valid input parameters for line to surface geometry pairs.
   *
   * This function adds all available input parameters for line to surface pairs, that are not yet
   * added in the SetValidParametersLineTo3D function.
   *
   * @param (out) Parameter list to add the line to surface parameters to.
   */
  void valid_parameters_line_to_surface(std::vector<Core::IO::InputSpec>& list);

}  // namespace GeometryPair


FOUR_C_NAMESPACE_CLOSE

#endif
