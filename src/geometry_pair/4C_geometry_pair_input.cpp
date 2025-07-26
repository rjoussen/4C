// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_geometry_pair_input.hpp"

#include "4C_io_input_spec_builders.hpp"

FOUR_C_NAMESPACE_OPEN


void GeometryPair::set_valid_parameters_line_to3_d(std::vector<Core::IO::InputSpec>& list)
{
  using namespace Core::IO::InputSpecBuilders;

  // Add the input parameters for line to 3D coupling.

  // Segmentation strategy.
  list.push_back(parameter<LineTo3DStrategy>(
      "GEOMETRY_PAIR_STRATEGY", {.description = "Type of employed segmentation strategy",
                                    .default_value = LineTo3DStrategy::segmentation}));

  // Number of search points for segmentation.
  list.push_back(parameter<int>("GEOMETRY_PAIR_SEGMENTATION_SEARCH_POINTS",
      {.description = "Number of search points for segmentation", .default_value = 6}));

  // What to do if not all Gauss points of a segment project valid
  list.push_back(parameter<NotAllGaussPointsProjectValidAction>(
      "GEOMETRY_PAIR_SEGMENTATION_NOT_ALL_GAUSS_POINTS_PROJECT_VALID_ACTION",
      {.description = "What to do if not all Gauss points of a segment project valid",
          .default_value = NotAllGaussPointsProjectValidAction::fail}));

  // Number of integration points on the line.
  list.push_back(parameter<int>("GAUSS_POINTS",
      {.description = "Number of Gauss Points for the integral evaluations", .default_value = 6}));

  // Number of integration along the circumference in cross section coupling.
  list.push_back(parameter<int>("INTEGRATION_POINTS_CIRCUMFERENCE",
      {.description = "Number of Integration points along the circumferential direction of the "
                      "beam. This is parameter is only used in beam to cylinder meshtying. No "
                      "gauss integration is used along the circumferential direction, equally "
                      "spaced integration points are used.",
          .default_value = 6}));
}

void GeometryPair::set_valid_parameters_line_to_surface(std::vector<Core::IO::InputSpec>& list)
{
  // Add the input parameters for line to surface coupling.

  // Add the surface normal option.
  list.push_back(Core::IO::InputSpecBuilders::parameter<GeometryPair::SurfaceNormals>(
      "GEOMETRY_PAIR_SURFACE_NORMALS",
      {.description = "How the surface normals are evaluated",
          .default_value = GeometryPair::SurfaceNormals::standard}));
}

Core::FE::GaussRule1D GeometryPair::int_to_gauss_rule1_d(const int n_gauss_points)
{
  switch (n_gauss_points)
  {
    case 1:
      return Core::FE::GaussRule1D::line_1point;
    case 2:
      return Core::FE::GaussRule1D::line_2point;
    case 3:
      return Core::FE::GaussRule1D::line_3point;
    case 4:
      return Core::FE::GaussRule1D::line_4point;
    case 5:
      return Core::FE::GaussRule1D::line_5point;
    case 6:
      return Core::FE::GaussRule1D::line_6point;
    case 7:
      return Core::FE::GaussRule1D::line_7point;
    case 8:
      return Core::FE::GaussRule1D::line_8point;
    case 9:
      return Core::FE::GaussRule1D::line_9point;
    case 10:
      return Core::FE::GaussRule1D::line_10point;
    case 20:
      return Core::FE::GaussRule1D::line_20point;
    case 32:
      return Core::FE::GaussRule1D::line_32point;
    case 50:
      return Core::FE::GaussRule1D::line_50point;
    default:
    {
      FOUR_C_THROW("No Gauss rule defined for {} points", n_gauss_points);
      return Core::FE::GaussRule1D::undefined;
    }
  }
};

FOUR_C_NAMESPACE_CLOSE
