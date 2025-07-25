// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_geometry_pair_line_projection.hpp"

#include "4C_geometry_pair_element.hpp"
#include "4C_geometry_pair_element_evaluation_functions.hpp"
#include "4C_geometry_pair_line_to_3D_evaluation_data.hpp"
#include "4C_geometry_pair_line_to_surface_evaluation_data.hpp"
#include "4C_geometry_pair_line_to_surface_gauss_point_projection.hpp"
#include "4C_geometry_pair_line_to_surface_segmentation.hpp"
#include "4C_geometry_pair_line_to_volume_gauss_point_projection.hpp"
#include "4C_geometry_pair_line_to_volume_segmentation.hpp"
#include "4C_geometry_pair_scalar_types.hpp"
#include "4C_geometry_pair_utility_functions.hpp"
#include "4C_utils_fad.hpp"

FOUR_C_NAMESPACE_OPEN

/**
 *
 */
template <typename PairType>
void GeometryPair::LineTo3DBase<PairType>::project_point_on_line_to_other(const PairType* pair,
    const ElementData<line, scalar_type>& element_data_line,
    const ElementData<other, scalar_type>& element_data_other, const scalar_type& eta,
    Core::LinAlg::Matrix<3, 1, scalar_type>& xi, ProjectionResult& projection_result)
{
  // Get the point on the line.
  Core::LinAlg::Matrix<3, 1, scalar_type> r_line;
  GeometryPair::evaluate_position<line>(eta, element_data_line, r_line);

  // Project the point to the solid.
  pair->project_point_to_other(r_line, element_data_other, xi, projection_result);
}

/**
 *
 */
template <typename PairType>
void GeometryPair::LineTo3DBase<PairType>::project_points_on_line_to_other(const PairType* pair,
    const ElementData<line, scalar_type>& element_data_line,
    const ElementData<other, scalar_type>& element_data_other,
    std::vector<ProjectionPoint1DTo3D<scalar_type>>& projection_points,
    unsigned int& n_projections_valid, unsigned int& n_projections)
{
  // Initialize counters.
  n_projections_valid = 0;
  n_projections = 0;

  // Loop over points and check if they project to this other geometry.
  for (auto& point : projection_points)
  {
    // Project the point.
    project_point_on_line_to_other(pair, element_data_line, element_data_other, point.get_eta(),
        point.get_xi(), point.get_projection_result());

    // Update the counters.
    if (point.get_projection_result() == ProjectionResult::projection_found_valid)
    {
      n_projections_valid++;
      n_projections++;
    }
    if (point.get_projection_result() == ProjectionResult::projection_found_not_valid)
      n_projections++;
  }
}

/**
 *
 */
template <typename PairType>
void GeometryPair::LineTo3DBase<PairType>::project_points_on_line_to_other(const PairType* pair,
    const ElementData<line, scalar_type>& element_data_line,
    const ElementData<other, scalar_type>& element_data_other,
    std::vector<ProjectionPoint1DTo3D<scalar_type>>& projection_points,
    unsigned int& n_projections_valid)
{
  // Initialize dummy variable.
  unsigned int n_projections_dummy;

  // Project the points.
  project_points_on_line_to_other(pair, element_data_line, element_data_other, projection_points,
      n_projections_valid, n_projections_dummy);
}

/**
 *
 */
template <typename PairType>
void GeometryPair::LineTo3DBase<PairType>::project_gauss_points_on_segment_to_other(
    const PairType* pair, const ElementData<line, scalar_type>& element_data_line,
    const ElementData<other, scalar_type>& element_data_other, LineSegment<scalar_type>& segment)
{
  const auto& evaluation_data = *(pair->get_evaluation_data());
  const Core::FE::IntegrationPoints1D& gauss_points = evaluation_data.get_gauss_points();

  // Set up the vector with the projection points.
  std::vector<ProjectionPoint1DTo3D<scalar_type>>& projection_points =
      segment.get_projection_points();
  projection_points.clear();
  projection_points.reserve(gauss_points.nquad);
  Core::LinAlg::Matrix<3, 1, scalar_type> xi_start;
  StartValues<other::geometry_type_>::set(xi_start);
  for (unsigned int i = 0; i < (unsigned int)gauss_points.nquad; i++)
  {
    scalar_type eta = segment.get_eta_a() + (segment.get_eta_b() - segment.get_eta_a()) * 0.5 *
                                                (gauss_points.qxg[i][0] + 1.);
    projection_points.push_back(
        ProjectionPoint1DTo3D<scalar_type>(eta, xi_start, gauss_points.qwgt[i]));
  }

  // Project the Gauss points to the other geometry.
  unsigned int n_valid_projections;
  unsigned int n_projections;
  project_points_on_line_to_other(pair, element_data_line, element_data_other, projection_points,
      n_valid_projections, n_projections);

  // Check if a warning or an error should be output
  const bool all_valid = n_valid_projections == static_cast<unsigned int>(gauss_points.nquad);
  const bool all_projected = n_projections == static_cast<unsigned int>(gauss_points.nquad);
  const bool is_warning = all_projected and (!all_valid) and
                          evaluation_data.get_not_all_gauss_points_project_valid_action() ==
                              GeometryPair::NotAllGaussPointsProjectValidAction::warning;
  const bool is_error =
      !all_projected or
      (!all_valid and evaluation_data.get_not_all_gauss_points_project_valid_action() !=
                          GeometryPair::NotAllGaussPointsProjectValidAction::warning);
  if (is_warning or is_error)
  {
    // Add detailed output that allows for a reconstruction of the failed projection
    std::stringstream error_message;

    // Get the geometry information of the line and other geometry
    print_pair_information(error_message, pair, element_data_line, element_data_other);

    // Print a projection point
    auto print_projection_point = [&error_message](const auto& projection_point)
    {
      error_message << "\n  line parameter coordinate: "
                    << Core::FADUtils::cast_to_double(projection_point.get_eta());
      error_message << "\n  other parameter coordinate: ";
      Core::FADUtils::cast_to_double(projection_point.get_xi()).print(error_message);
      error_message << "  projection result: " << (int)projection_point.get_projection_result();
    };

    // Print the segment information
    error_message << "\nSegment start point:";
    print_projection_point(segment.get_start_point());
    error_message << "\nSegment end point:";
    print_projection_point(segment.get_end_point());

    // Print the projection result for each Gauss point on the segment
    for (unsigned int i_point = 0; i_point < projection_points.size(); i_point++)
    {
      error_message << "\nGauss point " << i_point;
      print_projection_point(projection_points[i_point]);
    }

    // Depending on the input file, print a warning or fail
    if (is_warning)
    {
      std::cout << error_message.str();
    }
    else
    {
      FOUR_C_THROW(
          "{}\n\nError when projecting the Gauss points. There are {} Gauss points, {} "
          "could be projected and {} of the projected ones are valid.",
          error_message.str(), gauss_points.nquad, n_projections, n_valid_projections);
    }
  }
}


/**
 *
 */
template <typename PairType>
void GeometryPair::LineTo3DBase<PairType>::intersect_line_with_other(const PairType* pair,
    const ElementData<line, scalar_type>& element_data_line,
    const ElementData<other, scalar_type>& element_data_other,
    std::vector<ProjectionPoint1DTo3D<scalar_type>>& intersection_points)
{
  // Set default values for the parameter coordinates.
  scalar_type eta_start;
  Core::LinAlg::Matrix<3, 1, scalar_type> xi_start;
  StartValues<line::geometry_type_>::set(eta_start);
  StartValues<other::geometry_type_>::set(xi_start);

  // Call the intersect function.
  pair->intersect_line_with_other(
      element_data_line, element_data_other, intersection_points, eta_start, xi_start);
}

/**
 *
 */
template <typename PairType>
void GeometryPair::LineTo3DGaussPointProjection<PairType>::pre_evaluate(const PairType* pair,
    const ElementData<line, scalar_type>& element_data_line,
    const ElementData<other, scalar_type>& element_data_other,
    std::vector<LineSegment<scalar_type>>& segments)
{
  // Get the Gauss point projection tracker for this line element.
  std::vector<bool>& line_projection_tracker = get_line_projection_vector(pair);

  // Gauss rule.
  Core::FE::IntegrationPoints1D gauss_points = pair->get_evaluation_data()->get_gauss_points();

  // Initialize variables for the projection.
  scalar_type eta;
  Core::LinAlg::Matrix<3, 1, scalar_type> xi;
  ProjectionResult projection_result;
  LineSegment<scalar_type> line_segment;
  bool one_projects = false;

  // Loop over Gauss points and check if they project to this other geometry.
  for (unsigned int index_gp = 0; index_gp < line_projection_tracker.size(); index_gp++)
  {
    // Only check points that do not already have a valid projection.
    if (line_projection_tracker[index_gp] == false)
    {
      eta = gauss_points.qxg[index_gp][0];
      LineTo3DBase<PairType>::project_point_on_line_to_other(
          pair, element_data_line, element_data_other, eta, xi, projection_result);
      if (projection_result == ProjectionResult::projection_found_valid)
      {
        // Valid Gauss point was found, add to this segment and set tracking point to true.
        line_segment.add_projection_point(
            ProjectionPoint1DTo3D<scalar_type>(eta, xi, gauss_points.qwgt[index_gp]));
        line_projection_tracker[index_gp] = true;

        one_projects = true;
      }
    }
  }

  if (one_projects)
  {
    // Clear the segment vector and add the found segment for the current line-to-xxx pair.
    segments.clear();
    segments.push_back(line_segment);
  }
}

/**
 *
 */
template <typename PairType>
void GeometryPair::LineTo3DGaussPointProjection<PairType>::evaluate(const PairType* pair,
    const ElementData<line, scalar_type>& element_data_line,
    const ElementData<other, scalar_type>& element_data_other,
    std::vector<LineSegment<scalar_type>>& segments)
{
  // We only check for boundary segmentation if it is needed.
  switch (pair->get_evaluation_data()->get_strategy())
  {
    case GeometryPair::LineTo3DStrategy::gauss_point_projection_without_boundary_segmentation:
      return;
    case GeometryPair::LineTo3DStrategy::gauss_point_projection_boundary_segmentation:
      break;
    default:
      FOUR_C_THROW("Wrong LineTo3DStrategy in Evaluate of Gauss point projection pairs.");
  }

  // Only zero one segments are expected.
  if (segments.size() > 1)
    FOUR_C_THROW(
        "There should be zero or one segments for the Gauss point projection method. The actual "
        "value is {}!",
        segments.size());

  // Check if one point projected in pre_evaluate.
  if (segments.size() == 1 && segments[0].get_number_of_projection_points() > 0)
  {
    // Flag if segmentation is needed.
    bool need_segmentation = false;

    // Check if all Gauss points projected for this line.
    const std::vector<bool>& line_projection_tracker = get_line_projection_vector(pair);
    for (auto const& projects : line_projection_tracker)
      if (!projects) need_segmentation = true;

    if (need_segmentation)
    {
      // Segmentation is needed. First get the intersection points with the other geometry.
      std::vector<ProjectionPoint1DTo3D<scalar_type>> intersection_points;
      LineTo3DBase<PairType>::intersect_line_with_other(
          pair, element_data_line, element_data_other, intersection_points);

      // This algorithm only works if one intersection point was found.
      if (intersection_points.size() != 1)
        FOUR_C_THROW(
            "In the segmentation case we expect exactly one found intersection point. Got: {}!",
            intersection_points.size());

      // Get the limits of the segmented line.
      scalar_type eta_a, eta_b, eta_intersection_point, eta_first_gauss_point;
      eta_intersection_point = intersection_points[0].get_eta();
      eta_first_gauss_point = segments[0].get_projection_points()[0].get_eta();
      if (eta_intersection_point < eta_first_gauss_point)
      {
        eta_a = eta_intersection_point;
        eta_b = 1.;
      }
      else
      {
        eta_a = -1.;
        eta_b = eta_intersection_point;
      }

      // Reproject the Gauss points on the segmented line.
      segments[0] = LineSegment<scalar_type>(eta_a, eta_b);
      LineTo3DBase<PairType>::project_gauss_points_on_segment_to_other(
          pair, element_data_line, element_data_other, segments[0]);
    }
  }
}

/**
 *
 */
template <typename PairType>
std::vector<bool>& GeometryPair::LineTo3DGaussPointProjection<PairType>::get_line_projection_vector(
    const PairType* pair)
{
  // Get the Gauss point projection tracker for this line element.
  int line_element_id = pair->element1()->id();
  std::map<int, std::vector<bool>>& projection_tracker =
      pair->get_evaluation_data()->get_gauss_point_projection_tracker();
  auto find = projection_tracker.find(line_element_id);
  if (find == projection_tracker.end())
    FOUR_C_THROW("Could not find the projection tracker for line id {}.", line_element_id);
  return find->second;
}


/**
 *
 */
template <typename PairType>
void GeometryPair::LineTo3DSegmentation<PairType>::evaluate(const PairType* pair,
    const ElementData<line, scalar_type>& element_data_line,
    const ElementData<other, scalar_type>& element_data_other,
    std::vector<LineSegment<scalar_type>>& segments)
{
  // Only zero segments are expected.
  if (segments.size() > 0)
    FOUR_C_THROW(
        "There should be zero segments for the segmentation method. The actual value is {}!",
        segments.size());

  // Number of search points.
  unsigned int n_search_points = pair->get_evaluation_data()->get_number_of_search_points();

  // Set up vector with projection points for the search points.
  std::vector<ProjectionPoint1DTo3D<scalar_type>> search_points;
  search_points.reserve(n_search_points);
  Core::LinAlg::Matrix<3, 1, scalar_type> xi_start;
  StartValues<other::geometry_type_>::set(xi_start);
  scalar_type eta;
  for (unsigned int i_search_point = 0; i_search_point < n_search_points; i_search_point++)
  {
    eta = -1. + i_search_point * 2. / (n_search_points - 1.);
    search_points.push_back(ProjectionPoint1DTo3D<scalar_type>(eta, xi_start));
  }

  // Project all the search points.
  unsigned int dummy;
  unsigned int n_projections;
  LineTo3DBase<PairType>::project_points_on_line_to_other(
      pair, element_data_line, element_data_other, search_points, dummy, n_projections);

  // If no point could be projected return, as we assume that we wont find a surface projection.
  // This usually happens for higher order elements, more search points can be a solution to
  // that problem.
  if (n_projections != 0)
  {
    // Now we start from every search point that could be projected (not necessary valid) and look
    // for surface intersections of the line with the other geometry. The intersection points are
    // stored in a std::set that will keep them unique and in order.
    std::set<ProjectionPoint1DTo3D<scalar_type>> intersection_points;

    // If the start and/or end points (nodes of the line element) could be projected valid, add them
    // to the intersections set. This helps for some cases where the line just touches the surface
    // of the other geometry.
    if (search_points.front().get_projection_result() == ProjectionResult::projection_found_valid)
      intersection_points.insert(search_points.front());
    if (search_points.back().get_projection_result() == ProjectionResult::projection_found_valid)
      intersection_points.insert(search_points.back());

    // Vector for intersection point search.
    std::vector<ProjectionPoint1DTo3D<scalar_type>> search_intersection_points;

    // Starting from each search point, try to project to all surfaces of the other geometry.
    for (auto const& point : search_points)
    {
      // Only use search points that could be projected.
      if (point.get_projection_result() != ProjectionResult::projection_not_found)
      {
        // Get the intersections with the other geometry.
        pair->intersect_line_with_other(element_data_line, element_data_other,
            search_intersection_points, point.get_eta(), point.get_xi());

        // Add the found intersection points to the set.
        for (auto& found_point : search_intersection_points)
          intersection_points.insert(found_point);
      }
    }

    // In the case of zero and one intersection points, no segmentation is needed. One point only
    // occurs when the line just touches the other geometry. This is the reason why the start and/or
    // end points are added to the set of found intersection.
    if (intersection_points.size() > 1)
    {
      // The intersection points in intersection_points are in order. Now it is checked if a
      // point between the intersection points is inside or outside the other geometry and therefore
      // it can be decided if the segment is part of this pair. By doing so we can avoid
      // complications in cases when a line is exactly between two other geometries.

      // We reuse the vector created for the search points.
      search_points.clear();
      search_points.reserve(intersection_points.size() - 1);

      // Loop over the middle points and check if they are projected valid or not.
      std::set<GeometryPair::LineSegment<double>>& segment_tracker = get_segment_tracking_set(pair);
      ProjectionResult projection_result;
      bool last_segment_active = false;
      ProjectionPoint1DTo3D<scalar_type> segment_start;
      unsigned int counter = 0;
      for (typename std::set<ProjectionPoint1DTo3D<scalar_type>>::iterator set_iterator =
               intersection_points.begin();
          set_iterator != intersection_points.end(); ++set_iterator)
      {
        // Reference to this current point.
        const ProjectionPoint1DTo3D<scalar_type>& start_point = *set_iterator;

        // Get the next intersection point and calculate the projection. This can only be done if
        // the iterator is not on its last iteration.
        if (counter != intersection_points.size() - 1)
        {
          // Reference to the next point.
          const ProjectionPoint1DTo3D<scalar_type>& end_point = *std::next(set_iterator);

          // Get starting points for projection.
          eta = 0.5 * (start_point.get_eta() + end_point.get_eta());
          xi_start = start_point.get_xi();
          xi_start += end_point.get_xi();
          xi_start.scale(0.5);

          // Project and check result.
          LineTo3DBase<PairType>::project_point_on_line_to_other(
              pair, element_data_line, element_data_other, eta, xi_start, projection_result);
        }
        else
        {
          // In case the iterator is on the last intersection point, finish up the active segment.
          projection_result = ProjectionResult::projection_found_not_valid;
        }

        // Check if a new segment is found.
        if (projection_result == ProjectionResult::projection_found_valid)
        {
          // The line part between the current start_point and end_point is a segment for this pair.
          // Start its counter if no segment is currently active.
          if (!last_segment_active)
          {
            last_segment_active = true;
            segment_start = start_point;
          }
        }
        else
        {
          // This point is outside of the other geometry -> if current segment exists finish it.
          if (last_segment_active)
          {
            // Create a segment with double as the scalar type.
            LineSegment<double> new_segment_double(
                Core::FADUtils::cast_to_double(segment_start.get_eta()),
                Core::FADUtils::cast_to_double(start_point.get_eta()));

            // Check if the segment already exists for this line.
            if (segment_tracker.find(new_segment_double) == segment_tracker.end())
            {
              // Add the new segment to this pair and to the evaluation tracker.
              segments.push_back(LineSegment<scalar_type>(segment_start, start_point));
              segment_tracker.insert(new_segment_double);

              // Project the Gauss points on the segment.
              LineTo3DBase<PairType>::project_gauss_points_on_segment_to_other(
                  pair, element_data_line, element_data_other, segments.back());
            }

            // Deactivate the current segment.
            last_segment_active = false;
          }
        }

        // Advance the counter.
        counter++;
      }
    }
  }
}

/**
 *
 */
template <typename PairType>
std::set<GeometryPair::LineSegment<double>>&
GeometryPair::LineTo3DSegmentation<PairType>::get_segment_tracking_set(const PairType* pair)
{
  // Get the segment tracker for this line element.
  int line_element_id = pair->element1()->id();
  std::map<int, std::set<GeometryPair::LineSegment<double>>>& segment_tracker_map =
      pair->get_evaluation_data()->get_segment_tracker();
  auto find = segment_tracker_map.find(line_element_id);
  if (find == segment_tracker_map.end())
    FOUR_C_THROW("Could not find the segment tracker for line id {}.", line_element_id);
  return find->second;
}



/**
 * Explicit template initialization of template class.
 *
 * With the optional arguments the template initiations become quite long and unreadable, therefore
 * they are initialized with preprocessor macros in the following section.
 */
namespace GeometryPair
{
  // Define line-to-volume Gauss point projection pairs.
#define initialize_template_volume_gauss_point(a, b, c) \
  template class LineTo3DGaussPointProjection<          \
      GeometryPairLineToVolumeGaussPointProjection<a, b, c>>;


  initialize_template_volume_gauss_point(double, t_hermite, t_hex8);
  initialize_template_volume_gauss_point(double, t_hermite, t_hex20);
  initialize_template_volume_gauss_point(double, t_hermite, t_hex27);
  initialize_template_volume_gauss_point(double, t_hermite, t_tet4);
  initialize_template_volume_gauss_point(double, t_hermite, t_tet10);
  initialize_template_volume_gauss_point(double, t_hermite, t_nurbs27);
  initialize_template_volume_gauss_point(double, t_hermite, t_wedge6);

  // Define line-to-volume segmentation pairs.
#define initialize_template_volume_segmentation(a, b, c) \
  template class LineTo3DSegmentation<GeometryPairLineToVolumeSegmentation<a, b, c>>;

  initialize_template_volume_segmentation(double, t_hermite, t_hex8);
  initialize_template_volume_segmentation(double, t_hermite, t_hex20);
  initialize_template_volume_segmentation(double, t_hermite, t_hex27);
  initialize_template_volume_segmentation(double, t_hermite, t_tet4);
  initialize_template_volume_segmentation(double, t_hermite, t_tet10);
  initialize_template_volume_segmentation(double, t_hermite, t_nurbs27);
  initialize_template_volume_segmentation(double, t_hermite, t_wedge6);

  // Helper types for the macro initialization. The compiler has troubles inserting the templated
  // typenames into the macros.
  using line_to_surface_patch_scalar_type_fixed_size_1st_order_line2_nurbs9 =
      line_to_surface_patch_scalar_type_fixed_size_1st_order<t_line2, t_nurbs9>;
  using line_to_surface_patch_scalar_type_fixed_size_line2_nurbs9 =
      line_to_surface_patch_scalar_type_fixed_size<t_line2, t_nurbs9>;
  using line_to_surface_patch_scalar_type_fixed_size_1st_order_hermite_nurbs9 =
      line_to_surface_patch_scalar_type_fixed_size_1st_order<t_hermite, t_nurbs9>;
  using line_to_surface_patch_scalar_type_fixed_size_hermite_nurbs9 =
      line_to_surface_patch_scalar_type_fixed_size<t_hermite, t_nurbs9>;

  // Define line-to-surface Gauss point projection pairs.
#define initialize_template_surface_gauss_point(a, b, c) \
  template class LineTo3DGaussPointProjection<           \
      GeometryPairLineToSurfaceGaussPointProjection<a, b, c>>;

  initialize_template_surface_gauss_point(double, t_line2, t_quad4);
  initialize_template_surface_gauss_point(double, t_line2, t_quad8);
  initialize_template_surface_gauss_point(double, t_line2, t_quad9);
  initialize_template_surface_gauss_point(double, t_line2, t_nurbs9);
  initialize_template_surface_gauss_point(double, t_line2, t_tri3);
  initialize_template_surface_gauss_point(double, t_line2, t_tri6);

  initialize_template_surface_gauss_point(
      line_to_surface_patch_scalar_type_1st_order, t_line2, t_quad4);
  initialize_template_surface_gauss_point(
      line_to_surface_patch_scalar_type_1st_order, t_line2, t_quad8);
  initialize_template_surface_gauss_point(
      line_to_surface_patch_scalar_type_1st_order, t_line2, t_quad9);
  initialize_template_surface_gauss_point(
      line_to_surface_patch_scalar_type_fixed_size_1st_order_line2_nurbs9, t_line2, t_nurbs9);
  initialize_template_surface_gauss_point(
      line_to_surface_patch_scalar_type_1st_order, t_line2, t_tri3);
  initialize_template_surface_gauss_point(
      line_to_surface_patch_scalar_type_1st_order, t_line2, t_tri6);

  initialize_template_surface_gauss_point(line_to_surface_patch_scalar_type, t_line2, t_quad4);
  initialize_template_surface_gauss_point(line_to_surface_patch_scalar_type, t_line2, t_quad8);
  initialize_template_surface_gauss_point(line_to_surface_patch_scalar_type, t_line2, t_quad9);
  initialize_template_surface_gauss_point(
      line_to_surface_patch_scalar_type_fixed_size_line2_nurbs9, t_line2, t_nurbs9);
  initialize_template_surface_gauss_point(line_to_surface_patch_scalar_type, t_line2, t_tri3);
  initialize_template_surface_gauss_point(line_to_surface_patch_scalar_type, t_line2, t_tri6);

  initialize_template_surface_gauss_point(double, t_hermite, t_quad4);
  initialize_template_surface_gauss_point(double, t_hermite, t_quad8);
  initialize_template_surface_gauss_point(double, t_hermite, t_quad9);
  initialize_template_surface_gauss_point(double, t_hermite, t_nurbs9);
  initialize_template_surface_gauss_point(double, t_hermite, t_tri3);
  initialize_template_surface_gauss_point(double, t_hermite, t_tri6);

  initialize_template_surface_gauss_point(
      line_to_surface_patch_scalar_type_1st_order, t_hermite, t_quad4);
  initialize_template_surface_gauss_point(
      line_to_surface_patch_scalar_type_1st_order, t_hermite, t_quad8);
  initialize_template_surface_gauss_point(
      line_to_surface_patch_scalar_type_1st_order, t_hermite, t_quad9);
  initialize_template_surface_gauss_point(
      line_to_surface_patch_scalar_type_fixed_size_1st_order_hermite_nurbs9, t_hermite, t_nurbs9);
  initialize_template_surface_gauss_point(
      line_to_surface_patch_scalar_type_1st_order, t_hermite, t_tri3);
  initialize_template_surface_gauss_point(
      line_to_surface_patch_scalar_type_1st_order, t_hermite, t_tri6);

  initialize_template_surface_gauss_point(line_to_surface_patch_scalar_type, t_hermite, t_quad4);
  initialize_template_surface_gauss_point(line_to_surface_patch_scalar_type, t_hermite, t_quad8);
  initialize_template_surface_gauss_point(line_to_surface_patch_scalar_type, t_hermite, t_quad9);
  initialize_template_surface_gauss_point(
      line_to_surface_patch_scalar_type_fixed_size_hermite_nurbs9, t_hermite, t_nurbs9);
  initialize_template_surface_gauss_point(line_to_surface_patch_scalar_type, t_hermite, t_tri3);
  initialize_template_surface_gauss_point(line_to_surface_patch_scalar_type, t_hermite, t_tri6);

  // Define line-to-surface segmentation pairs.
#define initialize_template_surface_segmentation(a, b, c) \
  template class LineTo3DSegmentation<GeometryPairLineToSurfaceSegmentation<a, b, c>>;

  initialize_template_surface_segmentation(double, t_line2, t_quad4);
  initialize_template_surface_segmentation(double, t_line2, t_quad8);
  initialize_template_surface_segmentation(double, t_line2, t_quad9);
  initialize_template_surface_segmentation(double, t_line2, t_nurbs9);
  initialize_template_surface_segmentation(double, t_line2, t_tri3);
  initialize_template_surface_segmentation(double, t_line2, t_tri6);

  initialize_template_surface_segmentation(
      line_to_surface_patch_scalar_type_1st_order, t_line2, t_quad4);
  initialize_template_surface_segmentation(
      line_to_surface_patch_scalar_type_1st_order, t_line2, t_quad8);
  initialize_template_surface_segmentation(
      line_to_surface_patch_scalar_type_1st_order, t_line2, t_quad9);
  initialize_template_surface_segmentation(
      line_to_surface_patch_scalar_type_fixed_size_1st_order_line2_nurbs9, t_line2, t_nurbs9);
  initialize_template_surface_segmentation(
      line_to_surface_patch_scalar_type_1st_order, t_line2, t_tri3);
  initialize_template_surface_segmentation(
      line_to_surface_patch_scalar_type_1st_order, t_line2, t_tri6);

  initialize_template_surface_segmentation(line_to_surface_patch_scalar_type, t_line2, t_quad4);
  initialize_template_surface_segmentation(line_to_surface_patch_scalar_type, t_line2, t_quad8);
  initialize_template_surface_segmentation(line_to_surface_patch_scalar_type, t_line2, t_quad9);
  initialize_template_surface_segmentation(
      line_to_surface_patch_scalar_type_fixed_size_line2_nurbs9, t_line2, t_nurbs9);
  initialize_template_surface_segmentation(line_to_surface_patch_scalar_type, t_line2, t_tri3);
  initialize_template_surface_segmentation(line_to_surface_patch_scalar_type, t_line2, t_tri6);

  initialize_template_surface_segmentation(double, t_hermite, t_quad4);
  initialize_template_surface_segmentation(double, t_hermite, t_quad8);
  initialize_template_surface_segmentation(double, t_hermite, t_quad9);
  initialize_template_surface_segmentation(double, t_hermite, t_nurbs9);
  initialize_template_surface_segmentation(double, t_hermite, t_tri3);
  initialize_template_surface_segmentation(double, t_hermite, t_tri6);

  initialize_template_surface_segmentation(
      line_to_surface_patch_scalar_type_1st_order, t_hermite, t_quad4);
  initialize_template_surface_segmentation(
      line_to_surface_patch_scalar_type_1st_order, t_hermite, t_quad8);
  initialize_template_surface_segmentation(
      line_to_surface_patch_scalar_type_1st_order, t_hermite, t_quad9);
  initialize_template_surface_segmentation(
      line_to_surface_patch_scalar_type_fixed_size_1st_order_hermite_nurbs9, t_hermite, t_nurbs9);
  initialize_template_surface_segmentation(
      line_to_surface_patch_scalar_type_1st_order, t_hermite, t_tri3);
  initialize_template_surface_segmentation(
      line_to_surface_patch_scalar_type_1st_order, t_hermite, t_tri6);

  initialize_template_surface_segmentation(line_to_surface_patch_scalar_type, t_hermite, t_quad4);
  initialize_template_surface_segmentation(line_to_surface_patch_scalar_type, t_hermite, t_quad8);
  initialize_template_surface_segmentation(line_to_surface_patch_scalar_type, t_hermite, t_quad9);
  initialize_template_surface_segmentation(
      line_to_surface_patch_scalar_type_fixed_size_hermite_nurbs9, t_hermite, t_nurbs9);
  initialize_template_surface_segmentation(line_to_surface_patch_scalar_type, t_hermite, t_tri3);
  initialize_template_surface_segmentation(line_to_surface_patch_scalar_type, t_hermite, t_tri6);
}  // namespace GeometryPair

FOUR_C_NAMESPACE_CLOSE
