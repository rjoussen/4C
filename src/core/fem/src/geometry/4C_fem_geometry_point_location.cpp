// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_fem_geometry_point_location.hpp"

#include "4C_fem_discretization.hpp"
#include "4C_fem_general_element.hpp"
#include "4C_fem_general_node.hpp"
#include "4C_fem_geometry_element_coordtrafo.hpp"
#include "4C_geometric_search_bvh.hpp"
#include "4C_linalg_serialdensematrix.hpp"

#include <algorithm>
#include <limits>
#include <vector>

FOUR_C_NAMESPACE_OPEN

namespace
{
  /*!
   * @brief Distance of a point outside of the parameter space of an element measured in parameter
   * space coordinates (zero or negative if the point lies inside)
   */
  double parameter_space_violation(
      const Core::FE::CellType celltype, const Core::LinAlg::Matrix<3, 1>& xi)
  {
    switch (celltype)
    {
      case Core::FE::CellType::hex8:
      case Core::FE::CellType::hex20:
      case Core::FE::CellType::hex27:
        return std::max({std::abs(xi(0)), std::abs(xi(1)), std::abs(xi(2))}) - 1.0;
      case Core::FE::CellType::tet4:
      case Core::FE::CellType::tet10:
        return std::max({-xi(0), -xi(1), -xi(2), xi(0) + xi(1) + xi(2) - 1.0});
      case Core::FE::CellType::wedge6:
      case Core::FE::CellType::wedge15:
        return std::max({-xi(0), -xi(1), xi(0) + xi(1) - 1.0, std::abs(xi(2)) - 1.0});
      case Core::FE::CellType::pyramid5:
        return std::max({std::abs(xi(0)) - (1.0 - xi(2)), std::abs(xi(1)) - (1.0 - xi(2)), -xi(2),
            xi(2) - 1.0});
      default:
        FOUR_C_THROW("Point location is not implemented for cell type {}.",
            Core::FE::cell_type_to_string(celltype));
    }
  }

  /*!
   * @brief Compute the parameter space coordinates of a point and the violation of the parameter
   * space bounds for one element
   */
  std::pair<Core::LinAlg::Matrix<3, 1>, double> parameter_space_coordinates(
      const Core::Elements::Element& element, const Core::LinAlg::Matrix<3, 1>& point,
      const Core::Geo::NodePositionFunction& node_position, const double tolerance)
  {
    const int num_nodes = element.num_node();
    Core::LinAlg::SerialDenseMatrix xyze(3, num_nodes);
    const Core::Nodes::Node* const* nodes = element.nodes();
    for (int i = 0; i < num_nodes; ++i)
    {
      const Core::LinAlg::Matrix<3, 1> x = node_position(*nodes[i]);
      for (int d = 0; d < 3; ++d) xyze(d, i) = x(d);
    }

    // quick rejection based on the bounding box of the element nodes, enlarged such that points
    // within the tolerance are not rejected (the parameter space of all cell types spans at least 1
    // per direction)
    const double box_tolerance = 1.0e-8 + tolerance;
    for (int d = 0; d < 3; ++d)
    {
      double x_min = std::numeric_limits<double>::max();
      double x_max = std::numeric_limits<double>::lowest();
      for (int i = 0; i < num_nodes; ++i)
      {
        x_min = std::min(x_min, xyze(d, i));
        x_max = std::max(x_max, xyze(d, i));
      }
      const double margin = box_tolerance * (x_max - x_min) + 1.0e-8;
      if (point(d) < x_min - margin || point(d) > x_max + margin)
        return {Core::LinAlg::Matrix<3, 1>(Core::LinAlg::Initialization::zero),
            std::numeric_limits<double>::max()};
    }

    Core::LinAlg::Matrix<3, 1> xi(Core::LinAlg::Initialization::zero);
    Core::Geo::current_to_volume_element_coordinates(element.shape(), xyze, point, xi);

    if (!std::isfinite(xi(0)) || !std::isfinite(xi(1)) || !std::isfinite(xi(2)))
      return {xi, std::numeric_limits<double>::max()};

    return {xi, parameter_space_violation(element.shape(), xi)};
  }
}  // namespace

std::vector<std::optional<Core::Geo::PointInElement>> Core::Geo::locate_points_in_elements(
    std::span<const Core::Elements::Element* const> elements,
    std::span<const Core::LinAlg::Matrix<3, 1>> points, const NodePositionFunction& node_position,
    const double tolerance)
{
  // the bounding boxes of the elements, extended such that they contain all points within the
  // tolerance in parameter space
  std::vector<std::pair<int, Core::GeometricSearch::BoundingVolume>> element_boxes;
  element_boxes.reserve(elements.size());
  for (std::size_t e = 0; e < elements.size(); ++e)
  {
    const Core::Elements::Element& element = *elements[e];
    Core::GeometricSearch::BoundingVolume box;
    Core::LinAlg::Matrix<3, 1> lower(Core::LinAlg::Initialization::zero);
    Core::LinAlg::Matrix<3, 1> upper(Core::LinAlg::Initialization::zero);
    for (int i = 0; i < element.num_node(); ++i)
    {
      const Core::LinAlg::Matrix<3, 1> x = node_position(*element.nodes()[i]);
      box.add_point(x);
      for (int d = 0; d < 3; ++d)
      {
        lower(d) = i == 0 ? x(d) : std::min(lower(d), x(d));
        upper(d) = i == 0 ? x(d) : std::max(upper(d), x(d));
      }
    }
    double extent = 0.0;
    for (int d = 0; d < 3; ++d) extent = std::max(extent, upper(d) - lower(d));
    box.extend_boundaries((1.0e-8 + tolerance) * extent + 1.0e-8);
    element_boxes.emplace_back(static_cast<int>(e), box);
  }

  std::vector<std::pair<int, Core::GeometricSearch::BoundingVolume>> point_boxes;
  point_boxes.reserve(points.size());
  for (std::size_t p = 0; p < points.size(); ++p)
  {
    Core::GeometricSearch::BoundingVolume box;
    box.add_point(points[p]);
    point_boxes.emplace_back(static_cast<int>(p), box);
  }

  std::vector<std::vector<int>> candidates(points.size());
  for (const auto& [point_index, element_index] :
      Core::GeometricSearch::collision_search(element_boxes, point_boxes))
    candidates[point_index].push_back(element_index);

  // the element containing the point, or the one with the smallest violation of its parameter
  // space bounds (the candidates are sorted for a result independent of the search order)
  std::vector<std::optional<PointInElement>> locations(points.size());
  for (std::size_t p = 0; p < points.size(); ++p)
  {
    std::ranges::sort(candidates[p]);
    double best_violation = std::numeric_limits<double>::max();
    for (int element_index : candidates[p])
    {
      const Core::Elements::Element& element = *elements[element_index];
      const auto [xi, violation] =
          parameter_space_coordinates(element, points[p], node_position, tolerance);
      if (violation < best_violation)
      {
        best_violation = violation;
        locations[p] = PointInElement{&element, xi, violation};
      }
    }
    if (best_violation > tolerance) locations[p].reset();
  }
  return locations;
}

FOUR_C_NAMESPACE_CLOSE
