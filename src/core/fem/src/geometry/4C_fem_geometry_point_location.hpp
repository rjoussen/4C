// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_FEM_GEOMETRY_POINT_LOCATION_HPP
#define FOUR_C_FEM_GEOMETRY_POINT_LOCATION_HPP

#include "4C_config.hpp"

#include "4C_linalg_fixedsizematrix.hpp"

#include <functional>
#include <optional>
#include <span>
#include <vector>

FOUR_C_NAMESPACE_OPEN

namespace Core::Elements
{
  class Element;
}

namespace Core::Nodes
{
  class Node;
}

namespace Core::Geo
{
  /*!
   * @brief An element containing a point together with the parameter space coordinates of the
   * point within this element.
   */
  struct PointInElement
  {
    const Core::Elements::Element* element;
    Core::LinAlg::Matrix<3, 1> xi;
    //! distance of the point outside of the parameter space of the element, measured in parameter
    //! space coordinates (zero or negative if the point lies inside)
    double violation;
  };

  /*!
   * @brief Position of a node in the configuration in which points are located
   */
  using NodePositionFunction = std::function<Core::LinAlg::Matrix<3, 1>(const Core::Nodes::Node&)>;

  /*!
   * @brief Locate points in a set of 3D elements.
   *
   * The candidate elements of each point are found by a collision search of the bounding boxes of
   * the elements, extended by @p tolerance times their size. Of all candidates, the element
   * containing the point is returned, otherwise the element with the smallest violation of its
   * parameter space bounds if this violation is below @p tolerance.
   *
   * @param elements Elements in which the points are searched (e.g., the column elements of a
   * discretization)
   * @param points Positions of the points in the configuration defined by @p node_position
   * @param node_position Function returning the position of a node in the configuration in which
   * the points are searched (e.g., the current configuration)
   * @param tolerance Tolerance in parameter space for accepting a point outside of all elements
   * @return For each point, the element containing it and the parameter space coordinates of the
   * point, or std::nullopt if no element was found
   */
  std::vector<std::optional<PointInElement>> locate_points_in_elements(
      std::span<const Core::Elements::Element* const> elements,
      std::span<const Core::LinAlg::Matrix<3, 1>> points, const NodePositionFunction& node_position,
      double tolerance = 1.0e-6);
}  // namespace Core::Geo

FOUR_C_NAMESPACE_CLOSE

#endif
