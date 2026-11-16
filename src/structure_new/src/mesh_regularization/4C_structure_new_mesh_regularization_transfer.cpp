// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_structure_new_mesh_regularization_transfer.hpp"

#include "4C_comm_exporter.hpp"
#include "4C_comm_mpi_utils.hpp"
#include "4C_fem_discretization.hpp"
#include "4C_fem_general_node.hpp"
#include "4C_fem_general_utils_fem_shapefunctions.hpp"
#include "4C_fem_geometry_point_location.hpp"
#include "4C_geometric_search_bounding_volume.hpp"
#include "4C_geometric_search_distributed_tree.hpp"
#include "4C_linalg_serialdensevector.hpp"
#include "4C_linalg_utils_sparse_algebra_manipulation.hpp"
#include "4C_mat_history_transfer.hpp"
#include "4C_mat_so3_material.hpp"
#include "4C_material_base.hpp"
#include "4C_material_parameter_base.hpp"
#include "4C_solid_ele.hpp"
#include "4C_structure_new_mesh_regularization_remote_evaluation.hpp"
#include "4C_utils_exceptions.hpp"

#include <algorithm>
#include <cmath>
#include <functional>
#include <map>
#include <memory>
#include <optional>
#include <set>
#include <span>
#include <tuple>

FOUR_C_NAMESPACE_OPEN

namespace
{
  //! tolerance in parameter space for new node positions: the weighted gaps of the boundary
  //! sliding keep the boundary only in an averaged sense, at bends the new faces cut the old ones
  //! by a fraction of the element size, so that boundary nodes may lie slightly outside the old
  //! mesh
  constexpr double node_location_tolerance = 1.0;

  //! tolerance in parameter space for new Gauss points: where the boundary is concave (e.g.,
  //! wrapped around a tool edge), the straight edges between sliding boundary nodes cut across it,
  //! and Gauss points near them lie outside of the old mesh; their state is extrapolated from the
  //! patch of the nearest old element
  constexpr double gauss_point_location_tolerance = 1.0;

  const Discret::Elements::Solid<3>& cast_to_solid(const Core::Elements::Element& element)
  {
    const auto* solid = dynamic_cast<const Discret::Elements::Solid<3>*>(&element);
    FOUR_C_ASSERT_ALWAYS(solid != nullptr,
        "Mesh regularization only supports 3D solid elements, but element {} is not one.",
        element.id());
    return *solid;
  }

  //! the transferred quantities of an element at a Gauss point: the deformation gradient of the
  //! material and its history
  std::optional<std::vector<Mat::HistoryTransferQuantity>> history_quantities(
      const Discret::Elements::Solid<3>& element)
  {
    auto material_quantities = element.solid_material()->history_transfer_quantities();
    if (!material_quantities) return std::nullopt;

    std::vector<Mat::HistoryTransferQuantity> quantities{
        {"deformation_gradient", Mat::HistoryTransferKind::invertible_tensor}};
    quantities.insert(quantities.end(), material_quantities->begin(), material_quantities->end());
    return quantities;
  }

  //! number of values of the deformation gradient in the transferred data
  constexpr int deformation_gradient_size = 9;

  //! the state of a Gauss point from the transferred values (deformation gradient and history)
  Discret::Elements::TransferredGaussPointState gauss_point_state(std::span<const double> values)
  {
    Discret::Elements::TransferredGaussPointState state;
    for (int r = 0; r < 3; ++r)
      for (int c = 0; c < 3; ++c) state.deformation_gradient(r, c) = values[3 * r + c];
    state.material_history.assign(values.begin() + deformation_gradient_size, values.end());
    return state;
  }

  //! the state of the old mesh at the column nodes
  struct OldNodalState
  {
    std::vector<Core::LinAlg::Matrix<3, 1>> reference_position;
    std::vector<Core::LinAlg::Matrix<3, 1>> current_position;
    //! displacement followed by all nodal fields
    std::vector<std::vector<Core::LinAlg::Matrix<3, 1>>> fields;
  };

  Core::LinAlg::Matrix<3, 1> nodal_values(const Core::FE::Discretization& discretization,
      const Core::Nodes::Node& node, const Core::LinAlg::Vector<double>& col_vector)
  {
    const std::vector<int> dofs = discretization.dof(0, &node);
    Core::LinAlg::Matrix<3, 1> values;
    for (int d = 0; d < 3; ++d)
      values(d) = col_vector.local_values_as_span()[col_vector.get_map().lid(dofs[d])];
    return values;
  }

  OldNodalState gather_old_nodal_state(const Core::FE::Discretization& discretization,
      const std::vector<const Core::LinAlg::Vector<double>*>& fields)
  {
    std::vector<Core::LinAlg::Vector<double>> col_fields;
    for (const auto* field : fields)
    {
      col_fields.emplace_back(*discretization.dof_col_map());
      Core::LinAlg::export_to(*field, col_fields.back());
    }

    const int num_col_nodes = discretization.num_my_col_nodes();
    OldNodalState state;
    state.reference_position.resize(num_col_nodes);
    state.current_position.resize(num_col_nodes);
    state.fields.assign(fields.size(), std::vector<Core::LinAlg::Matrix<3, 1>>(num_col_nodes));

    for (int lid = 0; lid < num_col_nodes; ++lid)
    {
      const Core::Nodes::Node& node = *discretization.l_col_node(lid);
      for (int d = 0; d < 3; ++d) state.reference_position[lid](d) = node.x()[d];
      for (std::size_t f = 0; f < fields.size(); ++f)
        state.fields[f][lid] = nodal_values(discretization, node, col_fields[f]);

      // the first field is the displacement
      state.current_position[lid] = state.reference_position[lid];
      state.current_position[lid].update(1.0, state.fields[0][lid], 1.0);
    }
    return state;
  }

  /*!
   * @brief The converged state at all Gauss points of the column elements whose material supports
   * the transfer of its history
   *
   * The deformation is evaluated with the state "displacement" of the discretization.
   */
  std::map<int, std::vector<double>> gather_old_history(
      const Core::FE::Discretization& discretization)
  {
    std::map<int, std::vector<double>> history;
    for (int lid = 0; lid < discretization.num_my_row_elements(); ++lid)
    {
      const Core::Elements::Element* element = discretization.l_row_element(lid);
      auto& solid = const_cast<Discret::Elements::Solid<3>&>(cast_to_solid(*element));
      const auto quantities = history_quantities(solid);
      if (!quantities) continue;
      const int size = Mat::history_transfer_size(*quantities);
      const int num_gp = solid.get_gauss_rule().num_points();

      const std::vector<Core::LinAlg::Tensor<double, 3, 3>> deformation_gradients =
          solid.evaluate_material_deformation_gradients(discretization);

      std::vector<double>& data = history[element->id()];
      data.reserve(size * num_gp);
      for (int gp = 0; gp < num_gp; ++gp)
      {
        const std::size_t old_size = data.size();
        for (int r = 0; r < 3; ++r)
          for (int c = 0; c < 3; ++c) data.push_back(deformation_gradients[gp](r, c));
        solid.solid_material()->gather_history_for_transfer(gp, data);
        FOUR_C_ASSERT_ALWAYS(static_cast<int>(data.size() - old_size) == size,
            "The material of element {} gathered {} history values, but declared {}.",
            element->id(), data.size() - old_size, size);
      }
    }

    Core::Communication::Exporter exporter(*discretization.element_row_map(),
        *discretization.element_col_map(), discretization.get_comm());
    exporter.do_export(history);
    return history;
  }

  bool same_quantities(const std::optional<std::vector<Mat::HistoryTransferQuantity>>& a,
      const std::vector<Mat::HistoryTransferQuantity>& b)
  {
    return a.has_value() && a->size() == b.size() &&
           std::equal(a->begin(), a->end(), b.begin(),
               [](const auto& x, const auto& y) { return x.name == y.name && x.kind == y.kind; });
  }

  //! position of a Gauss point of an element in the configuration given by the node positions
  Core::LinAlg::Matrix<3, 1> gauss_point_position(const Discret::Elements::Solid<3>& element,
      const int gp, const Core::Geo::NodePositionFunction& node_position)
  {
    const double* xi = element.get_gauss_rule().point(gp);
    Core::LinAlg::SerialDenseVector N(element.num_node());
    Core::FE::shape_function_3d(N, xi[0], xi[1], xi[2], element.shape());
    Core::LinAlg::Matrix<3, 1> x(Core::LinAlg::Initialization::zero);
    for (int i = 0; i < element.num_node(); ++i)
      x.update(N(i), node_position(*element.nodes()[i]), 1.0);
    return x;
  }

  //! size of an element in the configuration given by the node positions (bounding box diagonal)
  double element_size(
      const Core::Elements::Element& element, const Core::Geo::NodePositionFunction& node_position)
  {
    Core::LinAlg::Matrix<3, 1> lower = node_position(*element.nodes()[0]);
    Core::LinAlg::Matrix<3, 1> upper = lower;
    for (int i = 1; i < element.num_node(); ++i)
    {
      const Core::LinAlg::Matrix<3, 1> x = node_position(*element.nodes()[i]);
      for (int d = 0; d < 3; ++d)
      {
        lower(d) = std::min(lower(d), x(d));
        upper(d) = std::max(upper(d), x(d));
      }
    }
    upper.update(-1.0, lower, 1.0);
    return upper.norm2() / std::sqrt(3.0);
  }

  int material_id(const Core::Elements::Element& element)
  {
    return element.material()->parameter()->id();
  }

  std::vector<double> shape_functions(
      const Core::Elements::Element& element, const Core::LinAlg::Matrix<3, 1>& xi)
  {
    Core::LinAlg::SerialDenseVector values(element.num_node());
    Core::FE::shape_function_3d(values, xi(0), xi(1), xi(2), element.shape());
    return {values.values(), values.values() + element.num_node()};
  }

  //! kind of a point at which the old state is evaluated
  enum class PointKind
  {
    node,
    gauss_point
  };

  //! a point at which the old state is evaluated, located among the old elements with its material
  struct OldStatePoint
  {
    PointKind kind;
    int material;
    Core::LinAlg::Matrix<3, 1> position;
  };

  //! the old state at a point, evaluated in the old element containing it
  struct OldStateAtPoint
  {
    //! violation of the parameter space bounds of the old element by the point
    double violation;
    int element;
    std::vector<double> values;
  };

  //! evaluation of the old state in the old element containing a point
  using OldStateEvaluator =
      std::function<std::vector<double>(const Core::Elements::Element& old_element,
          const Core::LinAlg::Matrix<3, 1>& xi, const OldStatePoint& point)>;

  /*!
   * @brief Evaluate the old state at points, which are located among the old elements of all
   * processors
   *
   * The candidate elements of each point are found by a global search over the bounding boxes of
   * the old row elements. The owners of the candidates locate the point among those with its
   * material and evaluate the old state in the element containing it, of which the one with the
   * smallest violation of its parameter space bounds is taken. Since the owner of an old element
   * evaluates its state, the movement of the mesh is not limited by the ghosting. Collective.
   *
   * @return For each point, the old state, or std::nullopt if the point could not be located within
   * the tolerance of its kind
   */
  std::vector<std::optional<OldStateAtPoint>> evaluate_old_state_at_points(
      const Core::FE::Discretization& discretization, const std::vector<OldStatePoint>& points,
      const Core::Geo::NodePositionFunction& old_position, const OldStateEvaluator& evaluate)
  {
    // the bounding boxes of the old row elements, extended such that they contain all points within
    // the larger tolerance in parameter space
    std::vector<std::pair<int, Core::GeometricSearch::BoundingVolume>> element_boxes;
    for (int lid = 0; lid < discretization.num_my_row_elements(); ++lid)
    {
      const Core::Elements::Element& element = *discretization.l_row_element(lid);
      Core::GeometricSearch::BoundingVolume box;
      for (int i = 0; i < element.num_node(); ++i) box.add_point(old_position(*element.nodes()[i]));
      const double extent = std::sqrt(3.0) * element_size(element, old_position);
      box.extend_boundaries((1.0e-8 + gauss_point_location_tolerance) * extent + 1.0e-8);
      element_boxes.emplace_back(element.id(), box);
    }

    std::vector<std::pair<int, Core::GeometricSearch::BoundingVolume>> point_boxes;
    for (std::size_t p = 0; p < points.size(); ++p)
    {
      Core::GeometricSearch::BoundingVolume box;
      box.add_point(points[p].position);
      point_boxes.emplace_back(static_cast<int>(p), box);
    }

    // the candidate elements of each point, by their processor
    std::vector<std::map<int, std::vector<int>>> candidates(points.size());
    for (const auto& collision : Core::GeometricSearch::global_collision_search(
             element_boxes, point_boxes, discretization.get_comm()))
    {
      candidates[collision.lid_predicate][collision.pid_primitive].push_back(
          collision.gid_primitive);
    }

    // a request: the kind and material of the point, its position and the candidate elements
    std::vector<Solid::MeshRegularization::RemoteRequest> requests;
    std::vector<std::size_t> request_points;
    for (std::size_t p = 0; p < points.size(); ++p)
    {
      for (auto& [processor, ids] : candidates[p])
      {
        std::ranges::sort(ids);
        std::vector<double> data{static_cast<double>(points[p].kind),
            static_cast<double>(points[p].material), points[p].position(0), points[p].position(1),
            points[p].position(2)};
        data.insert(data.end(), ids.begin(), ids.end());
        requests.push_back({processor, std::move(data)});
        request_points.push_back(p);
      }
    }

    // an answer: the violation, the id of the old element and the old state there
    const auto answer = [&](std::span<const double> request) -> std::vector<double>
    {
      OldStatePoint point{static_cast<PointKind>(std::lround(request[0])),
          static_cast<int>(std::lround(request[1])),
          Core::LinAlg::Matrix<3, 1>(Core::LinAlg::Initialization::zero)};
      for (int d = 0; d < 3; ++d) point.position(d) = request[2 + d];

      std::vector<const Core::Elements::Element*> elements;
      for (const double id : request.subspan(5))
      {
        const Core::Elements::Element* element =
            discretization.g_element(static_cast<int>(std::lround(id)));
        if (material_id(*element) == point.material) elements.push_back(element);
      }

      const double tolerance =
          point.kind == PointKind::node ? node_location_tolerance : gauss_point_location_tolerance;
      const auto location = Core::Geo::locate_points_in_elements(elements,
          std::span<const Core::LinAlg::Matrix<3, 1>>(&point.position, 1), old_position, tolerance)
                                .front();
      if (!location) return {};

      std::vector<double> result{location->violation, static_cast<double>(location->element->id())};
      const std::vector<double> values = evaluate(*location->element, location->xi, point);
      result.insert(result.end(), values.begin(), values.end());
      return result;
    };

    const std::vector<std::vector<double>> answers =
        Solid::MeshRegularization::answer_on_processors(
            discretization.get_comm(), requests, answer);

    // the element with the smallest violation (the smallest id among equal ones) contains the point
    std::vector<std::optional<OldStateAtPoint>> old_state(points.size());
    for (std::size_t r = 0; r < answers.size(); ++r)
    {
      if (answers[r].empty()) continue;
      OldStateAtPoint candidate{answers[r][0], static_cast<int>(std::lround(answers[r][1])),
          std::vector<double>(answers[r].begin() + 2, answers[r].end())};
      std::optional<OldStateAtPoint>& best = old_state[request_points[r]];
      if (!best || std::tie(candidate.violation, candidate.element) <
                       std::tie(best->violation, best->element))
        best = std::move(candidate);
    }
    return old_state;
  }
}  // namespace

void Solid::MeshRegularization::move_nodes_and_transfer_state(
    Core::FE::Discretization& discretization, const Core::LinAlg::Vector<double>& node_displacement,
    Core::LinAlg::Vector<double>& displacement,
    const std::vector<Core::LinAlg::Vector<double>*>& nodal_fields,
    const TransferParameters& parameters)
{
  // --- state of the old mesh -------------------------------------------------------------------
  std::vector<const Core::LinAlg::Vector<double>*> old_fields{&displacement};
  for (const auto* field : nodal_fields) old_fields.push_back(field);
  const OldNodalState old_state = gather_old_nodal_state(discretization, old_fields);
  discretization.clear_state();
  discretization.set_state(0, "displacement", displacement);
  const std::map<int, std::vector<double>> old_history = gather_old_history(discretization);

  const Core::LinAlg::Map& node_col_map = *discretization.node_col_map();
  const Core::Geo::NodePositionFunction old_position = [&](const Core::Nodes::Node& node)
  { return old_state.current_position[node_col_map.lid(node.id())]; };

  Core::LinAlg::Vector<double> node_displacement_col(*discretization.dof_col_map());
  Core::LinAlg::export_to(node_displacement, node_displacement_col);

  const auto moved = [&](const Core::Nodes::Node& node)
  { return nodal_values(discretization, node, node_displacement_col).norm_inf() != 0.0; };

  const auto element_moved = [&](const Core::Elements::Element& element)
  {
    return std::ranges::any_of(std::span(element.nodes(), element.num_node()),
        [&](const auto* node) { return moved(*node); });
  };

  const Core::Geo::NodePositionFunction new_node_position = [&](const Core::Nodes::Node& node)
  {
    Core::LinAlg::Matrix<3, 1> x = old_position(node);
    x.update(1.0, nodal_values(discretization, node, node_displacement_col), 1.0);
    return x;
  };

  // --- the points of the new mesh: the moved nodes and the Gauss points of the moved elements --
  std::vector<OldStatePoint> points;
  std::vector<const Core::Nodes::Node*> moved_nodes;
  for (int row_lid = 0; row_lid < discretization.num_my_row_nodes(); ++row_lid)
  {
    const Core::Nodes::Node& node = *discretization.l_row_node(row_lid);
    if (!moved(node)) continue;
    moved_nodes.push_back(&node);
    points.push_back({PointKind::node, material_id(*node.adjacent_elements()[0].user_element()),
        new_node_position(node)});
  }

  std::map<int, std::size_t> first_gauss_point;
  for (int lid = 0; lid < discretization.num_my_row_elements(); ++lid)
  {
    const Core::Elements::Element* element = discretization.l_row_element(lid);
    if (!element_moved(*element)) continue;
    const auto& solid = cast_to_solid(*element);
    FOUR_C_ASSERT_ALWAYS(history_quantities(solid).has_value(),
        "The material of element {} does not support the transfer of its history, which is "
        "required since the mesh regularization moves its nodes. Exclude the element from the "
        "mesh movement with ALE Dirichlet conditions.",
        solid.id());
    first_gauss_point[solid.id()] = points.size();
    for (int gp = 0; gp < solid.get_gauss_rule().num_points(); ++gp)
    {
      points.push_back({PointKind::gauss_point, material_id(solid),
          gauss_point_position(solid, gp, new_node_position)});
    }
  }

  // --- the old state at the points, evaluated by the owners of the old elements ----------------
  // at a node, the further nodal fields interpolated with the shape functions of the old element;
  // at a Gauss point, the deformation gradient and the history approximated from the Gauss points
  // of the patch of old elements with the same material sharing a node with the old element, in
  // the old current configuration relative to the point and scaled by the element size
  const OldStateEvaluator evaluate = [&](const Core::Elements::Element& old_element,
                                         const Core::LinAlg::Matrix<3, 1>& xi,
                                         const OldStatePoint& point) -> std::vector<double>
  {
    if (point.kind == PointKind::node)
    {
      const std::vector<double> N = shape_functions(old_element, xi);
      std::vector<double> values;
      for (std::size_t f = 1; f < old_state.fields.size(); ++f)
      {
        Core::LinAlg::Matrix<3, 1> value(Core::LinAlg::Initialization::zero);
        for (std::size_t i = 0; i < N.size(); ++i)
          value.update(
              N[i], old_state.fields[f][node_col_map.lid(old_element.nodes()[i]->id())], 1.0);
        for (int d = 0; d < 3; ++d) values.push_back(value(d));
      }
      return values;
    }

    const auto& old_solid = cast_to_solid(old_element);
    const auto quantities = history_quantities(old_solid);
    FOUR_C_ASSERT_ALWAYS(quantities.has_value(),
        "The material of element {} does not support the transfer of its history.", old_solid.id());
    const int size = Mat::history_transfer_size(*quantities);
    const double h = element_size(old_solid, old_position);

    std::set<int> patch;
    for (int i = 0; i < old_solid.num_node(); ++i)
      for (auto adjacent : old_solid.nodes()[i]->adjacent_elements())
        patch.insert(adjacent.user_element()->id());

    std::vector<std::vector<double>> gauss_point_data;
    std::vector<Core::LinAlg::Matrix<3, 1>> gauss_point_locations;
    for (int id : patch)
    {
      const auto& element = cast_to_solid(*discretization.g_element(id));
      if (!same_quantities(history_quantities(element), *quantities)) continue;

      const std::vector<double>& data = old_history.at(element.id());
      for (int gp = 0; gp < element.get_gauss_rule().num_points(); ++gp)
      {
        gauss_point_data.emplace_back(data.begin() + gp * size, data.begin() + (gp + 1) * size);
        Core::LinAlg::Matrix<3, 1> x = gauss_point_position(element, gp, old_position);
        x.update(-1.0, point.position, 1.0);
        x.scale(1.0 / h);
        gauss_point_locations.push_back(x);
      }
    }

    return Mat::interpolate_history(*quantities, gauss_point_data, gauss_point_locations,
        Core::LinAlg::Matrix<3, 1>(Core::LinAlg::Initialization::zero),
        parameters.history_weight_decay);
  };

  const std::vector<std::optional<OldStateAtPoint>> old_state_at_points =
      evaluate_old_state_at_points(discretization, points, old_position, evaluate);

  // --- nodes: displacement relative to the unchanged reference position and nodal fields -------
  std::vector<Core::LinAlg::Vector<double>*> new_fields{&displacement};
  for (auto* field : nodal_fields) new_fields.push_back(field);

  for (std::size_t n = 0; n < moved_nodes.size(); ++n)
  {
    const Core::Nodes::Node& node = *moved_nodes[n];
    const Core::LinAlg::Matrix<3, 1>& new_position = points[n].position;
    const auto& old_values = old_state_at_points[n];
    if (!old_values.has_value())
    {
      const Core::LinAlg::Matrix<3, 1> movement =
          nodal_values(discretization, node, node_displacement_col);
      FOUR_C_THROW(
          "The new position ({}, {}, {}) of node {} could not be located in the old mesh. The node "
          "moved by ({}, {}, {}) and left the old mesh.",
          new_position(0), new_position(1), new_position(2), node.id(), movement(0), movement(1),
          movement(2));
    }

    std::vector<Core::LinAlg::Matrix<3, 1>> new_values(new_fields.size());
    new_values[0].update(
        1.0, new_position, -1.0, old_state.reference_position[node_col_map.lid(node.id())], 0.0);
    for (std::size_t f = 1; f < new_fields.size(); ++f)
      for (int d = 0; d < 3; ++d) new_values[f](d) = old_values->values[3 * (f - 1) + d];

    const std::vector<int> dofs = discretization.dof(0, &node);
    for (int d = 0; d < 3; ++d)
    {
      const int dof_lid = discretization.dof_row_map()->lid(dofs[d]);
      for (std::size_t f = 0; f < new_fields.size(); ++f)
        new_fields[f]->get_values()[dof_lid] = new_values[f](d);
    }
  }

  // --- Gauss points: deformation gradient of the material and its history ---------------------
  discretization.clear_state();
  discretization.set_state(0, "displacement", displacement);

  // the owned elements receive their state, the values are recorded for the ghosted copies
  std::map<int, std::vector<double>> transferred_values;
  for (int lid = 0; lid < discretization.num_my_row_elements(); ++lid)
  {
    Core::Elements::Element* element = discretization.l_row_element(lid);
    if (!element_moved(*element)) continue;

    auto& solid = const_cast<Discret::Elements::Solid<3>&>(cast_to_solid(*element));
    std::vector<double>& element_values = transferred_values[solid.id()];

    solid.receive_transferred_state(discretization,
        [&](int gp, const Core::LinAlg::Tensor<double, 3>& new_position)
        {
          const auto& old_values = old_state_at_points[first_gauss_point.at(solid.id()) + gp];
          FOUR_C_ASSERT_ALWAYS(old_values.has_value(),
              "Gauss point {} of element {} at ({}, {}, {}) could not be located in the old mesh.",
              gp, solid.id(), new_position(0), new_position(1), new_position(2));

          element_values.insert(
              element_values.end(), old_values->values.begin(), old_values->values.end());
          return gauss_point_state(old_values->values);
        });
  }

  // the discretization also evaluates the ghosted elements, which receive the state of the owner
  Core::Communication::Exporter exporter(*discretization.element_row_map(),
      *discretization.element_col_map(), discretization.get_comm());
  exporter.do_export(transferred_values);
  const int my_rank = Core::Communication::my_mpi_rank(discretization.get_comm());
  for (int lid = 0; lid < discretization.num_my_col_elements(); ++lid)
  {
    Core::Elements::Element* element = discretization.l_col_element(lid);
    if (element->owner() == my_rank) continue;
    const auto it = transferred_values.find(element->id());
    if (it == transferred_values.end()) continue;

    auto& solid = const_cast<Discret::Elements::Solid<3>&>(cast_to_solid(*element));
    const std::span<const double> values = it->second;
    const std::size_t size = values.size() / solid.get_gauss_rule().num_points();
    solid.receive_transferred_state(discretization,
        [&](int gp, const Core::LinAlg::Tensor<double, 3>&)
        { return gauss_point_state(values.subspan(gp * size, size)); });
  }

  discretization.clear_state();
}

FOUR_C_NAMESPACE_CLOSE
