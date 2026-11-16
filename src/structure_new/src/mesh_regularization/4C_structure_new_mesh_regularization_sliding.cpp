// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_structure_new_mesh_regularization_sliding.hpp"

#include "4C_fem_discretization.hpp"
#include "4C_fem_general_element.hpp"
#include "4C_fem_general_node.hpp"
#include "4C_fem_general_utils_local_connectivity_matrices.hpp"
#include "4C_geometric_search_distributed_tree.hpp"
#include "4C_linalg_sparsematrix.hpp"
#include "4C_linalg_utils_sparse_algebra_math.hpp"
#include "4C_structure_new_mesh_regularization_remote_evaluation.hpp"
#include "4C_utils_exceptions.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <deque>
#include <functional>
#include <limits>
#include <numbers>
#include <numeric>
#include <set>

FOUR_C_NAMESPACE_OPEN

namespace
{
  Core::LinAlg::Matrix<3, 1> cross(
      const Core::LinAlg::Matrix<3, 1>& a, const Core::LinAlg::Matrix<3, 1>& b)
  {
    Core::LinAlg::Matrix<3, 1> c;
    c(0) = a(1) * b(2) - a(2) * b(1);
    c(1) = a(2) * b(0) - a(0) * b(2);
    c(2) = a(0) * b(1) - a(1) * b(0);
    return c;
  }

  //! number of corner nodes of a face with the given number of nodes
  int num_corners(const std::size_t num_face_nodes)
  {
    if (num_face_nodes == 3 || num_face_nodes == 6) return 3;
    if (num_face_nodes == 4 || num_face_nodes == 8 || num_face_nodes == 9) return 4;
    FOUR_C_THROW("Unsupported face with {} nodes.", num_face_nodes);
  }

  //! enlargement of the search boxes around the moved nodes, which covers the single precision of
  //! the bounding boxes
  constexpr double bounding_box_tolerance = 1.0e-6;

  //! number of iterations of the alternating projection onto the two patches of a feature edge
  constexpr int num_edge_projections = 5;

  //! largest number of faces of a 3D element, which numbers the faces globally
  constexpr int max_faces_per_element = 6;

  //! unit normal of a face with the given corner positions
  Core::LinAlg::Matrix<3, 1> face_normal(const std::vector<Core::LinAlg::Matrix<3, 1>>& x)
  {
    Core::LinAlg::Matrix<3, 1> a;
    Core::LinAlg::Matrix<3, 1> b;
    if (x.size() == 4)
    {
      a.update(1.0, x[2], -1.0, x[0], 0.0);
      b.update(1.0, x[3], -1.0, x[1], 0.0);
    }
    else
    {
      a.update(1.0, x[1], -1.0, x[0], 0.0);
      b.update(1.0, x[2], -1.0, x[0], 0.0);
    }
    Core::LinAlg::Matrix<3, 1> normal = cross(a, b);
    normal.scale(1.0 / normal.norm2());
    return normal;
  }

  //! the sorted ids of the two nodes of an edge
  std::pair<int, int> edge_key(const int a, const int b) { return {std::min(a, b), std::max(a, b)}; }
}  // namespace

Solid::MeshRegularization::BoundarySliding::BoundarySliding(
    const Core::FE::Discretization& discretization, const double feature_angle)
    : discretization_(discretization),
      feature_angle_cosine_(std::cos(feature_angle * std::numbers::pi / 180.0))
{
  // count the occurrences of the faces of the column elements: all elements adjacent to a row node
  // are column elements, so the faces with a row node are classified correctly, and the boundary
  // faces among them belong to one element only, which numbers them globally
  struct FaceCount
  {
    int count = 0;
    int id = -1;
    std::vector<int> ordered_ids;
  };
  std::map<std::vector<int>, FaceCount> face_count;
  for (int lid = 0; lid < discretization_.num_my_col_elements(); ++lid)
  {
    const Core::Elements::Element* element = discretization_.l_col_element(lid);
    const std::vector<std::vector<int>> surfaces =
        Core::FE::get_ele_node_numbering_surfaces(element->shape());
    FOUR_C_ASSERT_ALWAYS(static_cast<int>(surfaces.size()) <= max_faces_per_element,
        "Element {} has more than {} faces.", element->id(), max_faces_per_element);
    for (std::size_t s = 0; s < surfaces.size(); ++s)
    {
      std::vector<int> ids;
      for (int i = 0; i < num_corners(surfaces[s].size()); ++i)
        ids.push_back(element->node_ids()[surfaces[s][i]]);
      std::vector<int> key = ids;
      std::ranges::sort(key);
      FaceCount& entry = face_count[key];
      ++entry.count;
      entry.ordered_ids = ids;
      entry.id = element->id() * max_faces_per_element + static_cast<int>(s);
    }
  }

  // the boundary faces with a row node, with their normals in the initial configuration; the
  // processor owning the smallest node of a face owns the face
  const Core::LinAlg::Map& node_row_map = *discretization_.node_row_map();
  for (const auto& [key, entry] : face_count)
  {
    if (entry.count != 1) continue;
    if (std::ranges::none_of(key, [&](int id) { return node_row_map.lid(id) >= 0; })) continue;

    Face face;
    face.node_ids = entry.ordered_ids;
    face.corners = corner_positions(face.node_ids);
    face.normal = face_normal(face.corners);
    face.initial_normal = face.normal;
    row_node_faces_[entry.id] = face;
    if (node_row_map.lid(key.front()) >= 0) owned_faces_.push_back(entry.id);

    for (int id : entry.ordered_ids)
      if (node_row_map.lid(id) >= 0) node_faces_[id].push_back(entry.id);
  }
}

Core::LinAlg::Matrix<3, 1> Solid::MeshRegularization::BoundarySliding::position(
    const int node_id) const
{
  const Core::Nodes::Node* node = discretization_.g_node(node_id);
  Core::LinAlg::Matrix<3, 1> x;
  for (int d = 0; d < 3; ++d) x(d) = node->x()[d];
  return x;
}

std::vector<Core::LinAlg::Matrix<3, 1>>
Solid::MeshRegularization::BoundarySliding::corner_positions(
    const std::vector<int>& node_ids) const
{
  std::vector<Core::LinAlg::Matrix<3, 1>> x;
  for (int id : node_ids) x.push_back(position(id));
  return x;
}

std::vector<int> Solid::MeshRegularization::BoundarySliding::dof_lids(const int node_id) const
{
  const std::vector<int> dofs = discretization_.dof(0, discretization_.g_node(node_id));
  std::vector<int> lids;
  for (int d = 0; d < 3; ++d) lids.push_back(discretization_.dof_row_map()->lid(dofs[d]));
  return lids;
}

bool Solid::MeshRegularization::BoundarySliding::same_patch(const int a, const int b) const
{
  return std::abs(faces_.at(a).initial_normal.dot(faces_.at(b).initial_normal)) >=
         feature_angle_cosine_;
}

void Solid::MeshRegularization::BoundarySliding::add_face(const int id, Face face)
{
  const std::size_t n = face.node_ids.size();
  for (std::size_t i = 0; i < n; ++i)
    edge_faces_[edge_key(face.node_ids[i], face.node_ids[(i + 1) % n])].push_back(id);
  faces_[id] = std::move(face);
}

void Solid::MeshRegularization::BoundarySliding::setup()
{
  // the known faces are the boundary faces with a row node, in the current reference configuration
  faces_.clear();
  edge_faces_.clear();
  for (const auto& [id, row_node_face] : row_node_faces_)
  {
    Face face = row_node_face;
    face.corners = corner_positions(face.node_ids);
    face.normal = face_normal(face.corners);
    add_face(id, std::move(face));
  }

  // bounding boxes of the owned faces for the global search of the closest points
  face_boxes_.clear();
  for (int id : owned_faces_)
  {
    Core::GeometricSearch::BoundingVolume box;
    for (const auto& x : faces_.at(id).corners) box.add_point(x);
    face_boxes_.emplace_back(id, box);
  }

  // the faces at each row node grouped into patches: faces sharing an edge belong to the same patch
  // if their normals in the initial configuration enclose an angle below the feature angle
  constraints_.clear();
  for (const auto& [node_id, face_ids] : node_faces_)
  {
    std::vector<int> patch_of(face_ids.size());
    std::iota(patch_of.begin(), patch_of.end(), 0);
    const auto root = [&](int i)
    {
      while (patch_of[i] != i) i = patch_of[i];
      return i;
    };
    for (std::size_t i = 0; i < face_ids.size(); ++i)
    {
      for (std::size_t j = i + 1; j < face_ids.size(); ++j)
      {
        const auto& a = faces_.at(face_ids[i]).node_ids;
        const auto& b = faces_.at(face_ids[j]).node_ids;
        const auto shared = std::ranges::count_if(
            a, [&](int id) { return std::ranges::find(b, id) != b.end(); });
        if (shared >= 2 && same_patch(face_ids[i], face_ids[j]))
          patch_of[root(static_cast<int>(i))] = root(static_cast<int>(j));
      }
    }

    NodeConstraint constraint;
    std::vector<int> roots;
    for (std::size_t i = 0; i < face_ids.size(); ++i)
    {
      const auto it = std::ranges::find(roots, root(static_cast<int>(i)));
      if (it == roots.end())
      {
        roots.push_back(root(static_cast<int>(i)));
        constraint.patch_faces.push_back({face_ids[i]});
      }
      else
        constraint.patch_faces[it - roots.begin()].push_back(face_ids[i]);
    }

    // the normal of each patch at the node is the average of the normals of its faces
    for (const auto& faces : constraint.patch_faces)
    {
      Core::LinAlg::Matrix<3, 1> normal(Core::LinAlg::Initialization::zero);
      for (int face_id : faces)
        normal.update(faces_.at(face_id).normal.dot(faces_.at(faces.front()).normal) >= 0.0 ? 1.0
                                                                                            : -1.0,
            faces_.at(face_id).normal, 1.0);
      normal.scale(1.0 / normal.norm2());
      constraint.normals.push_back(normal);
    }
    constraint.active_normals = constraint.normals;

    constraints_[node_id] = constraint;
  }
}

void Solid::MeshRegularization::BoundarySliding::fetch_faces(
    const std::map<int, int>& face_processors)
{
  // a face is sent as its number of corners, the ids and positions of its corners and its normal
  // in the initial configuration
  std::vector<RemoteRequest> requests;
  std::vector<int> ids;
  for (const auto& [id, processor] : face_processors)
  {
    requests.push_back({processor, {static_cast<double>(id)}});
    ids.push_back(id);
  }

  const auto answer = [&](std::span<const double> request) -> std::vector<double>
  {
    const Face& face = faces_.at(static_cast<int>(std::lround(request[0])));
    std::vector<double> data{static_cast<double>(face.node_ids.size())};
    for (int id : face.node_ids) data.push_back(id);
    for (const auto& x : face.corners) data.insert(data.end(), {x(0), x(1), x(2)});
    for (int d = 0; d < 3; ++d) data.push_back(face.initial_normal(d));
    return data;
  };

  const std::vector<std::vector<double>> answers =
      answer_on_processors(discretization_.get_comm(), requests, answer);

  for (std::size_t r = 0; r < answers.size(); ++r)
  {
    const std::vector<double>& data = answers[r];
    const auto n = static_cast<std::size_t>(std::lround(data[0]));
    Face face;
    for (std::size_t i = 0; i < n; ++i)
      face.node_ids.push_back(static_cast<int>(std::lround(data[1 + i])));
    for (std::size_t i = 0; i < n; ++i)
    {
      Core::LinAlg::Matrix<3, 1> x;
      for (int d = 0; d < 3; ++d) x(d) = data[1 + n + 3 * i + d];
      face.corners.push_back(x);
    }
    for (int d = 0; d < 3; ++d) face.initial_normal(d) = data[1 + 4 * n + d];
    face.normal = face_normal(face.corners);
    add_face(ids[r], std::move(face));
  }
}

std::vector<int> Solid::MeshRegularization::BoundarySliding::patch_candidates(
    const std::vector<int>& patch_faces, const std::vector<int>& candidates) const
{
  // the candidates reachable from the faces of the patch at the node without crossing a feature
  // edge
  std::set<int> reached(patch_faces.begin(), patch_faces.end());
  std::deque<int> queue(patch_faces.begin(), patch_faces.end());
  while (!queue.empty())
  {
    const int face_id = queue.front();
    queue.pop_front();
    const std::vector<int>& node_ids = faces_.at(face_id).node_ids;
    for (std::size_t i = 0; i < node_ids.size(); ++i)
    {
      const auto it =
          edge_faces_.find(edge_key(node_ids[i], node_ids[(i + 1) % node_ids.size()]));
      for (int neighbor : it->second)
      {
        if (reached.contains(neighbor) || !std::ranges::binary_search(candidates, neighbor) ||
            !same_patch(face_id, neighbor))
          continue;
        reached.insert(neighbor);
        queue.push_back(neighbor);
      }
    }
  }
  return {reached.begin(), reached.end()};
}

Core::LinAlg::Matrix<3, 3> Solid::MeshRegularization::BoundarySliding::tangent_projection(
    const int node_id, const std::array<bool, 3>& fixed_directions) const
{
  Core::LinAlg::Matrix<3, 3> projection(Core::LinAlg::Initialization::zero);
  const auto it = constraints_.find(node_id);
  if (it != constraints_.end() && it->second.is_corner()) return projection;

  for (int d = 0; d < 3; ++d) projection(d, d) = 1.0;

  // remove the span of the fixed directions, the normals of the patches and the normals of all
  // faces at the current position (Gram-Schmidt): a node at a kink between faces of a patch stays
  // there, since its minimum is usually at the kink, which steps across it would only oscillate
  // around
  std::vector<Core::LinAlg::Matrix<3, 1>> directions;
  for (int d = 0; d < 3; ++d)
  {
    if (!fixed_directions[d]) continue;
    Core::LinAlg::Matrix<3, 1> unit(Core::LinAlg::Initialization::zero);
    unit(d) = 1.0;
    directions.push_back(unit);
  }
  if (it != constraints_.end())
  {
    directions.insert(directions.end(), it->second.normals.begin(), it->second.normals.end());
    directions.insert(
        directions.end(), it->second.active_normals.begin(), it->second.active_normals.end());
  }

  std::vector<Core::LinAlg::Matrix<3, 1>> basis;
  for (Core::LinAlg::Matrix<3, 1> direction : directions)
  {
    for (const auto& b : basis) direction.update(-direction.dot(b), b, 1.0);
    if (direction.norm2() < 1.0e-6) continue;
    direction.scale(1.0 / direction.norm2());
    basis.push_back(direction);
    projection.multiply_nt(-1.0, direction, direction, 1.0);
  }
  return projection;
}

void Solid::MeshRegularization::BoundarySliding::restrict_to_tangent_space(
    Core::LinAlg::SparseMatrix& matrix, Core::LinAlg::Vector<double>& residual,
    const Core::LinAlg::Map& dirichlet_dofs) const
{
  // the projection P and its complement I - P, which are the identity and zero at interior dofs
  const Core::LinAlg::Map& dof_row_map = *discretization_.dof_row_map();
  Core::LinAlg::SparseMatrix projection(dof_row_map, 3);
  Core::LinAlg::SparseMatrix complement(dof_row_map, 3);
  for (int node_lid = 0; node_lid < discretization_.num_my_row_nodes(); ++node_lid)
  {
    const Core::Nodes::Node* node = discretization_.l_row_node(node_lid);
    const std::vector<int> dofs = discretization_.dof(0, node);
    std::array<bool, 3> fixed{};
    for (int d = 0; d < 3; ++d) fixed[d] = dirichlet_dofs.lid(dofs[d]) >= 0;

    if (!constraints_.contains(node->id()) && !std::ranges::any_of(fixed, std::identity{}))
    {
      for (int d = 0; d < 3; ++d) projection.assemble(1.0, dofs[d], dofs[d]);
      continue;
    }

    const Core::LinAlg::Matrix<3, 3> p = tangent_projection(node->id(), fixed);
    for (int a = 0; a < 3; ++a)
    {
      for (int b = 0; b < 3; ++b)
      {
        projection.assemble(p(a, b), dofs[a], dofs[b]);
        complement.assemble((a == b ? 1.0 : 0.0) - p(a, b), dofs[a], dofs[b]);
      }
    }
  }
  projection.complete();
  complement.complete();

  // P^T K P + (I - P): the increments of the boundary nodes are tangential to the boundary
  const std::unique_ptr<Core::LinAlg::SparseMatrix> stiffness_projection =
      Core::LinAlg::matrix_multiply(matrix, false, projection, false);
  std::unique_ptr<Core::LinAlg::SparseMatrix> restricted = Core::LinAlg::matrix_multiply(
      projection, true, *stiffness_projection, false, false, false, false);
  restricted->add(complement, false, 1.0, 1.0);
  restricted->complete();
  matrix.assign(Core::LinAlg::DataAccess::Copy, *restricted);

  Core::LinAlg::Vector<double> projected_residual(residual);
  projection.multiply(true, residual, projected_residual);
  residual = projected_residual;
}

void Solid::MeshRegularization::BoundarySliding::remove_constrained_components(
    Core::LinAlg::Vector<double>& residual) const
{
  for (const auto& [node_id, constraint] : constraints_)
  {
    const std::vector<int> lids = dof_lids(node_id);
    Core::LinAlg::Matrix<3, 1> r;
    for (int i = 0; i < 3; ++i) r(i) = residual.local_values_as_span()[lids[i]];

    if (constraint.is_corner())
      r.clear();
    else
    {
      // remove the span of the active normals (Gram-Schmidt)
      std::vector<Core::LinAlg::Matrix<3, 1>> basis;
      for (const auto& normal : constraint.active_normals)
      {
        Core::LinAlg::Matrix<3, 1> direction(normal);
        for (const auto& b : basis) direction.update(-direction.dot(b), b, 1.0);
        if (direction.norm2() < 1.0e-6) continue;
        direction.scale(1.0 / direction.norm2());
        basis.push_back(direction);
      }
      for (const auto& b : basis) r.update(-r.dot(b), b, 1.0);
    }

    for (int i = 0; i < 3; ++i) residual.get_values()[lids[i]] = r(i);
  }
}


Solid::MeshRegularization::BoundarySliding::ClosestPoint
Solid::MeshRegularization::BoundarySliding::closest_point(
    const Core::LinAlg::Matrix<3, 1>& point, const std::vector<int>& faces) const
{
  Core::LinAlg::Matrix<3, 1> best = point;
  int best_face = -1;
  double best_distance = std::numeric_limits<double>::max();
  std::vector<std::pair<int, Core::LinAlg::Matrix<3, 1>>> candidates;

  for (int face_id : faces)
  {
    const std::vector<Core::LinAlg::Matrix<3, 1>>& x = faces_.at(face_id).corners;

    Core::LinAlg::Matrix<3, 1> candidate;
    if (x.size() == 4)
    {
      // projected Gauss-Newton iteration on the bilinear quadrilateral
      double s = 0.0;
      double t = 0.0;
      for (int iteration = 0; iteration < 20; ++iteration)
      {
        const double N[4] = {0.25 * (1 - s) * (1 - t), 0.25 * (1 + s) * (1 - t),
            0.25 * (1 + s) * (1 + t), 0.25 * (1 - s) * (1 + t)};
        const double dNds[4] = {-0.25 * (1 - t), 0.25 * (1 - t), 0.25 * (1 + t), -0.25 * (1 + t)};
        const double dNdt[4] = {-0.25 * (1 - s), -0.25 * (1 + s), 0.25 * (1 + s), 0.25 * (1 - s)};
        Core::LinAlg::Matrix<3, 1> r(Core::LinAlg::Initialization::zero);
        Core::LinAlg::Matrix<3, 1> xs(Core::LinAlg::Initialization::zero);
        Core::LinAlg::Matrix<3, 1> xt(Core::LinAlg::Initialization::zero);
        for (int i = 0; i < 4; ++i)
        {
          r.update(N[i], x[i], 1.0);
          xs.update(dNds[i], x[i], 1.0);
          xt.update(dNdt[i], x[i], 1.0);
        }
        r.update(-1.0, point, 1.0);
        const double a11 = xs.dot(xs), a12 = xs.dot(xt), a22 = xt.dot(xt);
        const double b1 = -xs.dot(r), b2 = -xt.dot(r);
        const double det = a11 * a22 - a12 * a12;
        if (det <= 0.0) break;
        const double ds = (a22 * b1 - a12 * b2) / det;
        const double dt = (a11 * b2 - a12 * b1) / det;
        s = std::clamp(s + ds, -1.0, 1.0);
        t = std::clamp(t + dt, -1.0, 1.0);
        if (std::abs(ds) + std::abs(dt) < 1.0e-12) break;
      }
      const double N[4] = {0.25 * (1 - s) * (1 - t), 0.25 * (1 + s) * (1 - t),
          0.25 * (1 + s) * (1 + t), 0.25 * (1 - s) * (1 + t)};
      candidate.clear();
      for (int i = 0; i < 4; ++i) candidate.update(N[i], x[i], 1.0);
    }
    else
    {
      // projection onto the plane of the triangle with clamped barycentric coordinates
      Core::LinAlg::Matrix<3, 1> e1, e2, r;
      e1.update(1.0, x[1], -1.0, x[0], 0.0);
      e2.update(1.0, x[2], -1.0, x[0], 0.0);
      r.update(1.0, point, -1.0, x[0], 0.0);
      const double a11 = e1.dot(e1), a12 = e1.dot(e2), a22 = e2.dot(e2);
      const double det = a11 * a22 - a12 * a12;
      double u = (a22 * e1.dot(r) - a12 * e2.dot(r)) / det;
      double v = (a11 * e2.dot(r) - a12 * e1.dot(r)) / det;
      u = std::max(u, 0.0);
      v = std::max(v, 0.0);
      if (u + v > 1.0)
      {
        u /= (u + v);
        v = 1.0 - u;
      }
      candidate = x[0];
      candidate.update(u, e1, v, e2, 1.0);
    }

    candidates.emplace_back(face_id, candidate);
    Core::LinAlg::Matrix<3, 1> difference(candidate);
    difference.update(-1.0, point, 1.0);
    if (difference.norm2() < best_distance)
    {
      best_distance = difference.norm2();
      best = candidate;
      best_face = face_id;
    }
  }

  // all faces whose closest point coincides with the closest point; where the boundary has a kink
  // between faces, the minimum of a sliding node is often at the kink, which the Newton iteration
  // only approaches asymptotically (steps across it are halved), so the tolerance has to be well
  // above round-off for both normals to balance the residual there
  ClosestPoint result{best, best_face, {}};
  if (best_face < 0) return result;
  Core::LinAlg::Matrix<3, 1> diagonal = faces_.at(best_face).corners[0];
  diagonal.update(-1.0, faces_.at(best_face).corners[2], 1.0);
  const double tolerance = 1.0e-6 * diagonal.norm2();
  for (const auto& [face_id, candidate] : candidates)
  {
    Core::LinAlg::Matrix<3, 1> difference(candidate);
    difference.update(-1.0, best, 1.0);
    if (difference.norm2() <= tolerance) result.touching_faces.push_back(face_id);
  }
  return result;
}


void Solid::MeshRegularization::BoundarySliding::project_onto_boundary(
    Core::LinAlg::Vector<double>& displacement, const Core::LinAlg::Map& dirichlet_dofs)
{
  // candidate faces near the moved nodes: a boundary node starts on the boundary, so its closest
  // point lies within the distance of its displacement, and the bounding boxes of all faces at the
  // closest point contain it; the faces of all processors are searched, those not known on this
  // processor are fetched from their owners
  std::vector<int> searched_nodes;
  std::vector<std::pair<int, Core::GeometricSearch::BoundingVolume>> node_boxes;
  for (const auto& [node_id, constraint] : constraints_)
  {
    if (constraint.is_corner()) continue;
    const std::vector<int> lids = dof_lids(node_id);
    Core::LinAlg::Matrix<3, 1> point = position(node_id);
    for (int d = 0; d < 3; ++d) point(d) += displacement.local_values_as_span()[lids[d]];
    Core::LinAlg::Matrix<3, 1> movement(point);
    movement.update(-1.0, position(node_id), 1.0);
    Core::GeometricSearch::BoundingVolume box;
    box.add_point(point);
    box.extend_boundaries(movement.norm2() + bounding_box_tolerance);
    node_boxes.emplace_back(static_cast<int>(searched_nodes.size()), box);
    searched_nodes.push_back(node_id);
  }

  std::map<int, std::vector<int>> node_candidates;
  std::map<int, int> missing_faces;
  for (const auto& collision : Core::GeometricSearch::global_collision_search(
           face_boxes_, node_boxes, discretization_.get_comm()))
  {
    node_candidates[searched_nodes[collision.lid_predicate]].push_back(collision.gid_primitive);
    if (!faces_.contains(collision.gid_primitive))
      missing_faces[collision.gid_primitive] = collision.pid_primitive;
  }
  fetch_faces(missing_faces);

  for (auto& [node_id, constraint] : constraints_)
  {
    const std::vector<int> lids = dof_lids(node_id);
    const std::vector<int> dofs = discretization_.dof(0, discretization_.g_node(node_id));
    const Core::LinAlg::Matrix<3, 1> reference = position(node_id);

    // the Dirichlet dofs keep their values exactly
    const auto set_displacement = [&](const Core::LinAlg::Matrix<3, 1>& value)
    {
      for (int d = 0; d < 3; ++d)
        if (dirichlet_dofs.lid(dofs[d]) < 0) displacement.get_values()[lids[d]] = value(d);
    };

    if (constraint.is_corner())
    {
      set_displacement(Core::LinAlg::Matrix<3, 1>(Core::LinAlg::Initialization::zero));
      continue;
    }

    // the candidate faces of each patch of the node
    std::vector<int>& node_faces = node_candidates[node_id];
    std::ranges::sort(node_faces);
    const std::size_t num_patches = constraint.patch_faces.size();
    std::vector<std::vector<int>> candidates(num_patches);
    for (std::size_t g = 0; g < num_patches; ++g)
      candidates[g] = patch_candidates(constraint.patch_faces[g], node_faces);

    Core::LinAlg::Matrix<3, 1> point(reference);
    for (int d = 0; d < 3; ++d) point(d) += displacement.local_values_as_span()[lids[d]];

    const int num_iterations = num_patches == 1 ? 1 : num_edge_projections;
    for (int iteration = 0; iteration < num_iterations; ++iteration)
    {
      for (std::size_t g = 0; g < num_patches; ++g)
      {
        const ClosestPoint closest = closest_point(point, candidates[g]);
        point = closest.point;
        // the next increment is tangential to the face at the projected position
        if (closest.face >= 0) constraint.normals[g] = faces_.at(closest.face).normal;
      }
    }

    // the normals of all faces at the projected position balance the potential at the minimum
    constraint.active_normals.clear();
    for (const auto& faces : candidates)
      for (int face_id : closest_point(point, faces).touching_faces)
        constraint.active_normals.push_back(faces_.at(face_id).normal);

    point.update(-1.0, reference, 1.0);
    set_displacement(point);
  }
}

FOUR_C_NAMESPACE_CLOSE
