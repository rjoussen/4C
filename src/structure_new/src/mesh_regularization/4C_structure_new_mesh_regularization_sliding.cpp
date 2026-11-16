// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_structure_new_mesh_regularization_sliding.hpp"

#include "4C_adapter_coupling_nonlin_mortar.hpp"
#include "4C_comm_mpi_utils.hpp"
#include "4C_contact_element.hpp"
#include "4C_contact_input.hpp"
#include "4C_contact_interface.hpp"
#include "4C_contact_node.hpp"
#include "4C_fem_discretization.hpp"
#include "4C_fem_general_element.hpp"
#include "4C_fem_general_node.hpp"
#include "4C_fem_general_utils_local_connectivity_matrices.hpp"
#include "4C_global_data.hpp"
#include "4C_linalg_sparsematrix.hpp"
#include "4C_linalg_utils_sparse_algebra_manipulation.hpp"
#include "4C_linalg_utils_sparse_algebra_math.hpp"
#include "4C_mortar_element.hpp"
#include "4C_mortar_input.hpp"
#include "4C_utils_exceptions.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numbers>
#include <set>
#include <string>

FOUR_C_NAMESPACE_OPEN

namespace
{
  //! largest number of faces of a 3D element, which numbers the faces globally
  constexpr int max_faces_per_element = 6;

  //! squared length below which the part of a normal in the free directions of a node vanishes
  constexpr double free_normal_tolerance = 1.0e-12;

  //! number of corner nodes of a face with the given shape
  int num_corners(const Core::FE::CellType shape)
  {
    switch (shape)
    {
      case Core::FE::CellType::tri3:
      case Core::FE::CellType::tri6:
        return 3;
      case Core::FE::CellType::quad4:
      case Core::FE::CellType::quad8:
      case Core::FE::CellType::quad9:
        return 4;
      default:
        FOUR_C_THROW(
            "Unsupported boundary face of shape {}.", Core::FE::cell_type_to_string(shape));
    }
  }

  //! node ids of a face with the opposite orientation (corners, then edge midpoints, then center)
  std::vector<int> reversed(const std::vector<int>& node_ids, const Core::FE::CellType shape)
  {
    const int n = num_corners(shape);
    std::vector<int> result;
    // corners in reversed order, starting at the first one
    for (int i = 0; i < n; ++i) result.push_back(node_ids[(n - i) % n]);
    // the midpoint of the edge from corner i to corner i + 1 becomes the one of edge n - 1 - i
    if (static_cast<int>(node_ids.size()) >= 2 * n)
      for (int i = 0; i < n; ++i) result.push_back(node_ids[n + (2 * n - 1 - i) % n]);
    if (static_cast<int>(node_ids.size()) == 2 * n + 1) result.push_back(node_ids.back());
    FOUR_C_ASSERT_ALWAYS(
        result.size() == node_ids.size(), "Unsupported face with {} nodes.", node_ids.size());
    return result;
  }

  //! closest point on a triangle to a point (Ericson, Real-Time Collision Detection, 5.1.5)
  std::array<double, 3> face_normal(const std::vector<std::array<double, 3>>& x)
  {
    std::array<double, 3> a;
    std::array<double, 3> b;
    for (int d = 0; d < 3; ++d)
    {
      a[d] = x.size() == 4 ? x[2][d] - x[0][d] : x[1][d] - x[0][d];
      b[d] = x.size() == 4 ? x[3][d] - x[1][d] : x[2][d] - x[0][d];
    }
    std::array<double, 3> n = {
        a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]};
    const double length = std::sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
    for (double& value : n) value /= length;
    return n;
  }

  /*!
   * @brief Parameters of the mortar couplings of the sliding, independent of the settings of other
   * contact or mortar problems: dual Lagrange multipliers, segment-based mortar integration, and
   * frictionless constraints in Cartesian directions, such that fixed directions can be removed
   * from the normals
   */
  Teuchos::ParameterList sliding_parameters(const Global::Problem& problem)
  {
    Teuchos::ParameterList parameters;
    parameters.setParameters(problem.mortar_coupling_params());
    parameters.setParameters(problem.contact_dynamic_params());
    parameters.setParameters(problem.wear_params());
    parameters.set<Mortar::ShapeFcn>("LM_SHAPEFCN", Mortar::shape_dual);
    parameters.set<Mortar::AlgorithmType>("ALGORITHM", Mortar::algorithm_mortar);
    parameters.set<Mortar::IntType>("INTTYPE", Mortar::inttype_segments);
    parameters.set<bool>("CROSSPOINTS", false);
    // the segments of matching slave and master faces are cocircular, which the Delaunay
    // triangulation does not handle robustly
    parameters.set<Mortar::Triangulation>("TRIANGULATION", Mortar::triangulation_center);
    parameters.set<CONTACT::ConstraintDirection>(
        "CONSTRAINT_DIRECTIONS", CONTACT::ConstraintDirection::xyz);
    parameters.set<CONTACT::FrictionType>("FRICTION", CONTACT::FrictionType::none);
    parameters.set<bool>("NONSMOOTH_GEOMETRIES", false);
    parameters.set<CONTACT::Problemtype>("PROBTYPE", CONTACT::Problemtype::other);
    parameters.set<bool>("NURBS", false);
    parameters.set<int>("DIMENSION", problem.n_dim());
    parameters.set<bool>("Two_half_pass", false);
    // the interfaces are small and rebuilt for each mesh movement
    parameters.sublist("PARALLEL REDISTRIBUTION")
        .set<Mortar::ParallelRedist>("PARALLEL_REDIST", Mortar::ParallelRedist::redist_none);
    return parameters;
  }

  /*!
   * @brief Mortar coupling of the nodes of a patch (slave) to the fixed copy of its faces (master)
   *
   * The coupled nodes and elements are given instead of being read from conditions, the parameters
   * are the ones of sliding_parameters().
   */
  class PatchCoupling : public Adapter::CouplingNonLinMortar
  {
   public:
    PatchCoupling(Global::Problem& problem, std::map<int, Core::Nodes::Node*> master_nodes,
        std::map<int, std::shared_ptr<Core::Elements::Element>> master_elements,
        std::map<int, Core::Nodes::Node*> slave_nodes,
        std::map<int, std::shared_ptr<Core::Elements::Element>> slave_elements)
        : Adapter::CouplingNonLinMortar(problem, problem.n_dim(), sliding_parameters(problem),
              sliding_parameters(problem), problem.spatial_approximation_type()),
          parameters_(sliding_parameters(problem)),
          master_nodes_(std::move(master_nodes)),
          master_elements_(std::move(master_elements)),
          slave_nodes_(std::move(slave_nodes)),
          slave_elements_(std::move(slave_elements))
    {
      print_parallel_distribution_ = false;
    }

   protected:
    void read_mortar_condition(std::shared_ptr<Core::FE::Discretization>,
        std::shared_ptr<Core::FE::Discretization>, std::vector<int>, const std::string&,
        Teuchos::ParameterList& input, std::map<int, Core::Nodes::Node*>& target_global_nodes,
        std::map<int, Core::Nodes::Node*>& source_global_nodes,
        std::map<int, std::shared_ptr<Core::Elements::Element>>& target_elements,
        std::map<int, std::shared_ptr<Core::Elements::Element>>& source_elements) override
    {
      input.setParameters(parameters_);
      target_global_nodes = master_nodes_;
      target_elements = master_elements_;
      source_global_nodes = slave_nodes_;
      source_elements = slave_elements_;
    }

   private:
    Teuchos::ParameterList parameters_;
    std::map<int, Core::Nodes::Node*> master_nodes_;
    std::map<int, std::shared_ptr<Core::Elements::Element>> master_elements_;
    std::map<int, Core::Nodes::Node*> slave_nodes_;
    std::map<int, std::shared_ptr<Core::Elements::Element>> slave_elements_;
  };

  //! add the entries of a matrix in the columns of the mesh movement (dropping the fixed boundary)
  void add_mesh_movement_columns(const Core::LinAlg::SparseMatrix& source,
      Core::LinAlg::SparseMatrix& target, const int max_mesh_movement_dof)
  {
    for (int lid = 0; lid < source.row_map().num_my_elements(); ++lid)
    {
      int num_entries = 0;
      double* values = nullptr;
      int* indices = nullptr;
      source.extract_my_row_view(lid, num_entries, values, indices);
      const int row = source.row_map().gid(lid);
      for (int i = 0; i < num_entries; ++i)
      {
        const int column = source.col_map().gid(indices[i]);
        if (column <= max_mesh_movement_dof) target.assemble(values[i], row, column);
      }
    }
  }

  //! add the entries of a vector at the dofs of the mesh movement
  void add_mesh_movement_values(
      const Core::LinAlg::Vector<double>& source, Core::LinAlg::Vector<double>& target)
  {
    for (int lid = 0; lid < source.get_map().num_my_elements(); ++lid)
    {
      const int target_lid = target.get_map().lid(source.get_map().gid(lid));
      FOUR_C_ASSERT_ALWAYS(target_lid >= 0, "Slave dof {} is not a row dof of the mesh movement.",
          source.get_map().gid(lid));
      target.get_values()[target_lid] += source.local_values_as_span()[lid];
    }
  }
}  // namespace

Solid::MeshRegularization::BoundarySliding::BoundarySliding(Global::Problem& problem,
    std::shared_ptr<Core::FE::Discretization> discretization, const double feature_angle)
    : problem_(problem),
      discretization_(std::move(discretization)),
      feature_angle_cosine_(std::cos(feature_angle * std::numbers::pi / 180.0))
{
  // count the occurrences of the faces of the column elements: all elements adjacent to a row node
  // are column elements, so the faces with a row node are classified correctly, and the boundary
  // faces among them belong to one element only, which numbers them globally
  struct FaceCount
  {
    int count = 0;
    int id = -1;
    Face face;
  };
  std::map<std::vector<int>, FaceCount> face_count;
  for (int lid = 0; lid < discretization_->num_my_col_elements(); ++lid)
  {
    const Core::Elements::Element* element = discretization_->l_col_element(lid);
    const std::vector<std::vector<int>> surfaces =
        Core::FE::get_ele_node_numbering_surfaces(element->shape());
    FOUR_C_ASSERT_ALWAYS(static_cast<int>(surfaces.size()) <= max_faces_per_element,
        "Element {} has more than {} faces.", element->id(), max_faces_per_element);
    for (std::size_t s = 0; s < surfaces.size(); ++s)
    {
      Face face{.node_ids = {},
          .shape = Core::FE::get_ele_face_shape_type(element->shape(), static_cast<unsigned>(s))};
      for (int i : surfaces[s]) face.node_ids.push_back(element->node_ids()[i]);
      std::vector<int> key = face.node_ids;
      std::ranges::sort(key);
      FaceCount& entry = face_count[key];
      ++entry.count;
      entry.id = element->id() * max_faces_per_element + static_cast<int>(s);
      entry.face = face;
    }
  }

  const Core::LinAlg::Map& node_row_map = *discretization_->node_row_map();
  for (const auto& [key, entry] : face_count)
  {
    if (entry.count != 1) continue;
    if (std::ranges::none_of(key, [&](int id) { return node_row_map.lid(id) >= 0; })) continue;
    faces_[entry.id] = entry.face;
    for (int id : entry.face.node_ids)
      if (node_row_map.lid(id) >= 0) node_faces_[id].push_back(entry.id);
  }

  group_into_patches();
}

void Solid::MeshRegularization::BoundarySliding::group_into_patches()
{
  // adjacent faces whose normals in the current configuration enclose an angle below the feature
  // angle belong to the same patch: each processor finds the pairs at its row nodes, the patches
  // are the connected components of all pairs
  std::vector<std::pair<int, int>> pairs;
  for (const auto& [node_id, face_ids] : node_faces_)
  {
    for (std::size_t i = 0; i < face_ids.size(); ++i)
    {
      for (std::size_t j = i + 1; j < face_ids.size(); ++j)
      {
        const Face& a = faces_.at(face_ids[i]);
        const Face& b = faces_.at(face_ids[j]);
        const auto shared = std::ranges::count_if(a.node_ids,
            [&](int id) { return std::ranges::find(b.node_ids, id) != b.node_ids.end(); });
        if (shared < 2) continue;
        const std::array<double, 3> n_a = normal(a);
        const std::array<double, 3> n_b = normal(b);
        if (std::abs(n_a[0] * n_b[0] + n_a[1] * n_b[1] + n_a[2] * n_b[2]) >= feature_angle_cosine_)
          pairs.emplace_back(face_ids[i], face_ids[j]);
      }
    }
  }
  pairs = Core::Communication::all_reduce(pairs, discretization_->get_comm());

  // union-find with the smallest face id as the root of a patch
  std::map<int, int> parent;
  const auto root = [&](int id)
  {
    while (parent.contains(id) && parent.at(id) != id) id = parent.at(id);
    return id;
  };
  for (const auto& [a, b] : pairs)
  {
    const int root_a = root(a);
    const int root_b = root(b);
    parent[std::max(root_a, root_b)] = std::min(root_a, root_b);
  }
  patch_of_face_.clear();
  for (const auto& [id, face] : faces_) patch_of_face_[id] = root(id);
}

Solid::MeshRegularization::BoundarySliding::~BoundarySliding() = default;

std::array<double, 3> Solid::MeshRegularization::BoundarySliding::normal(const Face& face) const
{
  std::vector<std::array<double, 3>> corners;
  for (int i = 0; i < num_corners(face.shape); ++i)
  {
    const Core::Nodes::Node* node = discretization_->g_node(face.node_ids[i]);
    corners.push_back({node->x()[0], node->x()[1], node->x()[2]});
  }
  return face_normal(corners);
}

std::vector<int> Solid::MeshRegularization::BoundarySliding::patches(const int node_id) const
{
  std::set<int> node_patches;
  if (node_faces_.contains(node_id))
    for (int face_id : node_faces_.at(node_id)) node_patches.insert(patch_of_face_.at(face_id));
  return {node_patches.begin(), node_patches.end()};
}

void Solid::MeshRegularization::BoundarySliding::setup(const Core::LinAlg::Map& dirichlet_dofs)
{
  couplings_.clear();
  fixed_boundary_.reset();

  MPI_Comm comm = discretization_->get_comm();
  const int my_rank = Core::Communication::my_mpi_rank(comm);
  const Core::LinAlg::Map& dof_row_map = *discretization_->dof_row_map();
  const auto is_fixed = [&](const Core::Nodes::Node* node)
  {
    const std::vector<int> dofs = discretization_->dof(0, node);
    std::array<bool, 3> fixed{};
    for (int d = 0; d < 3; ++d) fixed[d] = dirichlet_dofs.lid(dofs[d]) >= 0;
    return fixed;
  };
  const auto free_normal_part = [](const std::array<double, 3>& n, const std::array<bool, 3>& fixed)
  {
    double squared = 0.0;
    for (int d = 0; d < 3; ++d)
      if (!fixed[d]) squared += n[d] * n[d];
    return squared;
  };

  // a patch is coupled unless its normal directions are fixed at all of its nodes
  std::set<int> coupled_patches;
  for (const auto& [node_id, face_ids] : node_faces_)
  {
    const std::array<bool, 3> fixed = is_fixed(discretization_->g_node(node_id));
    for (int face_id : face_ids)
    {
      if (free_normal_part(normal(faces_.at(face_id)), fixed) > free_normal_tolerance)
        coupled_patches.insert(patch_of_face_.at(face_id));
    }
  }
  coupled_patches = Core::Communication::all_reduce(coupled_patches, comm);

  // corners on several coupled patches are fixed; at nodes on one coupled patch, the fixed
  // directions are removed from the normal unless it lies in the fixed directions or all directions
  // are fixed
  std::vector<int> fixed_dofs;
  Core::LinAlg::Vector<double> symmetry_dofs(dof_row_map, true);
  for (int lid = 0; lid < discretization_->num_my_row_nodes(); ++lid)
  {
    const Core::Nodes::Node* node = discretization_->l_row_node(lid);
    const std::vector<int> dofs = discretization_->dof(0, node);
    const std::array<bool, 3> fixed = is_fixed(node);
    std::vector<int> node_patches;
    for (int patch : patches(node->id()))
      if (coupled_patches.contains(patch)) node_patches.push_back(patch);

    if (node_patches.size() > 1)
    {
      fixed_dofs.insert(fixed_dofs.end(), dofs.begin(), dofs.begin() + 3);
      continue;
    }
    for (int d = 0; d < 3; ++d)
      if (fixed[d]) fixed_dofs.push_back(dofs[d]);

    const int num_fixed = static_cast<int>(std::ranges::count(fixed, true));
    if (node_patches.size() != 1 || num_fixed == 0 || num_fixed == 3) continue;
    std::array<double, 3> patch_normal{};
    std::array<double, 3> first{};
    for (int face_id : node_faces_.at(node->id()))
    {
      if (patch_of_face_.at(face_id) != node_patches.front()) continue;
      const std::array<double, 3> n = normal(faces_.at(face_id));
      if (first == std::array<double, 3>{}) first = n;
      const double sign = n[0] * first[0] + n[1] * first[1] + n[2] * first[2] >= 0.0 ? 1.0 : -1.0;
      for (int d = 0; d < 3; ++d) patch_normal[d] += sign * n[d];
    }
    if (free_normal_part(patch_normal, fixed) > free_normal_tolerance)
      for (int d = 0; d < 3; ++d)
        if (fixed[d]) symmetry_dofs.get_values()[dof_row_map.lid(dofs[d])] = 1.0;
  }
  {
    int num_corners = 0;
    for (int lid = 0; lid < discretization_->num_my_row_nodes(); ++lid)
    {
      std::set<int> coupled;
      for (int patch : patches(discretization_->l_row_node(lid)->id()))
        if (coupled_patches.contains(patch)) coupled.insert(patch);
      if (coupled.size() > 1) ++num_corners;
    }
    num_corners = Core::Communication::sum_all(num_corners, comm);
    if (my_rank == 0)
      std::cout << "Boundary sliding: " << num_corners << " corner nodes are fixed" << std::endl;
  }
  fixed_dofs_ = std::make_shared<Core::LinAlg::Map>(
      -1, static_cast<int>(fixed_dofs.size()), fixed_dofs.data(), 0, comm);

  // the fixed copy of the faces of the coupled patches: the faces with a row node, owned by the
  // processor owning their smallest node, with offset node and face ids; the copies face the
  // slave faces as in contact
  const int node_offset = discretization_->node_row_map()->max_all_gid() + 1;
  const int face_offset =
      max_faces_per_element * (discretization_->element_row_map()->max_all_gid() + 1);
  fixed_boundary_ = std::make_shared<Core::FE::Discretization>(
      "mesh regularization fixed boundary", comm, discretization_->n_dim());
  std::map<int, std::vector<int>> patch_faces;
  std::map<int, std::shared_ptr<Core::Elements::Element>> copied_faces;
  for (const auto& [face_id, face] : faces_)
  {
    if (!coupled_patches.contains(patch_of_face_.at(face_id))) continue;
    patch_faces[patch_of_face_.at(face_id)].push_back(face_id);

    std::vector<int> copied_ids;
    for (int id : reversed(face.node_ids, face.shape))
    {
      copied_ids.push_back(id + node_offset);
      if (fixed_boundary_->have_global_node(id + node_offset)) continue;
      const Core::Nodes::Node* node = discretization_->g_node(id);
      fixed_boundary_->add_node(node->x(), id + node_offset,
          std::make_shared<CONTACT::Node>(
              id + node_offset, node->x(), node->owner(), std::vector<int>(3, -1), false, false));
    }
    const int owner = discretization_->g_node(std::ranges::min(face.node_ids))->owner();
    copied_faces[face_id] = std::make_shared<CONTACT::Element>(face_id + face_offset, owner,
        face.shape, static_cast<int>(copied_ids.size()), copied_ids.data(), false);
    fixed_boundary_->add_element(copied_faces.at(face_id));
  }
  fixed_boundary_->fill_complete(
      {.assign_degrees_of_freedom = true, .init_elements = false, .do_boundary_conditions = false});
  const int max_mesh_movement_dof = dof_row_map.max_all_gid();
  FOUR_C_ASSERT_ALWAYS(fixed_boundary_->dof_row_map()->num_global_elements() == 0 ||
                           fixed_boundary_->dof_row_map()->min_all_gid() > max_mesh_movement_dof,
      "The dofs of the fixed boundary have to follow the dofs of the mesh movement.");

  // one mortar coupling per coupled patch, in the same order on all processors
  for (int patch : coupled_patches)
  {
    std::map<int, Core::Nodes::Node*> master_nodes;
    std::map<int, std::shared_ptr<Core::Elements::Element>> master_elements;
    std::map<int, Core::Nodes::Node*> slave_nodes;
    std::map<int, std::shared_ptr<Core::Elements::Element>> slave_elements;
    if (patch_faces.contains(patch))
    {
      for (int face_id : patch_faces.at(patch))
      {
        const Face& face = faces_.at(face_id);
        for (int id : face.node_ids)
        {
          slave_nodes[id] = discretization_->g_node(id);
          master_nodes[id + node_offset] = fixed_boundary_->g_node(id + node_offset);
        }
        const int owner = discretization_->g_node(std::ranges::min(face.node_ids))->owner();
        master_elements[face_id + face_offset] = copied_faces.at(face_id);
        slave_elements[face_id] = std::make_shared<CONTACT::Element>(face_id, owner, face.shape,
            static_cast<int>(face.node_ids.size()), face.node_ids.data(), true);
      }
    }

    auto coupling = std::make_unique<PatchCoupling>(problem_, std::move(master_nodes),
        std::move(master_elements), std::move(slave_nodes), std::move(slave_elements));
    coupling->setup(fixed_boundary_, discretization_, {1, 1, 1}, "");
    search_inflation_factor_ =
        coupling->interface()->geometric_search_params().bounding_volume_inflation_factor;

    // mortar symmetry conditions: the fixed directions are removed from the normals (the interface
    // may be distributed differently)
    Core::LinAlg::Vector<double> symmetry_dofs_col(*coupling->interface()->source_col_dofs());
    Core::LinAlg::export_to(symmetry_dofs, symmetry_dofs_col);
    const Core::FE::Discretization& interface_discretization = coupling->interface()->discret();
    for (int lid = 0; lid < interface_discretization.num_my_col_nodes(); ++lid)
    {
      auto* node = dynamic_cast<CONTACT::Node*>(interface_discretization.l_col_node(lid));
      if (!node->is_source()) continue;
      for (int d = 0; d < 3; ++d)
      {
        const int symmetry_lid = symmetry_dofs_col.get_map().lid(node->dofs()[d]);
        FOUR_C_ASSERT_ALWAYS(symmetry_lid >= 0,
            "Slave node {} of the boundary sliding has no column dofs on processor {}.", node->id(),
            my_rank);
        node->dbc_dofs()[d] = symmetry_dofs_col.local_values_as_span()[symmetry_lid] == 1.0;
      }
    }
    couplings_.push_back(std::move(coupling));
  }

  state_map_ = Core::LinAlg::merge_map(dof_row_map, *fixed_boundary_->dof_row_map(), false);
  initial_gap_.reset();
}

void Solid::MeshRegularization::BoundarySliding::condense(Core::LinAlg::SparseMatrix& matrix,
    Core::LinAlg::Vector<double>& residual, const Core::LinAlg::Vector<double>& displacement)
{
  const Core::LinAlg::Map& dof_row_map = *discretization_->dof_row_map();
  const int max_mesh_movement_dof = dof_row_map.max_all_gid();

  // the states on the dofs of the mesh movement and of the fixed boundary (zero): the multipliers
  // are the ones of the current displacement, D z = -r_s at the slave dofs, which balance the
  // forces of the potential (exact at the solution, where these forces are normal); the fixed
  // directions are included, since the tangential equilibrium -T r_s contains their reactions
  // where a tangent tilts into them
  auto displacement_state = std::make_shared<Core::LinAlg::Vector<double>>(*state_map_, true);
  auto multiplier_state = std::make_shared<Core::LinAlg::Vector<double>>(*state_map_, true);
  for (int lid = 0; lid < dof_row_map.num_my_elements(); ++lid)
  {
    const int state_lid = state_map_->lid(dof_row_map.gid(lid));
    displacement_state->get_values()[state_lid] = displacement.local_values_as_span()[lid];
    multiplier_state->get_values()[state_lid] = -residual.local_values_as_span()[lid];
  }

  // the mortar matrices of all couplings in the dofs of the mesh movement
  Core::LinAlg::SparseMatrix gap_linearization(dof_row_map, 81, true, true);
  Core::LinAlg::SparseMatrix tangents(dof_row_map, 3, true, true);
  Core::LinAlg::SparseMatrix tangent_linearization(dof_row_map, 81, true, true);
  Core::LinAlg::Vector<double> weights(dof_row_map, true);
  Core::LinAlg::Vector<double> gap(dof_row_map, true);
  // the slave faces may move far from the fixed boundary within the mesh movement, in particular
  // when sliding along a curved boundary: a slave face lies within the largest movement of the
  // slave nodes from the faces it overlaps, so the search boxes are enlarged by this movement
  double largest_movement = 0.0;
  for (const auto& [node_id, face_ids] : node_faces_)
  {
    const std::vector<int> dofs = discretization_->dof(0, discretization_->g_node(node_id));
    double squared = 0.0;
    for (int d = 0; d < 3; ++d)
    {
      const double u = displacement.local_values_as_span()[dof_row_map.lid(dofs[d])];
      squared += u * u;
    }
    largest_movement = std::max(largest_movement, std::sqrt(squared));
  }
  largest_movement = Core::Communication::max_all(largest_movement, discretization_->get_comm());

  // a slave face that does not overlap the fixed boundary has no constraint
  bool uncovered = false;
  for (const auto& coupling : couplings_)
  {
    // the search inflates the boxes by the factor times the shortest edge of the interface
    const Core::FE::Discretization& interface_discretization = coupling->interface()->discret();
    double shortest_edge = std::numeric_limits<double>::max();
    for (int lid = 0; lid < interface_discretization.num_my_col_elements(); ++lid)
    {
      shortest_edge = std::min(shortest_edge,
          dynamic_cast<const Mortar::Element*>(interface_discretization.l_col_element(lid))
              ->min_edge_size());
    }
    if (shortest_edge < std::numeric_limits<double>::max())
    {
      coupling->interface()->geometric_search_params().bounding_volume_inflation_factor =
          search_inflation_factor_ + largest_movement / shortest_edge;
    }

    coupling->evaluate_sliding("displacement", displacement_state, multiplier_state);
    add_mesh_movement_columns(*coupling->n_matrix(), gap_linearization, max_mesh_movement_dof);
    add_mesh_movement_columns(*coupling->t_matrix(), tangents, max_mesh_movement_dof);
    add_mesh_movement_columns(*coupling->h_matrix(), tangent_linearization, max_mesh_movement_dof);

    const Core::LinAlg::SparseMatrix& d = *coupling->get_mortar_matrix_d();
    Core::LinAlg::Vector<double> d_diagonal(d.row_map(), true);
    d.extract_diagonal_copy(d_diagonal);
    for (int lid = 0; lid < d_diagonal.local_length(); ++lid)
    {
      // fixed dofs need no constraint, e.g., at a corner
      if (d_diagonal.local_values_as_span()[lid] > 0.0 ||
          fixed_dofs_->lid(d_diagonal.get_map().gid(lid)) >= 0)
        continue;
      uncovered = true;
    }
    add_mesh_movement_values(d_diagonal, weights);
    add_mesh_movement_values(*coupling->gap(), gap);
  }
  // such a displacement is not admissible: the residual is not a number, so that a Newton step to
  // it is halved
  if (Core::Communication::sum_all(uncovered ? 1 : 0, discretization_->get_comm()) > 0)
  {
    residual.put_scalar(std::numeric_limits<double>::quiet_NaN());
    return;
  }
  gap_linearization.complete();
  tangents.complete();
  tangent_linearization.complete();

  // the slave rows are replaced by -T (r_s + K_s du) + H(D z) du + N du + g = 0, where H(D z) is
  // the linearization of the tangents T contracted with the multipliers D z = -r_s
  Core::LinAlg::SparseMatrix row_operator(dof_row_map, 3, true, true);
  for (int lid = 0; lid < dof_row_map.num_my_elements(); ++lid)
    if (weights.local_values_as_span()[lid] == 0.0)
      row_operator.assemble(1.0, dof_row_map.gid(lid), dof_row_map.gid(lid));
  row_operator.add(tangents, false, -1.0, 1.0);
  row_operator.complete();

  std::unique_ptr<Core::LinAlg::SparseMatrix> condensed =
      Core::LinAlg::matrix_multiply(row_operator, false, matrix, false, false, false, true);
  Core::LinAlg::SparseMatrix system(dof_row_map, 81, true, true);
  system.add(*condensed, false, 1.0, 0.0);
  system.add(tangent_linearization, false, 1.0, 1.0);
  system.add(gap_linearization, false, 1.0, 1.0);
  system.complete();
  system.apply_dirichlet(*fixed_dofs_, true);
  matrix.assign(Core::LinAlg::DataAccess::Copy, system);

  Core::LinAlg::Vector<double> condensed_residual(dof_row_map, true);
  row_operator.multiply(false, residual, condensed_residual);
  // the weighted gaps keep their values of the start of the mesh movement (zero up to round-off)
  if (!initial_gap_) initial_gap_ = std::make_unique<Core::LinAlg::Vector<double>>(gap);
  condensed_residual.update(1.0, gap, -1.0, *initial_gap_, 1.0);
  for (int lid = 0; lid < fixed_dofs_->num_my_elements(); ++lid)
  {
    const int row = dof_row_map.lid(fixed_dofs_->gid(lid));
    condensed_residual.get_values()[row] = displacement.local_values_as_span()[row];
  }
  residual = condensed_residual;
}

FOUR_C_NAMESPACE_CLOSE
