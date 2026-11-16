// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_STRUCTURE_NEW_MESH_REGULARIZATION_SLIDING_HPP
#define FOUR_C_STRUCTURE_NEW_MESH_REGULARIZATION_SLIDING_HPP

#include "4C_config.hpp"

#include "4C_geometric_search_bounding_volume.hpp"
#include "4C_linalg_fixedsizematrix.hpp"
#include "4C_linalg_vector.hpp"

#include <array>
#include <map>
#include <span>
#include <utility>
#include <vector>

FOUR_C_NAMESPACE_OPEN

namespace Core::FE
{
  class Discretization;
}

namespace Core::LinAlg
{
  class Map;
  class SparseMatrix;
}  // namespace Core::LinAlg

namespace Solid::MeshRegularization
{
  /*!
   * @brief Sliding of boundary nodes along the boundary during the mesh movement
   *
   * The boundary faces of the discretization are detected automatically and grouped into smooth
   * patches once, in the initial configuration: adjacent faces whose normals enclose an angle
   * below the feature angle belong to the same patch. Since the connectivity does not change, the
   * patches are kept for all later configurations, such that kinks created by the deformation are
   * not mistaken for features. A node on one patch slides on the surface, a node on two patches
   * slides along their common feature edge, and a node on more patches is a corner and is fixed.
   *
   * The Newton steps of the mesh movement are restricted to the tangent planes at the current
   * positions of the nodes on the boundary, and the nodes are projected onto the boundary of the
   * reference configuration after each step, such that the shape of the domain is preserved.
   *
   * Each processor knows the boundary faces with a row node and owns those whose smallest node it
   * owns. The faces near a moved node are searched on all processors and fetched from their owners
   * if unknown, so a node may slide arbitrarily far within a mesh movement.
   */
  class BoundarySliding
  {
   public:
    /*!
     * @brief Detect the boundary faces and their patches in the current configuration of the
     * discretization, which is the initial one
     *
     * @param discretization Discretization of the mesh movement, whose column elements contain all
     * elements adjacent to the row nodes
     * @param feature_angle Angle in degrees above which adjacent boundary faces belong to different
     * patches
     */
    BoundarySliding(const Core::FE::Discretization& discretization, double feature_angle);

    /*!
     * @brief Set up the constraints of the boundary nodes in the current reference configuration
     * of the discretization (before each mesh movement)
     */
    void setup();

    /*!
     * @brief Projection of the movement of a node onto the tangent space of the boundary
     *
     * The projection removes the normals of the patches of the node and the given fixed
     * directions. It is the identity at interior nodes and zero at corners.
     */
    [[nodiscard]] Core::LinAlg::Matrix<3, 3> tangent_projection(
        int node_id, const std::array<bool, 3>& fixed_directions = {}) const;

    /*!
     * @brief Restrict the Newton step to the tangent spaces of the boundary
     *
     * The system matrix K and the residual r are replaced by P^T K P + (I - P) and P^T r with the
     * projections P of the nodes, such that the increments of boundary nodes are tangential to the
     * boundary. This enforces the constraints exactly, the projection onto the boundary after the
     * update accounts for its curvature. The dofs in @p dirichlet_dofs are fixed directions.
     */
    void restrict_to_tangent_space(Core::LinAlg::SparseMatrix& matrix,
        Core::LinAlg::Vector<double>& residual, const Core::LinAlg::Map& dirichlet_dofs) const;

    /*!
     * @brief Project the moved boundary nodes onto the boundary of the reference configuration
     *
     * The normals of the next Newton step are updated to the faces at the projected positions. The
     * dofs in @p dirichlet_dofs keep their values.
     *
     * The closest point on each patch of a node is searched among the faces of the reference
     * configuration near the moved node (global geometric search over their bounding boxes), which
     * are reachable from the faces of the patch at the node without crossing a feature edge.
     * Collective.
     */
    void project_onto_boundary(
        Core::LinAlg::Vector<double>& displacement, const Core::LinAlg::Map& dirichlet_dofs);

    /*!
     * @brief Remove the components of a residual that are balanced by the boundary
     *
     * At the constrained minimum, the gradient of the potential at a boundary node is balanced by
     * the boundary in the directions of the normals of the boundary faces at the node: one normal
     * inside a face, several at the edges and vertices of the piecewise flat boundary. Only the
     * remaining components measure convergence.
     */
    void remove_constrained_components(Core::LinAlg::Vector<double>& residual) const;

   private:
    struct Face
    {
      //! ids of the corner nodes, ordered around the face
      std::vector<int> node_ids;
      //! positions of the corner nodes in the reference configuration of the mesh movement
      std::vector<Core::LinAlg::Matrix<3, 1>> corners;
      //! unit normal in the reference configuration of the mesh movement
      Core::LinAlg::Matrix<3, 1> normal;
      //! unit normal in the initial configuration, which defines the feature edges
      Core::LinAlg::Matrix<3, 1> initial_normal;
    };

    struct NodeConstraint
    {
      //! the boundary faces adjacent to the node, grouped by patch
      std::vector<std::vector<int>> patch_faces;
      //! normal of the tangent plane for each patch
      std::vector<Core::LinAlg::Matrix<3, 1>> normals;
      //! normals of all boundary faces at the current position, whose span balances the gradient
      //! of the potential at the constrained minimum (several at edges and vertices of the faces)
      std::vector<Core::LinAlg::Matrix<3, 1>> active_normals;

      [[nodiscard]] bool is_corner() const { return patch_faces.size() > 2; }
    };

    //! position of a node in the reference configuration
    [[nodiscard]] Core::LinAlg::Matrix<3, 1> position(int node_id) const;

    //! positions of nodes in the reference configuration
    [[nodiscard]] std::vector<Core::LinAlg::Matrix<3, 1>> corner_positions(
        const std::vector<int>& node_ids) const;

    //! dof row lids of a node
    [[nodiscard]] std::vector<int> dof_lids(int node_id) const;

    //! whether two known faces sharing an edge belong to the same patch
    [[nodiscard]] bool same_patch(int a, int b) const;

    //! add a face to the known faces
    void add_face(int id, Face face);

    //! fetch unknown faces (by id) from their processors (by id) and add them to the known faces
    void fetch_faces(const std::map<int, int>& face_processors);

    //! the (sorted) candidate faces reachable from the faces of a patch at a node without crossing
    //! a feature edge
    [[nodiscard]] std::vector<int> patch_candidates(
        const std::vector<int>& patch_faces, const std::vector<int>& candidates) const;

    //! closest point on a set of faces and the face containing it
    struct ClosestPoint
    {
      Core::LinAlg::Matrix<3, 1> point;
      //! the face containing the closest point
      int face;
      //! all faces containing the closest point (several at their edges and vertices)
      std::vector<int> touching_faces;
    };

    [[nodiscard]] ClosestPoint closest_point(
        const Core::LinAlg::Matrix<3, 1>& point, const std::vector<int>& faces) const;

    const Core::FE::Discretization& discretization_;

    double feature_angle_cosine_;

    //! boundary faces with a row node (by global id), with their initial normals
    std::map<int, Face> row_node_faces_;

    //! ids of the boundary faces owned by this processor
    std::vector<int> owned_faces_;

    //! boundary faces adjacent to each row node (by node id)
    std::map<int, std::vector<int>> node_faces_;

    //! known faces in the reference configuration of the current mesh movement (by global id): the
    //! faces with a row node and the faces fetched from other processors
    std::map<int, Face> faces_;

    //! bounding boxes of the owned faces in the reference configuration of the current mesh movement
    std::vector<std::pair<int, Core::GeometricSearch::BoundingVolume>> face_boxes_;

    //! known faces adjacent to each edge (by the sorted ids of its nodes)
    std::map<std::pair<int, int>, std::vector<int>> edge_faces_;

    //! constraints of the boundary row nodes (by node id)
    std::map<int, NodeConstraint> constraints_;
  };

}  // namespace Solid::MeshRegularization

FOUR_C_NAMESPACE_CLOSE

#endif
