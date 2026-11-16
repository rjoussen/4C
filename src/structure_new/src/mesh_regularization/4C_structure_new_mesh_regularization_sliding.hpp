// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_STRUCTURE_NEW_MESH_REGULARIZATION_SLIDING_HPP
#define FOUR_C_STRUCTURE_NEW_MESH_REGULARIZATION_SLIDING_HPP

#include "4C_config.hpp"

#include "4C_fem_general_cell_type.hpp"
#include "4C_linalg_vector.hpp"

#include <array>
#include <map>
#include <memory>
#include <set>
#include <vector>

FOUR_C_NAMESPACE_OPEN

namespace Adapter
{
  class CouplingNonLinMortar;
}

namespace Core::FE
{
  class Discretization;
}

namespace Core::LinAlg
{
  class Map;
  class SparseMatrix;
}  // namespace Core::LinAlg

namespace Global
{
  class Problem;
}

namespace Solid::MeshRegularization
{
  /*!
   * @brief Sliding of boundary nodes along the boundary during the mesh movement
   *
   * The boundary faces of the discretization are detected automatically and grouped into smooth
   * patches in the initial configuration: adjacent faces whose normals enclose an angle below the
   * feature angle belong to the same patch.
   *
   * Before each mesh movement, the boundary faces are copied into a fixed discretization, which
   * keeps the shape of the domain in the reference configuration of the mesh movement. The nodes
   * of each patch slide on the copy of the patch: a mortar coupling with dual Lagrange multipliers
   * keeps the weighted normal gaps, and the multipliers are condensed. A patch whose normal
   * directions are fixed at all of its nodes (e.g., a symmetry plane) needs no coupling. A node on
   * several coupled patches is a corner and is fixed. Fixed directions of the other nodes are
   * removed from their normals (mortar symmetry conditions).
   */
  class BoundarySliding
  {
   public:
    /*!
     * @brief Detect the boundary faces and their patches in the current configuration of the
     * discretization, which is the initial one (collective)
     *
     * @param problem Problem with the parameters of the mortar coupling
     * @param discretization Discretization of the mesh movement, whose column elements contain all
     * elements adjacent to the row nodes
     * @param feature_angle Angle in degrees above which adjacent boundary faces belong to different
     * patches
     */
    BoundarySliding(Global::Problem& problem,
        std::shared_ptr<Core::FE::Discretization> discretization, double feature_angle);

    ~BoundarySliding();

    /*!
     * @brief Set up the mortar couplings in the current reference configuration of the
     * discretization (before each mesh movement, collective)
     *
     * @param dirichlet_dofs Fixed dofs of the mesh movement
     */
    void setup(const Core::LinAlg::Map& dirichlet_dofs);

    /*!
     * @brief Add the sliding constraints to the linearized mesh movement at a displacement
     *
     * The multipliers z = -D^{-1} (r_s + K_s du) of the slave dofs are condensed. The rows of each
     * slave node are replaced by its tangential equilibrium -T (r_s + K_s du) + H du = 0 and its
     * linearized weighted normal gap N du + g = 0, where H is the linearization of the tangents T
     * contracted with the multipliers D z = -r_s of the displacement, which balance the forces of
     * the potential at the solution. The residual and the system depend on the displacement only.
     * The fixed dofs (Dirichlet dofs and corners) get unit rows and zero residuals. If a slave face
     * does not overlap the fixed boundary, the displacement is not admissible and the residual is
     * not a number. Collective.
     */
    void condense(Core::LinAlg::SparseMatrix& matrix, Core::LinAlg::Vector<double>& residual,
        const Core::LinAlg::Vector<double>& displacement);

    //! The patches of the boundary faces at a row node (sorted global ids)
    [[nodiscard]] std::vector<int> patches(int node_id) const;

   private:
    struct Face
    {
      //! ids of all nodes, in the order of the face
      std::vector<int> node_ids;
      //! shape of the face
      Core::FE::CellType shape;
    };

    //! unit normal of a face in the current reference configuration
    [[nodiscard]] std::array<double, 3> normal(const Face& face) const;

    //! group the faces into patches: adjacent faces whose normals enclose an angle below the
    //! feature angle
    void group_into_patches();

    Global::Problem& problem_;

    std::shared_ptr<Core::FE::Discretization> discretization_;

    //! boundary faces with a row node (by global id)
    std::map<int, Face> faces_;

    //! cosine of the feature angle
    double feature_angle_cosine_;

    //! patch of each face in faces_ (by face id)
    std::map<int, int> patch_of_face_;

    //! boundary faces adjacent to each row node (by node id)
    std::map<int, std::vector<int>> node_faces_;

    //! weighted gaps (times the normals) at the start of the current mesh movement
    std::unique_ptr<Core::LinAlg::Vector<double>> initial_gap_;

    //! copies of the faces of the coupled patches in the reference configuration of the current
    //! mesh movement, with offset ids and their own dofs
    std::shared_ptr<Core::FE::Discretization> fixed_boundary_;

    //! inflation factor of the search boxes of the mortar couplings for no movement
    double search_inflation_factor_ = 0.0;

    //! mortar coupling of each coupled patch
    std::vector<std::unique_ptr<Adapter::CouplingNonLinMortar>> couplings_;

    //! dofs of the mesh movement and of the fixed boundary, on which the mortar states are given
    std::shared_ptr<Core::LinAlg::Map> state_map_;

    //! fixed dofs of the mesh movement: Dirichlet dofs and corners
    std::shared_ptr<Core::LinAlg::Map> fixed_dofs_;
  };

}  // namespace Solid::MeshRegularization

FOUR_C_NAMESPACE_CLOSE

#endif
