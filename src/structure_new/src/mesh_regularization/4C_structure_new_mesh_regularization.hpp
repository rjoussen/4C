// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_STRUCTURE_NEW_MESH_REGULARIZATION_HPP
#define FOUR_C_STRUCTURE_NEW_MESH_REGULARIZATION_HPP

#include "4C_config.hpp"

#include "4C_ale_mesh_quality.hpp"
#include "4C_linalg_vector.hpp"
#include "4C_structure_new_mesh_regularization_sliding.hpp"
#include "4C_structure_new_timint_basedatasdyn.hpp"

#include <Teuchos_ParameterList.hpp>

#include <memory>
#include <optional>

FOUR_C_NAMESPACE_OPEN

namespace Adapter
{
  class Ale;
}

namespace Core::FE
{
  class Discretization;
}

namespace Global
{
  class Problem;
}

namespace Solid::MeshRegularization
{
  /*!
   * @brief Computes a regularizing movement of the nodes of a solid discretization
   *
   * The movement is computed with an ALE field, which is cloned from the solid discretization.
   * Before each regularization, the reference configuration of the ALE field is set to the
   * current configuration of the solid, such that the ALE displacement is the movement of the
   * nodes in the current configuration. The mesh quality potential of the ALE field, whose target
   * is the initial mesh with its element shapes and sizes, is minimized
   * with the Newton scheme of the ALE field, whose steps are halved as long as they would invert an
   * element or do not decrease the residual sufficiently. Boundary nodes slide along the boundary.
   * Dofs with ALE Dirichlet conditions and dofs with Dirichlet conditions of the solid are fixed.
   */
  class MeshMovement
  {
   public:
    /*!
     * @param problem Global problem providing the ALE parameters and linear solvers
     * @param structure Solid discretization
     * @param structural_dynamics_parameters Structural dynamics parameters
     * @param config Parameters of the mesh regularization
     */
    MeshMovement(Global::Problem& problem, Core::FE::Discretization& structure,
        const Teuchos::ParameterList& structural_dynamics_parameters,
        const Solid::TimeInt::MeshRegularizationConfig& config);

    /*!
     * @brief Compute the movement of the nodes in the current configuration
     *
     * @param displacement Displacement of the solid (dof row map of the solid)
     * @param dirichlet_dofs Dofs of the solid with Dirichlet conditions, which are not moved
     * @return Movement of the nodes (dof row map of the solid), or nothing if the minimization did
     * not converge
     */
    std::optional<Core::LinAlg::Vector<double>> compute(
        const Core::LinAlg::Vector<double>& displacement, const Core::LinAlg::Map& dirichlet_dofs);

   private:
    //! move the reference configuration of the ALE field to the current solid configuration
    void set_ale_reference_to_current_configuration(
        const Core::LinAlg::Vector<double>& displacement);

    //! solve the ALE problem by minimizing the mesh quality potential
    void solve();

    //! total mesh quality potential of the current ALE state
    [[nodiscard]] double total_energy() const;

    //! ALE dofs with ALE Dirichlet conditions or Dirichlet conditions of the solid
    [[nodiscard]] std::shared_ptr<Core::LinAlg::Map> fixed_ale_dofs(
        const Core::LinAlg::Map& solid_dirichlet_dofs) const;

    Core::FE::Discretization& structure_;

    std::shared_ptr<Core::FE::Discretization> ale_discretization_;

    std::shared_ptr<Adapter::Ale> ale_;

    ALE::MeshQualityParameters quality_parameters_;

    //! sliding of boundary nodes
    std::unique_ptr<BoundarySliding> sliding_;

    //! target configuration of the mesh quality (the initial mesh) as displacement from the
    //! reference configuration of the ALE field (dof row map of the ALE field)
    std::shared_ptr<Core::LinAlg::Vector<double>> target_;

    //! whether the minimization of the current mesh movement did not converge
    bool failed_ = false;
  };

}  // namespace Solid::MeshRegularization

FOUR_C_NAMESPACE_CLOSE

#endif
