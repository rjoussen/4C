// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_structure_new_mesh_regularization.hpp"

#include "4C_adapter_ale.hpp"
#include "4C_ale_input.hpp"
#include "4C_ale_mesh_quality.hpp"
#include "4C_ale_utils_clonestrategy.hpp"
#include "4C_ale_utils_mapextractor.hpp"
#include "4C_fem_discretization.hpp"
#include "4C_fem_general_node.hpp"
#include "4C_fem_general_utils_createdis.hpp"
#include "4C_global_data.hpp"
#include "4C_io.hpp"
#include "4C_io_control.hpp"
#include "4C_io_pstream.hpp"
#include "4C_linalg_sparsematrix.hpp"
#include "4C_linalg_utils_sparse_algebra_manipulation.hpp"
#include "4C_linear_solver_method_linalg.hpp"
#include "4C_utils_parameter_list.hpp"

#include <Teuchos_StandardParameterEntryValidators.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <string>

FOUR_C_NAMESPACE_OPEN

namespace
{
  //! the ALE elements of the mesh movement have no material, the mesh quality potential does not
  //! need one
  constexpr int no_material = -1;

  /*!
   * @brief Clone strategy for the ALE field of the mesh movement, which creates ALE elements
   * without material
   */
  class MeshMovementCloneStrategy : public ALE::Utils::AleCloneStrategy
  {
   protected:
    void check_material_type(const int matid) {}

    void set_element_data(std::shared_ptr<Core::Elements::Element> newele,
        Core::Elements::Element* oldele, const int matid, const bool nurbsdis)
    {
    }
  };

  //! copy nodal values between discretizations with matching nodes
  void copy_nodal_values(const Core::FE::Discretization& source_discretization,
      const Core::LinAlg::Vector<double>& source,
      const Core::FE::Discretization& target_discretization, Core::LinAlg::Vector<double>& target)
  {
    for (int lid = 0; lid < target_discretization.num_my_row_nodes(); ++lid)
    {
      const Core::Nodes::Node* target_node = target_discretization.l_row_node(lid);
      const Core::Nodes::Node* source_node = source_discretization.g_node(target_node->id());
      const std::vector<int> source_dofs = source_discretization.dof(0, source_node);
      const std::vector<int> target_dofs = target_discretization.dof(0, target_node);
      for (int d = 0; d < 3; ++d)
      {
        target.get_values()[target.get_map().lid(target_dofs[d])] =
            source.local_values_as_span()[source.get_map().lid(source_dofs[d])];
      }
    }
  }
}  // namespace

Solid::MeshRegularization::MeshMovement::MeshMovement(Global::Problem& problem,
    Core::FE::Discretization& structure,
    const Teuchos::ParameterList& structural_dynamics_parameters,
    const Solid::TimeInt::MeshRegularizationConfig& config)
    : structure_(structure),
      quality_parameters_{
          .size_weight = problem.ale_dynamic_params().get<double>("MESH_QUALITY_SIZE_WEIGHT"),
          .motion_weight = problem.ale_dynamic_params().get<double>("MESH_QUALITY_MOTION_WEIGHT"),
          .size_growth_only =
              problem.ale_dynamic_params().get<bool>("MESH_QUALITY_SIZE_GROWTH_ONLY"),
          .shape_exponent = problem.ale_dynamic_params().get<double>("MESH_QUALITY_SHAPE_EXPONENT")}
{
  FOUR_C_ASSERT_ALWAYS(Teuchos::getIntegralValue<ALE::AleDynamic>(
                           problem.ale_dynamic_params(), "ALE_TYPE") == ALE::mesh_quality,
      "The mesh regularization requires ALE_TYPE mesh_quality in the section ALE DYNAMIC.");
  FOUR_C_ASSERT_ALWAYS(quality_parameters_.shape_exponent >= 1.0,
      "MESH_QUALITY_SHAPE_EXPONENT in the section ALE DYNAMIC has to be at least 1.");
  FOUR_C_ASSERT_ALWAYS(
      quality_parameters_.size_weight >= 0.0 && quality_parameters_.motion_weight >= 0.0 &&
          quality_parameters_.size_weight + quality_parameters_.motion_weight > 0.0,
      "The mesh regularization requires a positive MESH_QUALITY_SIZE_WEIGHT or "
      "MESH_QUALITY_MOTION_WEIGHT in the section ALE DYNAMIC: the shape term alone does not "
      "determine the size of the elements, so that the mesh movement is not unique.");
  FOUR_C_ASSERT_ALWAYS(
      !quality_parameters_.size_growth_only || quality_parameters_.motion_weight > 0.0,
      "The mesh regularization with MESH_QUALITY_SIZE_GROWTH_ONLY requires a positive "
      "MESH_QUALITY_MOTION_WEIGHT in the section ALE DYNAMIC: the size term does not determine the "
      "size of shrinking elements, so that the mesh movement is not unique.");

  ale_discretization_ =
      std::make_shared<Core::FE::Discretization>("ale", structure.get_comm(), structure.n_dim());
  ale_discretization_->fill_complete();
  Core::FE::DiscretizationCreator<MeshMovementCloneStrategy> creator;
  creator.create_matching_discretization(structure, *ale_discretization_, no_material);
  ale_discretization_->fill_complete();
  ale_discretization_->set_writer(std::make_shared<Core::IO::DiscretizationWriter>(
      *ale_discretization_, *problem.output_control_file(), problem.spatial_approximation_type()));

  Adapter::AleBaseAlgorithm ale_algorithm(
      problem, structural_dynamics_parameters, ale_discretization_);
  ale_ = ale_algorithm.ale_field();

  sliding_ = std::make_unique<BoundarySliding>(problem, ale_discretization_, config.feature_angle);

  Adapter::AleNonlinearSolverHooks hooks;
  // the target of the mesh quality is the initial mesh
  hooks.set_states = [this](Core::FE::Discretization& discretization)
  { discretization.set_state("mesh quality target", *target_); };
  // the boundary nodes slide on the boundary (mortar constraints with condensed multipliers), fixed
  // dofs do not move
  hooks.modify_system = [this](Core::LinAlg::SparseOperator& matrix,
                            Core::LinAlg::Vector<double>& residual,
                            const Core::LinAlg::Vector<double>& dispnp)
  { sliding_->condense(dynamic_cast<Core::LinAlg::SparseMatrix&>(matrix), residual, dispnp); };
  // steps into inverted elements (infinite potential) are not admissible
  hooks.accept_update = [this](Core::LinAlg::Vector<double>& dispnp)
  { return std::isfinite(total_energy()); };
  // a mesh movement that does not converge may have moved the mesh arbitrarily far, it is discarded
  hooks.on_failure = [this]() { failed_ = true; };
  ale_->set_nonlinear_solver_hooks(std::move(hooks));
}

void Solid::MeshRegularization::MeshMovement::set_ale_reference_to_current_configuration(
    const Core::LinAlg::Vector<double>& displacement)
{
  Core::LinAlg::Vector<double> displacement_col(*structure_.dof_col_map());
  Core::LinAlg::export_to(displacement, displacement_col);

  for (int lid = 0; lid < ale_discretization_->num_my_col_nodes(); ++lid)
  {
    Core::Nodes::Node* ale_node = ale_discretization_->l_col_node(lid);
    const Core::Nodes::Node* structure_node = structure_.g_node(ale_node->id());
    const std::vector<int> dofs = structure_.dof(0, structure_node);

    std::vector<double> position(3);
    for (int d = 0; d < 3; ++d)
    {
      position[d] =
          structure_node->x()[d] +
          displacement_col.local_values_as_span()[displacement_col.get_map().lid(dofs[d])];
    }
    ale_node->set_pos(position);
  }
}

double Solid::MeshRegularization::MeshMovement::total_energy() const
{
  Core::LinAlg::Vector<double> displacement_col(*ale_discretization_->dof_col_map());
  Core::LinAlg::export_to(*ale_->dispnp(), displacement_col);
  Core::LinAlg::Vector<double> target_col(*ale_discretization_->dof_col_map());
  Core::LinAlg::export_to(*target_, target_col);

  double energy = 0.0;
  for (int lid = 0; lid < ale_discretization_->num_my_row_elements(); ++lid)
  {
    const Core::Elements::Element* element = ale_discretization_->l_row_element(lid);
    std::vector<double> reference_coordinates;
    std::vector<double> displacements;
    std::vector<double> target_coordinates;
    for (int i = 0; i < element->num_node(); ++i)
    {
      const Core::Nodes::Node* node = element->nodes()[i];
      const std::vector<int> dofs = ale_discretization_->dof(0, node);
      for (int d = 0; d < 3; ++d)
      {
        const int lid = displacement_col.get_map().lid(dofs[d]);
        reference_coordinates.push_back(node->x()[d]);
        displacements.push_back(displacement_col.local_values_as_span()[lid]);
        target_coordinates.push_back(node->x()[d] + target_col.local_values_as_span()[lid]);
      }
    }
    energy += ALE::evaluate_mesh_quality_potential(element->shape(), reference_coordinates,
        displacements, quality_parameters_, nullptr, nullptr, target_coordinates);
  }
  energy = Core::Communication::sum_all(energy, ale_discretization_->get_comm());
  return energy;
}

std::shared_ptr<Core::LinAlg::Map> Solid::MeshRegularization::MeshMovement::fixed_ale_dofs(
    const Core::LinAlg::Map& solid_dirichlet_dofs) const
{
  const Core::LinAlg::Map& ale_dirichlet_dofs =
      *ale_->get_dbc_map_extractor(ALE::Utils::MapExtractor::dbc_set_std)->cond_map();

  std::vector<int> fixed;
  for (int lid = 0; lid < ale_discretization_->num_my_row_nodes(); ++lid)
  {
    const Core::Nodes::Node* ale_node = ale_discretization_->l_row_node(lid);
    const std::vector<int> ale_dofs = ale_discretization_->dof(0, ale_node);
    const std::vector<int> solid_dofs = structure_.dof(0, structure_.g_node(ale_node->id()));
    for (int d = 0; d < 3; ++d)
    {
      if (ale_dirichlet_dofs.lid(ale_dofs[d]) >= 0 || solid_dirichlet_dofs.lid(solid_dofs[d]) >= 0)
        fixed.push_back(ale_dofs[d]);
    }
  }
  return std::make_shared<Core::LinAlg::Map>(
      -1, static_cast<int>(fixed.size()), fixed.data(), 0, ale_discretization_->get_comm());
}

void Solid::MeshRegularization::MeshMovement::solve()
{
  FOUR_C_ASSERT_ALWAYS(std::isfinite(total_energy()),
      "The mesh movement starts from a mesh with an element that is inverted at a corner. The "
      "mesh has to be regularized earlier (START_STEP) or more often (INTERVAL).");

  ale_->time_step();
}

std::optional<Core::LinAlg::Vector<double>> Solid::MeshRegularization::MeshMovement::compute(
    const Core::LinAlg::Vector<double>& displacement, const Core::LinAlg::Map& dirichlet_dofs)
{
  set_ale_reference_to_current_configuration(displacement);

  ale_->write_access_dispnp()->put_scalar(0.0);
  ale_->prepare_time_step();

  // the initial mesh is the current configuration minus the displacement
  target_ =
      std::make_shared<Core::LinAlg::Vector<double>>(*ale_discretization_->dof_row_map(), true);
  copy_nodal_values(structure_, displacement, *ale_discretization_, *target_);
  target_->scale(-1.0);

  sliding_->setup(*fixed_ale_dofs(dirichlet_dofs));
  failed_ = false;
  solve();
  if (failed_) return std::nullopt;
  Core::LinAlg::Vector<double> movement(*structure_.dof_row_map(), true);
  copy_nodal_values(*ale_discretization_, *ale_->dispnp(), structure_, movement);

  return movement;
}

FOUR_C_NAMESPACE_CLOSE
