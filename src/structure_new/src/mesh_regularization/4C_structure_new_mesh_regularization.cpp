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
#include "4C_utils_function.hpp"
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

  //! relative round-off of the potential summed over all elements
  constexpr double energy_round_off = 1.0e-12;

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

  /*!
   * @brief Target volumes h^3 of the elements for a size function h(x, t) at the current position
   * of a node
   *
   * The size h is the edge length of a cube with the target volume. Its second derivatives are
   * central differences of its first derivatives.
   */
  ALE::TargetVolumeField target_volume_field(
      const Core::Utils::FunctionOfSpaceTime& function, const double time)
  {
    return [&function, time](const std::array<double, 3>& position)
    {
      const double h = function.evaluate(position, time, 0);
      FOUR_C_ASSERT_ALWAYS(
          h > 0.0, "The size function MESH_QUALITY_SIZE_FUNCTION has to be positive.");
      const std::vector<double> gradient = function.evaluate_spatial_derivative(position, time, 0);

      std::array<std::array<double, 3>, 3> hessian{};
      for (int j = 0; j < 3; ++j)
      {
        const double step = 1.0e-6 * std::max(1.0, std::abs(position[j]));
        std::array<double, 3> forward = position;
        std::array<double, 3> backward = position;
        forward[j] += step;
        backward[j] -= step;
        const std::vector<double> gradient_forward =
            function.evaluate_spatial_derivative(forward, time, 0);
        const std::vector<double> gradient_backward =
            function.evaluate_spatial_derivative(backward, time, 0);
        for (int i = 0; i < 3; ++i)
          hessian[i][j] = (gradient_forward[i] - gradient_backward[i]) / (2.0 * step);
      }

      // V = h^3: dV = 3 h^2 dh, d2V = 6 h dh dh + 3 h^2 d2h
      ALE::TargetVolume volume{.value = h * h * h};
      for (int i = 0; i < 3; ++i)
      {
        volume.gradient[i] = 3.0 * h * h * gradient[i];
        for (int j = 0; j < 3; ++j)
          volume.hessian[i][j] = 6.0 * h * gradient[i] * gradient[j] + 3.0 * h * h * hessian[i][j];
      }
      return volume;
    };
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
          .target_volume = {}}
{
  FOUR_C_ASSERT_ALWAYS(Teuchos::getIntegralValue<ALE::AleDynamic>(
                           problem.ale_dynamic_params(), "ALE_TYPE") == ALE::mesh_quality,
      "The mesh regularization requires ALE_TYPE mesh_quality in the section ALE DYNAMIC.");
  FOUR_C_ASSERT_ALWAYS(quality_parameters_.size_weight >= 0.0 &&
                           quality_parameters_.motion_weight >= 0.0 &&
                           quality_parameters_.size_weight + quality_parameters_.motion_weight > 0.0,
      "The mesh regularization requires a positive MESH_QUALITY_SIZE_WEIGHT or "
      "MESH_QUALITY_MOTION_WEIGHT in the section ALE DYNAMIC: the shape term alone does not "
      "determine the size of the elements, so that the mesh movement is not unique.");

  ale_discretization_ =
      std::make_shared<Core::FE::Discretization>("ale", structure.get_comm(), structure.n_dim());
  ale_discretization_->fill_complete();
  Core::FE::DiscretizationCreator<MeshMovementCloneStrategy> creator;
  creator.create_matching_discretization(structure, *ale_discretization_, no_material);
  ale_discretization_->fill_complete();
  ale_discretization_->set_writer(std::make_shared<Core::IO::DiscretizationWriter>(
      *ale_discretization_, *problem.output_control_file(), problem.spatial_approximation_type()));

  // the target sizes of the elements change by the ratio of a function of the position
  if (Teuchos::getIntegralValue<ALE::MeshQualitySizeTarget>(problem.ale_dynamic_params(),
          "MESH_QUALITY_SIZE_TARGET") == ALE::MeshQualitySizeTarget::function)
  {
    const int function_id = problem.ale_dynamic_params().get<int>("MESH_QUALITY_SIZE_FUNCTION");
    FOUR_C_ASSERT_ALWAYS(function_id > 0,
        "MESH_QUALITY_SIZE_TARGET function requires MESH_QUALITY_SIZE_FUNCTION in the section ALE "
        "DYNAMIC.");
    size_function_ = &problem.function_by_id<Core::Utils::FunctionOfSpaceTime>(function_id);
  }

  Adapter::AleBaseAlgorithm ale_algorithm(
      problem, structural_dynamics_parameters, ale_discretization_);
  ale_ = ale_algorithm.ale_field();

  sliding_ = std::make_unique<BoundarySliding>(*ale_discretization_, config.feature_angle);

  Adapter::AleNonlinearSolverHooks hooks;
  // the target of the mesh quality is the initial mesh, with sizes scaled to the sizes of the
  // initial mesh at the positions of the elements if given
  hooks.set_states = [this](Core::FE::Discretization& discretization)
  { discretization.set_state("mesh quality target", *target_); };
  hooks.set_element_parameters = [this](Teuchos::ParameterList& element_parameters)
  {
    if (quality_parameters_.target_volume)
      element_parameters.set("mesh quality target volume", quality_parameters_.target_volume);
  };
  // the Newton steps of the boundary nodes are tangential to the boundary, fixed dofs do not move
  hooks.modify_system = [this](Core::LinAlg::SparseOperator& matrix,
                            Core::LinAlg::Vector<double>& residual,
                            const Core::LinAlg::Vector<double>&)
  {
    sliding_->restrict_to_tangent_space(
        dynamic_cast<Core::LinAlg::SparseMatrix&>(matrix), residual, *fixed_dofs_);
  };
  // boundary nodes are projected onto the boundary after each step, steps that do not decrease the
  // potential are halved (the potential is infinite for inverted elements)
  hooks.accept_update = [this](Core::LinAlg::Vector<double>& dispnp)
  {
    sliding_->project_onto_boundary(dispnp, *fixed_dofs_);
    // a step may increase the potential within its round-off, which allows the last steps of a
    // converging iteration
    const double energy = total_energy();
    if (!(energy <= accepted_energy_ + energy_round_off * std::abs(accepted_energy_))) return false;
    accepted_energy_ = energy;
    return true;
  };
  hooks.filter_residual = [this](Core::LinAlg::Vector<double>& residual)
  { sliding_->remove_constrained_components(residual); };
  // a mesh movement that does not converge may have moved the mesh arbitrarily far
  hooks.require_convergence = true;
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
  accepted_energy_ = total_energy();
  FOUR_C_ASSERT_ALWAYS(std::isfinite(accepted_energy_),
      "The mesh movement starts from a mesh with an element that is inverted at a corner. The "
      "mesh has to be regularized earlier (START_STEP) or more often (INTERVAL).");

  ale_->time_step();
}

Core::LinAlg::Vector<double> Solid::MeshRegularization::MeshMovement::compute(
    const Core::LinAlg::Vector<double>& displacement, const Core::LinAlg::Map& dirichlet_dofs,
    const double time)
{
  set_ale_reference_to_current_configuration(displacement);

  ale_->write_access_dispnp()->put_scalar(0.0);
  ale_->prepare_time_step();
  fixed_dofs_ = fixed_ale_dofs(dirichlet_dofs);

  // the initial mesh is the current configuration minus the displacement
  target_ =
      std::make_shared<Core::LinAlg::Vector<double>>(*ale_discretization_->dof_row_map(), true);
  copy_nodal_values(structure_, displacement, *ale_discretization_, *target_);
  target_->scale(-1.0);

  // the target volume of the elements at the current position of each node
  if (size_function_) quality_parameters_.target_volume = target_volume_field(*size_function_, time);
  sliding_->setup();
  solve();
  Core::LinAlg::Vector<double> movement(*structure_.dof_row_map(), true);
  copy_nodal_values(*ale_discretization_, *ale_->dispnp(), structure_, movement);

  return movement;
}

FOUR_C_NAMESPACE_CLOSE
