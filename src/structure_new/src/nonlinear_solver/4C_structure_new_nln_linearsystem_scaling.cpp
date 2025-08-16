// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_structure_new_nln_linearsystem_scaling.hpp"

#include "4C_fem_condition_utils.hpp"
#include "4C_fem_discretization.hpp"
#include "4C_fem_general_cell_type.hpp"
#include "4C_fem_general_cell_type_traits.hpp"
#include "4C_fem_general_extract_values.hpp"
#include "4C_inpar_structure.hpp"
#include "4C_linalg_sparsematrix.hpp"
#include "4C_linalg_symmetric_tensor.hpp"
#include "4C_linalg_tensor.hpp"
#include "4C_linalg_utils_sparse_algebra_create.hpp"
#include "4C_linalg_utils_sparse_algebra_manipulation.hpp"
#include "4C_linalg_utils_sparse_algebra_math.hpp"
#include "4C_linalg_vector.hpp"
#include "4C_solid_3D_ele_calc_lib.hpp"
#include "4C_solver_nonlin_nox_linearproblem.hpp"
#include "4C_structure_new_timint_basedataglobalstate.hpp"
#include "4C_structure_new_timint_basedatasdyn.hpp"
#include "4C_utils_exceptions.hpp"

FOUR_C_NAMESPACE_OPEN

namespace
{

  template <Core::FE::CellType celltype>
  double compute_aspect_ratio(
      const Core::Elements::Element& ele, const std::vector<double>& displacements)
  {
    static_assert(celltype == Core::FE::CellType::hex8,
        "Only Hex8 elements are supported for STC matrix calculation");

    Discret::Elements::ElementNodes<celltype> element_nodes =
        Discret::Elements::evaluate_element_nodes<celltype>(ele, displacements);

    Discret::Elements::JacobianMapping<celltype> jacobian_mapping =
        Discret::Elements::evaluate_jacobian_mapping_centroid(element_nodes);


    Core::LinAlg::SymmetricTensor<double, Core::FE::dim<celltype>, Core::FE::dim<celltype>>
        jac0stretch = Core::LinAlg::assume_symmetry(
            Core::LinAlg::transpose(jacobian_mapping.inverse_jacobian_) *
            jacobian_mapping.inverse_jacobian_);
    const double r_stretch = sqrt(jac0stretch(0, 0));
    const double s_stretch = sqrt(jac0stretch(1, 1));
    const double t_stretch = sqrt(jac0stretch(2, 2));

    // return an averaged aspect ratio
    if (r_stretch >= s_stretch and r_stretch >= t_stretch)
    {
      return 0.5 * (r_stretch / s_stretch + r_stretch / t_stretch);
    }
    else if (s_stretch > r_stretch and s_stretch >= t_stretch)
    {
      return 0.5 * (s_stretch / r_stretch + s_stretch / t_stretch);
    }
    else if (t_stretch > r_stretch and t_stretch > s_stretch)
    {
      return 0.5 * (t_stretch / s_stretch + t_stretch / r_stretch);
    }

    return 0.0;
  }

  template <Core::FE::CellType celltype>
  Core::LinAlg::SerialDenseMatrix evaluate_stc_matrix(const Core::Elements::Element& ele,
      const std::vector<double> displacements, const Inpar::Solid::StcScale stc_scaling,
      int stc_layer,
      const std::multimap<int, const Core::Conditions::Condition*>& stc_layer_conditions)
  {
    const Core::Nodes::Node* const* nodes = ele.nodes();

    const double aspect_ratio = compute_aspect_ratio<celltype>(ele, displacements);
    const double stc_factor =
        stc_scaling == Inpar::Solid::stc_currsym ? aspect_ratio : aspect_ratio * aspect_ratio;

    const double factor1 = (stc_factor + 1.0) / (2.0 * stc_factor);
    const double factor2 = (stc_factor - 1.0) / (2.0 * stc_factor);
    const double factor3 = (1.0 / stc_factor);
    const double factor4 = (1.0 - 1.0 / stc_factor);

    Core::LinAlg::SerialDenseMatrix stc_matrix(
        Core::FE::num_nodes(celltype) * Core::FE::dim<celltype>,
        Core::FE::num_nodes(celltype) * Core::FE::dim<celltype>);

    std::vector<const Core::Conditions::Condition*> cond0;
    {
      auto range = stc_layer_conditions.equal_range(nodes[0]->id());
      std::ranges::copy(std::ranges::subrange(range.first, range.second) | std::views::values,
          std::back_inserter(cond0));
    }
    int condnum0 = 1000;    // minimum STCid of layer with nodes 0..3
    bool current0 = false;  // layer with nodes 0..4 to be scaled

    std::vector<const Core::Conditions::Condition*> cond1;
    {
      auto range = stc_layer_conditions.equal_range(nodes[Core::FE::num_nodes(celltype) / 2]->id());
      std::ranges::copy(std::ranges::subrange(range.first, range.second) | std::views::values,
          std::back_inserter(cond1));
    }
    int condnum1 = 1000;    // minimum STCid of layer with nodes 4..7
    bool current1 = false;  // minimum STCid of layer with nodes 4..7

    for (auto& conu : cond0)
    {
      int tmp = conu->parameters().get<int>("ConditionID");
      if (tmp < condnum0) condnum0 = tmp;
    }
    if (condnum0 == stc_layer) current0 = true;


    for (auto& conu : cond1)
    {
      int tmp = conu->parameters().get<int>("ConditionID");
      if (tmp < condnum1) condnum1 = tmp;
    }
    if (condnum1 == stc_layer) current1 = true;

    std::array<int, Core::FE::num_nodes(celltype)> number_adj_elements{};
    std::transform(nodes, nodes + Core::FE::num_nodes(celltype), number_adj_elements.begin(),
        [](const Core::Nodes::Node* node) { return node->num_element(); });

    constexpr int numdof = Core::FE::num_nodes(celltype) * Core::FE::dim<celltype>;

    // both surfaces are to be scaled
    if (current0 and current1)
    {
      // only valid for first round
      if (condnum0 != 1)
        FOUR_C_THROW("STC error: non-initial layer is not connected to a smaller id");
      else
      {
        for (int ind1 = 0; ind1 < numdof / 2; ind1++)
        {
          stc_matrix(ind1, ind1) +=
              factor1 / number_adj_elements[(ind1) % Core::FE::dim<celltype>] * cond0.size();

          stc_matrix(ind1 + numdof / 2, ind1 + numdof / 2) +=
              factor1 / number_adj_elements[(ind1 + numdof / 2) % Core::FE::dim<celltype>] *
              cond1.size();
          stc_matrix(ind1, ind1 + numdof / 2) +=
              factor2 / number_adj_elements[(ind1) % Core::FE::dim<celltype>] * cond0.size();
          stc_matrix(ind1 + numdof / 2, ind1) +=
              factor2 / number_adj_elements[(ind1 + numdof / 2) % Core::FE::dim<celltype>] *
              cond1.size();
        }
      }
    }
    // surface with nodes 0..3 is to be scaled
    else if (current0)
    {
      // but not by this element
      if (condnum1 > condnum0)
      {
        for (int ind1 = numdof / 2; ind1 < numdof; ind1++)
        {
          stc_matrix(ind1, ind1) += 1.0 / number_adj_elements[(ind1) % Core::FE::dim<celltype>];
        }
      }
      // this element has to do the whole scaling
      else if (condnum1 <= condnum0)
      {
        for (int ind1 = 0; ind1 < numdof; ind1++)
        {
          if (ind1 < numdof / 2)
          {
            stc_matrix(ind1, ind1) +=
                factor3 / number_adj_elements[(ind1) % Core::FE::dim<celltype>] * cond0.size();
            stc_matrix(ind1, ind1 + numdof / 2) +=
                factor4 / number_adj_elements[(ind1) % Core::FE::dim<celltype>] * cond0.size();
          }
          else
          {
            stc_matrix(ind1, ind1) += 1.0 / number_adj_elements[(ind1) % Core::FE::dim<celltype>];
          }
        }
      }
    }
    // surface with nodes 4..7 is to be scaled
    else if (current1)
    {
      // but not by this element
      if (condnum0 > condnum1)
      {
        for (int ind1 = 0; ind1 < numdof / 2; ind1++)
        {
          stc_matrix(ind1, ind1) += 1.0 / number_adj_elements[(ind1) % Core::FE::dim<celltype>];
        }
      }
      // this element has to do the whole scaling
      else if (condnum0 <= condnum1)
      {
        for (int ind1 = 0; ind1 < numdof; ind1++)
        {
          if (ind1 >= numdof / 2)
          {
            stc_matrix(ind1, ind1) +=
                factor3 / number_adj_elements[(ind1) % Core::FE::dim<celltype>] * cond1.size();
            stc_matrix(ind1, -numdof / 2 + ind1) +=
                factor4 / number_adj_elements[(ind1) % Core::FE::dim<celltype>] * cond1.size();
          }
          else
          {
            stc_matrix(ind1, ind1) += 1.0 / number_adj_elements[(ind1) % Core::FE::dim<celltype>];
          }
        }
      }
    }
    else
    {
      for (int ind1 = 0; ind1 < numdof; ind1++)
      {
        stc_matrix(ind1, ind1) += 1.0 / number_adj_elements[(ind1) % Core::FE::dim<celltype>];
      }
    }

    return stc_matrix;
  }

  void evaluate_and_assemble_stc_matrix(Core::FE::Discretization& discret,
      Core::LinAlg::Vector<double> global_displacements, Inpar::Solid::StcScale stc_scale,
      int stc_layer, Core::LinAlg::SparseMatrix& global_matrix)
  {
    std::vector<const Core::Conditions::Condition*> stc_conditions;
    discret.get_condition("STC Layer", stc_conditions);
    auto stc_layer_conditions = Core::Conditions::find_conditioned_node_ids_and_conditions(
        discret, stc_conditions, Core::Conditions::LookFor::locally_owned_and_ghosted);

    discret.evaluate(
        [&](Core::Elements::Element& ele)
        {
          FOUR_C_ASSERT_ALWAYS(
              ele.shape() == Core::FE::CellType::hex8, "STC only supported for hex8 cell type");
          Core::Elements::LocationArray la(1);
          ele.location_vector(discret, la);

          const std::vector<double> displacements =
              Core::FE::extract_values(global_displacements, la[0].lm_);

          Core::LinAlg::SerialDenseMatrix stc_matrix =
              evaluate_stc_matrix<Core::FE::CellType::hex8>(
                  ele, displacements, stc_scale, stc_layer, stc_layer_conditions);

          global_matrix.assemble(ele.id(), stc_matrix, la[0].lm_, la[0].lmowner_);
        });

    global_matrix.complete();
  }
}  // namespace

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Solid::Nln::LinSystem::StcScaling::StcScaling(
    const Solid::TimeInt::BaseDataSDyn& DataSDyn, Solid::TimeInt::BaseDataGlobalState& GState)
    : stcscale_(DataSDyn.get_stc_algo_type()), stclayer_(DataSDyn.get_stc_layer()), stcmat_(nullptr)
{
  // prepare matrix for scaled thickness business of thin shell structures
  stcmat_ =
      std::make_shared<Core::LinAlg::SparseMatrix>(*GState.dof_row_map_view(), 81, true, true);
  stcmat_->zero();


  Core::LinAlg::Vector<double> disp_with_ghosted(*GState.get_discret()->dof_col_map(), true);
  Core::LinAlg::export_to(*GState.get_dis_np(), disp_with_ghosted);

  evaluate_and_assemble_stc_matrix(
      *GState.get_discret(), disp_with_ghosted, stcscale_, 1, *stcmat_);

  for (int lay = 2; lay <= stclayer_; ++lay)
  {
    Core::LinAlg::SparseMatrix tmpstcmat(*GState.dof_row_map_view(), 81, true, true);
    evaluate_and_assemble_stc_matrix(
        *GState.get_discret(), disp_with_ghosted, stcscale_, lay, tmpstcmat);
    stcmat_ = Core::LinAlg::matrix_multiply(tmpstcmat, false, *stcmat_, false, true, false, true);
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Solid::Nln::LinSystem::StcScaling::compute_scaling(const NOX::Nln::LinearProblem& problem)
{
  (void)problem;  // avoid unused parameter warning
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Solid::Nln::LinSystem::StcScaling::scale_linear_system(NOX::Nln::LinearProblem& problem)
{
  // get stiffness matrix
  auto stiff_linalg = std::dynamic_pointer_cast<Core::LinAlg::SparseMatrix>(problem.jac);

  // right multiplication of stiffness matrix
  stiff_scaled_ =
      Core::LinAlg::matrix_multiply(*stiff_linalg, false, *stcmat_, false, true, false, true);

  // left multiplication of stiffness matrix and rhs
  if (stcscale_ == Inpar::Solid::stc_currsym)
  {
    stiff_scaled_ =
        Core::LinAlg::matrix_multiply(*stcmat_, true, *stiff_scaled_, false, true, false, true);

    Core::LinAlg::Vector<double> rhs_scaled(problem.rhs->get_map(), true);
    stcmat_->multiply(true, *problem.rhs, rhs_scaled);
    problem.rhs->update(1.0, rhs_scaled, 0.0);
  }

  // set new stiffness matrix
  problem.jac = stiff_scaled_;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Solid::Nln::LinSystem::StcScaling::unscale_linear_system(NOX::Nln::LinearProblem& problem)
{
  Core::LinAlg::Vector<double> lhs_unscaled(problem.lhs->get_map(), true);

  stcmat_->multiply(false, *problem.lhs, lhs_unscaled);
  problem.lhs->update(1.0, lhs_unscaled, 0.0);
}

FOUR_C_NAMESPACE_CLOSE
