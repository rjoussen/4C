// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_scatra_timint_heterogeneous_reaction_strategy.hpp"

#include "4C_fem_discretization.hpp"
#include "4C_fem_dofset_gidbased_wrapper.hpp"
#include "4C_fem_dofset_merged_wrapper.hpp"
#include "4C_fem_general_utils_createdis.hpp"
#include "4C_global_data.hpp"
#include "4C_linalg_sparsematrix.hpp"
#include "4C_linear_solver_method_linalg.hpp"
#include "4C_rebalance_print.hpp"
#include "4C_scatra_ele.hpp"
#include "4C_scatra_ele_action.hpp"
#include "4C_scatra_timint_implicit.hpp"
#include "4C_scatra_utils_clonestrategy.hpp"
#include "4C_utils_parameter_list.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 | constructor                                               vuong 06/16 |
 *----------------------------------------------------------------------*/
ScaTra::HeterogeneousReactionStrategy::HeterogeneousReactionStrategy(
    ScaTra::ScaTraTimIntImpl* scatratimint)
    : MeshtyingStrategyStd(scatratimint), issetup_(false), isinit_(false)
{
  return;
}  // ScaTra::HeterogeneousReactionStrategy::HeterogeneousReactionStrategy


/*------------------------------------------------------------------------*
 | evaluate heterogeneous reactions (actually no mesh tying    vuong 06/16 |
 *------------------------------------------------------------------------*/
void ScaTra::HeterogeneousReactionStrategy::evaluate_meshtying()
{
  check_is_init();
  check_is_setup();

  // create parameter list
  Teuchos::ParameterList condparams;

  // action for elements
  Core::Utils::add_enum_class_to_parameter_list<ScaTra::Action>(
      "action", ScaTra::Action::calc_heteroreac_mat_and_rhs, condparams);

  // set global state vectors according to time-integration scheme
  discret_->set_state("phinp", *scatratimint_->phiafnp());
  discret_->set_state("hist", *scatratimint_->hist());

  // provide scatra discretization with convective velocity
  discret_->set_state(scatratimint_->nds_vel(), "convective velocity field",
      *scatratimint_->discretization()->get_state(
          scatratimint_->nds_vel(), "convective velocity field"));

  // provide scatra discretization with velocity
  discret_->set_state(scatratimint_->nds_vel(), "velocity field",
      *scatratimint_->discretization()->get_state(scatratimint_->nds_vel(), "velocity field"));

  if (scatratimint_->is_ale())
  {
    discret_->set_state(scatratimint_->nds_disp(), "dispnp",
        *scatratimint_->discretization()->get_state(scatratimint_->nds_disp(), "dispnp"));
  }

  discret_->evaluate(condparams, scatratimint_->system_matrix(), scatratimint_->residual());

  // now we clear all states.
  // it would be nicer to do this directly before all
  // states are set at the beginning of this method.
  // However, in this case we are not able to set states externally
  // before this method is called.
  // See the call hierarchy of HeterogeneousReactionStrategy::set_state()
  // to check in which algorithms states are set on discret_ .
  discret_->clear_state();
  return;
}  // ScaTra::HeterogeneousReactionStrategy::evaluate_meshtying


/*----------------------------------------------------------------------*
 | initialize meshtying objects                              rauch 09/16 |
 *----------------------------------------------------------------------*/
void ScaTra::HeterogeneousReactionStrategy::setup_meshtying()
{
  // call init() of base class
  ScaTra::MeshtyingStrategyStd::setup_meshtying();

  // make sure we set up everything properly
  heterogeneous_reaction_sanity_check();

  MPI_Comm com(scatratimint_->discretization()->get_comm());

  // standard case
  discret_ = std::make_shared<Core::FE::Discretization>(
      scatratimint_->discretization()->name(), com, Global::Problem::instance()->n_dim());

  // call complete without assigning degrees of freedom
  discret_->fill_complete({
      .assign_degrees_of_freedom = false,
      .init_elements = true,
      .do_boundary_conditions = false,
  });

  std::shared_ptr<Core::FE::Discretization> scatradis = scatratimint_->discretization();

  // create scatra elements if the scatra discretization is empty
  {
    // fill scatra discretization by cloning fluid discretization
    Core::FE::clone_discretization_from_condition<ScaTra::ScatraReactionCloneStrategy>(*scatradis,
        *discret_, "ScatraHeteroReactionSlave",
        Global::Problem::instance()->cloning_material_map());

    // set implementation type of cloned scatra elements
    for (int i = 0; i < discret_->num_my_col_elements(); ++i)
    {
      Discret::Elements::Transport* element =
          dynamic_cast<Discret::Elements::Transport*>(discret_->l_col_element(i));
      if (element == nullptr) FOUR_C_THROW("Invalid element type!");

      if (element->material()->material_type() == Core::Materials::m_matlist_reactions)
        element->set_impl_type(Inpar::ScaTra::impltype_advreac);
      else
        FOUR_C_THROW("Invalid material type for HeterogeneousReactionStrategy!");
    }  // loop over all column elements
  }

  {
    // build a dofset that merges the DOFs from both sides
    // convention: the order of the merged dofset will be
    //  _            _
    // | slave dofs   |
    // |              |
    // |_master dofs _|
    //
    // slave side is supposed to be the surface discretization
    //
    std::shared_ptr<Core::DOFSets::DofSetMergedWrapper> newdofset =
        std::make_shared<Core::DOFSets::DofSetMergedWrapper>(scatradis->get_dof_set_proxy(),
            scatradis, "ScatraHeteroReactionMaster", "ScatraHeteroReactionSlave");

    // assign the dofset to the reaction discretization
    discret_->replace_dof_set(newdofset, false);

    // add all secondary dofsets as proxies
    for (int ndofset = 1; ndofset < scatratimint_->discretization()->num_dof_sets(); ++ndofset)
    {
      std::shared_ptr<Core::DOFSets::DofSetGIDBasedWrapper> gidmatchingdofset =
          std::make_shared<Core::DOFSets::DofSetGIDBasedWrapper>(scatratimint_->discretization(),
              scatratimint_->discretization()->get_dof_set_proxy(ndofset));
      discret_->add_dof_set(gidmatchingdofset);
    }

    // done. Rebuild all maps and boundary condition geometries
    discret_->fill_complete();

    if (Core::Communication::my_mpi_rank(com) == 0 and Core::Communication::num_mpi_ranks(com) > 1)
      std::cout << "parallel distribution of auxiliary discr. with standard ghosting" << std::endl;
    Core::Rebalance::print_parallel_distribution(*discret_);
  }

  set_is_setup(true);
  return;
}


/*----------------------------------------------------------------------*
 | setup meshtying objects                                  vuong 06/16 |
 *----------------------------------------------------------------------*/
void ScaTra::HeterogeneousReactionStrategy::init_meshtying()
{
  set_is_setup(false);

  // call init() of base class
  ScaTra::MeshtyingStrategyStd::init_meshtying();

  set_is_init(true);
  return;
}


/*----------------------------------------------------------------------*
 | Evaluate conditioned elements                            rauch 08/16 |
 *----------------------------------------------------------------------*/
void ScaTra::HeterogeneousReactionStrategy::evaluate_condition(Teuchos::ParameterList& params,
    std::shared_ptr<Core::LinAlg::SparseOperator> systemmatrix1,
    std::shared_ptr<Core::LinAlg::SparseOperator> systemmatrix2,
    std::shared_ptr<Core::LinAlg::Vector<double>> systemvector1,
    std::shared_ptr<Core::LinAlg::Vector<double>> systemvector2,
    std::shared_ptr<Core::LinAlg::Vector<double>> systemvector3, const std::string& condstring,
    const int condid)
{
  check_is_init();
  check_is_setup();

  // Call evaluate_condition on auxiliary discretization.
  // This condition has all dofs, both from the volume-
  // bound scalars and from the surface-bound scalars.
  discret_->evaluate_condition(params, systemmatrix1, systemmatrix2, systemvector1, systemvector2,
      systemvector3, condstring, condid);

  return;
}


/*----------------------------------------------------------------------*
 | Set state on auxiliary discretization                    rauch 12/16 |
 *----------------------------------------------------------------------*/
void ScaTra::HeterogeneousReactionStrategy::set_state(unsigned nds, const std::string& name,
    std::shared_ptr<const Core::LinAlg::Vector<double>> state)
{
  discret_->set_state(nds, name, *state);
  return;
}


/*----------------------------------------------------------------------*
 | sanity check for some assumptions and conventions        rauch 06/17 |
 *----------------------------------------------------------------------*/
void ScaTra::HeterogeneousReactionStrategy::heterogeneous_reaction_sanity_check()
{
  bool valid_slave = false;

  MPI_Comm com = scatratimint_->discretization()->get_comm();

  if (Core::Communication::my_mpi_rank(com) == 0)
    std::cout << " Sanity check for HeterogeneousReactionStrategy ...";

  std::vector<const Core::Conditions::Condition*> slave_conds;
  scatratimint_->discretization()->get_condition("ScatraHeteroReactionSlave", slave_conds);
  FOUR_C_ASSERT_ALWAYS(slave_conds.size() == 1, "Expected exactly one slave condition!");
  const auto* slave_cond = slave_conds[0];

  const Core::LinAlg::Map* element_row_map = scatratimint_->discretization()->element_row_map();

  // loop over row elements
  for (int lid = 0; lid < element_row_map->num_my_elements(); lid++)
  {
    const int gid = element_row_map->gid(lid);

    Core::Elements::Element* ele = scatratimint_->discretization()->g_element(gid);
    Core::Nodes::Node** nodes = ele->nodes();
    if (ele->shape() == Core::FE::CellType::quad4 or ele->shape() == Core::FE::CellType::tri3)
    {
      for (int node = 0; node < ele->num_node(); node++)
      {
        const int node_gid = nodes[node]->id();

        if (not slave_cond->contains_node(node_gid))
        {
          FOUR_C_THROW(
              "Surface discretization for membrane transport is "
              "supposed to wear ScatraHeteroReactionSlave condition!");
        }
        else
        {
          valid_slave = true;
          break;
        }
      }  // loop over nodes of row ele

      if (valid_slave) break;
    }  // if surface transport element

    else if (ele->shape() == Core::FE::CellType::hex8 or ele->shape() == Core::FE::CellType::tet4)
    {
      // no check so far
    }  // if volume transport element

    else
    {
      FOUR_C_THROW(
          "please implement check for new combination of volume transport "
          "- surface transport elements.");
    }

  }  // loop over row elements


  Core::Communication::barrier(com);
  if (Core::Communication::my_mpi_rank(com) == 0) std::cout << " Passed." << std::endl;

  return;
}

FOUR_C_NAMESPACE_CLOSE
