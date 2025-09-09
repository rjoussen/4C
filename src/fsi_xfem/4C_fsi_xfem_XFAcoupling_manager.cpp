// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_fsi_xfem_XFAcoupling_manager.hpp"

#include "4C_adapter_ale_fpsi.hpp"
#include "4C_adapter_str_structure.hpp"
#include "4C_fluid_xfluid.hpp"
#include "4C_linalg_mapextractor.hpp"
#include "4C_xfem_condition_manager.hpp"

FOUR_C_NAMESPACE_OPEN

/*-----------------------------------------------------------------------------------------*
| Constructor                                                                 ager 06/2016 |
*-----------------------------------------------------------------------------------------*/
XFEM::XfaCouplingManager::XfaCouplingManager(std::shared_ptr<FLD::XFluid> xfluid,
    std::shared_ptr<Adapter::AleFpsiWrapper> ale, std::vector<int> idx,
    std::shared_ptr<Adapter::Structure> structure)
    : CouplingCommManager(xfluid->discretization(), ale->discretization(), "", 0, 3),
      ale_(ale),
      xfluid_(xfluid),
      idx_(idx),
      structure_(structure)
{
  if (idx_.size() != static_cast<std::size_t>(2 + (structure_ != nullptr)))
    FOUR_C_THROW("XFACoupling_Manager required (two + num coupled block) ( {} != {})",
        (2 + (structure_ != nullptr)), idx_.size());

  if (structure_ != nullptr)
  {
    if (ale_->discretization()->has_condition("StructAleCoupling") &&
        structure_->discretization()->has_condition("StructAleCoupling"))
    {
      std::vector<const Core::Conditions::Condition*> conds;
      ale_->discretization()->get_condition("StructAleCoupling", conds);
      const int num_nodes_ale = conds.front()->get_nodes()->size();
      structure_->discretization()->get_condition("StructAleCoupling", conds);
      const int num_nodes_struct = conds.front()->get_nodes()->size();

      if (num_nodes_ale != num_nodes_struct)
      {
        FOUR_C_THROW("XFACoupling_Manager: For StructAleCoupling NumNodes not equal! ({} != {})",
            num_nodes_ale, num_nodes_struct);
      }

      std::cout << "|== XFACoupling_Manager: Setup of Ale Structure Coupling! ==|" << std::endl;
      ale_struct_coupling_ = std::make_shared<XFEM::CouplingCommManager>(
          ale_->discretization(), structure_->discretization(), "StructAleCoupling");
    }
  }
}

void XFEM::XfaCouplingManager::predict_coupling_states()
{
  /*
    if (Ale_Struct_coupling_ != nullptr)
    {
      std::cout << "XFEM::XFACoupling_Manager::predict_coupling_states"<< std::endl;

      //-------------------------------------------
      // Perform a 1st predictor to the ALE field (required, otherwise relaxation solve is bad
    conditioned afterwards) ale_->DoPredictor(); // const vel predictor on dispnp_
      //-------------------------------------------

      //-------------------------------------------
      // Perform a 2nd predictor in the sense of a relaxation solve (required, otherwise the ALE
    displacements are too large at the interface)
      //-------------------------------------------

      // first manipulate the fluid-solid interface, then setup the ale system
      Ale_Struct_coupling_->insert_vector(1,structure_->Dispnp(),0,ale_->WriteAccessDispnp(),XFEM::Coupling_Comm_Manager::full_to_full);

      // apply inner Dirichlet conditions (don't forget to reset the time and the step!)
      ale_->prepare_time_step(); // applies DBCs to the current dispnp

      ale_->TimeStep(ALE::Utils::MapExtractor::dbc_set_structale);

      // Reset the time and dt to be called incremented again in the actual ale->prepare_time_step
      ale_->reset_time(ale_->Dt());
    }
  */
}
/*-----------------------------------------------------------------------------------------*
| Set required displacement & velocity states in the coupling object          ager 06/2016 |
*-----------------------------------------------------------------------------------------*/
void XFEM::XfaCouplingManager::set_coupling_states()
{
  // 1 Sets structural conditioned Dispnp onto Ale
  if (ale_struct_coupling_ != nullptr)
    ale_struct_coupling_->insert_vector(1, structure_->dispnp(), 0, ale_->write_access_dispnp(),
        XFEM::CouplingCommManager::full_to_full);

  // 2 Get AleDisplacements
  std::shared_ptr<Core::LinAlg::Vector<double>> aledisplacements =
      std::make_shared<Core::LinAlg::Vector<double>>(*get_map_extractor(0)->map(1), true);
  insert_vector(1, ale_->dispnp(), 0, aledisplacements, CouplingCommManager::partial_to_partial);
  // 3 Set Fluid Dispnp
  get_map_extractor(0)->insert_vector(*aledisplacements, 1, *xfluid_->write_access_dispnp());

  // 4 new grid velocity
  xfluid_->update_gridv();

  // update also ALE vectors w.r.t. current state
  xfluid_->update_ale_state_vectors();

  return;
}

/*-----------------------------------------------------------------------------------------*
| Add the coupling matrixes to the global systemmatrix                        ager 06/2016 |
*-----------------------------------------------------------------------------------------*/
void XFEM::XfaCouplingManager::add_coupling_matrix(
    Core::LinAlg::BlockSparseMatrixBase& systemmatrix, double scaling)
{
  // Get Idx of fluid and ale field map extractors
  const int& aidx_other = ALE::Utils::MapExtractor::cond_other;
  const std::shared_ptr<Core::LinAlg::BlockSparseMatrixBase> a = ale_->block_system_matrix();

  // ALE Condensation
  Core::LinAlg::SparseMatrix& aii = a->matrix(aidx_other, aidx_other);

  systemmatrix.assign(idx_[1], idx_[1], Core::LinAlg::DataAccess::Share, aii);

  if (ale_struct_coupling_ != nullptr)
  {
    const int& aidx_as = ALE::Utils::MapExtractor::cond_lung_asi;
    Core::LinAlg::SparseMatrix& ai_gau = a->matrix(aidx_other, aidx_as);
    ale_struct_coupling_->insert_matrix(0, 0, ai_gau, 1, systemmatrix.matrix(idx_[1], idx_[2]),
        XFEM::CouplingCommManager::col, 1.0, true, false);
  }

  //  //////////////////////////////////////////////
  //  //////                                  //////
  //  //////    Linearization of fluid_field   //////
  //  //////    with respect to ale mesh      //////
  //  //////             motion               //////
  //  //////                                  //////
  //  //////////////////////////////////////////////

  // TODO: THIS IS STILL MISSING, BUT USUALLY DOES NOT HAVE A BIG INFLUENCE ON THE CONVERGENCE!!!
  return;
}

/*-----------------------------------------------------------------------------------------*
| Add the coupling rhs                                                        ager 06/2016 |
*-----------------------------------------------------------------------------------------*/
void XFEM::XfaCouplingManager::add_coupling_rhs(std::shared_ptr<Core::LinAlg::Vector<double>> rhs,
    const Core::LinAlg::MultiMapExtractor& me, double scaling)
{
  std::shared_ptr<const Core::LinAlg::Vector<double>> av = ale_->rhs();
  std::shared_ptr<Core::LinAlg::Vector<double>> aov = ale_->interface()->extract_other_vector(*av);
  me.insert_vector(*aov, idx_[1], *rhs);  // add ALE contributions to 'rhs'
  return;
}

FOUR_C_NAMESPACE_CLOSE
