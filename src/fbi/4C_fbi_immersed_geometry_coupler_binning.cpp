// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_fbi_immersed_geometry_coupler_binning.hpp"

#include "4C_beam3_base.hpp"
#include "4C_beaminteraction_calc_utils.hpp"
#include "4C_binstrategy.hpp"
#include "4C_binstrategy_utils.hpp"
#include "4C_fem_discretization_faces.hpp"
#include "4C_fem_general_element.hpp"
#include "4C_fem_general_extract_values.hpp"
#include "4C_fem_general_node.hpp"
#include "4C_linalg_fixedsizematrix.hpp"
#include "4C_linalg_map.hpp"

#include <Teuchos_TimeMonitor.hpp>

FOUR_C_NAMESPACE_OPEN
/*----------------------------------------------------------------------*/

FBI::FBIBinningGeometryCoupler::FBIBinningGeometryCoupler()
    : FBIGeometryCoupler::FBIGeometryCoupler(), binstrategy_(), bintoelemap_(), binrowmap_(nullptr)
{
}
/*----------------------------------------------------------------------*/
void FBI::FBIBinningGeometryCoupler::setup_binning(
    std::vector<std::shared_ptr<Core::FE::Discretization>>& discretizations,
    std::shared_ptr<const Core::LinAlg::Vector<double>> structure_displacement)
{
  std::shared_ptr<const Core::LinAlg::Vector<double>> disp2 =
      std::make_shared<Core::LinAlg::Vector<double>>(*(discretizations[1]->dof_col_map()));

  std::vector<std::shared_ptr<const Core::LinAlg::Vector<double>>> disp_vec = {
      structure_displacement, disp2};

  partition_geometry(discretizations, structure_displacement);
}
/*----------------------------------------------------------------------*/
void FBI::FBIBinningGeometryCoupler::partition_geometry(
    std::vector<std::shared_ptr<Core::FE::Discretization>>& discretizations,
    std::shared_ptr<const Core::LinAlg::Vector<double>> structure_displacement)
{
  std::shared_ptr<const Core::LinAlg::Vector<double>> disp2 =
      std::make_shared<Core::LinAlg::Vector<double>>(*(discretizations[1]->dof_col_map()));

  std::vector<std::shared_ptr<const Core::LinAlg::Vector<double>>> disp_vec = {
      structure_displacement, disp2};

  // nodes, that are owned by a proc, are distributed to the bins of this proc
  std::vector<std::map<int, std::vector<int>>> nodesinbin(2);

  std::map<int, std::set<int>> bintorowelemap_fluid;

  binstrategy_->distribute_elements_to_bins_using_ele_aabb(*discretizations[0],
      discretizations[0]->my_row_element_range(), bintoelemap_, structure_displacement);

  binstrategy_->bin_discret()->fill_complete(Core::FE::OptionsFillComplete::none());

  std::set<int> colbins;

  // first, add default one layer ghosting

  std::vector<int> binvec(27);
  for (auto i = 0; i < binstrategy_->bin_discret()->element_row_map()->num_my_elements(); ++i)
  {
    auto currbin = binstrategy_->bin_discret()->l_row_element(i);
    int it = currbin->id();
    {
      binstrategy_->get_neighbor_and_own_bin_ids(it, binvec);
      colbins.insert(binvec.begin(), binvec.end());
      binvec.clear();
    }
  }


  // extend ghosting of bin discretization
  binstrategy_->extend_ghosting_of_binning_discretization(*binrowmap_, colbins, true);

  // assign Elements to bins
  binstrategy_->remove_all_eles_from_bins();
  binstrategy_->assign_eles_to_bins(*discretizations[0], bintoelemap_,
      BeamInteraction::Utils::convert_element_to_bin_content_type);
}
/*----------------------------------------------------------------------*/
void FBI::FBIBinningGeometryCoupler::update_binning(
    std::shared_ptr<Core::FE::Discretization>& structure_discretization,
    std::shared_ptr<const Core::LinAlg::Vector<double>> structure_column_displacement)
{
  binstrategy_->distribute_elements_to_bins_using_ele_aabb(*structure_discretization,
      structure_discretization->my_col_element_range(), bintoelemap_,
      structure_column_displacement);


  // assign Elements to bins
  binstrategy_->remove_all_eles_from_bins();
  binstrategy_->assign_eles_to_bins(*structure_discretization, bintoelemap_,
      BeamInteraction::Utils::convert_element_to_bin_content_type);
}
/*----------------------------------------------------------------------*/
void FBI::FBIBinningGeometryCoupler::setup(
    std::vector<std::shared_ptr<Core::FE::Discretization>>& discretizations,
    std::shared_ptr<const Core::LinAlg::Vector<double>> structure_displacement)
{
  auto t = Teuchos::TimeMonitor::getNewTimer("FBI::FBICoupler::Setup");
  Teuchos::TimeMonitor monitor(*t);

  setup_binning(discretizations, structure_displacement);

  FBI::FBIGeometryCoupler::setup(discretizations, structure_displacement);
}
/*----------------------------------------------------------------------*/

std::shared_ptr<std::map<int, std::vector<int>>> FBI::FBIBinningGeometryCoupler::search(
    std::vector<std::shared_ptr<Core::FE::Discretization>>& discretizations,
    std::shared_ptr<const Core::LinAlg::Vector<double>>& column_structure_displacement)
{
  auto t = Teuchos::TimeMonitor::getNewTimer("FBI::FBIBinningCoupler::Search");
  Teuchos::TimeMonitor monitor(*t);

  update_binning(discretizations[0], column_structure_displacement);
  // Vector to hand elements pointers to the bridge object
  std::shared_ptr<std::map<int, std::vector<int>>> pairids =
      std::make_shared<std::map<int, std::vector<int>>>();

  pairids = FBI::FBIGeometryCoupler::search(discretizations, column_structure_displacement);

  return pairids;
}

/*----------------------------------------------------------------------*/

void FBI::FBIBinningGeometryCoupler::compute_current_positions(Core::FE::Discretization& dis,
    std::shared_ptr<std::map<int, Core::LinAlg::Matrix<3, 1>>> positions,
    std::shared_ptr<const Core::LinAlg::Vector<double>> disp) const
{
  positions->clear();
  std::vector<int> src_dofs(
      9);  // todo this does not work for all possible elements, does it? Variable size?

  const Core::LinAlg::Map* bincolmap = binstrategy_->bin_discret()->element_col_map();
  std::vector<int> colbinvec;
  colbinvec.reserve(bincolmap->num_my_elements());

  for (int lid = 0; lid < bincolmap->num_my_elements(); ++lid)
  {
    Core::Elements::Element* currbin = binstrategy_->bin_discret()->l_col_element(lid);
    colbinvec.push_back(currbin->id());
  }

  std::set<Core::Elements::Element*> beam_element_list;

  binstrategy_->get_bin_content(
      beam_element_list, {Core::Binstrategy::Utils::BinContentType::Beam}, colbinvec, false);

  for (std::set<Core::Elements::Element*>::iterator element = beam_element_list.begin();
      element != beam_element_list.end(); element++)
  {
    Core::Nodes::Node** node_list = (*element)->nodes();
    unsigned int numnode = (*element)->num_node();
    for (unsigned int i = 0; i < numnode; i++)
    {
      const Core::Nodes::Node* node = node_list[i];
      if (disp != nullptr)
      {
        // get the DOF numbers of the current node
        dis.dof(node, 0, src_dofs);
        // get the current displacements
        std::vector<double> mydisp = Core::FE::extract_values(*disp, src_dofs);

        for (int d = 0; d < 3; ++d) (*positions)[node->id()](d) = node->x()[d] + mydisp.at(d);
      }
    }
  }
}

/*----------------------------------------------------------------------*/

void FBI::FBIBinningGeometryCoupler::set_binning(
    std::shared_ptr<Core::Binstrategy::BinningStrategy> binning)
{
  binstrategy_ = binning;
  binstrategy_->bin_discret()->fill_complete(Core::FE::OptionsFillComplete::none());
  binrowmap_ =
      std::make_shared<Core::LinAlg::Map>(*(binstrategy_->bin_discret()->element_row_map()));
};

FOUR_C_NAMESPACE_CLOSE
