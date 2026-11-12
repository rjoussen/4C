// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_coupling_volmortar.hpp"

#include "4C_coupling_volmortar_cell.hpp"
#include "4C_coupling_volmortar_defines.hpp"
#include "4C_coupling_volmortar_integrator.hpp"
#include "4C_coupling_volmortar_utils.hpp"
#include "4C_cut_cutwizard.hpp"
#include "4C_cut_elementhandle.hpp"
#include "4C_cut_volumecell.hpp"
#include "4C_fem_discretization.hpp"
#include "4C_fem_dofset_predefineddofnumber.hpp"
#include "4C_fem_geometry_delaunay_triangulation.hpp"
#include "4C_fem_geometry_searchtree.hpp"
#include "4C_fem_geometry_searchtree_service.hpp"
#include "4C_linalg_mapextractor.hpp"
#include "4C_linalg_sparsematrix.hpp"
#include "4C_linalg_utils_sparse_algebra_assemble.hpp"
#include "4C_linalg_utils_sparse_algebra_math.hpp"
#include "4C_linear_solver_method.hpp"
#include "4C_linear_solver_method_linalg.hpp"
#include "4C_mortar_calc_utils.hpp"
#include "4C_mortar_coupling3d.hpp"
#include "4C_mortar_element.hpp"
#include "4C_mortar_utils.hpp"
#include "4C_utils_exceptions.hpp"
#include "4C_utils_pairedvector.hpp"
#include "4C_utils_parameter_list.hpp"

#include <Teuchos_Time.hpp>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 |  ctor (public)                                            farah 10/13|
 *----------------------------------------------------------------------*/
Coupling::VolMortar::VolMortarCoupl::VolMortarCoupl(int dim,
    std::shared_ptr<Core::FE::Discretization> dis1,  // on Omega_1
    std::shared_ptr<Core::FE::Discretization> dis2,  // on Omega_2
    const Teuchos::ParameterList& volmortar_parameters,
    const Teuchos::ParameterList& cut_parameters,
    std::vector<int>* coupleddof12,  // 2-->1
    std::vector<int>* coupleddof21,  // 1-->2
    std::pair<int, int>* dofset12,   // 2-->1
    std::pair<int, int>* dofset21,   // 1-->2
    std::shared_ptr<Coupling::VolMortar::Utils::DefaultMaterialStrategy>
        materialstrategy  // strategy for element information transfer
    )
    : dim_(dim),
      cut_params_(cut_parameters),
      dis1_(dis1),
      dis2_(dis2),
      materialstrategy_(materialstrategy)
{
  // check
  if (not dis1_->filled() or not dis2_->filled())
    FOUR_C_THROW(
        "ERROR: fill_complete() has to be called on both discretizations before setup of "
        "VolMortarCoupl");

  // its the same communicator for all discr.
  comm_ = dis1->get_comm();
  myrank_ = Core::Communication::my_mpi_rank(comm_);

  // define dof sets
  if (dofset21 == nullptr)
    dofset21_ = std::pair<int, int>(1, 0);
  else
    dofset21_ = *dofset21;

  if (dofset12 == nullptr)
    dofset12_ = std::pair<int, int>(1, 0);
  else
    dofset12_ = *dofset12;

  const int* rownodes1 = dis1->node_row_map()->my_global_elements();
  const int numrownode1 = dis1->node_row_map()->num_my_elements();

  const int* rownodes2 = dis2->node_row_map()->my_global_elements();
  const int numrownode2 = dis2->node_row_map()->num_my_elements();

  // determine row map of projection operators
  build_maps(dis1, p12_dofrowmap_, coupleddof21, rownodes1, numrownode1, dofset12_.first);

  build_maps(dis2, p21_dofrowmap_, coupleddof12, rownodes2, numrownode2, dofset21_.first);

  // determine domain map of projection operators
  build_maps(dis1, p21_dofdomainmap_, coupleddof21, rownodes1, numrownode1, dofset21_.second);
  build_maps(dis2, p12_dofdomainmap_, coupleddof12, rownodes2, numrownode2, dofset12_.second);

  const int* colnodes1 = dis1->node_col_map()->my_global_elements();
  const int numcolnode1 = dis1->node_col_map()->num_my_elements();

  const int* colnodes2 = dis2->node_col_map()->my_global_elements();
  const int numcolnode2 = dis2->node_col_map()->num_my_elements();

  // determine column map of projection operators
  build_maps(dis1, p21_dofcolmap_, coupleddof21, colnodes1, numcolnode1, dofset21_.second);
  build_maps(dis2, p12_dofcolmap_, coupleddof12, colnodes2, numcolnode2, dofset12_.second);

  // get required parameter list
  read_and_check_input(volmortar_parameters);

  // init dop normals
  init_dop_normals();

  // init aux normal TODO: no fixed direction!!! ONLY FOR 2D CASE !!!
  auxn_[0] = 0.0;
  auxn_[1] = 0.0;
  auxn_[2] = 1.0;

  // initialize the counter
  polygoncounter_ = 0;
  cellcounter_ = 0;
  inteles_ = 0;
  volume_ = 0.0;
}

/*----------------------------------------------------------------------*
 |  Build maps based on coupling dofs                        farah 03/15|
 *----------------------------------------------------------------------*/
void Coupling::VolMortar::VolMortarCoupl::build_maps(std::shared_ptr<Core::FE::Discretization>& dis,
    std::shared_ptr<const Core::LinAlg::Map>& dofmap, const std::vector<int>* coupleddof,
    const int* nodes, int numnode, int dofset)
{
  std::vector<int> dofmapvec;
  std::map<int, std::vector<int>> dofs;

  const bool couplalldofs = (coupleddof == nullptr);

  unsigned numcoupleddofs = 0;
  if (not couplalldofs)
  {
    // count number of dofs per node that are coupled
    for (unsigned j = 0; j < coupleddof->size(); ++j)
      if ((*coupleddof)[j] == 1) numcoupleddofs++;
  }
  else
    numcoupleddofs = dis->num_dof(dofset, dis->l_row_node(0));

  for (int i = 0; i < numnode; ++i)
  {
    const Core::Nodes::Node* actnode = dis->g_node(nodes[i]);

    const std::vector<int> dof = dis->dof(dofset, actnode);
    if (numcoupleddofs > dof.size())
      FOUR_C_THROW("got just {} dofs at node {} (lid={}) but expected {}", dof.size(), nodes[i], i,
          numcoupleddofs);

    for (unsigned j = 0; j < dof.size(); ++j)
    {
      if (not couplalldofs)
        if ((*coupleddof)[j] == 0) continue;

      dofmapvec.push_back(dof[j]);
    }
  }
  // dof map is the original, unpermuted distribution of dofs
  dofmap = std::make_shared<Core::LinAlg::Map>(-1, dofmapvec.size(), dofmapvec.data(), 0, comm_);
}

/*----------------------------------------------------------------------*
 |  Evaluate (public)                                        farah 10/13|
 *----------------------------------------------------------------------*/
void Coupling::VolMortar::VolMortarCoupl::evaluate_volmortar()
{
  /***********************************************************
   * Welcome                                                 *
   ***********************************************************/
  if (myrank_ == 0)
  {
    std::cout << "**************************************************" << std::endl;
    std::cout << "*****     Welcome to Core::VOLMORTAR-Coupling!     *****" << std::endl;
    std::cout << "**************************************************" << std::endl;
  }

  // time measurement
  Core::Communication::barrier(comm_);
  const double t_start = Teuchos::Time::wallTime();

  /***********************************************************
   * Check initial residuum and perform mesh init             *
   ***********************************************************/
  // check_initial_residuum();
  // mesh initialization procedure
  if (params().get<bool>("MESH_INIT")) mesh_init();

  /***********************************************************
   * initialize global matrices                              *
   ***********************************************************/
  initialize();

  /***********************************************************
   * Segment-based integration                               *
   ***********************************************************/
  if (Teuchos::getIntegralValue<IntType>(params(), "INTTYPE") == inttype_segments)
    evaluate_segments();

  /***********************************************************
   * Element-based Integration                               *
   ***********************************************************/
  else if (Teuchos::getIntegralValue<IntType>(params(), "INTTYPE") == inttype_elements)
    evaluate_elements();

  // no other possibility
  else
    FOUR_C_THROW("ERROR: Chosen INTTYPE not provided");

  /***********************************************************
   * complete global matrices and create projection operator *
   ***********************************************************/
  complete();
  create_projection_operator();

  /**************************************************
   * Bye                                            *
   **************************************************/
  Core::Communication::barrier(comm_);
  const double inttime = Teuchos::Time::wallTime() - t_start;

  if (myrank_ == 0)
  {
    std::cout << "**************************************************" << std::endl;
    std::cout << "*****       Core::VOLMORTAR-Coupling Done!!!       *****" << std::endl;
    std::cout << "**************************************************" << std::endl;

    // output
    std::cout << "required evaluation time: " << inttime << std::endl;
    //  std::cout << "Polyogns/Polyhedra = " << polygoncounter_ << std::endl;
    //  std::cout << "Created Cells      = " << cellcounter_    << std::endl;
    //  std::cout << "Integr. Elements   = " << inteles_        << std::endl;
    std::cout << "Integrated Volume  = " << volume_ << "\n" << std::endl;
  }

  // reset counter
  polygoncounter_ = 0;
  cellcounter_ = 0;
  inteles_ = 0;
  volume_ = 0.0;
}

/*----------------------------------------------------------------------*
 |  Init normals for Dop calculation                         farah 05/14|
 *----------------------------------------------------------------------*/
void Coupling::VolMortar::VolMortarCoupl::init_dop_normals()
{
  dopnormals_(0, 0) = 1.0;
  dopnormals_(0, 1) = 0.0;
  dopnormals_(0, 2) = 0.0;

  dopnormals_(1, 0) = 0.0;
  dopnormals_(1, 1) = 1.0;
  dopnormals_(1, 2) = 0.0;

  dopnormals_(2, 0) = 0.0;
  dopnormals_(2, 1) = 0.0;
  dopnormals_(2, 2) = 1.0;

  dopnormals_(3, 0) = 1.0;
  dopnormals_(3, 1) = 1.0;
  dopnormals_(3, 2) = 0.0;

  dopnormals_(4, 0) = 1.0;
  dopnormals_(4, 1) = 0.0;
  dopnormals_(4, 2) = 1.0;

  dopnormals_(5, 0) = 0.0;
  dopnormals_(5, 1) = 1.0;
  dopnormals_(5, 2) = 1.0;

  dopnormals_(6, 0) = 1.0;
  dopnormals_(6, 1) = 0.0;
  dopnormals_(6, 2) = -1.0;

  dopnormals_(7, 0) = 1.0;
  dopnormals_(7, 1) = -1.0;
  dopnormals_(7, 2) = 0.0;

  dopnormals_(8, 0) = 0.0;
  dopnormals_(8, 1) = 1.0;
  dopnormals_(8, 2) = -1.0;
}
/*----------------------------------------------------------------------*
 |  Init search tree                                         farah 05/14|
 *----------------------------------------------------------------------*/
std::shared_ptr<Core::Geo::SearchTree> Coupling::VolMortar::VolMortarCoupl::init_search(
    std::shared_ptr<Core::FE::Discretization> searchdis)
{
  // init current positions
  std::map<int, Core::LinAlg::Matrix<3, 1>> currentpositions;

  for (int lid = 0; lid < searchdis->num_my_col_elements(); ++lid)
  {
    Core::Elements::Element* source_ele = searchdis->l_col_element(lid);

    // calculate slabs for every node on every element
    for (int k = 0; k < source_ele->num_node(); k++)
    {
      Core::Nodes::Node* node = source_ele->nodes()[k];
      Core::LinAlg::Matrix<3, 1> currpos;

      const auto x = node->x();
      currpos(0) = x[0];
      currpos(1) = x[1];
      if (x.size() > 2) currpos(2) = x[2];

      currentpositions[node->id()] = currpos;
    }
  }

  // init of 3D search tree
  std::shared_ptr<Core::Geo::SearchTree> searchTree = std::make_shared<Core::Geo::SearchTree>(5);

  // find the bounding box of the elements and initialize the search tree
  const Core::LinAlg::Matrix<3, 2> rootBox =
      Core::Geo::get_xaab_bof_dis(*searchdis, currentpositions);
  searchTree->initialize_tree(rootBox, *searchdis, Core::Geo::TreeType(Core::Geo::OCTTREE));

  return searchTree;
}

/*----------------------------------------------------------------------*
 |  Calculate Dops for background mesh                       farah 05/14|
 *----------------------------------------------------------------------*/
std::map<int, Core::LinAlg::Matrix<9, 2>> Coupling::VolMortar::VolMortarCoupl::calc_background_dops(
    std::shared_ptr<Core::FE::Discretization> searchdis)
{
  std::map<int, Core::LinAlg::Matrix<9, 2>> currentKDOPs;

  for (int lid = 0; lid < searchdis->num_my_col_elements(); ++lid)
  {
    Core::Elements::Element* source_ele = searchdis->l_col_element(lid);

    currentKDOPs[source_ele->id()] = calc_dop(*source_ele);
  }

  return currentKDOPs;
}

/*----------------------------------------------------------------------*
 |  Calculate Dop for one Element                            farah 05/14|
 *----------------------------------------------------------------------*/
Core::LinAlg::Matrix<9, 2> Coupling::VolMortar::VolMortarCoupl::calc_dop(
    Core::Elements::Element& ele)
{
  Core::LinAlg::Matrix<9, 2> dop;

  // calculate slabs
  for (int j = 0; j < 9; j++)
  {
    // initialize slabs
    dop(j, 0) = 1.0e12;
    dop(j, 1) = -1.0e12;
  }

  // calculate slabs for every node on every element
  for (int k = 0; k < ele.num_node(); k++)
  {
    Core::Nodes::Node* node = ele.nodes()[k];

    // get current node position
    std::array<double, 3> pos = {0.0, 0.0, 0.0};
    for (int j = 0; j < dim_; ++j) pos[j] = node->x()[j];

    // calculate slabs
    for (int j = 0; j < 9; j++)
    {
      //= ax+by+cz=d/sqrt(aa+bb+cc)
      double num =
          dopnormals_(j, 0) * pos[0] + dopnormals_(j, 1) * pos[1] + dopnormals_(j, 2) * pos[2];
      double denom =
          sqrt((dopnormals_(j, 0) * dopnormals_(j, 0)) + (dopnormals_(j, 1) * dopnormals_(j, 1)) +
               (dopnormals_(j, 2) * dopnormals_(j, 2)));
      double dcurrent = num / denom;

      if (dcurrent > dop(j, 1)) dop(j, 1) = dcurrent;
      if (dcurrent < dop(j, 0)) dop(j, 0) = dcurrent;
    }
  }

  return dop;
}

/*----------------------------------------------------------------------*
 |  Perform searching procedure                              farah 05/14|
 *----------------------------------------------------------------------*/
std::vector<int> Coupling::VolMortar::VolMortarCoupl::search(Core::Elements::Element& ele,
    std::shared_ptr<Core::Geo::SearchTree> SearchTree,
    std::map<int, Core::LinAlg::Matrix<9, 2>>& currentKDOPs)
{
  // vector of global ids of found elements
  std::vector<int> gids;
  gids.clear();
  std::set<int> gid;
  gid.clear();

  Core::LinAlg::Matrix<9, 2> queryKDOP;

  // calc dop for considered element
  queryKDOP = calc_dop(ele);

  //**********************************************************
  // search for near elements to the background node's coord
  SearchTree->search_collisions(currentKDOPs, queryKDOP, 0, gid);

  for (std::set<int>::iterator iter = gid.begin(); iter != gid.end(); ++iter) gids.push_back(*iter);

  return gids;
}

/*----------------------------------------------------------------------*
 |  Assign materials for both fields                         vuong 09/14|
 *----------------------------------------------------------------------*/
void Coupling::VolMortar::VolMortarCoupl::assign_materials()
{
  if (dis1_ == nullptr or dis2_ == nullptr)
    FOUR_C_THROW("no discretization for assigning materials!");

  // init search trees
  std::shared_ptr<Core::Geo::SearchTree> SearchTreeA = init_search(dis1_);
  std::shared_ptr<Core::Geo::SearchTree> SearchTreeB = init_search(dis2_);

  // calculate DOPs for search algorithm
  std::map<int, Core::LinAlg::Matrix<9, 2>> CurrentDOPsA = calc_background_dops(dis1_);
  std::map<int, Core::LinAlg::Matrix<9, 2>> CurrentDOPsB = calc_background_dops(dis2_);

  /**************************************************
   * loop over all Adis elements                    *
   **************************************************/
  for (int j = 0; j < dis1_->num_my_col_elements(); ++j)
  {
    // get target element
    Core::Elements::Element* Aele = dis1_->l_col_element(j);

    std::vector<int> found = search(*Aele, SearchTreeB, CurrentDOPsB);

    /***********************************************************
     * Assign materials                                        *
     ***********************************************************/
    materialstrategy_->assign_material2_to1(this, Aele, found, dis1_, dis2_);
  }

  /**************************************************
   * loop over all Bdis elements                    *
   **************************************************/
  for (int j = 0; j < dis2_->num_my_col_elements(); ++j)
  {
    // get target element
    Core::Elements::Element* Bele = dis2_->l_col_element(j);

    std::vector<int> found = search(*Bele, SearchTreeA, CurrentDOPsA);

    /***********************************************************
     * Assign materials                                        *
     ***********************************************************/
    materialstrategy_->assign_material1_to2(this, Bele, found, dis1_, dis2_);
  }
}

/*----------------------------------------------------------------------*
 |  Calculate trafo matrix for quadr. elements               farah 05/14|
 *----------------------------------------------------------------------*/
std::vector<int> Coupling::VolMortar::VolMortarCoupl::get_adjacent_nodes(
    Core::FE::CellType shape, int& lid)
{
  // vector of adjacent node ids
  std::vector<int> ids;

  switch (shape)
  {
    case Core::FE::CellType::hex20:
    {
      if (lid == 8)
      {
        ids.push_back(0);
        ids.push_back(1);
      }
      else if (lid == 9)
      {
        ids.push_back(1);
        ids.push_back(2);
      }
      else if (lid == 10)
      {
        ids.push_back(2);
        ids.push_back(3);
      }
      else if (lid == 11)
      {
        ids.push_back(3);
        ids.push_back(0);
      }
      else if (lid == 12)
      {
        ids.push_back(0);
        ids.push_back(4);
      }
      else if (lid == 13)
      {
        ids.push_back(1);
        ids.push_back(5);
      }
      else if (lid == 14)
      {
        ids.push_back(2);
        ids.push_back(6);
      }
      else if (lid == 15)
      {
        ids.push_back(3);
        ids.push_back(7);
      }
      else if (lid == 16)
      {
        ids.push_back(4);
        ids.push_back(5);
      }
      else if (lid == 17)
      {
        ids.push_back(5);
        ids.push_back(6);
      }
      else if (lid == 18)
      {
        ids.push_back(6);
        ids.push_back(7);
      }
      else if (lid == 19)
      {
        ids.push_back(4);
        ids.push_back(7);
      }
      else
        FOUR_C_THROW("ERROR: Given Id is wrong!!!");

      break;
    }
    case Core::FE::CellType::tet10:
    {
      if (lid == 4)
      {
        ids.push_back(0);
        ids.push_back(1);
      }
      else if (lid == 5)
      {
        ids.push_back(1);
        ids.push_back(2);
      }
      else if (lid == 6)
      {
        ids.push_back(0);
        ids.push_back(2);
      }
      else if (lid == 7)
      {
        ids.push_back(0);
        ids.push_back(3);
      }
      else if (lid == 8)
      {
        ids.push_back(1);
        ids.push_back(3);
      }
      else if (lid == 9)
      {
        ids.push_back(2);
        ids.push_back(3);
      }
      else
        FOUR_C_THROW("ERROR: Given Id is wrong!!!");

      break;
    }
    default:
      FOUR_C_THROW("ERROR: shape unknown\n");
      break;
  }

  return ids;
}
/*----------------------------------------------------------------------*
 |  Calculate trafo matrix for quadr. elements               farah 05/14|
 *----------------------------------------------------------------------*/
void Coupling::VolMortar::VolMortarCoupl::create_trafo_operator(Core::Elements::Element& ele,
    std::shared_ptr<Core::FE::Discretization> searchdis, bool dis, std::set<int>& donebefore)
{
  // trafo parameter
  const double alpha = 0.3;

  // ids for corner and edge nodes
  int corner_min = 0;
  int corner_max = 0;
  int edge_min = 0;
  int edge_max = 0;

  if (ele.shape() == Core::FE::CellType::hex20)
  {
    corner_min = 0;
    corner_max = 7;
    edge_min = 8;
    edge_max = 19;
  }
  else if (ele.shape() == Core::FE::CellType::tet10)
  {
    corner_min = 0;
    corner_max = 3;
    edge_min = 4;
    edge_max = 9;
  }
  else
    FOUR_C_THROW("ERROR: Unknown shape for trafo matrix");

  // loop over element nodes
  for (int i = 0; i < ele.num_node(); ++i)
  {
    Core::Nodes::Node* cnode = ele.nodes()[i];
    if (cnode->owner() != myrank_) continue;

    std::set<int>::iterator iter = donebefore.find(cnode->id());
    if (iter != donebefore.end()) continue;

    donebefore.insert(cnode->id());

    //************************************************
    // corner nodes
    if (i >= corner_min and i <= corner_max)
    {
      int nsource_dof = searchdis->num_dof(1, cnode);

      // loop over source dofs
      for (int jdof = 0; jdof < nsource_dof; ++jdof)
      {
        int row = searchdis->dof(1, cnode, jdof);

        if (dis)
          t1_->assemble(1.0, row, row);
        else
          t2_->assemble(1.0, row, row);
      }
    }
    //************************************************
    // edge nodes
    else if (i >= edge_min and i <= edge_max)
    {
      // get ids
      std::vector<int> ids = get_adjacent_nodes(ele.shape(), i);

      int nsource_dof = searchdis->num_dof(1, cnode);

      // loop over dofs
      for (int jdof = 0; jdof < nsource_dof; ++jdof)
      {
        int row = searchdis->dof(1, cnode, jdof);

        // assemble diagonal entries
        if (dis)
          t1_->assemble(1.0 - 2.0 * alpha, row, row);
        else
          t2_->assemble(1.0 - 2.0 * alpha, row, row);

        // found ids
        for (int id = 0; id < (int)ids.size(); ++id)
        {
          Core::Nodes::Node* fnode = ele.nodes()[ids[id]];
          int nfdof = searchdis->num_dof(1, fnode);

          for (int fdof = 0; fdof < nfdof; ++fdof)
          {
            int col = searchdis->dof(1, fnode, fdof);

            if (jdof == fdof)
            {
              // assemble off-diagonal entries
              if (dis)
                t1_->assemble(alpha, row, col);
              else
                t2_->assemble(alpha, row, col);
            }
          }
        }
      }
    }
  }
}

/*----------------------------------------------------------------------*
 |  Consistent interpolation routine                         farah 06/14|
 *----------------------------------------------------------------------*/
void Coupling::VolMortar::VolMortarCoupl::evaluate_consistent_interpolation()
{
  /***********************************************************
   * Welcome                                                 *
   ***********************************************************/
  if (myrank_ == 0)
  {
    std::cout << "*********************************************************" << std::endl;
    std::cout << "***** Welcome to Consistent-Interpolation-Coupling! *****" << std::endl;
    std::cout << "*********************************************************" << std::endl;
  }

  /***********************************************************
   * Init P-matrices                                         *
   ***********************************************************/
  p12_ = std::make_shared<Core::LinAlg::SparseMatrix>(*p12_dofrowmap_, 10);
  p21_ = std::make_shared<Core::LinAlg::SparseMatrix>(*p21_dofrowmap_, 100);

  /***********************************************************
   * create search tree and current dops                     *
   ***********************************************************/
  std::shared_ptr<Core::Geo::SearchTree> SearchTreeA = init_search(dis1_);
  std::shared_ptr<Core::Geo::SearchTree> SearchTreeB = init_search(dis2_);
  std::map<int, Core::LinAlg::Matrix<9, 2>> CurrentDOPsA = calc_background_dops(dis1_);
  std::map<int, Core::LinAlg::Matrix<9, 2>> CurrentDOPsB = calc_background_dops(dis2_);

  /***********************************************************
   * Create P operators                                      *
   ***********************************************************/
  for (int i = 0; i < dis1_->num_my_col_nodes(); ++i)
  {
    // 1 map node into bele
    int gid = dis1_->node_col_map()->gid(i);
    Core::Nodes::Node* anode = dis1_->g_node(gid);

    // get found elements from other discr.
    std::vector<int> found =
        search(*anode->adjacent_elements()[0].user_element(), SearchTreeB, CurrentDOPsB);

    assemble_consistent_interpolation_p12(anode, found);
  }  // end node loop

  //================================================
  for (int i = 0; i < dis2_->num_my_col_nodes(); ++i)
  {
    // 1 map node into bele
    int gid = dis2_->node_col_map()->gid(i);
    Core::Nodes::Node* bnode = dis2_->g_node(gid);

    // get found elements from other discr.
    std::vector<int> found =
        search(*bnode->adjacent_elements()[0].user_element(), SearchTreeA, CurrentDOPsA);

    assemble_consistent_interpolation_p21(bnode, found);
  }  // end node loop

  /***********************************************************
   * Complete                                                *
   ***********************************************************/
  p12_->complete(*p12_dofdomainmap_, *p12_dofrowmap_);
  p21_->complete(*p21_dofdomainmap_, *p21_dofrowmap_);
}

/*----------------------------------------------------------------------*
 |  Element-based routine                                    farah 04/14|
 *----------------------------------------------------------------------*/
void Coupling::VolMortar::VolMortarCoupl::evaluate_elements()
{
  // output
  if (myrank_ == 0)
  {
    std::cout << "*****       Element-based Integration        *****" << std::endl;
    std::cout << "*****       Calc First Projector:            *****" << std::endl;
  }

  // init search trees
  std::shared_ptr<Core::Geo::SearchTree> SearchTreeA = init_search(dis1_);
  std::shared_ptr<Core::Geo::SearchTree> SearchTreeB = init_search(dis2_);

  // calculate DOPs for search algorithm
  std::map<int, Core::LinAlg::Matrix<9, 2>> CurrentDOPsA = calc_background_dops(dis1_);
  std::map<int, Core::LinAlg::Matrix<9, 2>> CurrentDOPsB = calc_background_dops(dis2_);

  /**************************************************
   * loop over all Adis elements                    *
   **************************************************/
  std::set<int> donebeforea;
  for (int j = 0; j < dis1_->num_my_col_elements(); ++j)
  {
    // print_status(j,false);

    // get target element
    Core::Elements::Element* Aele = dis1_->l_col_element(j);

    std::vector<int> found = search(*Aele, SearchTreeB, CurrentDOPsB);
    integrate_3d_ele_based_p12(*Aele, found);

    // create trafo operator for quadr. modification
    if (dualquad_ != dualquad_no_mod) create_trafo_operator(*Aele, dis1_, true, donebeforea);
  }

  // half-time output
  if (myrank_ == 0)
  {
    std::cout << "**************************************************" << std::endl;
    std::cout << "*****       Calc Second Projector:           *****" << std::endl;
  }

  /**************************************************
   * loop over all Bdis elements                    *
   **************************************************/
  std::set<int> donebeforeb;
  for (int j = 0; j < dis2_->num_my_col_elements(); ++j)
  {
    // print_status(j,true);

    // get target element
    Core::Elements::Element* Bele = dis2_->l_col_element(j);

    std::vector<int> found = search(*Bele, SearchTreeA, CurrentDOPsA);
    integrate_3d_ele_based_p21(*Bele, found);

    // create trafo operator for quadr. modification
    if (dualquad_ != dualquad_no_mod) create_trafo_operator(*Bele, dis2_, false, donebeforeb);
  }

  // stats
  inteles_ += dis1_->num_global_elements();
  inteles_ += dis2_->num_global_elements();
}

/*----------------------------------------------------------------------*
 |  Segment-based routine                                    farah 04/14|
 *----------------------------------------------------------------------*/
void Coupling::VolMortar::VolMortarCoupl::evaluate_segments()
{
  // create search tree and current dops
  std::shared_ptr<Core::Geo::SearchTree> SearchTreeB = init_search(dis2_);
  std::map<int, Core::LinAlg::Matrix<9, 2>> CurrentDOPsB = calc_background_dops(dis2_);

  /**************************************************
   * loop over all source elements                   *
   **************************************************/
  for (int i = 0; i < dis1_->num_my_col_elements(); ++i)
  {
    // output of coupling status
    // print_status(i);

    // get source element
    Core::Elements::Element* Aele = dis1_->l_col_element(i);

    // get found elements from other discr.
    std::vector<int> found = search(*Aele, SearchTreeB, CurrentDOPsB);

    /***********************************************************
     * Assign materials                                        *
     ***********************************************************/
    materialstrategy_->assign_material2_to1(this, Aele, found, dis1_, dis2_);

    /**************************************************
     * loop over all target elements                  *
     **************************************************/
    for (int foundeles = 0; foundeles < (int)found.size(); ++foundeles)
    {
      // get b element
      Core::Elements::Element* Bele = dis2_->g_element(found[foundeles]);

      /***********************************************************
       * Assign materials                                        *
       ***********************************************************/
      materialstrategy_->assign_material1_to2(
          this, Bele, std::vector<int>(1, Aele->id()), dis1_, dis2_);

      /**************************************************
       *                    2D                          *
       **************************************************/
      if (dim_ == 2)
      {
        evaluate_segments_2d(*Aele, *Bele);
      }

      /**************************************************
       *                    3D                          *
       **************************************************/
      else if (dim_ == 3)
      {
        evaluate_segments_3d(Aele, Bele);
      }

      /**************************************************
       * Wrong Dimension !!!                            *
       **************************************************/
      else
        FOUR_C_THROW("ERROR: Problem dimension is not correct!");

    }  // end target element loop
  }  // end source element loop
}

/*----------------------------------------------------------------------*
 |  Segment-based routine 2D                                 farah 04/14|
 *----------------------------------------------------------------------*/
void Coupling::VolMortar::VolMortarCoupl::evaluate_segments_2d(
    Core::Elements::Element& Aele, Core::Elements::Element& Bele)
{
  // define polygon vertices
  static std::vector<Mortar::Vertex> source_vertices;
  static std::vector<Mortar::Vertex> target_vertices;
  static std::vector<Mortar::Vertex> ClippedPolygon;
  static std::vector<std::shared_ptr<Mortar::IntCell>> cells;

  // clear old polygons
  source_vertices.clear();
  target_vertices.clear();
  ClippedPolygon.clear();
  cells.clear();
  cells.resize(0);

  // build new polygons
  define_vertices_target(Bele, target_vertices);
  define_vertices_source(Aele, source_vertices);

  double tol = 1e-12;
  ClippedPolygon = Mortar::Utils::clip_polygons(source_vertices, target_vertices, auxn(), tol);
  int clipsize = (int)(ClippedPolygon.size());

  // proceed only if clipping polygon is at least a triangle
  if (clipsize < 3)
    return;
  else
    polygoncounter_ += 1;

  // triangulation
  delaunay_triangulation(cells, ClippedPolygon, tol);

  cellcounter_ += (int)cells.size();

  // integrate cells
  integrate_2d(Aele, Bele, cells);
}

/*----------------------------------------------------------------------*
 |  Segment-based routine 3D                                 farah 04/14|
 *----------------------------------------------------------------------*/
void Coupling::VolMortar::VolMortarCoupl::evaluate_segments_3d(
    Core::Elements::Element* Aele, Core::Elements::Element* Bele)
{
  // check need element-based integration over source_ele:
  bool integrateA = check_ele_integration(*Aele, *Bele);

  // check need element-based integration over target_ele:
  bool integrateB = check_ele_integration(*Bele, *Aele);

  // check need for cut:
  bool performcut = check_cut(*Aele, *Bele);

  /**************************************************
   * Integrate element-based source_ele                   *
   **************************************************/
  if (integrateA)
  {
    // integrate over Aele
    integrate_3d(*Aele, *Bele, 0);
  }
  /**************************************************
   * Integrate element-based target_ele                   *
   **************************************************/
  else if (integrateB)
  {
    // integrate over Bele
    integrate_3d(*Aele, *Bele, 1);
  }
  /**************************************************
   * Start Cut                                      *
   **************************************************/
  else if (performcut)
  {
    // call cut and create integration cells
    try
    {
      bool switched_conf = false;
      perform_cut(Aele, Bele, switched_conf);
    }
    catch (Core::Exception& err1)
    {
      try
      {
        bool switched_conf = true;
        perform_cut(Bele, Aele, switched_conf);
      }
      catch (Core::Exception& err2)
      {
        std::cout << "runtime error 2 = " << err2.what() << std::endl;
      }
    }
  }  // end performcut
  else
  {
    // do nothing...
  }
}

/*----------------------------------------------------------------------*
 |  get parameters and check for validity                    farah 04/14|
 *----------------------------------------------------------------------*/
void Coupling::VolMortar::VolMortarCoupl::read_and_check_input(
    const Teuchos::ParameterList& volmortar_parameters)
{
  // check validity
  if (Teuchos::getIntegralValue<IntType>(volmortar_parameters, "INTTYPE") == inttype_segments)
  {
    if (myrank_ == 0)
    {
      std::cout
          << "WARNING: The chosen integration type for volmortar coupling requires cut procedure !"
          << std::endl;
      std::cout
          << "WARNING: The cut is up to now not able to exactly calculate the required segments!"
          << std::endl;
    }
  }

  if (volmortar_parameters.get<bool>("MESH_INIT") and
      Teuchos::getIntegralValue<IntType>(volmortar_parameters, "INTTYPE") == inttype_segments)
  {
    FOUR_C_THROW("ERROR: mesh_init only for ele-based integration!!!");
  }

  if (Teuchos::getIntegralValue<Coupling::VolMortar::Shapefcn>(volmortar_parameters, "SHAPEFCN") ==
      shape_std)
  {
    std::cout << "WARNING: Standard shape functions are employed! D is lumped!" << std::endl;
  }

  // merge to global parameter list
  params_.setParameters(volmortar_parameters);

  // get specific and frequently reused parameters
  dualquad_ = Teuchos::getIntegralValue<DualQuad>(params_, "DUALQUAD");
}

/*----------------------------------------------------------------------*
 |  check initial residuum                                   farah 04/14|
 *----------------------------------------------------------------------*/
void Coupling::VolMortar::VolMortarCoupl::check_initial_residuum()
{
  // create vectors of initial primary variables
  std::shared_ptr<Core::LinAlg::Vector<double>> var_A =
      std::make_shared<Core::LinAlg::Vector<double>>(*discret1()->dof_row_map(0), true);
  std::shared_ptr<Core::LinAlg::Vector<double>> var_B =
      std::make_shared<Core::LinAlg::Vector<double>>(*discret2()->dof_row_map(1), true);

  // solution
  std::shared_ptr<Core::LinAlg::Vector<double>> result_A =
      std::make_shared<Core::LinAlg::Vector<double>>(*discret2()->dof_row_map(1), true);
  std::shared_ptr<Core::LinAlg::Vector<double>> result_B =
      std::make_shared<Core::LinAlg::Vector<double>>(*discret2()->dof_row_map(1), true);

  // node positions for Discr A
  for (int i = 0; i < discret1()->num_my_row_elements(); ++i)
  {
    Core::Elements::Element* Aele = discret1()->l_row_element(i);
    for (int j = 0; j < Aele->num_node(); ++j)
    {
      Core::Nodes::Node* cnode = Aele->nodes()[j];
      int nsource_dof = discret1()->num_dof(0, cnode);

      // loop over source dofs
      for (int jdof = 0; jdof < nsource_dof; ++jdof)
      {
        int id = discret1()->dof(0, cnode, jdof);
        double val = cnode->x()[jdof];

        var_A->replace_global_value(id, val);
      }
    }
  }

  // node positions for Discr B
  for (int i = 0; i < discret2()->num_my_row_elements(); ++i)
  {
    Core::Elements::Element* Bele = discret2()->l_row_element(i);
    for (int j = 0; j < Bele->num_node(); ++j)
    {
      Core::Nodes::Node* cnode = Bele->nodes()[j];
      int nsource_dof = discret2()->num_dof(1, cnode);

      // loop over source dofs
      for (int jdof = 0; jdof < nsource_dof; ++jdof)
      {
        int id = discret2()->dof(1, cnode, jdof);
        double val = cnode->x()[jdof];

        var_B->replace_global_value(id, val);
      }
    }
  }

  p21_->multiply(false, *var_A, *result_A);

  // subtract both results
  result_A->update(-1.0, *var_B, 1.0);

  std::cout << "Result of init check= " << std::endl;
  result_A->print(std::cout);
}

/*----------------------------------------------------------------------*
 |  check initial residuum                                   farah 04/14|
 *----------------------------------------------------------------------*/
void Coupling::VolMortar::VolMortarCoupl::mesh_init()
{
  // create merged map:
  std::shared_ptr<Core::DOFSets::DofSetInterface> dofsetaux;
  dofsetaux = std::make_shared<Core::DOFSets::DofSetPredefinedDoFNumber>(dim_, 0, 0, true);
  int dofseta = dis1_->add_dof_set(dofsetaux);
  dofsetaux = std::make_shared<Core::DOFSets::DofSetPredefinedDoFNumber>(dim_, 0, 0, true);
  int dofsetb = dis2_->add_dof_set(dofsetaux);
  dis1_->fill_complete({
      .assign_degrees_of_freedom = true,
      .init_elements = false,
      .do_boundary_conditions = false,
  });
  dis2_->fill_complete({
      .assign_degrees_of_freedom = true,
      .init_elements = false,
      .do_boundary_conditions = false,
  });

  double omega = 1.0;
  double nu = 0.0;
  int maxiter = 1;

  // init zero residuum vector for old iteration
  std::shared_ptr<Core::LinAlg::Vector<double>> ResoldA =
      std::make_shared<Core::LinAlg::Vector<double>>(*discret1()->dof_row_map(dofseta), true);
  std::shared_ptr<Core::LinAlg::Vector<double>> ResoldB =
      std::make_shared<Core::LinAlg::Vector<double>>(*discret2()->dof_row_map(dofsetb), true);

  // output
  if (myrank_ == 0)
  {
    std::cout << "*****       MESHINIT       *****" << std::endl;
  }

  // global mesh init loop
  for (int mi = 0; mi < maxiter; mi++)
  {
    // init
    dmatrix_xa_ =
        std::make_shared<Core::LinAlg::SparseMatrix>(*discret1()->dof_row_map(dofseta), 10);
    mmatrix_xa_ =
        std::make_shared<Core::LinAlg::SparseMatrix>(*discret1()->dof_row_map(dofseta), 100);

    dmatrix_xb_ =
        std::make_shared<Core::LinAlg::SparseMatrix>(*discret2()->dof_row_map(dofsetb), 10);
    mmatrix_xb_ =
        std::make_shared<Core::LinAlg::SparseMatrix>(*discret2()->dof_row_map(dofsetb), 100);

    // output
    if (myrank_ == 0) std::cout << "*****       step " << mi << std::endl;

    // init search trees
    std::shared_ptr<Core::Geo::SearchTree> SearchTreeA = init_search(dis1_);
    std::shared_ptr<Core::Geo::SearchTree> SearchTreeB = init_search(dis2_);

    // calculate DOPs for search algorithm
    std::map<int, Core::LinAlg::Matrix<9, 2>> CurrentDOPsA = calc_background_dops(dis1_);
    std::map<int, Core::LinAlg::Matrix<9, 2>> CurrentDOPsB = calc_background_dops(dis2_);

    /**************************************************
     * loop over all Adis elements                    *
     **************************************************/
    std::set<int> donebeforea;
    for (int j = 0; j < dis1_->num_my_col_elements(); ++j)
    {
      // get target element
      Core::Elements::Element* Aele = dis1_->l_col_element(j);

      std::vector<int> found = search(*Aele, SearchTreeB, CurrentDOPsB);
      integrate_3d_ele_based_a_dis_mesh_init(*Aele, found, dofseta, dofsetb);
    }

    /**************************************************
     * loop over all Bdis elements                    *
     **************************************************/
    std::set<int> donebeforeb;
    for (int j = 0; j < dis2_->num_my_col_elements(); ++j)
    {
      // get target element
      Core::Elements::Element* Bele = dis2_->l_col_element(j);

      std::vector<int> found = search(*Bele, SearchTreeA, CurrentDOPsA);
      integrate_3d_ele_based_b_dis_mesh_init(*Bele, found, dofseta, dofsetb);
    }

    // complete...
    dmatrix_xa_->complete();
    mmatrix_xa_->complete(*discret2()->dof_row_map(dofsetb), *discret1()->dof_row_map(dofseta));
    dmatrix_xb_->complete();
    mmatrix_xb_->complete(*discret1()->dof_row_map(dofseta), *discret2()->dof_row_map(dofsetb));

    mergedmap_ = Core::LinAlg::merge_map(
        *discret1()->dof_row_map(dofseta), *discret2()->dof_row_map(dofsetb), false);
    std::shared_ptr<Core::LinAlg::Vector<double>> mergedsol =
        std::make_shared<Core::LinAlg::Vector<double>>(*mergedmap_);
    std::shared_ptr<Core::LinAlg::Vector<double>> mergedX =
        std::make_shared<Core::LinAlg::Vector<double>>(*mergedmap_);
    std::shared_ptr<Core::LinAlg::Vector<double>> mergedXa =
        std::make_shared<Core::LinAlg::Vector<double>>(*discret1()->dof_row_map(dofseta));
    std::shared_ptr<Core::LinAlg::Vector<double>> mergedXb =
        std::make_shared<Core::LinAlg::Vector<double>>(*discret2()->dof_row_map(dofsetb));

    for (int n = 0; n < dis1_->node_row_map()->num_my_elements(); ++n)
    {
      int gid = dis1_->node_row_map()->gid(n);
      Core::Nodes::Node* node = dis1_->g_node(gid);

      Core::LinAlg::SerialDenseVector pos(3);
      pos(0) = node->x()[0];
      pos(1) = node->x()[1];
      pos(2) = node->x()[2];

      std::vector<int> id(3);
      id[0] = dis1_->dof(dofseta, node, 0);
      id[1] = dis1_->dof(dofseta, node, 1);
      id[2] = dis1_->dof(dofseta, node, 2);

      std::vector<int> owner(3);
      owner[0] = node->owner();
      owner[1] = node->owner();
      owner[2] = node->owner();

      Core::LinAlg::assemble(*mergedX, pos, id, owner);
      Core::LinAlg::assemble(*mergedXa, pos, id, owner);
    }

    for (int n = 0; n < dis2_->node_row_map()->num_my_elements(); ++n)
    {
      int gid = dis2_->node_row_map()->gid(n);
      Core::Nodes::Node* node = dis2_->g_node(gid);

      Core::LinAlg::SerialDenseVector pos(3);
      pos(0) = node->x()[0];
      pos(1) = node->x()[1];
      pos(2) = node->x()[2];

      std::vector<int> id(3);
      id[0] = dis2_->dof(dofsetb, node, 0);
      id[1] = dis2_->dof(dofsetb, node, 1);
      id[2] = dis2_->dof(dofsetb, node, 2);

      std::vector<int> owner(3);
      owner[0] = node->owner();
      owner[1] = node->owner();
      owner[2] = node->owner();

      Core::LinAlg::assemble(*mergedX, pos, id, owner);
      Core::LinAlg::assemble(*mergedXb, pos, id, owner);
    }

    //--------------------------------------------------------------
    //--------------------------------------------------------------
    // Check:
    std::shared_ptr<Core::LinAlg::Vector<double>> solDA =
        std::make_shared<Core::LinAlg::Vector<double>>(*discret1()->dof_row_map(dofseta));
    std::shared_ptr<Core::LinAlg::Vector<double>> solMA =
        std::make_shared<Core::LinAlg::Vector<double>>(*discret1()->dof_row_map(dofseta));

    dmatrix_xa_->multiply(false, *mergedXa, *solDA);
    mmatrix_xa_->multiply(false, *mergedXb, *solMA);

    std::shared_ptr<Core::LinAlg::Vector<double>> solDB =
        std::make_shared<Core::LinAlg::Vector<double>>(*discret2()->dof_row_map(dofsetb));
    std::shared_ptr<Core::LinAlg::Vector<double>> solMB =
        std::make_shared<Core::LinAlg::Vector<double>>(*discret2()->dof_row_map(dofsetb));

    dmatrix_xb_->multiply(false, *mergedXb, *solDB);
    mmatrix_xb_->multiply(false, *mergedXa, *solMB);

    solDA->update(-1.0, *solMA, 1.0);
    solDB->update(-1.0, *solMB, 1.0);

    double ra = 0.0;
    double rb = 0.0;

    // Residuum k+1
    solDA->norm_2(&ra);
    solDB->norm_2(&rb);

    std::cout << "ra= " << ra << "    rb= " << rb << std::endl;

    if (ra < 1e-10 and rb < 1e-10) return;

    solDA->scale(-1.0 * omega);
    solDB->scale(-1.0 * omega);

    // Calc relaxation parameter nu:
    if (mi > 0)
    {
      double fac = 0.0;
      std::shared_ptr<Core::LinAlg::Vector<double>> DiffA =
          std::make_shared<Core::LinAlg::Vector<double>>(*discret1()->dof_row_map(dofseta), true);
      std::shared_ptr<Core::LinAlg::Vector<double>> DiffB =
          std::make_shared<Core::LinAlg::Vector<double>>(*discret2()->dof_row_map(dofsetb), true);
      DiffA->update(1.0, *solDA, 0.0);
      DiffA->update(-1.0, *ResoldA, 1.0);
      DiffB->update(1.0, *solDB, 0.0);
      DiffB->update(-1.0, *ResoldB, 1.0);

      double topa = 0.0;
      DiffA->dot(*solDA, &topa);

      double topb = 0.0;
      DiffB->dot(*solDB, &topb);

      double top = -(topa + topb);

      double downa = 0.0;
      DiffA->dot(*DiffA, &downa);

      double downb = 0.0;
      DiffB->dot(*DiffB, &downb);
      double down = (downa + downb);

      fac = top / down;
      std::cout << "fac= " << fac << std::endl;
      nu = nu + (nu - 1.0) * fac;
    }

    ResoldA->update(1.0, *solDA, 0.0);
    ResoldB->update(1.0, *solDB, 0.0);
    //--------------------------------------------------------------
    //--------------------------------------------------------------
    // relaxation parameter
    //    omega = 1.0 - nu;
    omega = 1.0;

    Core::LinAlg::SparseMatrix k(*mergedmap_, 100, false, true);
    Core::LinAlg::matrix_add(*dmatrix_xa_, false, -1.0 * omega, k, 1.0);
    Core::LinAlg::matrix_add(*mmatrix_xa_, false, 1.0 * omega, k, 1.0);
    Core::LinAlg::matrix_add(*dmatrix_xb_, false, -1.0 * omega, k, 1.0);
    Core::LinAlg::matrix_add(*mmatrix_xb_, false, 1.0 * omega, k, 1.0);

    Core::LinAlg::Vector<double> ones(*mergedmap_);
    ones.put_scalar(1.0);
    Core::LinAlg::SparseMatrix onesdiag(ones);
    onesdiag.complete();
    Core::LinAlg::matrix_add(onesdiag, false, 1.0, k, 1.0);
    k.complete();

    // solve with default solver

    Teuchos::ParameterList solvparams;
    Core::Utils::add_enum_class_to_parameter_list<Core::LinearSolver::SolverType>(
        "SOLVER", Core::LinearSolver::SolverType::UMFPACK, solvparams);
    Core::LinAlg::Solver solver(solvparams, comm_, nullptr, Core::IO::Verbositylevel::standard);

    Core::LinAlg::SolverParams solver_params;
    solver_params.refactor = true;
    solver.solve(Core::Utils::shared_ptr_from_ref(k), mergedsol, mergedX, solver_params);

    std::shared_ptr<Core::LinAlg::Vector<double>> sola =
        std::make_shared<Core::LinAlg::Vector<double>>(*discret1()->dof_row_map(dofseta));
    std::shared_ptr<Core::LinAlg::Vector<double>> solb =
        std::make_shared<Core::LinAlg::Vector<double>>(*discret2()->dof_row_map(dofsetb));

    Core::LinAlg::MapExtractor mapext(*mergedmap_,
        std::make_shared<Core::LinAlg::Map>(*(discret1()->dof_row_map(dofseta))),
        std::make_shared<Core::LinAlg::Map>(*(discret2()->dof_row_map(dofsetb))));
    mapext.extract_cond_vector(*mergedsol, *sola);
    mapext.extract_other_vector(*mergedsol, *solb);

    //--------------------------------------------------------------
    // Post Processing... Node Relocation
    // node positions
    for (int i = 0; i < discret1()->num_my_row_elements(); ++i)
    {
      Core::Elements::Element* Aele = discret1()->l_row_element(i);
      for (int j = 0; j < Aele->num_node(); ++j)
      {
        Core::Nodes::Node* cnode = Aele->nodes()[j];
        int nsource_dof = discret1()->num_dof(dofseta, cnode);
        std::vector<double> nvector(3);

        // loop over source dofs
        for (int jdof = 0; jdof < nsource_dof; ++jdof)
        {
          const int lid = sola->get_map().lid(discret1()->dof(dofseta, cnode, jdof));
          nvector[jdof] = sola->local_values_as_span()[lid] - cnode->x()[jdof];
        }
        cnode->change_pos(nvector);
      }
    }

    for (int i = 0; i < discret2()->num_my_row_elements(); ++i)
    {
      Core::Elements::Element* Bele = discret2()->l_row_element(i);
      for (int j = 0; j < Bele->num_node(); ++j)
      {
        Core::Nodes::Node* cnode = Bele->nodes()[j];
        int nsource_dof = discret2()->num_dof(dofsetb, cnode);
        std::vector<double> nvector(3);

        // loop over source dofs
        for (int jdof = 0; jdof < nsource_dof; ++jdof)
        {
          const int lid = solb->get_map().lid(discret2()->dof(dofsetb, cnode, jdof));
          nvector[jdof] = solb->local_values_as_span()[lid] - cnode->x()[jdof];
        }
        cnode->change_pos(nvector);
      }
    }

    dis1_->fill_complete({
        .assign_degrees_of_freedom = false,
        .init_elements = true,
        .do_boundary_conditions = true,
    });
    dis2_->fill_complete({
        .assign_degrees_of_freedom = false,
        .init_elements = true,
        .do_boundary_conditions = true,
    });

    // last check:
    std::shared_ptr<Core::LinAlg::Vector<double>> checka =
        std::make_shared<Core::LinAlg::Vector<double>>(*discret1()->dof_row_map(dofseta));
    std::shared_ptr<Core::LinAlg::Vector<double>> checkb =
        std::make_shared<Core::LinAlg::Vector<double>>(*discret2()->dof_row_map(dofsetb));

    for (int n = 0; n < dis1_->node_row_map()->num_my_elements(); ++n)
    {
      int gid = dis1_->node_row_map()->gid(n);
      Core::Nodes::Node* node = dis1_->g_node(gid);

      Core::LinAlg::SerialDenseVector pos(3);
      pos(0) = node->x()[0];
      pos(1) = node->x()[1];
      pos(2) = node->x()[2];

      std::vector<int> id(3);
      id[0] = dis1_->dof(dofseta, node, 0);
      id[1] = dis1_->dof(dofseta, node, 1);
      id[2] = dis1_->dof(dofseta, node, 2);

      std::vector<int> owner(3);
      owner[0] = node->owner();
      owner[1] = node->owner();
      owner[2] = node->owner();

      Core::LinAlg::assemble(*checka, pos, id, owner);
    }

    for (int n = 0; n < dis2_->node_row_map()->num_my_elements(); ++n)
    {
      int gid = dis2_->node_row_map()->gid(n);
      Core::Nodes::Node* node = dis2_->g_node(gid);

      Core::LinAlg::SerialDenseVector pos(3);
      pos(0) = node->x()[0];
      pos(1) = node->x()[1];
      pos(2) = node->x()[2];

      std::vector<int> id(3);
      id[0] = dis2_->dof(dofsetb, node, 0);
      id[1] = dis2_->dof(dofsetb, node, 1);
      id[2] = dis2_->dof(dofsetb, node, 2);

      std::vector<int> owner(3);
      owner[0] = node->owner();
      owner[1] = node->owner();
      owner[2] = node->owner();

      Core::LinAlg::assemble(*checkb, pos, id, owner);
    }

    //--------------------------------------------------------------
    //--------------------------------------------------------------
    // Check:
    std::shared_ptr<Core::LinAlg::Vector<double>> finalDA =
        std::make_shared<Core::LinAlg::Vector<double>>(*discret1()->dof_row_map(dofseta));
    std::shared_ptr<Core::LinAlg::Vector<double>> finalMA =
        std::make_shared<Core::LinAlg::Vector<double>>(*discret1()->dof_row_map(dofseta));

    dmatrix_xa_->multiply(false, *checka, *finalDA);
    mmatrix_xa_->multiply(false, *checkb, *finalMA);

    std::shared_ptr<Core::LinAlg::Vector<double>> finalDB =
        std::make_shared<Core::LinAlg::Vector<double>>(*discret2()->dof_row_map(dofsetb));
    std::shared_ptr<Core::LinAlg::Vector<double>> finalMB =
        std::make_shared<Core::LinAlg::Vector<double>>(*discret2()->dof_row_map(dofsetb));

    dmatrix_xb_->multiply(false, *checkb, *finalDB);
    mmatrix_xb_->multiply(false, *checka, *finalMB);

    finalDA->update(-1.0, *finalMA, 1.0);

    finalDB->update(-1.0, *finalMB, 1.0);

    double finalra = 0.0;
    double finalrb = 0.0;

    // Residuum k+1
    finalDA->norm_2(&finalra);
    finalDB->norm_2(&finalrb);

    std::cout << "final ra= " << finalra << "   final rb= " << finalrb << std::endl;
  }
}
/*----------------------------------------------------------------------*
 |  print_status (public)                                     farah 02/14|
 *----------------------------------------------------------------------*/
void Coupling::VolMortar::VolMortarCoupl::print_status(int& i, bool dis_switch)
{
  static int percent_counter = 0;
  static int EleSum = 0;

  if (i == 0)
  {
    if (dis_switch)
      EleSum = dis2_->num_global_elements();
    else
      EleSum = dis1_->num_global_elements();

    percent_counter = 0;
  }

  if ((int)((i * 100) / EleSum) > (int)(10 * percent_counter))
  {
    std::cout << "---------------------------" << std::endl;
    std::cout << (int)((i * 100) / EleSum) - 1 << "% of Coupling Evaluations are done!"
              << std::endl;
    std::cout << "---------------------------" << std::endl;

    percent_counter++;
  }
}

/*----------------------------------------------------------------------*
 |  Start Cut routine                                        farah 01/14|
 *----------------------------------------------------------------------*/
void Coupling::VolMortar::VolMortarCoupl::perform_cut(
    Core::Elements::Element* source_ele, Core::Elements::Element* target_ele, bool switched_conf)
{
  // create empty vector of integration cells
  std::vector<std::shared_ptr<Cell>> IntCells;

  // the cut wizard wants discretizations to perform the cut. One is supposed to be the background
  // discretization and the other the interface discretization. As we want to cut only two 3D
  // elements, we build two auxiliary discretizations. The first holds only a copy of the target
  // element, acting as background mesh, and the other one is build of the surface elements of the
  // source element, being the interface discretization. We use temporary copies of all elements and
  // nodes, as we only need the geometry to perform the cut, but want to make sure that the gids
  // and dofs of the original elements are kept untouched.

  std::shared_ptr<Core::FE::Discretization> sourceauxdis =
      std::make_shared<Core::FE::Discretization>((std::string) "sourceauxdis", comm_, dim_);
  std::shared_ptr<Core::FE::Discretization> targetauxdis =
      std::make_shared<Core::FE::Discretization>((std::string) "targetauxdis", comm_, dim_);

  // build surface elements for all surfaces of source element
  std::vector<std::shared_ptr<Core::Elements::Element>> sele_surfs = source_ele->surfaces();
  const int numsurf = sele_surfs.size();
  for (int isurf = 0; isurf < numsurf; ++isurf)
  {
    // add all surface elements to auxiliary source discretization (no need for cloning, as surface
    // elements are
    // rebuild every time when calling Surfaces(), anyway)
    sourceauxdis->add_element(sele_surfs[isurf]);
  }

  // add clone of element to auxiliary discretization
  targetauxdis->add_element(std::shared_ptr<Core::Elements::Element>(target_ele->clone()));

  // add clones of nodes to auxiliary discretizations
  for (int node = 0; node < source_ele->num_node(); ++node)
  {
    auto cloned_node = std::shared_ptr<Core::Nodes::Node>(source_ele->nodes()[node]->clone());
    sourceauxdis->add_node(cloned_node->x(), cloned_node->id(), cloned_node);
  }
  for (int node = 0; node < target_ele->num_node(); ++node)
  {
    auto cloned_node = std::shared_ptr<Core::Nodes::Node>(target_ele->nodes()[node]->clone());
    targetauxdis->add_node(cloned_node->x(), cloned_node->id(), cloned_node);
  }

  // complete dis
  sourceauxdis->fill_complete({
      .assign_degrees_of_freedom = true,
      .init_elements = false,
      .do_boundary_conditions = false,
  });
  targetauxdis->fill_complete({
      .assign_degrees_of_freedom = true,
      .init_elements = false,
      .do_boundary_conditions = false,
  });

  //--------------------------------------------------------------------------------------
  // Initialize the cut wizard

  // create new cut wizard
  Cut::CutWizard wizard(targetauxdis);

  // *************************************
  // TESSELATION *************************
  if (Teuchos::getIntegralValue<CutType>(params(), "CUTTYPE") == cuttype_tessellation)
  {
    // Set options for the cut wizard
    wizard.set_options(cut_params_, Cut::NDS_Strategy_full,
        Cut::VCellGaussPts_Tessellation,  // how to create volume cell Gauss points?
        Cut::BCellGaussPts_Tessellation,  // how to create boundary cell Gauss points?
        "invalid_file",                   // no file needed
        false,                            // gmsh output for cut library
        true,                             // find point positions
        true,                             // generate only tet cells
        false                             // print screen output
    );

    // cut in reference configuration
    wizard.set_background_state(nullptr, nullptr, -1);
    wizard.add_cutter_state(0, sourceauxdis, nullptr);

    wizard.prepare();
    wizard.cut(true);  // include_inner

    Cut::plain_volumecell_set target_cells_out;
    Cut::plain_volumecell_set target_cells_in;
    Cut::ElementHandle* em = wizard.get_element(target_ele);

    // is target_ele in cut involved?
    if (em != nullptr)
    {
      em->collect_volume_cells(target_cells_in, target_cells_out);

      int count = 0;

      for (Cut::plain_volumecell_set::iterator u = target_cells_in.begin();
          u != target_cells_in.end(); u++)
      {
        Cut::VolumeCell* vc = *u;
        const Cut::plain_integrationcell_set& intcells = vc->integration_cells();

        for (Cut::plain_integrationcell_set::const_iterator z = intcells.begin();
            z != intcells.end(); z++)
        {
          Cut::IntegrationCell* ic = *z;

          IntCells.push_back(std::make_shared<Coupling::VolMortar::Cell>(
              count, 4, ic->coordinates(), ic->shape()));
          volume_ += IntCells[count]->vol();

          count++;
        }
      }
      // integrate found cells - tesselation
      if (!switched_conf)
        integrate_3d_cell(*source_ele, *target_ele, IntCells);
      else
        integrate_3d_cell(*target_ele, *source_ele, IntCells);

      // count the cells and the polygons/polyhedra
      polygoncounter_ += target_cells_in.size();
      cellcounter_ += count;
    }
  }

  // *******************************************
  // DIRECT DIVERGENCE *************************
  else if (Teuchos::getIntegralValue<CutType>(params(), "CUTTYPE") == cuttype_directdivergence)
  {
    // Set options for the cut wizard
    wizard.set_options(cut_params_, Cut::NDS_Strategy_full,
        Cut::VCellGaussPts_DirectDivergence,  // how to create volume cell Gauss points?
        Cut::BCellGaussPts_Tessellation,      // how to create boundary cell Gauss points?
        "invalid_file",                       // no file needed
        false,                                // gmsh output for cut library
        true,                                 // find point positions
        false,                                // generate only tet cells
        false                                 // print screen output
    );

    // cut in reference configuration
    wizard.set_background_state(nullptr, nullptr, -1);
    wizard.add_cutter_state(0, sourceauxdis, nullptr);

    wizard.cut(true);  // include_inner

    Cut::plain_volumecell_set target_cells_out;
    Cut::plain_volumecell_set target_cells_in;
    Cut::ElementHandle* em = wizard.get_element(target_ele);

    // for safety
    volcell_.clear();
    if (em != nullptr)
    {
      em->collect_volume_cells(volcell_, target_cells_out);

      // start integration
      if (switched_conf)
        integrate_3d_cell_direct_divergence(*target_ele, *source_ele, switched_conf);
      else
        integrate_3d_cell_direct_divergence(*source_ele, *target_ele, switched_conf);

      polygoncounter_ += volcell_.size();
    }
  }

  // *******************************************
  // DEFAULT           *************************
  else
    FOUR_C_THROW("ERROR: Chosen Cuttype for volmortar not supported!");
}

/*----------------------------------------------------------------------*
 |  Check need for element-based integration                 farah 01/14|
 *----------------------------------------------------------------------*/
bool Coupling::VolMortar::VolMortarCoupl::check_ele_integration(
    Core::Elements::Element& source_ele, Core::Elements::Element& target_ele)
{
  bool integrateele = true;
  bool converged = false;

  double xi[3] = {0.0, 0.0, 0.0};
  double xgl[3] = {0.0, 0.0, 0.0};

  //--------------------------------------------------------
  // 1. all source nodes within with target ele ?
  for (int u = 0; u < source_ele.num_node(); ++u)
  {
    xgl[0] = source_ele.nodes()[u]->x()[0];
    xgl[1] = source_ele.nodes()[u]->x()[1];
    xgl[2] = source_ele.nodes()[u]->x()[2];

    // global to local:
    if (target_ele.shape() == Core::FE::CellType::hex8)
      Mortar::Utils::global_to_local<Core::FE::CellType::hex8>(target_ele, xgl, xi, converged);
    else if (target_ele.shape() == Core::FE::CellType::hex20)
      Mortar::Utils::global_to_local<Core::FE::CellType::hex20>(target_ele, xgl, xi, converged);
    else if (target_ele.shape() == Core::FE::CellType::hex27)
      Mortar::Utils::global_to_local<Core::FE::CellType::hex27>(target_ele, xgl, xi, converged);
    else if (target_ele.shape() == Core::FE::CellType::tet4)
      Mortar::Utils::global_to_local<Core::FE::CellType::tet4>(target_ele, xgl, xi, converged);
    else if (target_ele.shape() == Core::FE::CellType::tet10)
      Mortar::Utils::global_to_local<Core::FE::CellType::tet10>(target_ele, xgl, xi, converged);
    else if (target_ele.shape() == Core::FE::CellType::pyramid5)
      Mortar::Utils::global_to_local<Core::FE::CellType::pyramid5>(target_ele, xgl, xi, converged);
    else
      FOUR_C_THROW("ERROR: Shape function not supported!");

    if (converged == true)
    {
      if (target_ele.shape() == Core::FE::CellType::hex8 or
          target_ele.shape() == Core::FE::CellType::hex20 or
          target_ele.shape() == Core::FE::CellType::hex27)
      {
        if (xi[0] > -1.0 - VOLMORTARELETOL and xi[0] < 1.0 + VOLMORTARELETOL and
            xi[1] > -1.0 - VOLMORTARELETOL and xi[1] < 1.0 + VOLMORTARELETOL and
            xi[2] > -1.0 - VOLMORTARELETOL and xi[2] < 1.0 + VOLMORTARELETOL)
          integrateele = true;
        else
          return false;
      }
      else if (target_ele.shape() == Core::FE::CellType::tet4 or
               target_ele.shape() == Core::FE::CellType::tet10)
      {
        if (xi[0] > 0.0 - VOLMORTARELETOL and xi[0] < 1.0 + VOLMORTARELETOL and
            xi[1] > 0.0 - VOLMORTARELETOL and xi[1] < 1.0 + VOLMORTARELETOL and
            xi[2] > 0.0 - VOLMORTARELETOL and xi[2] < 1.0 + VOLMORTARELETOL and
            (xi[0] + xi[1] + xi[2]) < 1.0 + 3 * VOLMORTARELETOL)
          integrateele = true;
        else
          return false;
      }
      else
        FOUR_C_THROW("ERROR: Element type not supported!");
    }
    else
    {
      //      std::cout << "!!! GLOBAL TO LOCAL NOT CONVERGED !!!" << std::endl;
      return false;
    }
  }

  return integrateele;
}

/*----------------------------------------------------------------------*
 |  Check need for cut and element-based integration         farah 01/14|
 *----------------------------------------------------------------------*/
bool Coupling::VolMortar::VolMortarCoupl::check_cut(
    Core::Elements::Element& source_ele, Core::Elements::Element& target_ele)
{
  double xi[3] = {0.0, 0.0, 0.0};
  double xgl[3] = {0.0, 0.0, 0.0};
  bool converged = false;

  {
    //--------------------------------------------------------
    // 1. check if all nodes are outside a parameter space surface
    // bools for tet
    bool xi0 = false;
    bool xi1 = false;
    bool xi2 = false;
    bool all = false;

    // bools for hex
    bool xi0n = false;
    bool xi1n = false;
    bool xi2n = false;

    for (int u = 0; u < target_ele.num_node(); ++u)
    {
      xgl[0] = target_ele.nodes()[u]->x()[0];
      xgl[1] = target_ele.nodes()[u]->x()[1];
      xgl[2] = target_ele.nodes()[u]->x()[2];

      // global to local:
      if (source_ele.shape() == Core::FE::CellType::hex8)
        Mortar::Utils::global_to_local<Core::FE::CellType::hex8>(source_ele, xgl, xi, converged);
      else if (source_ele.shape() == Core::FE::CellType::hex20)
        Mortar::Utils::global_to_local<Core::FE::CellType::hex20>(source_ele, xgl, xi, converged);
      else if (source_ele.shape() == Core::FE::CellType::hex27)
        Mortar::Utils::global_to_local<Core::FE::CellType::hex27>(source_ele, xgl, xi, converged);
      else if (source_ele.shape() == Core::FE::CellType::tet4)
        Mortar::Utils::global_to_local<Core::FE::CellType::tet4>(source_ele, xgl, xi, converged);
      else if (source_ele.shape() == Core::FE::CellType::tet10)
        Mortar::Utils::global_to_local<Core::FE::CellType::tet10>(source_ele, xgl, xi, converged);
      else if (source_ele.shape() == Core::FE::CellType::pyramid5)
        Mortar::Utils::global_to_local<Core::FE::CellType::pyramid5>(
            source_ele, xgl, xi, converged);
      else
        FOUR_C_THROW("ERROR: Shape function not supported!");

      if (converged == true)
      {
        if (source_ele.shape() == Core::FE::CellType::hex8 or
            source_ele.shape() == Core::FE::CellType::hex20 or
            source_ele.shape() == Core::FE::CellType::hex27)
        {
          if (xi[0] > -1.0 + VOLMORTARCUTTOL) xi0 = true;
          if (xi[1] > -1.0 + VOLMORTARCUTTOL) xi1 = true;
          if (xi[2] > -1.0 + VOLMORTARCUTTOL) xi2 = true;

          if (xi[0] < 1.0 - VOLMORTARCUTTOL) xi0n = true;
          if (xi[1] < 1.0 - VOLMORTARCUTTOL) xi1n = true;
          if (xi[2] < 1.0 - VOLMORTARCUTTOL) xi2n = true;
        }
        else if (source_ele.shape() == Core::FE::CellType::tet4 or
                 source_ele.shape() == Core::FE::CellType::tet10)
        {
          if (xi[0] > 0.0 + VOLMORTARCUTTOL) xi0 = true;
          if (xi[1] > 0.0 + VOLMORTARCUTTOL) xi1 = true;
          if (xi[2] > 0.0 + VOLMORTARCUTTOL) xi2 = true;
          if ((xi[0] + xi[1] + xi[2]) < 1.0 - 3.0 * VOLMORTARCUTTOL) all = true;
        }
        else
          FOUR_C_THROW("ERROR: Element not supported!");
      }
    }  // end node loop

    if (source_ele.shape() == Core::FE::CellType::tet4 or
        source_ele.shape() == Core::FE::CellType::tet10)
    {
      if (!xi0 or !xi1 or !xi2 or !all) return false;
    }
    else if (source_ele.shape() == Core::FE::CellType::hex8 or
             source_ele.shape() == Core::FE::CellType::hex20 or
             source_ele.shape() == Core::FE::CellType::hex27)
    {
      if (!xi0 or !xi1 or !xi2 or !xi0n or !xi1n or !xi2n) return false;
    }
    else
      FOUR_C_THROW("ERROR: Element not supported!");
  }

  {
    //--------------------------------------------------------
    // 2. check if all nodes are outside a parameter space surface (changed elements)
    bool xi0 = false;
    bool xi1 = false;
    bool xi2 = false;
    bool all = false;

    // bools for hex
    bool xi0n = false;
    bool xi1n = false;
    bool xi2n = false;

    for (int u = 0; u < source_ele.num_node(); ++u)
    {
      xgl[0] = source_ele.nodes()[u]->x()[0];
      xgl[1] = source_ele.nodes()[u]->x()[1];
      xgl[2] = source_ele.nodes()[u]->x()[2];

      // global to local:
      if (target_ele.shape() == Core::FE::CellType::hex8)
        Mortar::Utils::global_to_local<Core::FE::CellType::hex8>(target_ele, xgl, xi, converged);
      else if (target_ele.shape() == Core::FE::CellType::hex20)
        Mortar::Utils::global_to_local<Core::FE::CellType::hex20>(target_ele, xgl, xi, converged);
      else if (target_ele.shape() == Core::FE::CellType::hex27)
        Mortar::Utils::global_to_local<Core::FE::CellType::hex27>(target_ele, xgl, xi, converged);
      else if (target_ele.shape() == Core::FE::CellType::tet4)
        Mortar::Utils::global_to_local<Core::FE::CellType::tet4>(target_ele, xgl, xi, converged);
      else if (target_ele.shape() == Core::FE::CellType::tet10)
        Mortar::Utils::global_to_local<Core::FE::CellType::tet10>(target_ele, xgl, xi, converged);
      else if (target_ele.shape() == Core::FE::CellType::pyramid5)
        Mortar::Utils::global_to_local<Core::FE::CellType::pyramid5>(
            target_ele, xgl, xi, converged);
      else
        FOUR_C_THROW("ERROR: Shape function not supported!");

      if (converged == true)
      {
        if (target_ele.shape() == Core::FE::CellType::hex8 or
            target_ele.shape() == Core::FE::CellType::hex20 or
            target_ele.shape() == Core::FE::CellType::hex27)
        {
          if (xi[0] > -1.0 + VOLMORTARCUTTOL) xi0 = true;
          if (xi[1] > -1.0 + VOLMORTARCUTTOL) xi1 = true;
          if (xi[2] > -1.0 + VOLMORTARCUTTOL) xi2 = true;

          if (xi[0] < 1.0 - VOLMORTARCUTTOL) xi0n = true;
          if (xi[1] < 1.0 - VOLMORTARCUTTOL) xi1n = true;
          if (xi[2] < 1.0 - VOLMORTARCUTTOL) xi2n = true;
        }

        else if (target_ele.shape() == Core::FE::CellType::tet4 or
                 target_ele.shape() == Core::FE::CellType::tet10)
        {
          if (xi[0] > 0.0 + VOLMORTARCUTTOL) xi0 = true;
          if (xi[1] > 0.0 + VOLMORTARCUTTOL) xi1 = true;
          if (xi[2] > 0.0 + VOLMORTARCUTTOL) xi2 = true;
          if ((xi[0] + xi[1] + xi[2]) < 1.0 - 3.0 * VOLMORTARCUTTOL) all = true;
        }
        else
          FOUR_C_THROW("ERROR: Element not supported!");
      }
    }

    if (target_ele.shape() == Core::FE::CellType::tet4 or
        target_ele.shape() == Core::FE::CellType::tet10)
    {
      if (!xi0 or !xi1 or !xi2 or !all) return false;
    }
    else if (target_ele.shape() == Core::FE::CellType::hex8 or
             target_ele.shape() == Core::FE::CellType::hex20 or
             target_ele.shape() == Core::FE::CellType::hex27)
    {
      if (!xi0 or !xi1 or !xi2 or !xi0n or !xi1n or !xi2n) return false;
    }
    else
      FOUR_C_THROW("ERROR: Element not supported!");
  }

  //--------------------------------------------------------
  // 3. target nodes within source parameter space?
  for (int u = 0; u < target_ele.num_node(); ++u)
  {
    xgl[0] = target_ele.nodes()[u]->x()[0];
    xgl[1] = target_ele.nodes()[u]->x()[1];
    xgl[2] = target_ele.nodes()[u]->x()[2];

    // global to local:
    if (source_ele.shape() == Core::FE::CellType::hex8)
      Mortar::Utils::global_to_local<Core::FE::CellType::hex8>(source_ele, xgl, xi, converged);
    else if (source_ele.shape() == Core::FE::CellType::hex20)
      Mortar::Utils::global_to_local<Core::FE::CellType::hex20>(source_ele, xgl, xi, converged);
    else if (source_ele.shape() == Core::FE::CellType::hex27)
      Mortar::Utils::global_to_local<Core::FE::CellType::hex27>(source_ele, xgl, xi, converged);
    else if (source_ele.shape() == Core::FE::CellType::tet4)
      Mortar::Utils::global_to_local<Core::FE::CellType::tet4>(source_ele, xgl, xi, converged);
    else if (source_ele.shape() == Core::FE::CellType::tet10)
      Mortar::Utils::global_to_local<Core::FE::CellType::tet10>(source_ele, xgl, xi, converged);
    else if (source_ele.shape() == Core::FE::CellType::pyramid5)
      Mortar::Utils::global_to_local<Core::FE::CellType::pyramid5>(source_ele, xgl, xi, converged);
    else
      FOUR_C_THROW("ERROR: Shape function not supported!");

    if (converged == true)
    {
      if (source_ele.shape() == Core::FE::CellType::hex8 or
          source_ele.shape() == Core::FE::CellType::hex20 or
          source_ele.shape() == Core::FE::CellType::hex27)
      {
        if (abs(xi[0]) < 1.0 - VOLMORTARCUT2TOL and abs(xi[1]) < 1.0 - VOLMORTARCUT2TOL and
            abs(xi[2]) < 1.0 - VOLMORTARCUT2TOL)
          return true;
      }
      else if (source_ele.shape() == Core::FE::CellType::tet4 or
               source_ele.shape() == Core::FE::CellType::tet10)
      {
        if (xi[0] > 0.0 + VOLMORTARCUT2TOL and xi[0] < 1.0 - VOLMORTARCUT2TOL and
            xi[1] > 0.0 + VOLMORTARCUT2TOL and xi[1] < 1.0 - VOLMORTARCUT2TOL and
            xi[2] > 0.0 + VOLMORTARCUT2TOL and xi[2] < 1.0 - VOLMORTARCUT2TOL and
            (xi[0] + xi[1] + xi[2]) < 1.0 - 3.0 * VOLMORTARCUT2TOL)
          return true;
      }
      else
        FOUR_C_THROW("ERROR: Element not supported!");
    }
  }

  //--------------------------------------------------------
  // 4. source nodes within target parameter space?
  for (int u = 0; u < source_ele.num_node(); ++u)
  {
    xgl[0] = source_ele.nodes()[u]->x()[0];
    xgl[1] = source_ele.nodes()[u]->x()[1];
    xgl[2] = source_ele.nodes()[u]->x()[2];

    // global to local:
    if (target_ele.shape() == Core::FE::CellType::hex8)
      Mortar::Utils::global_to_local<Core::FE::CellType::hex8>(target_ele, xgl, xi, converged);
    else if (target_ele.shape() == Core::FE::CellType::hex20)
      Mortar::Utils::global_to_local<Core::FE::CellType::hex20>(target_ele, xgl, xi, converged);
    else if (target_ele.shape() == Core::FE::CellType::hex27)
      Mortar::Utils::global_to_local<Core::FE::CellType::hex27>(target_ele, xgl, xi, converged);
    else if (target_ele.shape() == Core::FE::CellType::tet4)
      Mortar::Utils::global_to_local<Core::FE::CellType::tet4>(target_ele, xgl, xi, converged);
    else if (target_ele.shape() == Core::FE::CellType::tet10)
      Mortar::Utils::global_to_local<Core::FE::CellType::tet10>(target_ele, xgl, xi, converged);
    else if (target_ele.shape() == Core::FE::CellType::pyramid5)
      Mortar::Utils::global_to_local<Core::FE::CellType::pyramid5>(target_ele, xgl, xi, converged);
    else
      FOUR_C_THROW("ERROR: Shape function not supported!");

    if (converged == true)
    {
      if (target_ele.shape() == Core::FE::CellType::hex8 or
          target_ele.shape() == Core::FE::CellType::hex20 or
          target_ele.shape() == Core::FE::CellType::hex27)
      {
        if (abs(xi[0]) < 1.0 - VOLMORTARCUT2TOL and abs(xi[1]) < 1.0 - VOLMORTARCUT2TOL and
            abs(xi[2]) < 1.0 - VOLMORTARCUT2TOL)
          return true;
      }
      else if (target_ele.shape() == Core::FE::CellType::tet4 or
               target_ele.shape() == Core::FE::CellType::tet10)
      {
        if (xi[0] > 0.0 + VOLMORTARCUT2TOL and xi[0] < 1.0 - VOLMORTARCUT2TOL and
            xi[1] > 0.0 + VOLMORTARCUT2TOL and xi[1] < 1.0 - VOLMORTARCUT2TOL and
            xi[2] > 0.0 + VOLMORTARCUT2TOL and xi[2] < 1.0 - VOLMORTARCUT2TOL and
            (xi[0] + xi[1] + xi[2]) < 1.0 - 3.0 * VOLMORTARCUT2TOL)
          return true;
      }
      else
        FOUR_C_THROW("ERROR: Element not supported!");
    }
  }

  return false;
}

/*----------------------------------------------------------------------*
 |  integrate_2d Cells                                        farah 01/14|
 *----------------------------------------------------------------------*/
void Coupling::VolMortar::VolMortarCoupl::integrate_2d(Core::Elements::Element& source_ele,
    Core::Elements::Element& target_ele, std::vector<std::shared_ptr<Mortar::IntCell>>& cells)
{
  //--------------------------------------------------------------------
  // loop over cells for A Field
  // loop over cells
  for (int q = 0; q < (int)cells.size(); ++q)
  {
    switch (source_ele.shape())
    {
      // 2D surface elements
      case Core::FE::CellType::quad4:
      {
        switch (target_ele.shape())
        {
          // 2D surface elements
          case Core::FE::CellType::quad4:
          {
            static VolMortarIntegrator<Core::FE::CellType::quad4, Core::FE::CellType::quad4>
                integrator(params());
            integrator.integrate_cells_2d(source_ele, target_ele, *cells[q], *d1_, *m12_, *dis1_,
                *dis2_, dofset12_.first, dofset12_.second);
            break;
          }
          case Core::FE::CellType::tri3:
          {
            static VolMortarIntegrator<Core::FE::CellType::quad4, Core::FE::CellType::tri3>
                integrator(params());
            integrator.integrate_cells_2d(source_ele, target_ele, *cells[q], *d1_, *m12_, *dis1_,
                *dis2_, dofset12_.first, dofset12_.second);
            break;
          }
          default:
          {
            FOUR_C_THROW("ERROR: unknown shape!");
            break;
          }
        }
        break;
      }
      case Core::FE::CellType::tri3:
      {
        switch (target_ele.shape())
        {
          // 2D surface elements
          case Core::FE::CellType::quad4:
          {
            static VolMortarIntegrator<Core::FE::CellType::tri3, Core::FE::CellType::quad4>
                integrator(params());
            integrator.integrate_cells_2d(source_ele, target_ele, *cells[q], *d1_, *m12_, *dis1_,
                *dis2_, dofset12_.first, dofset12_.second);
            break;
          }
          case Core::FE::CellType::tri3:
          {
            static VolMortarIntegrator<Core::FE::CellType::tri3, Core::FE::CellType::tri3>
                integrator(params());
            integrator.integrate_cells_2d(source_ele, target_ele, *cells[q], *d1_, *m12_, *dis1_,
                *dis2_, dofset12_.first, dofset12_.second);
            break;
          }
          default:
          {
            FOUR_C_THROW("ERROR: unknown shape!");
            break;
          }
        }
        break;
      }
      default:
      {
        FOUR_C_THROW("ERROR: unknown shape!");
        break;
      }
    }

    //--------------------------------------------------------------------
    // loop over cells for A Field
    switch (target_ele.shape())
    {
      // 2D surface elements
      case Core::FE::CellType::quad4:
      {
        switch (source_ele.shape())
        {
          // 2D surface elements
          case Core::FE::CellType::quad4:
          {
            static VolMortarIntegrator<Core::FE::CellType::quad4, Core::FE::CellType::quad4>
                integrator(params());
            integrator.integrate_cells_2d(target_ele, source_ele, *cells[q], *d2_, *m21_, *dis2_,
                *dis1_, dofset21_.first, dofset21_.second);
            break;
          }
          case Core::FE::CellType::tri3:
          {
            static VolMortarIntegrator<Core::FE::CellType::quad4, Core::FE::CellType::tri3>
                integrator(params());
            integrator.integrate_cells_2d(target_ele, source_ele, *cells[q], *d2_, *m21_, *dis2_,
                *dis1_, dofset21_.first, dofset21_.second);
            break;
          }
          default:
          {
            FOUR_C_THROW("ERROR: unknown shape!");
            break;
          }
        }
        break;
      }
      case Core::FE::CellType::tri3:
      {
        switch (source_ele.shape())
        {
          // 2D surface elements
          case Core::FE::CellType::quad4:
          {
            static VolMortarIntegrator<Core::FE::CellType::tri3, Core::FE::CellType::quad4>
                integrator(params());
            integrator.integrate_cells_2d(target_ele, source_ele, *cells[q], *d2_, *m21_, *dis2_,
                *dis1_, dofset21_.first, dofset21_.second);
            break;
          }
          case Core::FE::CellType::tri3:
          {
            static VolMortarIntegrator<Core::FE::CellType::tri3, Core::FE::CellType::tri3>
                integrator(params());
            integrator.integrate_cells_2d(target_ele, source_ele, *cells[q], *d2_, *m21_, *dis2_,
                *dis1_, dofset21_.first, dofset21_.second);
            break;
          }
          default:
          {
            FOUR_C_THROW("ERROR: unknown shape!");
            break;
          }
        }
        break;
      }
      default:
      {
        FOUR_C_THROW("ERROR: unknown shape!");
        break;
      }
    }
  }  // end loop
}

/*----------------------------------------------------------------------*
 |  integrate_3d Cells                                        farah 01/14|
 *----------------------------------------------------------------------*/
void Coupling::VolMortar::VolMortarCoupl::integrate_3d_cell(Core::Elements::Element& source_ele,
    Core::Elements::Element& target_ele, std::vector<std::shared_ptr<Cell>>& cells)
{
  //--------------------------------------------------------------------
  // loop over cells for A Field
  for (int q = 0; q < (int)cells.size(); ++q)
  {
    switch (source_ele.shape())
    {
      // 2D surface elements
      case Core::FE::CellType::hex8:
      {
        switch (target_ele.shape())
        {
          // 2D surface elements
          case Core::FE::CellType::hex8:
          {
            static VolMortarIntegrator<Core::FE::CellType::hex8, Core::FE::CellType::hex8>
                integrator(params());
            integrator.initialize_gp(false, 0, cells[q]->shape());
            integrator.integrate_cells_3d(source_ele, target_ele, *cells[q], *d1_, *m12_, *d2_,
                *m21_, *dis1_, *dis2_, dofset12_.first, dofset12_.second, dofset21_.first,
                dofset21_.second);
            break;
          }
          case Core::FE::CellType::hex20:
          {
            static VolMortarIntegrator<Core::FE::CellType::hex8, Core::FE::CellType::hex20>
                integrator(params());
            integrator.initialize_gp(false, 0, cells[q]->shape());
            integrator.integrate_cells_3d(source_ele, target_ele, *cells[q], *d1_, *m12_, *d2_,
                *m21_, *dis1_, *dis2_, dofset12_.first, dofset12_.second, dofset21_.first,
                dofset21_.second);
            break;
          }
          case Core::FE::CellType::hex27:
          {
            static VolMortarIntegrator<Core::FE::CellType::hex8, Core::FE::CellType::hex27>
                integrator(params());
            integrator.initialize_gp(false, 0, cells[q]->shape());
            integrator.integrate_cells_3d(source_ele, target_ele, *cells[q], *d1_, *m12_, *d2_,
                *m21_, *dis1_, *dis2_, dofset12_.first, dofset12_.second, dofset21_.first,
                dofset21_.second);
            break;
          }
          case Core::FE::CellType::tet4:
          {
            static VolMortarIntegrator<Core::FE::CellType::hex8, Core::FE::CellType::tet4>
                integrator(params());
            integrator.initialize_gp(false, 0, cells[q]->shape());
            integrator.integrate_cells_3d(source_ele, target_ele, *cells[q], *d1_, *m12_, *d2_,
                *m21_, *dis1_, *dis2_, dofset12_.first, dofset12_.second, dofset21_.first,
                dofset21_.second);
            break;
          }
          case Core::FE::CellType::tet10:
          {
            static VolMortarIntegrator<Core::FE::CellType::hex8, Core::FE::CellType::tet10>
                integrator(params());
            integrator.initialize_gp(false, 0, cells[q]->shape());
            integrator.integrate_cells_3d(source_ele, target_ele, *cells[q], *d1_, *m12_, *d2_,
                *m21_, *dis1_, *dis2_, dofset12_.first, dofset12_.second, dofset21_.first,
                dofset21_.second);
            break;
          }
          case Core::FE::CellType::pyramid5:
          {
            static VolMortarIntegrator<Core::FE::CellType::hex8, Core::FE::CellType::pyramid5>
                integrator(params());
            integrator.initialize_gp(false, 0, cells[q]->shape());
            integrator.integrate_cells_3d(source_ele, target_ele, *cells[q], *d1_, *m12_, *d2_,
                *m21_, *dis1_, *dis2_, dofset12_.first, dofset12_.second, dofset21_.first,
                dofset21_.second);
            break;
          }
          default:
          {
            FOUR_C_THROW("ERROR: unknown shape!");
            break;
          }
        }
        break;
      }
      // #######################################################
      case Core::FE::CellType::hex20:
      {
        switch (target_ele.shape())
        {
          // 2D surface elements
          case Core::FE::CellType::hex8:
          {
            static VolMortarIntegrator<Core::FE::CellType::hex20, Core::FE::CellType::hex8>
                integrator(params());
            integrator.initialize_gp(false, 0, cells[q]->shape());
            integrator.integrate_cells_3d(source_ele, target_ele, *cells[q], *d1_, *m12_, *d2_,
                *m21_, *dis1_, *dis2_, dofset12_.first, dofset12_.second, dofset21_.first,
                dofset21_.second);
            break;
          }
          case Core::FE::CellType::hex20:
          {
            static VolMortarIntegrator<Core::FE::CellType::hex20, Core::FE::CellType::hex20>
                integrator(params());
            integrator.initialize_gp(false, 0, cells[q]->shape());
            integrator.integrate_cells_3d(source_ele, target_ele, *cells[q], *d1_, *m12_, *d2_,
                *m21_, *dis1_, *dis2_, dofset12_.first, dofset12_.second, dofset21_.first,
                dofset21_.second);
            break;
          }
          case Core::FE::CellType::hex27:
          {
            static VolMortarIntegrator<Core::FE::CellType::hex20, Core::FE::CellType::hex27>
                integrator(params());
            integrator.initialize_gp(false, 0, cells[q]->shape());
            integrator.integrate_cells_3d(source_ele, target_ele, *cells[q], *d1_, *m12_, *d2_,
                *m21_, *dis1_, *dis2_, dofset12_.first, dofset12_.second, dofset21_.first,
                dofset21_.second);
            break;
          }
          case Core::FE::CellType::tet4:
          {
            static VolMortarIntegrator<Core::FE::CellType::hex20, Core::FE::CellType::tet4>
                integrator(params());
            integrator.initialize_gp(false, 0, cells[q]->shape());
            integrator.integrate_cells_3d(source_ele, target_ele, *cells[q], *d1_, *m12_, *d2_,
                *m21_, *dis1_, *dis2_, dofset12_.first, dofset12_.second, dofset21_.first,
                dofset21_.second);
            break;
          }
          case Core::FE::CellType::tet10:
          {
            static VolMortarIntegrator<Core::FE::CellType::hex20, Core::FE::CellType::tet10>
                integrator(params());
            integrator.initialize_gp(false, 0, cells[q]->shape());
            integrator.integrate_cells_3d(source_ele, target_ele, *cells[q], *d1_, *m12_, *d2_,
                *m21_, *dis1_, *dis2_, dofset12_.first, dofset12_.second, dofset21_.first,
                dofset21_.second);
            break;
          }
          case Core::FE::CellType::pyramid5:
          {
            static VolMortarIntegrator<Core::FE::CellType::hex20, Core::FE::CellType::pyramid5>
                integrator(params());
            integrator.initialize_gp(false, 0, cells[q]->shape());
            integrator.integrate_cells_3d(source_ele, target_ele, *cells[q], *d1_, *m12_, *d2_,
                *m21_, *dis1_, *dis2_, dofset12_.first, dofset12_.second, dofset21_.first,
                dofset21_.second);
            break;
          }
          default:
          {
            FOUR_C_THROW("ERROR: unknown shape!");
            break;
          }
        }
        break;
      }
      // #######################################################
      case Core::FE::CellType::hex27:
      {
        switch (target_ele.shape())
        {
          // 2D surface elements
          case Core::FE::CellType::hex8:
          {
            static VolMortarIntegrator<Core::FE::CellType::hex27, Core::FE::CellType::hex8>
                integrator(params());
            integrator.initialize_gp(false, 0, cells[q]->shape());
            integrator.integrate_cells_3d(source_ele, target_ele, *cells[q], *d1_, *m12_, *d2_,
                *m21_, *dis1_, *dis2_, dofset12_.first, dofset12_.second, dofset21_.first,
                dofset21_.second);
            break;
          }
          case Core::FE::CellType::hex20:
          {
            static VolMortarIntegrator<Core::FE::CellType::hex27, Core::FE::CellType::hex20>
                integrator(params());
            integrator.initialize_gp(false, 0, cells[q]->shape());
            integrator.integrate_cells_3d(source_ele, target_ele, *cells[q], *d1_, *m12_, *d2_,
                *m21_, *dis1_, *dis2_, dofset12_.first, dofset12_.second, dofset21_.first,
                dofset21_.second);
            break;
          }
          case Core::FE::CellType::hex27:
          {
            static VolMortarIntegrator<Core::FE::CellType::hex27, Core::FE::CellType::hex27>
                integrator(params());
            integrator.initialize_gp(false, 0, cells[q]->shape());
            integrator.integrate_cells_3d(source_ele, target_ele, *cells[q], *d1_, *m12_, *d2_,
                *m21_, *dis1_, *dis2_, dofset12_.first, dofset12_.second, dofset21_.first,
                dofset21_.second);
            break;
          }
          case Core::FE::CellType::tet4:
          {
            static VolMortarIntegrator<Core::FE::CellType::hex27, Core::FE::CellType::tet4>
                integrator(params());
            integrator.initialize_gp(false, 0, cells[q]->shape());
            integrator.integrate_cells_3d(source_ele, target_ele, *cells[q], *d1_, *m12_, *d2_,
                *m21_, *dis1_, *dis2_, dofset12_.first, dofset12_.second, dofset21_.first,
                dofset21_.second);
            break;
          }
          case Core::FE::CellType::tet10:
          {
            static VolMortarIntegrator<Core::FE::CellType::hex27, Core::FE::CellType::tet10>
                integrator(params());
            integrator.initialize_gp(false, 0, cells[q]->shape());
            integrator.integrate_cells_3d(source_ele, target_ele, *cells[q], *d1_, *m12_, *d2_,
                *m21_, *dis1_, *dis2_, dofset12_.first, dofset12_.second, dofset21_.first,
                dofset21_.second);
            break;
          }
          case Core::FE::CellType::pyramid5:
          {
            static VolMortarIntegrator<Core::FE::CellType::hex27, Core::FE::CellType::pyramid5>
                integrator(params());
            integrator.initialize_gp(false, 0, cells[q]->shape());
            integrator.integrate_cells_3d(source_ele, target_ele, *cells[q], *d1_, *m12_, *d2_,
                *m21_, *dis1_, *dis2_, dofset12_.first, dofset12_.second, dofset21_.first,
                dofset21_.second);
            break;
          }
          default:
          {
            FOUR_C_THROW("ERROR: unknown shape!");
            break;
          }
        }
        break;
      }
      // #######################################################
      case Core::FE::CellType::tet4:
      {
        switch (target_ele.shape())
        {
          // 2D surface elements
          case Core::FE::CellType::hex8:
          {
            static VolMortarIntegrator<Core::FE::CellType::tet4, Core::FE::CellType::hex8>
                integrator(params());
            integrator.initialize_gp(false, 0, cells[q]->shape());
            integrator.integrate_cells_3d(source_ele, target_ele, *cells[q], *d1_, *m12_, *d2_,
                *m21_, *dis1_, *dis2_, dofset12_.first, dofset12_.second, dofset21_.first,
                dofset21_.second);
            break;
          }
          case Core::FE::CellType::hex20:
          {
            static VolMortarIntegrator<Core::FE::CellType::tet4, Core::FE::CellType::hex20>
                integrator(params());
            integrator.initialize_gp(false, 0, cells[q]->shape());
            integrator.integrate_cells_3d(source_ele, target_ele, *cells[q], *d1_, *m12_, *d2_,
                *m21_, *dis1_, *dis2_, dofset12_.first, dofset12_.second, dofset21_.first,
                dofset21_.second);
            break;
          }
          case Core::FE::CellType::hex27:
          {
            static VolMortarIntegrator<Core::FE::CellType::tet4, Core::FE::CellType::hex27>
                integrator(params());
            integrator.initialize_gp(false, 0, cells[q]->shape());
            integrator.integrate_cells_3d(source_ele, target_ele, *cells[q], *d1_, *m12_, *d2_,
                *m21_, *dis1_, *dis2_, dofset12_.first, dofset12_.second, dofset21_.first,
                dofset21_.second);
            break;
          }
          case Core::FE::CellType::tet4:
          {
            static VolMortarIntegrator<Core::FE::CellType::tet4, Core::FE::CellType::tet4>
                integrator(params());
            integrator.initialize_gp(false, 0, cells[q]->shape());
            integrator.integrate_cells_3d(source_ele, target_ele, *cells[q], *d1_, *m12_, *d2_,
                *m21_, *dis1_, *dis2_, dofset12_.first, dofset12_.second, dofset21_.first,
                dofset21_.second);
            break;
          }
          case Core::FE::CellType::tet10:
          {
            static VolMortarIntegrator<Core::FE::CellType::tet4, Core::FE::CellType::tet10>
                integrator(params());
            integrator.initialize_gp(false, 0, cells[q]->shape());
            integrator.integrate_cells_3d(source_ele, target_ele, *cells[q], *d1_, *m12_, *d2_,
                *m21_, *dis1_, *dis2_, dofset12_.first, dofset12_.second, dofset21_.first,
                dofset21_.second);
            break;
          }
          case Core::FE::CellType::pyramid5:
          {
            static VolMortarIntegrator<Core::FE::CellType::tet4, Core::FE::CellType::pyramid5>
                integrator(params());
            integrator.initialize_gp(false, 0, cells[q]->shape());
            integrator.integrate_cells_3d(source_ele, target_ele, *cells[q], *d1_, *m12_, *d2_,
                *m21_, *dis1_, *dis2_, dofset12_.first, dofset12_.second, dofset21_.first,
                dofset21_.second);
            break;
          }
          default:
          {
            FOUR_C_THROW("ERROR: unknown shape!");
            break;
          }
        }
        break;
      }
      // #######################################################
      case Core::FE::CellType::tet10:
      {
        switch (target_ele.shape())
        {
          // 2D surface elements
          case Core::FE::CellType::hex8:
          {
            static VolMortarIntegrator<Core::FE::CellType::tet10, Core::FE::CellType::hex8>
                integrator(params());
            integrator.initialize_gp(false, 0, cells[q]->shape());
            integrator.integrate_cells_3d(source_ele, target_ele, *cells[q], *d1_, *m12_, *d2_,
                *m21_, *dis1_, *dis2_, dofset12_.first, dofset12_.second, dofset21_.first,
                dofset21_.second);
            break;
          }
          case Core::FE::CellType::hex20:
          {
            static VolMortarIntegrator<Core::FE::CellType::tet10, Core::FE::CellType::hex20>
                integrator(params());
            integrator.initialize_gp(false, 0, cells[q]->shape());
            integrator.integrate_cells_3d(source_ele, target_ele, *cells[q], *d1_, *m12_, *d2_,
                *m21_, *dis1_, *dis2_, dofset12_.first, dofset12_.second, dofset21_.first,
                dofset21_.second);
            break;
          }
          case Core::FE::CellType::hex27:
          {
            static VolMortarIntegrator<Core::FE::CellType::tet10, Core::FE::CellType::hex27>
                integrator(params());
            integrator.initialize_gp(false, 0, cells[q]->shape());
            integrator.integrate_cells_3d(source_ele, target_ele, *cells[q], *d1_, *m12_, *d2_,
                *m21_, *dis1_, *dis2_, dofset12_.first, dofset12_.second, dofset21_.first,
                dofset21_.second);
            break;
          }
          case Core::FE::CellType::tet4:
          {
            static VolMortarIntegrator<Core::FE::CellType::tet10, Core::FE::CellType::tet4>
                integrator(params());
            integrator.initialize_gp(false, 0, cells[q]->shape());
            integrator.integrate_cells_3d(source_ele, target_ele, *cells[q], *d1_, *m12_, *d2_,
                *m21_, *dis1_, *dis2_, dofset12_.first, dofset12_.second, dofset21_.first,
                dofset21_.second);
            break;
          }
          case Core::FE::CellType::tet10:
          {
            static VolMortarIntegrator<Core::FE::CellType::tet10, Core::FE::CellType::tet10>
                integrator(params());
            integrator.initialize_gp(false, 0, cells[q]->shape());
            integrator.integrate_cells_3d(source_ele, target_ele, *cells[q], *d1_, *m12_, *d2_,
                *m21_, *dis1_, *dis2_, dofset12_.first, dofset12_.second, dofset21_.first,
                dofset21_.second);
            break;
          }
          case Core::FE::CellType::pyramid5:
          {
            static VolMortarIntegrator<Core::FE::CellType::tet10, Core::FE::CellType::pyramid5>
                integrator(params());
            integrator.initialize_gp(false, 0, cells[q]->shape());
            integrator.integrate_cells_3d(source_ele, target_ele, *cells[q], *d1_, *m12_, *d2_,
                *m21_, *dis1_, *dis2_, dofset12_.first, dofset12_.second, dofset21_.first,
                dofset21_.second);
            break;
          }
          default:
          {
            FOUR_C_THROW("ERROR: unknown shape!");
            break;
          }
        }
        break;
      }
      // #######################################################
      case Core::FE::CellType::pyramid5:
      {
        switch (target_ele.shape())
        {
          // 2D surface elements
          case Core::FE::CellType::hex8:
          {
            static VolMortarIntegrator<Core::FE::CellType::pyramid5, Core::FE::CellType::hex8>
                integrator(params());
            integrator.initialize_gp(false, 0, cells[q]->shape());
            integrator.integrate_cells_3d(source_ele, target_ele, *cells[q], *d1_, *m12_, *d2_,
                *m21_, *dis1_, *dis2_, dofset12_.first, dofset12_.second, dofset21_.first,
                dofset21_.second);
            break;
          }
          case Core::FE::CellType::hex20:
          {
            static VolMortarIntegrator<Core::FE::CellType::pyramid5, Core::FE::CellType::hex20>
                integrator(params());
            integrator.initialize_gp(false, 0, cells[q]->shape());
            integrator.integrate_cells_3d(source_ele, target_ele, *cells[q], *d1_, *m12_, *d2_,
                *m21_, *dis1_, *dis2_, dofset12_.first, dofset12_.second, dofset21_.first,
                dofset21_.second);
            break;
          }
          case Core::FE::CellType::hex27:
          {
            static VolMortarIntegrator<Core::FE::CellType::pyramid5, Core::FE::CellType::hex27>
                integrator(params());
            integrator.initialize_gp(false, 0, cells[q]->shape());
            integrator.integrate_cells_3d(source_ele, target_ele, *cells[q], *d1_, *m12_, *d2_,
                *m21_, *dis1_, *dis2_, dofset12_.first, dofset12_.second, dofset21_.first,
                dofset21_.second);
            break;
          }
          case Core::FE::CellType::tet4:
          {
            static VolMortarIntegrator<Core::FE::CellType::pyramid5, Core::FE::CellType::tet4>
                integrator(params());
            integrator.initialize_gp(false, 0, cells[q]->shape());
            integrator.integrate_cells_3d(source_ele, target_ele, *cells[q], *d1_, *m12_, *d2_,
                *m21_, *dis1_, *dis2_, dofset12_.first, dofset12_.second, dofset21_.first,
                dofset21_.second);
            break;
          }
          case Core::FE::CellType::tet10:
          {
            static VolMortarIntegrator<Core::FE::CellType::pyramid5, Core::FE::CellType::tet10>
                integrator(params());
            integrator.initialize_gp(false, 0, cells[q]->shape());
            integrator.integrate_cells_3d(source_ele, target_ele, *cells[q], *d1_, *m12_, *d2_,
                *m21_, *dis1_, *dis2_, dofset12_.first, dofset12_.second, dofset21_.first,
                dofset21_.second);
            break;
          }
          case Core::FE::CellType::pyramid5:
          {
            static VolMortarIntegrator<Core::FE::CellType::pyramid5, Core::FE::CellType::pyramid5>
                integrator(params());
            integrator.initialize_gp(false, 0, cells[q]->shape());
            integrator.integrate_cells_3d(source_ele, target_ele, *cells[q], *d1_, *m12_, *d2_,
                *m21_, *dis1_, *dis2_, dofset12_.first, dofset12_.second, dofset21_.first,
                dofset21_.second);
            break;
          }
          default:
          {
            FOUR_C_THROW("ERROR: unknown shape!");
            break;
          }
        }
        break;
      }
      default:
      {
        FOUR_C_THROW("ERROR: unknown shape!");
        break;
      }
    }
  }  // end cell loop
}

/*----------------------------------------------------------------------*
 |  integrate_3d elebased                                     farah 04/14|
 *----------------------------------------------------------------------*/
void Coupling::VolMortar::VolMortarCoupl::integrate_3d_ele_based_p12(
    Core::Elements::Element& Aele, std::vector<int>& foundeles)
{
  switch (Aele.shape())
  {
    // 2D "volume" elements
    case Core::FE::CellType::quad4:
    {
      static VolMortarIntegratorEleBased<Core::FE::CellType::quad4> integrator(params());
      integrator.initialize_gp();
      integrator.integrate_ele_based_3d(Aele, foundeles, *d1_, *m12_, *dis1_, *dis2_,
          dofset12_.first, dofset12_.second, *p12_dofrowmap_, *p12_dofcolmap_);
      break;
    }
    case Core::FE::CellType::quad8:
    {
      static VolMortarIntegratorEleBased<Core::FE::CellType::quad8> integrator(params());
      integrator.initialize_gp();
      integrator.integrate_ele_based_3d(Aele, foundeles, *d1_, *m12_, *dis1_, *dis2_,
          dofset12_.first, dofset12_.second, *p12_dofrowmap_, *p12_dofcolmap_);
      break;
    }
    case Core::FE::CellType::quad9:
    {
      static VolMortarIntegratorEleBased<Core::FE::CellType::quad9> integrator(params());
      integrator.initialize_gp();
      integrator.integrate_ele_based_3d(Aele, foundeles, *d1_, *m12_, *dis1_, *dis2_,
          dofset12_.first, dofset12_.second, *p12_dofrowmap_, *p12_dofcolmap_);
      break;
    }
    case Core::FE::CellType::tri3:
    {
      static VolMortarIntegratorEleBased<Core::FE::CellType::tri3> integrator(params());
      integrator.initialize_gp();
      integrator.integrate_ele_based_3d(Aele, foundeles, *d1_, *m12_, *dis1_, *dis2_,
          dofset12_.first, dofset12_.second, *p12_dofrowmap_, *p12_dofcolmap_);
      break;
    }
    case Core::FE::CellType::tri6:
    {
      static VolMortarIntegratorEleBased<Core::FE::CellType::tri6> integrator(params());
      integrator.initialize_gp();
      integrator.integrate_ele_based_3d(Aele, foundeles, *d1_, *m12_, *dis1_, *dis2_,
          dofset12_.first, dofset12_.second, *p12_dofrowmap_, *p12_dofcolmap_);
      break;
    }
    // 3D volume elements
    case Core::FE::CellType::hex8:
    {
      static VolMortarIntegratorEleBased<Core::FE::CellType::hex8> integrator(params());
      integrator.initialize_gp();
      integrator.integrate_ele_based_3d(Aele, foundeles, *d1_, *m12_, *dis1_, *dis2_,
          dofset12_.first, dofset12_.second, *p12_dofrowmap_, *p12_dofcolmap_);
      break;
    }
    case Core::FE::CellType::hex27:
    {
      static VolMortarIntegratorEleBased<Core::FE::CellType::hex27> integrator(params());
      integrator.initialize_gp();
      integrator.integrate_ele_based_3d(Aele, foundeles, *d1_, *m12_, *dis1_, *dis2_,
          dofset12_.first, dofset12_.second, *p12_dofrowmap_, *p12_dofcolmap_);
      break;
    }
    case Core::FE::CellType::hex20:
    {
      static VolMortarIntegratorEleBased<Core::FE::CellType::hex20> integrator(params());
      integrator.initialize_gp();
      integrator.integrate_ele_based_3d(Aele, foundeles, *d1_, *m12_, *dis1_, *dis2_,
          dofset12_.first, dofset12_.second, *p12_dofrowmap_, *p12_dofcolmap_);
      break;
    }
    case Core::FE::CellType::tet4:
    {
      static VolMortarIntegratorEleBased<Core::FE::CellType::tet4> integrator(params());
      integrator.initialize_gp();
      integrator.integrate_ele_based_3d(Aele, foundeles, *d1_, *m12_, *dis1_, *dis2_,
          dofset12_.first, dofset12_.second, *p12_dofrowmap_, *p12_dofcolmap_);
      break;
    }
    case Core::FE::CellType::tet10:
    {
      static VolMortarIntegratorEleBased<Core::FE::CellType::tet10> integrator(params());
      integrator.initialize_gp();
      integrator.integrate_ele_based_3d(Aele, foundeles, *d1_, *m12_, *dis1_, *dis2_,
          dofset12_.first, dofset12_.second, *p12_dofrowmap_, *p12_dofcolmap_);
      break;
    }
    case Core::FE::CellType::pyramid5:
    {
      static VolMortarIntegratorEleBased<Core::FE::CellType::pyramid5> integrator(params());
      integrator.initialize_gp();
      integrator.integrate_ele_based_3d(Aele, foundeles, *d1_, *m12_, *dis1_, *dis2_,
          dofset12_.first, dofset12_.second, *p12_dofrowmap_, *p12_dofcolmap_);
      break;
    }
    default:
    {
      FOUR_C_THROW("ERROR: unknown shape!");
      break;
    }
  }
}

/*----------------------------------------------------------------------*
 |  integrate_3d element based for projector P21              farah 04/14|
 *----------------------------------------------------------------------*/
void Coupling::VolMortar::VolMortarCoupl::integrate_3d_ele_based_p21(
    Core::Elements::Element& Bele, std::vector<int>& foundeles)
{
  switch (Bele.shape())
  {
    // 2D "volume" elements
    case Core::FE::CellType::quad4:
    {
      static VolMortarIntegratorEleBased<Core::FE::CellType::quad4> integrator(params());
      integrator.initialize_gp();
      integrator.integrate_ele_based_3d(Bele, foundeles, *d2_, *m21_, *dis2_, *dis1_,
          dofset21_.first, dofset21_.second, *p21_dofrowmap_, *p21_dofcolmap_);
      break;
    }
    case Core::FE::CellType::quad8:
    {
      static VolMortarIntegratorEleBased<Core::FE::CellType::quad8> integrator(params());
      integrator.initialize_gp();
      integrator.integrate_ele_based_3d(Bele, foundeles, *d2_, *m21_, *dis2_, *dis1_,
          dofset21_.first, dofset21_.second, *p21_dofrowmap_, *p21_dofcolmap_);
      break;
    }
    case Core::FE::CellType::quad9:
    {
      static VolMortarIntegratorEleBased<Core::FE::CellType::quad9> integrator(params());
      integrator.initialize_gp();
      integrator.integrate_ele_based_3d(Bele, foundeles, *d2_, *m21_, *dis2_, *dis1_,
          dofset21_.first, dofset21_.second, *p21_dofrowmap_, *p21_dofcolmap_);
      break;
    }
    case Core::FE::CellType::tri3:
    {
      static VolMortarIntegratorEleBased<Core::FE::CellType::tri3> integrator(params());
      integrator.initialize_gp();
      integrator.integrate_ele_based_3d(Bele, foundeles, *d2_, *m21_, *dis2_, *dis1_,
          dofset21_.first, dofset21_.second, *p21_dofrowmap_, *p21_dofcolmap_);
      break;
    }
    case Core::FE::CellType::tri6:
    {
      static VolMortarIntegratorEleBased<Core::FE::CellType::tri6> integrator(params());
      integrator.initialize_gp();
      integrator.integrate_ele_based_3d(Bele, foundeles, *d2_, *m21_, *dis2_, *dis1_,
          dofset21_.first, dofset21_.second, *p21_dofrowmap_, *p21_dofcolmap_);
      break;
    }
    // 3D volume elements
    case Core::FE::CellType::hex8:
    {
      static VolMortarIntegratorEleBased<Core::FE::CellType::hex8> integrator(params());
      integrator.initialize_gp();
      integrator.integrate_ele_based_3d(Bele, foundeles, *d2_, *m21_, *dis2_, *dis1_,
          dofset21_.first, dofset21_.second, *p21_dofrowmap_, *p21_dofcolmap_);
      break;
    }
    case Core::FE::CellType::hex20:
    {
      static VolMortarIntegratorEleBased<Core::FE::CellType::hex20> integrator(params());
      integrator.initialize_gp();
      integrator.integrate_ele_based_3d(Bele, foundeles, *d2_, *m21_, *dis2_, *dis1_,
          dofset21_.first, dofset21_.second, *p21_dofrowmap_, *p21_dofcolmap_);
      break;
    }
    case Core::FE::CellType::hex27:
    {
      static VolMortarIntegratorEleBased<Core::FE::CellType::hex27> integrator(params());
      integrator.initialize_gp();
      integrator.integrate_ele_based_3d(Bele, foundeles, *d2_, *m21_, *dis2_, *dis1_,
          dofset21_.first, dofset21_.second, *p21_dofrowmap_, *p21_dofcolmap_);
      break;
    }
    case Core::FE::CellType::tet4:
    {
      static VolMortarIntegratorEleBased<Core::FE::CellType::tet4> integrator(params());
      integrator.initialize_gp();
      integrator.integrate_ele_based_3d(Bele, foundeles, *d2_, *m21_, *dis2_, *dis1_,
          dofset21_.first, dofset21_.second, *p21_dofrowmap_, *p21_dofcolmap_);
      break;
    }
    case Core::FE::CellType::tet10:
    {
      static VolMortarIntegratorEleBased<Core::FE::CellType::tet10> integrator(params());
      integrator.initialize_gp();
      integrator.integrate_ele_based_3d(Bele, foundeles, *d2_, *m21_, *dis2_, *dis1_,
          dofset21_.first, dofset21_.second, *p21_dofrowmap_, *p21_dofcolmap_);
      break;
    }
    case Core::FE::CellType::pyramid5:
    {
      static VolMortarIntegratorEleBased<Core::FE::CellType::pyramid5> integrator(params());
      integrator.initialize_gp();
      integrator.integrate_ele_based_3d(Bele, foundeles, *d2_, *m21_, *dis2_, *dis1_,
          dofset21_.first, dofset21_.second, *p21_dofrowmap_, *p21_dofcolmap_);
      break;
    }
    default:
    {
      FOUR_C_THROW("ERROR: unknown shape!");
      break;
    }
  }
}

/*----------------------------------------------------------------------*
 |  integrate_3d elebased                                     farah 03/15|
 *----------------------------------------------------------------------*/
void Coupling::VolMortar::VolMortarCoupl::integrate_3d_ele_based_a_dis_mesh_init(
    Core::Elements::Element& Aele, std::vector<int>& foundeles, int dofseta, int dofsetb)
{
  switch (Aele.shape())
  {
    // 2D surface elements
    case Core::FE::CellType::hex8:
    {
      static VolMortarIntegratorEleBased<Core::FE::CellType::hex8> integrator(params());
      integrator.initialize_gp();
      integrator.integrate_ele_based_3d(Aele, foundeles, *dmatrix_xa_, *mmatrix_xa_, *dis1_, *dis2_,
          dofseta, dofsetb, *p12_dofrowmap_, *p12_dofcolmap_);
      break;
    }
    case Core::FE::CellType::tet4:
    {
      static VolMortarIntegratorEleBased<Core::FE::CellType::tet4> integrator(params());
      integrator.initialize_gp();
      integrator.integrate_ele_based_3d(Aele, foundeles, *dmatrix_xa_, *mmatrix_xa_, *dis1_, *dis2_,
          dofseta, dofsetb, *p12_dofrowmap_, *p12_dofcolmap_);
      break;
    }
    case Core::FE::CellType::hex27:
    {
      static VolMortarIntegratorEleBased<Core::FE::CellType::hex27> integrator(params());
      integrator.initialize_gp();
      integrator.integrate_ele_based_3d(Aele, foundeles, *dmatrix_xa_, *mmatrix_xa_, *dis1_, *dis2_,
          dofseta, dofsetb, *p12_dofrowmap_, *p12_dofcolmap_);
      break;
    }
    case Core::FE::CellType::hex20:
    {
      static VolMortarIntegratorEleBased<Core::FE::CellType::hex20> integrator(params());
      integrator.initialize_gp();
      integrator.integrate_ele_based_3d(Aele, foundeles, *dmatrix_xa_, *mmatrix_xa_, *dis1_, *dis2_,
          dofseta, dofsetb, *p12_dofrowmap_, *p12_dofcolmap_);
      break;
    }
    case Core::FE::CellType::tet10:
    {
      static VolMortarIntegratorEleBased<Core::FE::CellType::tet10> integrator(params());
      integrator.initialize_gp();
      integrator.integrate_ele_based_3d(Aele, foundeles, *dmatrix_xa_, *mmatrix_xa_, *dis1_, *dis2_,
          dofseta, dofsetb, *p12_dofrowmap_, *p12_dofcolmap_);
      break;
    }
    case Core::FE::CellType::pyramid5:
    {
      static VolMortarIntegratorEleBased<Core::FE::CellType::pyramid5> integrator(params());
      integrator.initialize_gp();
      integrator.integrate_ele_based_3d(Aele, foundeles, *dmatrix_xa_, *mmatrix_xa_, *dis1_, *dis2_,
          dofseta, dofsetb, *p12_dofrowmap_, *p12_dofcolmap_);
      break;
    }
    default:
    {
      FOUR_C_THROW("ERROR: unknown shape!");
      break;
    }
  }
}

/*----------------------------------------------------------------------*
 |  integrate_3d Cells                                        farah 03/15|
 *----------------------------------------------------------------------*/
void Coupling::VolMortar::VolMortarCoupl::integrate_3d_ele_based_b_dis_mesh_init(
    Core::Elements::Element& Bele, std::vector<int>& foundeles, int dofsetb, int dofseta)
{
  switch (Bele.shape())
  {
    // 2D surface elements
    case Core::FE::CellType::hex8:
    {
      static VolMortarIntegratorEleBased<Core::FE::CellType::hex8> integrator(params());
      integrator.initialize_gp();
      integrator.integrate_ele_based_3d(Bele, foundeles, *dmatrix_xb_, *mmatrix_xb_, *dis2_, *dis1_,
          dofseta, dofsetb, *p21_dofrowmap_, *p21_dofcolmap_);
      break;
    }
    case Core::FE::CellType::tet4:
    {
      static VolMortarIntegratorEleBased<Core::FE::CellType::tet4> integrator(params());
      integrator.initialize_gp();
      integrator.integrate_ele_based_3d(Bele, foundeles, *dmatrix_xb_, *mmatrix_xb_, *dis2_, *dis1_,
          dofseta, dofsetb, *p21_dofrowmap_, *p21_dofcolmap_);
      break;
    }
    case Core::FE::CellType::hex27:
    {
      static VolMortarIntegratorEleBased<Core::FE::CellType::hex27> integrator(params());
      integrator.initialize_gp();
      integrator.integrate_ele_based_3d(Bele, foundeles, *dmatrix_xb_, *mmatrix_xb_, *dis2_, *dis1_,
          dofseta, dofsetb, *p21_dofrowmap_, *p21_dofcolmap_);
      break;
    }
    case Core::FE::CellType::hex20:
    {
      static VolMortarIntegratorEleBased<Core::FE::CellType::hex20> integrator(params());
      integrator.initialize_gp();
      integrator.integrate_ele_based_3d(Bele, foundeles, *dmatrix_xb_, *mmatrix_xb_, *dis2_, *dis1_,
          dofseta, dofsetb, *p21_dofrowmap_, *p21_dofcolmap_);
      break;
    }
    case Core::FE::CellType::tet10:
    {
      static VolMortarIntegratorEleBased<Core::FE::CellType::tet10> integrator(params());
      integrator.initialize_gp();
      integrator.integrate_ele_based_3d(Bele, foundeles, *dmatrix_xb_, *mmatrix_xb_, *dis2_, *dis1_,
          dofseta, dofsetb, *p21_dofrowmap_, *p21_dofcolmap_);
      break;
    }
    case Core::FE::CellType::pyramid5:
    {
      static VolMortarIntegratorEleBased<Core::FE::CellType::pyramid5> integrator(params());
      integrator.initialize_gp();
      integrator.integrate_ele_based_3d(Bele, foundeles, *dmatrix_xb_, *mmatrix_xb_, *dis2_, *dis1_,
          dofseta, dofsetb, *p21_dofrowmap_, *p21_dofcolmap_);
      break;
    }
    default:
    {
      FOUR_C_THROW("ERROR: unknown shape!");
      break;
    }
  }
}

/*----------------------------------------------------------------------*
 |  Assemble p matrix for cons. interpolation approach       farah 06/14|
 *----------------------------------------------------------------------*/
void Coupling::VolMortar::VolMortarCoupl::assemble_consistent_interpolation_p12(
    Core::Nodes::Node* node, std::vector<int>& foundeles)
{
  static ConsInterpolator interpolator;
  interpolator.interpolate(
      node, *p12_, *dis1_, *dis2_, foundeles, dofset12_, *p12_dofrowmap_, *p12_dofcolmap_);
}

/*----------------------------------------------------------------------*
 |  Assemble p matrix for cons. interpolation approach       farah 06/14|
 *----------------------------------------------------------------------*/
void Coupling::VolMortar::VolMortarCoupl::assemble_consistent_interpolation_p21(
    Core::Nodes::Node* node, std::vector<int>& foundeles)
{
  static ConsInterpolator interpolator;
  interpolator.interpolate(
      node, *p21_, *dis2_, *dis1_, foundeles, dofset21_, *p21_dofrowmap_, *p21_dofcolmap_);
}

/*----------------------------------------------------------------------*
 |  integrate_3d Cells for direct divergence approach         farah 04/14|
 *----------------------------------------------------------------------*/
void Coupling::VolMortar::VolMortarCoupl::integrate_3d_cell_direct_divergence(
    Core::Elements::Element& source_ele, Core::Elements::Element& target_ele, bool switched_conf)
{
  if ((int)(volcell_.size()) > 1)
    std::cout << "****************************   CELL SIZE > 1 ***************************"
              << std::endl;

  for (Cut::plain_volumecell_set::iterator i = volcell_.begin(); i != volcell_.end(); i++)
  {
    if (*i == nullptr) continue;

    Cut::VolumeCell* vc = *i;

    if (vc->is_negligibly_small()) continue;

    // main gp rule
    std::shared_ptr<Core::FE::GaussPoints> intpoints = vc->get_gauss_rule();

    //--------------------------------------------------------------------
    // loop over cells for A Field

    switch (source_ele.shape())
    {
      // 2D surface elements
      case Core::FE::CellType::hex8:
      {
        switch (target_ele.shape())
        {
          // 2D surface elements
          case Core::FE::CellType::hex8:
          {
            static VolMortarIntegrator<Core::FE::CellType::hex8, Core::FE::CellType::hex8>
                integrator(params());
            integrator.integrate_cells_3d_direct_diveregence(source_ele, target_ele, *vc, intpoints,
                switched_conf, *d1_, *m12_, *d2_, *m21_, *dis1_, *dis2_, dofset12_.first,
                dofset12_.second, dofset21_.first, dofset21_.second);
            break;
          }
          case Core::FE::CellType::tet4:
          {
            static VolMortarIntegrator<Core::FE::CellType::hex8, Core::FE::CellType::tet4>
                integrator(params());
            integrator.integrate_cells_3d_direct_diveregence(source_ele, target_ele, *vc, intpoints,
                switched_conf, *d1_, *m12_, *d2_, *m21_, *dis1_, *dis2_, dofset12_.first,
                dofset12_.second, dofset21_.first, dofset21_.second);
            break;
          }
          case Core::FE::CellType::pyramid5:
          {
            static VolMortarIntegrator<Core::FE::CellType::hex8, Core::FE::CellType::pyramid5>
                integrator(params());
            integrator.integrate_cells_3d_direct_diveregence(source_ele, target_ele, *vc, intpoints,
                switched_conf, *d1_, *m12_, *d2_, *m21_, *dis1_, *dis2_, dofset12_.first,
                dofset12_.second, dofset21_.first, dofset21_.second);
            break;
          }
          default:
          {
            FOUR_C_THROW("ERROR: unknown shape!");
            break;
          }
        }
        break;
      }
      case Core::FE::CellType::tet4:
      {
        switch (target_ele.shape())
        {
          // 2D surface elements
          case Core::FE::CellType::hex8:
          {
            static VolMortarIntegrator<Core::FE::CellType::tet4, Core::FE::CellType::hex8>
                integrator(params());
            integrator.integrate_cells_3d_direct_diveregence(source_ele, target_ele, *vc, intpoints,
                switched_conf, *d1_, *m12_, *d2_, *m21_, *dis1_, *dis2_, dofset12_.first,
                dofset12_.second, dofset21_.first, dofset21_.second);
            break;
          }
          case Core::FE::CellType::tet4:
          {
            static VolMortarIntegrator<Core::FE::CellType::tet4, Core::FE::CellType::tet4>
                integrator(params());
            integrator.integrate_cells_3d_direct_diveregence(source_ele, target_ele, *vc, intpoints,
                switched_conf, *d1_, *m12_, *d2_, *m21_, *dis1_, *dis2_, dofset12_.first,
                dofset12_.second, dofset21_.first, dofset21_.second);
            break;
          }
          case Core::FE::CellType::pyramid5:
          {
            static VolMortarIntegrator<Core::FE::CellType::tet4, Core::FE::CellType::pyramid5>
                integrator(params());
            integrator.integrate_cells_3d_direct_diveregence(source_ele, target_ele, *vc, intpoints,
                switched_conf, *d1_, *m12_, *d2_, *m21_, *dis1_, *dis2_, dofset12_.first,
                dofset12_.second, dofset21_.first, dofset21_.second);
            break;
          }
          default:
          {
            FOUR_C_THROW("ERROR: unknown shape!");
            break;
          }
        }
        break;
      }
      case Core::FE::CellType::pyramid5:
      {
        switch (target_ele.shape())
        {
          // 2D surface elements
          case Core::FE::CellType::hex8:
          {
            static VolMortarIntegrator<Core::FE::CellType::pyramid5, Core::FE::CellType::hex8>
                integrator(params());
            integrator.integrate_cells_3d_direct_diveregence(source_ele, target_ele, *vc, intpoints,
                switched_conf, *d1_, *m12_, *d2_, *m21_, *dis1_, *dis2_, dofset12_.first,
                dofset12_.second, dofset21_.first, dofset21_.second);
            break;
          }
          case Core::FE::CellType::tet4:
          {
            static VolMortarIntegrator<Core::FE::CellType::pyramid5, Core::FE::CellType::tet4>
                integrator(params());
            integrator.integrate_cells_3d_direct_diveregence(source_ele, target_ele, *vc, intpoints,
                switched_conf, *d1_, *m12_, *d2_, *m21_, *dis1_, *dis2_, dofset12_.first,
                dofset12_.second, dofset21_.first, dofset21_.second);
            break;
          }
          case Core::FE::CellType::pyramid5:
          {
            static VolMortarIntegrator<Core::FE::CellType::pyramid5, Core::FE::CellType::pyramid5>
                integrator(params());
            integrator.integrate_cells_3d_direct_diveregence(source_ele, target_ele, *vc, intpoints,
                switched_conf, *d1_, *m12_, *d2_, *m21_, *dis1_, *dis2_, dofset12_.first,
                dofset12_.second, dofset21_.first, dofset21_.second);
            break;
          }
          default:
          {
            FOUR_C_THROW("ERROR: unknown shape!");
            break;
          }
        }
        break;
      }
      default:
      {
        FOUR_C_THROW("ERROR: unknown shape!");
        break;
      }
    }
  }
}

/*----------------------------------------------------------------------*
 |  Integrate over A-element domain                          farah 01/14|
 *----------------------------------------------------------------------*/
void Coupling::VolMortar::VolMortarCoupl::integrate_3d(
    Core::Elements::Element& source_ele, Core::Elements::Element& target_ele, int domain)
{
  //--------------------------------------------------------------------
  // loop over cells for A Field
  switch (source_ele.shape())
  {
    // 2D surface elements
    case Core::FE::CellType::hex8:
    {
      switch (target_ele.shape())
      {
        // 2D surface elements
        case Core::FE::CellType::hex8:
        {
          static VolMortarIntegrator<Core::FE::CellType::hex8, Core::FE::CellType::hex8> integrator(
              params());
          integrator.initialize_gp(true, domain);
          integrator.integrate_ele_3d(domain, source_ele, target_ele, *d1_, *m12_, *d2_, *m21_,
              *dis1_, *dis2_, dofset12_.first, dofset12_.second, dofset21_.first, dofset21_.second);
          break;
        }
        case Core::FE::CellType::hex20:
        {
          static VolMortarIntegrator<Core::FE::CellType::hex8, Core::FE::CellType::hex20>
              integrator(params());
          integrator.initialize_gp(true, domain);
          integrator.integrate_ele_3d(domain, source_ele, target_ele, *d1_, *m12_, *d2_, *m21_,
              *dis1_, *dis2_, dofset12_.first, dofset12_.second, dofset21_.first, dofset21_.second);
          break;
        }
        case Core::FE::CellType::hex27:
        {
          static VolMortarIntegrator<Core::FE::CellType::hex8, Core::FE::CellType::hex27>
              integrator(params());
          integrator.initialize_gp(true, domain);
          integrator.integrate_ele_3d(domain, source_ele, target_ele, *d1_, *m12_, *d2_, *m21_,
              *dis1_, *dis2_, dofset12_.first, dofset12_.second, dofset21_.first, dofset21_.second);
          break;
        }
        case Core::FE::CellType::tet4:
        {
          static VolMortarIntegrator<Core::FE::CellType::hex8, Core::FE::CellType::tet4> integrator(
              params());
          integrator.initialize_gp(true, domain);
          integrator.integrate_ele_3d(domain, source_ele, target_ele, *d1_, *m12_, *d2_, *m21_,
              *dis1_, *dis2_, dofset12_.first, dofset12_.second, dofset21_.first, dofset21_.second);
          break;
        }
        case Core::FE::CellType::tet10:
        {
          static VolMortarIntegrator<Core::FE::CellType::hex8, Core::FE::CellType::tet10>
              integrator(params());
          integrator.initialize_gp(true, domain);
          integrator.integrate_ele_3d(domain, source_ele, target_ele, *d1_, *m12_, *d2_, *m21_,
              *dis1_, *dis2_, dofset12_.first, dofset12_.second, dofset21_.first, dofset21_.second);
          break;
        }
        case Core::FE::CellType::pyramid5:
        {
          static VolMortarIntegrator<Core::FE::CellType::hex8, Core::FE::CellType::pyramid5>
              integrator(params());
          integrator.initialize_gp(true, domain);
          integrator.integrate_ele_3d(domain, source_ele, target_ele, *d1_, *m12_, *d2_, *m21_,
              *dis1_, *dis2_, dofset12_.first, dofset12_.second, dofset21_.first, dofset21_.second);
          break;
        }
        default:
        {
          FOUR_C_THROW("ERROR: unknown shape!");
          break;
        }
      }
      break;
    }
    // ########################################################
    case Core::FE::CellType::hex20:
    {
      switch (target_ele.shape())
      {
        // 2D surface elements
        case Core::FE::CellType::hex8:
        {
          static VolMortarIntegrator<Core::FE::CellType::hex20, Core::FE::CellType::hex8>
              integrator(params());
          integrator.initialize_gp(true, domain);
          integrator.integrate_ele_3d(domain, source_ele, target_ele, *d1_, *m12_, *d2_, *m21_,
              *dis1_, *dis2_, dofset12_.first, dofset12_.second, dofset21_.first, dofset21_.second);
          break;
        }
        case Core::FE::CellType::hex20:
        {
          static VolMortarIntegrator<Core::FE::CellType::hex20, Core::FE::CellType::hex20>
              integrator(params());
          integrator.initialize_gp(true, domain);
          integrator.integrate_ele_3d(domain, source_ele, target_ele, *d1_, *m12_, *d2_, *m21_,
              *dis1_, *dis2_, dofset12_.first, dofset12_.second, dofset21_.first, dofset21_.second);
          break;
        }
        case Core::FE::CellType::hex27:
        {
          static VolMortarIntegrator<Core::FE::CellType::hex20, Core::FE::CellType::hex27>
              integrator(params());
          integrator.initialize_gp(true, domain);
          integrator.integrate_ele_3d(domain, source_ele, target_ele, *d1_, *m12_, *d2_, *m21_,
              *dis1_, *dis2_, dofset12_.first, dofset12_.second, dofset21_.first, dofset21_.second);
          break;
        }
        case Core::FE::CellType::tet4:
        {
          static VolMortarIntegrator<Core::FE::CellType::hex20, Core::FE::CellType::tet4>
              integrator(params());
          integrator.initialize_gp(true, domain);
          integrator.integrate_ele_3d(domain, source_ele, target_ele, *d1_, *m12_, *d2_, *m21_,
              *dis1_, *dis2_, dofset12_.first, dofset12_.second, dofset21_.first, dofset21_.second);
          break;
        }
        case Core::FE::CellType::tet10:
        {
          static VolMortarIntegrator<Core::FE::CellType::hex20, Core::FE::CellType::tet10>
              integrator(params());
          integrator.initialize_gp(true, domain);
          integrator.integrate_ele_3d(domain, source_ele, target_ele, *d1_, *m12_, *d2_, *m21_,
              *dis1_, *dis2_, dofset12_.first, dofset12_.second, dofset21_.first, dofset21_.second);
          break;
        }
        case Core::FE::CellType::pyramid5:
        {
          static VolMortarIntegrator<Core::FE::CellType::hex20, Core::FE::CellType::pyramid5>
              integrator(params());
          integrator.initialize_gp(true, domain);
          integrator.integrate_ele_3d(domain, source_ele, target_ele, *d1_, *m12_, *d2_, *m21_,
              *dis1_, *dis2_, dofset12_.first, dofset12_.second, dofset21_.first, dofset21_.second);
          break;
        }
        default:
        {
          FOUR_C_THROW("ERROR: unknown shape!");
          break;
        }
      }
      break;
    }
    // ########################################################
    case Core::FE::CellType::hex27:
    {
      switch (target_ele.shape())
      {
        // 2D surface elements
        case Core::FE::CellType::hex8:
        {
          static VolMortarIntegrator<Core::FE::CellType::hex27, Core::FE::CellType::hex8>
              integrator(params());
          integrator.initialize_gp(true, domain);
          integrator.integrate_ele_3d(domain, source_ele, target_ele, *d1_, *m12_, *d2_, *m21_,
              *dis1_, *dis2_, dofset12_.first, dofset12_.second, dofset21_.first, dofset21_.second);
          break;
        }
        case Core::FE::CellType::hex20:
        {
          static VolMortarIntegrator<Core::FE::CellType::hex27, Core::FE::CellType::hex20>
              integrator(params());
          integrator.initialize_gp(true, domain);
          integrator.integrate_ele_3d(domain, source_ele, target_ele, *d1_, *m12_, *d2_, *m21_,
              *dis1_, *dis2_, dofset12_.first, dofset12_.second, dofset21_.first, dofset21_.second);
          break;
        }
        case Core::FE::CellType::hex27:
        {
          static VolMortarIntegrator<Core::FE::CellType::hex27, Core::FE::CellType::hex27>
              integrator(params());
          integrator.initialize_gp(true, domain);
          integrator.integrate_ele_3d(domain, source_ele, target_ele, *d1_, *m12_, *d2_, *m21_,
              *dis1_, *dis2_, dofset12_.first, dofset12_.second, dofset21_.first, dofset21_.second);
          break;
        }
        case Core::FE::CellType::tet4:
        {
          static VolMortarIntegrator<Core::FE::CellType::hex27, Core::FE::CellType::tet4>
              integrator(params());
          integrator.initialize_gp(true, domain);
          integrator.integrate_ele_3d(domain, source_ele, target_ele, *d1_, *m12_, *d2_, *m21_,
              *dis1_, *dis2_, dofset12_.first, dofset12_.second, dofset21_.first, dofset21_.second);
          break;
        }
        case Core::FE::CellType::tet10:
        {
          static VolMortarIntegrator<Core::FE::CellType::hex27, Core::FE::CellType::tet10>
              integrator(params());
          integrator.initialize_gp(true, domain);
          integrator.integrate_ele_3d(domain, source_ele, target_ele, *d1_, *m12_, *d2_, *m21_,
              *dis1_, *dis2_, dofset12_.first, dofset12_.second, dofset21_.first, dofset21_.second);
          break;
        }
        case Core::FE::CellType::pyramid5:
        {
          static VolMortarIntegrator<Core::FE::CellType::hex27, Core::FE::CellType::pyramid5>
              integrator(params());
          integrator.initialize_gp(true, domain);
          integrator.integrate_ele_3d(domain, source_ele, target_ele, *d1_, *m12_, *d2_, *m21_,
              *dis1_, *dis2_, dofset12_.first, dofset12_.second, dofset21_.first, dofset21_.second);
          break;
        }
        default:
        {
          FOUR_C_THROW("ERROR: unknown shape!");
          break;
        }
      }
      break;
    }
    // ########################################################
    case Core::FE::CellType::tet4:
    {
      switch (target_ele.shape())
      {
        // 2D surface elements
        case Core::FE::CellType::hex8:
        {
          static VolMortarIntegrator<Core::FE::CellType::tet4, Core::FE::CellType::hex8> integrator(
              params());
          integrator.initialize_gp(true, domain);
          integrator.integrate_ele_3d(domain, source_ele, target_ele, *d1_, *m12_, *d2_, *m21_,
              *dis1_, *dis2_, dofset12_.first, dofset12_.second, dofset21_.first, dofset21_.second);
          break;
        }
        case Core::FE::CellType::hex20:
        {
          static VolMortarIntegrator<Core::FE::CellType::tet4, Core::FE::CellType::hex20>
              integrator(params());
          integrator.initialize_gp(true, domain);
          integrator.integrate_ele_3d(domain, source_ele, target_ele, *d1_, *m12_, *d2_, *m21_,
              *dis1_, *dis2_, dofset12_.first, dofset12_.second, dofset21_.first, dofset21_.second);
          break;
        }
        case Core::FE::CellType::hex27:
        {
          static VolMortarIntegrator<Core::FE::CellType::tet4, Core::FE::CellType::hex27>
              integrator(params());
          integrator.initialize_gp(true, domain);
          integrator.integrate_ele_3d(domain, source_ele, target_ele, *d1_, *m12_, *d2_, *m21_,
              *dis1_, *dis2_, dofset12_.first, dofset12_.second, dofset21_.first, dofset21_.second);
          break;
        }
        case Core::FE::CellType::tet4:
        {
          static VolMortarIntegrator<Core::FE::CellType::tet4, Core::FE::CellType::tet4> integrator(
              params());
          integrator.initialize_gp(true, domain);
          integrator.integrate_ele_3d(domain, source_ele, target_ele, *d1_, *m12_, *d2_, *m21_,
              *dis1_, *dis2_, dofset12_.first, dofset12_.second, dofset21_.first, dofset21_.second);
          break;
        }
        case Core::FE::CellType::tet10:
        {
          static VolMortarIntegrator<Core::FE::CellType::tet4, Core::FE::CellType::tet10>
              integrator(params());
          integrator.initialize_gp(true, domain);
          integrator.integrate_ele_3d(domain, source_ele, target_ele, *d1_, *m12_, *d2_, *m21_,
              *dis1_, *dis2_, dofset12_.first, dofset12_.second, dofset21_.first, dofset21_.second);
          break;
        }
        case Core::FE::CellType::pyramid5:
        {
          static VolMortarIntegrator<Core::FE::CellType::tet4, Core::FE::CellType::pyramid5>
              integrator(params());
          integrator.initialize_gp(true, domain);
          integrator.integrate_ele_3d(domain, source_ele, target_ele, *d1_, *m12_, *d2_, *m21_,
              *dis1_, *dis2_, dofset12_.first, dofset12_.second, dofset21_.first, dofset21_.second);
          break;
        }
        default:
        {
          FOUR_C_THROW("ERROR: unknown shape!");
          break;
        }
      }
      break;
    }
    // ########################################################
    case Core::FE::CellType::tet10:
    {
      switch (target_ele.shape())
      {
        // 2D surface elements
        case Core::FE::CellType::hex8:
        {
          static VolMortarIntegrator<Core::FE::CellType::tet10, Core::FE::CellType::hex8>
              integrator(params());
          integrator.initialize_gp(true, domain);
          integrator.integrate_ele_3d(domain, source_ele, target_ele, *d1_, *m12_, *d2_, *m21_,
              *dis1_, *dis2_, dofset12_.first, dofset12_.second, dofset21_.first, dofset21_.second);
          break;
        }
        case Core::FE::CellType::hex20:
        {
          static VolMortarIntegrator<Core::FE::CellType::tet10, Core::FE::CellType::hex20>
              integrator(params());
          integrator.initialize_gp(true, domain);
          integrator.integrate_ele_3d(domain, source_ele, target_ele, *d1_, *m12_, *d2_, *m21_,
              *dis1_, *dis2_, dofset12_.first, dofset12_.second, dofset21_.first, dofset21_.second);
          break;
        }
        case Core::FE::CellType::hex27:
        {
          static VolMortarIntegrator<Core::FE::CellType::tet10, Core::FE::CellType::hex27>
              integrator(params());
          integrator.initialize_gp(true, domain);
          integrator.integrate_ele_3d(domain, source_ele, target_ele, *d1_, *m12_, *d2_, *m21_,
              *dis1_, *dis2_, dofset12_.first, dofset12_.second, dofset21_.first, dofset21_.second);
          break;
        }
        case Core::FE::CellType::tet4:
        {
          static VolMortarIntegrator<Core::FE::CellType::tet10, Core::FE::CellType::tet4>
              integrator(params());
          integrator.initialize_gp(true, domain);
          integrator.integrate_ele_3d(domain, source_ele, target_ele, *d1_, *m12_, *d2_, *m21_,
              *dis1_, *dis2_, dofset12_.first, dofset12_.second, dofset21_.first, dofset21_.second);
          break;
        }
        case Core::FE::CellType::tet10:
        {
          static VolMortarIntegrator<Core::FE::CellType::tet10, Core::FE::CellType::tet10>
              integrator(params());
          integrator.initialize_gp(true, domain);
          integrator.integrate_ele_3d(domain, source_ele, target_ele, *d1_, *m12_, *d2_, *m21_,
              *dis1_, *dis2_, dofset12_.first, dofset12_.second, dofset21_.first, dofset21_.second);
          break;
        }
        case Core::FE::CellType::pyramid5:
        {
          static VolMortarIntegrator<Core::FE::CellType::tet10, Core::FE::CellType::pyramid5>
              integrator(params());
          integrator.initialize_gp(true, domain);
          integrator.integrate_ele_3d(domain, source_ele, target_ele, *d1_, *m12_, *d2_, *m21_,
              *dis1_, *dis2_, dofset12_.first, dofset12_.second, dofset21_.first, dofset21_.second);
          break;
        }
        default:
        {
          FOUR_C_THROW("ERROR: unknown shape!");
          break;
        }
      }
      break;
    }
    // ########################################################
    case Core::FE::CellType::pyramid5:
    {
      switch (target_ele.shape())
      {
        // 2D surface elements
        case Core::FE::CellType::hex8:
        {
          static VolMortarIntegrator<Core::FE::CellType::pyramid5, Core::FE::CellType::hex8>
              integrator(params());
          integrator.initialize_gp(true, domain);
          integrator.integrate_ele_3d(domain, source_ele, target_ele, *d1_, *m12_, *d2_, *m21_,
              *dis1_, *dis2_, dofset12_.first, dofset12_.second, dofset21_.first, dofset21_.second);
          break;
        }
        case Core::FE::CellType::hex20:
        {
          static VolMortarIntegrator<Core::FE::CellType::pyramid5, Core::FE::CellType::hex20>
              integrator(params());
          integrator.initialize_gp(true, domain);
          integrator.integrate_ele_3d(domain, source_ele, target_ele, *d1_, *m12_, *d2_, *m21_,
              *dis1_, *dis2_, dofset12_.first, dofset12_.second, dofset21_.first, dofset21_.second);
          break;
        }
        case Core::FE::CellType::hex27:
        {
          static VolMortarIntegrator<Core::FE::CellType::pyramid5, Core::FE::CellType::hex27>
              integrator(params());
          integrator.initialize_gp(true, domain);
          integrator.integrate_ele_3d(domain, source_ele, target_ele, *d1_, *m12_, *d2_, *m21_,
              *dis1_, *dis2_, dofset12_.first, dofset12_.second, dofset21_.first, dofset21_.second);
          break;
        }
        case Core::FE::CellType::tet4:
        {
          static VolMortarIntegrator<Core::FE::CellType::pyramid5, Core::FE::CellType::tet4>
              integrator(params());
          integrator.initialize_gp(true, domain);
          integrator.integrate_ele_3d(domain, source_ele, target_ele, *d1_, *m12_, *d2_, *m21_,
              *dis1_, *dis2_, dofset12_.first, dofset12_.second, dofset21_.first, dofset21_.second);
          break;
        }
        case Core::FE::CellType::tet10:
        {
          static VolMortarIntegrator<Core::FE::CellType::pyramid5, Core::FE::CellType::tet10>
              integrator(params());
          integrator.initialize_gp(true, domain);
          integrator.integrate_ele_3d(domain, source_ele, target_ele, *d1_, *m12_, *d2_, *m21_,
              *dis1_, *dis2_, dofset12_.first, dofset12_.second, dofset21_.first, dofset21_.second);
          break;
        }
        case Core::FE::CellType::pyramid5:
        {
          static VolMortarIntegrator<Core::FE::CellType::pyramid5, Core::FE::CellType::pyramid5>
              integrator(params());
          integrator.initialize_gp(true, domain);
          integrator.integrate_ele_3d(domain, source_ele, target_ele, *d1_, *m12_, *d2_, *m21_,
              *dis1_, *dis2_, dofset12_.first, dofset12_.second, dofset21_.first, dofset21_.second);
          break;
        }
        default:
        {
          FOUR_C_THROW("ERROR: unknown shape!");
          break;
        }
      }
      break;
    }
    default:
    {
      FOUR_C_THROW("ERROR: unknown shape!");
      break;
    }
  }

  // integration element counter
  inteles_++;
}

/*----------------------------------------------------------------------*
 |  init (public)                                            farah 10/13|
 *----------------------------------------------------------------------*/
void Coupling::VolMortar::VolMortarCoupl::initialize()
{
  /* ******************************************************************
   * (re)setup global Mortar Core::LinAlg::SparseMatrices                   *
   * unknowns which are going to be condensed are defined on the source*
   * side. Therefore, the rows are the auxiliary variables on the     *
   * source side!                                                      *
   * ******************************************************************/

  d1_ = std::make_shared<Core::LinAlg::SparseMatrix>(*p12_dofrowmap_, 10);
  m12_ = std::make_shared<Core::LinAlg::SparseMatrix>(*p12_dofrowmap_, 100);

  d2_ = std::make_shared<Core::LinAlg::SparseMatrix>(*p21_dofrowmap_, 10);
  m21_ = std::make_shared<Core::LinAlg::SparseMatrix>(*p21_dofrowmap_, 100);

  // initialize trafo operator for quadr. modification
  if (dualquad_ != dualquad_no_mod)
  {
    t1_ = std::make_shared<Core::LinAlg::SparseMatrix>(*p12_dofrowmap_, 10);
    t2_ = std::make_shared<Core::LinAlg::SparseMatrix>(*p21_dofrowmap_, 10);
  }
}

/*----------------------------------------------------------------------*
 |  Complete (public)                                        farah 01/14|
 *----------------------------------------------------------------------*/
void Coupling::VolMortar::VolMortarCoupl::complete()
{
  // complete...
  d1_->complete(*p12_dofrowmap_, *p12_dofrowmap_);
  m12_->complete(*p12_dofdomainmap_, *p12_dofrowmap_);

  d2_->complete(*p21_dofrowmap_, *p21_dofrowmap_);
  m21_->complete(*p21_dofdomainmap_, *p21_dofrowmap_);

  // complete trafo operator for quadr. modification
  if (dualquad_ != dualquad_no_mod)
  {
    t1_->complete(*p12_dofrowmap_, *p12_dofrowmap_);
    t2_->complete(*p21_dofrowmap_, *p21_dofrowmap_);
  }
}

/*----------------------------------------------------------------------*
 |  compute projection operator P                            farah 01/14|
 *----------------------------------------------------------------------*/
void Coupling::VolMortar::VolMortarCoupl::create_projection_operator()
{
  /********************************************************************/
  /* Multiply Mortar matrices: P = inv(D) * M         A               */
  /********************************************************************/
  Core::LinAlg::SparseMatrix invd1(*d1_);
  std::shared_ptr<Core::LinAlg::Vector<double>> diag1 =
      std::make_shared<Core::LinAlg::Vector<double>>(*p12_dofrowmap_, true);

  // extract diagonal of invd into diag
  invd1.extract_diagonal_copy(*diag1);

  // set zero diagonal values to dummy 1.0
  for (int i = 0; i < diag1->local_length(); ++i)
    if (abs(diag1->local_values_as_span()[i]) < 1e-12) (*diag1).get_values()[i] = 1.0;

  // scalar inversion of diagonal values
  diag1->reciprocal(*diag1);

  // re-insert inverted diagonal into invd
  invd1.replace_diagonal_values(*diag1);

  // do the multiplication P = inv(D) * M
  std::shared_ptr<Core::LinAlg::SparseMatrix> aux12 =
      Core::LinAlg::matrix_multiply(invd1, false, *m12_, false, false, false, true);

  /********************************************************************/
  /* Multiply Mortar matrices: P = inv(D) * M         B               */
  /********************************************************************/
  Core::LinAlg::SparseMatrix invd2(*d2_);
  std::shared_ptr<Core::LinAlg::Vector<double>> diag2 =
      std::make_shared<Core::LinAlg::Vector<double>>(*p21_dofrowmap_, true);

  // extract diagonal of invd into diag
  invd2.extract_diagonal_copy(*diag2);

  // set zero diagonal values to dummy 1.0
  for (int i = 0; i < diag2->local_length(); ++i)
    if (abs(diag2->local_values_as_span()[i]) < 1e-12) (*diag2).get_values()[i] = 1.0;

  // scalar inversion of diagonal values
  diag2->reciprocal(*diag2);

  // re-insert inverted diagonal into invd
  invd2.replace_diagonal_values(*diag2);

  // do the multiplication P = inv(D) * M
  std::shared_ptr<Core::LinAlg::SparseMatrix> aux21 =
      Core::LinAlg::matrix_multiply(invd2, false, *m21_, false, false, false, true);

  // initialize trafo operator for quadr. modification
  if (dualquad_ != dualquad_no_mod)
  {
    p12_ = Core::LinAlg::matrix_multiply(*t1_, false, *aux12, false, false, false, true);
    p21_ = Core::LinAlg::matrix_multiply(*t2_, false, *aux21, false, false, false, true);
  }
  else
  {
    p12_ = aux12;
    p21_ = aux21;
  }
}

/*----------------------------------------------------------------------*
 |  Define polygon of mortar vertices                        farah 01/14|
 *----------------------------------------------------------------------*/
void Coupling::VolMortar::VolMortarCoupl::define_vertices_source(
    Core::Elements::Element& ele, std::vector<Mortar::Vertex>& source_vertices)
{
  // project source nodes onto auxiliary plane
  int nnodes = ele.num_node();
  Core::Nodes::Node** mynodes = ele.nodes();
  if (!mynodes) FOUR_C_THROW("ERROR: project_source: Null pointer!");

  // initialize storage for source coords + their ids
  std::vector<double> vertices(3);
  std::vector<int> snodeids(1);

  for (int i = 0; i < nnodes; ++i)
  {
    // compute projection
    const auto x = mynodes[i]->x();
    for (size_t k = 0; k < x.size(); ++k) vertices[k] = x[k];

    // get node id, too
    snodeids[0] = mynodes[i]->id();

    // store into vertex data structure
    source_vertices.emplace_back(vertices, Mortar::Vertex::source, snodeids, nullptr);
  }
}

/*----------------------------------------------------------------------*
 |  Define polygon of mortar vertices                        farah 01/14|
 *----------------------------------------------------------------------*/
void Coupling::VolMortar::VolMortarCoupl::define_vertices_target(
    Core::Elements::Element& ele, std::vector<Mortar::Vertex>& target_vertices)
{
  // project target nodes onto auxiliary plane
  int nnodes = ele.num_node();
  Core::Nodes::Node** mynodes = ele.nodes();
  if (!mynodes) FOUR_C_THROW("ERROR: project_target: Null pointer!");

  // initialize storage for target coords + their ids
  std::vector<double> vertices(3);
  std::vector<int> tnodeids(1);

  for (int i = 0; i < nnodes; ++i)
  {
    // compute projection
    const auto x = mynodes[i]->x();
    for (size_t k = 0; k < x.size(); ++k) vertices[k] = x[k];

    // get node id, too
    tnodeids[0] = mynodes[i]->id();

    // store into vertex data structure
    target_vertices.emplace_back(vertices, Mortar::Vertex::projtarget, tnodeids, nullptr);
  }
}

/*----------------------------------------------------------------------*
 |  Triangulation of clip polygon (3D) - DELAUNAY             popp 08/11|
 *----------------------------------------------------------------------*/
void Coupling::VolMortar::VolMortarCoupl::delaunay_triangulation(
    std::vector<std::shared_ptr<Mortar::IntCell>>& cells, std::vector<Mortar::Vertex>& clip,
    double tol)
{
  cells.resize(0);

  // no linearization needed for volume mortar
  std::vector<Core::Gen::Pairedvector<int, double>> linvertex(3, 100);
  std::vector<Core::Gen::Pairedvector<int, double>> derivauxn;

  const Mortar::Utils::PlaneCoordinateSystem plane(clip, auxn());
  const auto triangles =
      Core::Geo::delaunay_triangulation_convex_polygon_2d(plane.to_plane(clip), tol);

  for (int t = 0; t < static_cast<int>(triangles.size()); ++t)
  {
    Core::LinAlg::Matrix<3, 3> coords;
    for (int i = 0; i < 3; ++i)
      for (int k = 0; k < 3; ++k) coords(k, i) = clip[triangles[t][i]].coord()[k];

    cells.push_back(std::make_shared<Mortar::IntCell>(t, 3, coords, auxn(),
        Core::FE::CellType::tri3, linvertex, linvertex, linvertex, derivauxn));
  }
}

FOUR_C_NAMESPACE_CLOSE
