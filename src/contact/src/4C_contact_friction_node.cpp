// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_contact_friction_node.hpp"

#include "4C_comm_pack_helpers.hpp"
#include "4C_contact_defines.hpp"
#include "4C_contact_element.hpp"
#include "4C_utils_exceptions.hpp"

FOUR_C_NAMESPACE_OPEN

CONTACT::FriNodeType CONTACT::FriNodeType::instance_;

Core::Communication::ParObject* CONTACT::FriNodeType::create(
    Core::Communication::UnpackBuffer& buffer)
{
  std::vector<double> x(3, 0.0);
  std::vector<int> dofs;

  // TODO: friplus = true for all nodes!!! change this with pack/unpack
  CONTACT::FriNode* node = new CONTACT::FriNode(0, x, 0, dofs, false, false, true);
  node->unpack(buffer);

  return node;
}

/*----------------------------------------------------------------------*
 |  ctor (public)                                             mgit 01/10|
 *----------------------------------------------------------------------*/
CONTACT::FriNodeDataContainer::FriNodeDataContainer()
    : slip_(false), slipold_(false), drowsold_(0), drowsoldLTL_(0)
{
  for (int i = 0; i < 3; ++i)
  {
    jump()[i] = 0.0;
    traction()[i] = 0.0;
    tractionold()[i] = 0.0;
    traction_ltl()[i] = 0.0;
    tractionold_ltl()[i] = 0.0;
  }
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                            mgit 01/10|
 *----------------------------------------------------------------------*/
void CONTACT::FriNodeDataContainer::pack(Core::Communication::PackBuffer& data) const
{
  // add jump_
  add_to_pack(data, jump_, 3 * sizeof(double));
  // add slip_
  add_to_pack(data, slip_);
  // add slip_
  add_to_pack(data, slipold_);
  // add traction_
  add_to_pack(data, traction_, 3 * sizeof(double));
  // add tractionold_
  add_to_pack(data, tractionold_, 3 * sizeof(double));
  // add traction_
  add_to_pack(data, tractionLTL_, 3 * sizeof(double));
  // add tractionold_
  add_to_pack(data, tractionoldLTL_, 3 * sizeof(double));

  // add drowsold_,mrowsold_,mnodesold_
  int hasdata = drowsold_.size();

  add_to_pack(data, hasdata);

  if (hasdata != 0)
  {
    int dentries = (int)drowsold_.size();
    add_to_pack(data, dentries);
    add_to_pack(data, drowsold_);
    add_to_pack(data, mrowsold_);
    add_to_pack(data, mnodesold_);
  }

  int hasdata2 = drowsoldLTL_.size();
  add_to_pack(data, hasdata2);
  if (hasdata2 != 0)
  {
    int dentries = (int)drowsoldLTL_.size();
    add_to_pack(data, dentries);
    add_to_pack(data, drowsoldLTL_);
    add_to_pack(data, mrowsoldLTL_);
  }
  // add derivjump
  int hasdataderivjump = derivjump_.size();
  add_to_pack(data, hasdataderivjump);

  if (hasdataderivjump != 0)
  {
    for (int i = 0; i < hasdataderivjump; i++) add_to_pack(data, (derivjump_[i]));
  }
}

/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                            mgit 01/10|
 *----------------------------------------------------------------------*/
void CONTACT::FriNodeDataContainer::unpack(Core::Communication::UnpackBuffer& buffer)
{
  // jump_
  extract_from_pack(buffer, jump_, 3 * sizeof(double));
  // slip_
  extract_from_pack(buffer, slip_);
  // slipold_
  extract_from_pack(buffer, slipold_);
  // traction_
  extract_from_pack(buffer, traction_, 3 * sizeof(double));
  // tractionold_
  extract_from_pack(buffer, tractionold_, 3 * sizeof(double));
  // traction_
  extract_from_pack(buffer, tractionLTL_, 3 * sizeof(double));
  // tractionold_
  extract_from_pack(buffer, tractionoldLTL_, 3 * sizeof(double));

  // drowsold_,mrowsold_,mnodesold_
  int hasdata;
  extract_from_pack(buffer, hasdata);

  if (hasdata != 0)
  {
    int dentries;
    extract_from_pack(buffer, dentries);

    drowsold_.resize(dentries);
    extract_from_pack(buffer, drowsold_);
    extract_from_pack(buffer, mrowsold_);
    extract_from_pack(buffer, mnodesold_);
  }

  // drowsold_,mrowsold_,mnodesold_
  int hasdata2;
  extract_from_pack(buffer, hasdata2);

  if (hasdata2 != 0)
  {
    int dentries;
    extract_from_pack(buffer, dentries);

    drowsoldLTL_.resize(dentries);
    extract_from_pack(buffer, drowsoldLTL_);
    extract_from_pack(buffer, mrowsoldLTL_);
  }

  // and derivjump_
  int hasdataderivjump;
  extract_from_pack(buffer, hasdataderivjump);

  if (hasdataderivjump != 0)
  {
    derivjump_.resize(hasdataderivjump);
    for (int i = 0; i < hasdataderivjump; i++)
    {
      extract_from_pack(buffer, derivjump_[i]);
    }
  }
}

/*----------------------------------------------------------------------*
 |  ctor (public)                                             mgit 07/11|
 *----------------------------------------------------------------------*/
CONTACT::FriNodeWearDataContainer::FriNodeWearDataContainer()
    : weightedwear_(0.0), deltaweightedwear_(0.0)
{
  wcurr_[0] = 0.0;
  wold_[0] = 0.0;
  waccu_[0] = 0.0;
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                           farah 10/13|
 *----------------------------------------------------------------------*/
void CONTACT::FriNodeWearDataContainer::pack(Core::Communication::PackBuffer& data) const
{
  add_to_pack(data, weightedwear_);
  add_to_pack(data, deltaweightedwear_);

  // add d2row
  int hasdata = d2rows_.size();
  add_to_pack(data, hasdata);

  if (hasdata != 0)
  {
    for (int i = 0; i < hasdata; i++)
    {
      add_to_pack(data, (d2rows_[i]));
    }
  }
}

/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                           farah 10/13|
 *----------------------------------------------------------------------*/
void CONTACT::FriNodeWearDataContainer::unpack(Core::Communication::UnpackBuffer& buffer)
{
  extract_from_pack(buffer, weightedwear_);
  extract_from_pack(buffer, deltaweightedwear_);

  // d2rows_
  int hasdata;
  extract_from_pack(buffer, hasdata);

  if (hasdata != 0)
  {
    d2rows_.resize(hasdata);
    for (int i = 0; i < hasdata; i++)
    {
      extract_from_pack(buffer, d2rows_[i]);
    }
  }
}

/*----------------------------------------------------------------------*
 |  ctor (public)                                             mgit 02/10|
 *----------------------------------------------------------------------*/
CONTACT::FriNode::FriNode(int id, std::span<const double> coords, const int owner,
    const std::vector<int>& dofs, const bool isslave, const bool initactive, const bool friplus)
    : CONTACT::Node(id, coords, owner, dofs, isslave, initactive), wear_(friplus)
{
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                        mgit 02/10|
 *----------------------------------------------------------------------*/
CONTACT::FriNode::FriNode(const CONTACT::FriNode& old) : CONTACT::Node(old), wear_(old.wear_)
{
  // not yet used and thus not necessarily consistent
  FOUR_C_THROW("FriNode copy-ctor not yet implemented");
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance of FriNode and return pointer to it (public)|
 |                                                            mgit 02/10|
 *----------------------------------------------------------------------*/
CONTACT::FriNode* CONTACT::FriNode::clone() const
{
  CONTACT::FriNode* newnode = new CONTACT::FriNode(*this);
  return newnode;
}

/*----------------------------------------------------------------------*
 |  << operator                                               mgit 02/10|
 *----------------------------------------------------------------------*/
std::ostream& operator<<(std::ostream& os, const CONTACT::FriNode& frinode)
{
  frinode.print(os);
  return os;
}

/*----------------------------------------------------------------------*
 |  print this element (public)                               mgit 02/10|
 *----------------------------------------------------------------------*/
void CONTACT::FriNode::print(std::ostream& os) const
{
  // Print id and coordinates
  os << "Contact ";
  CONTACT::Node::print(os);
  if (is_slave())
    if (is_init_active()) os << " InitActive ";
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                            mgit 02/10|
 *----------------------------------------------------------------------*/
void CONTACT::FriNode::pack(Core::Communication::PackBuffer& data) const
{
  // pack type of this instance of ParObject
  int type = unique_par_object_id();
  add_to_pack(data, type);

  // add base class Mortar::Node
  CONTACT::Node::pack(data);

  // add data_
  bool hasdata = (fridata_ != nullptr);
  add_to_pack(data, hasdata);
  if (hasdata) fridata_->pack(data);

  bool hasweardata = (weardata_ != nullptr);
  add_to_pack(data, hasweardata);
  if (hasweardata) weardata_->pack(data);
}

/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                            mgit 02/10|
 *----------------------------------------------------------------------*/
void CONTACT::FriNode::unpack(Core::Communication::UnpackBuffer& buffer)
{
  Core::Communication::extract_and_assert_id(buffer, unique_par_object_id());

  // extract base class CONTACT::Node
  CONTACT::Node::unpack(buffer);

  // **************************
  // FriData
  bool hasdata;
  extract_from_pack(buffer, hasdata);
  if (hasdata)
  {
    fridata_ = std::make_shared<CONTACT::FriNodeDataContainer>();
    fridata_->unpack(buffer);
  }
  else
    fridata_ = nullptr;

  // **************************
  // FriDataPlus
  bool hasdataplus;
  extract_from_pack(buffer, hasdataplus);
  if (hasdataplus)
  {
    weardata_ = std::make_shared<CONTACT::FriNodeWearDataContainer>();
    weardata_->unpack(buffer);
  }
  else
    weardata_ = nullptr;

  // Check
}

/*----------------------------------------------------------------------*
 | calculate the apparent coefficient of friction            seitz 11/15|
 *----------------------------------------------------------------------*/
double CONTACT::FriNode::fr_coeff(const double frcoeff_in) const
{
  // return the friction coefficient, if we do not have a TSI problem
  if (cTSIdata_ == nullptr) return frcoeff_in;

  // in TSI case, the friction coefficient is temperature dependent
  else
  {
    double maxT = std::max(tsi_data().temp(), tsi_data().temp_master());
    return frcoeff_in * (maxT - cTSIdata_->temp_dam()) * (maxT - cTSIdata_->temp_dam()) /
           ((cTSIdata_->temp_dam() - cTSIdata_->temp_ref()) *
               (cTSIdata_->temp_dam() - cTSIdata_->temp_ref()));
  }
}

/*----------------------------------------------------------------------*
 | calculate derivative of apparent coefficient of friction  seitz 11/15|
 *----------------------------------------------------------------------*/
void CONTACT::FriNode::deriv_fr_coeff_temp(
    const double frcoeff_in, std::map<int, double>& derivT, std::map<int, double>& derivDisp) const
{
  derivT.clear();
  derivDisp.clear();

  // if we do not have a TSI problem, the friction coefficient is constant
  if (cTSIdata_ == nullptr) return;

  double T_dam = cTSIdata_->temp_dam();
  double T_ref = cTSIdata_->temp_ref();
  if (cTSIdata_->temp() > cTSIdata_->temp_master())
  {
    double maxT = tsi_data().temp();
    derivT[dofs()[0]] += 2. * frcoeff_in * (maxT - T_dam) / ((T_dam - T_ref) * (T_dam - T_ref));
    derivDisp.clear();
  }
  else
  {
    double maxT = tsi_data().temp_master();
    for (std::map<int, double>::const_iterator i = tsi_data().deriv_temp_master_temp().begin();
        i != tsi_data().deriv_temp_master_temp().end(); ++i)
      derivT[i->first] +=
          2. * frcoeff_in * (maxT - T_dam) / ((T_dam - T_ref) * (T_dam - T_ref)) * i->second;
    for (std::map<int, double>::const_iterator i = tsi_data().deriv_temp_master_disp().begin();
        i != tsi_data().deriv_temp_master_disp().end(); ++i)
      derivDisp[i->first] +=
          2. * frcoeff_in * (maxT - T_dam) / ((T_dam - T_ref) * (T_dam - T_ref)) * i->second;
  }
}

/*----------------------------------------------------------------------*
 |  Add a value to the 'SNodes' set                        gitterle 11/09|
 *----------------------------------------------------------------------*/
void CONTACT::FriNode::add_s_node(int node) { fri_data().get_s_nodes().insert(node); }

/*----------------------------------------------------------------------*
 |  Add a value to the 'MNodes' set                        gitterle 11/09|
 *----------------------------------------------------------------------*/
void CONTACT::FriNode::add_m_node(int node) { fri_data().get_m_nodes().insert(node); }

/*----------------------------------------------------------------------*
 |  Add a value to the 'D2' map                              farah 06/13|
 *----------------------------------------------------------------------*/
void CONTACT::FriNode::add_d2_value(int row, int col, double val)
{
  // check if this is a master node or slave boundary node
  if (is_slave() == true) FOUR_C_THROW("AddD2Value: function called for slave node {}", id());

  // check if this has been called before
  if ((int)wear_data().get_d2().size() == 0) wear_data().get_d2().resize(num_dof());

  // check row index input
  if ((int)wear_data().get_d2().size() <= row)
    FOUR_C_THROW("AddD2Value: tried to access invalid row index!");

  // add the pair (col,val) to the given row
  std::map<int, double>& d2map = wear_data().get_d2()[row];
  d2map[col] += val;
}

/*----------------------------------------------------------------------*
 |  Add a value to the 'DerivJump' map                     gitterle 11/09|
 *----------------------------------------------------------------------*/
void CONTACT::FriNode::add_deriv_jump_value(int row, int col, double val)
{
  // check if this is a master node or slave boundary node
  if (is_slave() == false) FOUR_C_THROW("AddJumpValue: function called for master node {}", id());
  if (is_on_bound() == true)
    FOUR_C_THROW("AddJumpValue: function called for boundary node {}", id());

  // check if this has been called before
  if ((int)fri_data().get_deriv_jump().size() == 0) fri_data().get_deriv_jump().resize(num_dof());

  // check row index input
  if ((int)fri_data().get_deriv_jump().size() <= row)
    FOUR_C_THROW("AddDerivJumpValue: tried to access invalid row index!");

  // add the pair (col,val) to the given row
  std::map<int, double>& zmap = fri_data().get_deriv_jump()[row];
  zmap[col] += val;
}

/*----------------------------------------------------------------------*
|  Add jump value from gp-wise integration of slip          farah 08/13|
*----------------------------------------------------------------------*/
void CONTACT::FriNode::add_jump_value(double val, int k)
{
  // check if this is a master node or slave boundary node
  if (is_slave() == false) FOUR_C_THROW("AddJumpValue: function called for master node {}", id());
  if (is_on_bound() == true)
    FOUR_C_THROW("AddJumpValue: function called for boundary node {}", id());

  fri_data().jump_var()[k] += val;
}

/*----------------------------------------------------------------------*
 |  Add a value to the 'T' map                               farah 09/13|
 *----------------------------------------------------------------------*/
void CONTACT::FriNode::add_t_value(int row, int col, double val)
{
  // check if this is a master node or slave boundary node
  //  if (IsSlave()==false)
  //     FOUR_C_THROW("AddTValue: function called for master node {}", Id());

  // check if this has been called before
  if ((int)wear_data().get_t().size() == 0) wear_data().get_t().resize(num_dof());

  // check row index input
  if ((int)wear_data().get_t().size() <= row)
    FOUR_C_THROW("AddTValue: tried to access invalid row index!");

  // add the pair (col,val) to the given row
  std::map<int, double>& tmap = wear_data().get_t()[row];
  tmap[col] += val;
}

/*----------------------------------------------------------------------*
 |  Add a value to the 'E' map                               farah 09/13|
 *----------------------------------------------------------------------*/
void CONTACT::FriNode::add_e_value(int row, int col, double val)
{
  // check if this is a master node or slave boundary node
  //  if (IsSlave()==false)
  //     FOUR_C_THROW("AddEValue: function called for master node {}", Id());

  // check if this has been called before
  if ((int)wear_data().get_e().size() == 0) wear_data().get_e().resize(num_dof());

  // check row index input
  if ((int)wear_data().get_e().size() <= row)
    FOUR_C_THROW("AddEValue: tried to access invalid row index!");

  // add the pair (col,val) to the given row
  std::map<int, double>& emap = wear_data().get_e()[row];
  emap[col] += val;
}

/*----------------------------------------------------------------------*
 |  Store nodal entries of D and M to old ones             gitterle 12/08|
 *----------------------------------------------------------------------*/
void CONTACT::FriNode::store_dm_old()
{
  // clear and zero nodal vectors
  fri_data().get_d_old().clear();
  fri_data().get_m_old().clear();
  fri_data().get_d_old_ltl().clear();
  fri_data().get_m_old_ltl().clear();

  // write drows_ to drowsold_
  fri_data().get_d_old() = mo_data().get_d();
  fri_data().get_m_old() = mo_data().get_m();

  // write drows_ to drowsold_
  fri_data().get_d_old_ltl() = mo_data().get_dltl();
  fri_data().get_m_old_ltl() = mo_data().get_mltl();

  // also vectors containing the according master nodes
  fri_data().get_m_nodes_old().clear();
  fri_data().get_m_nodes_old() = fri_data().get_m_nodes();
}

/*----------------------------------------------------------------------*
 |  Store nodal entries penalty tractions to old ones      gitterle 10/09|
 *----------------------------------------------------------------------*/
void CONTACT::FriNode::store_trac_old()
{
  // write entries to old ones
  for (int j = 0; j < 3; ++j) fri_data().tractionold()[j] = fri_data().traction()[j];

  for (int j = 0; j < 3; ++j) fri_data().tractionold_ltl()[j] = fri_data().traction_ltl()[j];
}

/*-----------------------------------------------------------------------*
 |  Set the value of deltawear                             gitterle 12/10|
 *----------------------------------------------------------------------*/
void CONTACT::FriNode::add_delta_weighted_wear_value(double val)
{
  // add given value to deltawear_
  wear_data().delta_weighted_wear() += val;
}

/*----------------------------------------------------------------------*
 |  Initialize data container                             gitterle 10/09|
 *----------------------------------------------------------------------*/
void CONTACT::FriNode::initialize_data_container()
{
  // get maximum size of lin vectors
  linsize_ = 0;
  for (int i = 0; i < num_element(); ++i)
    for (int j = 0; j < elements()[i]->num_node(); ++j)
      linsize_ += elements()[i]->num_dof_per_node(*(elements()[i]->nodes()[j]));

  // get maximum size of lin vectors
  dentries_ = 0;
  for (int i = 0; i < num_element(); ++i)
    for (int j = 0; j < elements()[i]->num_node(); ++j)
      dentries_ += elements()[i]->num_dof_per_node(*(elements()[i]->nodes()[j]));

  // only initialize if not yet done
  if (modata_ == nullptr && codata_ == nullptr && fridata_ == nullptr)
  {
    modata_ = std::make_shared<Mortar::NodeDataContainer>();
    codata_ = std::make_shared<CONTACT::NodeDataContainer>();
    fridata_ = std::make_shared<CONTACT::FriNodeDataContainer>();
  }

  // initialize data container for wear and tsi problems
  if (wear_ == true)
  {
    if (weardata_ == nullptr) weardata_ = std::make_shared<CONTACT::FriNodeWearDataContainer>();
  }
}

/*----------------------------------------------------------------------*
 |  Reset data container                                      popp 09/10|
 *----------------------------------------------------------------------*/
void CONTACT::FriNode::reset_data_container()
{
  // reset to nullptr
  fridata_ = nullptr;
  weardata_ = nullptr;
  codata_ = nullptr;
  modata_ = nullptr;
}

FOUR_C_NAMESPACE_CLOSE
