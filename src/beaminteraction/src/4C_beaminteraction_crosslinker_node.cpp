// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_beaminteraction_crosslinker_node.hpp"

#include "4C_comm_pack_helpers.hpp"
#include "4C_comm_utils_factory.hpp"
#include "4C_mat_crosslinkermat.hpp"
#include "4C_utils_exceptions.hpp"

FOUR_C_NAMESPACE_OPEN

CrossLinking::CrosslinkerNodeType CrossLinking::CrosslinkerNodeType::instance_;



/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Core::Communication::ParObject* CrossLinking::CrosslinkerNodeType::create(
    Core::Communication::UnpackBuffer& buffer)
{
  std::vector<double> dummycoord(3, 999.0);
  auto* crosslinker = new CrossLinking::CrosslinkerNode(-1, dummycoord, -1);
  crosslinker->unpack(buffer);
  return crosslinker;
}

/*----------------------------------------------------------------------------*
 |  ctor (public)                                              eichinger 10/16|
 *----------------------------------------------------------------------------*/
CrossLinking::CrosslinkerNodeDataContainer::CrosslinkerNodeDataContainer() : numbond_(0)
{
  clbspots_.clear();
  std::pair<int, int> pair;
  pair.first = -1;
  pair.second = -1;
  clbspots_.push_back(pair);  // first binding spot of crosslinker
  clbspots_.push_back(pair);  // second binding spot of crosslinker

  return;
}

/*----------------------------------------------------------------------------*
 |  Pack data                                                        (public) |
 |                                                             eichinger 10/16|
 *----------------------------------------------------------------------------*/
void CrossLinking::CrosslinkerNodeDataContainer::pack(Core::Communication::PackBuffer& data) const
{
  // add numbond
  add_to_pack(data, numbond_);
  // add clbspots_
  add_to_pack(data, clbspots_);

  return;
}

/*----------------------------------------------------------------------------*
 |  Unpack data                                                      (public) |
 |                                                             eichinger 10/16|
 *----------------------------------------------------------------------------*/
void CrossLinking::CrosslinkerNodeDataContainer::unpack(Core::Communication::UnpackBuffer& buffer)
{
  // numbond
  extract_from_pack(buffer, numbond_);
  // clbspots_
  extract_from_pack(buffer, clbspots_);

  return;
}

/*----------------------------------------------------------------------------*
 *  ctor (public)                                              eichinger 10/16|
 *----------------------------------------------------------------------------*/
CrossLinking::CrosslinkerNode::CrosslinkerNode(
    int id, std::span<const double> coords, const int owner)
    : Core::Nodes::Node(id, coords, owner), mat_(nullptr)
{
  return;
}

/*----------------------------------------------------------------------------*
 |  copy-ctor (public)                                         eichinger 10/16|
 *----------------------------------------------------------------------------*/
CrossLinking::CrosslinkerNode::CrosslinkerNode(const CrossLinking::CrosslinkerNode& old)
    : Core::Nodes::Node(old)
{
  FOUR_C_THROW(
      "Copy constructor of CrosslinkerNodeDataContainer needs to "
      "implemented first");
  return;
}

/*----------------------------------------------------------------------------*
 |  Deep copy this instance and return pointer to it (public)                 |
 |                                                             eichinger 10/16|
 *----------------------------------------------------------------------------*/
CrossLinking::CrosslinkerNode* CrossLinking::CrosslinkerNode::clone() const
{
  CrossLinking::CrosslinkerNode* newnode = new CrossLinking::CrosslinkerNode(*this);
  return newnode;
}

/*----------------------------------------------------------------------------*
 |  << operator                                                eichinger 10/16|
 *----------------------------------------------------------------------------*/
std::ostream& operator<<(std::ostream& os, const CrossLinking::CrosslinkerNode& crosslinker_node)
{
  crosslinker_node.print(os);
  return os;
}

/*----------------------------------------------------------------------------*
 |  print this CrossslinkerNode (public)                       eichinger 10/16|
 *----------------------------------------------------------------------------*/
void CrossLinking::CrosslinkerNode::print(std::ostream& os) const
{
  // Print id and coordinates
  os << "Crosslinker ";
  Core::Nodes::Node::print(os);


  return;
}

/*----------------------------------------------------------------------------*
 |  Pack data                                                        (public) |
 |                                                             eichinger 10/16|
 *----------------------------------------------------------------------------*/
void CrossLinking::CrosslinkerNode::pack(Core::Communication::PackBuffer& data) const
{
  // pack type of this instance of ParObject
  int type = unique_par_object_id();
  add_to_pack(data, type);
  // add base class Core::Nodes::Node
  Core::Nodes::Node::pack(data);

  // add material
  bool hasmat = (mat_ != nullptr);
  add_to_pack(data, hasmat);
  if (hasmat) mat_->pack(data);

  return;
}

/*----------------------------------------------------------------------------*
 |  Unpack data                                                      (public) |
 |                                                             eichinger 10/16|
 *----------------------------------------------------------------------------*/
void CrossLinking::CrosslinkerNode::unpack(Core::Communication::UnpackBuffer& buffer)
{
  Core::Communication::extract_and_assert_id(buffer, unique_par_object_id());

  // extract base class Core::Nodes::Node
  Core::Nodes::Node::unpack(buffer);

  // mat
  bool hasmat;
  extract_from_pack(buffer, hasmat);
  if (hasmat)
  {
    Core::Communication::ParObject* o = Core::Communication::factory(buffer);
    Mat::CrosslinkerMat* mat = dynamic_cast<Mat::CrosslinkerMat*>(o);
    if (mat == nullptr) FOUR_C_THROW("failed to unpack material");
    // unpack material
    mat_ = std::shared_ptr<Mat::CrosslinkerMat>(mat);
  }
  else
  {
    mat_ = nullptr;
  }


  return;
}

///*----------------------------------------------------------------------------*
// |  Initialize data container                                  eichinger 10/16|
// *----------------------------------------------------------------------------*/
// void CrossLinking::CrosslinkerNode::initialize_data_container()
//{
//  // only initialize if not yet done
//  if (cldata_ ==nullptr)
//    cldata_ = Teuchos::rcp(new CrossLinking::CrosslinkerNodeDataContainer());
//
//  return;
//}

/*----------------------------------------------------------------------------*
 |  create material class (public)                             eichinger 10/16|
 *---------------------------------------------------------------------------*/
void CrossLinking::CrosslinkerNode::set_material(int const matnum)
{
  std::shared_ptr<Mat::CrosslinkerMat> mat =
      std::dynamic_pointer_cast<Mat::CrosslinkerMat>(Mat::factory(matnum));
  if (mat == nullptr) FOUR_C_THROW("Invalid material given to crosslinker node. \n");
  mat_ = mat;
}

/*----------------------------------------------------------------------------*
 |  create material class (public)                             eichinger 10/16|
 *---------------------------------------------------------------------------*/
void CrossLinking::CrosslinkerNode::set_material(std::shared_ptr<Core::Mat::Material> material)
{
  std::shared_ptr<Mat::CrosslinkerMat> mat =
      std::dynamic_pointer_cast<Mat::CrosslinkerMat>(material);
  if (mat == nullptr) FOUR_C_THROW("Invalid material given to crosslinker node. \n");
  mat_ = mat;
}

///*----------------------------------------------------------------------------*
// |  Reset data container                                       eichinger 10/16|
// *----------------------------------------------------------------------------*/
// void CrossLinking::CrosslinkerNode::ResetDataContainer()
//{
//  // reset to nullptr
//  cldata_  = nullptr;
//
//  return;
//}

FOUR_C_NAMESPACE_CLOSE
