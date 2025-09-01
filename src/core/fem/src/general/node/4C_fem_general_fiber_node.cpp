// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_fem_general_fiber_node.hpp"

#include "4C_comm_pack_helpers.hpp"

#include <utility>

FOUR_C_NAMESPACE_OPEN

Core::Nodes::FiberNodeType Core::Nodes::FiberNodeType::instance_;


Core::Communication::ParObject* Core::Nodes::FiberNodeType::create(
    Core::Communication::UnpackBuffer& buffer)
{
  std::vector<double> dummy_coords(3, 999.0);
  std::map<CoordinateSystemDirection, std::array<double, 3>> coordinateSystemDirections;
  std::vector<std::array<double, 3>> fibers;
  std::map<AngleType, double> angles;
  auto* object = new FiberNode(-1, dummy_coords, coordinateSystemDirections, fibers, angles, -1);
  object->unpack(buffer);
  return object;
}

Core::Nodes::FiberNode::FiberNode(int id, std::span<const double> coords,
    std::map<CoordinateSystemDirection, std::array<double, 3>> coordinateSystemDirections,
    std::vector<std::array<double, 3>> fibers, std::map<AngleType, double> angles, const int owner)
    : Core::Nodes::Node(id, coords, owner),
      coordinateSystemDirections_(std::move(coordinateSystemDirections)),
      fibers_(std::move(fibers)),
      angles_(std::move(angles))
{
}

/*
  Deep copy the derived class and return pointer to it
*/
Core::Nodes::FiberNode* Core::Nodes::FiberNode::clone() const
{
  auto* newfn = new Core::Nodes::FiberNode(*this);

  return newfn;
}

/*
  Pack this class so it can be communicated

  Pack and Unpack are used to communicate this fiber node

*/
void Core::Nodes::FiberNode::pack(Core::Communication::PackBuffer& data) const
{
  // pack type of this instance of ParObject
  int type = unique_par_object_id();
  add_to_pack(data, type);
  // add base class of fiber node
  Core::Nodes::Node::pack(data);

  // Add fiber data
  add_to_pack(data, fibers_);
  add_to_pack(data, coordinateSystemDirections_);
  add_to_pack(data, angles_);
}

/*
  Unpack data from a char vector into this class

  Pack and Unpack are used to communicate this fiber node
*/
void Core::Nodes::FiberNode::unpack(Core::Communication::UnpackBuffer& buffer)
{
  Core::Communication::extract_and_assert_id(buffer, unique_par_object_id());

  Node::unpack(buffer);

  // extract fiber data
  extract_from_pack(buffer, fibers_);
  extract_from_pack(buffer, coordinateSystemDirections_);
  extract_from_pack(buffer, angles_);
}

/*
  Print this fiber node
*/
void Core::Nodes::FiberNode::print(std::ostream& os) const
{
  os << "Fiber Node :";
  Core::Nodes::Node::print(os);
  os << "(" << fibers_.size() << " fibers, " << angles_.size() << " angles)";
}

FOUR_C_NAMESPACE_CLOSE
