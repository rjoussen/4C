// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_fem_general_node.hpp"

#include "4C_comm_pack_helpers.hpp"
#include "4C_utils_exceptions.hpp"

FOUR_C_NAMESPACE_OPEN


Core::Nodes::NodeType Core::Nodes::NodeType::instance_;


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Core::Communication::ParObject* Core::Nodes::NodeType::create(
    Core::Communication::UnpackBuffer& buffer)
{
  std::vector<double> dummycoord(3, 999.0);
  auto* object = new Core::Nodes::Node(-1, dummycoord, -1);
  object->unpack(buffer);
  return object;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Core::Nodes::Node::Node(const int id, std::span<const double> coords, const int owner)
    : ParObject(), id_(id), lid_(-1), owner_(owner), x_(coords.begin(), coords.end())
{
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Core::Nodes::Node* Core::Nodes::Node::clone() const
{
  auto* newnode = new Core::Nodes::Node(*this);
  return newnode;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
std::ostream& operator<<(std::ostream& os, const Core::Nodes::Node& node)
{
  node.print(os);
  return os;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Core::Nodes::Node::print(std::ostream& os) const
{
  // Print id and coordinates
  os << "Node " << std::setw(12) << id() << " Owner " << std::setw(4) << owner() << " Coords "
     << std::setw(12) << x()[0] << " " << std::setw(12) << x()[1] << " " << std::setw(12) << x()[2]
     << " ";
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Core::Nodes::Node::pack(Core::Communication::PackBuffer& data) const
{
  // pack type of this instance of ParObject
  int type = unique_par_object_id();
  add_to_pack(data, type);
  // add id
  add_to_pack(data, id());
  // add owner
  add_to_pack(data, owner());
  // x_
  add_to_pack(data, x_);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Core::Nodes::Node::unpack(Core::Communication::UnpackBuffer& buffer)
{
  Core::Communication::extract_and_assert_id(buffer, unique_par_object_id());

  // id_
  extract_from_pack(buffer, id_);
  // owner_
  extract_from_pack(buffer, owner_);
  // x_
  extract_from_pack(buffer, x_);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Core::Nodes::Node::get_condition(
    const std::string& name, std::vector<Core::Conditions::Condition*>& out) const
{
  const int num = condition_.count(name);
  out.resize(num);
  auto startit = condition_.lower_bound(name);
  auto endit = condition_.upper_bound(name);
  int count = 0;
  std::multimap<std::string, std::shared_ptr<Core::Conditions::Condition>>::const_iterator curr;
  for (curr = startit; curr != endit; ++curr) out[count++] = curr->second.get();
  if (count != num) FOUR_C_THROW("Mismatch in number of conditions found");
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Core::Nodes::Node::change_pos(std::vector<double> nvector)
{
  FOUR_C_ASSERT(x_.size() == nvector.size(),
      "Mismatch in size of the nodal coordinates vector and the vector to change the nodal "
      "position");
  for (std::size_t i = 0; i < x_.size(); ++i) x_[i] = x_[i] + nvector[i];
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Core::Nodes::Node::set_pos(std::vector<double> nvector)
{
  FOUR_C_ASSERT(x_.size() == nvector.size(),
      "Mismatch in size of the nodal coordinates vector and the vector to set the new nodal "
      "position");
  for (std::size_t i = 0; i < x_.size(); ++i) x_[i] = nvector[i];
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool Core::Nodes::Node::vis_data(const std::string& name, std::vector<double>& data)
{
  if (name == "Nodeowner")
  {
    if (static_cast<int>(data.size()) < 1) FOUR_C_THROW("Size mismatch");
    data[0] = owner();
    return true;
  }
  return false;
}

FOUR_C_NAMESPACE_CLOSE
