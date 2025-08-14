// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_io_input_parameter_container.templates.hpp"

#include <Teuchos_ParameterList.hpp>

FOUR_C_NAMESPACE_OPEN

void Core::IO::InputParameterContainer::print(std::ostream& os) const
{
  for (const auto& [key, entry] : entries_)
  {
    os << key << " : ";
    get_type_actions().at(entry.data.type()).print(os, entry.data);
  }

  for (const auto& [key, group] : groups_)
  {
    os << key << " : ";
    group.print(os);
  }
}


Core::IO::InputParameterContainer& Core::IO::InputParameterContainer::group(const std::string& name)
{
  FOUR_C_ASSERT_ALWAYS(!name.empty(), "Group name must not be empty.");
  return groups_[name];
}


const Core::IO::InputParameterContainer& Core::IO::InputParameterContainer::group(
    const std::string& name) const
{
  FOUR_C_ASSERT_ALWAYS(groups_.contains(name), "Group '{}' not found in container.", name);
  return groups_.at(name);
}

void Core::IO::InputParameterContainer::add_list(const std::string& name, List&& list)
{
  FOUR_C_ASSERT_ALWAYS(!entries_.contains(name), "Key '{}' already exists in container.", name);

  entries_[name] = {std::move(list)};

  auto& type_actions = Core::IO::InputParameterContainer::get_type_actions();
  if (!type_actions.contains(typeid(List)))
  {
    type_actions[typeid(List)] = {
        .print =
            [](std::ostream& os, const std::any& data)
        {
          const auto& list = std::any_cast<List>(data);
          for (const auto& container : list)
          {
            container.print(os);
          }
        },
        .write_to_pl =
            [](Teuchos::ParameterList& pl, const std::string& key, const std::any& data)
        {
          const auto& list = std::any_cast<List>(data);
          std::vector<Teuchos::ParameterList> sublists;
          for (const auto& container : list)
          {
            container.to_teuchos_parameter_list(sublists.emplace_back());
          }
          pl.set(key, sublists);
        },
    };
  }
}

auto Core::IO::InputParameterContainer::get_list(const std::string& name) const -> const List&
{
  return get<List>(name);
}

bool Core::IO::InputParameterContainer::has_group(const std::string& name) const
{
  return groups_.contains(name);
}

void Core::IO::InputParameterContainer::merge(const Core::IO::InputParameterContainer& other)
{
  const auto combine_maps = [](auto& map1, const auto& map2)
  {
    for (const auto& [key, value] : map2)
    {
      if (map1.contains(key))
      {
        FOUR_C_THROW("Duplicate key '{}' encountered while merging two containers.", key);
      }
      map1[key] = value;
    }
  };

  combine_maps(entries_, other.entries_);
  combine_maps(groups_, other.groups_);
}


void Core::IO::InputParameterContainer::clear()
{
  entries_.clear();
  groups_.clear();
}

void Core::IO::InputParameterContainer::to_teuchos_parameter_list(
    Teuchos::ParameterList& list) const
{
  for (const auto& [key, entry] : entries_)
  {
    get_type_actions().at(entry.data.type()).write_to_pl(list, key, entry.data);
  }

  for (const auto& [key, group] : groups_)
  {
    group.to_teuchos_parameter_list(list.sublist(key));
  }
}

std::map<std::type_index, Core::IO::InputParameterContainer::TypeActions>&
Core::IO::InputParameterContainer::get_type_actions()
{
  static std::map<std::type_index, TypeActions> type_actions;
  return type_actions;
}

// Explicit instantiation of the template functions.

template void Core::IO::InputParameterContainer::add(const std::string& name, const int& data);
template void Core::IO::InputParameterContainer::add(const std::string& name, const double& data);
template void Core::IO::InputParameterContainer::add(
    const std::string& name, const std::string& data);
template void Core::IO::InputParameterContainer::add(const std::string& name, const bool& data);

template void Core::IO::InputParameterContainer::add(
    const std::string& name, const std::vector<int>& data);
template void Core::IO::InputParameterContainer::add(
    const std::string& name, const std::vector<double>& data);
template void Core::IO::InputParameterContainer::add(
    const std::string& name, const std::vector<std::string>& data);
template void Core::IO::InputParameterContainer::add(
    const std::string& name, const std::vector<bool>& data);

FOUR_C_NAMESPACE_CLOSE
