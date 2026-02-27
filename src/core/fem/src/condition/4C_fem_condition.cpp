// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_fem_condition.hpp"

#include "4C_fem_general_element.hpp"

#include <string>
#include <utility>
#include <variant>

FOUR_C_NAMESPACE_OPEN


Core::Conditions::Condition::Condition(const std::variant<int, std::string>& id,
    const Core::Conditions::ConditionType type, const bool buildgeometry,
    const Core::Conditions::GeometryType gtype, EntityType entity_type)
    : buildgeometry_(buildgeometry), type_(type), gtype_(gtype), entity_type_(entity_type)
{
  if (entity_type_ == EntityType::node_set_name)
  {
    FOUR_C_ASSERT(std::holds_alternative<std::string>(id),
        "Condition of ENTITY_TYPE 'node_set_name' requires a string id referring to the node set "
        "name.");
    node_set_name_ = std::get<std::string>(id);
  }
  else
  {
    FOUR_C_ASSERT(std::holds_alternative<int>(id),
        "Condition with ID-based identification requires an integer id referring to the respective "
        "entity.");
    id_ = std::get<int>(id);
  }
}

std::ostream& operator<<(std::ostream& os, const Core::Conditions::Condition& cond)
{
  cond.print(os);
  return os;
}


void Core::Conditions::Condition::print(std::ostream& os) const
{
  os << "Condition "
     << (node_set_name_.has_value() ? node_set_name_.value() : std::to_string(id_.value())) << " "
     << to_string(type_) << ": ";
  container_.print(os);
  os << std::endl;
  if (nodes_.size() != 0)
  {
    os << "Nodes of this condition:";
    for (const auto& node_gid : nodes_) os << " " << node_gid;
    os << std::endl;
  }
  if (!geometry_.empty())
  {
    os << "Elements of this condition:";
    for (const auto& [ele_id, ele] : geometry_) os << " " << ele_id;
    os << std::endl;
  }
}

std::shared_ptr<Core::Conditions::Condition> Core::Conditions::Condition::copy_without_geometry()
    const
{
  std::shared_ptr<Core::Conditions::Condition> copy(new Condition(*this));
  copy->clear_geometry();
  return copy;
}

void Core::Conditions::Condition::resolve_geometry(const GeometryContext& geometry_context)
{
  switch (entity_type_)
  {
    case Core::Conditions::EntityType::legacy_id:
    {
      if (geometry_context.dnode_fenode.size() == 0 && geometry_context.dline_fenode.size() == 0 &&
          geometry_context.dsurf_fenode.size() == 0 && geometry_context.dvol_fenode.size() == 0)
      {
        FOUR_C_THROW(
            "{} condition {} uses legacy_id entity type but no legacy entities were defined in "
            "the input file.\n"
            "This is probably because the geometry is handled in an external file.\n"
            "If this is the case, you must specify a specific entity type (node_set_id or "
            "element_block_id) or identify the node set via its name.\n",
            type_, id_.value());
      }
      switch (gtype_)
      {
        case Core::Conditions::geometry_type_point:
          if (id_ < 0 or static_cast<unsigned>(id_.value()) >= geometry_context.dnode_fenode.size())
          {
            FOUR_C_THROW(
                "DPoint {} not in range [0:{}[\n"
                "DPoint condition on non existent DPoint?"
                "Could not read set from entity type.",
                id_.value(), geometry_context.dnode_fenode.size());
          }
          set_nodes(geometry_context.dnode_fenode[id_.value()]);
          break;
        case Core::Conditions::geometry_type_line:
          if (id_ < 0 or static_cast<unsigned>(id_.value()) >= geometry_context.dline_fenode.size())
          {
            FOUR_C_THROW(
                "DLine {} not in range [0:{}[\n"
                "DLine condition on non existent DLine?"
                "Could not read set from entity type.",
                id_.value(), geometry_context.dline_fenode.size());
          }
          set_nodes(geometry_context.dline_fenode[id_.value()]);
          break;
        case Core::Conditions::geometry_type_surface:
          if (id_ < 0 or static_cast<unsigned>(id_.value()) >= geometry_context.dsurf_fenode.size())
          {
            FOUR_C_THROW(
                "DSurface {} not in range [0:{}[\n"
                "DSurface condition on non existent DSurface?"
                "Could not read set from entity type.",
                id_.value(), geometry_context.dsurf_fenode.size());
          }
          set_nodes(geometry_context.dsurf_fenode[id_.value()]);
          break;
        case Core::Conditions::geometry_type_volume:
          if (id_ < 0 or static_cast<unsigned>(id_.value()) >= geometry_context.dvol_fenode.size())
          {
            FOUR_C_THROW(
                "DVolume {} not in range [0:{}[\n"
                "DVolume condition on non existent DVolume?",
                id_.value(), geometry_context.dvol_fenode.size());
          }
          set_nodes(geometry_context.dvol_fenode[id_.value()]);
          break;
        default:
          FOUR_C_THROW("geometry type unspecified");
          break;
      }

      break;
    }
    case Core::Conditions::EntityType::node_set_id:
    {
      FOUR_C_ASSERT_ALWAYS(geometry_context.node_sets.contains(id_.value()),
          "Cannot apply condition '{}' to node set {} which is not specified in the mesh file.",
          type_, id_.value());
      set_nodes(geometry_context.node_sets.at(id_.value()));
      break;
    }
    case Core::Conditions::EntityType::node_set_name:
    {
      FOUR_C_ASSERT_ALWAYS(geometry_context.node_sets_names.contains(node_set_name_.value()),
          "NODE_SET_NAME '{}' could not be found in the meshfile.", node_set_name_.value());

      const auto& ids = geometry_context.node_sets_names.at(node_set_name_.value());
      FOUR_C_ASSERT_ALWAYS(ids.size() == 1,
          "NODE_SET_NAME '{}' is not unique in the meshfile ({} occurrences).",
          node_set_name_.value(), ids.size());

      set_nodes(geometry_context.node_sets.at(ids[0]));
      break;
    }
    case Core::Conditions::EntityType::element_block_id:
    {
      FOUR_C_ASSERT_ALWAYS(geometry_context.element_block_nodes.contains(id_.value()),
          "Cannot apply condition '{}' to element block {} which is not specified in the mesh "
          "file.",
          type_, id_.value());
      set_nodes(geometry_context.element_block_nodes.at(id_.value()));
      break;
    }
  }
}

int Core::Conditions::Condition::id() const
{
  FOUR_C_ASSERT(entity_type_ != EntityType::node_set_name,
      "Condition {} of ENTITY_TYPE 'node_set_name' has no id assigned!", node_set_name_.value());
  return id_.value();
}

std::string Core::Conditions::Condition::node_set_name() const
{
  FOUR_C_ASSERT(entity_type_ == EntityType::node_set_name,
      "Condition of ENTITY_TYPE 'node_set_name' has no node_set_name assigned!");
  return node_set_name_.value();
}


FOUR_C_NAMESPACE_CLOSE
