// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_fem_general_element_definition.hpp"

#include "4C_comm_parobjectfactory.hpp"
#include "4C_io_input_parameter_container.hpp"
#include "4C_io_input_spec_builders.hpp"

FOUR_C_NAMESPACE_OPEN


Core::Elements::ElementDefinition::ElementDefinition()
{
  Core::Communication::ParObjectFactory::instance().setup_element_definition(definitions);
}

const Core::IO::InputSpec& Core::Elements::ElementDefinition::get(
    const std::string& element_name, Core::FE::CellType cell_type) const
{
  auto it = definitions.find(element_name);
  if (it == definitions.end()) FOUR_C_THROW("No element '{}' found.", element_name);
  auto it2 = it->second.find(cell_type);
  if (it2 == it->second.end())
    FOUR_C_THROW("Element '{}' does not seem to know cell type '{}'.", element_name, cell_type);
  return it2->second;
}

std::tuple<std::string, Core::FE::CellType, Core::IO::InputParameterContainer>
Core::Elements::ElementDefinition::unpack_element_data(
    const Core::IO::InputParameterContainer& data) const
{
  const auto& [element_name, element_group] =
      data.exactly_one_group(definitions | std::views::keys);

  const auto& cell_specs = definitions.at(element_name);
  const auto& [cell_type_name, specific_data] = element_group.exactly_one_group(
      cell_specs | std::views::keys | std::views::transform(Core::FE::cell_type_to_string));

  return {element_name, Core::FE::string_to_cell_type(std::string(cell_type_name)), specific_data};
}


Core::IO::InputSpec Core::Elements::ElementDefinition::element_data_spec() const
{
  using namespace Core::IO::InputSpecBuilders;
  std::vector<Core::IO::InputSpec> element_choices;
  for (const auto& [element, cell_specs] : definitions)
  {
    std::vector<Core::IO::InputSpec> cell_specs_choices;
    for (const auto& [cell_type, spec] : cell_specs)
    {
      cell_specs_choices.emplace_back(group(Core::FE::cell_type_to_string(cell_type), {spec}));
    }
    element_choices.emplace_back(group(element, {one_of(std::move(cell_specs_choices))}));
  }
  return one_of(std::move(element_choices));
}

FOUR_C_NAMESPACE_CLOSE
