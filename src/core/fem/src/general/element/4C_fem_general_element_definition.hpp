// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_FEM_GENERAL_ELEMENT_DEFINITION_HPP
#define FOUR_C_FEM_GENERAL_ELEMENT_DEFINITION_HPP

#include "4C_config.hpp"

#include "4C_fem_general_cell_type.hpp"
#include "4C_io_input_spec.hpp"

#include <map>
#include <string>

FOUR_C_NAMESPACE_OPEN


namespace Core::Elements
{
  /**
   * Collection of valid element input.
   *
   * This glass gathers all valid element definitions from the global state. This is possible
   * since the ElementType subclass register themselves in the global Factory.
   */
  struct ElementDefinition
  {
    //! Gather all valid element definitions from global state
    ElementDefinition();

    /**
     * Convenience access with check for existence of element definition.
     */
    const Core::IO::InputSpec& get(
        const std::string& element_name, Core::FE::CellType cell_type) const;

    /**
     * Get an InputSpec that describes all valid element definitions.
     */
    [[nodiscard]] Core::IO::InputSpec element_data_spec() const;

    /**
     * Given a @p data container that matches the element_data_spec(), unpack the information
     * into a tuple of (element_name, cell_type, specific_data), where specific_data is a
     * container that matches the spec from get(element_name, cell_type).
     */
    [[nodiscard]] std::tuple<std::string, Core::FE::CellType, Core::IO::InputParameterContainer>
    unpack_element_data(const Core::IO::InputParameterContainer& data) const;

    //! Map from physics to cell type to InputSpec.
    std::map<std::string, std::map<Core::FE::CellType, Core::IO::InputSpec>> definitions;
  };

}  // namespace Core::Elements


FOUR_C_NAMESPACE_CLOSE

#endif
