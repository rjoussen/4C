// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_io_input_spec_validators.hpp"

#include <regex>

FOUR_C_NAMESPACE_OPEN

Core::IO::InputSpecBuilders::Validators::Validator<std::string>
Core::IO::InputSpecBuilders::Validators::pattern(std::string pattern)
{
  return Validator<std::string>([re = std::regex(pattern, std::regex::ECMAScript)](
                                    const std::string& v) { return std::regex_search(v, re); },
      [pattern](std::ostream& os) { os << "pattern{" << pattern << "}"; },
      [pattern](YamlNodeRef yaml)
      {
        auto& node = yaml.node;
        node |= ryml::MAP;
        node["pattern"] |= ryml::MAP;
        node["pattern"]["pattern"] << pattern;
      });
}

FOUR_C_NAMESPACE_CLOSE
