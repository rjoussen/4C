// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_FEM_GEOMETRY_UPDATE_REFERENCE_CONFIG_HPP
#define FOUR_C_FEM_GEOMETRY_UPDATE_REFERENCE_CONFIG_HPP


#include "4C_config.hpp"

#include "4C_fem_general_element.hpp"
#include "4C_fem_general_node.hpp"

#include <memory>

FOUR_C_NAMESPACE_OPEN

namespace Core::Geo
{
  //! Update material configuration of @p dis with @p disp
  void update_reference_config_with_disp(
      Core::FE::Discretization& dis, const Core::LinAlg::Vector<double>& disp);
}  // namespace Core::Geo

FOUR_C_NAMESPACE_CLOSE

#endif
