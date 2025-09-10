// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_config.hpp"

#include "4C_linalg_tensor.hpp"
#include "4C_solid_3D_ele_fibers.hpp"
#include "4C_utils_exceptions.hpp"

#include <algorithm>

FOUR_C_NAMESPACE_OPEN


Discret::Elements::Fibers Discret::Elements::read_fibers(
    const Core::IO::InputParameterContainer& input_data)
{
  Fibers fibers{};

  // for now only support the old style fiber input via FIBER1, FIBER2, ... keywords
  for (int i = 1;; ++i)
  {
    const std::string fiber_name = "FIBER" + std::to_string(i);

    const auto* fiber_ptr = input_data.get_if<std::optional<std::vector<double>>>(fiber_name);
    if (!fiber_ptr || !fiber_ptr->has_value())
    {
      break;
    }
    Core::LinAlg::Tensor<double, 3> tensor{};
    std::ranges::copy(fiber_ptr->value(), tensor.data());
    // Normalize the input vectors
    tensor /= Core::LinAlg::norm2(tensor);

    fibers.element_fibers.emplace_back(tensor);
  }

  return fibers;
}

std::optional<Discret::Elements::CoordinateSystem> Discret::Elements::read_coordinate_system(
    const Core::IO::InputParameterContainer& input_data)
{
  // for now only support the old style fiber input via RAD, AXI, CIR keywords
  const std::array rad_axi_cir = {input_data.get_if<std::optional<std::vector<double>>>("RAD"),
      input_data.get_if<std::optional<std::vector<double>>>("AXI"),
      input_data.get_if<std::optional<std::vector<double>>>("CIR")};

  if (std::ranges::none_of(rad_axi_cir, [](const auto* p) { return p && p->has_value(); }))
  {
    // no local coordinate system defined
    return std::nullopt;
  }

  FOUR_C_ASSERT_ALWAYS(
      std::ranges::all_of(rad_axi_cir, [](const auto* p) { return p && p->has_value(); }),
      "If you specify a coordinate system, you need to define all of RAD, AXI and CIR!");

  CoordinateSystem coord_sys{};

  for (std::size_t i = 0; i < 3; ++i)
  {
    std::ranges::copy(rad_axi_cir[i]->value(), coord_sys.element_system[i].data());

    // Normalize the input vectors
    coord_sys.element_system[i] /= Core::LinAlg::norm2(coord_sys.element_system[i]);
  }

  return coord_sys;
}

FOUR_C_NAMESPACE_CLOSE