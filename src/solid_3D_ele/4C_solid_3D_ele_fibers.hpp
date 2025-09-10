// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_SOLID_3D_ELE_FIBERS_HPP
#define FOUR_C_SOLID_3D_ELE_FIBERS_HPP

#include "4C_config.hpp"

#include "4C_io_input_parameter_container.hpp"
#include "4C_linalg_tensor.hpp"

#include <optional>



FOUR_C_NAMESPACE_OPEN


namespace Discret::Elements
{

  /*!
   * @brief A container holding fiber directions within an element
   */
  struct Fibers
  {
    /// List of fibers defined on the element level
    std::vector<Core::LinAlg::Tensor<double, 3>> element_fibers;
  };

  /*!
   * @brief A container holding a local coordinate system within an element
   */
  struct CoordinateSystem
  {
    /// Local coordinate system defined on the element level
    std::array<Core::LinAlg::Tensor<double, 3>, 3> element_system;
  };

  /*!
   * @brief Read fiber directions from the input definition of the element (either from the legacy
   * FIBER1, ... keywords or from input fields FIBERS)
   */
  Fibers read_fibers(const Core::IO::InputParameterContainer& input_data);

  /*!
   * @brief Read local coordinate system from the input definition of the element (either from
   * RAD,AXI,CIR or from input field COORDINATE_SYSTEM), if present.
   */
  std::optional<CoordinateSystem> read_coordinate_system(
      const Core::IO::InputParameterContainer& input_data);

}  // namespace Discret::Elements

FOUR_C_NAMESPACE_CLOSE

#endif
