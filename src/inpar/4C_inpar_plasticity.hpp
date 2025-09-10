// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_INPAR_PLASTICITY_HPP
#define FOUR_C_INPAR_PLASTICITY_HPP

#include "4C_config.hpp"

#include "4C_io_input_spec.hpp"


FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
namespace Inpar
{
  namespace Plasticity
  {
    /// plasticity parameters
    Core::IO::InputSpec valid_parameters();

  }  // namespace Plasticity

}  // namespace Inpar

/*----------------------------------------------------------------------*/
FOUR_C_NAMESPACE_CLOSE

#endif
