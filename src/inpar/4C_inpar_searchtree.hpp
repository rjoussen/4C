// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_INPAR_SEARCHTREE_HPP
#define FOUR_C_INPAR_SEARCHTREE_HPP

#include "4C_config.hpp"

#include "4C_io_input_spec.hpp"


FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
namespace Inpar
{
  namespace Geo
  {
    /// specify tree type
    enum TreeType
    {
      Notree,
      Octree3D,
      Quadtree3D,
      Quadtree2D
    };

    /// searchtree parameters
    Core::IO::InputSpec valid_parameters();

  }  // namespace Geo

}  // namespace Inpar

/*----------------------------------------------------------------------*/
FOUR_C_NAMESPACE_CLOSE

#endif
