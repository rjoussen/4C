// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_FLUID_ELE_NULLSPACE_HPP
#define FOUR_C_FLUID_ELE_NULLSPACE_HPP

#include "4C_config.hpp"

#include "4C_linalg_serialdensematrix.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Core::Nodes
{
  class Node;
}

namespace FLD
{
  template <unsigned dim>
  Core::LinAlg::SerialDenseMatrix compute_fluid_null_space()
  {
    /* !\brief Helper function for the nodal nullspace of fluid elements
      The rigid body modes for fluids are:

                xtrans   ytrans  ztrans   pressure
                mode[0]  mode[1] mode[2]  mode[3]
          ----------------------------------------
          x   |    1       0       0       0
          y   |    0       1       0       0
          z   |    0       0       1       0
          p   |    0       0       0       1

          valid element types: fluid3, xfluid3
      */
    Core::LinAlg::SerialDenseMatrix nullspace(dim, dim, true);

    for (unsigned int i = 0; i < dim; i++)
    {
      for (unsigned int j = 0; j < dim; j++)
      {
        if (i == j) nullspace(i, j) = 1.0;
      }
    }

    return nullspace;
  }
}  // namespace FLD

FOUR_C_NAMESPACE_CLOSE

#endif
