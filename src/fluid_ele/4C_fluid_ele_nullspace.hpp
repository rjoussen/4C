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
  /*!
  \brief Helper function for the nodal nullspace of fluid elements

    \param numdof (in):  number of degrees of freedom
    \param dimnsp (in):  dimension of the nullspace
                         */
  Core::LinAlg::SerialDenseMatrix compute_fluid_null_space(const int numdof, const int dimnsp);
}  // namespace FLD

FOUR_C_NAMESPACE_CLOSE

#endif
