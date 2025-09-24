// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_FEM_DISCRETIZATION_NULLSPACE_HPP
#define FOUR_C_FEM_DISCRETIZATION_NULLSPACE_HPP

#include "4C_config.hpp"

#include "4C_linalg_multi_vector.hpp"

#include <memory>

FOUR_C_NAMESPACE_OPEN

namespace Core::FE
{
  class Discretization;

  /*!
   \brief Calculate the nullspace based on a given discretization

  The nullspace is build by looping over all nodes of a discretization and stored
          in the respective variable.

     \param dis (in): discretization
     \param dimns (in): nullspace dimension
     \param map (in): nullspace map
      */
  std::shared_ptr<Core::LinAlg::MultiVector<double>> compute_null_space(
      const Core::FE::Discretization& dis, const int dimns, const Core::LinAlg::Map& dofmap);
}  // namespace Core::FE

FOUR_C_NAMESPACE_CLOSE

#endif
