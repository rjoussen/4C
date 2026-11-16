// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_STRUCTURE_NEW_MESH_REGULARIZATION_REMOTE_EVALUATION_HPP
#define FOUR_C_STRUCTURE_NEW_MESH_REGULARIZATION_REMOTE_EVALUATION_HPP

#include "4C_config.hpp"

#include <mpi.h>

#include <functional>
#include <span>
#include <utility>
#include <vector>

FOUR_C_NAMESPACE_OPEN

namespace Solid::MeshRegularization
{
  /*!
   * @brief A request to another processor: the processor and the data of the request
   */
  struct RemoteRequest
  {
    int processor;
    std::vector<double> data;
  };

  /*!
   * @brief Answer requests on the processors they are addressed to
   *
   * Each request is sent to its processor, which answers it with @p answer, and the answer is sent
   * back. This lets a processor evaluate data owned by other processors (e.g., the state of their
   * elements at a point) without changing the distribution of the discretization. Collective.
   *
   * @param requests The requests of this processor
   * @param answer The answer to a request on the processor it is addressed to
   * @return The answers to the requests of this processor, in their order
   */
  std::vector<std::vector<double>> answer_on_processors(MPI_Comm comm,
      const std::vector<RemoteRequest>& requests,
      const std::function<std::vector<double>(std::span<const double> request)>& answer);
}  // namespace Solid::MeshRegularization

FOUR_C_NAMESPACE_CLOSE

#endif
