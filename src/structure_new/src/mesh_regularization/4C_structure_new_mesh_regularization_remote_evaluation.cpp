// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_structure_new_mesh_regularization_remote_evaluation.hpp"

#include "4C_comm_exporter.hpp"
#include "4C_comm_mpi_utils.hpp"
#include "4C_linalg_map.hpp"
#include "4C_linalg_utils_densematrix_communication.hpp"

#include <map>
#include <numeric>

FOUR_C_NAMESPACE_OPEN

std::vector<std::vector<double>> Solid::MeshRegularization::answer_on_processors(MPI_Comm comm,
    const std::vector<RemoteRequest>& requests,
    const std::function<std::vector<double>(std::span<const double> request)>& answer)
{
  const int num_processors = Core::Communication::num_mpi_ranks(comm);
  const int my_rank = Core::Communication::my_mpi_rank(comm);

  // unique ids of the requests: numbered consecutively over the processors
  const std::vector<int> counts =
      Core::Communication::all_gather(static_cast<int>(requests.size()), comm);
  const int offset = std::accumulate(counts.begin(), counts.begin() + my_rank, 0);

  std::vector<int> request_ids(requests.size());
  std::vector<std::vector<int>> ids_to_processors(num_processors);
  std::map<int, std::vector<double>> data;
  for (std::size_t i = 0; i < requests.size(); ++i)
  {
    request_ids[i] = offset + static_cast<int>(i);
    ids_to_processors[requests[i].processor].push_back(request_ids[i]);
    data[request_ids[i]] = requests[i].data;
  }

  // each processor learns the ids of the requests addressed to it
  std::vector<int> received_ids;
  Core::LinAlg::all_to_all_communication(comm, ids_to_processors, received_ids);

  const Core::LinAlg::Map request_map(
      -1, static_cast<int>(request_ids.size()), request_ids.data(), 0, comm);
  const Core::LinAlg::Map received_map(
      -1, static_cast<int>(received_ids.size()), received_ids.data(), 0, comm);

  Core::Communication::Exporter(request_map, received_map, comm).do_export(data);

  std::map<int, std::vector<double>> answers;
  for (const auto& [id, request] : data) answers[id] = answer(request);

  Core::Communication::Exporter(received_map, request_map, comm).do_export(answers);

  std::vector<std::vector<double>> result(requests.size());
  for (std::size_t i = 0; i < requests.size(); ++i) result[i] = std::move(answers.at(request_ids[i]));
  return result;
}

FOUR_C_NAMESPACE_CLOSE
