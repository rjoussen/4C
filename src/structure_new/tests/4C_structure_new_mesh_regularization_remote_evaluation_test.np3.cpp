// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <gtest/gtest.h>

#include "4C_structure_new_mesh_regularization_remote_evaluation.hpp"

#include "4C_comm_mpi_utils.hpp"

#include <vector>

namespace
{
  using namespace FourC;

  TEST(AnswerOnProcessorsTest, AnswersEachRequestOnItsProcessor)
  {
    const MPI_Comm comm = MPI_COMM_WORLD;
    const int num_processors = Core::Communication::num_mpi_ranks(comm);
    const int my_rank = Core::Communication::my_mpi_rank(comm);

    // rank r sends r + 1 requests, the i-th to rank (r + i) mod n, with data (r, i); rank 0 sends
    // one to itself
    std::vector<Solid::MeshRegularization::RemoteRequest> requests;
    for (int i = 0; i <= my_rank; ++i)
    {
      requests.push_back({(my_rank + i) % num_processors,
          {static_cast<double>(my_rank), static_cast<double>(i)}});
    }

    // the answer names the answering processor and doubles the data
    const auto answers = Solid::MeshRegularization::answer_on_processors(comm, requests,
        [&](std::span<const double> request)
        {
          return std::vector<double>{
              static_cast<double>(my_rank), 2.0 * request[0], 2.0 * request[1]};
        });

    ASSERT_EQ(answers.size(), requests.size());
    for (int i = 0; i <= my_rank; ++i)
    {
      EXPECT_EQ(answers[i], (std::vector<double>{static_cast<double>((my_rank + i) % num_processors),
                                2.0 * my_rank, 2.0 * i}));
    }
  }
}  // namespace
