// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_utils_flat_vector_vector.hpp"

#include <benchmark/benchmark.h>

namespace
{

  void fill_data(std::vector<std::vector<int>>& data, size_t n_outer)
  {
    for (size_t i = 0; i < n_outer; ++i)
    {
      const auto n_inner = n_outer % 10 + 1;
      data[i].resize(n_inner);
      for (size_t j = 0; j < n_inner; ++j)
      {
        data[i][j] = (i + j) % 4;
      }
    }
  }
  void flat_vector_vector_naive(benchmark::State& state)
  {
    size_t n_outer = 100;
    std::vector<std::vector<int>> data(n_outer);
    fill_data(data, n_outer);

    int result = 0;
    for (auto _ : state)
    {
      for (size_t i = 0; i < n_outer; ++i)
        for (size_t j = 0; j < data[i].size(); ++j) result += data[i][j];
      benchmark::DoNotOptimize(result);
    }
  }
  BENCHMARK(flat_vector_vector_naive);

  void flat_vector_vector_opt(benchmark::State& state)
  {
    size_t n_outer = 100;
    std::vector<std::vector<int>> data(n_outer);
    fill_data(data, n_outer);

    FourC::Core::Utils::FlatVectorVector<int> flat_data(data);
    int result = 0;
    for (auto _ : state)
    {
      for (size_t i = 0; i < n_outer; ++i)
        for (size_t j = 0; j < flat_data[i].size(); ++j) result += flat_data[i][j];
      benchmark::DoNotOptimize(result);
    }
  }
  BENCHMARK(flat_vector_vector_opt);
}  // namespace
