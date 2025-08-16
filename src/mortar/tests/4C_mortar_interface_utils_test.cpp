// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <gtest/gtest.h>

#include "4C_mortar_interface_utils.hpp"

namespace
{
  using namespace FourC;

  TEST(InterfaceUtilsTestSuite, TestComputeParallelDistributionStatistics)
  {
    std::vector<int> quantityAcrossAllRanks = {2, 1, 4, 3};

    int minElement = 0;
    int maxElement = 0;
    double meanOfAllElements = 0;

    Mortar::InterfaceUtils::compute_parallel_distribution_statistics(
        quantityAcrossAllRanks, minElement, maxElement, meanOfAllElements);

    EXPECT_EQ(minElement, 1);
    EXPECT_EQ(maxElement, 4);
    EXPECT_EQ(meanOfAllElements, 2.5);
  }
}  // namespace
