// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_MORTAR_INTERFACE_UTILS_HPP
#define FOUR_C_MORTAR_INTERFACE_UTILS_HPP

#include "4C_config.hpp"

#include <string>
#include <vector>

FOUR_C_NAMESPACE_OPEN

namespace Mortar
{
  namespace InterfaceUtils
  {
    /*!
    \brief Compute statistics on parallel distribution

    \param[in] quantityAcrossAllRanks Vector with number of nodes/elements per rank
    \param[out] minOverAllRanks Smallest element in input vector
    \param[out] maxOverAllRanks Largest element in input vector
    \param[out] meanOverAllRanks Mean of all elements in input vector
    */
    void compute_parallel_distribution_statistics(const std::vector<int>& quantityAcrossAllRanks,
        int& minOverAllRanks, int& maxOverAllRanks, double& meanOverAllRanks);

    /*!
    \brief Compute and print one row of statistics on parallel distribution table

    \param[in] nameOfQuantity Name describing the quantity to be printed
    \param[in] quantityAcrossAllRanks Vector with number of nodes/elements per rank
    \param[in] printOnThisRank Flag to en-/disable screen output on this rank
    */
    void compute_and_print_row_of_parallel_distribution_statisctics(
        const std::string& nameOfQuantity, const std::vector<int>& quantityAcrossAllRanks,
        const bool printOnThisRank);

  }  // namespace InterfaceUtils
}  // namespace Mortar

FOUR_C_NAMESPACE_CLOSE

#endif
