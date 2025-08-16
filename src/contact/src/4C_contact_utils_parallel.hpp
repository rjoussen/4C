// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_CONTACT_UTILS_PARALLEL_HPP
#define FOUR_C_CONTACT_UTILS_PARALLEL_HPP

#include "4C_config.hpp"

#include "4C_utils_parameter_list.fwd.hpp"

FOUR_C_NAMESPACE_OPEN

namespace CONTACT
{
  namespace Utils
  {
    /*!
    \brief Decide whether to use the new code path that performs ghosting in a safe way or not

    The new code path performing redistribution and ghosting in a safe way, i.e. such that ghosting
    is extended often and far enough, is not working for all contact scenarios, yet.
    Use this function to check, whether the scenario given in the input file can use the new path or
    has to stick to the old path (with bugs in the extension of the interface ghosting).

    @param[in] contactParams Parameter list with all contact-relevant input parameters
    @return True, if new path is chosen. False otherwise.
    */
    bool use_safe_redistribute_and_ghosting(const Teuchos::ParameterList& contactParams);

  }  // namespace Utils
}  // namespace CONTACT

FOUR_C_NAMESPACE_CLOSE

#endif
