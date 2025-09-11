// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_INPAR_FS3I_HPP
#define FOUR_C_INPAR_FS3I_HPP

#include "4C_config.hpp"

#include "4C_io_input_spec.hpp"

#include <vector>


FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
namespace Inpar
{
  namespace FS3I
  {
    /// Type of coupling strategy for SSI problems
    enum SolutionSchemeOverFields
    {
      fs3i_SequStagg,
      //    fs3i_SequStagg_SolidToScatra,
      fs3i_IterStagg,
      //    fs3i_IterStaggFixedRel_ScatraToSolid,
      //    fs3i_IterStaggFixedRel_SolidToScatra,
      //    fs3i_IterStaggAitken_ScatraToSolid,
      //    fs3i_IterStaggAitken_SolidToScatra
      //    fs3i_IterStaggAitkenIrons,
      //    fs3i_Monolithic
    };

    /// Type of coupling strategy between the two fields of the SSI problems
    enum VolumeCoupling
    {
      coupling_match,
      coupling_nonmatch
    };

    /// fs3i parameters
    std::vector<Core::IO::InputSpec> valid_parameters();

  }  // namespace FS3I

}  // namespace Inpar

FOUR_C_NAMESPACE_CLOSE

#endif
