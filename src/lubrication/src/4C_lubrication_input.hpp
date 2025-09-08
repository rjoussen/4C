// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_LUBRICATION_INPUT_HPP
#define FOUR_C_LUBRICATION_INPUT_HPP

#include "4C_config.hpp"

#include "4C_io_input_spec.hpp"
#include "4C_utils_exceptions.hpp"
#include "4C_utils_parameter_list.fwd.hpp"


FOUR_C_NAMESPACE_OPEN

namespace Lubrication
{
  /// compute error compared to analytical solution
  enum CalcError
  {
    calcerror_no,
    calcerror_byfunction
  };

  /// compute velocity by function
  enum VelocityField
  {
    velocity_zero,
    velocity_function,
    velocity_EHL
  };

  /// compute height by function
  enum HeightField
  {
    height_zero,
    height_function,
    height_EHL
  };

  //! type of norm to check for convergence
  enum ConvNorm
  {
    convnorm_abs,  //!< absolute norm
    convnorm_rel,  //!< relative norm
    convnorm_mix   //!< mixed absolute-relative norm
  };

  //! type of vector norm used for error/residual vectors
  enum VectorNorm
  {
    norm_vague = 0,  //!< undetermined norm
    norm_l1,         //!< L1/linear norm
    norm_l2,         //!< L2/Euclidean norm
    norm_rms,        //!< root mean square (RMS) norm
    norm_inf         //!< Maximum/infinity norm
  };

  /// set the lubrication parameters
  Core::IO::InputSpec set_valid_parameters();

}  // namespace Lubrication

FOUR_C_NAMESPACE_CLOSE

#endif
