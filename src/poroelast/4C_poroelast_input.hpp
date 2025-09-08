// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_POROELAST_INPUT_HPP
#define FOUR_C_POROELAST_INPUT_HPP

#include "4C_config.hpp"

#include "4C_io_input_spec.hpp"
#include "4C_utils_exceptions.hpp"


FOUR_C_NAMESPACE_OPEN

namespace PoroElast
{
  /// Type of coupling strategy for poroelasticity problems
  enum class SolutionSchemeOverFields
  {
    //    OneWay,
    //   SequStagg,
    //   IterStagg,
    Partitioned,
    Monolithic,
    Monolithic_structuresplit,
    Monolithic_fluidsplit,
    Monolithic_nopenetrationsplit,
    Monolithic_meshtying
  };

  /// flag for control, in which equation the transient terms are included
  enum TransientEquationsOfPoroFluid
  {
    transient_none,
    transient_momentum_only,
    transient_continuity_only,
    transient_all
  };

  /// @name Solution technique and related

  /// type of norm to check for convergence
  enum ConvNorm
  {
    convnorm_undefined,
    convnorm_abs_global,       ///< absolute norm of global solution vectors
    convnorm_abs_singlefields  ///< absolute norm of single field solution vectors
    //    convnorm_rel_global,       ///< absolute norm of global solution vectors
    //    convnorm_rel_singlefields  ///< absolute norm of single field solution vectors
  };

  /// type of norm to be calculated
  enum VectorNorm
  {
    norm_undefined,
    norm_l1,         //!< L1/linear norm
    norm_l1_scaled,  //!< L1/linear norm scaled by length of vector
    norm_l2,         //!< L2/Euclidean norm
    norm_rms,        //!< root mean square (RMS) norm
    norm_inf         //!< Maximum/infinity norm
  };

  /// type of norm to check for convergence
  enum BinaryOp
  {
    bop_undefined,
    bop_and,  ///<  and
    bop_or    ///<  or
  };

  /// type of initial field for poroelasticity problem
  enum InitialField
  {
    //  initfield_zero_field,
    initfield_field_by_function
    //   initfield_field_by_condition
  };

  //@}


  /// set the poroelast parameters
  Core::IO::InputSpec set_valid_parameters();

}  // namespace PoroElast

FOUR_C_NAMESPACE_CLOSE

#endif
