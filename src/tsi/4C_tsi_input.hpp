// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_TSI_INPUT_HPP
#define FOUR_C_TSI_INPUT_HPP

#include "4C_config.hpp"

#include "4C_io_input_spec.hpp"

#include <optional>
#include <vector>


FOUR_C_NAMESPACE_OPEN

namespace TSI
{
  //! Type of coupling strategy for TSI problems
  enum class SolutionSchemeOverFields
  {
    OneWay,
    SequStagg,
    IterStagg,
    IterStaggAitken,
    IterStaggAitkenIrons,
    IterStaggFixedRel,
    Monolithic
  };

  //! Type of coupling variable for TSI problems
  enum class CouplingVariable
  {
    Displacement,
    Temperature
  };

  //! @name Solution technique and related stuff

  //! type of norm to check for convergence
  enum class ConvNorm
  {
    Abs,  //!< absolute norm
    Rel,  //!< relative norm of TSI problem with initial TSI rhs
    Mix   //!< mixed absolute-relative norm
  };

  //! type of norm to check for convergence
  enum class BinaryOp
  {
    bop_and,               //!< and
    bop_or,                //!< or
    bop_coupl_or_single,   //!< either TSI problem or single field problems converged
    bop_coupl_and_single,  //!< either TSI problem or single field problems converged
    bop_and_single,        //!< and in single field problems
    bop_or_single          //!< or in single field problems
  };

  //! type of solution techniques
  enum class NlnSolTech
  {
    fullnewton,  //!< full Newton-Raphson iteration
    ptc,         //!< pseudo transient continuation nonlinear iteration
  };

  //! type of line-search strategy
  enum class LineSearch
  {
    LS_none = 0,   //!< no line search
    LS_structure,  //!< line-search based on structural residual
    LS_thermo,     //!< line-search based on thermal residual
    LS_or,         //!< line-search based on structural or thermal residual
    LS_and,        //!< line-search based on structural and thermal residual
  };

  //! events that trigger monolithic TSI time step reduction
  enum class TimeStepReductionReason
  {
    none,
    material,
    nonlinear_solver,
    both,
  };

  //! settings for monolithic TSI time-step reduction and recovery
  struct TimeStepControlSettings
  {
    TimeStepReductionReason reduction_reason = TimeStepReductionReason::both;
    double decrease_factor = 0.5;
    double min_timestep = 1.0e-12;
    int steps_to_increase = 4;
    std::optional<double> max_timestep = std::nullopt;
    std::optional<double> increase_factor = std::nullopt;
  };

  //@}

  //! @name General
  //@{

  //! type of vector norm used for error/residual vectors
  enum class VectorNorm
  {
    L1,         //!< L1/linear norm
    L1_Scaled,  //!< L1/linear norm scaled by length of vector
    L2,         //!< L2/Euclidean norm
    Rms,        //!< root mean square (RMS) norm
    Inf         //!< Maximum/infinity norm
  };

  //! Method used to calculate plastic dissipation
  enum class DissipationMode
  {
    pl_multiplier,  //!< Dissipation = yield stress times plastic multiplier
    pl_flow,        //!< Dissipation = Mandel stress : sym(L^p)
    Taylor_Quinney  //!< Dissipation based on Taylor Quinney factor
  };

  //@}

  /// tsi parameters
  std::vector<Core::IO::InputSpec> valid_parameters();

}  // namespace TSI

/*----------------------------------------------------------------------*/
FOUR_C_NAMESPACE_CLOSE

#endif
