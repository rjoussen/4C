// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_REDUCED_LUNG_TERMINAL_UNIT_RECRUITMENT_HPP
#define FOUR_C_REDUCED_LUNG_TERMINAL_UNIT_RECRUITMENT_HPP

#include "4C_config.hpp"

#include "4C_reduced_lung_terminal_unit_common.hpp"
#include "4C_utils_exceptions.hpp"

#include <utility>
#include <variant>
#include <vector>

FOUR_C_NAMESPACE_OPEN

namespace ReducedLung::TerminalUnits
{
  // Everything here is named after recruitment alone, but models both directions of the process:
  // terminal units recruit on the opening branch of the hysteresis loop and derecruit on the
  // closing one. The same holds for the input keys and the output fields.

  /**
   * @brief Enum types used in reduced-lung input for recruitment model selection.
   *
   * The recruitment data structs below are built on them, so unlike the rheology and elasticity
   * selection enums these aliases live next to the data rather than in the model namespace.
   */
  using PressureLawType =
      ReducedLungParameters::LungTree::TerminalUnits::RecruitmentModel::PressureLawType;
  using TimeLawType = ReducedLungParameters::LungTree::TerminalUnits::RecruitmentModel::TimeLawType;
  using HysteresisPath =
      ReducedLungParameters::LungTree::TerminalUnits::RecruitmentModel::HysteresisPath;
  using ReferenceVolumeLinearization = ReducedLungParameters::LungTree::TerminalUnits::
      RecruitmentModel::ReferenceVolumeLinearization;

  /**
   * @brief Terminal units whose reference volume is the geometry-derived constant.
   */
  struct NoRecruitment
  {
    ///< Reference volume of each element, derived from the geometry and never changed.
    std::vector<double> v0;
  };

  /**
   * @brief Element-wise relaxation of the reference volume towards its quasi-static target.
   *
   * The time law composes with whichever pressure law drives the target, so it is kept separate
   * from the pressure law parameters.
   */
  struct TimeLaw
  {
    ///< Time law delaying the reference volume response, applied on top of the pressure law.
    std::vector<TimeLawType> type;
    ///< Relaxation time of the exponential relaxation time law; 0 without that time law.
    std::vector<double> tau;
  };

  /**
   * @brief Piecewise-linear pressure law on an opening and a closing hysteresis branch.
   */
  struct LinearPressureLaw
  {
    ///< Hysteresis branch of the last converged time step, initialized from the input field.
    std::vector<HysteresisPath> active_path_n;

    ///< Reference volume of the fully derecruited and the fully recruited element.
    std::vector<double> v0_min;
    std::vector<double> v0_max;
    ///< Transpulmonary pressure at which recruitment sets in on the closing and the opening branch.
    std::vector<double> p_closing_min;
    std::vector<double> p_opening_min;
    ///< Pressure interval over which the reference volume sweeps from v0_min to v0_max.
    std::vector<double> delta_p_minmax;
    ///< Relative distance to v0_min/v0_max at which the hysteresis branch is switched.
    std::vector<double> epsilon_v0_switch;
  };

  /**
   * @brief Recruitment data for a block of terminal units on the piecewise-linear pressure law.
   *
   * Recruitment turns the reference volume into a state variable driven by the transpulmonary
   * pressure: the pressure law sets a quasi-static target, the time law relaxes towards it. The
   * state vectors hold the last converged time step, denoted by the suffix `_n`; the value of the
   * new time step only exists within update_recruitment_state().
   */
  struct LinearPressureRecruitment
  {
    LinearPressureLaw pressure_law;
    TimeLaw time_law;

    ///< Whether the reference volume derivative enters the Jacobian or is dropped.
    std::vector<ReferenceVolumeLinearization> reference_volume_linearization;

    ///< Reference volume of the last converged time step.
    std::vector<double> v0_n;
    ///< Quasi-static reference volume the time law relaxes towards. Output only; seeded with the
    ///< initial reference volume until the first time step advances it.
    std::vector<double> v0_target;
  };

  /**
   * @brief Variant containing all supported terminal-unit recruitment model data structs.
   */
  using RecruitmentModel = std::variant<NoRecruitment, LinearPressureRecruitment>;
}  // namespace ReducedLung::TerminalUnits

namespace ReducedLung::TerminalUnits::Recruitment
{
  /**
   * @brief Human-readable name for recruitment pressure-law enum values.
   */
  inline const char* pressure_law_name(const PressureLawType pressure_law_type)
  {
    switch (pressure_law_type)
    {
      case PressureLawType::None:
        return "None";
      case PressureLawType::LinearPressure:
        return "LinearPressure";
    }
    FOUR_C_THROW("Unknown recruitment pressure-law type enum value.");
  }

  /**
   * @brief Dispatch a recruitment pressure-law enum value to its concrete C++ model type.
   *
   * The callable must provide templated overloads via
   * `callable.template operator()<NoRecruitment>()` and
   * `callable.template operator()<LinearPressureRecruitment>()`.
   */
  template <typename Callable>
  void dispatch_pressure_law_type(const PressureLawType pressure_law_type, Callable&& callable)
  {
    switch (pressure_law_type)
    {
      case PressureLawType::None:
        std::forward<Callable>(callable).template operator()<NoRecruitment>();
        return;
      case PressureLawType::LinearPressure:
        std::forward<Callable>(callable).template operator()<LinearPressureRecruitment>();
        return;
    }
    FOUR_C_THROW("Unknown recruitment pressure-law type enum value.");
  }

  /**
   * @brief Reference volume of one element at the last converged time step.
   *
   * Blocks without a recruitment law report their constant geometry-derived value.
   */
  [[nodiscard]] double reference_volume_n(
      const RecruitmentModel& recruitment_model, size_t element_index);

  /**
   * @brief Build internal-state updater callback for the concrete recruitment variant.
   *
   * The returned callback refreshes TerminalUnitData::reference_volume_context of all elements in
   * one model block, so that elasticity and rheology can consume the reference volume without
   * knowing which recruitment law produced it. Everything reading the reference volume runs after
   * it. Under the frozen linearization the reference volume stays at the last converged value and
   * the derivatives vanish; under the coupled one it follows the recruitment law within the Newton
   * step, exactly as the end-of-timestep routine will advance it. Blocks without a recruitment law
   * report the constant geometry-derived reference volume.
   */
  InternalStateUpdater make_internal_state_updater(const RecruitmentModel& recruitment_model);

  /**
   * @brief Build end-of-timestep callback for the concrete recruitment variant.
   *
   * The returned callback advances reference volume and hysteresis branch of all elements in one
   * model block by one time step. Blocks without a recruitment law have no state to advance.
   */
  EndOfTimestepRoutine make_end_of_timestep_routine(RecruitmentModel& recruitment_model);

  /**
   * @brief Append element parameters and initialize the recruitment state vectors.
   */
  void append_model_parameters(RecruitmentModel& recruitment_model, int global_element_id,
      double geometry_volume,
      const ReducedLungParameters::LungTree::TerminalUnits::RecruitmentModel& parameters);

  /**
   * @brief Build output evaluator callback for the concrete recruitment variant.
   *
   * Blocks without a recruitment law contribute no output fields.
   */
  OutputEvaluator make_output_evaluator(const RecruitmentModel& recruitment_model);
}  // namespace ReducedLung::TerminalUnits::Recruitment

FOUR_C_NAMESPACE_CLOSE

#endif
