// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_config.hpp"

#include "4C_linalg_symmetric_tensor.hpp"
#include "4C_linalg_tensor.hpp"
#include "4C_linalg_tensor_generators.hpp"
#include "4C_mat_monolithic_solid_scalar_material.hpp"
#include "4C_mat_so3_material.hpp"
#include "4C_utils_exceptions.hpp"

#include <optional>

#ifndef FOUR_C_MAT_TRAIT_THERMO_SOLID_HPP
#define FOUR_C_MAT_TRAIT_THERMO_SOLID_HPP

FOUR_C_NAMESPACE_OPEN

namespace Mat
{

  /*!
   * @brief Mechanical heat source contribution produced by the material and the linearizations
   * needed by thermomechanical coupling. Default constructed with zero values.
   */
  struct HeatSource
  {
    //! heat produced per unit reference volume and time, i.e. positive values heat the material
    double value = 0.0;
    double derivative_wrt_temperature = 0.0;
    //! derivative w.r.t. the (Green-Lagrange) strain
    Core::LinAlg::SymmetricTensor<double, 3, 3> derivative_wrt_strain{};
    //! derivative w.r.t. the (Green-Lagrange) strain rate
    Core::LinAlg::SymmetricTensor<double, 3, 3> derivative_wrt_strain_rate{};
  };

  /*!
   * @brief Stress-temperature modulus, i.e. the partial derivative of the stress w.r.t. the
   * temperature at fixed strain and internal variables, and its derivatives. Default constructed
   * with zero values.
   */
  struct StressTemperatureModulus
  {
    Core::LinAlg::SymmetricTensor<double, 3, 3> value{};
    Core::LinAlg::SymmetricTensor<double, 3, 3> derivative_wrt_temperature{};
    //! derivative w.r.t. the strain
    Core::LinAlg::SymmetricTensor<double, 3, 3, 3, 3> derivative_wrt_strain{};
  };

  /*!
   * @brief Kinematic quantities needed to evaluate the mechanical heat source
   */
  struct KinematicState
  {
    //! strain, i.e. the linear strain for small and the Green-Lagrange strain for finite strains
    Core::LinAlg::SymmetricTensor<double, 3, 3> strain;
    //! rate of the strain
    Core::LinAlg::SymmetricTensor<double, 3, 3> strain_rate;
    //! deformation gradient, only available for finite strains
    std::optional<Core::LinAlg::Tensor<double, 3, 3>> defgrad;

    [[nodiscard]] static KinematicState from_linear_strain(
        const Core::LinAlg::SymmetricTensor<double, 3, 3>& strain,
        const Core::LinAlg::SymmetricTensor<double, 3, 3>& strain_rate)
    {
      return {.strain = strain, .strain_rate = strain_rate, .defgrad = std::nullopt};
    }

    [[nodiscard]] static KinematicState from_deformation_gradient(
        const Core::LinAlg::Tensor<double, 3, 3>& defgrad,
        const Core::LinAlg::Tensor<double, 3, 3>& defgrad_rate)
    {
      const Core::LinAlg::SymmetricTensor<double, 3, 3> strain =
          0.5 * (Core::LinAlg::assume_symmetry(Core::LinAlg::transpose(defgrad) * defgrad) -
                    Core::LinAlg::TensorGenerators::identity<double, 3, 3>);
      const auto strain_rate =
          0.5 * Core::LinAlg::assume_symmetry(Core::LinAlg::transpose(defgrad) * defgrad_rate +
                                              Core::LinAlg::transpose(defgrad_rate) * defgrad);
      return {.strain = strain, .strain_rate = strain_rate, .defgrad = defgrad};
    }

    [[nodiscard]] const Core::LinAlg::Tensor<double, 3, 3>& deformation_gradient() const
    {
      if (!defgrad.has_value())
        FOUR_C_THROW("The deformation gradient is only available for finite strains.");
      return *defgrad;
    }
  };

  namespace Trait
  {
    class ThermoSolid : public So3Material, public MonolithicSolidScalarMaterial
    {
     public:
      /*!
       * Evaluate the mechanical heat source, i.e., the thermoelastic heating
       * T . stm : dE/dt plus the additional heat sources of the material, and its derivatives.
       */
      [[nodiscard]] HeatSource evaluate_mechanical_heat_source(const double temperature,
          const KinematicState& kinematic_state, const EvaluationContext<3>& context, const int gp,
          const int eleGID)
      {
        const StressTemperatureModulus stress_temperature_modulus =
            evaluate_stress_temperature_modulus(temperature, kinematic_state, gp);

        const auto& strain_rate = kinematic_state.strain_rate;
        const double stress_temperature_modulus_power =
            Core::LinAlg::ddot(stress_temperature_modulus.value, strain_rate);
        const double stress_temperature_modulus_derivative_power =
            Core::LinAlg::ddot(stress_temperature_modulus.derivative_wrt_temperature, strain_rate);

        HeatSource source =
            evaluate_additional_heat_source(temperature, kinematic_state, context, gp, eleGID);
        source.value += temperature * stress_temperature_modulus_power;
        source.derivative_wrt_temperature +=
            stress_temperature_modulus_power +
            temperature * stress_temperature_modulus_derivative_power;
        source.derivative_wrt_strain_rate += temperature * stress_temperature_modulus.value;
        source.derivative_wrt_strain +=
            temperature *
            Core::LinAlg::ddot(strain_rate, stress_temperature_modulus.derivative_wrt_strain);
        return source;
      }

     protected:
      /*!
       * Evaluate the stress-temperature modulus and its derivatives w.r.t. temperature and strain.
       *
       * @param temperature current temperature
       * @param kinematic_state current material kinematics
       * @param gp Gauss-point index
       */
      [[nodiscard]] virtual StressTemperatureModulus evaluate_stress_temperature_modulus(
          double temperature, const KinematicState& kinematic_state, int gp) = 0;

      /*!
       * Evaluate heat sources in addition to the thermoelastic heating, e.g., plastic dissipation.
       *
       * @param temperature current temperature
       * @param kinematic_state current material kinematics
       * @param context evaluation context
       * @param gp Gauss-point index
       * @param eleGID global element id
       */
      [[nodiscard]] virtual HeatSource evaluate_additional_heat_source(const double temperature,
          const KinematicState& kinematic_state, const EvaluationContext<3>& context, const int gp,
          const int eleGID)
      {
        return {};
      }
    };
  }  // namespace Trait
}  // namespace Mat

FOUR_C_NAMESPACE_CLOSE

#endif