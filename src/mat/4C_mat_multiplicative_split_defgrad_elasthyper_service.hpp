// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_MAT_MULTIPLICATIVE_SPLIT_DEFGRAD_ELASTHYPER_SERVICE_HPP
#define FOUR_C_MAT_MULTIPLICATIVE_SPLIT_DEFGRAD_ELASTHYPER_SERVICE_HPP
#include "4C_config.hpp"

#include "4C_comm_pack_helpers.hpp"
#include "4C_linalg_fixedsizematrix.hpp"
#include "4C_linalg_fixedsizematrix_tensor_products.hpp"
#include "4C_linalg_fixedsizematrix_voigt_notation.hpp"
#include "4C_linalg_symmetric_tensor.hpp"
#include "4C_linalg_tensor.hpp"
#include "4C_linalg_tensor_conversion.hpp"
#include "4C_linalg_tensor_einstein.hpp"
#include "4C_linalg_tensor_generators.hpp"
#include "4C_linalg_utils_densematrix_determinant.hpp"
#include "4C_mat_elasthyper_service.hpp"
#include "4C_mat_service.hpp"
#include "4C_mat_so3_material.hpp"
#include "4C_utils_exceptions.hpp"

#include <Teuchos_ParameterList.hpp>

#include <cmath>
#include <limits>
#include <optional>
#include <vector>

FOUR_C_NAMESPACE_OPEN

namespace Mat
{
  /*!
   * @brief Mechanical heat source contribution produced by the material and the linearizations
   * needed by thermomechanical coupling. Default constructed with zero values.
   */
  struct HeatSource
  {
    double value = 0.0;
    //! derivative of the mechanical dissipation heat source w.r.t. temperature
    double derivative_wrt_temperature = 0.0;
    //! derivative of the mechanical dissipation heat source w.r.t. the right Cauchy-Green tensor
    //! (stress-form)
    Core::LinAlg::Matrix<1, 6> derivative_wrt_cauchy_green{Core::LinAlg::Initialization::zero};

    //! pack a vector of MechanicalDissipation
    static void pack(
        Core::Communication::PackBuffer& data, const std::vector<HeatSource>& heat_source)
    {
      Core::Communication::add_to_pack(data, static_cast<int>(heat_source.size()));

      for (const auto& md : heat_source)
      {
        Core::Communication::add_to_pack(data, md.value);
        Core::Communication::add_to_pack(data, md.derivative_wrt_temperature);
        Core::Communication::add_to_pack(data, md.derivative_wrt_cauchy_green);
      }
    }
  };


  /// thermal expansion stress \f$ \mathbf{S}_T(T) \f$ and its derivatives w.r.t. temperature. It
  /// only depends on the temperature, i.e., it can be evaluated once per temperature.
  struct ThermalExpansionStress
  {
    /// \f$ \mathbf{S}_T \f$
    Core::LinAlg::SymmetricTensor<double, 3, 3> value{};
    /// \f$ \partial \mathbf{S}_T / \partial T \f$
    Core::LinAlg::SymmetricTensor<double, 3, 3> temperature_derivative{};
    /// \f$ \partial^2 \mathbf{S}_T / \partial T^2 \f$, only evaluated if requested
    std::optional<Core::LinAlg::SymmetricTensor<double, 3, 3>> second_temperature_derivative;
  };

  /*!
   * @brief Evaluate the thermal expansion stress of a hyperelastic material due to isotropic
   * thermal expansion and its derivatives w.r.t. temperature
   *
   * The thermal expansion is described by the thermal right Cauchy-Green tensor
   * \f[
   *   \mathbf{C}_T = \left( 1 + 2 \alpha_T \Delta T \right) \mathbf{I}
   *   \quad \Leftrightarrow \quad
   *   \mathbf{E}_T = \alpha_T \Delta T \, \mathbf{I}
   * \f]
   * at which the hyperelastic material with the stress \f$ \mathbf{S}_\text{he} \f$ and the
   * stiffness \f$ \mathbb{C}_\text{he} = 2 \partial \mathbf{S}_\text{he} / \partial
   * \mathbf{C}_\text{he} \f$ is evaluated:
   * \f[
   *   \mathbf{S}_T = \mathbf{S}_\text{he}(\mathbf{C}_T), \qquad
   *   \frac{\partial \mathbf{S}_T}{\partial T}
   *   = \frac{1}{2} \mathbb{C}_\text{he}(\mathbf{C}_T) : \frac{\partial \mathbf{C}_T}{\partial T}
   *   = \alpha_T \, \mathbb{C}_\text{he}(\mathbf{C}_T) : \mathbf{I}
   * \f]
   * The stress only depends on \f$ \mathbf{E}_T \f$; the thermal deformation gradient
   * \f$ \mathbf{F}_T = \sqrt{1 + 2 \alpha_T \Delta T} \, \mathbf{I} \f$ is only used for the
   * polyconvexity check of MAT_ElastHyper. The hyperelastic material is evaluated
   * without parameters (only the Gauss point and the element ID enter, e.g., for fibers and
   * element-wise parameters).
   *
   * The second derivative
   * \f[
   *   \frac{\partial^2 \mathbf{S}_T}{\partial T^2}
   *   = 2 \alpha_T^2 \frac{\partial \mathbb{C}_\text{he}}{\partial
   * \mathbf{C}_\text{he}}(\mathbf{C}_T) : \mathbf{I} : \mathbf{I}
   * \f]
   * requires the derivative of the elasticity tensor, which the hyperelastic material does not
   * provide. Hence, it is approximated by central differences of the analytic first derivative. The
   * step size is chosen relative to the thermal stretch \f$ \lambda_T^2 = 1 + 2 \alpha_T \Delta T
   * \f$, whose natural scale is one, independent of the units of the temperature.
   *
   * @param[in] material                          Hyperelastic material
   * @param[in] thermal_expansion_coefficient     \f$ \alpha_T \f$ (no thermal stress if zero)
   * @param[in] delta_temperature                 \f$ \Delta T \f$, temperature difference to the
   *                                              reference temperature
   * @param[in] context                           Evaluation context
   * @param[in] gp                                Gauss point
   * @param[in] eleGID                            Element ID
   * @param[in] max_temperature_derivative_order  Highest derivative w.r.t. temperature to
   *                                              evaluate (1 or 2)
   */
  inline ThermalExpansionStress evaluate_thermal_expansion_stress(Mat::So3Material& material,
      const double thermal_expansion_coefficient, const double delta_temperature,
      const EvaluationContext<3>& context, const int gp, const int eleGID,
      const int max_temperature_derivative_order = 1)
  {
    FOUR_C_ASSERT(max_temperature_derivative_order == 1 || max_temperature_derivative_order == 2,
        "Only the first and second derivative w.r.t. temperature are available.");

    ThermalExpansionStress thermal_stress{};
    if (max_temperature_derivative_order == 2)
      thermal_stress.second_temperature_derivative = Core::LinAlg::SymmetricTensor<double, 3, 3>{};
    if (thermal_expansion_coefficient == 0.0) return thermal_stress;

    const auto& id = Core::LinAlg::TensorGenerators::identity<double, 3, 3>;
    const double alpha = thermal_expansion_coefficient;

    // S_T and dS_T/dT at the given temperature difference
    const auto evaluate_value_and_temperature_derivative = [&](const double temperature_difference)
    {
      const double thermal_stretch_squared = 1.0 + 2.0 * alpha * temperature_difference;
      FOUR_C_ASSERT(thermal_stretch_squared > 0.0,
          "Thermal expansion must result in a positive thermal stretch, got lambda_T^2 = {}.",
          thermal_stretch_squared);

      // thermal Green-Lagrange strain E_T and deformation gradient F_T
      const Core::LinAlg::SymmetricTensor<double, 3, 3> ET = alpha * temperature_difference * id;
      const Core::LinAlg::Tensor<double, 3, 3> FT =
          std::sqrt(thermal_stretch_squared) * Core::LinAlg::get_full(id);

      // S_T and the hyperelastic stiffness at C_T
      ThermalExpansionStress result{};
      Core::LinAlg::SymmetricTensor<double, 3, 3, 3, 3> cmat_T{};
      material.evaluate(&FT, ET, {}, context, result.value, cmat_T, gp, eleGID);

      // dS_T/dT = alpha_T cmat_he(C_T) : I
      result.temperature_derivative = alpha * Core::LinAlg::ddot(cmat_T, id);
      return result;
    };

    const ThermalExpansionStress at_temperature =
        evaluate_value_and_temperature_derivative(delta_temperature);
    thermal_stress.value = at_temperature.value;
    thermal_stress.temperature_derivative = at_temperature.temperature_derivative;

    if (max_temperature_derivative_order == 2)
    {
      const double thermal_stretch_squared = 1.0 + 2.0 * alpha * delta_temperature;
      const double step = std::cbrt(std::numeric_limits<double>::epsilon()) *
                          thermal_stretch_squared / (2.0 * std::abs(alpha));
      *thermal_stress.second_temperature_derivative =
          (evaluate_value_and_temperature_derivative(delta_temperature + step)
                  .temperature_derivative -
              evaluate_value_and_temperature_derivative(delta_temperature - step)
                  .temperature_derivative) /
          (2.0 * step);
    }
    return thermal_stress;
  }

  /// stress and its derivatives at fixed inelastic deformation gradient, obtained by pulling back
  /// the response of the elastic material (see pull_back_stress_and_derivatives())
  struct StressAndDerivatives
  {
    /// complete 2nd Piola-Kirchhoff stress \f$ \mathbf{S} \f$
    Core::LinAlg::SymmetricTensor<double, 3, 3> stress{};
    /// \f$ 2 \frac{\partial \mathbf{S}}{\partial \mathbf{C}} \f$ at fixed \f$
    /// \mathbf{F}_\text{in}^{-1} \f$
    Core::LinAlg::SymmetricTensor<double, 3, 3, 3, 3> cmat_at_fixed_iFin{};
    /// \f$ \frac{\partial \mathbf{S}}{\partial \mathbf{F}_\text{in}^{-1}} \f$, where the first
    /// index is in stress-like Voigt notation.
    Core::LinAlg::Matrix<6, 9> dstress_diFin{Core::LinAlg::Initialization::zero};
    /// \f$ \frac{\partial \mathbf{S}}{\partial T} \f$ at fixed \f$ \mathbf{F}_\text{in}^{-1} \f$.
    Core::LinAlg::SymmetricTensor<double, 3, 3> dstress_dT{};
  };

  /*!
   * @brief Pull back the stress of the intermediate configuration and its derivatives to the
   * reference configuration, while keeping the inelastic deformation gradient fixed.
   *
   * Given the stress \f$ \mathbf{S}_\theta \f$ of the intermediate configuration at the elastic
   * right Cauchy-Green tensor \f$ \mathbf{C}_\text{e} = \mathbf{F}_\text{in}^{-T} \cdot \mathbf{C}
   * \cdot \mathbf{F}_\text{in}^{-1} \f$ and its stiffness \f$ \mathbb{C}_\text{e} = 2 \partial
   * \mathbf{S}_\theta / \partial \mathbf{C}_\text{e} \f$, the complete 2nd Piola-Kirchhoff stress
   * is
   * \f[
   *   \mathbf{S} = \det(\mathbf{F}_\text{in}) \mathbf{F}_\text{in}^{-1} \cdot \mathbf{S}_\theta
   *   \cdot \mathbf{F}_\text{in}^{-T}
   * \f]
   * and its temperature derivative at fixed \f$ \mathbf{F}_\text{in}^{-1} \f$ is
   * \f[
   * \frac{\partial \mathbf{S}}{\partial T} = \det(\mathbf{F}_\text{in})
   *   \mathbf{F}_\text{in}^{-1} \cdot \frac{\partial \mathbf{S}_\theta}{\partial T}
   *   \cdot \mathbf{F}_\text{in}^{-T}
   * \f]
   * The tangent stiffness tensor at fixed \f$ \mathbf{F}_\text{in}^{-1} \f$ is given by the
   * stiffness tensor pulled back in its first index pair and applying the chain rule:
   * \f[
   *   \mathbb{C}_{ijkl} = 2 \frac{\partial S_{ij}}{\partial C_{kl}}
   *    = \det(\mathbf{F}_\text{in}) \underbrace{F^{-1}_{\text{in},ia}
   *   F^{-1}_{\text{in},jb} \mathbb{C}_{\text{e},abcd} }_{=: \mathbb{X}_{ijcd}}
   *   \frac{\partial C_{\text{e},cd}}{\partial C_{kl}}
   * \f]
   * and the derivative of the stress wrt. \f$ \mathbf{F}_\text{in}^{-1} \f$ is given by
   * \f[
   *   \frac{\partial S_{ij}}{\partial F^{-1}_{\text{in},kl}}
   *   = \det(\mathbf{F}_\text{in}) \left(
   *     \delta_{ik} \left[ \mathbf{F}_\text{in}^{-1} \cdot \mathbf{S}_\theta \right]_{jl}
   *   + \delta_{jk} \left[ \mathbf{F}_\text{in}^{-1} \cdot \mathbf{S}_\theta \right]_{il} \right)
   *   + \det(\mathbf{F}_\text{in}) \frac{1}{2}  \mathbb{X}_{ijab}
   *     \frac{\partial C_{\text{e},ab}}{\partial F^{-1}_{\text{in},kl}}
   *   - S_{ij} F^T_{\text{in},kl}
   * \f]
   * where the last term stems from the derivative of the determinant.
   *
   * The derivatives of the elastic Cauchy-Green tensor \f$ \mathbf{C}_\text{e} \f$ w.r.t. the
   * right Cauchy-Green tensor \f$ \mathbf{C} \f$ and the inverse inelastic deformation gradient
   * \f$ \mathbf{F}_\text{in}^{-1} \f$ are given by
   * \f[ \frac{\partial C_{\text{e},cd}}{\partial C_{kl}} = F^{-1}_{\text{in},kc}
   *   F^{-1}_{\text{in},ld} \f]
   *
   * and
   * \f[ \frac{\partial C_{\text{e},ab}}{\partial F^{-1}_{\text{in},kl}}
   *   = \delta_{al} \left[ \mathbf{C} \cdot \mathbf{F}_\text{in}^{-1} \right]_{kb}
   *   + \delta_{bl} \left[ \mathbf{C} \cdot \mathbf{F}_\text{in}^{-1} \right]_{ka} \f]
   *
   * @param[in] defgrad     Deformation gradient \f$ \mathbf{F} \f$
   * @param[in] iFin        Inverse inelastic deformation gradient \f$ \mathbf{F}_\text{in}^{-1} \f$
   * @param[in] Stheta      Stress of the intermediate configuration \f$ \mathbf{S}_\theta \f$
   * @param[in] cmat_e      Stiffness of the intermediate configuration \f$ \mathbb{C}_\text{e} \f$
   * @param[in] dStheta_dT  Temperature derivative \f$ \partial \mathbf{S}_\theta / \partial T \f$
   */
  inline StressAndDerivatives pull_back_stress_and_derivatives(
      const Core::LinAlg::Tensor<double, 3, 3>& defgrad,
      const Core::LinAlg::Tensor<double, 3, 3>& iFin,
      const Core::LinAlg::SymmetricTensor<double, 3, 3>& Stheta,
      const Core::LinAlg::SymmetricTensor<double, 3, 3, 3, 3>& cmat_e,
      const Core::LinAlg::SymmetricTensor<double, 3, 3>& dStheta_dT)
  {
    const auto& id = Core::LinAlg::TensorGenerators::identity<double, 3, 3>;
    const double detFin = 1.0 / Core::LinAlg::det(iFin);
    const Core::LinAlg::Tensor<double, 3, 3> iFinT = Core::LinAlg::transpose(iFin);

    StressAndDerivatives stress_and_derivatives;

    stress_and_derivatives.stress = detFin * Core::LinAlg::assume_symmetry(Core::LinAlg::dot(
                                                 Core::LinAlg::dot(iFin, Stheta), iFinT));
    stress_and_derivatives.dstress_dT = detFin * Core::LinAlg::assume_symmetry(Core::LinAlg::dot(
                                                     Core::LinAlg::dot(iFin, dStheta_dT), iFinT));

    // pull back the cmat_e, one index at a time (cheaper than pulling back all/multiple indices at
    // once)
    const Core::LinAlg::Tensor<double, 3, 3, 3, 3> cmat_e_index_1_pulled_back =
        Core::LinAlg::einsum<"ia", "ajkl">(iFin, cmat_e);
    const Core::LinAlg::Tensor<double, 3, 3, 3, 3> cmat_e_index_12_pulled_back =
        Core::LinAlg::einsum<"jb", "ibkl">(iFin, cmat_e_index_1_pulled_back);

    const Core::LinAlg::Tensor<double, 3, 3, 3, 3> cmat_e_index_123_pulled_back =
        Core::LinAlg::einsum<"kc", "ijcl">(iFin, cmat_e_index_12_pulled_back);

    // all indices pulled back, now multiply with det(F_in) and assume symmetry
    stress_and_derivatives.cmat_at_fixed_iFin =
        detFin * Core::LinAlg::assume_symmetry(
                     Core::LinAlg::einsum<"ld", "ijkd">(iFin, cmat_e_index_123_pulled_back));

    // helpers
    const Core::LinAlg::Tensor<double, 3, 3> iFin_Stheta = Core::LinAlg::dot(iFin, Stheta);
    const Core::LinAlg::Tensor<double, 3, 3> C_iFin =
        Core::LinAlg::dot(Core::LinAlg::dot(Core::LinAlg::transpose(defgrad), defgrad), iFin);

    // assemble dS_diFin (note that cmat_e_index_12_pulled_back is symmetric in the last two
    // indices, which allows for simplification)
    const Core::LinAlg::Tensor<double, 3, 3, 3, 3> dS_diFin =
        detFin * (Core::LinAlg::einsum<"ik", "jl">(id, iFin_Stheta) +
                     Core::LinAlg::einsum<"jk", "il">(id, iFin_Stheta) +
                     Core::LinAlg::einsum<"ijlb", "kb">(cmat_e_index_12_pulled_back, C_iFin)) -
        Core::LinAlg::einsum<"ij", "kl">(stress_and_derivatives.stress, Core::LinAlg::inv(iFinT));
    stress_and_derivatives.dstress_diFin =
        Core::LinAlg::make_6x9_voigt_matrix_from_tensor(dS_diFin);

    return stress_and_derivatives;
  }

  inline void evaluate_ce(const Core::LinAlg::Matrix<3, 3>& F,
      const Core::LinAlg::Matrix<3, 3>& iFin, Core::LinAlg::Matrix<3, 3>& Ce)
  {
    static Core::LinAlg::Matrix<3, 3> FiFin(Core::LinAlg::Initialization::uninitialized);
    FiFin.multiply_nn(F, iFin);
    Ce.multiply_tn(FiFin, FiFin);
  }

  inline void evaluatei_cin_ci_cin(const Core::LinAlg::Matrix<3, 3>& C,
      const Core::LinAlg::Matrix<3, 3>& iCin, Core::LinAlg::Matrix<3, 3>& iCinCiCin)
  {
    static Core::LinAlg::Matrix<3, 3> CiCin(Core::LinAlg::Initialization::uninitialized);
    CiCin.multiply_nn(C, iCin);
    iCinCiCin.multiply_nn(iCin, CiCin);
  }

  inline void elast_hyper_evaluate_elastic_part(const Core::LinAlg::Matrix<3, 3>& F,
      const Core::LinAlg::Matrix<3, 3>& iFin, Core::LinAlg::Matrix<6, 1>& S_stress,
      Core::LinAlg::Matrix<6, 6>& cmat,
      const std::vector<std::shared_ptr<Mat::Elastic::Summand>>& potsum,
      Mat::SummandProperties summandProperties, const int gp, const int eleGID)
  {
    if (summandProperties.anisomod or summandProperties.anisoprinc)
    {
      FOUR_C_THROW(
          "An additional inelastic part is not yet implemented for anisotropic materials.");
    }

    S_stress.clear();
    cmat.clear();

    // Variables needed for the computation of the stress resultants
    static Core::LinAlg::Matrix<3, 3> C(Core::LinAlg::Initialization::zero);
    static Core::LinAlg::Matrix<3, 3> Ce(Core::LinAlg::Initialization::zero);
    static Core::LinAlg::Matrix<3, 3> iC(Core::LinAlg::Initialization::zero);
    static Core::LinAlg::Matrix<3, 3> iCin(Core::LinAlg::Initialization::zero);
    static Core::LinAlg::Matrix<3, 3> iCinCiCin(Core::LinAlg::Initialization::zero);

    static Core::LinAlg::Matrix<6, 1> iCinv(Core::LinAlg::Initialization::zero);
    static Core::LinAlg::Matrix<6, 1> iCinCiCinv(Core::LinAlg::Initialization::zero);
    static Core::LinAlg::Matrix<6, 1> iCv(Core::LinAlg::Initialization::zero);
    static Core::LinAlg::Matrix<3, 1> principleInvariantsCe(Core::LinAlg::Initialization::zero);

    // Compute right Cauchy-Green tensor C=F^TF
    C.multiply_tn(F, F);

    // Compute inverse right Cauchy-Green tensor C^-1
    iC.invert(C);

    // Compute inverse inelastic right Cauchy-Green Tensor
    iCin.multiply_nt(iFin, iFin);

    // Compute iCin * C * iCin
    Mat::evaluatei_cin_ci_cin(C, iCin, iCinCiCin);

    // Compute Ce
    Mat::evaluate_ce(F, iFin, Ce);

    // Compute principal invariants
    Mat::invariants_principal(principleInvariantsCe, Ce);

    Core::LinAlg::Matrix<3, 1> dPIe(Core::LinAlg::Initialization::zero);
    Core::LinAlg::Matrix<6, 1> ddPIIe(Core::LinAlg::Initialization::zero);

    Mat::elast_hyper_evaluate_invariant_derivatives(
        principleInvariantsCe, dPIe, ddPIIe, potsum, summandProperties, gp, eleGID);

    // 2nd Piola Kirchhoff stress factors (according to Holzapfel-Nonlinear Solid Mechanics p. 216)
    static Core::LinAlg::Matrix<3, 1> gamma(Core::LinAlg::Initialization::zero);
    // constitutive tensor factors (according to Holzapfel-Nonlinear Solid Mechanics p. 261)
    static Core::LinAlg::Matrix<8, 1> delta(Core::LinAlg::Initialization::zero);

    Mat::calculate_gamma_delta(gamma, delta, principleInvariantsCe, dPIe, ddPIIe);

    // Convert necessary tensors to stress-like Voigt-Notation
    Core::LinAlg::Voigt::Stresses::matrix_to_vector(iCin, iCinv);
    Core::LinAlg::Voigt::Stresses::matrix_to_vector(iCinCiCin, iCinCiCinv);
    Core::LinAlg::Voigt::Stresses::matrix_to_vector(iC, iCv);

    // Contribution to 2nd Piola-Kirchhoff stress tensor
    S_stress.update(gamma(0), iCinv, 1.0);
    S_stress.update(gamma(1), iCinCiCinv, 1.0);
    S_stress.update(gamma(2), iCv, 1.0);

    // Contribution to the linearization
    cmat.multiply_nt(delta(0), iCinv, iCinv, 1.);
    cmat.multiply_nt(delta(1), iCinCiCinv, iCinv, 1.);
    cmat.multiply_nt(delta(1), iCinv, iCinCiCinv, 1.);
    cmat.multiply_nt(delta(2), iCinv, iCv, 1.);
    cmat.multiply_nt(delta(2), iCv, iCinv, 1.);
    cmat.multiply_nt(delta(3), iCinCiCinv, iCinCiCinv, 1.);
    cmat.multiply_nt(delta(4), iCinCiCinv, iCv, 1.);
    cmat.multiply_nt(delta(4), iCv, iCinCiCinv, 1.);
    cmat.multiply_nt(delta(5), iCv, iCv, 1.);
    Core::LinAlg::FourTensorOperations::add_holzapfel_product(cmat, iCv, delta(6));
    Core::LinAlg::FourTensorOperations::add_holzapfel_product(cmat, iCinv, delta(7));
  }

}  // namespace Mat
FOUR_C_NAMESPACE_CLOSE

#endif
