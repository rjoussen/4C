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
#include "4C_linalg_tensor_generators.hpp"
#include "4C_linalg_utils_densematrix_determinant.hpp"
#include "4C_mat_elasthyper_service.hpp"
#include "4C_mat_service.hpp"
#include "4C_mat_so3_material.hpp"
#include "4C_utils_exceptions.hpp"

#include <Teuchos_ParameterList.hpp>

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


  /// helper to store a second order tensor and its derivative w.r.t. temperature
  struct TensorAndTemperatureDerivative
  {
    /// value of the second order tensor
    Core::LinAlg::SymmetricTensor<double, 3, 3> value{};
    /// derivative of the second order tensor w.r.t. temperature
    Core::LinAlg::SymmetricTensor<double, 3, 3> temperature_derivative{};
  };

  /*!
   * @brief Evaluate the thermo-elastic stress of a hyperelastic material due to isotropic thermal
   * expansion and its derivative w.r.t. temperature
   *
   * The thermal expansion is described by the thermal right Cauchy-Green tensor
   * \f[
   *   \mathbf{C}_T = \left( 1 + 2 \alpha_T \Delta T \right) \mathbf{I}
   *   \quad \Leftrightarrow \quad
   *   \mathbf{E}_T = \alpha_T \Delta T \, \mathbf{I}
   * \f]
   * at which the hyperelastic material is evaluated:
   * \f[
   *   \mathbf{S}_T = \mathbf{S}_\text{e}(\mathbf{C}_T), \qquad
   *   \frac{\partial \mathbf{S}_T}{\partial T}
   *   = \frac{1}{2} \mathbb{C}_\text{e}(\mathbf{C}_T) : \frac{\partial \mathbf{C}_T}{\partial T}
   *   = \alpha_T \, \mathbb{C}_\text{e}(\mathbf{C}_T) : \mathbf{I}
   * \f]
   * The material only depends on \f$ \mathbf{E}_T \f$; the thermal deformation gradient
   * \f$ \mathbf{F}_T = \sqrt{1 + 2 \alpha_T \Delta T} \, \mathbf{I} \f$ is only passed since the
   * So3Material interface requires it.
   *
   * @param[in] material                       Hyperelastic material
   * @param[in] thermal_expansion_coefficient  \f$ \alpha_T \f$ (no thermal stress if zero)
   * @param[in] delta_temperature              \f$ \Delta T \f$, temperature difference to the
   *                                           reference temperature
   * @param[in] params                         Parameter list as handed in from the element
   * @param[in] context                        Evaluation context
   * @param[in] gp                             Gauss point
   * @param[in] eleGID                         Element ID
   * @return \f$ \mathbf{S}_T \f$ and \f$ \partial \mathbf{S}_T / \partial T \f$
   */
  inline TensorAndTemperatureDerivative evaluate_thermal_expansion_stress(
      Mat::So3Material& material, const double thermal_expansion_coefficient,
      const double delta_temperature, const Teuchos::ParameterList& params,
      const EvaluationContext<3>& context, const int gp, const int eleGID)
  {
    TensorAndTemperatureDerivative thermal_stress{};
    if (thermal_expansion_coefficient == 0.0) return thermal_stress;

    const auto& id = Core::LinAlg::TensorGenerators::identity<double, 3, 3>;
    const double thermal_stretch_squared =
        1.0 + 2.0 * thermal_expansion_coefficient * delta_temperature;
    FOUR_C_ASSERT(thermal_stretch_squared > 0.0,
        "Thermal expansion must result in a positive thermal stretch, got lambda_T^2 = {}.",
        thermal_stretch_squared);

    // thermal Green-Lagrange strain E_T and deformation gradient F_T
    const Core::LinAlg::SymmetricTensor<double, 3, 3> ET =
        thermal_expansion_coefficient * delta_temperature * id;
    const Core::LinAlg::Tensor<double, 3, 3> FT =
        std::sqrt(thermal_stretch_squared) * Core::LinAlg::get_full(id);

    // S_T and cmat_e at C_T
    Core::LinAlg::SymmetricTensor<double, 3, 3, 3, 3> cmat_T{};
    material.evaluate(&FT, ET, params, context, thermal_stress.value, cmat_T, gp, eleGID);

    // dS_T/dT = alpha_T cmat_e(C_T) : I
    thermal_stress.temperature_derivative =
        thermal_expansion_coefficient * Core::LinAlg::ddot(cmat_T, id);
    return thermal_stress;
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
