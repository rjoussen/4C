// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_MULTIPLICATIVE_SPLIT_DEFGRAD_ELASTHYPER_TEST_UTILS_HPP
#define FOUR_C_MULTIPLICATIVE_SPLIT_DEFGRAD_ELASTHYPER_TEST_UTILS_HPP

#include "4C_global_data.hpp"
#include "4C_io_input_field.hpp"
#include "4C_linalg_symmetric_tensor.hpp"
#include "4C_linalg_tensor.hpp"
#include "4C_linalg_tensor_einstein.hpp"
#include "4C_linalg_tensor_generators.hpp"
#include "4C_mat_material_factory.hpp"
#include "4C_mat_par_bundle.hpp"
#include "4C_mat_so3_material.hpp"
#include "4C_solid_ele_fibers.hpp"
#include "4C_ssi_input.hpp"
#include "4C_utils_exceptions.hpp"

#include <algorithm>
#include <cmath>
#include <functional>
#include <memory>
#include <optional>
#include <string>
#include <vector>

namespace FourC::MultiplicativeSplitDefgradElastHyperTestUtils
{
  /// elastic summands of the ElastHyper material used as elastic part of the multiplicative split
  enum class ElasticSummands
  {
    /// formulated in principal invariants
    coup_neo_hooke,
    /// formulated in modified invariants (isochoric/volumetric split)
    iso_neo_hooke_vol_sussman_bathe,
    /// adds an anisotropic summand
    coup_neo_hooke_transv_iso
  };

  /// material ids of the elastic part
  constexpr int matid_coup_neo_hooke = 1;
  constexpr int matid_iso_neo_hooke = 2;
  constexpr int matid_vol_sussman_bathe = 3;
  constexpr int matid_transv_iso = 4;
  constexpr int matid_structural_tensor_strategy = 5;
  constexpr int matid_elast_hyper = 6;
  /// first material id that is free for the remaining materials of a test
  constexpr int first_free_matid = 7;

  /// element fibers of the materials created by create_material(): the transversely isotropic
  /// summand is aligned with FIBER1, FIBER2 is not used here and free for a test, e.g., for a
  /// second anisotropy direction that differs from the one of the summand
  constexpr int elastic_fiber_id = 1;
  constexpr int free_fiber_id = 2;

  constexpr int gp = 0;
  constexpr int ele_gid = 0;

  /// parameters of a transversely isotropic summand along the element fiber @p fiber_id, also used
  /// for the fiber reader
  inline Core::IO::InputParameterContainer transv_iso_data(const double alpha, const double beta,
      const double gamma, const int structural_tensor_strategy_id, const int fiber_id)
  {
    Core::IO::InputParameterContainer data;
    data.add("ALPHA", alpha);
    data.add("BETA", beta);
    data.add("GAMMA", gamma);
    data.add("ANGLE", 0.0);
    data.add("STR_TENS_ID", structural_tensor_strategy_id);
    data.add("FIBER", fiber_id);
    data.add("INIT", 1);
    return data;
  }

  /// material ids of the given elastic summands
  inline std::vector<int> elastic_summand_ids(const ElasticSummands elastic_summands)
  {
    switch (elastic_summands)
    {
      case ElasticSummands::coup_neo_hooke:
        return {matid_coup_neo_hooke};
      case ElasticSummands::iso_neo_hooke_vol_sussman_bathe:
        return {matid_iso_neo_hooke, matid_vol_sussman_bathe};
      case ElasticSummands::coup_neo_hooke_transv_iso:
        return {matid_coup_neo_hooke, matid_transv_iso};
      default:
        FOUR_C_THROW("Unexpected elastic summands.");
    }
  }

  /// insert all elastic summands and an ElastHyper material with the given summands
  inline void set_up_elast_hyper(const ElasticSummands elastic_summands)
  {
    Core::IO::InputParameterContainer coup_neo_hooke_data;
    coup_neo_hooke_data.add("YOUNG", 2.0e5);
    coup_neo_hooke_data.add("NUE", 0.29);
    Global::Problem::instance()->materials()->insert(matid_coup_neo_hooke,
        Mat::make_parameter(matid_coup_neo_hooke, Core::Materials::MaterialType::mes_coupneohooke,
            coup_neo_hooke_data));

    // shear and bulk modulus of the coupled Neo-Hooke summand
    Core::IO::InputParameterContainer iso_neo_hooke_data;
    iso_neo_hooke_data.add("MUE", Core::IO::InputField<double>(7.75e4));
    Global::Problem::instance()->materials()->insert(matid_iso_neo_hooke,
        Mat::make_parameter(matid_iso_neo_hooke, Core::Materials::MaterialType::mes_isoneohooke,
            iso_neo_hooke_data));
    Core::IO::InputParameterContainer vol_sussman_bathe_data;
    vol_sussman_bathe_data.add("KAPPA", 1.59e5);
    Global::Problem::instance()->materials()->insert(matid_vol_sussman_bathe,
        Mat::make_parameter(matid_vol_sussman_bathe,
            Core::Materials::MaterialType::mes_volsussmanbathe, vol_sussman_bathe_data));

    Core::IO::InputParameterContainer structural_tensor_strategy_data;
    structural_tensor_strategy_data.add<std::string>("STRATEGY", "Standard");
    structural_tensor_strategy_data.add<std::string>("DISTR", "none");
    structural_tensor_strategy_data.add("C1", 1.0);
    structural_tensor_strategy_data.add("C2", 0.0);
    structural_tensor_strategy_data.add("C3", 0.0);
    structural_tensor_strategy_data.add("C4", 1.0e16);
    Global::Problem::instance()->materials()->insert(matid_structural_tensor_strategy,
        Mat::make_parameter(matid_structural_tensor_strategy,
            Core::Materials::MaterialType::mes_structuraltensorstratgy,
            structural_tensor_strategy_data));
    Global::Problem::instance()->materials()->insert(
        matid_transv_iso, Mat::make_parameter(matid_transv_iso,
                              Core::Materials::MaterialType::mes_couptransverselyisotropic,
                              transv_iso_data(1.0e4, 5.0e3, 2.0e4, matid_structural_tensor_strategy,
                                  elastic_fiber_id)));

    const std::vector<int> summand_ids = elastic_summand_ids(elastic_summands);
    Core::IO::InputParameterContainer elast_hyper_data;
    elast_hyper_data.add("NUMMAT", static_cast<int>(summand_ids.size()));
    elast_hyper_data.add("MATIDS", summand_ids);
    elast_hyper_data.add("DENS", 1.0);
    elast_hyper_data.add("POLYCONVEX", 0);
    Global::Problem::instance()->materials()->insert(
        matid_elast_hyper, Mat::make_parameter(matid_elast_hyper,
                               Core::Materials::MaterialType::m_elasthyper, elast_hyper_data));
  }

  /// insert the multiplicative split with id @p matid, the elastic part with the given elastic
  /// summands, the inelastic factors with the given ids and the given thermal expansion
  inline void set_up_multiplicative_split(const int matid, const ElasticSummands elastic_summands,
      const std::vector<int>& inelastic_factor_ids, const double ref_temperature,
      const double thermal_expansion_coefficient)
  {
    Core::IO::InputParameterContainer multiplicative_split_data;
    const std::vector<int> summand_ids = elastic_summand_ids(elastic_summands);
    multiplicative_split_data.add("NUMMATEL", static_cast<int>(summand_ids.size()));
    multiplicative_split_data.add("MATIDSEL", summand_ids);
    multiplicative_split_data.add("NUMFACINEL", static_cast<int>(inelastic_factor_ids.size()));
    multiplicative_split_data.add("INELDEFGRADFACIDS", inelastic_factor_ids);
    multiplicative_split_data.add("DENS", 1.0);
    multiplicative_split_data.add("REF_TEMPERATURE", ref_temperature);
    multiplicative_split_data.add("THERMAL_EXPANSION_COEFFICIENT", thermal_expansion_coefficient);
    Global::Problem::instance()->materials()->insert(
        matid, Mat::make_parameter(matid,
                   Core::Materials::MaterialType::m_multiplicative_split_defgrad_elasthyper,
                   multiplicative_split_data));
  }

  /// some generic parameters required for safety checks in the material setup
  inline void set_up_generic_parameters()
  {
    auto parameter_list = std::make_shared<Teuchos::ParameterList>();
    parameter_list->sublist("STRUCTURAL DYNAMIC", false).set("MASSLIN", Solid::MassLin::ml_none);
    parameter_list->sublist("SSI CONTROL")
        .set("COUPALGO", SSI::SolutionSchemeOverFields::ssi_IterStagg);
    Global::Problem::instance()->set_parameter_list(parameter_list);
  }

  /// create the material with id @p matid with a single Gauss point and two element fibers
  template <typename Material>
  std::shared_ptr<Material> create_material(const int matid)
  {
    auto material = std::dynamic_pointer_cast<Material>(Mat::factory(matid));
    FOUR_C_ASSERT_ALWAYS(material != nullptr, "Unexpected type of material {}.", matid);
    Discret::Elements::Fibers fibers;
    // FIBER1 (elastic_fiber_id)
    fibers.element_fibers.push_back(
        Core::LinAlg::Tensor<double, 3>{{1.0 / 3.0, 2.0 / 3.0, 2.0 / 3.0}});
    // FIBER2 (free_fiber_id)
    fibers.element_fibers.push_back(Core::LinAlg::Tensor<double, 3>{{0.0, 0.6, 0.8}});
    material->setup(1, fibers, std::nullopt);
    return material;
  }

  /// Green-Lagrange strain \f$ \boldsymbol{E} = \frac{1}{2} (\boldsymbol{F}^\mathrm{T}
  /// \boldsymbol{F} - \boldsymbol{I}) \f$
  inline Core::LinAlg::SymmetricTensor<double, 3, 3> green_lagrange(
      const Core::LinAlg::Tensor<double, 3, 3>& F)
  {
    return 0.5 * (Core::LinAlg::assume_symmetry(Core::LinAlg::dot(Core::LinAlg::transpose(F), F)) -
                     Core::LinAlg::TensorGenerators::identity<double, 3, 3>);
  }

  /// stress and stiffness tensor of a material evaluation
  struct StressResponse
  {
    Core::LinAlg::SymmetricTensor<double, 3, 3> stress;
    Core::LinAlg::SymmetricTensor<double, 3, 3, 3, 3> cmat;
  };

  /// evaluate the stress and stiffness tensor of @p material at @p F.
  /// The material is evaluated with the given @p params and @p context.
  inline StressResponse evaluate(Mat::So3Material& material,
      const Core::LinAlg::Tensor<double, 3, 3>& F, const Teuchos::ParameterList& params,
      const Mat::EvaluationContext<3>& context)
  {
    StressResponse response{};
    material.evaluate(
        &F, green_lagrange(F), params, context, response.stress, response.cmat, gp, ele_gid);
    return response;
  }

  /// largest absolute entry of @p tensor, e.g., to scale the tolerance of a comparison against
  /// finite differences, whose errors scale with the largest entry
  inline double max_abs(const auto& tensor)
  {
    return std::abs(std::ranges::max(
        tensor.container(), {}, [](const double value) { return std::abs(value); }));
  }

  /// Apply the chain rule to compute the derivative of the stress w.r.t. the deformation gradient
  /// from the given stiffness tensor
  /// \f$ \frac{\partial S_{ij}}{\partial F_{kl}} = \frac{1}{2} \mathbb{C}_{ijmn}
  /// \frac{\partial C_{mn}}{\partial F_{kl}} = \mathbb{C}_{ijml} F_{km} \f$
  inline Core::LinAlg::Tensor<double, 3, 3, 3, 3> compute_stress_deriv_wrt_defgrad(
      const Core::LinAlg::SymmetricTensor<double, 3, 3, 3, 3>& cmat,
      const Core::LinAlg::Tensor<double, 3, 3>& F)
  {
    return Core::LinAlg::einsum<"ijml", "km">(cmat, F);
  }

  /**
   * @brief Compute the derivative of the stress w.r.t. the deformation gradient with central finite
   * differences
   *
   * \f$ \frac{\partial S_{ij}}{\partial F_{kl}} \approx \frac{S_{ij}(\boldsymbol{F} + h
   * \boldsymbol{e}_k \otimes \boldsymbol{e}_l) - S_{ij}(\boldsymbol{F} - h \boldsymbol{e}_k
   * \otimes \boldsymbol{e}_l)}{2 h} \f$
   *
   * @param evaluate_stress function that computes the stress at a given deformation gradient
   * @param F deformation gradient at which the derivative is computed
   * @param h finite difference step size
   */
  inline Core::LinAlg::Tensor<double, 3, 3, 3, 3> finite_difference_stress_deriv_wrt_defgrad(
      const std::function<Core::LinAlg::SymmetricTensor<double, 3, 3>(
          const Core::LinAlg::Tensor<double, 3, 3>&)>& evaluate_stress,
      const Core::LinAlg::Tensor<double, 3, 3>& F, const double h = 1.0e-6)
  {
    Core::LinAlg::Tensor<double, 3, 3, 3, 3> dSdF_fd{};
    for (std::size_t k = 0; k < 3; ++k)
    {
      for (std::size_t l = 0; l < 3; ++l)
      {
        Core::LinAlg::Tensor<double, 3, 3> F_plus = F;
        Core::LinAlg::Tensor<double, 3, 3> F_minus = F;
        F_plus(k, l) += h;
        F_minus(k, l) -= h;

        const Core::LinAlg::SymmetricTensor<double, 3, 3> dS =
            (evaluate_stress(F_plus) - evaluate_stress(F_minus)) / (2.0 * h);
        for (std::size_t i = 0; i < 3; ++i)
          for (std::size_t j = 0; j < 3; ++j) dSdF_fd(i, j, k, l) = dS(i, j);
      }
    }
    return dSdF_fd;
  }
}  // namespace FourC::MultiplicativeSplitDefgradElastHyperTestUtils

#endif
