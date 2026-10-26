// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_mat_multiplicative_split_defgrad_elasthyper.hpp"

#include "4C_comm_pack_helpers.hpp"
#include "4C_global_data.hpp"
#include "4C_legacy_enum_definitions_materials.hpp"
#include "4C_linalg_fixedsizematrix.hpp"
#include "4C_linalg_fixedsizematrix_generators.hpp"
#include "4C_linalg_fixedsizematrix_tensor_products.hpp"
#include "4C_linalg_fixedsizematrix_voigt_notation.hpp"
#include "4C_linalg_symmetric_tensor.hpp"
#include "4C_linalg_tensor.hpp"
#include "4C_linalg_tensor_conversion.hpp"
#include "4C_linalg_tensor_einstein.hpp"
#include "4C_linalg_tensor_generators.hpp"
#include "4C_mat_inelastic_defgrad_factors.hpp"
#include "4C_mat_material_factory.hpp"
#include "4C_mat_multiplicative_split_defgrad_elasthyper_service.hpp"
#include "4C_mat_par_bundle.hpp"
#include "4C_mat_service.hpp"
#include "4C_mat_so3_material.hpp"
#include "4C_ssi_input.hpp"
#include "4C_structure_new_enum_lists.hpp"
#include "4C_utils_enum.hpp"
#include "4C_utils_exceptions.hpp"

#include <Teuchos_StandardParameterEntryValidators.hpp>

#include <memory>

FOUR_C_NAMESPACE_OPEN

namespace
{
  std::shared_ptr<Mat::So3Material> create_elastic_material(const int matid)
  {
    auto elastic_material = std::dynamic_pointer_cast<Mat::So3Material>(Mat::factory(matid));
    FOUR_C_ASSERT_ALWAYS(elastic_material != nullptr &&
                             elastic_material->material_type() == Core::Materials::m_elasthyper,
        "The elastic material (MATID_ELAST = {}) of MAT_MultiplicativeSplitDefgradElastHyper must "
        "be of type MAT_ElastHyper.",
        matid);
    return elastic_material;
  }
}  // namespace


/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
Mat::PAR::MultiplicativeSplitDefgradElastHyper::MultiplicativeSplitDefgradElastHyper(
    const Core::Mat::PAR::Parameter::Data& matdata)
    : Parameter(matdata),
      matid_elast_(matdata.parameters.get<int>("MATID_ELAST")),
      numfac_inel_(matdata.parameters.get<int>("NUMFACINEL")),
      inel_defgradfacids_(matdata.parameters.get<std::vector<int>>("INELDEFGRADFACIDS")),
      ref_temperature_(matdata.parameters.get<double>("REF_TEMPERATURE")),
      thermal_expansion_coefficient_(
          matdata.parameters.get<double>("THERMAL_EXPANSION_COEFFICIENT"))
{
  if (numfac_inel_ != static_cast<int>(inel_defgradfacids_.size()))
  {
    FOUR_C_THROW(
        "number of inelastic deformation gradient factors {} does not fit to size of inelastic "
        "deformation gradient ID vector {}",
        numfac_inel_, inel_defgradfacids_.size());
  }
}

std::shared_ptr<Core::Mat::Material>
Mat::PAR::MultiplicativeSplitDefgradElastHyper::create_material()
{
  return std::make_shared<Mat::MultiplicativeSplitDefgradElastHyper>(this);
}

Mat::MultiplicativeSplitDefgradElastHyperType
    Mat::MultiplicativeSplitDefgradElastHyperType::instance_;

Core::Communication::ParObject* Mat::MultiplicativeSplitDefgradElastHyperType::create(
    Core::Communication::UnpackBuffer& buffer)
{
  auto* splitdefgrad_elhy = new Mat::MultiplicativeSplitDefgradElastHyper();
  splitdefgrad_elhy->unpack(buffer);

  return splitdefgrad_elhy;
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
Mat::MultiplicativeSplitDefgradElastHyper::MultiplicativeSplitDefgradElastHyper()
    : inelastic_(std::make_shared<Mat::InelasticFactorsHandler>()),
      params_(nullptr),
      elastic_material_(nullptr)
{
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
Mat::MultiplicativeSplitDefgradElastHyper::MultiplicativeSplitDefgradElastHyper(
    Mat::PAR::MultiplicativeSplitDefgradElastHyper* params)
    : inelastic_(std::make_shared<Mat::InelasticFactorsHandler>()),
      params_(params),
      elastic_material_(create_elastic_material(params->matid_elast_))
{
  inelastic_->assign_to_source(params, *elastic_material_);
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
void Mat::MultiplicativeSplitDefgradElastHyper::pack(Core::Communication::PackBuffer& data) const
{
  // pack type of this instance of ParObject
  int type = unique_par_object_id();
  add_to_pack(data, type);
  // matid
  int matid = -1;
  if (params_ != nullptr) matid = params_->id();  // in case we are in post-process mode
  add_to_pack(data, matid);

  Core::Communication::PotentiallyUnusedBufferScope material_scope{data};
  if (params_ != nullptr)  // materials are not accessible in postprocessing mode
  {
    // pack elastic material first, the inelastic factors need it on construction when unpacking
    elastic_material_->pack(data);

    // pack inelastic factors handler
    inelastic_->pack_inelastic(data);
  }
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
void Mat::MultiplicativeSplitDefgradElastHyper::unpack(Core::Communication::UnpackBuffer& buffer)
{
  // make sure we have a pristine material
  params_ = nullptr;
  elastic_material_ = nullptr;

  Core::Communication::extract_and_assert_id(buffer, unique_par_object_id());

  // matid and recover params_
  int matid;
  extract_from_pack(buffer, matid);
  if (Global::Problem::instance()->materials() != nullptr)
  {
    if (Global::Problem::instance()->materials()->num() != 0)
    {
      const int probinst = Global::Problem::instance()->materials()->get_read_from_problem();
      auto* mat = Global::Problem::instance(probinst)->materials()->parameter_by_id(matid);
      if (mat->type() == material_type())
        params_ = dynamic_cast<Mat::PAR::MultiplicativeSplitDefgradElastHyper*>(mat);
      else
        FOUR_C_THROW("Type of parameter material {} does not fit to calling type {}", mat->type(),
            material_type());
    }
  }

  Core::Communication::PotentiallyUnusedBufferScope material_scope{buffer};
  if (params_ != nullptr)  // materials are not accessible in postprocessing mode
  {
    // unpack elastic material
    elastic_material_ = create_elastic_material(params_->matid_elast_);
    elastic_material_->unpack(buffer);

    // unpack inelastic deformation gradient factors
    inelastic_->assign_to_source(params_, *elastic_material_);
    inelastic_->unpack_inelastic(buffer);
  }
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
void Mat::MultiplicativeSplitDefgradElastHyper::evaluate(
    const Core::LinAlg::Tensor<double, 3, 3>* defgrad,
    const Core::LinAlg::SymmetricTensor<double, 3, 3>& glstrain,
    const Teuchos::ParameterList& params, const EvaluationContext<3>& context,
    Core::LinAlg::SymmetricTensor<double, 3, 3>& stress,
    Core::LinAlg::SymmetricTensor<double, 3, 3, 3, 3>& cmat, int gp, int eleGID)
{
  // do all stuff that only has to be done once per evaluate() call
  pre_evaluate(params, context, gp, eleGID);

  const Core::LinAlg::Matrix<3, 3> defgrd_mat = Core::LinAlg::make_matrix_view(*defgrad);
  Core::LinAlg::Tensor<double, 3, 3> iFin{};
  Core::LinAlg::Matrix<3, 3> iFinM = Core::LinAlg::make_matrix_view(iFin);
  inelastic_->evaluate_inverse_inelastic_def_grad(&defgrd_mat, iFinM);

  const StressAndDerivatives stress_and_derivatives =
      evaluate_stress_and_derivatives(*defgrad, iFin, params, context, gp, eleGID);
  stress = stress_and_derivatives.stress;
  cmat = stress_and_derivatives.cmat_at_fixed_iFin;

  // evaluate additional terms for the elasticity tensor
  // cmatadd = 2 \frac{\partial S}{\partial F_{in}^{-1}} : \frac{\partial F_{in}^{-1}}{\partial
  // C}, where F_{in}^{-1} can be multiplicatively composed of several inelastic contributions
  const Core::LinAlg::SymmetricTensor<double, 3, 3> iC =
      Core::LinAlg::inv(Core::LinAlg::assume_symmetry(
          Core::LinAlg::dot(Core::LinAlg::transpose(*defgrad), *defgrad)));
  const Core::LinAlg::Matrix<6, 1> iCV = Core::LinAlg::make_stress_like_voigt_view(iC);
  Core::LinAlg::Matrix<6, 6> cmat_view = Core::LinAlg::make_stress_like_voigt_view(cmat);
  cmat_view.update(
      1.0, evaluate_additional_cmat(&defgrd_mat, iCV, stress_and_derivatives.dstress_diFin), 1.0);
}

Core::LinAlg::SymmetricTensor<double, 3, 3>
Mat::MultiplicativeSplitDefgradElastHyper::evaluate_d_stress_d_scalar(
    const Core::LinAlg::Tensor<double, 3, 3>& defgrad,
    const Core::LinAlg::SymmetricTensor<double, 3, 3>& glstrain,
    const Teuchos::ParameterList& params, const EvaluationContext<3>& context, int gp, int eleGID)
{
  Core::LinAlg::Matrix<3, 3> defgrad_mat = Core::LinAlg::make_matrix_view(defgrad);
  // do all stuff that only has to be done once per evaluate() call
  pre_evaluate(params, context, gp, eleGID);

  // get source of deformation for this OD block depending on the differentiation type
  auto source(PAR::InelasticSource::none);
  const int differentiationtype = get_or<int>(
      params, "differentiationtype", static_cast<int>(Solid::DifferentiationType::none));
  if (differentiationtype == static_cast<int>(Solid::DifferentiationType::elch))
    source = PAR::InelasticSource::concentration;
  else if (differentiationtype == static_cast<int>(Solid::DifferentiationType::temp))
    source = PAR::InelasticSource::temperature;
  else
    FOUR_C_THROW("unknown scalaratype");

  Core::LinAlg::Tensor<double, 3, 3> iFin{};
  Core::LinAlg::Matrix<3, 3> iFinM = Core::LinAlg::make_matrix_view(iFin);
  inelastic_->evaluate_inverse_inelastic_def_grad(&defgrad_mat, iFinM);

  const StressAndDerivatives stress_and_derivatives =
      evaluate_stress_and_derivatives(defgrad, iFin, params, context, gp, eleGID);

  Core::LinAlg::SymmetricTensor<double, 3, 3> d_stress_d_scalar =
      make_symmetric_tensor_from_stress_like_voigt_matrix(
          evaluate_od_stiff_mat(source, &defgrad_mat, stress_and_derivatives.dstress_diFin));

  if (source == PAR::InelasticSource::temperature)
  {
    // stress depends on temperature not only through the inelastic deformation gradient but also
    // through thermal expansion
    d_stress_d_scalar += stress_and_derivatives.dstress_dT;
  }

  return d_stress_d_scalar;
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
double Mat::MultiplicativeSplitDefgradElastHyper::evaluate_cauchy_n_dir_and_derivatives(
    const Core::LinAlg::Tensor<double, 3, 3>& defgrd, const Core::LinAlg::Tensor<double, 3>& n,
    const Core::LinAlg::Tensor<double, 3>& dir, Core::LinAlg::Matrix<3, 1>* d_cauchyndir_dn,
    Core::LinAlg::Matrix<3, 1>* d_cauchyndir_ddir, Core::LinAlg::Matrix<9, 1>* d_cauchyndir_dF,
    Core::LinAlg::Matrix<9, 9>* d2_cauchyndir_dF2, Core::LinAlg::Matrix<9, 3>* d2_cauchyndir_dF_dn,
    Core::LinAlg::Matrix<9, 3>* d2_cauchyndir_dF_ddir, const EvaluationContext<3>& context,
    int eleGID, const double* concentration, const double* temp, double* d_cauchyndir_dT,
    Core::LinAlg::Matrix<9, 1>* d2_cauchyndir_dF_dT)
{
  // the thermal expansion stress S_T is not considered here/is deprecated
  FOUR_C_ASSERT_ALWAYS(
      temp == nullptr && d_cauchyndir_dT == nullptr && d2_cauchyndir_dF_dT == nullptr,
      "Temperature dependence of the Cauchy stress is not implemented for "
      "MAT_MultiplicativeSplitDefgradElastHyper.");

  if (concentration != nullptr) set_concentration_gp(*concentration);

  const Core::LinAlg::Matrix<3, 3> defgrd_mat = Core::LinAlg::make_matrix_view(defgrd);
  Core::LinAlg::Tensor<double, 3, 3> iFin{};
  Core::LinAlg::Matrix<3, 3> iFinM = Core::LinAlg::make_matrix_view(iFin);
  inelastic_->evaluate_inverse_inelastic_def_grad(&defgrd_mat, iFinM);
  const Core::LinAlg::Tensor<double, 3, 3> Fe = Core::LinAlg::dot(defgrd, iFin);

  // the Cauchy stress equals the Cauchy stress of the elastic material evaluated at F_e, the
  // derivatives w.r.t. F follow from dF_e/dF at constant F_in.
  // Note: this neglects the dependence of F_in on F (e.g., via the material concentration c det F
  // of scalar-dependent growth or via the local Newton of viscoplasticity). Including it would
  // require dF_in^{-1}/dC from the inelastic factors and, for the second derivatives,
  // d^2F_in^{-1}/dC^2, which the factors do not provide.
  Core::LinAlg::Matrix<9, 1> d_cauchyndir_dFe(Core::LinAlg::Initialization::zero);
  Core::LinAlg::Matrix<9, 9> d2_cauchyndir_dFe2(Core::LinAlg::Initialization::zero);
  Core::LinAlg::Matrix<9, 3> d2_cauchyndir_dFe_dn(Core::LinAlg::Initialization::zero);
  Core::LinAlg::Matrix<9, 3> d2_cauchyndir_dFe_ddir(Core::LinAlg::Initialization::zero);
  const double cauchy_n_dir = elastic_material_->evaluate_cauchy_n_dir_and_derivatives(Fe, n, dir,
      d_cauchyndir_dn, d_cauchyndir_ddir, d_cauchyndir_dF ? &d_cauchyndir_dFe : nullptr,
      d2_cauchyndir_dF2 ? &d2_cauchyndir_dFe2 : nullptr,
      d2_cauchyndir_dF_dn ? &d2_cauchyndir_dFe_dn : nullptr,
      d2_cauchyndir_dF_ddir ? &d2_cauchyndir_dFe_ddir : nullptr, context, eleGID, nullptr, nullptr,
      nullptr, nullptr);

  Core::LinAlg::Matrix<9, 9> d_Fe_dF(Core::LinAlg::Initialization::zero);
  Core::LinAlg::FourTensorOperations::add_non_symmetric_product(
      1.0, Core::LinAlg::identity_matrix<3>(), iFinM, d_Fe_dF);

  if (d_cauchyndir_dF) d_cauchyndir_dF->multiply_tn(1.0, d_Fe_dF, d_cauchyndir_dFe, 0.0);
  if (d2_cauchyndir_dF2)
  {
    Core::LinAlg::Matrix<9, 9> d2_cauchyndir_dFe2_dFe_dF(Core::LinAlg::Initialization::zero);
    d2_cauchyndir_dFe2_dFe_dF.multiply_nn(1.0, d2_cauchyndir_dFe2, d_Fe_dF, 0.0);
    d2_cauchyndir_dF2->multiply_tn(1.0, d_Fe_dF, d2_cauchyndir_dFe2_dFe_dF, 0.0);
  }
  if (d2_cauchyndir_dF_dn)
    d2_cauchyndir_dF_dn->multiply_tn(1.0, d_Fe_dF, d2_cauchyndir_dFe_dn, 0.0);
  if (d2_cauchyndir_dF_ddir)
    d2_cauchyndir_dF_ddir->multiply_tn(1.0, d_Fe_dF, d2_cauchyndir_dFe_ddir, 0.0);

  return cauchy_n_dir;
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
Mat::MultiplicativeSplitDefgradElastHyper::StressAndDerivatives
Mat::MultiplicativeSplitDefgradElastHyper::evaluate_stress_and_derivatives(
    const Core::LinAlg::Tensor<double, 3, 3>& defgrad,
    const Core::LinAlg::Tensor<double, 3, 3>& iFin, const Teuchos::ParameterList& params,
    const EvaluationContext<3>& context, const int gp, const int eleGID)
{
  const auto& id = Core::LinAlg::TensorGenerators::identity<double, 3, 3>;

  // kinematics: F_e = F iFin and det(F_in)
  const Core::LinAlg::Tensor<double, 3, 3> Fe = Core::LinAlg::dot(defgrad, iFin);
  const double detFin = 1.0 / Core::LinAlg::det(iFin);
  const Core::LinAlg::Tensor<double, 3, 3> iFinT = Core::LinAlg::transpose(iFin);

  // elastic stress and elasticity tensor at E_e = 1/2 (F_e^T F_e - I); F_e is only passed since the
  // So3Material interface requires it
  const Core::LinAlg::SymmetricTensor<double, 3, 3> Ee =
      0.5 *
      (Core::LinAlg::assume_symmetry(Core::LinAlg::dot(Core::LinAlg::transpose(Fe), Fe)) - id);
  Core::LinAlg::SymmetricTensor<double, 3, 3> Se{};
  Core::LinAlg::SymmetricTensor<double, 3, 3, 3, 3> cmat_e{};
  elastic_material_->evaluate(&Fe, Ee, params, context, Se, cmat_e, gp, eleGID);

  // thermal stress S_T (zero without temperature)
  TensorAndTemperatureDerivative ST{};
  if (params.isParameter("temperature"))
  {
    ST = evaluate_thermal_expansion_stress(*elastic_material_,
        params_->thermal_expansion_coefficient_,
        params.get<double>("temperature") - params_->ref_temperature_, params, context, gp, eleGID);
  }
  const Core::LinAlg::SymmetricTensor<double, 3, 3> Se_minus_ST = Se - ST.value;

  StressAndDerivatives stress_and_derivatives;

  stress_and_derivatives.stress = detFin * Core::LinAlg::assume_symmetry(Core::LinAlg::dot(
                                               Core::LinAlg::dot(iFin, Se_minus_ST), iFinT));
  stress_and_derivatives.dstress_dT =
      -detFin * Core::LinAlg::assume_symmetry(
                    Core::LinAlg::dot(Core::LinAlg::dot(iFin, ST.temperature_derivative), iFinT));

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
  const Core::LinAlg::Tensor<double, 3, 3> iFin_Se_minus_ST = Core::LinAlg::dot(iFin, Se_minus_ST);
  const Core::LinAlg::Tensor<double, 3, 3> C_iFin =
      Core::LinAlg::dot(Core::LinAlg::dot(Core::LinAlg::transpose(defgrad), defgrad), iFin);

  // assemble dS_diFin (note that cmat_e_index_12_pulled_back is symmetric in the last two indices,
  // which allows for simplification)
  const Core::LinAlg::Tensor<double, 3, 3, 3, 3> dS_diFin =
      detFin * (Core::LinAlg::einsum<"ik", "jl">(id, iFin_Se_minus_ST) +
                   Core::LinAlg::einsum<"jk", "il">(id, iFin_Se_minus_ST) +
                   Core::LinAlg::einsum<"ijlb", "kb">(cmat_e_index_12_pulled_back, C_iFin)) -
      Core::LinAlg::einsum<"ij", "kl">(stress_and_derivatives.stress, Core::LinAlg::inv(iFinT));
  stress_and_derivatives.dstress_diFin = Core::LinAlg::make_6x9_voigt_matrix_from_tensor(dS_diFin);

  return stress_and_derivatives;
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
void Mat::MultiplicativeSplitDefgradElastHyper::evaluate_linearization_od(
    const Core::LinAlg::Tensor<double, 3, 3>& defgrd, const double concentration,
    Core::LinAlg::Matrix<9, 1>& d_F_dx)
{
  set_concentration_gp(concentration);
  const Core::LinAlg::Matrix<3, 3> defgrd_mat = Core::LinAlg::make_matrix_view(defgrd);

  // References to vector of inelastic contributions and inelastic deformation gradients
  auto facdefgradin = inelastic_->fac_def_grad_in();

  // number of contributions for this source
  const int num_contributions = inelastic_->num_inelastic_def_grad();

  // build inverse inelastic deformation gradient
  Core::LinAlg::Tensor<double, 3, 3> iFin{};
  Core::LinAlg::Matrix<3, 3> iFinM = Core::LinAlg::make_matrix_view(iFin);

  inelastic_->evaluate_inverse_inelastic_def_grad(&defgrd_mat, iFinM);

  Core::LinAlg::Tensor<double, 3, 3> Fe = defgrd * iFin;

  // calculate the derivative of the deformation gradient w.r.t. the inelastic deformation gradient
  Core::LinAlg::Tensor<double, 3, 3, 3, 3> d_F_dFin =
      Core::LinAlg::reorder_axis<0, 2, 3, 1>(Core::LinAlg::dyadic(
          Fe, Core::LinAlg::get_full(Core::LinAlg::TensorGenerators::identity<double, 3, 3>)));

  // remains zero for inelastic factors that do not depend on the scalar
  Core::LinAlg::Tensor<double, 3, 3> d_Fin_dx{};

  // check number of factors the inelastic deformation gradient consists of and choose
  // implementation accordingly
  if (num_contributions == 1)
  {
    facdefgradin[0].second->evaluate_inelastic_def_grad_derivative(
        Core::LinAlg::det(defgrd), d_Fin_dx);
  }
  else
    FOUR_C_THROW("NOT YET IMPLEMENTED");

  auto dFdc = Core::LinAlg::ddot(d_F_dFin, d_Fin_dx);

  Core::LinAlg::Voigt::matrix_3x3_to_9x1(Core::LinAlg::make_matrix_view(dFdc), d_F_dx);
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
Core::LinAlg::Matrix<6, 6> Mat::MultiplicativeSplitDefgradElastHyper::evaluate_additional_cmat(
    const Core::LinAlg::Matrix<3, 3>* const defgrad, const Core::LinAlg::Matrix<6, 1>& iCV,
    const Core::LinAlg::Matrix<6, 9>& dSdiFin)
{
  // declare output variable
  Core::LinAlg::Matrix<6, 6> cmatadd;

  const auto& facdefgradin = inelastic_->fac_def_grad_in();
  const auto& iFinjM = inelastic_->geti_finj();
  const int num_contributions = inelastic_->num_inelastic_def_grad();

  // check amount of factors the inelastic deformation gradient consists of and choose
  // implementation accordingly
  if (num_contributions == 1)
  {
    Core::LinAlg::Matrix<3, 3> id3x3(Core::LinAlg::Initialization::zero);
    for (int i = 0; i < 3; ++i) id3x3(i, i) = 1.0;
    facdefgradin[0].second->evaluate_additional_cmat(
        defgrad, id3x3, iFinjM[0].second, iCV, dSdiFin, cmatadd);
  }
  else if (num_contributions > 1)
  {
    // dSdiFinj = dSdiFin : diFindiFinj
    // diFindiFinj = \Pi_(k=num_contributions_part)^(j+1) iFin_k : dFinjdFinj : \Pi_(k=j-1)^(0). The
    // double contraction and derivative in index notation for three contributions: diFindiFinj_abcd
    // = iFin_(j-1)_ae \delta_ec \delta_df iFin_(j+1)_fb = iFin_(j-1)_ac iFin_(j+1)_db. This is
    // performed by nonsymmetric product \Pi_(k=num_contributions_part)^(j+1) iFin_k (x)
    // \Pi_(k=j-1)^(0) iFin_k, where (x) denots the nonsymmetric product and \Pi is the
    // multiplication operator.
    Core::LinAlg::Matrix<6, 9> dSdiFinj(Core::LinAlg::Initialization::zero);
    Core::LinAlg::Matrix<9, 9> diFindiFinj(Core::LinAlg::Initialization::zero);
    Core::LinAlg::Matrix<3, 3> id(Core::LinAlg::Initialization::zero);
    for (int i = 0; i < 3; ++i) id(i, i) = 1.0;

    // product of all iFinj, except for the one, that is currently evaluated. In case of inner
    // (neither first or last) we have two products
    Core::LinAlg::Matrix<3, 3> producta(Core::LinAlg::Initialization::zero);
    Core::LinAlg::Matrix<3, 3> productb(Core::LinAlg::Initialization::zero);
    Core::LinAlg::Matrix<3, 3> producta_temp(Core::LinAlg::Initialization::zero);
    Core::LinAlg::Matrix<3, 3> productb_temp(Core::LinAlg::Initialization::zero);

    for (int i = 0; i < num_contributions; ++i)
    {
      // reset the product of the previous contribution
      diFindiFinj.clear();

      // multiply all inelastic deformation gradients except for range between first and current
      producta = id;
      for (int j = num_contributions - 1; j > i; --j)
      {
        producta_temp.multiply(1.0, producta, iFinjM[j].second, 0.0);
        producta.update(1.0, producta_temp, 0.0);
      }

      // multiply all inelastic deformation gradients except for range between last and current
      productb = id;
      if (i > 0)
      {
        for (int j = i - 1; j >= 0; --j)
        {
          productb_temp.multiply(1.0, productb, iFinjM[j].second, 0.0);
          productb.update(1.0, productb_temp, 0.0);
        }
      }

      // evaluate additional contribution to C by applying chain rule
      Core::LinAlg::FourTensorOperations::add_non_symmetric_product(
          1.0, producta, productb, diFindiFinj);
      dSdiFinj.multiply(1.0, dSdiFin, diFindiFinj, 0.0);
      facdefgradin[i].second->evaluate_additional_cmat(
          defgrad, productb, iFinjM[i].second, iCV, dSdiFinj, cmatadd);
    }
  }
  else
    FOUR_C_THROW("You should not be here");

  return cmatadd;
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
void Mat::MultiplicativeSplitDefgradElastHyper::setup(const int numgp,
    const Discret::Elements::Fibers& fibers,
    const std::optional<Discret::Elements::CoordinateSystem>& coord_system)
{
  elastic_material_->setup(numgp, fibers, coord_system);

  // setup inelastic materials
  inelastic_->setup(numgp, fibers, coord_system);
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
void Mat::MultiplicativeSplitDefgradElastHyper::post_setup(
    const Teuchos::ParameterList& params, const int eleGID)
{
  elastic_material_->post_setup(params, eleGID);
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
void Mat::MultiplicativeSplitDefgradElastHyper::update()
{
  elastic_material_->update();

  // update inelastic materials
  inelastic_->update();
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
Core::LinAlg::Matrix<6, 1> Mat::MultiplicativeSplitDefgradElastHyper::evaluate_od_stiff_mat(
    const PAR::InelasticSource source, const Core::LinAlg::Matrix<3, 3>* const defgrad,
    const Core::LinAlg::Matrix<6, 9>& dSdiFin)
{
  Core::LinAlg::Matrix<6, 1> dstressdx{Core::LinAlg::Initialization::zero};

  // References to vector of inelastic contributions and inelastic deformation gradients
  const auto& facdefgradin = inelastic_->fac_def_grad_in();
  const auto& iFinjM = inelastic_->geti_finj();

  // number of contributions for this source
  const int num_contributions = inelastic_->num_inelastic_def_grad();

  // check number of factors the inelastic deformation gradient consists of and choose
  // implementation accordingly
  if (num_contributions == 1)
  {
    if (facdefgradin[0].first == source)
    {
      facdefgradin[0].second->evaluate_od_stiff_mat(
          defgrad, Core::LinAlg::identity_matrix<3>(), iFinjM[0].second, dSdiFin, dstressdx);
    }
  }
  else if (num_contributions > 1)
  {
    // dSdiFinj = dSdiFin : diFindiFinj
    // diFindiFinj = \Pi_(k=num_contributions_part)^(j+1) iFin_k : dFinjdFinj : \Pi_(k=j-1)^(0). The
    // double contraction and derivative in index notation for three contributions: diFindiFinj_abcd
    // = iFin_(j-1)_ae \delta_ec \delta_df iFin_(j+1)_fb = iFin_(j-1)_ac iFin_(j+1)_db. This is
    // performed by nonsymmetric product \Pi_(k=num_contributions_part)^(j+1) iFin_k (x)
    // \Pi_(k=j-1)^(0) iFin_k, where (x) denots the nonsymmetric product and \Pi is the
    // multiplication operator.
    Core::LinAlg::Matrix<6, 9> dSdiFinj(Core::LinAlg::Initialization::zero);
    Core::LinAlg::Matrix<9, 9> diFindiFinj(Core::LinAlg::Initialization::zero);
    Core::LinAlg::Matrix<3, 3> id(Core::LinAlg::Initialization::zero);
    for (int i = 0; i < 3; ++i) id(i, i) = 1.0;

    // product of all iFinj, except for the one, that is currently evaluated. In case of inner
    // (neither first or last) we have two products
    Core::LinAlg::Matrix<3, 3> producta(Core::LinAlg::Initialization::zero);
    Core::LinAlg::Matrix<3, 3> productb(Core::LinAlg::Initialization::zero);
    Core::LinAlg::Matrix<3, 3> producta_temp(Core::LinAlg::Initialization::zero);
    Core::LinAlg::Matrix<3, 3> productb_temp(Core::LinAlg::Initialization::zero);

    for (int i = 0; i < num_contributions; ++i)
    {
      // only if the contribution is from this source, the derivative is non-zero
      if (facdefgradin[i].first == source)
      {
        // reset the product of the previous contribution
        diFindiFinj.clear();

        // multiply all inelastic deformation gradients except for range between first and current
        // in reverse order
        producta = id;
        for (int j = num_contributions - 1; j > i; --j)
        {
          producta_temp.multiply(1.0, producta, iFinjM[j].second, 0.0);
          producta.update(1.0, producta_temp, 0.0);
        }

        // multiply all inelastic deformation gradients except for range between last and current
        productb = id;
        if (i > 0)
        {
          for (int j = i - 1; j >= 0; --j)
          {
            productb_temp.multiply(1.0, productb, iFinjM[j].second, 0.0);
            productb.update(1.0, productb_temp, 0.0);
          }
        }

        // evaluate additional contribution to OD block by applying chain rule
        Core::LinAlg::FourTensorOperations::add_non_symmetric_product(
            1.0, producta, productb, diFindiFinj);
        dSdiFinj.multiply(1.0, dSdiFin, diFindiFinj, 0.0);
        facdefgradin[i].second->evaluate_od_stiff_mat(
            defgrad, productb, iFinjM[i].second, dSdiFinj, dstressdx);
      }
    }
  }
  else
    FOUR_C_THROW("You should not be here");

  return dstressdx;
}

void Mat::MultiplicativeSplitDefgradElastHyper::stress_temperature_modulus_and_deriv(
    Core::LinAlg::SymmetricTensor<double, 3, 3>& stm,
    Core::LinAlg::SymmetricTensor<double, 3, 3>& stm_dT, const int gp)
{
  // this contribution is intentionally neglected so far.
  stm.fill(0.0);
  stm_dT.fill(0.0);
}

Mat::HeatSource Mat::MultiplicativeSplitDefgradElastHyper::evaluate_additional_heat_source(
    const EvaluationContext<3>& context, const int gp, const int eleGID,
    const Core::LinAlg::Matrix<3, 3>* current_defgrad, const double current_temperature)
{
  HeatSource heat_source;

  // References to vector of inelastic contributions and inelastic deformation gradients
  const auto& facdefgradin = inelastic_->fac_def_grad_in();

  // number of contributions
  const int num_contributions = inelastic_->num_inelastic_def_grad();

  // check number of factors the inelastic deformation gradient consists of and choose
  // implementation accordingly
  if (num_contributions == 1)
  {
    if (const auto viscoplastic_factor =
            std::dynamic_pointer_cast<Mat::InelasticDefgradTransvIsotropElastViscoplast>(
                facdefgradin[0].second))
    {
      heat_source = viscoplastic_factor->evaluate_taylor_quinney_heat_source(context, gp, eleGID,
          current_defgrad, Core::LinAlg::identity_matrix<3>(), current_temperature);
    }
  }
  else if (num_contributions > 1)
  {
    FOUR_C_THROW(
        "Evaluation of additional heat source is only implemented for a single factor so far.");
  }
  else
    FOUR_C_THROW("You should not be here");

  return heat_source;
}
/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
void Mat::MultiplicativeSplitDefgradElastHyper::pre_evaluate(const Teuchos::ParameterList& params,
    const EvaluationContext<3>& context, const int gp, const int eleGID) const
{
  // loop over all inelastic contributions
  for (int p = 0; p < inelastic_->num_inelastic_def_grad(); ++p)
    inelastic_->fac_def_grad_in()[p].second->pre_evaluate(params, context, gp, eleGID);
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
void Mat::MultiplicativeSplitDefgradElastHyper::set_concentration_gp(const double concentration)
{
  for (int p = 0; p < inelastic_->num_inelastic_def_grad(); ++p)
    inelastic_->fac_def_grad_in()[p].second->set_concentration_gp(concentration);
}

void Mat::MultiplicativeSplitDefgradElastHyper::register_output_data_names(
    std::unordered_map<std::string, int>& names_and_size) const
{
  // register output for the inelastic defgrad factors
  inelastic_->register_output_data_names(names_and_size);
}

bool Mat::MultiplicativeSplitDefgradElastHyper::evaluate_output_data(
    const std::string& name, Core::LinAlg::SerialDenseMatrix& data) const
{
  // evaluate GP output for the inelastic defgrad factors
  return inelastic_->evaluate_output_data(name, data);
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
void Mat::InelasticFactorsHandler::assign_to_source(
    Mat::PAR::MultiplicativeSplitDefgradElastHyper* params, Mat::So3Material& elastic_material)
{
  facdefgradin_.clear();
  i_finj_.clear();

  // get inelastic deformation gradient factors and assign them to their source
  for (int inelastic_matnum : params->inel_defgradfacids_)
  {
    auto inelastic_factor =
        Mat::InelasticDefgradFactors::factory(inelastic_matnum, elastic_material);
    if (inelastic_factor == nullptr) FOUR_C_THROW("Failed to allocate!");
    std::pair<PAR::InelasticSource, std::shared_ptr<Mat::InelasticDefgradFactors>> temppair(
        inelastic_factor->get_inelastic_source(), inelastic_factor);
    facdefgradin_.push_back(temppair);
  }

  i_finj_.resize(facdefgradin_.size());

  // safety checks
  // get the scatra structure control parameter list
  const auto& ssicontrol = Global::Problem::instance()->ssi_control_params();
  if (Teuchos::getIntegralValue<SSI::SolutionSchemeOverFields>(ssicontrol, "COUPALGO") ==
      SSI::SolutionSchemeOverFields::ssi_Monolithic)
  {
    for (const auto& inelasitc_factor : facdefgradin_)
    {
      const auto materialtype = inelasitc_factor.second->material_type();
      if ((materialtype != Core::Materials::mfi_lin_scalar_aniso) and
          (materialtype != Core::Materials::mfi_lin_scalar_iso) and
          (materialtype != Core::Materials::mfi_lin_temp_iso) and
          (materialtype != Core::Materials::mfi_no_growth) and
          (materialtype != Core::Materials::mfi_time_funct_aniso) and
          (materialtype != Core::Materials::mfi_time_funct_iso) and
          (materialtype != Core::Materials::mfi_poly_intercal_frac_aniso) and
          (materialtype != Core::Materials::mfi_poly_intercal_frac_iso) and
          (materialtype != Core::Materials::mfi_transv_isotrop_elast_viscoplast))
      {
        FOUR_C_THROW(
            "When you use the 'COUPALGO' 'ssi_Monolithic' from the 'SSI CONTROL' section, you need "
            "to use one of the materials derived from 'Mat::InelasticDefgradFactors'!"
            " If you want to use a different material, feel free to implement it! ;-)");
      }
    }
  }
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
void Mat::InelasticFactorsHandler::evaluate_inverse_inelastic_def_grad(
    const Core::LinAlg::Matrix<3, 3>* const defgrad, Core::LinAlg::Matrix<3, 3>& iFinM)
{
  // temporary variables
  Core::LinAlg::Matrix<3, 3> iFinp(Core::LinAlg::Initialization::zero);
  Core::LinAlg::Matrix<3, 3> iFin_init_store(Core::LinAlg::Initialization::zero);

  // clear variables
  iFinM.clear();
  iFin_init_store.clear();

  for (int i = 0; i < 3; ++i) iFin_init_store(i, i) = 1.0;

  for (int i = 0; i < num_inelastic_def_grad(); ++i)
  {
    // clear tmp variable
    iFinp.clear();

    // calculate inelastic deformation gradient and its inverse
    facdefgradin_[i].second->evaluate_inverse_inelastic_def_grad(defgrad, iFin_init_store, iFinp);

    // store inelastic deformation gradient of p-th inelastic contribution
    i_finj_[i].second = iFinp;

    // update inverse inelastic deformation gradient
    iFinM.multiply(iFin_init_store, iFinp);

    // store result for next evaluation
    iFin_init_store.update(1.0, iFinM, 0.0);
  }
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
void Mat::InelasticFactorsHandler::update()
{
  // loop over all inelastic contributions
  for (const auto& [_, inelastic_defgrad_factor] : facdefgradin_)
  {
    inelastic_defgrad_factor->update();
  }
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
void Mat::InelasticFactorsHandler::setup(const int numgp, const Discret::Elements::Fibers& fibers,
    const std::optional<Discret::Elements::CoordinateSystem>& coord_system)
{
  // loop over all inelastic contributions
  for (int i = 0; i < num_inelastic_def_grad(); ++i)
  {
    facdefgradin_[i].second->setup(numgp, fibers, coord_system);
  }
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
void Mat::InelasticFactorsHandler::pack_inelastic(Core::Communication::PackBuffer& data) const
{
  // loop over all inelastic contributions
  for (int i = 0; i < num_inelastic_def_grad(); ++i)
  {
    facdefgradin_[i].second->pack_inelastic(data);
  }
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/

void Mat::InelasticFactorsHandler::unpack_inelastic(Core::Communication::UnpackBuffer& buffer)
{
  // loop over all inelastic contributions
  for (int i = 0; i < num_inelastic_def_grad(); ++i)
  {
    facdefgradin_[i].second->unpack_inelastic(buffer);
  }
}


void Mat::InelasticFactorsHandler::register_output_data_names(
    std::unordered_map<std::string, int>& names_and_size) const
{
  // loop over all inelastic contributions
  for (int i = 0; i < num_inelastic_def_grad(); ++i)
  {
    facdefgradin_[i].second->register_output_data_names(names_and_size);
  }
}

bool Mat::InelasticFactorsHandler::evaluate_output_data(
    const std::string& name, Core::LinAlg::SerialDenseMatrix& data) const
{
  // declare the number of inelastic defgrad factors for which we use Gauss
  // point data
  int num_inel_fac_GP_data = 0;

  // loop over all inelastic contributions
  for (int i = 0; i < num_inelastic_def_grad(); ++i)
  {
    if (facdefgradin_[i].second->evaluate_output_data(name, data))
    {
      ++num_inel_fac_GP_data;
    };
  }

  return (num_inel_fac_GP_data > 0);
}


FOUR_C_NAMESPACE_CLOSE
