// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <gtest/gtest.h>

#include "4C_global_data.hpp"
#include "4C_linalg_fixedsizematrix.hpp"
#include "4C_linalg_fixedsizematrix_voigt_notation.hpp"
#include "4C_linalg_symmetric_tensor.hpp"
#include "4C_linalg_tensor.hpp"
#include "4C_mat_material_factory.hpp"
#include "4C_mat_multiplicative_split_defgrad_elasthyper.hpp"
#include "4C_mat_multiplicative_split_defgrad_elasthyper_service.hpp"
#include "4C_mat_par_bundle.hpp"
#include "4C_material_parameter_base.hpp"
#include "4C_solid_ele_fibers.hpp"
#include "4C_ssi_input.hpp"
#include "4C_unittest_utils_assertions_test.hpp"
#include "4C_utils_singleton_owner.hpp"

namespace
{
  using namespace FourC;

  class MultiplicativeSplitDefgradElastHyperTest : public ::testing::Test
  {
   protected:
    void SetUp() override
    {
      F_ = {{{1.1, 0.01, 0.03}, {0.04, 1.2, 0.02}, {0.06, 0.05, 1.3}}};

      // set up the deformation gradient
      FM_(0, 0) = 1.1;
      FM_(0, 1) = 0.01;
      FM_(0, 2) = 0.03;
      FM_(1, 0) = 0.04;
      FM_(1, 1) = 1.2;
      FM_(1, 2) = 0.02;
      FM_(2, 0) = 0.06;
      FM_(2, 1) = 0.05;
      FM_(2, 2) = 1.3;

      set_ref_values_inverse_right_cauchy_green();

      set_ref_values_evaluated_sdi_fin();

      set_up_multiplicative_split_object_with_requirements();
    }


    // note: requirement of the generation of an object of the multiplicative split material is that
    // the global problem is equipped with the corresponding elastic and inelastic material. This is
    // set up within this method.
    void set_up_multiplicative_split_object_with_requirements()
    {
      // do problem instance specific stuff
      const int problemid(0);
      Global::Problem& problem = (*Global::Problem::instance());
      problem.materials()->set_read_from_problem(problemid);

      // set up elastic material to be added to problem instance
      const int matid_elastic(1);
      Core::IO::InputParameterContainer mat_elastic_neo_hooke_data;
      mat_elastic_neo_hooke_data.add("YOUNG", 1.5e2);
      mat_elastic_neo_hooke_data.add("NUE", 0.3);

      // add elastic material to problem instance
      problem.materials()->insert(
          matid_elastic, Mat::make_parameter(1, Core::Materials::MaterialType::mes_coupneohooke,
                             mat_elastic_neo_hooke_data));

      // set up hyperelastic material describing the elastic response
      const int matid_elast_hyper(4);
      Core::IO::InputParameterContainer mat_elast_hyper_data;
      mat_elast_hyper_data.add("NUMMAT", 1);
      mat_elast_hyper_data.add("MATIDS", std::vector<int>{matid_elastic});
      mat_elast_hyper_data.add("DENS", 1.32e1);
      mat_elast_hyper_data.add("POLYCONVEX", 0);
      problem.materials()->insert(matid_elast_hyper,
          Mat::make_parameter(matid_elast_hyper, Core::Materials::MaterialType::m_elasthyper,
              mat_elast_hyper_data));

      // set up inelastic material to be added to problem instance
      const int inelastic_defgrad_id(2);
      Core::IO::InputParameterContainer mat_inelastic_data;
      mat_inelastic_data.add("SCALAR1", 1);
      mat_inelastic_data.add("SCALAR1_MolarGrowthFac", 1.1);
      mat_inelastic_data.add("SCALAR1_RefConc", 1.2);

      // add inelastic material to problem instance
      problem.materials()->insert(inelastic_defgrad_id,
          Mat::make_parameter(
              1, Core::Materials::MaterialType::mfi_lin_scalar_iso, mat_inelastic_data));

      // set parameter list
      auto parameter_list_pointer = std::make_shared<Teuchos::ParameterList>();
      parameter_list_pointer->sublist("STRUCTURAL DYNAMIC", false)
          .set("MASSLIN", Solid::MassLin::ml_none);
      parameter_list_pointer->sublist("SSI CONTROL")
          .set("COUPALGO", SSI::SolutionSchemeOverFields::ssi_IterStagg);

      // set the parameter list in the global problem
      problem.set_parameter_list(parameter_list_pointer);

      // create MultiplicativeSplitDefgrad_ElastHyper object;
      // initialize container for material parameters first
      Core::IO::InputParameterContainer multiplicativeSplitDefgradData;

      multiplicativeSplitDefgradData.add("MATID_ELAST", matid_elast_hyper);
      multiplicativeSplitDefgradData.add("NUMFACINEL", 1);
      std::vector<int> inelastic_defgrad_factor_ids = {inelastic_defgrad_id};
      multiplicativeSplitDefgradData.add("INELDEFGRADFACIDS", inelastic_defgrad_factor_ids);
      multiplicativeSplitDefgradData.add("REF_TEMPERATURE", 293.0);
      multiplicativeSplitDefgradData.add("THERMAL_EXPANSION_COEFFICIENT", 0.1);

      // get pointer to parameter class
      parameters_multiplicative_split_defgrad_ =
          std::make_shared<Mat::PAR::MultiplicativeSplitDefgradElastHyper>(
              Core::Mat::PAR::Parameter::Data{.parameters = multiplicativeSplitDefgradData});

      // setup pointer to MultiplicativeSplitDefgrad_ElastHyper object
      multiplicative_split_defgrad_ = std::make_shared<Mat::MultiplicativeSplitDefgradElastHyper>(
          parameters_multiplicative_split_defgrad_.get());
      multiplicative_split_defgrad_->setup(1, Discret::Elements::Fibers{}, std::nullopt);
    }

    void set_ref_values_inverse_right_cauchy_green()
    {
      iCV_ref_(0) = 0.82944796303475;
      iCV_ref_(1) = 0.696728380339461;
      iCV_ref_(2) = 0.596601975917312;
      iCV_ref_(3) = -0.032359712004713;
      iCV_ref_(4) = -0.034647920866899;
      iCV_ref_(5) = -0.053120416239144;
    }

    void set_concentration_to_inelastic_material(const double concentration) const
    {
      Teuchos::ParameterList params;

      // set up a concentration vector and store it to the parameter list
      auto gpconc_lin =
          std::shared_ptr<std::vector<double>>(new std::vector<double>({concentration}));
      params.set<std::shared_ptr<std::vector<double>>>("scalars", gpconc_lin);



      // call pre evaluate to have concentration during actual call
      multiplicative_split_defgrad_->pre_evaluate(params, {}, 0, 0);
    }

    void set_ref_values_evaluated_sdi_fin()
    {
      dSdiFin_ref_(0, 0) = 114.1660435509188;
      dSdiFin_ref_(0, 1) = -1.661105324599284;
      dSdiFin_ref_(0, 2) = -1.657789966811478;
      dSdiFin_ref_(0, 3) = 0.6950689065562556;
      dSdiFin_ref_(0, 4) = -0.00329516271498935;
      dSdiFin_ref_(0, 5) = 0.9329432831258889;
      dSdiFin_ref_(0, 6) = 0.01011436563696388;
      dSdiFin_ref_(0, 7) = -0.006583593739068971;
      dSdiFin_ref_(0, 8) = 0.013463381879464;
      dSdiFin_ref_(1, 0) = -11.73283847873412;
      dSdiFin_ref_(1, 1) = 105.2701678289236;
      dSdiFin_ref_(1, 2) = -11.6040946680382;
      dSdiFin_ref_(1, 3) = -0.05898254736522643;
      dSdiFin_ref_(1, 4) = -0.4920621412911008;
      dSdiFin_ref_(1, 5) = -0.03535183236330931;
      dSdiFin_ref_(1, 6) = -0.5154481790422178;
      dSdiFin_ref_(1, 7) = -0.04608342825898787;
      dSdiFin_ref_(1, 8) = 0.0942401395889132;
      dSdiFin_ref_(2, 0) = -18.75794747350803;
      dSdiFin_ref_(2, 1) = -18.58921926562027;
      dSdiFin_ref_(2, 2) = 98.5798500483063;
      dSdiFin_ref_(2, 3) = -0.09429870933073022;
      dSdiFin_ref_(2, 4) = -0.03687574852596564;
      dSdiFin_ref_(2, 5) = -0.05651895879799482;
      dSdiFin_ref_(2, 6) = 0.1131885845970979;
      dSdiFin_ref_(2, 7) = -0.3081745971147502;
      dSdiFin_ref_(2, 8) = -0.201080649217473;
      dSdiFin_ref_(3, 0) = -2.547471706027534;
      dSdiFin_ref_(3, 1) = -1.882323090966196;
      dSdiFin_ref_(3, 2) = -2.22961181644894;
      dSdiFin_ref_(3, 3) = 58.4374016535939;
      dSdiFin_ref_(3, 4) = 0.4645651032280521;
      dSdiFin_ref_(3, 5) = -0.241290937884049;
      dSdiFin_ref_(3, 6) = 57.93471620610067;
      dSdiFin_ref_(3, 7) = -0.008854474142796649;
      dSdiFin_ref_(3, 8) = 0.01810730908547461;
      dSdiFin_ref_(4, 0) = -1.995076152570434;
      dSdiFin_ref_(4, 1) = -2.094379609332357;
      dSdiFin_ref_(4, 2) = -2.207682721946282;
      dSdiFin_ref_(4, 3) = -0.01002951450149832;
      dSdiFin_ref_(4, 4) = 58.56206171546444;
      dSdiFin_ref_(4, 5) = -0.006011298785499145;
      dSdiFin_ref_(4, 6) = -0.1638352028130229;
      dSdiFin_ref_(4, 7) = 58.44089844409201;
      dSdiFin_ref_(4, 8) = -0.2770982534557008;
      dSdiFin_ref_(5, 0) = -4.075765888689157;
      dSdiFin_ref_(5, 1) = -3.864812436613845;
      dSdiFin_ref_(5, 2) = -3.388101888830967;
      dSdiFin_ref_(5, 3) = -0.1368544954849143;
      dSdiFin_ref_(5, 4) = -0.007666693768910916;
      dSdiFin_ref_(5, 5) = 58.5542331450641;
      dSdiFin_ref_(5, 6) = 0.02353259936218022;
      dSdiFin_ref_(5, 7) = 0.3364299251029703;
      dSdiFin_ref_(5, 8) = 57.95243768164728;
    }

    // defined input quantities
    Core::LinAlg::Tensor<double, 3, 3> F_;
    Core::LinAlg::Matrix<3, 3> FM_;
    // reference solutions
    Core::LinAlg::Matrix<6, 1> iCV_ref_;
    Core::LinAlg::Matrix<6, 9> dSdiFin_ref_;

    // pointer to material parameters
    std::shared_ptr<Mat::PAR::MultiplicativeSplitDefgradElastHyper>
        parameters_multiplicative_split_defgrad_;

    // pointer to material
    std::shared_ptr<Mat::MultiplicativeSplitDefgradElastHyper> multiplicative_split_defgrad_;

    Core::Utils::SingletonOwnerRegistry::ScopeGuard guard;
  };

  TEST_F(MultiplicativeSplitDefgradElastHyperTest, TestEvaluateAdditionalCmat)
  {
    const double concentration(44327.362);
    set_concentration_to_inelastic_material(concentration);


    // reference solution
    Core::LinAlg::Matrix<6, 6> cMatAdd_ref;
    cMatAdd_ref(0, 0) = -0.701208493301168;
    cMatAdd_ref(0, 1) = -0.5890084484992926;
    cMatAdd_ref(0, 2) = -0.5043624088277512;
    cMatAdd_ref(0, 3) = 0.02735663466513797;
    cMatAdd_ref(0, 4) = 0.029291067637571;
    cMatAdd_ref(0, 5) = 0.04490756345738416;
    cMatAdd_ref(1, 0) = -0.5183018322515506;
    cMatAdd_ref(1, 1) = -0.4353685971936842;
    cMatAdd_ref(1, 2) = -0.3728020454851408;
    cMatAdd_ref(1, 3) = 0.02022079596387239;
    cMatAdd_ref(1, 4) = 0.02165064195626716;
    cMatAdd_ref(1, 5) = 0.03319365444696364;
    cMatAdd_ref(2, 0) = -0.3873521175618817;
    cMatAdd_ref(2, 1) = -0.3253720854320112;
    cMatAdd_ref(2, 2) = -0.2786130643659098;
    cMatAdd_ref(2, 3) = 0.01511198233926229;
    cMatAdd_ref(2, 4) = 0.01618057565396367;
    cMatAdd_ref(2, 5) = 0.02480722918495967;
    cMatAdd_ref(3, 0) = 0.04212677142365208;
    cMatAdd_ref(3, 1) = 0.03538608632607147;
    cMatAdd_ref(3, 2) = 0.03030077375609256;
    cMatAdd_ref(3, 3) = -0.00164351502651245;
    cMatAdd_ref(3, 4) = -0.001759730697661004;
    cMatAdd_ref(3, 5) = -0.00269792890279474;
    cMatAdd_ref(4, 0) = 0.0398350977693272;
    cMatAdd_ref(4, 1) = 0.03346110230705895;
    cMatAdd_ref(4, 2) = 0.02865242799932495;
    cMatAdd_ref(4, 3) = -0.001554108695111724;
    cMatAdd_ref(4, 4) = -0.001664002296403264;
    cMatAdd_ref(4, 5) = -0.002551163024973279;
    cMatAdd_ref(5, 0) = 0.07166415111957784;
    cMatAdd_ref(5, 1) = 0.06019720363802254;
    cMatAdd_ref(5, 2) = 0.05154630075158284;
    cMatAdd_ref(5, 3) = -0.002795873152556784;
    cMatAdd_ref(5, 4) = -0.002993573976479339;
    cMatAdd_ref(5, 5) = -0.004589594171729136;

    // actual call that is tested
    Core::LinAlg::Matrix<6, 6> cMatAdd =
        multiplicative_split_defgrad_->evaluate_additional_cmat(&FM_, iCV_ref_, dSdiFin_ref_);

    FOUR_C_EXPECT_NEAR(cMatAdd, cMatAdd_ref, 1.0e-10);
  }

  TEST_F(MultiplicativeSplitDefgradElastHyperTest, TestEvaluateCauchyNDirAndDerivatives)
  {
    // we need to first set this dummy concentration. otherwise the call below results in a
    // segfault
    const double dummy_conc(-1.0);
    set_concentration_to_inelastic_material(dummy_conc);

    // input variables
    const Core::LinAlg::Tensor<double, 3> n = {
        {1.0 / std::numbers::sqrt2, 0, -1.0 / std::numbers::sqrt2}};
    const Core::LinAlg::Tensor<double, 3> dir = {
        {std::numbers::inv_sqrt3, -1.0 / std::numbers::sqrt3, -1.0 / std::numbers::sqrt3}};
    const double concentration(1.0);

    // output variables
    Core::LinAlg::Matrix<3, 1> d_cauchyndir_dn(Core::LinAlg::Initialization::zero);
    Core::LinAlg::Matrix<3, 1> d_cauchyndir_ddir(Core::LinAlg::Initialization::zero);
    Core::LinAlg::Matrix<9, 1> d_cauchyndir_dF(Core::LinAlg::Initialization::zero);


    double total_time = 0.0;
    double time_step_size = 1.0;
    Mat::EvaluationContext<3> context{.total_time = &total_time,
        .time_step_size = &time_step_size,
        .xi = {},
        .ref_coords = nullptr};
    const double cauchy_n_dir =
        multiplicative_split_defgrad_->evaluate_cauchy_n_dir_and_derivatives(F_, n, dir,
            &d_cauchyndir_dn, &d_cauchyndir_ddir, &d_cauchyndir_dF, nullptr, nullptr, nullptr,
            context, 0, &concentration, nullptr, nullptr, nullptr);

    const double cauchy_n_dir_ref(6.019860168755);
    Core::LinAlg::Matrix<3, 1> d_cauchyndir_dn_ref(Core::LinAlg::Initialization::zero);
    Core::LinAlg::Matrix<3, 1> d_cauchyndir_ddir_ref(Core::LinAlg::Initialization::zero);
    Core::LinAlg::Matrix<9, 1> d_cauchyndir_dF_ref(Core::LinAlg::Initialization::zero);
    d_cauchyndir_dn_ref(0) = -2.856437080521;
    d_cauchyndir_dn_ref(1) = -6.736850094992;
    d_cauchyndir_dn_ref(2) = -1.136980497476e+01;

    d_cauchyndir_ddir_ref(0) = -1.933304201727;
    d_cauchyndir_ddir_ref(1) = -8.793331859785e-01;
    d_cauchyndir_ddir_ref(2) = -1.148067468249e+01;

    d_cauchyndir_dF_ref(0) = 7.834365817988e+01;
    d_cauchyndir_dF_ref(1) = 4.257334930649e+01;
    d_cauchyndir_dF_ref(2) = 7.961614320095e+01;
    d_cauchyndir_dF_ref(3) = -2.192663678809e+01;
    d_cauchyndir_dF_ref(4) = 1.867810957230e+01;
    d_cauchyndir_dF_ref(5) = -4.292543557935e+01;
    d_cauchyndir_dF_ref(6) = -1.694113997401e+01;
    d_cauchyndir_dF_ref(7) = 1.982432031426e+01;
    d_cauchyndir_dF_ref(8) = -3.361951777558e+01;

    EXPECT_NEAR(cauchy_n_dir, cauchy_n_dir_ref, 1.0e-10);
    FOUR_C_EXPECT_NEAR(d_cauchyndir_dn, d_cauchyndir_dn_ref, 1.0e-10);
    FOUR_C_EXPECT_NEAR(d_cauchyndir_ddir, d_cauchyndir_ddir_ref, 1.0e-10);
    FOUR_C_EXPECT_NEAR(d_cauchyndir_dF, d_cauchyndir_dF_ref, 1.0e-10);
  }

  TEST_F(MultiplicativeSplitDefgradElastHyperTest, TestEvaluateLinearizationOD)
  {
    // we need to first set this dummy concentration. otherwise the call below results in a
    // segfault
    const double dummy_conc(-1.0);
    set_concentration_to_inelastic_material(dummy_conc);

    // actual material call
    const double concentration(1.0);
    Core::LinAlg::Matrix<9, 1> DFDx(Core::LinAlg::Initialization::zero);
    multiplicative_split_defgrad_->evaluate_linearization_od(F_, concentration, DFDx);

    // define the reference solution
    Core::LinAlg::Matrix<9, 1> DFDx_ref;
    DFDx_ref(0) = 4.417109534556e-01;
    DFDx_ref(1) = 4.818664946788e-01;
    DFDx_ref(2) = 5.220220359020e-01;
    DFDx_ref(3) = 4.015554122323e-03;
    DFDx_ref(4) = 8.031108244647e-03;
    DFDx_ref(5) = 1.204666236697e-02;
    DFDx_ref(6) = 1.606221648929e-02;
    DFDx_ref(7) = 2.007777061162e-02;
    DFDx_ref(8) = 2.409332473394e-02;

    FOUR_C_EXPECT_NEAR(DFDx, DFDx_ref, 1.0e-10);
  }

  TEST_F(MultiplicativeSplitDefgradElastHyperTest, TestEvaluateLinearizationODForSeveralStates)
  {
    // we need to first set this dummy concentration. otherwise the call below results in a
    // segfault
    const double dummy_conc(-1.0);
    set_concentration_to_inelastic_material(dummy_conc);

    const double concentration(1.0);
    Core::LinAlg::Matrix<9, 1> DFDx(Core::LinAlg::Initialization::zero);
    multiplicative_split_defgrad_->evaluate_linearization_od(F_, concentration, DFDx);
    Core::LinAlg::Matrix<9, 1> DFDx_transposed(Core::LinAlg::Initialization::zero);
    multiplicative_split_defgrad_->evaluate_linearization_od(
        Core::LinAlg::transpose(F_), concentration, DFDx_transposed);

    // the isotropic growth yields dF/dc = s(det F) F, i.e., the linearization evaluated with F^T
    // is the transpose of the linearization evaluated with F
    Core::LinAlg::Matrix<3, 3> DFDx_matrix;
    Core::LinAlg::Voigt::matrix_9x1_to_3x3(DFDx, DFDx_matrix);
    Core::LinAlg::Matrix<3, 3> DFDx_matrix_transposed;
    DFDx_matrix_transposed.update_t(1.0, DFDx_matrix, 0.0);
    Core::LinAlg::Matrix<9, 1> DFDx_transposed_ref;
    Core::LinAlg::Voigt::matrix_3x3_to_9x1(DFDx_matrix_transposed, DFDx_transposed_ref);

    FOUR_C_EXPECT_NEAR(DFDx_transposed, DFDx_transposed_ref, 1.0e-10);
  }

  TEST_F(MultiplicativeSplitDefgradElastHyperTest, TestEvaluateLinearizationODWithoutGrowth)
  {
    // multiplicative split material whose inelastic deformation gradient does not depend on the
    // scalar
    const int no_growth_id(3);
    Global::Problem::instance()->materials()->insert(
        no_growth_id, Mat::make_parameter(1, Core::Materials::MaterialType::mfi_no_growth, {}));

    Core::IO::InputParameterContainer no_growth_material_data;
    no_growth_material_data.add("MATID_ELAST", 4);
    no_growth_material_data.add("NUMFACINEL", 1);
    no_growth_material_data.add("INELDEFGRADFACIDS", std::vector<int>{no_growth_id});
    no_growth_material_data.add("REF_TEMPERATURE", 293.0);
    no_growth_material_data.add("THERMAL_EXPANSION_COEFFICIENT", 0.1);
    Mat::PAR::MultiplicativeSplitDefgradElastHyper no_growth_parameters(
        Core::Mat::PAR::Parameter::Data{.parameters = no_growth_material_data});
    Mat::MultiplicativeSplitDefgradElastHyper no_growth_material(&no_growth_parameters);

    // evaluate the material with growth first, its linearization must not affect the one of the
    // material without growth
    const double dummy_conc(-1.0);
    set_concentration_to_inelastic_material(dummy_conc);
    const double concentration(1.0);
    Core::LinAlg::Matrix<9, 1> DFDx_growth(Core::LinAlg::Initialization::zero);
    multiplicative_split_defgrad_->evaluate_linearization_od(F_, concentration, DFDx_growth);

    Core::LinAlg::Matrix<9, 1> DFDx(Core::LinAlg::Initialization::zero);
    no_growth_material.evaluate_linearization_od(F_, concentration, DFDx);

    const Core::LinAlg::Matrix<9, 1> DFDx_ref(Core::LinAlg::Initialization::zero);
    FOUR_C_EXPECT_NEAR(DFDx, DFDx_ref, 1.0e-14);
  }

  TEST_F(MultiplicativeSplitDefgradElastHyperTest, TestEvaluateOdStiffMat)
  {
    const double concentration(44327.362);
    set_concentration_to_inelastic_material(concentration);

    // do the actual call that is tested
    auto source(Mat::PAR::InelasticSource::concentration);

    // reference solution
    Core::LinAlg::Matrix<6, 1> dSdx_ref;
    dSdx_ref(0) = -1.907155639254611e-05;
    dSdx_ref(1) = -1.409683812529051e-05;
    dSdx_ref(2) = -1.05352513901749e-05;
    dSdx_ref(3) = 1.14576920347856e-06;
    dSdx_ref(4) = 1.08343997650926e-06;
    dSdx_ref(5) = 1.949130554546719e-06;

    Core::LinAlg::Matrix<6, 1> dSdx =
        multiplicative_split_defgrad_->evaluate_od_stiff_mat(source, &FM_, dSdiFin_ref_);

    FOUR_C_EXPECT_NEAR(dSdx, dSdx_ref, 1.0e-10);
  }
}  // namespace
