// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_MAT_MULTIPLICATIVE_SPLIT_DEFGRAD_ELASTHYPER_HPP
#define FOUR_C_MAT_MULTIPLICATIVE_SPLIT_DEFGRAD_ELASTHYPER_HPP

#include "4C_config.hpp"

#include "4C_comm_parobjectfactory.hpp"
#include "4C_linalg_fixedsizematrix.hpp"
#include "4C_linalg_symmetric_tensor.hpp"
#include "4C_linalg_tensor.hpp"
#include "4C_mat_multiplicative_split_defgrad_elasthyper_service.hpp"
#include "4C_mat_so3_material.hpp"
#include "4C_mat_trait_thermo_solid.hpp"
#include "4C_material_parameter_base.hpp"
#include "4C_utils_exceptions.hpp"

#include <memory>
#include <vector>

FOUR_C_NAMESPACE_OPEN


// forward declaration
namespace Mat
{
  class InelasticDefgradFactors;

  namespace PAR
  {
    enum class InelasticSource
    {
      none,
      concentration,
      temperature
    };

    class MultiplicativeSplitDefgradElastHyper : public Core::Mat::PAR::Parameter
    {
     public:
      /// standard constructor
      explicit MultiplicativeSplitDefgradElastHyper(const Core::Mat::PAR::Parameter::Data& matdata);

      std::shared_ptr<Core::Mat::Material> create_material() override;

      /// ID of the hyperelastic material describing the elastic response
      const int matid_elast_;

      /// number of factors of inelastic deformation gradient F_{in} = F_{in,1} . F_{in,2}. ... .
      /// F_{in,n} (n factors)
      const int numfac_inel_;

      /// IDs of inelastic deformation gradient factors (i-th ID specifies calculation of F_{in,i})
      const std::vector<int> inel_defgradfacids_;

      /// reference temperature for thermal expansion
      const double ref_temperature_;

      /// thermal expansion coefficient
      const double thermal_expansion_coefficient_;

    };  // class MultiplicativeSplitDefgrad_ElastHyper

  }  // namespace PAR

  class MultiplicativeSplitDefgradElastHyperType : public Core::Communication::ParObjectType
  {
   public:
    std::string name() const override { return "MultiplicativeSplitDefgrad_ElastHyperType"; }

    static MultiplicativeSplitDefgradElastHyperType& instance() { return instance_; };

    Core::Communication::ParObject* create(Core::Communication::UnpackBuffer& buffer) override;

   private:
    static MultiplicativeSplitDefgradElastHyperType instance_;
  };

  /*----------------------------------------------------------------------*/
  /*! \class InelasticFactorsHandler

     All factors (class InelasticDefgradFactors) contributing to inelastic deformation are stored
     within this class. They are initially classified according to the source of inelastic
     deformation (e.g. concentration). All contributions from one source can be returned.
  */
  class InelasticFactorsHandler
  {
   public:
    /*!
     * @brief Evaluate inelastic deformation gradient from all factors and assign to different
     * sources
     *
     * @param[in]  defgrad  Deformation gradient
     * @param[out] iFinM    inverse inelastic deformation gradient
     */
    void evaluate_inverse_inelastic_def_grad(
        const Core::LinAlg::Matrix<3, 3>* defgrad, Core::LinAlg::Matrix<3, 3>& iFinM);

    /// Returns all inelastic factors as a vector
    const std::vector<
        std::pair<PAR::InelasticSource, std::shared_ptr<Mat::InelasticDefgradFactors>>>&
    fac_def_grad_in() const
    {
      return facdefgradin_;
    }

    /// Return vector of all inelastic deformation gradients
    const std::vector<std::pair<PAR::InelasticSource, Core::LinAlg::Matrix<3, 3>>>& geti_finj()
        const
    {
      return i_finj_;
    }

    /// total number of inelastic contributions
    int num_inelastic_def_grad() const { return static_cast<int>(facdefgradin_.size()); }

    /// Assigns the different inelastic factors to different sources
    void assign_to_source(
        Mat::PAR::MultiplicativeSplitDefgradElastHyper* params, Mat::So3Material& elastic_material);

    /// Update the history variables for the different inelastic factors
    void update();

    /// Set initial conditions for the history variables of all Gauss points
    void setup(const int numgp, const Discret::Elements::Fibers& fibers,
        const std::optional<Discret::Elements::CoordinateSystem>& coord_system);

    /// Pack
    void pack_inelastic(Core::Communication::PackBuffer& data) const;

    /// Unpack
    void unpack_inelastic(Core::Communication::UnpackBuffer& buffer);

    /*!
     * @brief Register names of the internal data that should be saved during runtime output
     *
     * @param[out] name_and_size Unordered map of names of the data with the respective vector size
     */
    void register_output_data_names(std::unordered_map<std::string, int>& names_and_size) const;

    /*!
     * @brief Evaluate internal data for every Gauss point saved for output during runtime
     * output
     *
     * @param[in] name  Name of the data to export
     * @param[out] data NUMGPxNUMDATA Matrix holding the data
     *
     * @return true if data is set by the material, otherwise false
     */
    bool evaluate_output_data(const std::string& name, Core::LinAlg::SerialDenseMatrix& data) const;

   private:
    /// vector that holds pairs of inelastic contribution and respective source
    std::vector<std::pair<PAR::InelasticSource, std::shared_ptr<Mat::InelasticDefgradFactors>>>
        facdefgradin_;

    /// vector that holds pairs of inelastic deformation gradients and respective source
    std::vector<std::pair<PAR::InelasticSource, Core::LinAlg::Matrix<3, 3>>> i_finj_;
  };

  /*----------------------------------------------------------------------*/
  /*! @class MultiplicativeSplitDefgrad_ElastHyper

    In this class the deformation gradient is modeled to be split multiplicatively
    in elastic and inelastic deformation gradients (F = F_{el} * F_{in}).
    The elastic response is evaluated by a hyperelastic material (MAT_ElastHyper) at the elastic
    deformation gradient F_{el} = F * F_{in}^{-1} and pulled back to the reference configuration.
    Only elastic strains cause stresses!
    The inelastic deformation gradient itself can be a product of different inelastic
    deformation gradients, i.e. F_{in} = F_{in,1} * F_{in,2} * ... * F_{in,n}.
    The inverse inelastic deformation gradients and their derivative w.r.t. the primary variables
    that are needed to set up the system to be solved are evaluated in the derived classes
    of the interface class 'InelasticDefgradFactors'.
*/
  class MultiplicativeSplitDefgradElastHyper : public Trait::ThermoSolid
  {
   public:
    /// construct empty material object
    MultiplicativeSplitDefgradElastHyper();

    /// construct the material object given material parameters
    explicit MultiplicativeSplitDefgradElastHyper(
        Mat::PAR::MultiplicativeSplitDefgradElastHyper* params);

    /// stress and its derivatives, obtained by pulling back the response of the elastic material
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

    int unique_par_object_id() const override
    {
      return MultiplicativeSplitDefgradElastHyperType::instance().unique_par_object_id();
    }

    void pack(Core::Communication::PackBuffer& data) const override;

    void unpack(Core::Communication::UnpackBuffer& buffer) override;

    void valid_kinematics(Solid::KinemType kinem) override
    {
      if (kinem != Solid::KinemType::nonlinearTotLag)
        FOUR_C_THROW("element and material kinematics are not compatible");
    }

    Core::Materials::MaterialType material_type() const override
    {
      return Core::Materials::m_multiplicative_split_defgrad_elasthyper;
    }

    std::shared_ptr<Core::Mat::Material> clone() const override
    {
      return std::make_shared<MultiplicativeSplitDefgradElastHyper>(*this);
    }

    /// material mass density (defined by the elastic material)
    double density() const override { return elastic_material_->density(); }

    Core::Mat::PAR::Parameter* parameter() const override { return params_; }

    void evaluate(const Core::LinAlg::Tensor<double, 3, 3>* defgrad,
        const Core::LinAlg::SymmetricTensor<double, 3, 3>& glstrain,
        const Teuchos::ParameterList& params, const EvaluationContext<3>& context,
        Core::LinAlg::SymmetricTensor<double, 3, 3>& stress,
        Core::LinAlg::SymmetricTensor<double, 3, 3, 3, 3>& cmat, int gp, int eleGID) override;

    Core::LinAlg::SymmetricTensor<double, 3, 3> evaluate_d_stress_d_scalar(
        const Core::LinAlg::Tensor<double, 3, 3>& defgrad,
        const Core::LinAlg::SymmetricTensor<double, 3, 3>& glstrain,
        const Teuchos::ParameterList& params, const EvaluationContext<3>& context, int gp,
        int eleGID) override;

    void reinit(double temperature, unsigned gp) override { /* do nothing */ };

    void stress_temperature_modulus_and_deriv(Core::LinAlg::SymmetricTensor<double, 3, 3>& stm,
        Core::LinAlg::SymmetricTensor<double, 3, 3>& stm_dT, const int gp) override;

    /*!
     * @brief Evaluate the heat source produced by this material. Currently, this is only
     * enabled for the viscoplastic factor.
     *
     * @param[in] context Evaluation context, providing access to the current timestep and total
     * time
     * @param[in] gp Gauss point
     * @param[in] eleGID global element ID
     * @param[in] defgrad Deformation gradient
     * @param[in] current_temperature Absolute current temperature
     * @return mechanical dissipation heat source and derivatives w.r.t. temperature and the right
     *         Cauchy-Green tensor
     */
    [[nodiscard]] HeatSource evaluate_additional_heat_source(const EvaluationContext<3>& context,
        const int gp, const int eleGID, const Core::LinAlg::Matrix<3, 3>* current_defgrad,
        const double current_temperature);

    /*!
     * @brief Evaluate the Cauchy stress contracted with n and dir and its derivatives
     *
     * @note The derivatives w.r.t. the deformation gradient are taken at fixed inelastic
     * deformation gradient, i.e., its dependence on the deformation gradient is neglected.
     */
    double evaluate_cauchy_n_dir_and_derivatives(const Core::LinAlg::Tensor<double, 3, 3>& defgrd,
        const Core::LinAlg::Tensor<double, 3>& n, const Core::LinAlg::Tensor<double, 3>& dir,
        Core::LinAlg::Matrix<3, 1>* d_cauchyndir_dn, Core::LinAlg::Matrix<3, 1>* d_cauchyndir_ddir,
        Core::LinAlg::Matrix<9, 1>* d_cauchyndir_dF, Core::LinAlg::Matrix<9, 9>* d2_cauchyndir_dF2,
        Core::LinAlg::Matrix<9, 3>* d2_cauchyndir_dF_dn,
        Core::LinAlg::Matrix<9, 3>* d2_cauchyndir_dF_ddir, const EvaluationContext<3>& context,
        int eleGID, const double* concentration, const double* temp, double* d_cauchyndir_dT,
        Core::LinAlg::Matrix<9, 1>* d2_cauchyndir_dF_dT) override;

    void evaluate_linearization_od(const Core::LinAlg::Tensor<double, 3, 3>& defgrd,
        double concentration, Core::LinAlg::Matrix<9, 1>& d_F_dx) override;

    void setup(int numgp, const Discret::Elements::Fibers& fibers,
        const std::optional<Discret::Elements::CoordinateSystem>& coord_system) override;

    void post_setup(const Teuchos::ParameterList& params, int eleGID) override;

    void update() override;

    void register_output_data_names(
        std::unordered_map<std::string, int>& names_and_size) const override;

    bool evaluate_output_data(
        const std::string& name, Core::LinAlg::SerialDenseMatrix& data) const override;

    /*!
     * @brief Evaluate off-diagonal stiffness matrix contributions due to the inelastic deformation,
     * i.e.:
     * \f[
     * \frac{\partial \mathbf{S}}{\partial \mathbf{F}_\text{in}^{-1}}
     * :\frac{\mathrm{d}\mathbf{F}_\text{in}^{-1}}{\mathrm{d}x}
     * \f]
     * where \f[x\f] is the scalar of the secondary field corresponding to \param source
     *
     * @param[in] source   Source of inelastic deformation
     * @param[in] defgrad  Deformation gradient
     * @param[in] dSdiFin  Already computed derivative of 2nd Piola Kirchhoff stresses w.r.t. the
     * inverse inelastic deformation gradient \f[ \frac{\partial \mathbf{S}}{\partial
     * \mathbf{F}_\text{in}^{-1}} \f]
     * @return off-diagonal stiffness matrix contribution in stress-like Voigt notation
     */
    Core::LinAlg::Matrix<6, 1> evaluate_od_stiff_mat(PAR::InelasticSource source,
        const Core::LinAlg::Matrix<3, 3>* defgrad, const Core::LinAlg::Matrix<6, 9>& dSdiFin);

    /*!
     * @brief Evaluate additional terms of the elasticity tensor
     *
     * \f[
     *   cmat_{add} = \frac{\partial\mathbf{S}}{\partial \mathbf{F}_\text{in}^{-1}}
     *             : \frac{\partial\mathbf{F}_\text{in}^{-1}}{\partial \mathbf{C}}
     * \f]
     *
     * @param[in] defgrad  Deformation gradient
     * @param[in] iCV      Inverse right Cauchy-Green tensor stored as 6x1 vectors
     * @param[in] dSdiFin  Derivative of 2nd Piola Kirchhoff stress w.r.t. the inverse inelastic
     *                     deformation gradient
     * @return Additional elasticity tensor \f$ cmat_{add} \f$
     */
    Core::LinAlg::Matrix<6, 6> evaluate_additional_cmat(const Core::LinAlg::Matrix<3, 3>* defgrad,
        const Core::LinAlg::Matrix<6, 1>& iCV, const Core::LinAlg::Matrix<6, 9>& dSdiFin);

    /*!
     * @brief pre-evaluation, intended to be used for stuff that have to be done only once
     *
     * @param[in] params  parameter list as handed in from the element
     * @param[in] gp      current gauss point
     * @param[in] eleGID  Element ID
     */
    void pre_evaluate(const Teuchos::ParameterList& params, const EvaluationContext<3>& context,
        int gp, int eleGID) const;

    /*!
     * @brief set the gauss point concentration to the respective parameter class of the inelastic
     * material
     *
     * @param[in]  concentration concentration at gauss point
     */
    void set_concentration_gp(double concentration);

   private:
    /*!
     * @brief Evaluate the 2nd Piola-Kirchhoff stress and its derivatives, while keeping the
     * inelastic deformation gradient fixed.
     *
     * The complete 2nd Piola Kirchhoff stress is
     * \f[
     *   \mathbf{S} = \det(\mathbf{F}_\text{in}) \mathbf{F}_\text{in}^{-1} \cdot
     *   (\mathbf{S}_\text{e} - \mathbf{S}_T) \cdot \mathbf{F}_\text{in}^{-T}
     * \f]
     * where \f$ \mathbf{S}_\text{e} \f$ is the hyperelastic 2nd Piola-Kirchhoff stress, obtained by
     * evaluating the hyperelastic material at the elastic right Cauchy-Green tensor
     * \f$ \mathbf{C}_\text{e} = \mathbf{F}_\text{in}^{-T} \cdot \mathbf{C} \cdot
     * \mathbf{F}_\text{in}^{-1} \f$. Its hyperelastic stiffness
     * \f$ \mathbb{C}_\text{e} = 2 \frac{\partial \mathbf{S}_\text{e}}{\partial \mathbf{C}_\text{e}}
     * \f$ is also provided by the hyperelastic material.
     * \f$ \mathbf{S}_T\f$ is the thermo-elastic stress, obtained by evaluating the hyperelastic
     * material at the thermal right Cauchy-Green tensor
     * \f$ \mathbf{C}_T = (1 + 2 \alpha_T \Delta T) \mathbf{I} \f$.
     *
     * This yields the temperature derivative at fixed \f$ \mathbf{F}_\text{in}^{-1} \f$:
     * \f[
     * \frac{\partial \mathbf{S}}{\partial T} = -\det(\mathbf{F}_\text{in})
     *   \mathbf{F}_\text{in}^{-1} \cdot \frac{\partial \mathbf{S}_T}{\partial T}
     *   \cdot \mathbf{F}_\text{in}^{-T}
     * \f]
     * The tangent stiffness tensor at fixed \f$ \mathbf{F}_\text{in}^{-1} \f$ is given by the
     * hyperelastic stiffness tensor pulled back in its first index pair and applying the chain
     * rule:
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
     *     \delta_{ik} \left[ \mathbf{F}_\text{in}^{-1} \cdot
     *     \left(\mathbf{S}_\text{e} - \mathbf{S}_T \right) \right]_{jl}
     *   + \delta_{jk} \left[ \mathbf{F}_\text{in}^{-1} \cdot
     *     \left(\mathbf{S}_\text{e} - \mathbf{S}_T \right) \right]_{il} \right)
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
     * @param[in] defgrad  Deformation gradient
     * @param[in] iFin     Inverse inelastic deformation gradient
     * @param[in] params   Parameter list as handed in from the element
     * @param[in] context  Evaluation context as handed in from the element
     * @param[in] gp       Gauss point
     * @param[in] eleGID   Element ID
     */
    [[nodiscard]] StressAndDerivatives evaluate_stress_and_derivatives(
        const Core::LinAlg::Tensor<double, 3, 3>& defgrad,
        const Core::LinAlg::Tensor<double, 3, 3>& iFin, const Teuchos::ParameterList& params,
        const EvaluationContext<3>& context, int gp, int eleGID);

    /// Holds and classifies all inelastic factors
    std::shared_ptr<InelasticFactorsHandler> inelastic_;

    /// My material parameters
    Mat::PAR::MultiplicativeSplitDefgradElastHyper* params_;

    /// hyperelastic material describing the elastic response
    std::shared_ptr<Mat::So3Material> elastic_material_;
  };

}  // namespace Mat
FOUR_C_NAMESPACE_CLOSE

#endif
