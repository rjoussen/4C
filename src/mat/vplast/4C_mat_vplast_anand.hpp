// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_MAT_VPLAST_ANAND_HPP
#define FOUR_C_MAT_VPLAST_ANAND_HPP
#include "4C_config.hpp"

#include "4C_comm_parobject.hpp"
#include "4C_comm_parobjectfactory.hpp"
#include "4C_linalg_fixedsizematrix.hpp"
#include "4C_mat_inelastic_defgrad_factors_service.hpp"
#include "4C_mat_vplast_law.hpp"
#include "4C_material_parameter_base.hpp"
#include "4C_utils_exceptions.hpp"
#include "4C_utils_parameter_list.fwd.hpp"

#include <Teuchos_RCP.hpp>

#include <cmath>
#include <memory>


FOUR_C_NAMESPACE_OPEN


namespace Mat
{

  namespace Viscoplastic
  {

    namespace PAR
    {
      /*----------------------------------------------------------------------*/
      /*! \class Anand
       *
       * Parameter class of Anand
       */
      class Anand : public Core::Mat::PAR::Parameter
      {
       public:
        explicit Anand(const Core::Mat::PAR::Parameter::Data& matdata);

        std::shared_ptr<Core::Mat::Material> create_material() override { return nullptr; };

        // getter methods
        //! get strain rate prefactor \f$ A \exp\left( - \frac{Q}{R T} \right) \f$
        [[nodiscard]] double strain_rate_pre_fac() const { return strain_rate_pre_fac_; }
        //! get strain rate sensitivity \f$ m \f$
        [[nodiscard]] double strain_rate_sensitivity() const { return strain_rate_sensitivity_; }
        //! get initial flow resistance \f$ S(t_0) \f$
        [[nodiscard]] double init_flow_res() const { return init_flow_res_; }
        //! get hardening rate sensitivity \f$ a \f$
        [[nodiscard]] double harden_rate_sensitivity() const { return harden_rate_sensitivity_; }
        //! get hardening rate prefactor \f$ H_0 \f$
        [[nodiscard]] double harden_rate_pre_fac() const { return harden_rate_pre_fac_; }
        //! get flow resistance saturation prefactor \f$ S_* \f$
        [[nodiscard]] double flow_res_sat_fac() const { return flow_res_sat_fac_; }
        //! get flow resistance saturation exponent \f$ N \f$
        [[nodiscard]] double flow_res_sat_exp() const { return flow_res_sat_exp_; }


       private:
        //! strain rate prefactor (computed from the preexponential factor \f$ A \f$ and the
        //! exponential term \f$ \exp\left( - \frac{Q}{R T}\right) \f$)
        const double strain_rate_pre_fac_;
        //! strain rate sensitivity \f$ m \f$
        const double strain_rate_sensitivity_;
        //! initial flow resistance \f$ S(t_0) \f$
        const double init_flow_res_;
        //! hardening rate sensitivity \f$ a \f$
        const double harden_rate_sensitivity_;
        //! hardening rate pre-factor \f$ H_0 \f$
        const double harden_rate_pre_fac_;
        //! flow resistance saturation prefactor \f$ S_* \f$
        const double flow_res_sat_fac_;
        //! flow resistance saturation exponent / saturation sensitivity \f$ N \f$
        const double flow_res_sat_exp_;
      };
    }  // namespace PAR

    /*----------------------------------------------------------------------*/
    /*! \class Anand
     *
     *  Viscoplasticity law associated to the Anand model presented in
     *  -# Anand et al., An Elastic-Viscoplastic Model for Lithium, J. Electrochem. Soc. 166, 2019
     * Further details on the implementation in the numerical framework of the viscoplasticity laws
     * are presented in:
     * -# Master's Thesis : Dragos-Corneliu Ana, Continuum Modeling and Calibration of
     * Viscoplasticity in the Context of the Lithium Anode in Solid State Batteries,
     * Supervisor: Christoph Schmidt, 2024
     */
    class Anand : public Law
    {
     public:
      explicit Anand(Core::Mat::PAR::Parameter* params);

      Mat::Viscoplastic::PAR::Anand* parameter() const override
      {
        return dynamic_cast<Mat::Viscoplastic::PAR::Anand*>(Mat::Viscoplastic::Law::parameter());
      }

      Core::Materials::MaterialType material_type() const override
      {
        return Core::Materials::mvl_Anand;
      };

      double evaluate_stress_ratio(
          const double equiv_stress, const double equiv_plastic_strain) override;

      double evaluate_plastic_strain_rate(const double equiv_stress,
          const double equiv_plastic_strain, const double dt, const double max_plastic_strain_incr,
          Mat::ViscoplastErrorType& err_status, const bool update_hist_var) override;

      Core::LinAlg::Matrix<2, 1> evaluate_derivatives_of_plastic_strain_rate(
          const double equiv_stress, const double equiv_plastic_strain, const double dt,
          const double max_plastic_strain_deriv, Mat::ViscoplastErrorType& err_status,
          const bool update_hist_var) override;

      void setup(const int numgp, const Core::IO::InputParameterContainer& container) override;

      void pre_evaluate(int gp) override;

      void update() override;

      void update_gp_state(int gp) override;

      void pack_viscoplastic_law(Core::Communication::PackBuffer& data) const override;

      void unpack_viscoplastic_law(Core::Communication::UnpackBuffer& buffer) override;

      std::string debug_get_error_info(int gp) override;

     private:
      /// current Gauss point
      int gp_;

      /// struct containing constant parameters to be evaluated only once
      struct ConstPars
      {
        /// prefactor \f$ P \f$ of the plastic strain rate
        double p;

        /// prefactor logarithm \f$ \log(P) \f$
        double log_p;

        /// exponent \f$ E \f$ of the plastic strain rate
        double e;

        /// \f$ \log(P*E) \f$
        double log_p_e;

        /// hardening rate prefactor \f$ H_0 \f$
        double H_0;

        /// hardening rate sensitivity \f$ a \f$
        double a;

        /// \f$ a H_0 \f$
        double aH_0;

        /// flow resistance saturation prefactor \f$ S_* \f$
        double S_star;

        /// \f$ 1.0 / S_*  \f$
        double inv_S_star;

        /// flow resistance saturation exponent \f$ N \f$
        double N;

        /// \f$ e N \f$
        double eN;

        /// initial flow resistance
        double S_0;

        /// constructor
        ConstPars(const double prefac, const double expon, const double init_flow_resistance,
            const double harden_rate_prefac, const double harden_rate_sensitive,
            const double flow_res_sat_prefac, const double flow_res_sat_exp)
            : p(prefac),
              log_p(std::log(p)),
              e(expon),
              log_p_e(std::log(prefac * expon)),
              H_0(harden_rate_prefac),
              a(harden_rate_sensitive),
              aH_0(harden_rate_prefac * harden_rate_sensitive),
              S_star(flow_res_sat_prefac),
              inv_S_star(1.0 / flow_res_sat_prefac),
              N(flow_res_sat_exp),
              eN(expon * flow_res_sat_exp),
              S_0(init_flow_resistance)
        {
        }
      };
      /// instance of ConstPars struct
      const ConstPars const_pars_;

      /*! @brief Calculate flow resistance \f$ S \f$ from the equivalent stress \f$
       * \overline{\sigma} \f$ and the equivalent plastic strain \f$
       * \varepsilon^{\text{p}} \f$ (along with the known and already saved flow resistance, plastic
       * strain and equivalent stress values at the previous time instant)
       *
       * @note The computation is performed by discretizing the hardening
       * rate equation \f$ \dot{S}( \overline{\sigma}, \varepsilon^{\text{p}},
       * \dot{\varepsilon}^{\text{p}(\overline{\sigma}, \varepsilon^{\text{p}})) \f$ in time.
       * Currently, we only use the implicit Euler method for consistency with the discretization
       * within InelasticDefgradTransvIsotropElastViscoplast.
       *
       * @param[in] equiv_stress equivalent tensile stress \f$ \overline{\sigma} \f$ at the new time
       * instant
       * @param[in] equiv_plastic_strain equivalent plastic strain \f$
       * \varepsilon^{\text{p}} \f$ at the new time instant
       * @param[in|out] err_status: error status of the evaluation
       * return flow resistance at the new time instant
       */
      double compute_flow_resistance(const double equiv_stress, const double equiv_plastic_strain,
          Mat::ViscoplastErrorType& err_status);

      /*! @brief Calculate the partial derivatives of the flow resistance \f$ S \f$ with respect to
       * the equivalent stress \f$ \overline{\sigma}  \f$ and the plastic strain \f$
       * \varepsilon^{\text{p}} \f$ (assuming constant equivalent stress and flow resistance at the
       * previous time instant)
       *
       *
       * @param[in] equiv_stress equivalent tensile stress \f$ \overline{\sigma} \f$
       * @param[in] flow_resistance flow resistance \f$ S \f$
       * @param[in] delta_plastic_strain plastic strain increment \f$ \varepsilon^{\text{p}} -
       * \varepsilon^{\text{p}}_{\text{last}}  \f$, with \f$ \varepsilon^{\text{p}}_{\text{last}}
       * \f$ the plastic strain at the previous time instant (after the last converged substep, when
       * substepping)
       * @param[in|out] err_status: error status of the evaluation
       * @return derivatives of the flow resistance (component 0: w.r.t. equivalent stress,
       * component 1: w.r.t. plastic strain)
       */
      Core::LinAlg::Matrix<2, 1> compute_derivatives_of_flow_resistance(const double equiv_stress,
          const double flow_resistance, const double delta_plastic_strain,
          Mat::ViscoplastErrorType& err_status);

      /*! @brief Calculate the hardening tangent \f$ \frac{\text{d} S}{ \text{d}
       * \varepsilon^{\text{p}}}\f$, derived from the hardening rate equation
       *
       * @param[in] equiv_stress equivalent tensile stress \f$ \overline{\sigma} \f$
       * @param[in] flow_resistance flow resistance \f$ S \f$
       * @param[in|out] err_status: error status of the evaluation
       * @return hardening tangent value
       */
      double compute_hardening_tangent(const double equiv_stress, const double flow_resistance,
          Mat::ViscoplastErrorType& err_status);

      /*! @brief Calculate the partial derivatives of the hardening tangent \f$ \frac{\text{d} S}{
       * \text{d} \varepsilon^{\text{p}}}\f$ with respect to the equivalent stress \f$
       * \overline{\sigma}  \f$ and the flow resistance \f$ S \f$
       *
       * @param[in] equiv_stress equivalent tensile stress \f$ \overline{\sigma} \f$
       * @param[in] flow_resistance flow resistance \f$ S \f$
       * @param[in|out] err_status: error status of the evaluation
       * @return derivatives of the hardening tangent (component 0: w.r.t. equivalent stress,
       * component 1: w.r.t. flow resistance)
       */
      Core::LinAlg::Matrix<2, 1> compute_derivatives_of_hardening_tangent(const double equiv_stress,
          const double flow_resistance, Mat::ViscoplastErrorType& err_status);


      //! struct containing quantities at the last and current time points (i.e., at \f[ t_n \f] and
      //! \f[ t_{n+1} \f], respectively). The quantities are tracked at all Gauss points, in order
      //! to update them simultaneously during the update method call
      struct TimeStepQuantities
      {
        //! flow resistance at the last time step (for all Gauss points)
        std::vector<double> last_flow_resistance_;

        //! currently evaluated flow resistance (for all Gauss points)
        std::vector<double> current_flow_resistance_;

        //! flow resistance at the last computed time instant (after the last converged substep when
        //! using a substepping procedure, for all Gauss points)
        std::vector<double> last_substep_flow_resistance_;

        //! plastic strain at the last time step (for all Gauss points)
        std::vector<double> last_plastic_strain_;

        //! currently evaluated plastic strain (for all Gauss points)
        std::vector<double> current_plastic_strain_;

        //! plastic strain at the last computed time instant (after the last converged substep when
        //! using a substepping procedure, for all Gauss points)
        std::vector<double> last_substep_plastic_strain_;
      };
      TimeStepQuantities time_step_quantities_;
    };

  }  // namespace Viscoplastic

}  // namespace Mat

FOUR_C_NAMESPACE_CLOSE

#endif
