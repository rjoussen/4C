// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_mat_vplast_anand.hpp"

#include "4C_comm_pack_buffer.hpp"
#include "4C_global_data.hpp"
#include "4C_linalg_fixedsizematrix.hpp"
#include "4C_mat_inelastic_defgrad_factors_service.hpp"
#include "4C_mat_par_bundle.hpp"
#include "4C_mat_vplast_law.hpp"
#include "4C_utils_exceptions.hpp"
#include "4C_utils_fad.hpp"

#include <cmath>
#include <string>
#include <unordered_map>
#include <utility>

FOUR_C_NAMESPACE_OPEN

namespace
{

}  // namespace

using namespace Mat::InelasticDefgradTransvIsotropElastViscoplastUtils;

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
Mat::Viscoplastic::PAR::Anand::Anand(const Core::Mat::PAR::Parameter::Data& matdata)
    : Parameter(matdata),
      strain_rate_pre_fac_(matdata.parameters.get<double>("STRAIN_RATE_PREFAC")),
      strain_rate_sensitivity_(matdata.parameters.get<double>("STRAIN_RATE_SENS")),
      init_flow_res_(matdata.parameters.get<double>("INIT_FLOW_RES")),
      harden_rate_sensitivity_(matdata.parameters.get<double>("HARDEN_RATE_SENS")),
      harden_rate_pre_fac_(matdata.parameters.get<double>("HARDEN_RATE_PREFAC")),
      flow_res_sat_fac_(matdata.parameters.get<double>("FLOW_RES_SAT_FAC")),
      flow_res_sat_exp_(matdata.parameters.get<double>("FLOW_RES_SAT_EXP"))
{
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
Mat::Viscoplastic::Anand::Anand(Core::Mat::PAR::Parameter* params)
    : Mat::Viscoplastic::Law(params),
      const_pars_(parameter()->strain_rate_pre_fac(), 1.0 / parameter()->strain_rate_sensitivity(),
          parameter()->init_flow_res(), parameter()->harden_rate_pre_fac(),
          parameter()->harden_rate_sensitivity(), parameter()->flow_res_sat_fac(),
          parameter()->flow_res_sat_exp())
{
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
void Mat::Viscoplastic::Anand::setup(
    const int numgp, const Core::IO::InputParameterContainer& container)
{
  // update last_ and current_ values of the flow resistance
  time_step_quantities_.last_flow_resistance_.resize(numgp, parameter()->init_flow_res());
  time_step_quantities_.current_flow_resistance_.resize(
      numgp, parameter()->init_flow_res());  // value irrelevant at this point in time
  time_step_quantities_.last_substep_flow_resistance_.resize(numgp, parameter()->init_flow_res());


  // update last_ and current_ values of the plastic strain
  time_step_quantities_.last_plastic_strain_.resize(numgp, 0.0);
  time_step_quantities_.current_plastic_strain_.resize(
      numgp, 0.0);  // value irrelevant at this point in time
  time_step_quantities_.last_substep_plastic_strain_.resize(numgp, 0.0);
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
void Mat::Viscoplastic::Anand::pre_evaluate(int gp)
{  // set current Gauss point
  gp_ = gp;

  // set the last_substep_ values to the last_ values (because the preevaluation method is called
  // before substepping)
  time_step_quantities_.last_substep_flow_resistance_[gp] =
      time_step_quantities_.last_flow_resistance_[gp];
  time_step_quantities_.last_substep_plastic_strain_[gp] =
      time_step_quantities_.last_plastic_strain_[gp];
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
void Mat::Viscoplastic::Anand::update()
{
  // update last_ values with the current_ values
  time_step_quantities_.last_flow_resistance_ = time_step_quantities_.current_flow_resistance_;
  time_step_quantities_.last_substep_flow_resistance_ =
      time_step_quantities_.current_flow_resistance_;
  time_step_quantities_.last_plastic_strain_ = time_step_quantities_.current_plastic_strain_;
  time_step_quantities_.last_substep_plastic_strain_ =
      time_step_quantities_.current_plastic_strain_;
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
void Mat::Viscoplastic::Anand::update_gp_state(int gp)
{
  // the last_substep_ values are updated with the
  // currently evaluated values after a successfully converged substep (at the currently evaluated
  // GP only)
  time_step_quantities_.last_substep_flow_resistance_[gp] =
      time_step_quantities_.current_flow_resistance_[gp];
  time_step_quantities_.last_substep_plastic_strain_[gp] =
      time_step_quantities_.current_plastic_strain_[gp];
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
void Mat::Viscoplastic::Anand::pack_viscoplastic_law(Core::Communication::PackBuffer& data) const
{
  // pack history variables of this specific inelastic factor
  if (parameter() != nullptr)
  {
    // pack last_ values inside time_step_quantities_
    add_to_pack(data, time_step_quantities_.last_flow_resistance_);
    add_to_pack(data, time_step_quantities_.last_plastic_strain_);
  }
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
void Mat::Viscoplastic::Anand::unpack_viscoplastic_law(Core::Communication::UnpackBuffer& buffer)
{
  // NOTE: factory method is called during assign_to_source in the unpack method of the
  // multiplicative split framework --> the param class of the viscoplastic law is created, we only
  // need to unpack the history variables
  if (parameter() != nullptr)
  {
    // unpack last_ values inside time_step_quantities_
    extract_from_pack(buffer, time_step_quantities_.last_flow_resistance_);
    extract_from_pack(buffer, time_step_quantities_.last_plastic_strain_);
  }

  // fill current_ and last_substep_ values with the last_ values
  time_step_quantities_.current_flow_resistance_.resize(
      time_step_quantities_.last_flow_resistance_.size(),
      time_step_quantities_.last_flow_resistance_[0]);  // value irrelevant
  time_step_quantities_.last_substep_flow_resistance_.resize(
      time_step_quantities_.last_flow_resistance_.size(),
      time_step_quantities_.last_flow_resistance_[0]);  // value irrelevant
  time_step_quantities_.current_plastic_strain_.resize(
      time_step_quantities_.last_plastic_strain_.size(),
      time_step_quantities_.last_plastic_strain_[0]);  // value irrelevant
  time_step_quantities_.last_substep_plastic_strain_.resize(
      time_step_quantities_.last_plastic_strain_.size(),
      time_step_quantities_.last_plastic_strain_[0]);  // value irrelevant
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
double Mat::Viscoplastic::Anand::evaluate_stress_ratio(
    const double equiv_stress, const double equiv_plastic_strain)
{
  // declare error status (set to NoErrors initially)
  ErrorType err_status = ErrorType::NoErrors;

  // compute flow resistance
  const double flow_resistance =
      compute_flow_resistance(equiv_stress, equiv_plastic_strain, err_status);
  if (err_status != ErrorType::NoErrors)
    FOUR_C_THROW(
        "This method should not be used and it also failed while computing the flow resistance!");

  return equiv_stress / flow_resistance;
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
double Mat::Viscoplastic::Anand::evaluate_plastic_strain_rate(const double equiv_stress,
    const double equiv_plastic_strain, const double dt, const double max_plastic_strain_incr,
    ErrorType& err_status, const bool update_hist_var)
{
  // first set error status to "no errors"
  err_status = ErrorType::NoErrors;

  // Check if plastic strain is negative and throw error (handled by the parent material,
  // substepping)
  if (equiv_plastic_strain < 0.0)
  {
    err_status = ErrorType::NegativePlasticStrain;
    return -1;
  }

  // compute the viscoplastic strain rate; first we set it to 0
  double equiv_plastic_strain_rate = 0.0;

  // compute flow resistance (and save its value)
  double flow_resistance = time_step_quantities_.last_substep_flow_resistance_[gp_];

  // save the current values
  if (update_hist_var)
  {
    time_step_quantities_.current_flow_resistance_[gp_] = flow_resistance;
    time_step_quantities_.current_plastic_strain_[gp_] = equiv_plastic_strain;
  }

  if (equiv_stress > 0.0)
  {
    flow_resistance = compute_flow_resistance(equiv_stress, equiv_plastic_strain, err_status);
    if (err_status != ErrorType::NoErrors) return -1.0;


    // save the currently computed flow resistance
    if (update_hist_var) time_step_quantities_.current_flow_resistance_[gp_] = flow_resistance;

    // compute stress ratio
    const double log_stress_ratio = std::log(equiv_stress) - std::log(flow_resistance);


    // compute logarithm \f$ \log (P \exp(E \left[\frac{\overline{\sigma}}{S}}
    // - 1.0]) ) \f$
    const double log_temp = const_pars_.log_p + const_pars_.e * log_stress_ratio;

    // check if characteristic term too large, throw error overflow
    // error if so
    if (std::log(dt) + log_temp > std::log(max_plastic_strain_incr))
    {
      err_status = ErrorType::OverflowError;
      return -1;
    }


    equiv_plastic_strain_rate = std::exp(log_temp);
  }

  return equiv_plastic_strain_rate;
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
Core::LinAlg::Matrix<2, 1> Mat::Viscoplastic::Anand::evaluate_derivatives_of_plastic_strain_rate(
    const double equiv_stress, const double equiv_plastic_strain, const double dt,
    const double max_plastic_strain_deriv, ErrorType& err_status, const bool update_hist_var)
{
  // first set error status to "no errors"
  err_status = ErrorType::NoErrors;

  // Check if plastic strain is negative and throw error (handled by the parent material,
  // substepping)
  if (equiv_plastic_strain < 0.0)
  {
    err_status = ErrorType::NegativePlasticStrain;
    return Core::LinAlg::Matrix<2, 1>{true};
  }

  // computation of derivatives

  // first we set derivatives to 0
  Core::LinAlg::Matrix<2, 1> equiv_plastic_strain_rate_ders(true);

  // then we check the yield condition
  if (equiv_stress > 0.0)
  {
    // compute delta in plastic strain w.r.t. last substep value (at the last converged time
    // instant)
    const double delta_equiv_plastic_strain =
        equiv_plastic_strain - time_step_quantities_.last_substep_plastic_strain_[gp_];

    // compute the inverse of the flow resistance (we have already computed and stored the current
    // flow resistance previously)
    const double inv_flow_resistance = 1.0 / time_step_quantities_.current_flow_resistance_[gp_];

    // compute stress ratio
    double stress_ratio = equiv_stress * inv_flow_resistance;


    // compute derivatives of the flow resistance w.r.t. equivalent stress and the plastic strain
    Core::LinAlg::Matrix<2, 1> derivs_of_flow_resistance = compute_derivatives_of_flow_resistance(
        equiv_stress, time_step_quantities_.current_flow_resistance_[gp_],
        delta_equiv_plastic_strain, err_status);
    if (err_status != ErrorType::NoErrors)
    {
      err_status = ErrorType::FailedComputationFlowResistanceDerivs;
      return Core::LinAlg::Matrix<2, 1>{true};
    }

    // get signs of the derivatives
    const int signum_d_flow_res_d_plastic_strain =
        Core::FADUtils::signum(derivs_of_flow_resistance(0, 0));


    // compute the logarithm of the flow resistance derivative w.r.t. plastic strain
    const double log_abs_d_flow_res_d_plastic_strain =
        std::log(std::abs(derivs_of_flow_resistance(1, 0)));


    // logarithm of the stress ratio
    const double log_stress_ratio = std::log(stress_ratio);

    // logarithm of the inverse flow resistance
    const double log_inv_flow_res = std::log(inv_flow_resistance);

    // logarithm of the time step
    const double log_dt = std::log(dt);

    // auxiliary terms
    const double aux_term = const_pars_.log_p_e + (const_pars_.e - 1.0) * log_stress_ratio;
    double temp = std::abs(
        inv_flow_resistance - stress_ratio * inv_flow_resistance * derivs_of_flow_resistance(0, 0));
    const int signum_temp = Core::FADUtils::signum(
        inv_flow_resistance - stress_ratio * inv_flow_resistance * derivs_of_flow_resistance(0, 0));
    temp = std::log(temp);

    // compute first the logarithms of our derivatives (try to avoid overflow!)
    double log_deriv_sigma = aux_term + temp;
    double log_deriv_eps =
        aux_term + log_inv_flow_res + log_stress_ratio + log_abs_d_flow_res_d_plastic_strain;


    // check overflow error using these logarithms
    const double log_max_plastic_strain_deriv = std::log(max_plastic_strain_deriv);
    if ((log_dt + log_deriv_sigma > log_max_plastic_strain_deriv) ||
        (log_dt + log_deriv_eps > log_max_plastic_strain_deriv))
    {
      err_status = ErrorType::OverflowError;
      return Core::LinAlg::Matrix<2, 1>{true};
    }

    // compute the exact derivatives using these logarithms
    equiv_plastic_strain_rate_ders(0, 0) = signum_temp * std::exp(log_deriv_sigma);
    equiv_plastic_strain_rate_ders(1, 0) =
        -signum_d_flow_res_d_plastic_strain * std::exp(log_deriv_eps);
  }

  return equiv_plastic_strain_rate_ders;
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
double Mat::Viscoplastic::Anand::compute_flow_resistance(
    const double equiv_stress, const double equiv_plastic_strain, ErrorType& err_status)
{
  // make sure the error status is on NoErrors
  err_status = ErrorType::NoErrors;

  // integration performed using the Newton-Raphson method

  // purely elastic predictor
  double flow_resistance = time_step_quantities_.last_substep_flow_resistance_[gp_];

  // config of Newton-Raphson loop
  int iter = 0;
  int max_iter = 50;
  double conv_tol = 1.0e-10;
  double res_S = 0;  // residual of the nonlinear equation
  double inv_J = 0;  // inverse Jacobian of the nonlinear equation

  // hardening tangent
  double harden_tang = 0.0;

  // derivatives of the hardening tangent w.r.t. equiv_stress and flow resistance
  Core::LinAlg::Matrix<2, 1> derivs_harden_tang{true};

  // 1.0 / equiv_stress
  const double inv_equiv_stress = 1.0 / equiv_stress;

  // delta of the plastic strain w.r.t. the last substep value
  const double delta_plastic_strain =
      equiv_plastic_strain - time_step_quantities_.last_substep_plastic_strain_[gp_];

  // auxiliary double
  double temp = 0.0;

  // Newton-Raphson loop
  while (true)
  {
    ++iter;  // increment iterator

    // check for negative flow resistance values
    if (flow_resistance < 0.0)
    {
      err_status = ErrorType::FailedComputationFlowResistance;
      return -1;
    }

    // compute hardening tangent
    harden_tang = compute_hardening_tangent(equiv_stress, flow_resistance, err_status);
    if (err_status != ErrorType::NoErrors) return -1.0;

    // compute residual
    res_S = inv_equiv_stress *
            (flow_resistance - time_step_quantities_.last_substep_flow_resistance_[gp_] -
                harden_tang * delta_plastic_strain);

    // check convergence
    if (std::abs(res_S) < conv_tol) return flow_resistance;

    // check whether the maximum number of iterations has been exceeded
    if (iter > max_iter)
    {
      err_status = ErrorType::FailedComputationFlowResistance;
      return -1.0;
    }

    // compute derivatives of the hardening tangent
    derivs_harden_tang =
        compute_derivatives_of_hardening_tangent(equiv_stress, flow_resistance, err_status);
    if (err_status != ErrorType::NoErrors) return -1.0;


    temp = (1.0 - derivs_harden_tang(1) * delta_plastic_strain);
    // compute inverse Jacobian
    inv_J = equiv_stress / temp;
    if (std::abs(inv_J) > 1.0e10)
    {
      err_status = ErrorType::FailedComputationFlowResistance;
      return -1.0;
    }

    // update flow resistance
    flow_resistance += -res_S * inv_J;
  }
}


/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
double Mat::Viscoplastic::Anand::compute_hardening_tangent(
    const double equiv_stress, const double flow_resistance, ErrorType& err_status)
{
  // ensure the error status is clean
  err_status = ErrorType::NoErrors;

  // compute inverse stress ratio
  const double inv_stress_ratio = flow_resistance / equiv_stress;

  // compute saturation level \f$ S_{level} \f$
  const double S_level =
      const_pars_.inv_S_star * std::pow(inv_stress_ratio, const_pars_.eN + 1.0) * equiv_stress;

  // compute hardening tangent depending on the sign of \f$ 1.0 - S_{level} \f$
  if (1.0 - S_level < 0.0)
  {
    return -const_pars_.H_0 * std::pow(S_level - 1.0, const_pars_.a);
  }
  return const_pars_.H_0 * std::pow(1.0 - S_level, const_pars_.a);
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
Core::LinAlg::Matrix<2, 1> Mat::Viscoplastic::Anand::compute_derivatives_of_hardening_tangent(
    const double equiv_stress, const double flow_resistance, ErrorType& err_status)
{
  // declare output variable
  Core::LinAlg::Matrix<2, 1> derivs_of_harden_tang{true};

  // compute inverse stress ratio
  const double inv_stress_ratio = flow_resistance / equiv_stress;

  // compute saturation level \f$ S_{level} \f$
  const double S_level =
      const_pars_.inv_S_star * std::pow(inv_stress_ratio, const_pars_.eN + 1.0) * equiv_stress;

  // compute derivative of the saturation level w.r.t. equivalent stress
  const double dS_level_d_equiv_stress =
      -const_pars_.eN * const_pars_.inv_S_star * std::pow(inv_stress_ratio, const_pars_.eN + 1.0);
  // compute derivative of the saturation level w.r.t. flow resistance
  const double dS_level_d_flow_res =
      (const_pars_.eN + 1.0) * const_pars_.inv_S_star * std::pow(inv_stress_ratio, const_pars_.eN);

  // exponential term
  double temp = 0.0;
  if (1.0 - S_level < 0.0)
  {
    temp = std::pow(S_level - 1.0, const_pars_.a - 1.0);
  }
  else
  {
    temp = std::pow(1.0 - S_level, const_pars_.a - 1.0);
  }


  // compute derivatives of the hardening tangent
  // w.r.t. equiv. stress
  derivs_of_harden_tang(0, 0) = -const_pars_.aH_0 * temp * dS_level_d_equiv_stress;
  // w.r.t. flow resistance
  derivs_of_harden_tang(1, 0) = -const_pars_.aH_0 * temp * dS_level_d_flow_res;


  return derivs_of_harden_tang;
}

Core::LinAlg::Matrix<2, 1> Mat::Viscoplastic::Anand::compute_derivatives_of_flow_resistance(
    const double equiv_stress, const double flow_resistance, const double delta_plastic_strain,
    ErrorType& err_status)
{
  // declare the output derivatives
  Core::LinAlg::Matrix<2, 1> derivs_of_flow_resistance{true};

  // ensure the error status is on NoErrors
  err_status = ErrorType::NoErrors;

  // compute the hardening tangent
  const double harden_tang = compute_hardening_tangent(equiv_stress, flow_resistance, err_status);
  if (err_status != ErrorType::NoErrors) return Core::LinAlg::Matrix<2, 1>{true};


  // compute the derivatives of the hardening tangent w.r.t. equivalent stress and the flow
  // resistance
  Core::LinAlg::Matrix<2, 1> derivs_of_harden_tang =
      compute_derivatives_of_hardening_tangent(equiv_stress, flow_resistance, err_status);
  if (err_status != ErrorType::NoErrors) return Core::LinAlg::Matrix<2, 1>{true};
  double d_harden_tang_d_equiv_stress = derivs_of_harden_tang(0);
  double d_harden_tang_d_flow_res = derivs_of_harden_tang(1);

  // \f$ (1.0 - \frac{\partial T_{\text{harden}}}{\partial S}) (\Delta
  // \varepsilon)^{\text{p}} \f$
  double temp = 1.0 - d_harden_tang_d_flow_res * delta_plastic_strain;
  // \f$ 1.0  / (1.0 - \frac{\partial T_{\text{harden}}}{\partial S}) (\Delta
  // \varepsilon)^{\text{p}} \f$
  double inv_temp = 1.0 / temp;

  // build the derivatives of the flow resistance
  derivs_of_flow_resistance(0, 0) = inv_temp * d_harden_tang_d_equiv_stress * delta_plastic_strain;
  derivs_of_flow_resistance(1, 0) = inv_temp * harden_tang;

  return derivs_of_flow_resistance;
}

std::string Mat::Viscoplastic::Anand::debug_get_error_info(const int gp)
{
  // auxiliaries
  std::ostringstream temp_ostream{};
  // set output format for the numbers -> we can set it here for the
  // entire error message
  std::cout << std::fixed << std::setprecision(16);
  temp_ostream << std::fixed << std::setprecision(16);

  // declare error message for output
  std::string extended_error_message{""};

  extended_error_message += "---> VISCOPLASTIC LAW: Anand \n";
  extended_error_message += "last_flow_resistance: \n";
  extended_error_message += "Double<1,1> \n";
  temp_ostream << time_step_quantities_.last_flow_resistance_[gp_] << std::endl;
  extended_error_message += temp_ostream.str();
  temp_ostream.str("");
  extended_error_message += "last_plastic_strain: \n";
  extended_error_message += "Double<1,1> \n";
  temp_ostream << time_step_quantities_.last_plastic_strain_[gp_] << std::endl;
  extended_error_message += temp_ostream.str();
  temp_ostream.str("");

  return extended_error_message;
}

void Mat::Viscoplastic::Anand::debug_set_last_values(
    const int gp, const double last_flow_resistance, const double last_plastic_strain)
{
  time_step_quantities_.last_flow_resistance_[gp] = last_flow_resistance;
  time_step_quantities_.last_substep_flow_resistance_[gp] = last_flow_resistance;
  time_step_quantities_.last_plastic_strain_[gp] = last_plastic_strain;
  time_step_quantities_.last_substep_plastic_strain_[gp] = last_plastic_strain;
}

void Mat::Viscoplastic::Anand::register_output_data_names(
    std::unordered_map<std::string, int>& names_and_size) const
{
  names_and_size["flow_resistance"] = 1;
}

bool Mat::Viscoplastic::Anand::evaluate_output_data(
    const std::string& name, Core::LinAlg::SerialDenseMatrix& data) const
{
  if (name == "flow_resistance")
  {
    for (int gp = 0; gp < static_cast<int>(time_step_quantities_.current_flow_resistance_.size());
        ++gp)
    {
      data(gp, 0) = time_step_quantities_.current_flow_resistance_[gp];
    }
    return true;
  }

  return false;
}

FOUR_C_NAMESPACE_CLOSE
