// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_REDUCED_LUNG_1D_PIPE_FLOW_STRUCTURED_TREE_HPP
#define FOUR_C_REDUCED_LUNG_1D_PIPE_FLOW_STRUCTURED_TREE_HPP

#include "4C_config.hpp"

#include "4C_reduced_lung_1d_pipe_flow_input.hpp"

#include <complex>
#include <functional>
#include <vector>

FOUR_C_NAMESPACE_OPEN
namespace ReducedLung1DPipe
{
  namespace TerminalUnit
  {
    /**
     * @brief Structured tree outflow condition following Olufsen et al.
     * systemic circulation: 10.1152/ajpheart.1999.276.1.H257,
     * pulmonary circulation: 10.1007/s10237-014-0563-y
     *
     * Models the downstream vasculature as a self-similar bifurcating tree.
     * The root impedance Z(omega) is computed using Womersley theory and converted
     * to an impulse response z(t) via IFFT. The boundary pressure is obtained
     * using z(t) and the flow history Q(t): p_ST(t) = integral_0^t z(tau) Q(t-tau) dtau
     */
    struct StructuredTreeTerminalUnit
    {
      double root_radius;                   ///< Root radius of the tree at the terminal node
      double exponent_xi;                   ///< Power-law exponent xi in r_p^xi = r_d1^xi + r_d2^xi
      double asymmetry_ratio_gamma;         ///< Asymmetry ratio gamma = (r_small / r_large)^2
      double termination_radius_r_min;      ///< Minimum radius where the tree terminates
      double peripheral_resistance_R_peri;  ///< Peripheral resistance at terminal arterioles
      double length_coefficient;            ///< Coefficient k_l in l = k_l * r^k_exp
      double length_exponent;               ///< Exponent k_exp in l = k_l * r^k_exp

      // ---- Pre-computed at init time ----
      std::vector<std::complex<double>> impedance_spectrum;  ///< Root impedance Z[k] = Z(k*domega)
      std::vector<double> impulse_response;  ///< Impulse response z[k] = z(k*Deltat)
      int n_cycle;                           ///< Steps per cardiac period (= T / dt)
      double dt;                             ///< Time step of the spectrum grid

      // ---- Runtime state ----
      std::vector<double> flow_history;  ///< Circular buffer of past flow values Q
      int flow_history_write_index;      ///< Next write index in the circular buffer
      int flow_history_count;            ///< Number of valid entries in flow_history

      /// Convolution of the impulse response with the flow history, i.e.
      /// sum_{k=1}^{n_history} impulse_response[k] * flow_history[write_index - k].
      double convolved_pressure;
    };

    // Functions for structured tree outlet condition built similar to terminal_unit
    /**
     * Initializes the structured tree model for each terminal unit.
     * @param st_model struct with information from input file
     * @param global_id of terminal node
     * @param root_radius of terminal node
     * @param density_rho fluid parameter
     * @param viscosity_mu fluid parameter
     * @param dt timestep-size
     * @param period cardiac period or final time if no cardiac period defined
     * @return filled structured tree container
     */
    StructuredTreeTerminalUnit create_structured_tree_model(
        const ReducedLung1dPipeFlow::Parameters::TerminalUnits::StructuredTreeModel& st_model,
        int global_id, double root_radius, double density_rho, double viscosity_mu, double dt,
        double period);

    /**
     * Compute the root impedance spectrum Z[k] = Z(k*domega) of a structured tree.
     *
     * The spectrum is stored on the model as StructuredTreeTerminalUnit::impedance_spectrum, which
     * is used by 1D BC to assemble the impulse response.
     * @param st_model structured-tree model
     * @param density_rho fluid parameter
     * @param viscosity_mu fluid parameter
     * @param evaluate_eh_over_r wall stiffness function Eh/r(r)
     * @return root impedance spectrum, Hermitian, of size n_cycle
     */
    std::vector<std::complex<double>> compute_structured_tree_impedance(
        const StructuredTreeTerminalUnit& st_model, double density_rho, double viscosity_mu,
        const std::function<double(double)>& evaluate_eh_over_r);

    /**
     * Evaluate residual when terminal unit is structured-tree model.
     * @param area_A A^n at terminal node coupling to terminal unit
     * @param flow_Q Q^n at terminal node coupling to terminal unit
     * @param reference_area_A0 A0 at terminal node coupling to terminal unit
     * @param beta vessel wall stiffness coefficient
     * @param Pext external pressure given in input file
     * @param st_model updated structured-tree struct
     * @return residual
     */
    double evaluate_structured_tree_residual(double area_A, double flow_Q, double reference_area_A0,
        double beta, double Pext, const StructuredTreeTerminalUnit& st_model);

    /**
     * Evaluate jacobian when terminal unit is structured-tree model.
     * @param area_A A^n at terminal node coupling to terminal unit
     * @param dQdA dQdA^n at terminal node coupling to terminal unit
     * @param beta vessel wall stiffness coefficient
     * @param st_model updated structured-tree struct
     * @return jacobian
     */
    double evaluate_structured_tree_jacobian(
        double area_A, double dQdA, double beta, const StructuredTreeTerminalUnit& st_model);

    /**
     * Update structured tree container: Push new flow_Q into circular buffer containing flow
     * history, then sum z[k].
     * @param st_model structured-tree struct from last time step
     * @param flow_Q new flow to update structured-tree struct with
     */
    void update_structured_tree_data(StructuredTreeTerminalUnit& st_model, double flow_Q);

  }  // namespace TerminalUnit
}  // namespace ReducedLung1DPipe
FOUR_C_NAMESPACE_CLOSE
#endif  // FOUR_C_REDUCED_LUNG_1D_PIPE_FLOW_STRUCTURED_TREE_HPP
