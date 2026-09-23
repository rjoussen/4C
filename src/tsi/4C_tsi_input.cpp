// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_tsi_input.hpp"

#include "4C_contact_input.hpp"
#include "4C_io_input_spec_builders.hpp"

FOUR_C_NAMESPACE_OPEN

std::vector<Core::IO::InputSpec> TSI::valid_parameters()
{
  using namespace Core::IO::InputSpecBuilders;

  std::vector<Core::IO::InputSpec> specs;
  specs.push_back(group("TSI DYNAMIC",
      {

          // coupling strategy for (partitioned and monolithic) TSI solvers
          deprecated_selection<SolutionSchemeOverFields>("COUPALGO",
              {
                  {"tsi_oneway", SolutionSchemeOverFields::OneWay},
                  {"tsi_sequstagg", SolutionSchemeOverFields::SequStagg},
                  {"tsi_iterstagg", SolutionSchemeOverFields::IterStagg},
                  {"tsi_iterstagg_aitken", SolutionSchemeOverFields::IterStaggAitken},
                  {"tsi_iterstagg_aitkenirons", SolutionSchemeOverFields::IterStaggAitkenIrons},
                  {"tsi_iterstagg_fixedrelax", SolutionSchemeOverFields::IterStaggFixedRel},
                  {"tsi_monolithic", SolutionSchemeOverFields::Monolithic},
              },
              {.description = "Coupling strategies for TSI solvers",
                  .default_value = SolutionSchemeOverFields::Monolithic}),


          parameter<bool>(
              "MATCHINGGRID", {.description = "is matching grid", .default_value = true}),

          // output type
          parameter<int>(
              "RESTARTEVERY", {.description = "write restart possibility every RESTARTEVERY steps",
                                  .default_value = 1}),

          // time loop control
          parameter<int>(
              "NUMSTEP", {.description = "maximum number of Timesteps", .default_value = 200}),
          parameter<double>(
              "MAXTIME", {.description = "total simulation time", .default_value = 1000.0}),

          parameter<double>(
              "TIMESTEP", {.description = "time step size dt", .default_value = 0.05}),
          parameter<int>("ITEMAX",
              {.description = "maximum number of iterations over fields", .default_value = 10}),
          parameter<int>("ITEMIN",
              {.description = "minimal number of iterations over fields", .default_value = 1}),
          parameter<int>("RESULTSEVERY",
              {.description = "increment for writing solution", .default_value = 1}),

          parameter<ConvNorm>("NORM_INC",
              {.description = "type of norm for convergence check of primary variables in TSI",
                  .default_value = ConvNorm::Abs})},
      {.required = false}));

  /*----------------------------------------------------------------------*/
  /* parameters for monolithic TSI */
  specs.push_back(group("TSI DYNAMIC/MONOLITHIC",
      {

          // convergence tolerance of tsi residual
          parameter<double>("CONVTOL",
              {.description = "tolerance for convergence check of TSI", .default_value = 1e-6}),
          // Iterationparameters
          parameter<double>("TOLINC",
              {.description = "tolerance for convergence check of TSI-increment in monolithic TSI",
                  .default_value = 1.0e-6}),

          parameter<ConvNorm>(
              "NORM_RESF", {.description = "type of norm for residual convergence check",
                               .default_value = ConvNorm::Abs}),

          deprecated_selection<BinaryOp>("NORMCOMBI_RESFINC",
              {
                  {"And", BinaryOp::bop_and},
                  {"Or", BinaryOp::bop_or},
                  {"Coupl_Or_Single", BinaryOp::bop_coupl_or_single},
                  {"Coupl_And_Single", BinaryOp::bop_coupl_and_single},
                  {"And_Single", BinaryOp::bop_and_single},
                  {"Or_Single", BinaryOp::bop_or_single},
              },
              {.description =
                      "binary operator to combine primary variables and residual force values",
                  .default_value = BinaryOp::bop_coupl_and_single}),

          parameter<VectorNorm>(
              "ITERNORM", {.description = "type of norm to be applied to residuals",
                              .default_value = VectorNorm::Rms}),

          parameter<NlnSolTech>("NLNSOL", {.description = "Nonlinear solution technique",
                                              .default_value = NlnSolTech::fullnewton}),


          parameter<double>(
              "PTCDT", {.description = "pseudo time step for pseudo-transient "
                                       "continuation (PTC) stabilised Newton procedure",
                           .default_value = 0.1}),

          // number of linear solver used for monolithic TSI
          parameter<int>("LINEAR_SOLVER",
              {.description = "number of linear solver used for monolithic TSI problems",
                  .default_value = -1}),

          // convergence criteria adaptivity of monolithic TSI solver
          parameter<bool>("ADAPTCONV", {.description = "Switch on adaptive control of linear "
                                                       "solver tolerance for nonlinear solution",
                                           .default_value = false}),
          parameter<double>("ADAPTCONV_BETTER",
              {.description =
                      "The linear solver shall be this much better than the current nonlinear "
                      "residual in the nonlinear convergence limit",
                  .default_value = 0.1}),

          parameter<bool>("INFNORMSCALING",
              {.description = "Scale blocks of matrix with row infnorm?", .default_value = true}),

          // merge TSI block matrix to enable use of direct solver in monolithic TSI
          // default: "No", i.e. use block matrix
          parameter<bool>("MERGE_TSI_BLOCK_MATRIX",
              {.description = "Merge TSI block matrix", .default_value = false}),

          deprecated_selection<LineSearch>("TSI_LINE_SEARCH",
              {
                  {"none", LineSearch::LS_none},
                  {"structure", LineSearch::LS_structure},
                  {"thermo", LineSearch::LS_thermo},
                  {"and", LineSearch::LS_and},
                  {"or", LineSearch::LS_or},
              },
              {.description = "line-search strategy", .default_value = LineSearch::LS_none})},
      {.required = false}));

  /*----------------------------------------------------------------------*/
  /* parameters for partitioned TSI */
  specs.push_back(group("TSI DYNAMIC/PARTITIONED",
      {

          parameter<CouplingVariable>(
              "COUPVARIABLE", {.description = "Coupling variable",
                                  .default_value = CouplingVariable::Displacement}),


          // Solver parameter for relaxation of iterative staggered partitioned TSI
          parameter<double>("MAXOMEGA",
              {.description =
                      "largest omega allowed for Aitken relaxation (0.0 means no constraint)",
                  .default_value = 0.0}),
          parameter<double>(
              "FIXEDOMEGA", {.description = "fixed relaxation parameter", .default_value = 1.0}),

          // convergence tolerance of outer iteration loop
          parameter<double>("CONVTOL",
              {.description =
                      "tolerance for convergence check of outer iteraiton within partitioned TSI",
                  .default_value = 1e-6}),
      },
      {.required = false}));

  /*----------------------------------------------------------------------*/
  /**
   * Constitutive model for condensed Lagrange-multiplier TSI contact.
   *
   * Notation follows Seitz (2019), Computational Methods for Thermo-Elasto-Plastic Contact,
   * Section 2.7, especially Eqs. (2.125)--(2.130). Side (1) is the source ("slave" in the
   * input), side (2) the target ("master"), and n is the outward source normal.
   * Compressive contact pressure is denoted by p_n <= 0.
   * The heat-transfer inputs are the pressure-independent constants in Eq. (2.127):
   * \f[
   *   \bar\gamma^{(1)}=\text{HEATTRANSSLAVE},\qquad
   *   \bar\gamma^{(2)}=\text{HEATTRANSMASTER},\qquad
   *   \gamma^{(i)}=|p_n|\bar\gamma^{(i)} .
   * \f]
   * They define the contact heat conductivity and dissipation split ratio of Eq. (2.130):
   * \f[
   *   \beta_c=\frac{\bar\gamma^{(1)}\bar\gamma^{(2)}}
   *                 {\bar\gamma^{(1)}+\bar\gamma^{(2)}},\qquad
   *   \delta_c=\frac{\bar\gamma^{(1)}}{\bar\gamma^{(1)}+\bar\gamma^{(2)}} .
   * \f]
   * Both constants must be nonnegative and their sum positive. The effective conductance
   * is beta_c |p_n|; equal constants split frictional heat equally between the bodies.
   *
   * With the projection chi_t onto the target surface, the temperature jump and outward
   * contact heat fluxes in Eqs. (2.128)--(2.129) are
   * \f[
   *   [\![T]\!]=T^{(1)}-(T^{(2)}\circ\chi_t),\qquad
   *   q_c^{(1)}=\beta_c|p_n|[\![T]\!]-\delta_c\,t_\tau\cdot v_\tau,\qquad
   *   q_c^{(2)}=-\beta_c|p_n|[\![T]\!]-(1-\delta_c)t_\tau\cdot v_\tau,\qquad
   *   q_c^{(1)}+q_c^{(2)}=-t_\tau\cdot v_\tau .
   * \f]
   * Thus negative q_c supplies heat to a body.
   * The frictional power t_tau . v_tau is only the mechanical part of the total contact
   * dissipation D_c in Eq. (2.120), which also includes heat-transfer terms.
   *
   * Coulomb friction uses mu_0 = FrCoeffOrBound from the contact condition, T_0 = TEMP_REF,
   * and T_d = TEMP_DAMAGE. Eqs. (2.124)--(2.125) give
   * \f[
   *   \vartheta_c=\max(T^{(1)},T^{(2)}\circ\chi_t),\qquad
   *   \mu(\vartheta_c)=\mu_0\frac{(\vartheta_c-T_d)^2}{(T_d-T_0)^2},\qquad
   *   \|t_\tau\|\leq\mu(\vartheta_c)|p_n| .
   * \f]
   * TEMP_DAMAGE must exceed TEMP_REF. The coefficient equals mu_0 at TEMP_REF and
   * vanishes at TEMP_DAMAGE. This is an unclamped quadratic: it increases again above
   * TEMP_DAMAGE and exceeds mu_0 below TEMP_REF. The default large damage temperature
   * makes friction approximately temperature independent at ordinary temperatures.
   *
   * The Nitsche parameters below are legacy registrations; the condensed multiplier
   * implementation does not use them. PENALTYPARAM_THERMO does not enable penalty TSI contact.
   */
  specs.push_back(group("TSI CONTACT",
      {parameter<double>("HEATTRANSSLAVE",
           {.description = "Source-side heat-transfer constant gamma_bar^(1) (Seitz, Eq. 2.127). "
                           "Conductance is beta_c*abs(p_n), with beta_c = "
                           "gamma_bar^(1)*gamma_bar^(2)/(gamma_bar^(1)+gamma_bar^(2)). "
                           "The source receives delta_c = "
                           "gamma_bar^(1)/(gamma_bar^(1)+gamma_bar^(2)) of the frictional heat. "
                           "Both coefficients must be nonnegative with a positive sum.",
               .default_value = 0.0}),
          parameter<double>("HEATTRANSMASTER",
              {.description =
                      "Target-side heat-transfer constant gamma_bar^(2) (Seitz, Eq. 2.127). "
                      "Conductance is beta_c*abs(p_n), with beta_c = "
                      "gamma_bar^(1)*gamma_bar^(2)/(gamma_bar^(1)+gamma_bar^(2)). "
                      "The target receives 1-delta_c = "
                      "gamma_bar^(2)/(gamma_bar^(1)+gamma_bar^(2)) of the frictional heat. "
                      "Both coefficients must be nonnegative with a positive sum.",
                  .default_value = 0.0}),
          parameter<double>("TEMP_DAMAGE",
              {.description =
                      "T_d in mu(vartheta_c) = mu_0*(vartheta_c-T_d)^2/(T_d-T_0)^2 "
                      "(Seitz, Eq. 2.125); vartheta_c is the maximum contact-side temperature. "
                      "Must exceed TEMP_REF. Friction vanishes here but is not clamped: "
                      "it increases again above this temperature.",
                  .default_value = 1.0e12}),

          parameter<double>(
              "TEMP_REF", {.description = "T_0 in the quadratic temperature-dependent friction "
                                          "law (Seitz, Eq. 2.125). When vartheta_c equals T_0, "
                                          "mu equals mu_0 = "
                                          "FrCoeffOrBound from the Coulomb contact condition.",
                              .default_value = 0.0}),

          parameter<bool>("CONDENSED_LM_INCREMENTS",
              {.description =
                      "Condense and recover mechanical and thermal contact Lagrange multiplier "
                      "increments. If false, retain the legacy absolute-multiplier condensation.",
                  .default_value = false}),

          parameter<double>("NITSCHE_THETA_TSI",
              {.description = "Legacy Nitsche option: +1 symmetric, 0 non-symmetric, "
                              "-1 skew-symmetric. Unused by condensed TSI contact.",
                  .default_value = 0.0}),

          parameter<CONTACT::NitscheWeighting>("NITSCHE_WEIGHTING_TSI",
              {.description = "Legacy weighting of Nitsche consistency terms. "
                              "Unused by condensed TSI contact.",
                  .default_value = CONTACT::NitscheWeighting::harmonic}),

          parameter<bool>("NITSCHE_PENALTY_ADAPTIVE_TSI",
              {.description = "Legacy Nitsche penalty adaptation after each converged time step. "
                              "Unused by condensed TSI contact.",
                  .default_value = true}),

          parameter<double>("PENALTYPARAM_THERMO",
              {.description = "Legacy thermal Nitsche penalty parameter. Unused by condensed "
                              "TSI contact; setting it does not enable penalty TSI contact.",
                  .default_value = 0.0})},
      {.required = false}));
  return specs;
}

FOUR_C_NAMESPACE_CLOSE
