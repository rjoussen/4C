// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_FLUID_ELE_PARAMETER_HPP
#define FOUR_C_FLUID_ELE_PARAMETER_HPP

#include "4C_config.hpp"

#include "4C_fluid_ele_parameter_timint.hpp"
#include "4C_fluid_turbulence_input.hpp"
#include "4C_inpar_fluid.hpp"
#include "4C_utils_exceptions.hpp"

#include <Teuchos_StandardParameterEntryValidators.hpp>

FOUR_C_NAMESPACE_OPEN



namespace Discret
{
  namespace Elements
  {
    class FluidEleParameter
    {
     public:
      virtual ~FluidEleParameter() = default;

      /*========================================================================*/
      //! @name set-routines
      /*========================================================================*/

      //! general fluid parameter are set
      void set_element_general_fluid_parameter(Teuchos::ParameterList& params,  //> parameter list
          int myrank);                                                          //> proc id

      //! turbulence parameters are set
      void set_element_turbulence_parameters(Teuchos::ParameterList& params);  //> parameter list

      /// set loma parameters
      void set_element_loma_parameter(Teuchos::ParameterList& params);  //> parameter list

      //! set two-phase parameters
      void set_element_two_phase_parameter(Teuchos::ParameterList& params);  //> parameter list

      /*========================================================================*/
      //! @name access-routines
      /*========================================================================*/

      /*----------------------------------------------------*/
      //! @name general parameters
      /*----------------------------------------------------*/

      //! Flag for physical type of the fluid flow (incompressible, loma, varying_density,
      //! Boussinesq, poro)
      Inpar::FLUID::PhysicalType physical_type() const { return physicaltype_; };
      //! flag to (de)activate conservative formulation
      bool is_conservative() const { return is_conservative_; };
      //! flag to (de)activate Newton linearization
      bool is_newton() const { return is_newton_; };
      //! flag to (de)activate second derivatives
      bool is_inconsistent() const { return is_inconsistent_; };
      //! flag to (de)activate potential reactive terms
      bool reaction() const { return reaction_; };
      //! Return function number of Oseen advective field
      int oseen_field_func_no() const { return oseenfieldfuncno_; };
      //! flag to activate consistent reconstruction of second derivatives
      bool is_reconstruct_der() const { return is_reconstructder_; };

      /*----------------------------------------------------*/
      //! @name stabilization parameters
      /*----------------------------------------------------*/

      //! get the stabtype
      Inpar::FLUID::StabType stab_type() const { return stabtype_; };
      /// parameter for residual stabilization
      //! Flag to (de)activate time-dependent subgrid stabilization
      Inpar::FLUID::SubscalesTD tds() const { return tds_; };
      //! Flag to (de)activate time-dependent term in large-scale momentum equation
      Inpar::FLUID::Transient transient() const { return transient_; };
      //! Flag to (de)activate PSPG stabilization
      bool pspg() const { return pspg_; };
      //! Flag to (de)activate SUPG stabilization
      bool supg() const { return supg_; };
      //! Flag to (de)activate viscous term in residual-based stabilization
      Inpar::FLUID::VStab v_stab() const { return vstab_; };
      //! Flag to (de)activate reactive term in residual-based stabilization
      Inpar::FLUID::RStab r_stab() const { return rstab_; };
      //! Flag to (de)activate least-squares stabilization of continuity equation
      bool c_stab() const { return graddiv_; };
      //! Flag to (de)activate cross-stress term -> residual-based VMM
      Inpar::FLUID::CrossStress cross() const { return cross_; };
      //! Flag to (de)activate Reynolds-stress term -> residual-based VMM
      Inpar::FLUID::ReynoldsStress reynolds() const { return reynolds_; };
      //! Flag to define tau
      Inpar::FLUID::TauType which_tau() const { return whichtau_; };
      //! Flag to define characteristic element length for tau_Mu
      Inpar::FLUID::CharEleLengthU char_ele_length_u() const { return charelelengthu_; };
      //! Flag to define characteristic element length for tau_Mp and tau_C
      Inpar::FLUID::CharEleLengthPC char_ele_length_pc() const { return charelelengthpc_; };
      //! (sign) factor for viscous and reactive stabilization terms
      double visc_rea_stab_fac() const { return viscreastabfac_; };

      //! Flag to (de)activate polynomial pressure projection stabilization
      bool ppp() const { return ppp_; };

      //! flag for material evaluation at Gaussian integration points
      bool mat_gp() const { return mat_gp_; };
      //! flag for stabilization parameter evaluation at Gaussian integration points
      bool tau_gp() const { return tau_gp_; };

      /*----------------------------------------------------*/
      //! @name two phase parameters
      /*----------------------------------------------------*/
      double get_interface_thickness() const { return interface_thickness_; };
      bool get_enhanced_gauss_rule_in_interface() const { return enhanced_gaussrule_; }
      bool get_include_surface_tension() const { return include_surface_tension_; };


      /*----------------------------------------------------*/
      //! @name turbulence model
      /*----------------------------------------------------*/

      /// constant parameters for the turbulence formulation
      /// subgrid viscosity models
      //! flag to define turbulence model
      Inpar::FLUID::TurbModelAction turb_mod_action() const { return turb_mod_action_; };
      double cs() const { return Cs_; };
      bool cs_averaged() const { return Cs_averaged_; };
      double ci() const { return Ci_; };
      void setvan_driestdamping(double damping)
      {
        van_Driest_damping_ = damping;
        return;
      };
      double van_driestdamping() const { return van_Driest_damping_; };
      bool include_ci() const { return include_Ci_; };
      double ltau() const { return l_tau_; };
      //! flag to (de)activate fine-scale subgrid viscosity
      Inpar::FLUID::FineSubgridVisc fssgv() const { return fssgv_; };
      // Flag to Vreman filter method
      Inpar::FLUID::VremanFiMethod vrfi() const { return vrfi_; };
      /// multifractal subgrid-scales
      double csgs() const { return Csgs_; };
      double csgs_phi() const
      {
        double tmp = 0.0;
        if (not adapt_Csgs_phi_)
          tmp = Csgs_phi_;
        else
          tmp = Csgs_phi_ * meanCai_;
        return tmp;
      };
      double alpha() const { return alpha_; };
      bool calc_n() const { return CalcN_; };
      double n() const { return N_; };
      FLUID::RefVelocity ref_vel() const { return refvel_; };
      FLUID::RefLength ref_length() const { return reflength_; };
      double c_nu() const { return c_nu_; };
      double c_diff() const { return c_diff_; };
      bool near_wall_limit() const { return near_wall_limit_; };
      bool near_wall_limit_scatra() const { return near_wall_limit_scatra_; };
      bool b_gp() const { return B_gp_; };
      double beta() const { return beta_; };
      double mfs_is_conservative() const { return mfs_is_conservative_; };
      double adapt_csgs_phi() const { return adapt_Csgs_phi_; };
      void set_csgs_phi(double meanCai)
      {
        meanCai_ = meanCai;
        return;
      };
      bool consistent_mfs_residual() const { return consistent_mfs_residual_; };

      /*----------------------------------------------------*/
      //! @name loma parameters
      /*----------------------------------------------------*/

      //! flag for material update
      virtual bool update_mat() const { return update_mat_; };
      //! flag to (de)activate continuity SUPG term
      virtual bool conti_supg() const { return conti_supg_; };
      //! flag to (de)activate continuity cross-stress term -> residual-based VMM
      virtual Inpar::FLUID::CrossStress conti_cross() const { return conti_cross_; };
      //! flag to (de)activate continuity Reynolds-stress term -> residual-based VMM
      virtual Inpar::FLUID::ReynoldsStress conti_reynolds() const { return conti_reynolds_; };
      //! flag to (de)activate cross- and Reynolds-stress terms in loma continuity equation
      virtual bool multi_frac_loma_conti() const { return multifrac_loma_conti_; };

     protected:
      /*----------------------------------------------------*/
      //! @name general parameters
      /*----------------------------------------------------*/

      //! Flag SetGeneralParameter was called
      bool set_general_fluid_parameter_;

      //! Flag for physical type of the fluid flow (incompressible, loma, varying_density,
      //! Boussinesq, Poro)
      Inpar::FLUID::PhysicalType physicaltype_;
      //! parameter to switch the stabilization
      Inpar::FLUID::StabType stabtype_;
      /// Flags to switch on/off the different fluid formulations
      //! flag to (de)activate conservative formulation
      bool is_conservative_;
      //! flag to (de)activate Newton linearization
      bool is_newton_;
      //! flag to (de)activate second derivatives
      bool is_inconsistent_;
      //! flag to (de)activate potential reactive terms
      bool reaction_;
      //! function number of Oseen advective field
      int oseenfieldfuncno_;
      //! flag to activate consistent reconstruction of second derivatives
      bool is_reconstructder_;

      /*----------------------------------------------------*/
      //! @name stabilization parameters
      /*----------------------------------------------------*/

      /// parameter for residual based stabilizations
      //! Flag to (de)activate time-dependent subgrid stabilization
      Inpar::FLUID::SubscalesTD tds_;
      //! Flag to (de)activate time-dependent term in large-scale momentum equation
      Inpar::FLUID::Transient transient_;
      //! Flag to (de)activate PSPG stabilization
      bool pspg_;
      //! Flag to (de)activate SUPG stabilization
      bool supg_;
      //! Flag to (de)activate viscous term in residual-based stabilization
      Inpar::FLUID::VStab vstab_;
      //! Flag to (de)activate reactive term in residual-based stabilization
      Inpar::FLUID::RStab rstab_;
      //! Flag to (de)activate least-squares stabilization of continuity equation
      bool graddiv_;
      //! Flag to (de)activate cross-stress term -> residual-based VMM
      Inpar::FLUID::CrossStress cross_;
      //! Flag to (de)activate Reynolds-stress term -> residual-based VMM
      Inpar::FLUID::ReynoldsStress reynolds_;
      //! Flag to define tau
      Inpar::FLUID::TauType whichtau_;
      //! Flag to define characteristic element length for tau_Mu
      Inpar::FLUID::CharEleLengthU charelelengthu_;
      //! Flag to define characteristic element length for tau_Mp and tau_C
      Inpar::FLUID::CharEleLengthPC charelelengthpc_;
      //! (sign) factor for viscous and reactive stabilization terms
      double viscreastabfac_;

      //! Flag to (de)activate PPP (polynomial pressure projection) stabilization
      bool ppp_;

      /// parameter for evaluation of material and stabilization parameter
      //! flag for material evaluation at Gaussian integration points
      bool mat_gp_;
      //! flag for stabilization parameter evaluation at Gaussian integration points
      bool tau_gp_;

      /*----------------------------------------------------*/
      //! @name two phase parameters
      /*----------------------------------------------------*/

      double interface_thickness_;
      bool enhanced_gaussrule_;
      bool include_surface_tension_;

      /*----------------------------------------------------*/
      //! @name turbulence model
      /*----------------------------------------------------*/

      /// constant parameters for the turbulence formulation
      /// subgrid viscosity models
      //! flag to define turbulence model
      Inpar::FLUID::TurbModelAction turb_mod_action_;
      //! smagorinsky constant
      double Cs_;
      bool Cs_averaged_;
      //! loma: constant for isotropic part of subgrid-stress tensor (according to Yoshizawa 1986)
      double Ci_;
      bool include_Ci_;
      //! van Driest damping for channel flow
      double van_Driest_damping_;
      //! channel length to normalize the normal wall distance
      double l_tau_;
      //! flag to (de)activate fine-scale subgrid viscosity
      Inpar::FLUID::FineSubgridVisc fssgv_;
      // flag to Vreman filter method
      Inpar::FLUID::VremanFiMethod vrfi_;
      /// multifractal subgrid-scales
      double Csgs_;
      double Csgs_phi_;
      double alpha_;
      bool CalcN_;
      double N_;
      FLUID::RefVelocity refvel_;
      FLUID::RefLength reflength_;
      double c_nu_;
      double c_diff_;
      bool near_wall_limit_;
      bool near_wall_limit_scatra_;
      bool B_gp_;
      double beta_;
      bool mfs_is_conservative_;
      bool adapt_Csgs_phi_;
      double meanCai_;
      bool consistent_mfs_residual_;

      /*----------------------------------------------------*/
      //! @name loma parameters
      /*----------------------------------------------------*/

      //! flag for material update
      bool update_mat_;
      //! flag to (de)activate continuity SUPG term
      bool conti_supg_;
      //! flag to (de)activate continuity cross-stress term -> residual-based VMM
      Inpar::FLUID::CrossStress conti_cross_;
      //! flag to (de)activate continuity Reynolds-stress term -> residual-based VMM
      Inpar::FLUID::ReynoldsStress conti_reynolds_;
      //! flag to (de)activate cross- and Reynolds-stress terms in loma continuity equation
      bool multifrac_loma_conti_;

      /*-----------------------------------------------------*/

      /// private Constructor since we are a Singleton.
      FluidEleParameter();

     private:
      //! access time-integration parameters
      Discret::Elements::FluidEleParameterTimInt* fldparatimint_;
    };

  }  // namespace Elements
}  // namespace Discret

FOUR_C_NAMESPACE_CLOSE

#endif
