// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_SCATRA_ELE_CALC_POROFLUID_PRESSURE_BASED_HPP
#define FOUR_C_SCATRA_ELE_CALC_POROFLUID_PRESSURE_BASED_HPP


#include "4C_config.hpp"

#include "4C_fem_discretization.hpp"
#include "4C_global_data.hpp"
#include "4C_mat_fluidporo_multiphase.hpp"
#include "4C_mat_fluidporo_multiphase_reactions.hpp"
#include "4C_mat_scatra_multiporo.hpp"
#include "4C_porofluid_pressure_based_ele_action.hpp"
#include "4C_porofluid_pressure_based_ele_calc.hpp"
#include "4C_porofluid_pressure_based_ele_calc_utils.hpp"
#include "4C_porofluid_pressure_based_ele_evaluator.hpp"
#include "4C_porofluid_pressure_based_ele_parameter.hpp"
#include "4C_porofluid_pressure_based_ele_phasemanager.hpp"
#include "4C_porofluid_pressure_based_ele_variablemanager.hpp"
#include "4C_scatra_ele_calc_poro_reac.hpp"
#include "4C_utils_enum.hpp"
#include "4C_utils_function.hpp"

FOUR_C_NAMESPACE_OPEN

// forward declarations
namespace Mat
{
  class ScatraMat;
}

namespace Discret::Elements
{
  template <int nsd, int nen>
  class ScaTraEleInternalVariableManagerPorofluidPressureBased;

  template <Core::FE::CellType discretization_type>
  class ScaTraEleCalcPorofluidPressureBased final
      : public ScaTraEleCalcPoroReac<discretization_type>
  {
    /// private constructor, since we are a Singleton.
    ScaTraEleCalcPorofluidPressureBased(
        int number_of_dofs_per_node, int number_of_scalars, const std::string& discretization_name);

    using scatra_ele_calc = ScaTraEleCalc<discretization_type>;
    using poro_reaction = ScaTraEleCalcPoroReac<discretization_type>;
    using poro = ScaTraEleCalcPoro<discretization_type>;
    using advanced_reaction = ScaTraEleCalcAdvReac<discretization_type>;
    using scatra_ele_calc::nen_;
    using scatra_ele_calc::nsd_;

   public:
    /// Singleton access method
    static ScaTraEleCalcPorofluidPressureBased* instance(
        int number_of_dofs_per_node, int number_of_scalars, const std::string& discretization_name);

    /// Setup element evaluation
    int setup_calc(
        Core::Elements::Element* element, Core::FE::Discretization& discretization) override;

   protected:
    //! extract element based or nodal values
    void extract_element_and_node_values(Core::Elements::Element* element,
        Teuchos::ParameterList& parameters, Core::FE::Discretization& discretization,
        Core::Elements::LocationArray& location_array) override;

    //! extract element based or nodal values
    //! --> L2-projection case: called within extract_element_and_node_values
    void extract_nodal_flux(const Core::FE::Discretization& discretization,
        Core::Elements::LocationArray& location_array, int number_of_fluid_phases);

    //! set internal variables
    void set_internal_variables_for_mat_and_rhs() override;

    //! evaluate material
    void materials(
        std::shared_ptr<const Core::Mat::Material> material,  //!< pointer to current material
        int scalar_id,                                        //!< id of current scalar
        double& densn,                                        //!< density at t_(n)
        double& densnp,          //!< density at t_(n+1) or t_(n+alpha_F)
        double& densam,          //!< density at t_(n+alpha_M)
        double& viscosity,       //!< fluid viscosity
        int gauss_point_id = -1  //!< id of current gauss point (default = -1)
        ) override;

    //! material for scalar transport in a fluid phase
    void material_scatra_in_fluid(
        std::shared_ptr<const Core::Mat::Material> material,  //!< pointer to current material
        int scalar_id,                                        //!< id of current scalar
        double& densn,                                        //!< density at t_(n)
        double& densnp,          //!< density at t_(n+1) or t_(n+alpha_F)
        double& densam,          //!< density at t_(n+alpha_M)
        int gauss_point_id = -1  //!< id of current gauss point (default = -1)
    );

    //! material for scalar transport in a volume fraction
    void material_scatra_in_volfrac(
        std::shared_ptr<const Core::Mat::Material> material,  //!< pointer to current material
        int scalar_id,                                        //!< id of current scalar
        double& densn,                                        //!< density at t_(n)
        double& densnp,          //!< density at t_(n+1) or t_(n+alpha_F)
        double& densam,          //!< density at t_(n+alpha_M)
        int gauss_point_id = -1  //!< id of current gauss point (default = -1)
    );

    //! material for scalar transport in the solid phase
    void material_scatra_in_solid(
        std::shared_ptr<const Core::Mat::Material> material,  //!< pointer to current material
        int scalar_id,                                        //!< id of current scalar
        double& densn,                                        //!< density at t_(n)
        double& densnp,          //!< density at t_(n+1) or t_(n+alpha_F)
        double& densam,          //!< density at t_(n+alpha_M)
        int gauss_point_id = -1  //!< id of current gauss point (default = -1)
    );

    //! material for scalar transport of temperature (energy balance)
    void material_scatra_as_temperature(
        std::shared_ptr<const Core::Mat::Material> material,  //!< pointer to current material
        int scalar_id,                                        //!< id of current scalar
        double& densn,                                        //!< density at t_(n)
        double& densnp,          //!< density at t_(n+1) or t_(n+alpha_F)
        double& densam,          //!< density at t_(n+alpha_M)
        int gauss_point_id = -1  //!< id of current gauss point (default = -1)
    );

    //! Set advanced reaction terms and derivatives
    void set_advanced_reaction_terms(int scalar_id,                     //!< index of current scalar
        std::shared_ptr<Mat::MatListReactions> material_reaction_list,  //!< material reaction list
        const double* gauss_point_coordinates  //!< current Gauss-point coordinates
        ) override;

    //! compute pore pressure
    double compute_pore_pressure() override;

    //! get internal variable manager for pressure-based porofluid
    std::shared_ptr<ScaTraEleInternalVariableManagerPorofluidPressureBased<nsd_, nen_>>
    variable_manager()
    {
      return std::static_pointer_cast<
          ScaTraEleInternalVariableManagerPorofluidPressureBased<nsd_, nen_>>(
          scatra_ele_calc::scatravarmanager_);
    }

    //! calculation of convective element matrix in convective form
    //! the only difference to the base class version is, that there is no scaling with the
    //! density
    void calc_mat_conv(
        Core::LinAlg::SerialDenseMatrix& element_matrix,  //!< element matrix to be filled
        int scalar_id,                                    //!< index of current scalar
        double timefacfac,  //!< domain-integration factor times time-integration factor
        double densnp,      //!< density at time_(n+1)
        const Core::LinAlg::Matrix<nen_, 1>& sgconv  //!< subgrid-scale convective operator
        ) override;

    //! adaption of convective term for rhs
    //! the only difference to the base class version is, that there is no scaling with the
    //! density
    void recompute_conv_phi_for_rhs(int scalar_id,      //!< index of current scalar
        const Core::LinAlg::Matrix<nsd_, 1>& sgvelint,  //!< subgrid-scale velocity at Gauss point
        double densnp,                                  //!< density at time_(n+1)
        double densn,                                   //!< density at time_(n)
        double vdiv                                     //!< velocity divergence
        ) override;

    //! calculation of convective element matrix: add conservative contributions
    //! the only difference to the base class version is, that there is no scaling with the
    //! density
    void calc_mat_conv_add_cons(
        Core::LinAlg::SerialDenseMatrix& element_matrix,  //!< element matrix to be filled
        int scalar_id,                                    //!< index of current scalar
        double timefacfac,  //!< domain-integration factor times time-integration factor
        double vdiv,        //!< velocity divergence
        double densnp       //!< density at time_(n+1)
        ) override;

    //! calculation of mass element matrix (standard shape functions)
    void calc_mat_mass(
        Core::LinAlg::SerialDenseMatrix& element_matrix,  //!< element matrix to be filled
        const int& scalar_id,                             //!< index of current scalar
        const double& fac,                                //!< domain-integration factor
        const double& densam                              //!< density at time_(n+am)
        ) override;

    //! calculation of convective element matrix (OD term structure coupling)
    //! difference to base class: linearization of mesh motion + shapederivatives pressure
    //! gradient have to be included
    void calc_conv_od_mesh(Core::LinAlg::SerialDenseMatrix& element_matrix, int scalar_id,
        int number_of_dofs_per_node_mesh, double fac, double rhsfac, double densnp, double J,
        const Core::LinAlg::Matrix<nsd_, 1>& gradphi,
        const Core::LinAlg::Matrix<nsd_, 1>& convelint) override;

    //! calculation of linearized mass (off-diagonal/shapederivative term mesh)
    //! difference to base class: linearization of porosity is included
    void calc_lin_mass_od_mesh(
        Core::LinAlg::SerialDenseMatrix& element_matrix,  //!< element matrix to be filled
        int scalar_id,                                    //!< index of current scalar
        int number_of_dofs_per_node_mesh,  //!< number of dofs per node of ale element
        double rhsfac,  //!< time-integration factor for rhs times domain-integration factor
        double fac,     //!< domain-integration factor
        double densam,  //!< density at time_(n+am)
        double densnp,  //!< density at time_(n+1)
        double phinp,   //!< scalar at time_(n+1)
        double hist,    //!< history of time integration
        double J,       //!< determinant of Jacobian det(dx/ds)
        const Core::LinAlg::Matrix<1, nsd_ * nen_>&
            dJ_dmesh  //!< derivative of det(dx/ds) w.r.t. mesh displacement
        ) override;

    //! standard Galerkin transient, old part of rhs and source term (off-diagonal/shapederivative
    //! term mesh) difference to base class: linearization of porosity and advanced reaction terms
    //! are included
    void calc_hist_and_source_od_mesh(
        Core::LinAlg::SerialDenseMatrix& element_matrix,  //!< element current to be filled
        int scalar_id,                                    //!< index of current scalar
        int number_of_dofs_per_node_mesh,  //!< number of dofs per node of ale element
        double fac,                        //!< domain-integration factor
        double rhsint,                     //!< rhs at Gauss point
        double J,                          //!< determinant of Jacobian det(dx/ds)
        const Core::LinAlg::Matrix<1, nsd_ * nen_>&
            dJ_dmesh,  //!< derivative of det(dx/ds) w.r.t. mesh displacement
        double densnp  //!< density
        ) override;

    //! standard Galerkin diffusive term (off-diagonal/shapederivative term mesh)
    //! difference to base class: linearization of porosity and effective diffusivity
    void calc_diff_od_mesh(
        Core::LinAlg::SerialDenseMatrix& element_matrix,  //!< element current to be filled
        int scalar_id,                                    //!< index of current scalar
        int number_of_dofs_per_node_mesh,  //!< number of dofs per node of ale element
        double diffcoeff,                  //!< diffusion coefficient
        double fac,                        //!< domain-integration factor
        double rhsfac,  //!< time-integration factor for rhs times domain-integration factor
        double J,       //!< determinant of Jacobian det(dx/ds)
        const Core::LinAlg::Matrix<nsd_, 1>& gradphi,    //!< scalar gradient at Gauss point
        const Core::LinAlg::Matrix<nsd_, 1>& convelint,  //!< convective velocity
        const Core::LinAlg::Matrix<1, nsd_ * nen_>&
            dJ_dmesh  //!< derivative of det(dx/ds) w.r.t. mesh displacement
        ) override;

    //! reactive terms (standard Galerkin) (off-diagonal/shapederivative term mesh)
    //! difference to base class: linearization of porosity
    void calc_react_od_mesh(
        Core::LinAlg::SerialDenseMatrix& element_matrix,  //!< element current to be filled
        int scalar_id,                                    //!< index of current scalar
        int number_of_dofs_per_node_mesh,  //!< number of dofs per node of ale element
        double rhsfac,   //!< time-integration factor for rhs times domain-integration factor
        double rea_phi,  //!< reactive term
        double J,        //!< determinant of Jacobian det(dx/ds)
        const Core::LinAlg::Matrix<1, nsd_ * nen_>&
            dJ_dmesh  //!< derivative of det(dx/ds) w.r.t. mesh displacement
        ) override;

    //! calculation of convective element matrix in convective form (off-diagonal term fluid)
    void calc_mat_conv_od_fluid(
        Core::LinAlg::SerialDenseMatrix& element_matrix,  //!< element matrix to be filled
        int scalar_id,                                    //!< index of current scalar
        int number_of_dofs_per_node_mesh,  //!< number of dofs per node of fluid element // only a
                                           //!< dummy variable
        double rhsfac,  //!< domain-integration factor times time-integration factor
        double densnp,  //!< density at time_(n+1)
        const Core::LinAlg::Matrix<nsd_, 1>& gradphi  //!< scalar gradient
        ) override;

    //! calculation of convective element matrix in convective form -- additional conservative
    //! contributions (off-diagonal term fluid) not yet implemented --> FOUR_C_THROW
    void calc_mat_conv_add_cons_od_fluid(
        Core::LinAlg::SerialDenseMatrix& element_matrix,  //!< element matrix to be filled
        int scalar_id,                                    //!< index of current scalar
        int number_of_dofs_per_node_mesh,  //!< number of dofs per node of fluid element
        double timefacfac,  //!< domain-integration factor times time-integration factor
        double densnp,      //!< density at time_(n+1)
        double phinp        //!< scalar at time_(n+1)
        ) override;

    //! calculation of linearized mass (off-diagonal porofluid terms)
    //! linearization of porosity * saturation * vtrans
    void calc_lin_mass_od_fluid(
        Core::LinAlg::SerialDenseMatrix& element_matrix,  //!< element matrix to be filled
        int scalar_id,                                    //!< index of current scalar
        int number_of_dofs_per_node_mesh,  //!< number of dofs per node of fluid element // only a
                                           //!< dummy variable
        double rhsfac,  //!< time-integration factor for rhs times domain-integration factor
        double fac,     //!< domain-integration factor
        double densam,  //!< density at time_(n+am)
        double densnp,  //!< density at time_(n+1)
        double phinp,   //!< scalar at time_(n+1)
        double hist     //!< history of time integration
        ) override;

    //! calculate linearization of a mass matrix type matrix (off-diagonal porofluid terms)
    void calc_mass_matrix_type_od_porofluid(Core::LinAlg::SerialDenseMatrix& element_matrix,
        int scalar_id, const std::vector<double>* per_factor_mass_matrix_OD_porofluid,
        int total_number_of_porofluid_dofs, double pre_factor);

    //! standard Galerkin transient, old part of rhs and source term (off-diagonal porofluid terms)
    //! linearization of porosity * saturation * vtrans + advanced reaction terms
    void calc_hist_and_source_od_fluid(
        Core::LinAlg::SerialDenseMatrix& element_matrix,  //!< element current to be filled
        int scalar_id,                                    //!< index of current scalar
        int number_of_dofs_per_node_mesh,  //!< number of dofs per node of ale element
        double fac,                        //!< domain-integration factor
        double rhsint,                     //!< rhs at Gauss point
        double densnp                      //!< density
        ) override;

    //! standard Galerkin reactive term (off-diagonal porofluid terms)
    //! linearization of porosity * saturation * vreact
    void calc_react_od_fluid(
        Core::LinAlg::SerialDenseMatrix& element_matrix,  //!< element current to be filled
        int scalar_id,                                    //!< index of current scalar
        int number_of_dofs_per_node_mesh,  //!< number of dofs per node of ale element
        double rhsfac,  //!< time-integration factor for rhs times domain-integration factor
        double rea_phi  //!< rhs at Gauss point
        ) override;

    //! standard Galerkin diffusive term (off-diagonal porofluid terms)
    //! linearization of porosity*saturation*d_eff
    void calc_diff_od_fluid(
        Core::LinAlg::SerialDenseMatrix& element_matrix,  //!< element current to be filled
        int scalar_id,                                    //!< index of current scalar
        int number_of_dofs_per_node_mesh,  //!< number of dofs per node of ale element
        double rhsfac,  //!< time-integration factor for rhs times domain-integration factor
        const Core::LinAlg::Matrix<nsd_, 1>& gradphi  //!< scalar gradient at Gauss point
        ) override;

    //! fill coupling vector and add variables to reaction in order to compute reaction values and
    //! derivatives
    void fill_coupling_vector_and_add_variables(Mat::MatListReactions& material_list_reactions,
        ScaTraEleReaManagerAdvReac& reaction_manager);

    //! Get right hand side including reaction body force term
    void get_rhs_int(double& rhsint,  //!< rhs containing body force at Gauss point
        double densnp,                //!< density at t_(n+1)
        int scalar_id                 //!< index of current scalar
        ) override;

    //! calculation of reactive element matrix
    void calc_mat_react(
        Core::LinAlg::SerialDenseMatrix& element_matrix,  //!< element matrix to be filled
        int scalar_id,                                    //!< index of current scalar
        double timefacfac,  //!< domain-integration factor times time-integration factor
        double timetaufac,  //!< domain-integration factor times time-integration factor times tau
        double taufac,      //!< domain-integration factor times tau
        double densnp,      //!< density at time_(n+1)
        const Core::LinAlg::Matrix<nen_, 1>& sgconv,  //!< subgrid-scale convective operator
        const Core::LinAlg::Matrix<nen_, 1>& diff     //!< laplace term
        ) override;

    //! "shape derivatives" pressure gradient
    void apply_shape_derivs_pressure_grad(Core::LinAlg::SerialDenseMatrix& element_matrix,
        int scalar_id, double vrhs, const Core::LinAlg::Matrix<nsd_, 1>& gradient_phi,
        Core::LinAlg::Matrix<nsd_, 1> pressure_gradient_reference);

    //! nodal flux values at t_(n+1)
    std::vector<Core::LinAlg::Matrix<nsd_, nen_>> element_flux_np_;

    //! a vector containing all quantities, the equation is coupled with
    //! (i.e. pressures, saturations and porosity)
    std::vector<std::pair<std::string, double>> coupling_values_;

    //! penalty factor to avoid very small "densities"
    const double penalty_ = 1.0;

    //! do we use L2-projection or evaluation at GP
    bool enable_L2_projection_ = false;
  };


  template <int nsd, int nen>
  class ScaTraEleInternalVariableManagerPorofluidPressureBased final
      : public ScaTraEleInternalVariableManager<nsd, nen>
  {
    using scatra_variable_manager = ScaTraEleInternalVariableManager<nsd, nen>;

   public:
    explicit ScaTraEleInternalVariableManagerPorofluidPressureBased(int number_of_scalars)
        : ScaTraEleInternalVariableManager<nsd, nen>(number_of_scalars),
          pressure_(0),
          saturation_(0),
          density_(0),
          solid_pressure_(0.0),
          delta_(number_of_scalars, 0.0),
          relative_mobility_funct_id_(number_of_scalars),
          heat_capacity_(0),
          thermal_diffusivity_(0),
          effective_heat_capacity_(0.0),
          minimal_value_of_phase_(number_of_scalars, 0.0),
          evaluate_scalar_(number_of_scalars, true),
          scalar_to_phase_map_(
              number_of_scalars, {-1, Mat::ScaTraMatMultiPoro::SpeciesType::species_undefined}),
          material_is_set_(false),
          action_(ScaTra::Action::calc_mat_and_rhs)
    {
    }

    //! compute and set internal variables -- no L2-projection but evaluation at GP
    void set_internal_variables_porofluid_pressure_based(
        const Core::LinAlg::Matrix<nen, 1>& shape_functions,  //! array for shape functions
        const Core::LinAlg::Matrix<nsd, nen>&
            shape_functions_deriv_xyz,  //! global derivatives of shape functions w.r.t x,y,z
        const Core::LinAlg::Matrix<nsd, nen>&
            shape_functions_deriv,  //! global derivatives of shape functions w.r.t r,s,t
        const Core::LinAlg::Matrix<nsd, nsd>& jacobian_matrix,
        const Core::LinAlg::Matrix<nsd, nen>& node_coordinates_reference,
        const std::vector<Core::LinAlg::Matrix<nen, 1>>&
            element_phi_np,  //! scalar at t_(n+1) or t_(n+alpha_F)
        const std::vector<Core::LinAlg::Matrix<nen, 1>>& element_phi_n,  //! scalar at t_(n)
        const std::vector<Core::LinAlg::Matrix<nen, 1>>&
            element_hist,  //! history vector of transported scalars
        const Core::LinAlg::Matrix<nsd, nen>&
            element_force_velocity  //! nodal velocity due to external force
    )
    {
      // call base class (scatra) with dummy variable
      const Core::LinAlg::Matrix<nsd, nen> dummy_element_conv;
      scatra_variable_manager::set_internal_variables(shape_functions, shape_functions_deriv_xyz,
          element_phi_np, element_phi_n, dummy_element_conv, element_hist, dummy_element_conv);

      // velocity due to the external force
      Core::LinAlg::Matrix<nsd, 1> force_velocity;
      force_velocity.multiply(element_force_velocity, shape_functions);

      // get determinant of Jacobian dX / ds
      // transposed jacobian "dX/ds"
      Core::LinAlg::Matrix<nsd, nsd> jacobian_matrix_reference;
      jacobian_matrix_reference.multiply_nt(shape_functions_deriv, node_coordinates_reference);

      // inverse of transposed jacobian "ds/dX"
      const double jacobian_determinant_reference = jacobian_matrix_reference.determinant();
      const double jacobian_determinant = jacobian_matrix.determinant();

      // determinant of deformation gradient
      // det F = det ( d x / d X ) = det (dx/ds) * ( det(dX/ds))^-1
      const double determinant_deformation_gradient =
          jacobian_determinant / jacobian_determinant_reference;

      // clear current gauss point data for safety
      phase_manager_->clear_gp_state();
      variable_manager_->evaluate_gp_variables(shape_functions, shape_functions_deriv_xyz);
      // access from outside to the phase manager:
      // scatra-discretization has fluid-discretization on dofset 2
      phase_manager_->evaluate_gp_state(
          determinant_deformation_gradient, *variable_manager_, nds_scatra_porofluid_);

      const int number_of_fluid_phases = phase_manager_->num_fluid_phases();
      const int number_of_volume_fractions = phase_manager_->num_vol_frac();

      // resize all phase related vectors
      pressure_.resize(number_of_fluid_phases);
      saturation_.resize(number_of_fluid_phases);
      density_.resize(number_of_fluid_phases + number_of_volume_fractions);
      heat_capacity_.resize(number_of_fluid_phases + number_of_volume_fractions + 1);
      thermal_diffusivity_.resize(number_of_fluid_phases + number_of_volume_fractions + 1);
      pressure_gradient_.resize(number_of_fluid_phases + number_of_volume_fractions);
      diffusion_tensor_porofluid_.resize(number_of_fluid_phases + number_of_volume_fractions);
      pressure_gradient_norm_.resize(number_of_fluid_phases);
      volume_fractions_.resize(number_of_volume_fractions);
      volume_fractions_pressure_.resize(number_of_volume_fractions);
      relative_mobility_funct_id_.resize(scatra_variable_manager::numscal_);

      const std::vector<Core::LinAlg::Matrix<nsd, 1>>& fluid_gradient_phi =
          *(variable_manager_->grad_phinp());

      volume_fractions_ = phase_manager_->vol_frac();
      volume_fractions_pressure_ = phase_manager_->vol_frac_pressure();

      //! convective velocity
      std::vector<Core::LinAlg::Matrix<nsd, 1>> convective_phase_velocity(0.0);
      convective_phase_velocity.resize(number_of_fluid_phases + number_of_volume_fractions);
      //! convective part in convective form: (u_x * N,x) + (u_y * N,y)
      std::vector<Core::LinAlg::Matrix<nen, 1>> convective_phase_velocity_convective_form(0.0);
      convective_phase_velocity_convective_form.resize(
          number_of_fluid_phases + number_of_volume_fractions);

      //! temperature convective velocity
      Core::LinAlg::Matrix<nsd, 1> convective_temperature_velocity;
      //! temperature convective part in convective form
      Core::LinAlg::Matrix<nen, 1> convective_temperature_velocity_convective_form;

      for (int phase_id = 0; phase_id < number_of_fluid_phases; ++phase_id)
      {
        // current pressure gradient
        pressure_gradient_[phase_id].clear();

        // phase density
        density_[phase_id] = phase_manager_->density(phase_id);

        // compute the pressure gradient from the phi gradients
        for (int dof_id = 0; dof_id < number_of_fluid_phases; ++dof_id)
          pressure_gradient_[phase_id].update(
              phase_manager_->pressure_deriv(phase_id, dof_id), fluid_gradient_phi[dof_id], 1.0);

        // compute the absolute value of the pressure gradient from the phi gradients
        pressure_gradient_norm_[phase_id] = 0.0;
        for (int i = 0; i < nsd; i++)
          pressure_gradient_norm_[phase_id] +=
              pressure_gradient_[phase_id](i) * pressure_gradient_[phase_id](i);
        pressure_gradient_norm_[phase_id] = sqrt(pressure_gradient_norm_[phase_id]);

        // diffusion tensor
        diffusion_tensor_porofluid_[phase_id].clear();
        phase_manager_->permeability_tensor(phase_id, diffusion_tensor_porofluid_[phase_id]);
        diffusion_tensor_porofluid_[phase_id].scale(
            phase_manager_->rel_permeability(phase_id) / phase_manager_->dyn_viscosity(phase_id,
                                                             pressure_gradient_norm_[phase_id],
                                                             nds_scatra_porofluid_));

        // Insert Darcy's law:
        // porosity * S^phase * (v^phase - v^s) = - permeability/viscosity * grad p
        convective_phase_velocity[phase_id].multiply(
            -1.0, diffusion_tensor_porofluid_[phase_id], pressure_gradient_[phase_id]);
        convective_temperature_velocity.update(heat_capacity_[phase_id] * density_[phase_id],
            convective_phase_velocity[phase_id], 1.0);
        // in convective form: (u_x * N,x) + (u_y * N,y)
        convective_phase_velocity_convective_form[phase_id].multiply_tn(
            shape_functions_deriv_xyz, convective_phase_velocity[phase_id]);
        convective_temperature_velocity_convective_form.update(
            heat_capacity_[phase_id] * density_[phase_id],
            convective_phase_velocity_convective_form[phase_id], 1.0);

        // phase pressure
        pressure_[phase_id] = phase_manager_->pressure(phase_id);
        // phase saturation
        saturation_[phase_id] = phase_manager_->saturation(phase_id);
      }

      for (int volume_fraction_id = number_of_fluid_phases;
          volume_fraction_id < number_of_fluid_phases + number_of_volume_fractions;
          ++volume_fraction_id)
      {
        // current pressure gradient
        pressure_gradient_[volume_fraction_id].update(
            1.0, fluid_gradient_phi[volume_fraction_id + number_of_volume_fractions], 0.0);

        // density of the volume fraction
        density_[volume_fraction_id] =
            phase_manager_->vol_frac_density(volume_fraction_id - number_of_fluid_phases);

        // diffusion tensor
        diffusion_tensor_porofluid_[volume_fraction_id].clear();
        phase_manager_->permeability_tensor_vol_frac_pressure(
            volume_fraction_id - number_of_fluid_phases,
            diffusion_tensor_porofluid_[volume_fraction_id]);
        // -1.0 --> don't need pressure gradient norm
        diffusion_tensor_porofluid_[volume_fraction_id].scale(
            1.0 / phase_manager_->dyn_viscosity_vol_frac_pressure(
                      volume_fraction_id - number_of_fluid_phases, -1.0, nds_scatra_porofluid_));

        // Insert Darcy's law: volume fraction * (v^phase - v^s) = - permeability/viscosity * grad p
        convective_phase_velocity[volume_fraction_id].multiply(-1.0,
            diffusion_tensor_porofluid_[volume_fraction_id],
            pressure_gradient_[volume_fraction_id]);
        convective_temperature_velocity.update(
            heat_capacity_[volume_fraction_id] * density_[volume_fraction_id],
            convective_phase_velocity[volume_fraction_id], 1.0);
        // in convective form: (u_x * N,x) + (u_y * N,y)
        convective_phase_velocity_convective_form[volume_fraction_id].multiply_tn(
            shape_functions_deriv_xyz, convective_phase_velocity[volume_fraction_id]);
        convective_temperature_velocity_convective_form.update(
            heat_capacity_[volume_fraction_id] * density_[volume_fraction_id],
            convective_phase_velocity_convective_form[volume_fraction_id], 1.0);
      }

      // solid pressure
      solid_pressure_ = phase_manager_->solid_pressure();

      for (int scalar_id = 0; scalar_id < scatra_variable_manager::numscal_; ++scalar_id)
      {
        // overwrite convective term: - rho * k/\mu*grad p * grad phi
        switch (scalar_to_phase_map_[scalar_id].species_type)
        {
          case Mat::ScaTraMatMultiPoro::SpeciesType::species_in_fluid:
          case Mat::ScaTraMatMultiPoro::SpeciesType::species_in_volfrac:
            scatra_variable_manager::convelint_[scalar_id] =
                convective_phase_velocity[scalar_to_phase_map_[scalar_id].phaseID];
            // if the scalar reacts to the external force, add the velocity due to the external
            // force scaled with the relative mobility and the porosity * saturation
            if (scatra_variable_manager::reacts_to_force_[scalar_id])
            {
              const auto prefactor =
                  evaluate_relative_mobility(scalar_id) * phase_manager_->porosity() *
                  phase_manager_->saturation(scalar_to_phase_map_[scalar_id].phaseID);
              scatra_variable_manager::convelint_[scalar_id].update(prefactor, force_velocity, 1.0);
            }
            scatra_variable_manager::conv_[scalar_id].multiply_tn(
                shape_functions_deriv_xyz, scatra_variable_manager::convelint_[scalar_id]);
            scatra_variable_manager::conv_phi_[scalar_id] =
                scatra_variable_manager::convelint_[scalar_id].dot(
                    scatra_variable_manager::gradphi_[scalar_id]);
            break;

          case Mat::ScaTraMatMultiPoro::SpeciesType::species_in_solid:
            scatra_variable_manager::convelint_[scalar_id] =
                Core::LinAlg::Matrix<nsd, 1>(Core::LinAlg::Initialization::zero);
            scatra_variable_manager::conv_[scalar_id] =
                Core::LinAlg::Matrix<nen, 1>(Core::LinAlg::Initialization::zero);
            scatra_variable_manager::conv_phi_[scalar_id] = 0;
            break;

          case Mat::ScaTraMatMultiPoro::SpeciesType::species_temperature:
            scatra_variable_manager::convelint_[scalar_id] = convective_temperature_velocity;
            scatra_variable_manager::conv_[scalar_id] =
                convective_temperature_velocity_convective_form;
            scatra_variable_manager::conv_phi_[scalar_id] =
                convective_temperature_velocity.dot(scatra_variable_manager::gradphi_[scalar_id]);
            break;

          default:
            FOUR_C_THROW("Unknown species type {} for species {}.",
                scalar_to_phase_map_[scalar_id].species_type, scalar_id);
        }
        // set flag if we actually have to evaluate the species
        set_evaluate_scalar_flag(scalar_id);
      }
    }

    //! adapt convective term in case of L2-projection
    void adapt_convective_term_for_l2(
        const Core::LinAlg::Matrix<nen, 1>& shape_functions,  //! array for shape functions
        const Core::LinAlg::Matrix<nsd, nen>&
            shape_functions_deriv_xyz,  //! global derivatives of shape functions w.r.t x,y,z
        const std::vector<Core::LinAlg::Matrix<nsd, nen>>&
            element_flux_np  //! nodal flux values at t_(n+1) or t_(n+alpha_F)
    )
    {
      const int number_of_fluid_phases = element_flux_np.size();

      std::vector<Core::LinAlg::Matrix<nsd, 1>> flux(0.0);
      flux.resize(number_of_fluid_phases);

      // in convective form: q_x*N,x + q_y*N,y
      std::vector<Core::LinAlg::Matrix<nen, 1>> flux_conv(0.0);
      flux_conv.resize(number_of_fluid_phases);

      for (int phase_id = 0; phase_id < number_of_fluid_phases; ++phase_id)
      {
        flux[phase_id].multiply(1.0, element_flux_np[phase_id], shape_functions);
        flux_conv[phase_id].multiply_tn(shape_functions_deriv_xyz, flux[phase_id]);
      }

      for (const auto& [phase_id, _] : scalar_to_phase_map_)
      {
        if (phase_id < 0 || phase_id >= number_of_fluid_phases)
          FOUR_C_THROW("Invalid phase ID {}", phase_id);
      }

      // set convective term
      for (int scalar_id = 0; scalar_id < scatra_variable_manager::numscal_; ++scalar_id)
      {
        switch (scalar_to_phase_map_[scalar_id].species_type)
        {
          case Mat::ScaTraMatMultiPoro::SpeciesType::species_in_fluid:
          case Mat::ScaTraMatMultiPoro::SpeciesType::species_in_volfrac:
            scatra_variable_manager::conv_phi_[scalar_id] =
                flux[scalar_to_phase_map_[scalar_id].phaseID].dot(
                    scatra_variable_manager::gradphi_[scalar_id]);
            break;

          case Mat::ScaTraMatMultiPoro::SpeciesType::species_in_solid:
          case Mat::ScaTraMatMultiPoro::SpeciesType::species_temperature:
            scatra_variable_manager::conv_phi_[scalar_id] = 0;
            break;

          default:
            FOUR_C_THROW("Unknown species type {} for species {}.",
                scalar_to_phase_map_[scalar_id].species_type, scalar_id);
        }
      }
    }

    // Set the porofluid-material in the scatra variable manager
    void set_porofluid_material(const Core::Elements::Element* element)
    {
      // check if we actually have three materials
      if (element->num_material() < 3) FOUR_C_THROW("no third material available");

      // requires PoroMultiPhase material to be the third material
      porofluid_material_ = std::dynamic_pointer_cast<Mat::FluidPoroMultiPhase>(
          element->material(nds_scatra_porofluid_));
      if (porofluid_material_ == nullptr) FOUR_C_THROW("cast to Mat::FluidPoroMultiPhase failed!");

      material_is_set_ = true;
    }

    //! get pressure associated with scalar
    [[nodiscard]] double get_pressure(const int scalar_id) const
    {
      switch (scalar_to_phase_map_[scalar_id].species_type)
      {
        case Mat::ScaTraMatMultiPoro::SpeciesType::species_in_fluid:
          return pressure_[scalar_to_phase_map_[scalar_id].phaseID];

        default:
          FOUR_C_THROW("scalar_to_phase_id = {} for species {}",
              scalar_to_phase_map_[scalar_id].phaseID, scalar_id);
      }
    }

    //! get pressures
    [[nodiscard]] const std::vector<double>& get_pressure() const { return pressure_; }

    //! get saturation associated with scalar
    [[nodiscard]] double get_saturation(const int scalar_id) const
    {
      switch (scalar_to_phase_map_[scalar_id].species_type)
      {
        case Mat::ScaTraMatMultiPoro::SpeciesType::species_in_fluid:
          return saturation_[scalar_to_phase_map_[scalar_id].phaseID];

        default:
          FOUR_C_THROW("scalar_to_phase_id = {} for species {}",
              scalar_to_phase_map_[scalar_id].phaseID, scalar_id);
      }
    }

    //! get saturation associated with scalar
    [[nodiscard]] const std::vector<double>& get_saturation() const { return saturation_; }

    //! get density associated with scalar
    [[nodiscard]] double get_density(const int scalar_id) const
    {
      switch (scalar_to_phase_map_[scalar_id].species_type)
      {
        case Mat::ScaTraMatMultiPoro::SpeciesType::species_in_fluid:
        case Mat::ScaTraMatMultiPoro::SpeciesType::species_in_volfrac:
          return density_[scalar_to_phase_map_[scalar_id].phaseID];

        case Mat::ScaTraMatMultiPoro::SpeciesType::species_in_solid:
          return phase_manager_->solid_density();

        case Mat::ScaTraMatMultiPoro::SpeciesType::species_temperature:
          // set to 1.0 because densities are included in the effective heat capacity
          return 1.0;

        default:
          FOUR_C_THROW("Unknown species type {} for species {}.",
              scalar_to_phase_map_[scalar_id].species_type, scalar_id);
      }
    }

    //! return density vector
    [[nodiscard]] const std::vector<double>& get_density() const { return density_; }

    //! return volfrac vector
    [[nodiscard]] const std::vector<double>& get_volume_fraction() const
    {
      return volume_fractions_;
    }

    //! get volume fraction associated with scalar
    [[nodiscard]] double get_volume_fraction(const int scalar_id) const
    {
      switch (scalar_to_phase_map_[scalar_id].species_type)
      {
        // case Mat::ScaTraMatMultiPoro::SpeciesType::species_in_fluid:
        case Mat::ScaTraMatMultiPoro::SpeciesType::species_in_volfrac:
          return volume_fractions_[get_phase_id(scalar_id) - phase_manager_->num_fluid_phases()];

        default:
          FOUR_C_THROW("scalar_to_phase_id = {} for species {}",
              scalar_to_phase_map_[scalar_id].phaseID, scalar_id);
      }
    }

    //! get volfrac pressure vector
    [[nodiscard]] const std::vector<double>& get_volume_fractions_pressure() const
    {
      return volume_fractions_pressure_;
    }

    //! get solid pressure
    [[nodiscard]] double get_solid_pressure() const { return solid_pressure_; };

    //! set scalar ID to phase ID mapping and species type
    void set_phase_id_and_species_type(const int scalar_id, const int phase_id,
        const Mat::ScaTraMatMultiPoro::SpeciesType& species_type)
    {
      scalar_to_phase_map_[scalar_id].phaseID = phase_id;
      scalar_to_phase_map_[scalar_id].species_type = species_type;
    }

    //! get phase ID from scalar ID
    [[nodiscard]] int get_phase_id(const int scalar_id) const
    {
      if (scalar_to_phase_map_[scalar_id].phaseID < 0)
        FOUR_C_THROW("scalar_to_phase_id = {} for species {}",
            scalar_to_phase_map_[scalar_id].phaseID, scalar_id);

      return scalar_to_phase_map_[scalar_id].phaseID;
    }

    //! get species type of scalar
    [[nodiscard]] Mat::ScaTraMatMultiPoro::SpeciesType get_species_type(const int scalar_id) const
    {
      return scalar_to_phase_map_[scalar_id].species_type;
    }

    //! set delta for evaluation of effective diffusivity
    void set_delta(const double delta, const int scalar_id) { delta_[scalar_id] = delta; }

    //! set relative mobility function ID
    void set_relative_mobility_function_id(
        const int relative_mobility_funct_id, const int scalar_id)
    {
      relative_mobility_funct_id_[scalar_id] = relative_mobility_funct_id;
    }

    //! set heat capacity
    void set_heat_capacity(const std::vector<double>& heat_capacity)
    {
      heat_capacity_ = heat_capacity;
    }

    //! set effective heat capacity
    void set_effective_heat_capacity(const double effective_heat_capacity)
    {
      effective_heat_capacity_ = effective_heat_capacity;
    }

    //! set thermal diffusivity
    void set_thermal_diffusivity(const std::vector<double>& kappa) { thermal_diffusivity_ = kappa; }

    //! set minimum value of corresponding phase under which we assume that mass fraction is zero
    void set_minimal_value_of_phase(const double minimal_value_of_phase, const int scalar_id)
    {
      minimal_value_of_phase_[scalar_id] = minimal_value_of_phase;
    }

    //! set action
    void set_action(const ScaTra::Action action) { action_ = action; };

    //! get delta
    [[nodiscard]] double get_delta(const int phase_id) const { return delta_[phase_id]; }

    //! get heat capacity
    //! order [ <fluid>  <volfrac>  <solid> ]
    [[nodiscard]] double get_heat_capacity(const int phase_id) const
    {
      return heat_capacity_[phase_id];
    }

    //! get effective heat capacity
    [[nodiscard]] double get_effective_heat_capacity() const { return effective_heat_capacity_; }

    //! get thermal diffusivity
    //! order [ <fluid>  <volfrac>  <solid> ]
    [[nodiscard]] double get_thermal_diffusivity(const int phase_id) const
    {
      return thermal_diffusivity_[phase_id];
    }

    //! get evaluate scalar flag
    bool evaluate_scalar(const int scalar_id) { return evaluate_scalar_[scalar_id]; }

    //! get minimum value of corresponding phase under which we assume that mass fraction is zero
    [[nodiscard]] double get_minimal_value_of_phase(const int scalar_id) const
    {
      return minimal_value_of_phase_[scalar_id];
    }

    //! get pre-factor needed for off-diagonal mesh-linearization of mass matrix
    double get_pre_factor_mass_matrix_od_mesh(const int scalar_id, const double fac)
    {
      double pre_factor = fac;

      switch (scalar_to_phase_map_[scalar_id].species_type)
      {
        case Mat::ScaTraMatMultiPoro::SpeciesType::species_in_fluid:
        {
          // linearization of porosity
          //
          // dporosity/dd = dporosity/dJ * dJ/dd = dporosity/dJ * J * N_x
          // J denotes the determinant of the deformation gradient, i.e.,
          // det F = det ( d x / d X ) = det (dx/ds) * ( det(dX/ds) )^-1
          // in our case: pre_factor is scaled with density, i.e., porosity
          // --> scale it with 1.0/porosity here
          //
          // We will pass fac/porosity * J * dporosity/dJ into CalcMatMassODMesh, where it will
          // internally be scaled correctly.
          if (phase_manager_->porosity_depends_on_struct())
            pre_factor += fac / (phase_manager_->porosity()) * phase_manager_->jacobian_def_grad() *
                          phase_manager_->porosity_deriv_wrt_jacobian_def_grad();
          break;
        }
        case Mat::ScaTraMatMultiPoro::SpeciesType::species_in_volfrac:
        {
          // do nothing: correct pre_factor = fac
          break;
        }
        case Mat::ScaTraMatMultiPoro::SpeciesType::species_in_solid:
        {
          if (phase_manager_->porosity_depends_on_struct())
          {
            pre_factor -= fac /
                          (1 - phase_manager_->porosity() - phase_manager_->sum_add_vol_frac()) *
                          phase_manager_->jacobian_def_grad() *
                          phase_manager_->porosity_deriv_wrt_jacobian_def_grad();
          }
          break;
        }
        case Mat::ScaTraMatMultiPoro::SpeciesType::species_temperature:
        {
          if (phase_manager_->porosity_depends_on_struct())
          {
            // evaluate the effective heat capacity factor
            const int total_number_of_porofluid_phases =
                phase_manager_->num_fluid_phases() + phase_manager_->num_vol_frac();
            const int number_of_fluid_phases = phase_manager_->num_fluid_phases();

            const auto effective_heat_capacity = get_effective_heat_capacity();

            pre_factor += fac / effective_heat_capacity * phase_manager_->solid_density() *
                          get_heat_capacity(total_number_of_porofluid_phases) * (-1.0) *
                          phase_manager_->jacobian_def_grad() *
                          phase_manager_->porosity_deriv_wrt_jacobian_def_grad();

            for (int phase = 0; phase < number_of_fluid_phases; ++phase)
            {
              pre_factor += fac / effective_heat_capacity * phase_manager_->density(phase) *
                            get_heat_capacity(phase) * phase_manager_->jacobian_def_grad() *
                            phase_manager_->saturation(phase) *
                            phase_manager_->porosity_deriv_wrt_jacobian_def_grad();
            }
          }
          break;
        }
        default:
          FOUR_C_THROW("Unknown species type {} for species {}.",
              scalar_to_phase_map_[scalar_id].species_type, scalar_id);
      }
      return pre_factor;
    }

    //! get pre-factor needed for off-diagonal porofluid-linearization of mass matrix
    void get_pre_factor_mass_matrix_od_porofluid(
        const int scalar_id, std::vector<double>* pre_factor_mass_matrix_od_porofluid)
    {
      // reset to zero
      std::ranges::fill(pre_factor_mass_matrix_od_porofluid->begin(),
          pre_factor_mass_matrix_od_porofluid->end(), 0.0);

      const int number_of_fluid_phases = phase_manager_->num_fluid_phases();

      switch (scalar_to_phase_map_[scalar_id].species_type)
      {
        case Mat::ScaTraMatMultiPoro::SpeciesType::species_in_fluid:
        {
          const int phase_id = get_phase_id(scalar_id);
          const int total_number_of_porofluid_dofs = porofluid_material()->num_mat();

          // d poro / d S_j = SaturationDeriv
          // --> scaling with 1.0/S since term is re-scaled with densnp = porosity * S * rho
          for (int dof_id = 0; dof_id < number_of_fluid_phases; ++dof_id)
          {
            (*pre_factor_mass_matrix_od_porofluid)[dof_id] +=
                phase_manager_->saturation_deriv(phase_id, dof_id) /
                phase_manager_->saturation(phase_id);
          }

          // linearization of porosity (only if porosity is pressure-dependent)
          // d poro / d psi_j = porosityDeriv
          // --> scaling with 1.0/porosity  since term is re-scaled with densnp = porosity * S *rho
          if (phase_manager_->porosity_depends_on_fluid())
          {
            for (int dof_id = 0; dof_id < total_number_of_porofluid_dofs; ++dof_id)
              (*pre_factor_mass_matrix_od_porofluid)[dof_id] +=
                  phase_manager_->porosity_deriv(dof_id) / phase_manager_->porosity();
          }
          break;
        }
        case Mat::ScaTraMatMultiPoro::SpeciesType::species_in_volfrac:
        {
          const int phase_id = get_phase_id(scalar_id);
          // d volfrac_j / d volfrac_j = 1.0
          // --> scaling with 1.0/volfrac since term is re-scaled with densnp = rho * volfrac
          (*pre_factor_mass_matrix_od_porofluid)[phase_id] += 1.0 / get_volume_fraction(scalar_id);

          break;
        }
        case Mat::ScaTraMatMultiPoro::SpeciesType::species_in_solid:
        {
          // porosity of solid epsilon_s does not depend on porofluid primary variables
          break;
        }
        case Mat::ScaTraMatMultiPoro::SpeciesType::species_temperature:
        {
          // evaluate the effective heat capacity factor
          const int total_number_of_porofluid_phases =
              phase_manager_->num_fluid_phases() + phase_manager_->num_vol_frac();
          // total number of porofluid dofs additionally includes volume fraction pressures
          const int total_number_of_porofluid_dofs = porofluid_material()->num_mat();

          const auto effective_heat_capacity = get_effective_heat_capacity();

          for (int phase = 0; phase < number_of_fluid_phases; ++phase)
          {
            for (int dof_id = 0; dof_id < number_of_fluid_phases; ++dof_id)
            {
              (*pre_factor_mass_matrix_od_porofluid)[dof_id] +=
                  phase_manager_->saturation_deriv(phase, dof_id) / effective_heat_capacity *
                  phase_manager_->density(phase) * get_heat_capacity(phase) *
                  phase_manager_->porosity();
            }

            if (phase_manager_->porosity_depends_on_fluid())
            {
              for (int dof_id = 0; dof_id < total_number_of_porofluid_dofs; ++dof_id)
              {
                (*pre_factor_mass_matrix_od_porofluid)[dof_id] +=
                    phase_manager_->porosity_deriv(dof_id) / effective_heat_capacity *
                    phase_manager_->density(phase) * get_heat_capacity(phase) *
                    phase_manager_->saturation(phase);
              }
            }
          }

          for (int phase = number_of_fluid_phases; phase < total_number_of_porofluid_phases;
              ++phase)
          {
            (*pre_factor_mass_matrix_od_porofluid)[phase] +=
                get_heat_capacity(phase) *
                phase_manager_->vol_frac_density(phase - number_of_fluid_phases) /
                effective_heat_capacity;
          }
          break;
        }
        default:
          FOUR_C_THROW("Unknown species type {} for species {}.",
              scalar_to_phase_map_[scalar_id].species_type, scalar_id);
      }
    }

    //! evaluate relative mobility
    auto evaluate_relative_mobility(const int current_scalar) const
    {
      std::vector<std::pair<std::string, double>> varfunction_variables;
      std::vector<std::pair<std::string, double>> varfunction_constants;
      const auto number_of_fluid_phases = phase_manager_->num_fluid_phases();

      for (int phase_id = 0; phase_id < number_of_fluid_phases; phase_id++)
      {
        std::ostringstream temp;
        temp << phase_id + 1;

        varfunction_variables.emplace_back("S" + temp.str(), phase_manager_->saturation(phase_id));
        varfunction_variables.emplace_back("p" + temp.str(), phase_manager_->pressure(phase_id));
      }
      varfunction_variables.emplace_back("porosity", phase_manager_->porosity());

      const auto relative_mobility = Global::Problem::instance()
                                         ->function_by_id<Core::Utils::FunctionOfAnything>(
                                             relative_mobility_funct_id_[current_scalar])
                                         .evaluate(varfunction_variables, varfunction_constants, 0);

      return relative_mobility;
    }

    //! get pre-factor needed for off-diagonal fluid-linearization of convective term
    void get_pre_factor_convection_od_porofluid(const unsigned ui,
        std::vector<double>* pre_factor_convection_od_porofluid,
        const Core::LinAlg::Matrix<nsd, 1>& grad_phi,
        const Core::LinAlg::Matrix<1, nsd>& grad_phi_transpose_diffusion_tensor,
        const Core::LinAlg::Matrix<nen, 1>& shape_functions,
        const Core::LinAlg::Matrix<nsd, nen>& shape_functions_deriv_xyz, const int phase_id)
    {
      // reset to zero
      std::ranges::fill(pre_factor_convection_od_porofluid->begin(),
          pre_factor_convection_od_porofluid->end(), 0.0);

      // get correct factor
      double laplacian(0.0);
      for (int dim = 0; dim < nsd; dim++)
        laplacian +=
            shape_functions_deriv_xyz(dim, ui) * grad_phi_transpose_diffusion_tensor(0, dim);

      const int number_of_fluid_phases = phase_manager_->num_fluid_phases();
      const int number_of_volume_fractions = phase_manager_->num_vol_frac();

      // fluid phase
      if (phase_id >= 0 && phase_id < number_of_fluid_phases)
      {
        // derivative after fluid pressures
        for (int dof_id = 0; dof_id < number_of_fluid_phases; ++dof_id)
          (*pre_factor_convection_od_porofluid)[dof_id] +=
              laplacian * phase_manager_->pressure_deriv(phase_id, dof_id);

        // derivative after relative permeability
        if (not phase_manager_->has_constant_rel_permeability(phase_id))
        {
          const Core::LinAlg::Matrix<nsd, 1> pressure_gradient = get_pressure_gradient(phase_id);
          const double pressure_gradient_norm = get_pressure_gradient_norm(phase_id);

          static Core::LinAlg::Matrix<nsd, nsd> diffusion_tensor;
          phase_manager_->permeability_tensor(phase_id, diffusion_tensor);
          diffusion_tensor.scale(
              phase_manager_->rel_permeability_deriv(phase_id) /
              phase_manager_->dyn_viscosity(phase_id, pressure_gradient_norm, 2));

          static Core::LinAlg::Matrix<1, nsd> phase_grad_phi_transpose_diffusion_tensor;
          phase_grad_phi_transpose_diffusion_tensor.multiply_tn(grad_phi, diffusion_tensor);

          double phase_laplacian(0.0);
          for (unsigned j = 0; j < nsd; j++)
            phase_laplacian +=
                pressure_gradient(j) * phase_grad_phi_transpose_diffusion_tensor(0, j);

          for (int dof_id = 0; dof_id < number_of_fluid_phases; ++dof_id)
          {
            (*pre_factor_convection_od_porofluid)[dof_id] +=
                phase_laplacian * shape_functions(ui) *
                phase_manager_->saturation_deriv(phase_id, dof_id);
          }
        }
      }
      // volume fraction (in additional porous network)
      else if (phase_id < number_of_fluid_phases + number_of_volume_fractions)
      {
        // derivative after volfrac-pressure
        (*pre_factor_convection_od_porofluid)[phase_id + number_of_volume_fractions] += laplacian;
      }
      else
      {
        FOUR_C_THROW(
            "get_pre_factor_convection_od_porofluid has been called with phase {}", phase_id);
      }
    }

    //! get pre-factor needed for off-diagonal mesh-linearization of diffusive term
    double get_pre_factor_diffusion_od_mesh(
        const int scalar_id, const double rhsfac, const double diffusion_coefficient)
    {
      double pre_factor = 0.0;

      switch (scalar_to_phase_map_[scalar_id].species_type)
      {
        case Mat::ScaTraMatMultiPoro::SpeciesType::species_in_fluid:
        {
          if (phase_manager_->porosity_depends_on_struct())
          {
            const double delta = get_delta(scalar_id);

            // linearization of porosity
            // dporosity/dd = dporosity/dJ * dJ/dd = dporosity/dJ * J * N_x
            // J denotes the determinant of the deformation gradient, i.e.,
            // det F = det ( d x / d X ) = det (dx/ds) * ( det(dX/ds) )^-1
            //
            // in our case: diffusivity is scaled with porosity^(delta + 1)
            // --> scale it with 1.0/porosity^(delta + 1) here and build derivative
            // d diff/d porosity = (delta + 1) * porosity^delta
            pre_factor = rhsfac * diffusion_coefficient /
                         std::pow(phase_manager_->porosity(), delta + 1.0) * (delta + 1.0) *
                         std::pow(phase_manager_->porosity(), delta) *
                         phase_manager_->jacobian_def_grad() *
                         phase_manager_->porosity_deriv_wrt_jacobian_def_grad();
          }
          break;
        }
        case Mat::ScaTraMatMultiPoro::SpeciesType::species_in_volfrac:
        {
          // do nothing: linearization performed in base class
          break;
        }
        case Mat::ScaTraMatMultiPoro::SpeciesType::species_in_solid:
        {
          if (phase_manager_->porosity_depends_on_struct())
          {
            const double delta = get_delta(scalar_id);

            // linearization of porosity
            // dporosity/dd = dporosity/dJ * dJ/dd = dporosity/dJ * J * N_x
            // J denotes the determinant of the deformation gradient, i.e.,
            // det F = det ( d x / d X ) = det (dx/ds) * ( det(dX/ds) )^-1
            //
            // in our case: diffusivity is scaled with (1-porosity-sum_volfrac)^(delta + 1)
            // --> scale it with 1.0/(1-porosity-sum_volfrac)^(delta + 1) here and build
            // derivative
            // d diff/d porosity = (delta + 1) * (1-porosity-sum_volfrac)^delta
            pre_factor =
                -1.0 * rhsfac * diffusion_coefficient /
                std::pow(1.0 - phase_manager_->porosity() - phase_manager_->sum_add_vol_frac(),
                    delta + 1.0) *
                (delta + 1.0) *
                std::pow(
                    1.0 - phase_manager_->porosity() - phase_manager_->sum_add_vol_frac(), delta) *
                phase_manager_->jacobian_def_grad() *
                phase_manager_->porosity_deriv_wrt_jacobian_def_grad();
          }
          break;
        }
        case Mat::ScaTraMatMultiPoro::SpeciesType::species_temperature:
        {
          if (phase_manager_->porosity_depends_on_struct())
          {
            const int total_number_of_porofluid_phases =
                phase_manager_->num_fluid_phases() + phase_manager_->num_vol_frac();
            const int number_of_fluid_phases = phase_manager_->num_fluid_phases();

            pre_factor = rhsfac * get_thermal_diffusivity(total_number_of_porofluid_phases) *
                         (-1.0) * phase_manager_->jacobian_def_grad() *
                         phase_manager_->porosity_deriv_wrt_jacobian_def_grad();

            for (int phase = 0; phase < number_of_fluid_phases; ++phase)
            {
              pre_factor += rhsfac * get_thermal_diffusivity(phase) *
                            phase_manager_->jacobian_def_grad() *
                            phase_manager_->saturation(phase) *
                            phase_manager_->porosity_deriv_wrt_jacobian_def_grad();
            }
          }
          break;
        }
        default:
          FOUR_C_THROW("Unknown species type {} for species {}.",
              scalar_to_phase_map_[scalar_id].species_type, scalar_id);
      }
      return pre_factor;
    }

    //! get pre-factor needed for off-diagonal porofluid linearization of diffusive term
    void get_pre_factor_diffusion_od_porofluid(const int scalar_id, const double rhsfac,
        const double diffusion_coefficient, std::vector<double>* pre_factor_diffusion_OD_porofluid)
    {
      // reset to zero
      std::ranges::fill(pre_factor_diffusion_OD_porofluid->begin(),
          pre_factor_diffusion_OD_porofluid->end(), 0.0);

      switch (scalar_to_phase_map_[scalar_id].species_type)
      {
        case Mat::ScaTraMatMultiPoro::SpeciesType::species_in_fluid:
        {
          const int phase_id = get_phase_id(scalar_id);
          const int number_of_fluid_phases = phase_manager_->num_fluid_phases();
          const int total_number_of_porofluid_dofs = porofluid_material()->num_mat();

          const double delta = get_delta(scalar_id);

          // linearization of saturation * porosity * d_eff = (saturation * porosity)^(delta+1)
          // w.r.t saturation
          //
          // in our case: diffusivity is scaled with saturation^(delta + 1)
          // --> scale it with 1.0/saturation^(delta + 1) here and build derivative
          // d diff/d saturation = (delta + 1) * saturation^delta
          const double vrhs_sat = rhsfac * diffusion_coefficient /
                                  std::pow(phase_manager_->saturation(phase_id), delta + 1.0) *
                                  (delta + 1.0) *
                                  std::pow(phase_manager_->saturation(phase_id), delta);

          // linearization of saturation * porosity * d_eff = (saturation * porosity)^(delta+1)
          // w.r.t porosity
          //
          // in our case: diffusivity is scaled with porosity^(delta + 1)
          // --> scale it with 1.0/porosity^(delta + 1) here and build derivative
          // d diff/d porosity = (delta + 1) * porosity^delta
          const double vrhs_poro = rhsfac * diffusion_coefficient /
                                   std::pow(phase_manager_->porosity(), delta + 1.0) *
                                   (delta + 1.0) * std::pow(phase_manager_->porosity(), delta);

          for (int dof_id = 0; dof_id < number_of_fluid_phases; ++dof_id)
            (*pre_factor_diffusion_OD_porofluid)[dof_id] +=
                vrhs_sat * phase_manager_->saturation_deriv(phase_id, dof_id);

          // linearization of porosity (only if porosity is pressure-dependent)
          if (phase_manager_->porosity_depends_on_fluid())
          {
            for (int dof_id = 0; dof_id < total_number_of_porofluid_dofs; ++dof_id)
              (*pre_factor_diffusion_OD_porofluid)[dof_id] +=
                  vrhs_poro * phase_manager_->porosity_deriv(dof_id);
          }
          break;
        }
        case Mat::ScaTraMatMultiPoro::SpeciesType::species_in_volfrac:
        {
          const int phase_id = get_phase_id(scalar_id);
          const double delta = get_delta(scalar_id);
          // diffusivity is scaled with volfrac^(delta + 1)
          // --> scale it with 1.0/volfrac^(delta + 1) here and build derivative
          // d diff/d volfrac = (delta + 1) * volfrac^delta
          (*pre_factor_diffusion_OD_porofluid)[phase_id] +=
              rhsfac * diffusion_coefficient /
              std::pow(get_volume_fraction(scalar_id), delta + 1.0) * (delta + 1.0) *
              std::pow(get_volume_fraction(scalar_id), delta);
          break;
        }
        case Mat::ScaTraMatMultiPoro::SpeciesType::species_in_solid:
        {
          // porosity of solid epsilon_s does not depend on fluid primary variables
          break;
        }
        case Mat::ScaTraMatMultiPoro::SpeciesType::species_temperature:
        {
          const int number_of_fluid_phases = phase_manager_->num_fluid_phases();
          const int total_number_of_porofluid_phases =
              number_of_fluid_phases + phase_manager_->num_vol_frac();
          // total number of porofluid dofs additionally includes volume fraction pressures
          const int total_number_of_porofluid_dofs = porofluid_material()->num_mat();

          for (int phase_id = 0; phase_id < number_of_fluid_phases; ++phase_id)
          {
            for (int dof_id = 0; dof_id < number_of_fluid_phases; ++dof_id)
            {
              (*pre_factor_diffusion_OD_porofluid)[dof_id] +=
                  rhsfac * phase_manager_->saturation_deriv(phase_id, dof_id) *
                  get_thermal_diffusivity(phase_id) * phase_manager_->porosity();
            }

            if (phase_manager_->porosity_depends_on_fluid())
            {
              for (int dof_id = 0; dof_id < total_number_of_porofluid_dofs; ++dof_id)
              {
                (*pre_factor_diffusion_OD_porofluid)[dof_id] +=
                    rhsfac * phase_manager_->porosity_deriv(dof_id) *
                    get_thermal_diffusivity(phase_id) * phase_manager_->saturation(phase_id);
              }
            }
          }

          for (int phase_id = number_of_fluid_phases; phase_id < total_number_of_porofluid_phases;
              ++phase_id)
            (*pre_factor_diffusion_OD_porofluid)[phase_id] +=
                rhsfac * get_thermal_diffusivity(phase_id);
          break;
        }
        default:
          FOUR_C_THROW("Unknown species type {} for species {}.",
              scalar_to_phase_map_[scalar_id].species_type, scalar_id);
      }
    }

    //! get diffusion tensor of fluid phases (permeability * relative permeability)/viscosity
    void get_diffusion_tensor_fluid(
        const int scalar_id, Core::LinAlg::Matrix<nsd, nsd>& diffusion_tensor, const int phase_id)
    {
      switch (scalar_to_phase_map_[scalar_id].species_type)
      {
        case Mat::ScaTraMatMultiPoro::SpeciesType::species_in_fluid:
        case Mat::ScaTraMatMultiPoro::SpeciesType::species_in_volfrac:
        case Mat::ScaTraMatMultiPoro::SpeciesType::species_temperature:
        {
          diffusion_tensor = diffusion_tensor_porofluid_[phase_id];
          break;
        }
        case Mat::ScaTraMatMultiPoro::SpeciesType::species_in_solid:
        default:
          FOUR_C_THROW("Unknown species type {} for species {}.",
              scalar_to_phase_map_[scalar_id].species_type, scalar_id);
      }
    }

    //! compute the pressure gradient in reference configuration
    void get_pressure_gradient_reference(const Core::LinAlg::Matrix<nsd, nsd>& jacobian_matrix,
        Core::LinAlg::Matrix<nsd, 1>& pressure_gradient_reference, const int phase_id)
    {
      pressure_gradient_reference.clear();

      const int number_of_fluid_phases = phase_manager_->num_fluid_phases();
      const int number_of_volume_fractions = phase_manager_->num_vol_frac();

      // fluid phase
      if (phase_id >= 0 && phase_id < number_of_fluid_phases)
      {
        const std::vector<Core::LinAlg::Matrix<nsd, 1>>& gradient_fluid_phi =
            *(variable_manager_->grad_phinp());

        // gradient of phi w.r.t. reference coordinates
        std::vector<Core::LinAlg::Matrix<nsd, 1>> gradient_fluid_phi_reference(
            number_of_fluid_phases,
            Core::LinAlg::Matrix<nsd, 1>(Core::LinAlg::Initialization::zero));
        for (int dof_id = 0; dof_id < number_of_fluid_phases; ++dof_id)
          gradient_fluid_phi_reference[dof_id].multiply(
              jacobian_matrix, gradient_fluid_phi[dof_id]);

        // compute the pressure gradient from the phi gradients
        for (int dof_id = 0; dof_id < number_of_fluid_phases; ++dof_id)
          pressure_gradient_reference.update(phase_manager_->pressure_deriv(phase_id, dof_id),
              gradient_fluid_phi_reference[dof_id], 1.0);
      }

      // volume fraction (in additional porous network)
      else if (phase_id < number_of_fluid_phases + number_of_volume_fractions)
      {
        pressure_gradient_reference.multiply(jacobian_matrix, pressure_gradient_[phase_id]);
      }

      else
      {
        FOUR_C_THROW("GetRefGradPres has been called with phase {}", phase_id);
      }
    }

    // set flag if we actually evaluate the scalar or set it to zero
    void set_evaluate_scalar_flag(const int scalar_id)
    {
      const int number_of_fluid_phases = phase_manager_->num_fluid_phases();

      switch (scalar_to_phase_map_[scalar_id].species_type)
      {
        case Mat::ScaTraMatMultiPoro::SpeciesType::species_in_fluid:
        {
          // we do not evaluate if smaller than threshold (at GP)
          evaluate_scalar_[scalar_id] =
              (fabs(saturation_[scalar_to_phase_map_[scalar_id].phaseID]) >
                  minimal_value_of_phase_[scalar_id]);
          break;
        }
        case Mat::ScaTraMatMultiPoro::SpeciesType::species_in_volfrac:
        {
          // we do not evaluate if smaller than minimum nodal volume fraction in element
          evaluate_scalar_[scalar_id] = variable_manager_->element_has_valid_vol_frac_species(
              scalar_to_phase_map_[scalar_id].phaseID - number_of_fluid_phases);
          break;
        }
        case Mat::ScaTraMatMultiPoro::SpeciesType::species_in_solid:
        case Mat::ScaTraMatMultiPoro::SpeciesType::species_temperature:
        {
          evaluate_scalar_[scalar_id] = true;
          break;
        }
        default:
          FOUR_C_THROW("Unknown species type {} for species {}.",
              scalar_to_phase_map_[scalar_id].species_type, scalar_id);
      }
    }

    //! get pre-factor needed for off-diagonal mesh-linearization of hist and source term
    double get_pre_factor_hist_and_source_od_mesh(const int scalar_id, const double fac,
        const double densnp, const double hist, const double rhsint)
    {
      // linearization of mesh motion: call base class with fac * rhsint
      double pre_factor = fac * rhsint;

      switch (scalar_to_phase_map_[scalar_id].species_type)
      {
        case Mat::ScaTraMatMultiPoro::SpeciesType::species_in_fluid:
        {
          // linearization of history: prefactor is densnp = porosity * rho * S
          // --> porosity has to be linearized with linearization of porosity
          //
          // dporosity/dd = dporosity/dJ * dJ/dd = dporosity/dJ * J * N_x
          // J denotes the determinant of the deformation gradient, i.e.,
          // det F = det ( d x / d X ) = det (dx/ds) * ( det(dX/ds) )^-1
          // in our case: pre_factor is scaled with density, i.e., porosity
          // --> scale it with 1.0/porosity here

          if (phase_manager_->porosity_depends_on_struct())
          {
            pre_factor += 1.0 * fac * hist * densnp / (phase_manager_->porosity()) *
                          phase_manager_->jacobian_def_grad() *
                          phase_manager_->porosity_deriv_wrt_jacobian_def_grad();
          }
          break;
        }
        case Mat::ScaTraMatMultiPoro::SpeciesType::species_in_volfrac:
        {
          // do nothing: correct pre_factor = fac * rhsint
          break;
        }
        case Mat::ScaTraMatMultiPoro::SpeciesType::species_in_solid:
        {
          if (phase_manager_->porosity_depends_on_struct())
          {
            pre_factor -= fac * hist * densnp /
                          (1 - phase_manager_->porosity() - phase_manager_->sum_add_vol_frac()) *
                          phase_manager_->jacobian_def_grad() *
                          phase_manager_->porosity_deriv_wrt_jacobian_def_grad();
          }
          break;
        }
        case Mat::ScaTraMatMultiPoro::SpeciesType::species_temperature:
        {
          // evaluate the effective heat capacity factor
          const int number_of_dofs =
              phase_manager_->num_fluid_phases() + phase_manager_->num_vol_frac();
          const int number_of_fluid_phases = phase_manager_->num_fluid_phases();

          const auto effective_heat_capacity = get_effective_heat_capacity();

          pre_factor += fac * hist * densnp / effective_heat_capacity *
                        phase_manager_->solid_density() * get_heat_capacity(number_of_dofs) *
                        (-1.0) * phase_manager_->jacobian_def_grad() *
                        phase_manager_->porosity_deriv_wrt_jacobian_def_grad();

          for (int phase_id = 0; phase_id < number_of_fluid_phases; ++phase_id)
          {
            pre_factor += fac * hist * densnp / effective_heat_capacity *
                          phase_manager_->density(phase_id) * get_heat_capacity(phase_id) *
                          phase_manager_->saturation(phase_id) *
                          phase_manager_->porosity_deriv_wrt_jacobian_def_grad() *
                          phase_manager_->jacobian_def_grad();
          }
          break;
        }

        default:
          FOUR_C_THROW("Unknown species type {} for species {}.",
              scalar_to_phase_map_[scalar_id].species_type, scalar_id);
      }

      return pre_factor;
    }

    //! get action
    [[nodiscard]] ScaTra::Action get_action() const { return action_; };

    //! return pressure gradient of current phase
    const Core::LinAlg::Matrix<nsd, 1>& get_pressure_gradient(const int current_phase) const
    {
      return pressure_gradient_[current_phase];
    }

    //! return absolute value of pressure gradient of current phase
    double get_pressure_gradient_norm(const int current_phase) const
    {
      return pressure_gradient_norm_[current_phase];
    }

    //! return porofluid material
    std::shared_ptr<Mat::FluidPoroMultiPhase> porofluid_material()
    {
      if (!material_is_set_)
        FOUR_C_THROW(
            "Pressure-based porofluid material has not yet been set in the variable manager.");

      return porofluid_material_;
    }

    //! Setup phase manager and variable manager of fluid
    void setup_porofluid_managers(const Core::Elements::Element* element,
        const Core::FE::Discretization& discretization, const int number_of_fluid_phases,
        const int total_number_of_porofluid_dofs)
    {
      PoroFluidMultiPhaseEleParameter* parameters =
          PoroFluidMultiPhaseEleParameter::instance(discretization.name());

      phase_manager_ = PoroFluidManager::PhaseManagerInterface::create_phase_manager(*parameters,
          nsd, porofluid_material()->material_type(),
          PoroPressureBased::Action::get_access_from_scatra, total_number_of_porofluid_dofs,
          number_of_fluid_phases);

      // access from outside to the phase manager:
      // scatra-discretization has fluid-discretization on dofset 2
      phase_manager_->setup(element, nds_scatra_porofluid_);

      variable_manager_ =
          PoroFluidManager::VariableManagerInterface<nsd, nen>::create_variable_manager(*parameters,
              PoroPressureBased::Action::get_access_from_scatra, porofluid_material(),
              total_number_of_porofluid_dofs, number_of_fluid_phases);
    }

    // extract the element and node values of the poro-fluid
    // --> extract them from its variable manager
    void extract_element_and_node_values_of_porofluid(Core::Elements::Element* element,
        Core::FE::Discretization& discretization, Core::Elements::LocationArray& location_array,
        Core::LinAlg::Matrix<nsd, nen>& node_coordinates)
    {
      // access from outside to the variable manager:
      // scatra-discretization has fluid-discretization on dofset 2
      variable_manager_->extract_element_and_node_values(
          *element, discretization, location_array, node_coordinates, nds_scatra_porofluid_);
    }

    // get the phase manager of the fluid
    std::shared_ptr<PoroFluidManager::PhaseManagerInterface> porofluid_phase_manager()
    {
      return phase_manager_;
    }

   private:
    //! phase pressure
    std::vector<double> pressure_;
    //! phase saturation
    std::vector<double> saturation_;
    //! phase density
    std::vector<double> density_;
    //! solid pressure
    double solid_pressure_;
    //! pressure gradient
    std::vector<Core::LinAlg::Matrix<nsd, 1>> pressure_gradient_;
    //! diffusion tensor of porofluid
    std::vector<Core::LinAlg::Matrix<nsd, nsd>> diffusion_tensor_porofluid_;
    //! norm of pressure-gradient
    std::vector<double> pressure_gradient_norm_;
    //! volume fraction
    std::vector<double> volume_fractions_;
    //! volume fraction pressure
    std::vector<double> volume_fractions_pressure_;

    //! pressure-based porofluid material
    std::shared_ptr<Mat::FluidPoroMultiPhase> porofluid_material_;

    //! delta for effective diffusivity
    std::vector<double> delta_;

    //! function IDs of relative mobility functions
    std::vector<int> relative_mobility_funct_id_;

    //! heat capacity: order [ <fluid>  <volfrac>  <solid> ]
    std::vector<double> heat_capacity_;
    //! thermal diffusivity: order [ <fluid>  <volfrac>  <solid> ]
    std::vector<double> thermal_diffusivity_;
    //! effective heat capacity
    double effective_heat_capacity_;

    //! minimum saturation under which the corresponding mass fraction is assumed to be zero
    std::vector<double> minimal_value_of_phase_;
    //! flag to check if we have to evaluate the species equation in this element
    std::vector<bool> evaluate_scalar_;
    //! scalar to phase map
    std::vector<Mat::ScaTraMatMultiPoro::ScalarToPhaseMap> scalar_to_phase_map_;
    //! check if multiphase material has been set
    bool material_is_set_;

    ScaTra::Action action_;

    //! phase manager of the fluid
    std::shared_ptr<PoroFluidManager::PhaseManagerInterface> phase_manager_;

    //! variable manager of the fluid
    std::shared_ptr<PoroFluidManager::VariableManagerInterface<nsd, nen>> variable_manager_;

    //! dofset of fluid field on scatra discretization
    // TODO: find a better way to do this
    const int nds_scatra_porofluid_ = 2;
  };

}  // namespace Discret::Elements


FOUR_C_NAMESPACE_CLOSE

#endif
