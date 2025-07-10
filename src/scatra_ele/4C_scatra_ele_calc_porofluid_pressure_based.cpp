// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_scatra_ele_calc_porofluid_pressure_based.hpp"

#include "4C_fem_discretization.hpp"
#include "4C_fem_general_extract_values.hpp"
#include "4C_mat_list.hpp"
#include "4C_mat_list_reactions.hpp"
#include "4C_mat_scatra.hpp"
#include "4C_mat_scatra_multiporo.hpp"
#include "4C_mat_structporo.hpp"
#include "4C_scatra_ele_parameter_std.hpp"
#include "4C_scatra_ele_parameter_timint.hpp"
#include "4C_utils_singleton_owner.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <Core::FE::CellType discretization_type>
Discret::Elements::ScaTraEleCalcPorofluidPressureBased<
    discretization_type>::ScaTraEleCalcPorofluidPressureBased(const int number_of_dofs_per_node,
    const int number_of_scalars, const std::string& discretization_name)
    : ScaTraEleCalc<discretization_type>::ScaTraEleCalc(
          number_of_dofs_per_node, number_of_scalars, discretization_name),
      ScaTraEleCalcPoro<discretization_type>::ScaTraEleCalcPoro(
          number_of_dofs_per_node, number_of_scalars, discretization_name),
      ScaTraEleCalcAdvReac<discretization_type>::ScaTraEleCalcAdvReac(
          number_of_dofs_per_node, number_of_scalars, discretization_name),
      ScaTraEleCalcPoroReac<discretization_type>::ScaTraEleCalcPoroReac(
          number_of_dofs_per_node, number_of_scalars, discretization_name),
      element_flux_np_(0)
{
  // replace internal variable manager by internal variable manager for pressure-based porofluid
  scatra_ele_calc::scatravarmanager_ =
      std::make_shared<ScaTraEleInternalVariableManagerPorofluidPressureBased<nsd_, nen_>>(
          scatra_ele_calc::numscal_);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <Core::FE::CellType discretization_type>
Discret::Elements::ScaTraEleCalcPorofluidPressureBased<discretization_type>*
Discret::Elements::ScaTraEleCalcPorofluidPressureBased<discretization_type>::instance(
    const int number_of_dofs_per_node, const int number_of_scalars,
    const std::string& discretization_name)
{
  static auto singleton_map = Core::Utils::make_singleton_map<std::string>(
      [number_of_dofs_per_node, number_of_scalars, &discretization_name]()
      {
        return std::unique_ptr<ScaTraEleCalcPorofluidPressureBased>(
            new ScaTraEleCalcPorofluidPressureBased(
                number_of_dofs_per_node, number_of_scalars, discretization_name));
      });

  return singleton_map[discretization_name].instance(Core::Utils::SingletonAction::create);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <Core::FE::CellType discretization_type>
int Discret::Elements::ScaTraEleCalcPorofluidPressureBased<discretization_type>::setup_calc(
    Core::Elements::Element* element, Core::FE::Discretization& discretization)
{
  poro_reaction::setup_calc(element, discretization);

  // get the material
  std::shared_ptr<Core::Mat::Material> material = element->material();

  // set the porofluid material in the element
  variable_manager()->set_porofluid_material(element);

  if (material->material_type() == Core::Materials::m_matlist or
      material->material_type() == Core::Materials::m_matlist_reactions)
  {
    const std::shared_ptr<const Mat::MatList>& material_list =
        std::dynamic_pointer_cast<const Mat::MatList>(material);
    if (material_list->num_mat() < scatra_ele_calc::numdofpernode_)
      FOUR_C_THROW("Not enough materials in MatList.");

    for (int scalar_id = 0; scalar_id < scatra_ele_calc::numdofpernode_; ++scalar_id)
    {
      int material_id = material_list->mat_id(scalar_id);
      std::shared_ptr<Core::Mat::Material> single_material =
          material_list->material_by_id(material_id);

      switch (single_material->material_type())
      {
        case Core::Materials::m_scatra_in_fluid_porofluid_pressure_based:
        {
          const std::shared_ptr<const Mat::ScatraMatMultiPoroFluid>& scatra_material =
              std::dynamic_pointer_cast<const Mat::ScatraMatMultiPoroFluid>(single_material);

          // smaller zero or greater equal to the number of fluid phases
          if (scatra_material->phase_id() < 0 or
              scatra_material->phase_id() >=
                  variable_manager()->porofluid_material()->num_fluid_phases())
          {
            FOUR_C_THROW(
                "Invalid phase ID {} for scalar {} (species in fluid = MAT_scatra_multiporo_fluid)",
                scatra_material->phase_id(), scalar_id);
          }

          const int porofluid_material_id =
              variable_manager()->porofluid_material()->mat_id(scatra_material->phase_id());
          std::shared_ptr<Core::Mat::Material> porofluid_material =
              variable_manager()->porofluid_material()->material_by_id(porofluid_material_id);

          if (porofluid_material->material_type() != Core::Materials::m_fluidporo_singlephase)
          {
            FOUR_C_THROW(
                "Invalid phase ID for scalar {} (species in fluid = MAT_scatra_multiporo_fluid)",
                scalar_id);
          }

          variable_manager()->set_phase_id_and_species_type(scalar_id, scatra_material->phase_id(),
              Mat::ScaTraMatMultiPoro::SpeciesType::species_in_fluid);
          // set delta in the variablemanager
          variable_manager()->set_delta(scatra_material->delta(), scalar_id);
          // set minimum saturation in the variablemanager
          variable_manager()->set_minimal_value_of_phase(scatra_material->min_sat(), scalar_id);
          // set reacts to external force
          variable_manager()->set_reacts_to_force(
              scatra_material->reacts_to_external_force(), scalar_id);
          // set relative mobility function ID
          variable_manager()->set_relative_mobility_function_id(
              scatra_material->relative_mobility_funct_id(), scalar_id);
          break;
        }
        case Core::Materials::m_scatra_in_volfrac_porofluid_pressure_based:
        {
          const std::shared_ptr<const Mat::ScatraMatMultiPoroVolFrac>& scatra_material =
              std::dynamic_pointer_cast<const Mat::ScatraMatMultiPoroVolFrac>(single_material);

          // smaller zero or greater equal the number of fluid phases + number of volume fractions
          if (scatra_material->phase_id() <
                  variable_manager()->porofluid_material()->num_fluid_phases() or
              scatra_material->phase_id() >=
                  variable_manager()->porofluid_material()->num_fluid_phases() +
                      variable_manager()->porofluid_material()->num_vol_frac())
          {
            FOUR_C_THROW(
                "Invalid phase ID {} for scalar {} (species in volume fraction = "
                "MAT_scatra_multiporo_volfrac)",
                scatra_material->phase_id(), scalar_id);
          }

          const int porofluid_material_id =
              variable_manager()->porofluid_material()->mat_id(scatra_material->phase_id());
          std::shared_ptr<Core::Mat::Material> porofluid_material =
              variable_manager()->porofluid_material()->material_by_id(porofluid_material_id);

          if (porofluid_material->material_type() != Core::Materials::m_fluidporo_singlevolfrac)
          {
            FOUR_C_THROW(
                "Invalid phase ID for scalar {} (species in volume fraction = "
                "MAT_scatra_multiporo_volfrac)",
                scalar_id);
          }

          variable_manager()->set_phase_id_and_species_type(scalar_id, scatra_material->phase_id(),
              Mat::ScaTraMatMultiPoro::SpeciesType::species_in_volfrac);
          // set delta in the variablemanager
          variable_manager()->set_delta(scatra_material->delta(), scalar_id);
          // set reacts to external force
          variable_manager()->set_reacts_to_force(
              scatra_material->reacts_to_external_force(), scalar_id);
          // set relative mobility function ID
          variable_manager()->set_relative_mobility_function_id(
              scatra_material->relative_mobility_funct_id(), scalar_id);

          break;
        }
        case Core::Materials::m_scatra_in_solid_porofluid_pressure_based:
        {
          const std::shared_ptr<const Mat::ScatraMatMultiPoroSolid>& scatra_material =
              std::dynamic_pointer_cast<const Mat::ScatraMatMultiPoroSolid>(single_material);

          // set delta in the variablemanager
          variable_manager()->set_delta(scatra_material->delta(), scalar_id);

          // dummy value -1000 for phaseID because species in solid do not have a phaseID
          variable_manager()->set_phase_id_and_species_type(
              scalar_id, -1000, Mat::ScaTraMatMultiPoro::SpeciesType::species_in_solid);

          break;
        }
        case Core::Materials::m_scatra_as_temperature_porofluid_pressure_based:
        {
          const std::shared_ptr<const Mat::ScatraMatMultiPoroTemperature>& scatra_material =
              std::dynamic_pointer_cast<const Mat::ScatraMatMultiPoroTemperature>(single_material);

          // assemble heat capacities of fluid phases, volume fractions and solid phase
          // cp order [ <fluid>  <volfrac>  <solid> ]
          std::vector<double> heat_capacity;
          std::vector heat_capacity_fluid(scatra_material->cp_fluid());
          std::vector heat_capacity_volfrac(scatra_material->cp_volfrac());

          heat_capacity.insert(
              heat_capacity.begin(), heat_capacity_fluid.begin(), heat_capacity_fluid.end());
          heat_capacity.insert(
              heat_capacity.end(), heat_capacity_volfrac.begin(), heat_capacity_volfrac.end());
          heat_capacity.insert(heat_capacity.end(), scatra_material->cp_solid());

          variable_manager()->set_heat_capacity(heat_capacity);

          // assemble thermal diffusivity of fluid phases, volume fractions and solid phase
          // kappa order [ <fluid>  <volfrac>  <solid> ]
          std::vector<double> thermal_diffusivity;
          std::vector thermal_diffusivity_fluid(scatra_material->kappa_fluid());
          std::vector thermal_diffusivity_volfrac(scatra_material->kappa_volfrac());

          thermal_diffusivity.insert(thermal_diffusivity.begin(), thermal_diffusivity_fluid.begin(),
              thermal_diffusivity_fluid.end());
          thermal_diffusivity.insert(thermal_diffusivity.end(), thermal_diffusivity_volfrac.begin(),
              thermal_diffusivity_volfrac.end());
          thermal_diffusivity.insert(thermal_diffusivity.end(), scatra_material->kappa_solid());

          variable_manager()->set_thermal_diffusivity(thermal_diffusivity);

          // dummy value -1000 for phaseID because temperature does not have a phaseID
          variable_manager()->set_phase_id_and_species_type(
              scalar_id, -1000, Mat::ScaTraMatMultiPoro::SpeciesType::species_temperature);

          break;
        }
        default:
          FOUR_C_THROW(
              "Material type {} is not supported for multiphase flow through porous media!",
              single_material->material_type());
      }
    }
  }
  else
  {
    switch (material->material_type())
    {
      case Core::Materials::m_scatra_in_fluid_porofluid_pressure_based:
      {
        const std::shared_ptr<const Mat::ScatraMatMultiPoroFluid>& scatra_material =
            std::dynamic_pointer_cast<const Mat::ScatraMatMultiPoroFluid>(material);

        // smaller zero or greater equal the number of fluid phases
        if (scatra_material->phase_id() < 0 or
            scatra_material->phase_id() >=
                variable_manager()->porofluid_material()->num_fluid_phases())
        {
          FOUR_C_THROW(
              "Invalid phase ID {} for scalar {} (species in fluid = MAT_scatra_multiporo_fluid)",
              scatra_material->phase_id(), 0);
        }

        const int porofluid_material_id =
            variable_manager()->porofluid_material()->mat_id(scatra_material->phase_id());
        std::shared_ptr<Core::Mat::Material> porofluid_material =
            variable_manager()->porofluid_material()->material_by_id(porofluid_material_id);

        if (porofluid_material->material_type() != Core::Materials::m_fluidporo_singlephase)
          FOUR_C_THROW(
              "Invalid phase ID for scalar {} (species in fluid = MAT_scatra_multiporo_fluid)", 0);

        variable_manager()->set_phase_id_and_species_type(
            0, scatra_material->phase_id(), Mat::ScaTraMatMultiPoro::SpeciesType::species_in_fluid);
        // set delta in the variablemanager
        variable_manager()->set_delta(scatra_material->delta(), 0);
        // set minimum saturation in the variablemanager
        variable_manager()->set_minimal_value_of_phase(scatra_material->min_sat(), 0);
        // set reacts to external force
        variable_manager()->set_reacts_to_force(scatra_material->reacts_to_external_force(), 0);
        // set relative mobility function ID
        variable_manager()->set_relative_mobility_function_id(
            scatra_material->relative_mobility_funct_id(), 0);
        break;
      }
      case Core::Materials::m_scatra_in_volfrac_porofluid_pressure_based:
      {
        const std::shared_ptr<const Mat::ScatraMatMultiPoroVolFrac>& poromat =
            std::dynamic_pointer_cast<const Mat::ScatraMatMultiPoroVolFrac>(material);

        // smaller zero or greater equal the number of fluid phases + number of volume fractions
        if (poromat->phase_id() < 0 or
            poromat->phase_id() >= variable_manager()->porofluid_material()->num_fluid_phases() +
                                       variable_manager()->porofluid_material()->num_vol_frac())
        {
          FOUR_C_THROW(
              "Invalid phase ID {} for scalar {} (species in volume fraction = "
              "MAT_scatra_multiporo_volfrac)",
              poromat->phase_id(), 0);
        }

        const int porofluid_material_id =
            variable_manager()->porofluid_material()->mat_id(poromat->phase_id());
        std::shared_ptr<Core::Mat::Material> porofluid_material =
            variable_manager()->porofluid_material()->material_by_id(porofluid_material_id);

        if (porofluid_material->material_type() != Core::Materials::m_fluidporo_singlevolfrac)
        {
          FOUR_C_THROW(
              "Invalid phase ID for scalar {} (species in volume fraction = "
              "MAT_scatra_multiporo_volfrac)",
              0);
        }

        variable_manager()->set_phase_id_and_species_type(
            0, poromat->phase_id(), Mat::ScaTraMatMultiPoro::SpeciesType::species_in_volfrac);
        // set delta in the variablemanager
        variable_manager()->set_delta(poromat->delta(), 0);
        // set reacts to external force
        variable_manager()->set_reacts_to_force(poromat->reacts_to_external_force(), 0);
        // set relative mobility function ID
        variable_manager()->set_relative_mobility_function_id(
            poromat->relative_mobility_funct_id(), 0);

        break;
      }
      case Core::Materials::m_scatra_in_solid_porofluid_pressure_based:
      {
        const std::shared_ptr<const Mat::ScatraMatMultiPoroSolid>& scatra_material =
            std::dynamic_pointer_cast<const Mat::ScatraMatMultiPoroSolid>(material);

        // set delta in the variablemanager
        variable_manager()->set_delta(scatra_material->delta(), 0);

        // dummy value -1000 for phaseID because species in solid do not have a phaseID
        variable_manager()->set_phase_id_and_species_type(
            0, -1000, Mat::ScaTraMatMultiPoro::SpeciesType::species_in_solid);

        break;
      }
      case Core::Materials::m_scatra_as_temperature_porofluid_pressure_based:
      {
        const std::shared_ptr<const Mat::ScatraMatMultiPoroTemperature>& scatra_material =
            std::dynamic_pointer_cast<const Mat::ScatraMatMultiPoroTemperature>(material);

        // assemble heat capacities of fluid phases, volume fractions and solid phase
        std::vector<double> heat_capacity;
        std::vector heat_capacity_fluid(scatra_material->cp_fluid());
        std::vector heat_capacity_volfrac(scatra_material->cp_volfrac());

        heat_capacity.insert(
            heat_capacity.begin(), heat_capacity_fluid.begin(), heat_capacity_fluid.end());
        heat_capacity.insert(
            heat_capacity.end(), heat_capacity_volfrac.begin(), heat_capacity_volfrac.end());
        heat_capacity.insert(heat_capacity.end(), scatra_material->cp_solid());

        variable_manager()->set_heat_capacity(heat_capacity);

        // assemble thermal diffusivity of fluid phases, volume fractions and solid phase
        std::vector<double> thermal_diffusivity;
        std::vector thermal_diffusivity_fluid(scatra_material->kappa_fluid());
        std::vector thermal_diffusivity_volfrac(scatra_material->kappa_volfrac());

        thermal_diffusivity.insert(thermal_diffusivity.begin(), thermal_diffusivity_fluid.begin(),
            thermal_diffusivity_fluid.end());
        thermal_diffusivity.insert(thermal_diffusivity.end(), thermal_diffusivity_volfrac.begin(),
            thermal_diffusivity_volfrac.end());
        thermal_diffusivity.insert(thermal_diffusivity.end(), scatra_material->kappa_solid());

        variable_manager()->set_thermal_diffusivity(thermal_diffusivity);

        // dummy value -1000 for phaseID because temperature does not have a phaseID
        variable_manager()->set_phase_id_and_species_type(
            0, -1000, Mat::ScaTraMatMultiPoro::SpeciesType::species_temperature);

        break;
      }
      default:
        FOUR_C_THROW("Material type {} is not supported for multiphase flow through porous media!",
            material->material_type());
    }
  }

  return 0;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <Core::FE::CellType discretization_type>
void Discret::Elements::ScaTraEleCalcPorofluidPressureBased<
    discretization_type>::extract_element_and_node_values(Core::Elements::Element* element,
    Teuchos::ParameterList& parameters, Core::FE::Discretization& discretization,
    Core::Elements::LocationArray& location_array)
{
  // extract action parameter
  const auto action = Teuchos::getIntegralValue<ScaTra::Action>(parameters, "action");
  variable_manager()->set_action(action);

  //---------------------------------------------------------------------------------------------
  //                                 STRUCTURE
  //---------------------------------------------------------------------------------------------

  // get additional state vector for ALE case: grid displacement
  if (scatra_ele_calc::scatrapara_->is_ale())
  {
    // get number of dofset associated with displacement related dofs
    const int nds_displacement = scatra_ele_calc::scatrapara_->nds_disp();

    const std::shared_ptr<const Core::LinAlg::Vector<double>> displacement_np =
        discretization.get_state(nds_displacement, "dispnp");
    if (displacement_np == nullptr) FOUR_C_THROW("Cannot get state vector 'dispnp'");

    // determine number of displacement related dofs per node
    const int number_of_displacement_dofs_per_node =
        location_array[nds_displacement].lm_.size() / nen_;

    // construct location vector for displacement related dofs
    std::vector<int> location_vector(nsd_ * nen_, -1);
    for (unsigned inode = 0; inode < nen_; ++inode)
    {
      for (unsigned idim = 0; idim < nsd_; ++idim)
      {
        location_vector[inode * nsd_ + idim] =
            location_array[nds_displacement]
                .lm_[inode * number_of_displacement_dofs_per_node + idim];
      }
    }

    // extract local values of displacement field from global state vector
    Core::FE::extract_my_values<Core::LinAlg::Matrix<nsd_, nen_>>(
        *displacement_np, scatra_ele_calc::edispnp_, location_vector);

    // add nodal displacements to point coordinates
    scatra_ele_calc::update_node_coordinates();
  }
  else
  {
    scatra_ele_calc::edispnp_.clear();
  }

  //---------------------------------------------------------------------------------------------
  //                                 SCATRA
  //---------------------------------------------------------------------------------------------

  // extract local values from the global vectors
  const std::shared_ptr<const Core::LinAlg::Vector<double>> hist = discretization.get_state("hist");
  const std::shared_ptr<const Core::LinAlg::Vector<double>> phi_np =
      discretization.get_state("phinp");
  if (hist == nullptr || phi_np == nullptr)
    FOUR_C_THROW("Cannot get state vector 'hist' and/or 'phinp'");

  // values of scatra field are always in first dofset
  const std::vector<int>& global_dof_ids = location_array[0].lm_;
  Core::FE::extract_my_values<Core::LinAlg::Matrix<nen_, 1>>(
      *hist, scatra_ele_calc::ehist_, global_dof_ids);
  Core::FE::extract_my_values<Core::LinAlg::Matrix<nen_, 1>>(
      *phi_np, scatra_ele_calc::ephinp_, global_dof_ids);

  if (scatra_ele_calc::scatraparatimint_->is_gen_alpha() and
      not scatra_ele_calc::scatraparatimint_->is_incremental())
  {
    // extract additional local values from global vector
    const std::shared_ptr<const Core::LinAlg::Vector<double>> phi_n =
        discretization.get_state("phin");
    if (phi_n == nullptr) FOUR_C_THROW("Cannot get state vector 'phin'");
    Core::FE::extract_my_values<Core::LinAlg::Matrix<nen_, 1>>(
        *phi_n, scatra_ele_calc::ephin_, global_dof_ids);
  }

  if (scatra_ele_calc::scatrapara_->has_external_force())
  {
    // get number of dofset associated with velocity related dofs
    const auto nds_velocity = scatra_ele_calc::scatrapara_->nds_vel();
    const auto force_velocity = discretization.get_state(nds_velocity, "force_velocity");

    const int number_of_dofs_per_node =
        location_array[nds_velocity].lm_.size() / scatra_ele_calc::nen_;
    std::vector<int> location_vector(scatra_ele_calc::nsd_ * scatra_ele_calc::nen_, -1);
    for (unsigned inode = 0; inode < scatra_ele_calc::nen_; ++inode)
    {
      for (unsigned idim = 0; idim < scatra_ele_calc::nsd_; ++idim)
        location_vector[inode * scatra_ele_calc::nsd_ + idim] =
            location_array[nds_velocity].lm_[inode * number_of_dofs_per_node + idim];
    }

    Core::FE::extract_my_values<Core::LinAlg::Matrix<scatra_ele_calc::nsd_, scatra_ele_calc::nen_>>(
        *force_velocity, scatra_ele_calc::eforcevelocity_, location_vector);
  }

  //---------------------------------------------------------------------------------------------
  //                                 FLUID
  //---------------------------------------------------------------------------------------------

  // get number of dofset associated with pressure/fluid related dofs
  const int nds_pressure = scatra_ele_calc::scatrapara_->nds_pres();

  // determine number of velocity related dofs per node (= number of phases)
  const int number_of_fluid_phases = variable_manager()->porofluid_material()->num_fluid_phases();
  const int total_number_of_porofluid_dofs = variable_manager()->porofluid_material()->num_mat();

  // extract element and node values of the porofluid
  if (discretization.has_state(nds_pressure, "phinp_fluid"))
  {
    variable_manager()->setup_porofluid_managers(
        element, discretization, number_of_fluid_phases, total_number_of_porofluid_dofs);
    variable_manager()->extract_element_and_node_values_of_porofluid(
        element, discretization, location_array, scatra_ele_calc::xyze_);
    enable_L2_projection_ = parameters.get<bool>("L2-projection");
    // extract the nodal flux
    if (enable_L2_projection_)
    {
      extract_nodal_flux(discretization, location_array, number_of_fluid_phases);
    }
  }
  else
  {
    FOUR_C_THROW("Something went wrong here, scatra-dis does not have fluid primary variable");
  }

  // ---------------------------------------------------------------------
  // call routine for calculation of body force in element nodes
  // (time n+alpha_F for generalized-alpha scheme, at time n+1 otherwise)
  // ---------------------------------------------------------------------
  scatra_ele_calc::body_force(element);
  //--------------------------------------------------------------------------------
  // further node-based source terms not given via Neumann volume condition
  // i.e., set special body force for homogeneous isotropic turbulence
  //--------------------------------------------------------------------------------
  scatra_ele_calc::other_node_based_source_terms(global_dof_ids, discretization, parameters);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <Core::FE::CellType discretization_type>
void Discret::Elements::ScaTraEleCalcPorofluidPressureBased<
    discretization_type>::extract_nodal_flux(const Core::FE::Discretization& discretization,
    Core::Elements::LocationArray& location_array, const int number_of_fluid_phases)
{
  // resize state vectors based on number of phases
  element_flux_np_.resize(number_of_fluid_phases);

  // get number of dofset associated with velocity related dofs
  const int nds_velocity = scatra_ele_calc::scatrapara_->nds_vel();

  const std::string state_prefix = "flux";
  for (int phase_id = 0; phase_id < number_of_fluid_phases; phase_id++)
  {
    std::stringstream state_name;
    state_name << state_prefix << phase_id;

    // get convective (velocity - mesh displacement) velocity at nodes
    std::shared_ptr<const Core::LinAlg::Vector<double>> nodal_convective_velocity =
        discretization.get_state(nds_velocity, state_name.str());
    if (nodal_convective_velocity == nullptr)
      FOUR_C_THROW("Cannot get state vector {}", state_name.str());

    // extract local values of convective velocity field from global state vector
    Core::FE::extract_my_values<Core::LinAlg::Matrix<nsd_, nen_>>(
        *nodal_convective_velocity, element_flux_np_[phase_id], location_array[nds_velocity].lm_);
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <Core::FE::CellType discretization_type>
double
Discret::Elements::ScaTraEleCalcPorofluidPressureBased<discretization_type>::compute_pore_pressure()
{
  return variable_manager()->get_solid_pressure();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <Core::FE::CellType discretization_type>
void Discret::Elements::ScaTraEleCalcPorofluidPressureBased<discretization_type>::materials(
    const std::shared_ptr<const Core::Mat::Material> material, const int scalar_id, double& densn,
    double& densnp, double& densam, double& viscosity, const int gauss_point_id)
{
  switch (material->material_type())
  {
    case Core::Materials::m_scatra_in_fluid_porofluid_pressure_based:
    {
      material_scatra_in_fluid(material, scalar_id, densn, densnp, densam, gauss_point_id);
      break;
    }
    case Core::Materials::m_scatra_in_volfrac_porofluid_pressure_based:
    {
      material_scatra_in_volfrac(material, scalar_id, densn, densnp, densam, gauss_point_id);
      break;
    }
    case Core::Materials::m_scatra_in_solid_porofluid_pressure_based:
    {
      material_scatra_in_solid(material, scalar_id, densn, densnp, densam, gauss_point_id);
      break;
    }
    case Core::Materials::m_scatra_as_temperature_porofluid_pressure_based:
    {
      material_scatra_as_temperature(material, scalar_id, densn, densnp, densam, gauss_point_id);
      break;
    }
    default:
      FOUR_C_THROW("Material type {} is not supported for multiphase flow through porous media!",
          material->material_type());
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <Core::FE::CellType discretization_type>
void Discret::Elements::ScaTraEleCalcPorofluidPressureBased<discretization_type>::
    material_scatra_in_fluid(const std::shared_ptr<const Core::Mat::Material> material,
        const int scalar_id, double& densn, double& densnp, double& densam,
        const int gauss_point_id)
{
  if (gauss_point_id == -1)
    FOUR_C_THROW(
        "No gauss point given for evaluation of MatMultiPoro material. Check your input file.");

  const std::shared_ptr<const Mat::ScatraMatMultiPoroFluid>& scatra_material =
      std::dynamic_pointer_cast<const Mat::ScatraMatMultiPoroFluid>(material);

  // volume fraction of fluid phase: volfrac_fluid = porosity * saturation_fluid
  double volfrac_fluid = 0.0;
  // d_eff = d_0 * (porosity * saturation(scalar_id))^delta
  double d_eff = 0.0;

  if (variable_manager()->evaluate_scalar(scalar_id))
  {
    volfrac_fluid = variable_manager()->porofluid_phase_manager()->porosity() *
                    variable_manager()->get_saturation(scalar_id);
    d_eff = std::pow(variable_manager()->porofluid_phase_manager()->porosity() *
                         variable_manager()->get_saturation(scalar_id),
        scatra_material->delta());
  }

  {
    // set diffusivity (scaled with volfrac_fluid)
    poro::set_diffusivity(scatra_material, scalar_id, volfrac_fluid * d_eff);

    // set densities (scaled with volfrac_fluid)
    poro::set_densities(volfrac_fluid, densn, densnp, densam);
  }
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <Core::FE::CellType discretization_type>
void Discret::Elements::ScaTraEleCalcPorofluidPressureBased<discretization_type>::
    material_scatra_in_volfrac(const std::shared_ptr<const Core::Mat::Material> material,
        const int scalar_id, double& densn, double& densnp, double& densam,
        const int gauss_point_id)
{
  if (gauss_point_id == -1)
    FOUR_C_THROW(
        "No gauss point given for evaluation of MatMultiPoro material. Check your input file.");

  const std::shared_ptr<const Mat::ScatraMatMultiPoroVolFrac>& scatra_material =
      std::dynamic_pointer_cast<const Mat::ScatraMatMultiPoroVolFrac>(material);

  // volume fraction
  double volfrac = 0.0;
  // d_eff = d_0 * (porosity * saturation(scalar_id))^delta
  double d_eff = 0.0;

  if (variable_manager()->evaluate_scalar(scalar_id))
  {
    volfrac = variable_manager()->get_volume_fraction(scalar_id);
    d_eff = std::pow(variable_manager()->get_volume_fraction(scalar_id), scatra_material->delta());
  }

  {
    // set diffusivity (scaled with volfrac)
    poro::set_diffusivity(scatra_material, scalar_id, volfrac * d_eff);

    // set densities (scaled with volfrac)
    poro::set_densities(volfrac, densn, densnp, densam);
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <Core::FE::CellType discretization_type>
void Discret::Elements::ScaTraEleCalcPorofluidPressureBased<discretization_type>::
    material_scatra_in_solid(const std::shared_ptr<const Core::Mat::Material> material,
        const int scalar_id, double& densn, double& densnp, double& densam,
        const int gauss_point_id)
{
  if (gauss_point_id == -1)
    FOUR_C_THROW(
        "No gauss point given for evaluation of MatMultiPoro material. Check your input file.");

  const std::shared_ptr<const Mat::ScatraMatMultiPoroSolid>& scatra_material =
      std::dynamic_pointer_cast<const Mat::ScatraMatMultiPoroSolid>(material);

  // volume fraction of solid phase: volfrac_solid_phase = (1 - porosity - sumaddvolfrac)
  double volfrac_solid_phase = 0.0;
  // d_eff = d_0 * (porosity * saturation(scalar_id))^delta
  double d_eff = 0.0;

  if (variable_manager()->evaluate_scalar(scalar_id))
  {
    volfrac_solid_phase = (1 - variable_manager()->porofluid_phase_manager()->porosity() -
                           variable_manager()->porofluid_phase_manager()->sum_add_vol_frac());
    d_eff = std::pow((1 - variable_manager()->porofluid_phase_manager()->porosity() -
                         variable_manager()->porofluid_phase_manager()->sum_add_vol_frac()),
        scatra_material->delta());
  }

  {
    // set diffusivity (scaled with volfrac_solid_phase)
    poro::set_diffusivity(scatra_material, scalar_id, volfrac_solid_phase * d_eff);

    // set densities (scaled with volfrac_solid_phase)
    poro::set_densities(volfrac_solid_phase, densn, densnp, densam);
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <Core::FE::CellType discretization_type>
void Discret::Elements::ScaTraEleCalcPorofluidPressureBased<discretization_type>::
    material_scatra_as_temperature(const std::shared_ptr<const Core::Mat::Material> material,
        const int scalar_id, double& densn, double& densnp, double& densam,
        const int gauss_point_id)
{
  if (gauss_point_id == -1)
    FOUR_C_THROW(
        "No gauss point given for evaluation of MatMultiPoro material. Check your input file.");

  const std::shared_ptr<const Mat::ScatraMatMultiPoroTemperature>& scatra_material =
      std::dynamic_pointer_cast<const Mat::ScatraMatMultiPoroTemperature>(material);

  double effective_heat_capacity = 0.0;
  double effective_thermal_conductivity = 0.0;

  if (variable_manager()->evaluate_scalar(scalar_id))
  {
    const int number_of_fluid_phases =
        variable_manager()->porofluid_phase_manager()->num_fluid_phases();
    const int number_of_volume_fractions =
        variable_manager()->porofluid_phase_manager()->num_vol_frac();

    effective_heat_capacity =
        (1 - variable_manager()->porofluid_phase_manager()->porosity() -
            variable_manager()->porofluid_phase_manager()->sum_add_vol_frac()) *
        variable_manager()->porofluid_phase_manager()->solid_density() *
        scatra_material->cp_solid();

    // kappa_eff = kappa_s*poro_s + kappa_fluids*poro*saturation_fluids + kappa_volfrac*poro_volfrac
    effective_thermal_conductivity =
        (1 - variable_manager()->porofluid_phase_manager()->porosity() -
            variable_manager()->porofluid_phase_manager()->sum_add_vol_frac()) *
        scatra_material->kappa_solid();

    for (int phase = 0; phase < number_of_fluid_phases; ++phase)
    {
      effective_heat_capacity += scatra_material->cp_fluid(phase) *
                                 variable_manager()->porofluid_phase_manager()->porosity() *
                                 variable_manager()->porofluid_phase_manager()->saturation(phase) *
                                 variable_manager()->porofluid_phase_manager()->density(phase);

      effective_thermal_conductivity +=
          scatra_material->kappa_fluid(phase) *
          variable_manager()->porofluid_phase_manager()->porosity() *
          variable_manager()->porofluid_phase_manager()->saturation(phase);
    }

    for (int phase = 0; phase < number_of_volume_fractions; ++phase)
    {
      effective_heat_capacity +=
          scatra_material->cp_volfrac(phase) *
          variable_manager()->porofluid_phase_manager()->vol_frac(phase) *
          variable_manager()->porofluid_phase_manager()->vol_frac_density(phase);

      effective_thermal_conductivity +=
          scatra_material->kappa_volfrac(phase) *
          variable_manager()->porofluid_phase_manager()->vol_frac(phase);
    }
  }

  {
    variable_manager()->set_effective_heat_capacity(effective_heat_capacity);

    // set diffusivity
    poro::set_diffusivity(scatra_material, scalar_id, effective_thermal_conductivity);

    // set densities
    poro::set_densities(effective_heat_capacity, densn, densnp, densam);
  }
}

/*------------------------------------------------------------------------------*
 *------------------------------------------------------------------------------*/
template <Core::FE::CellType discretization_type>
void Discret::Elements::ScaTraEleCalcPorofluidPressureBased<
    discretization_type>::set_internal_variables_for_mat_and_rhs()
{
  variable_manager()->set_internal_variables_porofluid_pressure_based(scatra_ele_calc::funct_,
      scatra_ele_calc::derxy_, scatra_ele_calc::deriv_, scatra_ele_calc::xjm_,
      poro_reaction::xyze0_, scatra_ele_calc::ephinp_, scatra_ele_calc::ephin_,
      scatra_ele_calc::ehist_, scatra_ele_calc::eforcevelocity_);

  if (enable_L2_projection_)
  {
    variable_manager()->adapt_convective_term_for_l2(
        scatra_ele_calc::funct_, scatra_ele_calc::derxy_, element_flux_np_);
  }
}

/*-------------------------------------------------------------------------------*
 *-------------------------------------------------------------------------------*/
template <Core::FE::CellType discretization_type>
void Discret::Elements::ScaTraEleCalcPorofluidPressureBased<
    discretization_type>::set_advanced_reaction_terms(const int scalar_id,
    const std::shared_ptr<Mat::MatListReactions> material_reaction_list,
    const double* gauss_point_coordinates)
{
  const std::shared_ptr<ScaTraEleReaManagerAdvReac> reaction_manager =
      advanced_reaction::rea_manager();

  fill_coupling_vector_and_add_variables(*material_reaction_list, *reaction_manager);

  const ScaTra::Action action = variable_manager()->get_action();

  auto time = scatra_ele_calc::scatraparatimint_->time();
  // note: we always need the reaction term to calculate rhsint, which is needed also for OD-terms
  reaction_manager->add_to_rea_body_force(
      material_reaction_list->calc_rea_body_force_term(scalar_id,
          scatra_ele_calc::scatravarmanager_->phinp(), coupling_values_, gauss_point_coordinates,
          time),
      scalar_id);

  std::vector<std::pair<std::string, double>> empty_constants;

  switch (action)
  {
    case ScaTra::Action::calc_mat_and_rhs:
    case ScaTra::Action::calc_initial_time_deriv:
    {
      material_reaction_list->calc_rea_body_force_deriv_matrix(scalar_id,
          reaction_manager->get_rea_body_force_deriv_vector(scalar_id),
          scatra_ele_calc::scatravarmanager_->phinp(), coupling_values_, gauss_point_coordinates,
          time);

      break;
    }
    case ScaTra::Action::calc_scatra_mono_odblock_fluid:
    {
      material_reaction_list->calc_rea_body_force_deriv_matrix_add_variables(scalar_id,
          reaction_manager->get_rea_body_force_deriv_vector_add_variables(scalar_id),
          scatra_ele_calc::scatravarmanager_->phinp(), coupling_values_, empty_constants,
          gauss_point_coordinates, time);

      break;
    }
    case ScaTra::Action::calc_scatra_mono_odblock_mesh:
    {
      if (variable_manager()->porofluid_phase_manager()->porosity_depends_on_struct())
      {
        material_reaction_list->calc_rea_body_force_deriv_matrix_add_variables(scalar_id,
            reaction_manager->get_rea_body_force_deriv_vector_add_variables(scalar_id),
            scatra_ele_calc::scatravarmanager_->phinp(), coupling_values_, empty_constants,
            gauss_point_coordinates, time);
      }
      break;
    }
    default:
      FOUR_C_THROW("Wrong action type in variable_manager(), action type is {}", action);
  }
}

/*-------------------------------------------------------------------------------*
 *-------------------------------------------------------------------------------*/
template <Core::FE::CellType discretization_type>
void Discret::Elements::ScaTraEleCalcPorofluidPressureBased<discretization_type>::get_rhs_int(
    double& rhsint, const double densnp, const int scalar_id)
{
  // only difference is the inverse scaling with density
  advanced_reaction::get_rhs_int(
      rhsint, 1.0 / variable_manager()->get_density(scalar_id), scalar_id);
}

/*------------------------------------------------------------------------------ *
 *-------------------------------------------------------------------------------*/
template <Core::FE::CellType discretization_type>
void Discret::Elements::ScaTraEleCalcPorofluidPressureBased<discretization_type>::calc_mat_react(
    Core::LinAlg::SerialDenseMatrix& element_matrix, const int scalar_id, const double timefacfac,
    const double timetaufac, const double taufac, const double densnp,
    const Core::LinAlg::Matrix<nen_, 1>& sgconv, const Core::LinAlg::Matrix<nen_, 1>& diff)
{
  // only difference is the inverse scaling with density
  advanced_reaction::calc_mat_react(element_matrix, scalar_id, timefacfac, timetaufac, taufac,
      1.0 / variable_manager()->get_density(scalar_id), sgconv, diff);
}

/*-------------------------------------------------------------------------------*
 *-------------------------------------------------------------------------------*/
template <Core::FE::CellType discretization_type>
void Discret::Elements::ScaTraEleCalcPorofluidPressureBased<discretization_type>::
    fill_coupling_vector_and_add_variables(Mat::MatListReactions& material_list_reactions,
        ScaTraEleReaManagerAdvReac& reaction_manager)
{
  // if it is empty rebuilt it
  if (coupling_values_.empty())
  {
    // pressures
    const std::vector<double>& pressures = variable_manager()->get_pressure();
    const int number_of_fluid_phases = variable_manager()->porofluid_material()->num_fluid_phases();
    for (int phase_id = 0; phase_id < number_of_fluid_phases; ++phase_id)
    {
      coupling_values_.emplace_back("p" + std::to_string(phase_id + 1), pressures[phase_id]);
    }

    // saturation
    const std::vector<double>& saturations = variable_manager()->get_saturation();
    for (int phase_id = 0; phase_id < number_of_fluid_phases; ++phase_id)
    {
      coupling_values_.emplace_back("S" + std::to_string(phase_id + 1), saturations[phase_id]);
    }
    // porosity
    coupling_values_.push_back(std::pair<std::string, double>(
        "porosity", variable_manager()->porofluid_phase_manager()->porosity()));

    // additional volume fractions
    const std::vector<double>& volume_fractions = variable_manager()->get_volume_fraction();
    const int number_of_volume_fractions =
        variable_manager()->porofluid_phase_manager()->num_vol_frac();
    for (int phase_id = 0; phase_id < number_of_volume_fractions; ++phase_id)
    {
      coupling_values_.emplace_back(
          "VF" + std::to_string(phase_id + 1), volume_fractions[phase_id]);
    }

    // additional volume fraction pressures
    const std::vector<double>& volume_fraction_pressures =
        variable_manager()->get_volume_fractions_pressure();
    for (int phase_id = 0; phase_id < number_of_volume_fractions; phase_id++)
    {
      coupling_values_.emplace_back(
          "VFP" + std::to_string(phase_id + 1), volume_fraction_pressures[phase_id]);
    }

    // initialize and add the variables to the reaction manager --> has to be done only once
    reaction_manager.initialize_rea_body_force_deriv_vector_add_variables(
        scatra_ele_calc::numdofpernode_, coupling_values_.size());
    // error will be thrown if reaction != by-reaction coupling is chosen
    for (int dof_id = 0; dof_id < scatra_ele_calc::numdofpernode_; dof_id++)
      material_list_reactions.add_additional_variables(dof_id, coupling_values_);
  }
  // directly copy values (rely on order for performance reasons)
  else
  {
    // pressures
    const std::vector<double>& pressures = variable_manager()->get_pressure();
    const int number_of_fluid_phases = variable_manager()->porofluid_material()->num_fluid_phases();
    for (int phase_id = 0; phase_id < number_of_fluid_phases; phase_id++)
    {
      coupling_values_[phase_id].second = pressures[phase_id];
    }
    // saturation
    const std::vector<double>& saturations = variable_manager()->get_saturation();
    for (int phase_id = 0; phase_id < number_of_fluid_phases; phase_id++)
    {
      coupling_values_[number_of_fluid_phases + phase_id].second = saturations[phase_id];
    }
    // porosity
    coupling_values_[2 * number_of_fluid_phases].second =
        variable_manager()->porofluid_phase_manager()->porosity();

    // additional volume fractions
    const std::vector<double>& volume_fractions = variable_manager()->get_volume_fraction();
    const std::vector<double>& volume_fraction_pressures =
        variable_manager()->get_volume_fractions_pressure();
    const int number_of_volume_fractions =
        variable_manager()->porofluid_phase_manager()->num_vol_frac();
    for (int i = 0; i < number_of_volume_fractions; i++)
    {
      coupling_values_[2 * number_of_fluid_phases + 1 + i].second = volume_fractions[i];
      coupling_values_[2 * number_of_fluid_phases + number_of_volume_fractions + 1 + i].second =
          volume_fraction_pressures[i];
    }
  }
}

/*-----------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------*/
template <Core::FE::CellType discretization_type>
void Discret::Elements::ScaTraEleCalcPorofluidPressureBased<discretization_type>::calc_mat_conv(
    Core::LinAlg::SerialDenseMatrix& element_matrix, const int scalar_id, const double timefacfac,
    const double densnp, const Core::LinAlg::Matrix<nen_, 1>& sgconv)
{
  // case of zero saturation/volfrac
  // no convective term for species in solid
  if (variable_manager()->evaluate_scalar(scalar_id) &&
      variable_manager()->get_species_type(scalar_id) !=
          Mat::ScaTraMatMultiPoro::SpeciesType::species_in_solid)
  {
    // the only difference to the base class version is that there is no scaling with the density
    poro_reaction::calc_mat_conv(element_matrix, scalar_id, timefacfac, 1.0, sgconv);
  }
}

/*-----------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------*/
template <Core::FE::CellType discretization_type>
void Discret::Elements::ScaTraEleCalcPorofluidPressureBased<discretization_type>::calc_mat_mass(
    Core::LinAlg::SerialDenseMatrix& element_matrix, const int& scalar_id, const double& fac,
    const double& densam)
{
  if (variable_manager()->evaluate_scalar(scalar_id))
  {
    // the only difference to the base class version is that there is no scaling with the density
    poro_reaction::calc_mat_mass(element_matrix, scalar_id, fac, densam);
  }
  else
  {
    if (variable_manager()->get_species_type(scalar_id) ==
        Mat::ScaTraMatMultiPoro::SpeciesType::species_in_fluid)
    {
      // If we have zero "densities" (porosity * saturation(scalar_id)), which mostly happens for
      // tumor cells, the whole equation will be equal to zero since it is scaled with the density.
      // In that case, the mass fraction of the species (necrotic tumor cells) also has to be zero.
      // --> Here, we explicitly force it to be zero through a "Dirichlet" boundary condition.
      for (unsigned vi = 0; vi < nen_; ++vi)
      {
        const int fvi = vi * scatra_ele_calc::numdofpernode_ + scalar_id;
        element_matrix(fvi, fvi) += penalty_;
      }
    }
  }
}


/*------------------------------------------------------------------------------------------*
 *------------------------------------------------------------------------------------------*/
template <Core::FE::CellType discretization_type>
void Discret::Elements::ScaTraEleCalcPorofluidPressureBased<
    discretization_type>::calc_mat_conv_add_cons(Core::LinAlg::SerialDenseMatrix& element_matrix,
    const int scalar_id, const double timefacfac, const double vdiv, const double densnp)
{
  // the only difference to the base class version is that there is no scaling with the density
  poro_reaction::calc_mat_conv_add_cons(element_matrix, scalar_id, timefacfac, vdiv, 1.0);
}

/*------------------------------------------------------------------- *
 *--------------------------------------------------------------------*/
template <Core::FE::CellType discretization_type>
void Discret::Elements::ScaTraEleCalcPorofluidPressureBased<
    discretization_type>::recompute_conv_phi_for_rhs(const int scalar_id,
    const Core::LinAlg::Matrix<nsd_, 1>& sgvelint, const double densnp, const double densn,
    const double vdiv)
{
  // the only difference to the base class version is that there is no scaling with the density
  poro_reaction::recompute_conv_phi_for_rhs(scalar_id, sgvelint, 1.0, 1.0, vdiv);
}

/*-------------------------------------------------------------------- *
 *---------------------------------------------------------------------*/
template <Core::FE::CellType discretization_type>
void Discret::Elements::ScaTraEleCalcPorofluidPressureBased<discretization_type>::calc_conv_od_mesh(
    Core::LinAlg::SerialDenseMatrix& element_matrix, const int scalar_id,
    const int number_of_dofs_per_node_mesh, const double fac, const double rhsfac,
    const double densnp, const double J, const Core::LinAlg::Matrix<nsd_, 1>& gradphi,
    const Core::LinAlg::Matrix<nsd_, 1>& convelint)
{
  // case of zero saturation/volfrac
  // no convective term for species in solid
  if (variable_manager()->evaluate_scalar(scalar_id) &&
      variable_manager()->get_species_type(scalar_id) !=
          Mat::ScaTraMatMultiPoro::SpeciesType::species_in_solid)
  {
    static Core::LinAlg::Matrix<nsd_, nsd_> diffusion_tensor;
    static Core::LinAlg::Matrix<nsd_, 1> pressure_gradient_reference;

    // linearization of mesh motion
    // dJ/dd = dJ/dF : dF/dd = J * F^-T . N_{\psi} = J * N_x
    // J denotes the determinant of the Jacobian of the mapping between current and parameter space,
    // i.e., det(dx/ds) in our case: rhsfac = J * dt * theta --> d(rhsfac)/dd = rhsfac * N_x
    for (unsigned vi = 0; vi < nen_; ++vi)
    {
      const int fvi = vi * scatra_ele_calc::numdofpernode_ + scalar_id;
      const double v =
          rhsfac * scatra_ele_calc::funct_(vi) * (-1.0) * variable_manager()->conv_phi(scalar_id);

      for (unsigned ui = 0; ui < nen_; ++ui)
      {
        for (unsigned idim = 0; idim < nsd_; ++idim)
        {
          const int fui = ui * nsd_ + idim;
          element_matrix(fvi, fui) += v * scatra_ele_calc::derxy_(idim, ui);
        }
      }
    }

    if (variable_manager()->get_species_type(scalar_id) ==
            Mat::ScaTraMatMultiPoro::SpeciesType::species_in_fluid ||
        variable_manager()->get_species_type(scalar_id) ==
            Mat::ScaTraMatMultiPoro::SpeciesType::species_in_volfrac)
    {
      variable_manager()->get_diffusion_tensor_fluid(
          scalar_id, diffusion_tensor, variable_manager()->get_phase_id(scalar_id));
      variable_manager()->get_pressure_gradient_reference(scatra_ele_calc::xjm_,
          pressure_gradient_reference, variable_manager()->get_phase_id(scalar_id));
      const double vrhs = rhsfac * 1.0 / J * diffusion_tensor(0, 0) * (-1.0);

      // linearization of pressure gradient
      // standard Galerkin terms  -- "shapederivatives" pressure gradient
      apply_shape_derivs_pressure_grad(
          element_matrix, scalar_id, vrhs, gradphi, pressure_gradient_reference);

      // linearization of gradphi
      // standard Galerkin terms  -- "shapederivatives" gradphi
      poro_reaction::apply_shape_derivs_conv(
          element_matrix, scalar_id, rhsfac, 1.0, J, gradphi, convelint);
    }

    else if (variable_manager()->get_species_type(scalar_id) ==
             Mat::ScaTraMatMultiPoro::SpeciesType::species_temperature)
    {
      const int number_of_fluid_phases =
          variable_manager()->porofluid_phase_manager()->num_fluid_phases();
      const int number_of_volume_fractions =
          variable_manager()->porofluid_phase_manager()->num_vol_frac();

      for (int phase = 0; phase < number_of_fluid_phases; ++phase)
      {
        variable_manager()->get_diffusion_tensor_fluid(scalar_id, diffusion_tensor, phase);
        variable_manager()->get_pressure_gradient_reference(
            scatra_ele_calc::xjm_, pressure_gradient_reference, phase);

        const double vrhs = rhsfac * 1.0 / J * diffusion_tensor(0, 0) * (-1.0) *
                            variable_manager()->porofluid_phase_manager()->density(phase) *
                            variable_manager()->get_heat_capacity(phase);

        // linearization of pressure gradient
        // standard Galerkin terms  -- "shapederivatives" pressure gradient
        apply_shape_derivs_pressure_grad(
            element_matrix, scalar_id, vrhs, gradphi, pressure_gradient_reference);
      }

      for (int phase_id = number_of_fluid_phases;
          phase_id < number_of_fluid_phases + number_of_volume_fractions; ++phase_id)
      {
        variable_manager()->get_diffusion_tensor_fluid(scalar_id, diffusion_tensor, phase_id);
        variable_manager()->get_pressure_gradient_reference(
            scatra_ele_calc::xjm_, pressure_gradient_reference, phase_id);

        const double vrhs = rhsfac * 1.0 / J * diffusion_tensor(0, 0) * (-1.0) *
                            variable_manager()->porofluid_phase_manager()->vol_frac_density(
                                phase_id - number_of_fluid_phases) *
                            variable_manager()->get_heat_capacity(phase_id);

        // linearization of pressure gradient
        // standard Galerkin terms  -- "shapederivatives" pressure gradient
        apply_shape_derivs_pressure_grad(
            element_matrix, scalar_id, vrhs, gradphi, pressure_gradient_reference);
      }

      // linearization of gradphi
      // standard Galerkin terms  -- "shapederivatives" gradphi
      poro_reaction::apply_shape_derivs_conv(
          element_matrix, scalar_id, rhsfac, 1.0, J, gradphi, convelint);
    }
    else
    {
      FOUR_C_THROW("Species type no valid!");
    }
  }
}

/*-------------------------------------------------------------------- *
 *---------------------------------------------------------------------*/
template <Core::FE::CellType discretization_type>
void Discret::Elements::ScaTraEleCalcPorofluidPressureBased<
    discretization_type>::calc_lin_mass_od_mesh(Core::LinAlg::SerialDenseMatrix& element_matrix,
    const int scalar_id, const int number_of_dofs_per_node_mesh, const double rhsfac,
    const double fac, const double densam, const double densnp, const double phinp,
    const double hist, const double J, const Core::LinAlg::Matrix<1, nsd_ * nen_>& dJ_dmesh)
{
  // case of zero saturation/volfrac
  if (variable_manager()->evaluate_scalar(scalar_id))
  {
    // get pre-factor for this scalar
    const double pre_factor =
        variable_manager()->get_pre_factor_mass_matrix_od_mesh(scalar_id, fac);

    scatra_ele_calc::calc_lin_mass_od_mesh(element_matrix, scalar_id, number_of_dofs_per_node_mesh,
        rhsfac, pre_factor, densam, densnp, phinp, hist, J, dJ_dmesh);
  }
}

/*-------------------------------------------------------------------- *
 *---------------------------------------------------------------------*/
template <Core::FE::CellType discretization_type>
void Discret::Elements::ScaTraEleCalcPorofluidPressureBased<discretization_type>::
    calc_hist_and_source_od_mesh(Core::LinAlg::SerialDenseMatrix& element_matrix,
        const int scalar_id, const int number_of_dofs_per_node_mesh, const double fac,
        const double rhsint, const double J, const Core::LinAlg::Matrix<1, nsd_ * nen_>& dJ_dmesh,
        const double densnp)
{
  // case of zero saturation/volfrac
  if (variable_manager()->evaluate_scalar(scalar_id))
  {
    const int number_of_fluid_phases =
        variable_manager()->porofluid_phase_manager()->num_fluid_phases();

    // get pre-factor for this scalar
    const double pre_factor = variable_manager()->get_pre_factor_hist_and_source_od_mesh(
        scalar_id, fac, densnp, scatra_ele_calc::scatravarmanager_->hist(scalar_id), rhsint);

    // linearization of mesh motion: call base class with correct pre-factor
    scatra_ele_calc::calc_hist_and_source_od_mesh(element_matrix, scalar_id,
        number_of_dofs_per_node_mesh, pre_factor, 1.0, J, dJ_dmesh, densnp);

    // linearization of advanced reaction terms
    const std::shared_ptr<ScaTraEleReaManagerAdvReac> reaction_manager =
        advanced_reaction::rea_manager();
    if (reaction_manager->active() &&
        variable_manager()->porofluid_phase_manager()->porosity_depends_on_struct())
    {
      const std::vector<double> reaction_derivatives =
          reaction_manager->get_rea_body_force_deriv_vector_add_variables(scalar_id);

      // porosity deriv at [2* number_of_fluid_phases]:
      // d reac / d d = d reac / d poro * d poro / d d
      // with
      // dporo/dd = dporo/dJ * dJ/dd = dporosity/dJ * J * N_x
      // J denotes the determinant of the deformation gradient, i.e.
      // det F = det ( d x / d X ) = det(dx/ds) * ( det(dX/ds) )^-1

      const double porosity_derivative =
          reaction_derivatives[2 * number_of_fluid_phases] *
          variable_manager()->porofluid_phase_manager()->jacobian_def_grad() *
          variable_manager()->porofluid_phase_manager()->porosity_deriv_wrt_jacobian_def_grad();

      for (unsigned vi = 0; vi < nen_; ++vi)
      {
        const int fvi = vi * scatra_ele_calc::numdofpernode_ + scalar_id;
        // TODO: gen-alpha
        const double v = scatra_ele_calc::funct_(vi) * porosity_derivative *
                         scatra_ele_calc::scatraparatimint_->time_fac() * fac * (-1.0) /
                         variable_manager()->get_density(scalar_id);

        for (unsigned ui = 0; ui < nen_; ++ui)
        {
          for (unsigned idim = 0; idim < nsd_; ++idim)
          {
            const int fui = ui * nsd_ + idim;

            element_matrix(fvi, fui) += v * scatra_ele_calc::derxy_(idim, ui);
          }
        }
      }
    }
  }
}

/*-------------------------------------------------------------------- *
 *---------------------------------------------------------------------*/
template <Core::FE::CellType discretization_type>
void Discret::Elements::ScaTraEleCalcPorofluidPressureBased<discretization_type>::calc_diff_od_mesh(
    Core::LinAlg::SerialDenseMatrix& element_matrix, const int scalar_id,
    const int number_of_dofs_per_node_mesh, const double diffcoeff, const double fac,
    const double rhsfac, const double J, const Core::LinAlg::Matrix<nsd_, 1>& gradphi,
    const Core::LinAlg::Matrix<nsd_, 1>& convelint,
    const Core::LinAlg::Matrix<1, nsd_ * nen_>& dJ_dmesh)
{
  // case of zero saturation/volfrac
  if (variable_manager()->evaluate_scalar(scalar_id))
  {
    // call base class
    scatra_ele_calc::calc_diff_od_mesh(element_matrix, scalar_id, number_of_dofs_per_node_mesh,
        diffcoeff, fac, rhsfac, J, gradphi, convelint, dJ_dmesh);

    // get pre-factor for this scalar
    const double pre_factor =
        variable_manager()->get_pre_factor_diffusion_od_mesh(scalar_id, rhsfac, diffcoeff);

    if (fabs(pre_factor) > 1.0e-12)
    {
      for (unsigned vi = 0; vi < nen_; ++vi)
      {
        const int fvi = vi * scatra_ele_calc::numdofpernode_ + scalar_id;

        double laplacian(0.0);
        scatra_ele_calc::get_laplacian_weak_form_rhs(laplacian, gradphi, vi);
        const double v = pre_factor * laplacian;
        for (unsigned ui = 0; ui < nen_; ++ui)
        {
          for (unsigned idim = 0; idim < nsd_; ++idim)
          {
            const int fui = ui * nsd_ + idim;

            element_matrix(fvi, fui) += v * scatra_ele_calc::derxy_(idim, ui);
          }
        }
      }
    }
  }
}

/*-------------------------------------------------------------------- *
 *---------------------------------------------------------------------*/
template <Core::FE::CellType discretization_type>
void Discret::Elements::ScaTraEleCalcPorofluidPressureBased<
    discretization_type>::calc_react_od_mesh(Core::LinAlg::SerialDenseMatrix& element_matrix,
    const int scalar_id, const int number_of_dofs_per_node_mesh, const double rhsfac,
    const double rea_phi, const double J, const Core::LinAlg::Matrix<1, nsd_ * nen_>& dJ_dmesh)
{
  // case of zero saturation/volfrac
  if (variable_manager()->evaluate_scalar(scalar_id))
  {
    // get pre-factor for this scalar
    const double pre_factor =
        variable_manager()->get_pre_factor_mass_matrix_od_mesh(scalar_id, rhsfac);

    scatra_ele_calc::calc_react_od_mesh(
        element_matrix, scalar_id, number_of_dofs_per_node_mesh, pre_factor, rea_phi, J, dJ_dmesh);
  }
}

/*-------------------------------------------------------------------- *
 *---------------------------------------------------------------------*/
template <Core::FE::CellType discretization_type>
void Discret::Elements::ScaTraEleCalcPorofluidPressureBased<
    discretization_type>::calc_mat_conv_od_fluid(Core::LinAlg::SerialDenseMatrix& element_matrix,
    const int scalar_id, const int number_of_dofs_per_node_mesh, const double rhsfac,
    const double densnp, const Core::LinAlg::Matrix<nsd_, 1>& gradphi)
{
  // case of zero saturation/volfrac
  if (variable_manager()->evaluate_scalar(scalar_id))
  {
    if (variable_manager()->get_species_type(scalar_id) ==
            Mat::ScaTraMatMultiPoro::SpeciesType::species_in_fluid ||
        variable_manager()->get_species_type(scalar_id) ==
            Mat::ScaTraMatMultiPoro::SpeciesType::species_in_volfrac)
    {
      const int total_number_of_porofluid_dofs =
          variable_manager()->porofluid_material()->num_mat();

      static Core::LinAlg::Matrix<nsd_, nsd_> diffusion_tensor(Core::LinAlg::Initialization::zero);
      variable_manager()->get_diffusion_tensor_fluid(
          scalar_id, diffusion_tensor, variable_manager()->get_phase_id(scalar_id));

      // gradphi^T * diffusion_tensor
      // TODO: not sure if this works for anisotropic fluid diffusion tensor
      static Core::LinAlg::Matrix<1, nsd_> gradient_phi_transpose_diffusion_tensor;
      gradient_phi_transpose_diffusion_tensor.multiply_tn(gradphi, diffusion_tensor);

      for (unsigned vi = 0; vi < nen_; ++vi)
      {
        const int fvi = vi * scatra_ele_calc::numdofpernode_ + scalar_id;
        const double v = rhsfac * scatra_ele_calc::funct_(vi) * (-1.0);

        for (unsigned ui = 0; ui < nen_; ++ui)
        {
          // get pre-factor vector for this scalar
          static std::vector pre_factor_convection_OD_porofluid(
              total_number_of_porofluid_dofs, 0.0);
          variable_manager()->get_pre_factor_convection_od_porofluid(ui,
              &pre_factor_convection_OD_porofluid, gradphi, gradient_phi_transpose_diffusion_tensor,
              scatra_ele_calc::funct_, scatra_ele_calc::derxy_,
              variable_manager()->get_phase_id(scalar_id));

          for (int dof_id = 0; dof_id < total_number_of_porofluid_dofs; ++dof_id)
          {
            const int fui = ui * total_number_of_porofluid_dofs + dof_id;
            element_matrix(fvi, fui) += v * pre_factor_convection_OD_porofluid[dof_id];
          }
        }
      }

      // TODO: linearization of dynamic viscosity w.r.t. dof, necessary? FD check does not fail
    }

    else if (variable_manager()->get_species_type(scalar_id) ==
             Mat::ScaTraMatMultiPoro::SpeciesType::species_temperature)
    {
      const int number_of_fluid_phases =
          variable_manager()->porofluid_phase_manager()->num_fluid_phases();
      const int number_of_volume_fractions =
          variable_manager()->porofluid_phase_manager()->num_vol_frac();

      for (int phase = 0; phase < number_of_fluid_phases + number_of_volume_fractions; ++phase)
      {
        const int total_number_of_porofluid_dofs =
            variable_manager()->porofluid_material()->num_mat();

        static Core::LinAlg::Matrix<nsd_, nsd_> diffusion_tensor;
        variable_manager()->get_diffusion_tensor_fluid(scalar_id, diffusion_tensor, phase);

        // gradphi^T * diffusion_tensor
        static Core::LinAlg::Matrix<1, nsd_> gradient_phi_transpose_diffusion_tensor;
        gradient_phi_transpose_diffusion_tensor.multiply_tn(gradphi, diffusion_tensor);

        // calculate density * heat capacity
        double density_heat_capacity =
            variable_manager()->get_heat_capacity(phase) * variable_manager()->get_density()[phase];

        for (unsigned vi = 0; vi < nen_; ++vi)
        {
          const int fvi = vi * scatra_ele_calc::numdofpernode_ + scalar_id;
          const double v = rhsfac * scatra_ele_calc::funct_(vi) * (-1.0) * density_heat_capacity;

          for (unsigned ui = 0; ui < nen_; ++ui)
          {
            // get pre-fac vector for this scalar
            static std::vector pre_factor_convection_OD_porofluid(
                total_number_of_porofluid_dofs, 0.0);
            variable_manager()->get_pre_factor_convection_od_porofluid(ui,
                &pre_factor_convection_OD_porofluid, gradphi,
                gradient_phi_transpose_diffusion_tensor, scatra_ele_calc::funct_,
                scatra_ele_calc::derxy_, phase);

            for (int dof_id = 0; dof_id < total_number_of_porofluid_dofs; ++dof_id)
            {
              const int fui = ui * total_number_of_porofluid_dofs + dof_id;
              element_matrix(fvi, fui) += v * pre_factor_convection_OD_porofluid[dof_id];
            }
          }
        }
      }
    }
  }
}

/*-------------------------------------------------------------------------- *
 *---------------------------------------------------------------------------*/
template <Core::FE::CellType discretization_type>
void Discret::Elements::ScaTraEleCalcPorofluidPressureBased<discretization_type>::
    calc_mat_conv_add_cons_od_fluid(Core::LinAlg::SerialDenseMatrix& element_matrix,
        const int scalar_id, const int number_of_dofs_per_node_mesh, const double timefacfac,
        const double densnp, const double phinp)
{
  FOUR_C_THROW(
      "calc_mat_conv_add_cons_od_fluid not yet available for "
      "scatra_ele_calc_porofluid_pressure_based.");
}

/*-------------------------------------------------------------------- *
 *---------------------------------------------------------------------*/
template <Core::FE::CellType discretization_type>
void Discret::Elements::ScaTraEleCalcPorofluidPressureBased<
    discretization_type>::calc_lin_mass_od_fluid(Core::LinAlg::SerialDenseMatrix& element_matrix,
    const int scalar_id, const int number_of_dofs_per_node_mesh, const double rhsfac,
    const double fac, const double densam, const double densnp, const double phinp,
    const double hist)
{
  // case of zero saturation/volfrac
  if (variable_manager()->evaluate_scalar(scalar_id))
  {
    if (scatra_ele_calc::scatraparatimint_->is_gen_alpha())
      FOUR_C_THROW("GenAlpha is not implemented for scatra_ele_calc_porofluid_pressure_based.");

    double vtrans = 0.0;
    vtrans = fac * densnp * phinp;
    const int total_number_of_porofluid_dofs = variable_manager()->porofluid_material()->num_mat();

    // get pre-fac vector for this scalar
    std::vector pre_factor_mass_matrix_OD_porofluid(total_number_of_porofluid_dofs, 0.0);
    variable_manager()->get_pre_factor_mass_matrix_od_porofluid(
        scalar_id, &pre_factor_mass_matrix_OD_porofluid);

    calc_mass_matrix_type_od_porofluid(element_matrix, scalar_id,
        &pre_factor_mass_matrix_OD_porofluid, total_number_of_porofluid_dofs, vtrans);
  }
}

/*---------------------------------------------------------------------- *
 *-----------------------------------------------------------------------*/
template <Core::FE::CellType discretization_type>
void Discret::Elements::ScaTraEleCalcPorofluidPressureBased<discretization_type>::
    calc_mass_matrix_type_od_porofluid(Core::LinAlg::SerialDenseMatrix& element_matrix,
        const int scalar_id, const std::vector<double>* per_factor_mass_matrix_OD_porofluid,
        const int total_number_of_porofluid_dofs, double pre_factor)
{
  for (unsigned vi = 0; vi < nen_; ++vi)
  {
    const double v = pre_factor * scatra_ele_calc::funct_(vi);
    const int fvi = vi * scatra_ele_calc::numdofpernode_ + scalar_id;

    for (unsigned ui = 0; ui < nen_; ++ui)
    {
      const double vfunct = v * scatra_ele_calc::funct_(ui);
      for (int dof_id = 0; dof_id < total_number_of_porofluid_dofs; ++dof_id)
      {
        const int fui = ui * total_number_of_porofluid_dofs + dof_id;

        element_matrix(fvi, fui) += vfunct * (*per_factor_mass_matrix_OD_porofluid)[dof_id];
      }
    }
  }
}

/*----------------------------------------------------------------------------- *
 *------------------------------------------------------------------------------*/
template <Core::FE::CellType discretization_type>
void Discret::Elements::ScaTraEleCalcPorofluidPressureBased<discretization_type>::
    calc_hist_and_source_od_fluid(Core::LinAlg::SerialDenseMatrix& element_matrix,
        const int scalar_id, const int number_of_dofs_per_node_mesh, const double fac,
        const double rhsint, const double densnp)
{
  // case of zero saturation/volfrac
  if (variable_manager()->evaluate_scalar(scalar_id))
  {
    const int number_of_fluid_phases =
        variable_manager()->porofluid_phase_manager()->num_fluid_phases();
    const int number_of_volume_fractions =
        variable_manager()->porofluid_phase_manager()->num_vol_frac();
    // total number of porofluid dofs additionally includes volume fraction pressures
    const int total_number_of_porofluid_dofs = variable_manager()->porofluid_material()->num_mat();

    // linearization of history: pre-factor is densnp = porosity * rho * S
    // --> porosity and saturation have to be linearized
    const double vrhs = -1.0 * fac * scatra_ele_calc::scatravarmanager_->hist(scalar_id) * densnp;

    // get pre-factor vector for this scalar
    static std::vector pre_factor_mass_matrix_OD_porofluid(total_number_of_porofluid_dofs, 0.0);
    variable_manager()->get_pre_factor_mass_matrix_od_porofluid(
        scalar_id, &pre_factor_mass_matrix_OD_porofluid);

    calc_mass_matrix_type_od_porofluid(element_matrix, scalar_id,
        &pre_factor_mass_matrix_OD_porofluid, total_number_of_porofluid_dofs, vrhs);

    // linearization of advanced reaction terms
    const std::shared_ptr<ScaTraEleReaManagerAdvReac> reaction_manager =
        advanced_reaction::rea_manager();
    if (reaction_manager->active())
    {
      const std::vector<double> reaction_derivatives =
          reaction_manager->get_rea_body_force_deriv_vector_add_variables(scalar_id);

      // derivatives after primary variables of fluid
      std::vector phi_derivatives(total_number_of_porofluid_dofs, 0.0);

      for (int phase_id = 0; phase_id < number_of_fluid_phases; phase_id++)
      {
        // porosity deriv at 2 * number_of_fluid_phases:
        // d reac / d phi_i = (d reac / d poro) * (d poro / d phi_i)
        phi_derivatives[phase_id] +=
            reaction_derivatives[2 * number_of_fluid_phases] *
            variable_manager()->porofluid_phase_manager()->porosity_deriv(phase_id);
        for (int j = 0; j < number_of_fluid_phases; j++)
        {
          // pressure derivs at [0 ... number_of_fluid_phases]:
          // d reac / d phi_i = (d reac / d pressure_j) * (d pressure_j / d phi_i)
          // saturation derivs at [number_of_fluid_phases ... (2*number_of_fluid_phases-1)]:
          // d reac / d phi_i = (d reac / d sat_j) * (d sat_j / d phi_i)
          phi_derivatives[phase_id] +=
              reaction_derivatives[j + number_of_fluid_phases] *
                  variable_manager()->porofluid_phase_manager()->saturation_deriv(j, phase_id) +
              reaction_derivatives[j] *
                  variable_manager()->porofluid_phase_manager()->pressure_deriv(j, phase_id);
        }
      }

      for (int volume_fraction_id = 0; volume_fraction_id < number_of_volume_fractions;
          volume_fraction_id++)
      {
        // derivatives after volume fractions at [2*number_of_fluid_phases + 1 + volume_fraction_id]
        phi_derivatives[volume_fraction_id + number_of_fluid_phases] +=
            reaction_derivatives[2 * number_of_fluid_phases + 1 + volume_fraction_id] +
            reaction_derivatives[2 * number_of_fluid_phases] *
                variable_manager()->porofluid_phase_manager()->porosity_deriv(
                    volume_fraction_id + number_of_fluid_phases);
        // derivatives after volume fraction pressures at
        // [2*number_of_fluid_phases + number_of_volume_fractions+ 1 + volume_fraction_id]
        phi_derivatives[volume_fraction_id + number_of_fluid_phases + number_of_volume_fractions] +=
            reaction_derivatives[2 * number_of_fluid_phases + number_of_volume_fractions + 1 +
                                 volume_fraction_id];
      }

      // fill matrix
      for (unsigned vi = 0; vi < nen_; ++vi)
      {
        const int fvi = vi * scatra_ele_calc::numdofpernode_ + scalar_id;
        // TODO: gen-alpha?
        const double v = scatra_ele_calc::funct_(vi) *
                         scatra_ele_calc::scatraparatimint_->time_fac() * fac * (-1.0) /
                         variable_manager()->get_density(scalar_id);

        for (unsigned ui = 0; ui < nen_; ++ui)
        {
          const double vfunct = v * scatra_ele_calc::funct_(ui);

          for (int dof_id = 0; dof_id < total_number_of_porofluid_dofs; ++dof_id)
          {
            const int fui = ui * total_number_of_porofluid_dofs + dof_id;

            element_matrix(fvi, fui) += vfunct * phi_derivatives[dof_id];
          }
        }
      }
    }
  }
}

/*-------------------------------------------------------------------- *
 *---------------------------------------------------------------------*/
template <Core::FE::CellType discretization_type>
void Discret::Elements::ScaTraEleCalcPorofluidPressureBased<
    discretization_type>::calc_react_od_fluid(Core::LinAlg::SerialDenseMatrix& element_matrix,
    const int scalar_id, const int number_of_dofs_per_node_mesh, const double rhsfac,
    const double rea_phi)
{
  if (scatra_ele_calc::reamanager_->active() && variable_manager()->evaluate_scalar(scalar_id))
  {
    // total number of porofluid dofs additionally includes volume fraction pressures
    const int total_number_of_porofluid_dofs = variable_manager()->porofluid_material()->num_mat();

    const double vrhs = rhsfac * rea_phi;

    // get pre-factor vector for this scalar
    static std::vector pre_factor_mass_matrix_OD_porofluid(total_number_of_porofluid_dofs, 0.0);
    variable_manager()->get_pre_factor_mass_matrix_od_porofluid(
        scalar_id, &pre_factor_mass_matrix_OD_porofluid);

    calc_mass_matrix_type_od_porofluid(element_matrix, scalar_id,
        &pre_factor_mass_matrix_OD_porofluid, total_number_of_porofluid_dofs, vrhs);
  }
}

/*------------------------------------------------------------------ *
 *----------------------------------------------------------------   */
template <Core::FE::CellType discretization_type>
void Discret::Elements::ScaTraEleCalcPorofluidPressureBased<
    discretization_type>::calc_diff_od_fluid(Core::LinAlg::SerialDenseMatrix& element_matrix,
    const int scalar_id, const int number_of_dofs_per_node_mesh, const double rhsfac,
    const Core::LinAlg::Matrix<nsd_, 1>& gradphi)
{
  // case of zero saturation/volfrac
  if (variable_manager()->evaluate_scalar(scalar_id))
  {
    // total number of porofluid dofs additionally includes volume fraction pressures
    const int total_number_of_porofluid_dofs = variable_manager()->porofluid_material()->num_mat();

    // get pre-fac vector for this scalar
    static std::vector pre_factor_diffusion_OD_porofluid(total_number_of_porofluid_dofs, 0.0);
    variable_manager()->get_pre_factor_diffusion_od_porofluid(scalar_id, rhsfac,
        scatra_ele_calc::diffmanager_->get_isotropic_diff(scalar_id),
        &pre_factor_diffusion_OD_porofluid);

    for (unsigned vi = 0; vi < nen_; ++vi)
    {
      const int fvi = vi * scatra_ele_calc::numdofpernode_ + scalar_id;

      double laplacian(0.0);
      scatra_ele_calc::get_laplacian_weak_form_rhs(laplacian, gradphi, vi);

      for (unsigned ui = 0; ui < nen_; ++ui)
      {
        const double laplacian_shape_function = laplacian * scatra_ele_calc::funct_(ui);

        // derivative w.r.t. porofluid variables
        for (int dof_id = 0; dof_id < total_number_of_porofluid_dofs; ++dof_id)
        {
          const int fui = ui * total_number_of_porofluid_dofs + dof_id;

          element_matrix(fvi, fui) +=
              laplacian_shape_function * pre_factor_diffusion_OD_porofluid[dof_id];
        }
      }
    }
  }
}

/*---------------------------------------------------------------------*
 *---------------------------------------------------------------------*/
template <Core::FE::CellType discretization_type>
void Discret::Elements::ScaTraEleCalcPorofluidPressureBased<discretization_type>::
    apply_shape_derivs_pressure_grad(Core::LinAlg::SerialDenseMatrix& element_matrix,
        const int scalar_id, const double vrhs, const Core::LinAlg::Matrix<nsd_, 1>& gradient_phi,
        const Core::LinAlg::Matrix<nsd_, 1> pressure_gradient_reference)
{
  if (nsd_ == 3)
  {
    const double xjm_0_0 = scatra_ele_calc::xjm_(0, 0);
    const double xjm_0_1 = scatra_ele_calc::xjm_(0, 1);
    const double xjm_0_2 = scatra_ele_calc::xjm_(0, 2);
    const double xjm_1_0 = scatra_ele_calc::xjm_(1, 0);
    const double xjm_1_1 = scatra_ele_calc::xjm_(1, 1);
    const double xjm_1_2 = scatra_ele_calc::xjm_(1, 2);
    const double xjm_2_0 = scatra_ele_calc::xjm_(2, 0);
    const double xjm_2_1 = scatra_ele_calc::xjm_(2, 1);
    const double xjm_2_2 = scatra_ele_calc::xjm_(2, 2);

    const double refgradpres_0 = pressure_gradient_reference(0);
    const double refgradpres_1 = pressure_gradient_reference(1);
    const double refgradpres_2 = pressure_gradient_reference(2);

    const double gradphi_0 = gradient_phi(0);
    const double gradphi_1 = gradient_phi(1);
    const double gradphi_2 = gradient_phi(2);

    for (unsigned ui = 0; ui < nen_; ++ui)
    {
      const double v00 =
          +gradphi_1 * (refgradpres_0 * (scatra_ele_calc::deriv_(2, ui) * xjm_1_2 -
                                            scatra_ele_calc::deriv_(1, ui) * xjm_2_2) +
                           refgradpres_1 * (scatra_ele_calc::deriv_(0, ui) * xjm_2_2 -
                                               scatra_ele_calc::deriv_(2, ui) * xjm_0_2) +
                           refgradpres_2 * (scatra_ele_calc::deriv_(1, ui) * xjm_0_2 -
                                               scatra_ele_calc::deriv_(0, ui) * xjm_1_2)) +
          gradphi_2 * (refgradpres_0 * (scatra_ele_calc::deriv_(1, ui) * xjm_2_1 -
                                           scatra_ele_calc::deriv_(2, ui) * xjm_1_1) +
                          refgradpres_1 * (scatra_ele_calc::deriv_(2, ui) * xjm_0_1 -
                                              scatra_ele_calc::deriv_(0, ui) * xjm_2_1) +
                          refgradpres_2 * (scatra_ele_calc::deriv_(0, ui) * xjm_1_1 -
                                              scatra_ele_calc::deriv_(1, ui) * xjm_0_1));
      const double v01 =
          +gradphi_0 * (refgradpres_0 * (scatra_ele_calc::deriv_(1, ui) * xjm_2_2 -
                                            scatra_ele_calc::deriv_(2, ui) * xjm_1_2) +
                           refgradpres_1 * (scatra_ele_calc::deriv_(2, ui) * xjm_0_2 -
                                               scatra_ele_calc::deriv_(0, ui) * xjm_2_2) +
                           refgradpres_2 * (scatra_ele_calc::deriv_(0, ui) * xjm_1_2 -
                                               scatra_ele_calc::deriv_(1, ui) * xjm_0_2)) +
          gradphi_2 * (refgradpres_0 * (scatra_ele_calc::deriv_(2, ui) * xjm_1_0 -
                                           scatra_ele_calc::deriv_(1, ui) * xjm_2_0) +
                          refgradpres_1 * (scatra_ele_calc::deriv_(0, ui) * xjm_2_0 -
                                              scatra_ele_calc::deriv_(2, ui) * xjm_0_0) +
                          refgradpres_2 * (scatra_ele_calc::deriv_(1, ui) * xjm_0_0 -
                                              scatra_ele_calc::deriv_(0, ui) * xjm_1_0));
      const double v02 =
          +gradphi_0 * (refgradpres_0 * (scatra_ele_calc::deriv_(2, ui) * xjm_1_1 -
                                            scatra_ele_calc::deriv_(1, ui) * xjm_2_1) +
                           refgradpres_1 * (scatra_ele_calc::deriv_(0, ui) * xjm_2_1 -
                                               scatra_ele_calc::deriv_(2, ui) * xjm_0_1) +
                           refgradpres_2 * (scatra_ele_calc::deriv_(1, ui) * xjm_0_1 -
                                               scatra_ele_calc::deriv_(0, ui) * xjm_1_1)) +
          gradphi_1 * (refgradpres_0 * (scatra_ele_calc::deriv_(1, ui) * xjm_2_0 -
                                           scatra_ele_calc::deriv_(2, ui) * xjm_1_0) +
                          refgradpres_1 * (scatra_ele_calc::deriv_(2, ui) * xjm_0_0 -
                                              scatra_ele_calc::deriv_(0, ui) * xjm_2_0) +
                          refgradpres_2 * (scatra_ele_calc::deriv_(0, ui) * xjm_1_0 -
                                              scatra_ele_calc::deriv_(1, ui) * xjm_0_0));

      for (unsigned vi = 0; vi < nen_; ++vi)
      {
        const int fvi = vi * scatra_ele_calc::numdofpernode_ + scalar_id;
        const double v = vrhs * scatra_ele_calc::funct_(vi);

        element_matrix(fvi, ui * 3 + 0) += v * v00;
        element_matrix(fvi, ui * 3 + 1) += v * v01;
        element_matrix(fvi, ui * 3 + 2) += v * v02;
      }
    }
  }
  else if (nsd_ == 2)
  {
    const double refgradpres_0 = pressure_gradient_reference(0);
    const double refgradpres_1 = pressure_gradient_reference(1);

    const double gradphi_0 = gradient_phi(0);
    const double gradphi_1 = gradient_phi(1);

    for (unsigned ui = 0; ui < nen_; ++ui)
    {
      const double v00 = +gradphi_1 * (-refgradpres_0 * scatra_ele_calc::deriv_(1, ui) +
                                          refgradpres_1 * scatra_ele_calc::deriv_(0, ui));
      const double v01 = +gradphi_0 * (refgradpres_0 * scatra_ele_calc::deriv_(1, ui) -
                                          refgradpres_1 * scatra_ele_calc::deriv_(0, ui));

      for (unsigned vi = 0; vi < nen_; ++vi)
      {
        const int fvi = vi * scatra_ele_calc::numdofpernode_ + scalar_id;
        const double v = vrhs * scatra_ele_calc::funct_(vi);

        element_matrix(fvi, ui * 2 + 0) += v * v00;
        element_matrix(fvi, ui * 2 + 1) += v * v01;
      }
    }
  }
  else
  {
    FOUR_C_THROW("shapederivatives not implemented for 1D!");
  }
}

/*------------------------------------------------------------------------------*
 *------------------------------------------------------------------------------*/
// template classes

// 1D elements
template class Discret::Elements::ScaTraEleCalcPorofluidPressureBased<Core::FE::CellType::line2>;
template class Discret::Elements::ScaTraEleCalcPorofluidPressureBased<Core::FE::CellType::line3>;

// 2D elements
template class Discret::Elements::ScaTraEleCalcPorofluidPressureBased<Core::FE::CellType::tri3>;
template class Discret::Elements::ScaTraEleCalcPorofluidPressureBased<Core::FE::CellType::tri6>;
template class Discret::Elements::ScaTraEleCalcPorofluidPressureBased<Core::FE::CellType::quad4>;
// template class
// Discret::Elements::ScaTraEleCalcPorofluidPressureBased<Core::FE::CellType::quad8>;
template class Discret::Elements::ScaTraEleCalcPorofluidPressureBased<Core::FE::CellType::quad9>;

// 3D elements
template class Discret::Elements::ScaTraEleCalcPorofluidPressureBased<Core::FE::CellType::hex8>;
// template class
// Discret::Elements::ScaTraEleCalcPorofluidPressureBased<Core::FE::CellType::hex20>;
template class Discret::Elements::ScaTraEleCalcPorofluidPressureBased<Core::FE::CellType::hex27>;
template class Discret::Elements::ScaTraEleCalcPorofluidPressureBased<Core::FE::CellType::tet4>;
template class Discret::Elements::ScaTraEleCalcPorofluidPressureBased<Core::FE::CellType::tet10>;
// template class
// Discret::Elements::ScaTraEleCalcPorofluidPressureBased<Core::FE::CellType::wedge6>;
template class Discret::Elements::ScaTraEleCalcPorofluidPressureBased<Core::FE::CellType::pyramid5>;
template class Discret::Elements::ScaTraEleCalcPorofluidPressureBased<Core::FE::CellType::nurbs9>;
// template class
// Discret::Elements::ScaTraEleCalcPorofluidPressureBased<Core::FE::CellType::nurbs27>;

FOUR_C_NAMESPACE_CLOSE
