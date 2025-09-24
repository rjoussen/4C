// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_particle_interaction_sph_heatsource.hpp"

#include "4C_global_data.hpp"
#include "4C_particle_engine_container.hpp"
#include "4C_particle_engine_interface.hpp"
#include "4C_particle_interaction_material_handler.hpp"
#include "4C_particle_interaction_sph_neighbor_pairs.hpp"
#include "4C_particle_interaction_utils.hpp"
#include "4C_utils_exceptions.hpp"
#include "4C_utils_function.hpp"

#include <Teuchos_StandardParameterEntryValidators.hpp>
#include <Teuchos_TimeMonitor.hpp>

FOUR_C_NAMESPACE_OPEN

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
ParticleInteraction::SPHHeatSourceBase::SPHHeatSourceBase(const Teuchos::ParameterList& params)
    : params_sph_(params), heatsourcefctnumber_(params.get<int>("HEATSOURCE_FUNCT"))
{
  // empty constructor
}

void ParticleInteraction::SPHHeatSourceBase::init()
{
  // nothing to do
}

void ParticleInteraction::SPHHeatSourceBase::setup(
    const std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface,
    const std::shared_ptr<ParticleInteraction::MaterialHandler> particlematerial,
    const std::shared_ptr<ParticleInteraction::SPHNeighborPairs> neighborpairs)
{
  // set interface to particle engine
  particleengineinterface_ = particleengineinterface;

  // set particle container bundle
  particlecontainerbundle_ = particleengineinterface_->get_particle_container_bundle();

  // set particle material handler
  particlematerial_ = particlematerial;

  // set neighbor pair handler
  neighborpairs_ = neighborpairs;

  // determine size of vectors indexed by particle types
  const int typevectorsize = *(--particlecontainerbundle_->get_particle_types().end()) + 1;

  // allocate memory to hold particle types
  thermomaterial_.resize(typevectorsize);

  // iterate over particle types
  for (const auto& type_i : particlecontainerbundle_->get_particle_types())
    thermomaterial_[type_i] = dynamic_cast<const Mat::PAR::ParticleMaterialThermo*>(
        particlematerial_->get_ptr_to_particle_mat_parameter(type_i));

  // set of potential absorbing particle types
  std::set<PARTICLEENGINE::TypeEnum> potentialabsorbingtypes = {
      PARTICLEENGINE::Phase1, PARTICLEENGINE::Phase2, PARTICLEENGINE::RigidPhase};

  // iterate over particle types
  for (const auto& type_i : particlecontainerbundle_->get_particle_types())
  {
    // determine absorbing particle types
    if (thermomaterial_[type_i]->thermalAbsorptivity_ > 0.0)
    {
      // safety check
      if (not potentialabsorbingtypes.count(type_i))
        FOUR_C_THROW("thermal absorptivity for particles of type '{}' not possible!",
            PARTICLEENGINE::enum_to_type_name(type_i));

      absorbingtypes_.insert(type_i);
    }
    // determine non-absorbing particle types
    else if (potentialabsorbingtypes.count(type_i))
    {
      nonabsorbingtypes_.insert(type_i);
    }
  }
}

ParticleInteraction::SPHHeatSourceVolume::SPHHeatSourceVolume(const Teuchos::ParameterList& params)
    : ParticleInteraction::SPHHeatSourceBase(params)
{
  // empty constructor
}

void ParticleInteraction::SPHHeatSourceVolume::evaluate_heat_source(const double& evaltime) const
{
  TEUCHOS_FUNC_TIME_MONITOR("ParticleInteraction::SPHHeatSourceVolume::EvaluateHeatSource");

  // init vector containing evaluated function
  std::vector<double> funct(1);

  // get reference to function
  const auto& function =
      Global::Problem::instance()->function_by_id<Core::Utils::FunctionOfSpaceTime>(
          heatsourcefctnumber_);

  // safety check
  if (function.number_components() != 1)
    FOUR_C_THROW("dimension of function defining heat source is not one!");

  // iterate over absorbing particle types
  for (const auto& type_i : absorbingtypes_)
  {
    // get container of owned particles of current particle type
    PARTICLEENGINE::ParticleContainer* container_i =
        particlecontainerbundle_->get_specific_container(type_i, PARTICLEENGINE::Owned);

    // get material for current particle type
    const Mat::PAR::ParticleMaterialBase* basematerial_i =
        particlematerial_->get_ptr_to_particle_mat_parameter(type_i);

    const Mat::PAR::ParticleMaterialThermo* thermomaterial_i = thermomaterial_[type_i];

    // iterate over particles in container
    for (int particle_i = 0; particle_i < container_i->particles_stored(); ++particle_i)
    {
      // get pointer to particle states
      const double* dens_i =
          (container_i->have_stored_state(PARTICLEENGINE::Density))
              ? container_i->get_ptr_to_state(PARTICLEENGINE::Density, particle_i)
              : &(basematerial_i->initDensity_);

      const double* pos_i = container_i->get_ptr_to_state(PARTICLEENGINE::Position, particle_i);
      double* tempdot_i = container_i->get_ptr_to_state(PARTICLEENGINE::TemperatureDot, particle_i);

      // evaluate function defining heat source
      funct = function.evaluate_time_derivative(std::span(pos_i, 3), evaltime, 0, 0);

      // add contribution of heat source
      tempdot_i[0] += thermomaterial_i->thermalAbsorptivity_ * funct[0] *
                      thermomaterial_i->invThermalCapacity_ / dens_i[0];
    }
  }
}

ParticleInteraction::SPHHeatSourceSurface::SPHHeatSourceSurface(
    const Teuchos::ParameterList& params)
    : ParticleInteraction::SPHHeatSourceBase(params), eval_direction_(false)
{
  // empty constructor
}

void ParticleInteraction::SPHHeatSourceSurface::init()
{
  // call base class init
  SPHHeatSourceBase::init();

  // init heat source direction vector
  double value;
  std::istringstream directionstream(
      Teuchos::getNumericStringParameter(params_sph_, "HEATSOURCE_DIRECTION"));

  while (directionstream >> value) direction_.push_back(value);

  // safety check
  if (static_cast<int>(direction_.size()) != 3)
    FOUR_C_THROW("dimension (dim = {}) of heat source direction vector is wrong!",
        static_cast<int>(direction_.size()));

  // normalize heat source direction vector
  const double direction_norm = Utils::vec_norm_two(direction_.data());
  if (direction_norm > 0.0)
  {
    eval_direction_ = true;
    Utils::vec_set_scale(direction_.data(), 1.0 / direction_norm, direction_.data());
  }
}

void ParticleInteraction::SPHHeatSourceSurface::evaluate_heat_source(const double& evaltime) const
{
  TEUCHOS_FUNC_TIME_MONITOR("ParticleInteraction::SPHHeatSourceSurface::EvaluateHeatSource");

  // determine size of vectors indexed by particle types
  const int typevectorsize = *(--absorbingtypes_.end()) + 1;

  // colorfield gradient of absorbing interface particles
  std::vector<std::vector<std::vector<double>>> cfg_i(typevectorsize);

  // iterate over absorbing particle types
  for (const auto& type_i : absorbingtypes_)
  {
    // get container of owned particles of current particle type
    PARTICLEENGINE::ParticleContainer* container_i =
        particlecontainerbundle_->get_specific_container(type_i, PARTICLEENGINE::Owned);

    // get number of particles stored in container
    const int particlestored = container_i->particles_stored();

    // allocate memory
    cfg_i[type_i].assign(particlestored, std::vector<double>(3, 0.0));
  }

  // get relevant particle pair indices
  std::vector<int> relindices;
  neighborpairs_->get_relevant_particle_pair_indices_for_disjoint_combination(
      absorbingtypes_, nonabsorbingtypes_, relindices);

  // iterate over relevant particle pairs
  for (const int particlepairindex : relindices)
  {
    const SPHParticlePair& particlepair =
        neighborpairs_->get_ref_to_particle_pair_data()[particlepairindex];

    // access values of local index tuples of particle i and j
    PARTICLEENGINE::TypeEnum type_i;
    PARTICLEENGINE::StatusEnum status_i;
    int particle_i;
    std::tie(type_i, status_i, particle_i) = particlepair.tuple_i_;

    PARTICLEENGINE::TypeEnum type_j;
    PARTICLEENGINE::StatusEnum status_j;
    int particle_j;
    std::tie(type_j, status_j, particle_j) = particlepair.tuple_j_;

    // get corresponding particle containers
    PARTICLEENGINE::ParticleContainer* container_i =
        particlecontainerbundle_->get_specific_container(type_i, status_i);

    PARTICLEENGINE::ParticleContainer* container_j =
        particlecontainerbundle_->get_specific_container(type_j, status_j);

    // get material for particle types
    const Mat::PAR::ParticleMaterialBase* material_i =
        particlematerial_->get_ptr_to_particle_mat_parameter(type_i);

    const Mat::PAR::ParticleMaterialBase* material_j =
        particlematerial_->get_ptr_to_particle_mat_parameter(type_j);

    // get pointer to particle states
    const double* mass_i = container_i->get_ptr_to_state(PARTICLEENGINE::Mass, particle_i);

    const double* dens_i = container_i->have_stored_state(PARTICLEENGINE::Density)
                               ? container_i->get_ptr_to_state(PARTICLEENGINE::Density, particle_i)
                               : &(material_i->initDensity_);

    const double* mass_j = container_j->get_ptr_to_state(PARTICLEENGINE::Mass, particle_j);

    const double* dens_j = container_j->have_stored_state(PARTICLEENGINE::Density)
                               ? container_j->get_ptr_to_state(PARTICLEENGINE::Density, particle_j)
                               : &(material_j->initDensity_);

    // (current) volume of particle i and j
    const double V_i = mass_i[0] / dens_i[0];
    const double V_j = mass_j[0] / dens_j[0];

    const double fac = (Utils::pow<2>(V_i) + Utils::pow<2>(V_j)) / (dens_i[0] + dens_j[0]);

    // evaluate contribution of neighboring particle j
    if (absorbingtypes_.count(type_i))
    {
      // sum contribution of neighboring particle j
      Utils::vec_add_scale(cfg_i[type_i][particle_i].data(),
          dens_i[0] / V_i * fac * particlepair.dWdrij_, particlepair.e_ij_);
    }

    // evaluate contribution of neighboring particle i
    if (absorbingtypes_.count(type_j) and status_j == PARTICLEENGINE::Owned)
    {
      // sum contribution of neighboring particle i
      Utils::vec_add_scale(cfg_i[type_j][particle_j].data(),
          -dens_j[0] / V_j * fac * particlepair.dWdrji_, particlepair.e_ij_);
    }
  }

  // init vector containing evaluated function
  std::vector<double> funct(1);

  // get reference to function
  const auto& function =
      Global::Problem::instance()->function_by_id<Core::Utils::FunctionOfSpaceTime>(
          heatsourcefctnumber_);

  // safety check
  if (function.number_components() != 1)
    FOUR_C_THROW("dimension of function defining heat source is not one!");

  // iterate over absorbing particle types
  for (const auto& type_i : absorbingtypes_)
  {
    // get container of owned particles of current particle type
    PARTICLEENGINE::ParticleContainer* container_i =
        particlecontainerbundle_->get_specific_container(type_i, PARTICLEENGINE::Owned);

    // get material for current particle type
    const Mat::PAR::ParticleMaterialBase* basematerial_i =
        particlematerial_->get_ptr_to_particle_mat_parameter(type_i);

    const Mat::PAR::ParticleMaterialThermo* thermomaterial_i = thermomaterial_[type_i];

    // iterate over particles in container
    for (int particle_i = 0; particle_i < container_i->particles_stored(); ++particle_i)
    {
      // norm of colorfield gradient of absorbing interface particles
      const double f_i = Utils::vec_norm_two(cfg_i[type_i][particle_i].data());

      // no heat source contribution to current particle
      if (not(f_i > 0.0)) continue;

      // projection of colorfield gradient with heat source direction
      const double f_i_proj =
          eval_direction_ ? -Utils::vec_dot(direction_.data(), cfg_i[type_i][particle_i].data())
                          : f_i;

      // heat source contribution only for surface opposing heat source
      if (f_i_proj < 0.0) continue;

      // get pointer to particle states
      const double* dens_i =
          (container_i->have_stored_state(PARTICLEENGINE::Density))
              ? container_i->get_ptr_to_state(PARTICLEENGINE::Density, particle_i)
              : &(basematerial_i->initDensity_);

      const double* pos_i = container_i->get_ptr_to_state(PARTICLEENGINE::Position, particle_i);
      double* tempdot_i = container_i->get_ptr_to_state(PARTICLEENGINE::TemperatureDot, particle_i);

      // evaluate function defining heat source
      funct = function.evaluate_time_derivative(std::span(pos_i, 3), evaltime, 0, 0);

      // add contribution of heat source
      tempdot_i[0] += f_i_proj * thermomaterial_i->thermalAbsorptivity_ * funct[0] *
                      thermomaterial_i->invThermalCapacity_ / dens_i[0];
    }
  }
}

FOUR_C_NAMESPACE_CLOSE
