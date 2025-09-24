// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_particle_algorithm_dirichlet_bc.hpp"

#include "4C_global_data.hpp"
#include "4C_particle_algorithm_utils.hpp"
#include "4C_particle_engine_container.hpp"
#include "4C_particle_engine_container_bundle.hpp"
#include "4C_particle_engine_enums.hpp"
#include "4C_particle_engine_interface.hpp"
#include "4C_utils_function.hpp"

FOUR_C_NAMESPACE_OPEN

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
PARTICLEALGORITHM::DirichletBoundaryConditionHandler::DirichletBoundaryConditionHandler(
    const Teuchos::ParameterList& params)
    : params_(params)
{
  // empty constructor
}

void PARTICLEALGORITHM::DirichletBoundaryConditionHandler::init()
{
  // get control parameters for initial/boundary conditions
  const Teuchos::ParameterList& params_conditions =
      params_.sublist("INITIAL AND BOUNDARY CONDITIONS");

  // read parameters relating particle types to values
  PARTICLEALGORITHM::Utils::read_params_types_related_to_values(
      params_conditions, "DIRICHLET_BOUNDARY_CONDITION", dirichletbctypetofunctid_);

  // iterate over particle types and insert into set
  for (auto& typeIt : dirichletbctypetofunctid_) typessubjectedtodirichletbc_.insert(typeIt.first);
}

void PARTICLEALGORITHM::DirichletBoundaryConditionHandler::setup(
    const std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface)
{
  // set interface to particle engine
  particleengineinterface_ = particleengineinterface;
}

void PARTICLEALGORITHM::DirichletBoundaryConditionHandler::insert_particle_states_of_particle_types(
    std::map<PARTICLEENGINE::TypeEnum, std::set<PARTICLEENGINE::StateEnum>>& particlestatestotypes)
    const
{
  // iterate over particle types subjected to dirichlet boundary conditions
  for (auto& particleType : typessubjectedtodirichletbc_)
  {
    // insert states for types subjected to dirichlet boundary conditions
    particlestatestotypes[particleType].insert(PARTICLEENGINE::ReferencePosition);
  }
}

void PARTICLEALGORITHM::DirichletBoundaryConditionHandler::set_particle_reference_position() const
{
  // get particle container bundle
  PARTICLEENGINE::ParticleContainerBundleShrdPtr particlecontainerbundle =
      particleengineinterface_->get_particle_container_bundle();

  // iterate over particle types subjected to dirichlet boundary conditions
  for (auto& particleType : typessubjectedtodirichletbc_)
  {
    // get container of owned particles of current particle type
    PARTICLEENGINE::ParticleContainer* container =
        particlecontainerbundle->get_specific_container(particleType, PARTICLEENGINE::Owned);

    // set particle reference position
    container->update_state(0.0, PARTICLEENGINE::ReferencePosition, 1.0, PARTICLEENGINE::Position);
  }
}

void PARTICLEALGORITHM::DirichletBoundaryConditionHandler::evaluate_dirichlet_boundary_condition(
    const double& evaltime, const bool evalpos, const bool evalvel, const bool evalacc) const
{
  // degree of maximal function derivative
  int deg = 0;
  if (evalacc)
    deg = 2;
  else if (evalvel)
    deg = 1;

  // get bounding box dimensions
  Core::LinAlg::Matrix<3, 2> boundingbox =
      particleengineinterface_->domain_bounding_box_corner_positions();

  // get bin size
  const std::array<double, 3> binsize = particleengineinterface_->bin_size();

  // init vector containing evaluated function and derivatives
  std::vector<double> functtimederiv(deg + 1);

  // get particle container bundle
  PARTICLEENGINE::ParticleContainerBundleShrdPtr particlecontainerbundle =
      particleengineinterface_->get_particle_container_bundle();

  // iterate over particle types subjected to dirichlet boundary conditions
  for (auto& typeIt : dirichletbctypetofunctid_)
  {
    // get type of particles
    PARTICLEENGINE::TypeEnum particleType = typeIt.first;

    // get container of owned particles of current particle type
    PARTICLEENGINE::ParticleContainer* container =
        particlecontainerbundle->get_specific_container(particleType, PARTICLEENGINE::Owned);

    // get number of particles stored in container
    const int particlestored = container->particles_stored();

    // no owned particles of current particle type
    if (particlestored <= 0) continue;

    // get id of function
    const int functid = typeIt.second;

    // get reference to function
    const auto& function =
        Global::Problem::instance()->function_by_id<Core::Utils::FunctionOfSpaceTime>(functid);

    // get pointer to particle states
    const double* refpos = container->get_ptr_to_state(PARTICLEENGINE::ReferencePosition, 0);
    double* pos = container->get_ptr_to_state(PARTICLEENGINE::Position, 0);
    double* vel = container->get_ptr_to_state(PARTICLEENGINE::Velocity, 0);
    double* acc = container->get_ptr_to_state(PARTICLEENGINE::Acceleration, 0);

    // get particle state dimension
    int statedim = container->get_state_dim(PARTICLEENGINE::Position);

    // safety check
    if (static_cast<std::size_t>(statedim) != function.number_components())
      FOUR_C_THROW("dimension of function defining dirichlet boundary condition not correct!");

    // iterate over owned particles of current type
    for (int i = 0; i < particlestored; ++i)
    {
      // iterate over spatial dimension
      for (int dim = 0; dim < statedim; ++dim)
      {
        // evaluate function, first and second time derivative
        functtimederiv = function.evaluate_time_derivative(
            std::span(&(refpos[statedim * i]), 3), evaltime, deg, dim);

        // set position state
        if (evalpos)
        {
          // check for periodic boundary condition in current spatial direction
          if (particleengineinterface_->have_periodic_boundary_conditions_in_spatial_direction(dim))
          {
            // length of binning domain in a spatial direction
            const double binningdomainlength =
                particleengineinterface_->length_of_binning_domain_in_a_spatial_direction(dim);

            // get displacement from reference position canceling out multiples of periodic length
            // in current spatial direction
            double displacement = std::fmod(functtimederiv[0], binningdomainlength);

            // get new position
            double newpos = refpos[statedim * i + dim] + displacement;

            // shift by periodic length if new position is close to the periodic boundary and old
            // position is on other end domain
            if ((newpos > (boundingbox(dim, 1) - binsize[dim])) and
                (std::abs(newpos - pos[statedim * i + dim]) > binsize[dim]))
              pos[statedim * i + dim] = newpos - binningdomainlength;
            else if ((newpos < (boundingbox(dim, 0) + binsize[dim])) and
                     (std::abs(newpos - pos[statedim * i + dim]) > binsize[dim]))
              pos[statedim * i + dim] = newpos + binningdomainlength;
            else
              pos[statedim * i + dim] = newpos;
          }
          // no periodic boundary conditions in current spatial direction
          else
            pos[statedim * i + dim] = refpos[statedim * i + dim] + functtimederiv[0];
        }

        // set velocity state
        if (evalvel) vel[statedim * i + dim] = functtimederiv[1];

        // set acceleration state
        if (evalacc) acc[statedim * i + dim] = functtimederiv[2];
      }
    }
  }
}

FOUR_C_NAMESPACE_CLOSE
