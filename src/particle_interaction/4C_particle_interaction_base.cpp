// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_particle_interaction_base.hpp"

#include "4C_comm_mpi_utils.hpp"
#include "4C_particle_engine_container.hpp"
#include "4C_particle_engine_interface.hpp"
#include "4C_particle_interaction_material_handler.hpp"
#include "4C_particle_interaction_runtime_writer.hpp"

FOUR_C_NAMESPACE_OPEN

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
ParticleInteraction::ParticleInteractionBase::ParticleInteractionBase(
    MPI_Comm comm, const Teuchos::ParameterList& params)
    : comm_(comm),
      myrank_(Core::Communication::my_mpi_rank(comm)),
      params_(params),
      time_(0.0),
      dt_(0.0)
{
  // empty constructor
}

void ParticleInteraction::ParticleInteractionBase::init()
{
  // init particle material handler
  init_particle_material_handler();

  // init particle interaction writer
  init_particle_interaction_writer();
}

void ParticleInteraction::ParticleInteractionBase::setup(
    const std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface,
    const std::shared_ptr<PARTICLEWALL::WallHandlerInterface> particlewallinterface)
{
  // set interface to particle engine
  particleengineinterface_ = particleengineinterface;

  // set particle container bundle
  particlecontainerbundle_ = particleengineinterface_->get_particle_container_bundle();

  // set interface to particle wall handler
  particlewallinterface_ = particlewallinterface;

  // setup particle material handler
  particlematerial_->setup();

  // setup particle interaction writer
  particleinteractionwriter_->setup();

  // init vector
  gravity_.resize(3, 0.0);
}

void ParticleInteraction::ParticleInteractionBase::write_restart() const
{
  // nothing to do
}

void ParticleInteraction::ParticleInteractionBase::read_restart(
    const std::shared_ptr<Core::IO::DiscretizationReader> reader)
{
  // read restart of particle interaction writer
  particleinteractionwriter_->read_restart(reader);
}

void ParticleInteraction::ParticleInteractionBase::
    check_particle_interaction_distance_concerning_bin_size() const
{
  // get maximum particle interaction distance
  double allprocmaxinteractiondistance = 0.0;
  double maxinteractiondistance = max_interaction_distance();
  allprocmaxinteractiondistance = Core::Communication::max_all(maxinteractiondistance, comm_);

  // bin size safety check
  if (allprocmaxinteractiondistance > particleengineinterface_->min_bin_size())
    FOUR_C_THROW("the particle interaction distance is larger than the minimal bin size ({} > {})!",
        allprocmaxinteractiondistance, particleengineinterface_->min_bin_size());

  // periodic length safety check
  if (particleengineinterface_->have_periodic_boundary_conditions())
  {
    // loop over all spatial directions
    for (int dim = 0; dim < 3; ++dim)
    {
      // check for periodic boundary condition in current spatial direction
      if (not particleengineinterface_->have_periodic_boundary_conditions_in_spatial_direction(dim))
        continue;

      // check periodic length in current spatial direction
      if ((2.0 * allprocmaxinteractiondistance) >
          particleengineinterface_->length_of_binning_domain_in_a_spatial_direction(dim))
        FOUR_C_THROW(
            "particles are not allowed to interact directly and across the periodic boundary!");
    }
  }
}

void ParticleInteraction::ParticleInteractionBase::set_current_time(const double currenttime)
{
  time_ = currenttime;
}

void ParticleInteraction::ParticleInteractionBase::set_current_step_size(
    const double currentstepsize)
{
  dt_ = currentstepsize;
}

void ParticleInteraction::ParticleInteractionBase::set_current_write_result_flag(
    bool writeresultsthisstep)
{
  // set current write result flag in particle interaction writer
  particleinteractionwriter_->set_current_write_result_flag(writeresultsthisstep);
}

void ParticleInteraction::ParticleInteractionBase::set_gravity(std::vector<double>& gravity)
{
  gravity_ = gravity;
}

void ParticleInteraction::ParticleInteractionBase::write_interaction_runtime_output(
    const int step, const double time)
{
  // write particle interaction runtime output
  particleinteractionwriter_->write_particle_interaction_runtime_output(step, time);
}

void ParticleInteraction::ParticleInteractionBase::init_particle_material_handler()
{
  // create particle material handler
  particlematerial_ = std::make_shared<ParticleInteraction::MaterialHandler>(params_);

  // init particle material handler
  particlematerial_->init();
}

void ParticleInteraction::ParticleInteractionBase::init_particle_interaction_writer()
{
  // create particle interaction writer
  particleinteractionwriter_ =
      std::make_shared<ParticleInteraction::InteractionWriter>(comm_, params_);

  // init particle interaction writer
  particleinteractionwriter_->init();
}

double ParticleInteraction::ParticleInteractionBase::max_particle_radius() const
{
  // init value of maximum radius
  double maxrad = 0.0;

  // iterate over particle types
  for (const auto& type_i : particlecontainerbundle_->get_particle_types())
  {
    // get container of owned particles of current particle type
    PARTICLEENGINE::ParticleContainer* container =
        particlecontainerbundle_->get_specific_container(type_i, PARTICLEENGINE::Owned);

    // get maximum stored value of state
    double currmaxrad = container->get_max_value_of_state(PARTICLEENGINE::Radius);

    // compare to current maximum
    maxrad = std::max(maxrad, currmaxrad);
  }

  return maxrad;
}

FOUR_C_NAMESPACE_CLOSE
