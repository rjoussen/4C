// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_adapter_fld_poro.hpp"

#include "4C_fem_condition_utils.hpp"
#include "4C_fem_general_assemblestrategy.hpp"
#include "4C_fluid_ele.hpp"
#include "4C_fluid_ele_action.hpp"
#include "4C_fluid_implicit_integration.hpp"
#include "4C_fluid_utils_mapextractor.hpp"
#include "4C_global_data.hpp"
#include "4C_io.hpp"
#include "4C_linalg_utils_sparse_algebra_math.hpp"

FOUR_C_NAMESPACE_OPEN

/*======================================================================*/
/* constructor */
Adapter::FluidPoro::FluidPoro(std::shared_ptr<Fluid> fluid,
    std::shared_ptr<Core::FE::Discretization> dis, std::shared_ptr<Core::LinAlg::Solver> solver,
    std::shared_ptr<Teuchos::ParameterList> params,
    std::shared_ptr<Core::IO::DiscretizationWriter> output, bool isale, bool dirichletcond)
    : Adapter::FluidFPSI::FluidFPSI(fluid, dis, solver, params, output, isale, dirichletcond)
{
  // make sure

  if (fluid_ == nullptr) FOUR_C_THROW("Failed to create the underlying fluid adapter");

  discretization()->get_condition("no_penetration", nopencond_);
}

/*======================================================================*/
/* evaluate poroelasticity specific constraint*/
void Adapter::FluidPoro::evaluate_no_penetration_cond(
    std::shared_ptr<Core::LinAlg::Vector<double>> Cond_RHS,
    std::shared_ptr<Core::LinAlg::SparseMatrix> ConstraintMatrix,
    std::shared_ptr<Core::LinAlg::SparseMatrix> struct_vel_constraint_matrix,
    std::shared_ptr<Core::LinAlg::Vector<double>> condVector, std::set<int>& condIDs,
    PoroElast::Coupltype coupltype)
{
  if (!(discretization()->filled())) FOUR_C_THROW("fill_complete() was not called");
  if (!discretization()->have_dofs()) FOUR_C_THROW("assign_degrees_of_freedom() was not called");

  discretization()->set_state(0, "dispnp", *dispnp());
  discretization()->set_state(0, "scaaf", *scaaf());

  Teuchos::ParameterList params;

  params.set("timescale", time_scaling());

  if (coupltype == PoroElast::fluidfluid)
  {
    // first, find out which dofs will be constraint
    params.set<FLD::BoundaryAction>("action", FLD::no_penetrationIDs);
    params.set<Inpar::FLUID::PhysicalType>("Physical Type", Inpar::FLUID::poro);
    discretization()->evaluate_condition(params, condVector, "no_penetration");

    // write global IDs of dofs on which the no penetration condition is applied (can vary in time
    // and iteration)
    {
      const int ndim = Global::Problem::instance()->n_dim();
      const int ndof = ndim + 1;
      const int length = condVector->local_length();
      const int nnod = length / ndof;
      const Core::LinAlg::Map& map = condVector->get_map();
      bool isset = false;
      for (int i = 0; i < nnod; i++)
      {
        isset = false;
        for (int j = 0; j < ndof; j++)
        {
          if ((*condVector).get_values()[i * ndof + j] != 0.0 and isset == false)
          {
            condIDs.insert(map.gid(i * ndof + j));
            isset = true;
            // break;
          }
          else
            (*condVector).get_values()[i * ndof + j] = 0.0;
        }
      }
    }

    // set action for elements
    params.set<FLD::BoundaryAction>("action", FLD::no_penetration);
    // params.set<std::shared_ptr< std::set<int> > >("condIDs",condIDs);
    params.set<PoroElast::Coupltype>("coupling", PoroElast::fluidfluid);
    params.set<Inpar::FLUID::PhysicalType>("Physical Type", Inpar::FLUID::poro);

    Core::FE::AssembleStrategy fluidstrategy(0,  // fluiddofset for row
        0,                                       // fluiddofset for column
        ConstraintMatrix,                        // fluid matrix
        nullptr, nullptr, nullptr, nullptr);

    discretization()->set_state(0, "condVector", *condVector);

    discretization()->evaluate_condition(params, fluidstrategy, "no_penetration");
  }
  else if (coupltype == PoroElast::fluidstructure)
  {
    discretization()->set_state(0, "velnp", *velnp());
    discretization()->set_state(0, "gridv", *grid_vel());

    discretization()->set_state(0, "condVector", *condVector);

    // set action for elements
    params.set<FLD::BoundaryAction>("action", FLD::no_penetration);
    params.set<PoroElast::Coupltype>("coupling", PoroElast::fluidstructure);
    params.set<Inpar::FLUID::PhysicalType>("Physical Type", Inpar::FLUID::poro);

    // build specific assemble strategy for the fluid-mechanical system matrix
    // from the point of view of fluid_field:
    // fluiddofset = 0, structdofset = 1
    Core::FE::AssembleStrategy couplstrategy(0,  // fluiddofset for row
        1,                                       // structdofset for column
        ConstraintMatrix,
        struct_vel_constraint_matrix,  // fluid-mechanical matrix
        Cond_RHS, nullptr, nullptr);

    // evaluate the fluid-mechanical system matrix on the fluid element
    discretization()->evaluate_condition(params, couplstrategy, "no_penetration");
  }
  else
    FOUR_C_THROW("unknown coupling type for no penetration BC");

  discretization()->clear_state();

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
std::shared_ptr<Core::LinAlg::MapExtractor> Adapter::FluidPoro::vel_pres_splitter()
{
  return fluidimpl_->vel_pres_splitter();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Adapter::FluidPoro::output(const int step, const double time)
{
  // set variables that allow the forced
  // output at this point.

  // we surely want to write now
  int upres_ = 1;
  // poro is always ALE
  bool alefluid_ = true;
  // for immersed we want to write eledata in every step
  bool write_eledata_every_step_ = true;

  // write standard output if no arguments are provided (default -1)
  if (step == -1 and time == -1.0) fluid_field()->output();
  // write extra output for specified step and time
  else
  {
    // print info to screen
    if (Core::Communication::my_mpi_rank(fluid_field()->discretization()->get_comm()) == 0)
      std::cout << "\n   Write EXTRA FLUID Output Step=" << step << " Time=" << time << " ...   \n"
                << std::endl;

    // step number and time
    fluid_field()->disc_writer()->new_step(step, time);

    // time step, especially necessary for adaptive dt
    fluid_field()->disc_writer()->write_double("timestep", fluid_field()->dt());

    // velocity/pressure vector
    fluid_field()->disc_writer()->write_vector("velnp", fluid_field()->velnp());
    // (hydrodynamic) pressure
    std::shared_ptr<Core::LinAlg::Vector<double>> pressure =
        fluid_field()->get_vel_press_splitter()->extract_cond_vector(*fluid_field()->velnp());
    fluid_field()->disc_writer()->write_vector("pressure", pressure);

    if (alefluid_) fluid_field()->disc_writer()->write_vector("dispnp", fluid_field()->dispnp());

    // write domain decomposition for visualization (only once!)
    if ((fluid_field()->step() == upres_ or fluid_field()->step() == 0) and
        !write_eledata_every_step_)
      fluid_field()->disc_writer()->write_element_data(true);
    else
      fluid_field()->disc_writer()->write_element_data(true);

    return;

  }  // write extra output
}

FOUR_C_NAMESPACE_CLOSE
