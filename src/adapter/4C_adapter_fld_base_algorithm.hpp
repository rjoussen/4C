// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_ADAPTER_FLD_BASE_ALGORITHM_HPP
#define FOUR_C_ADAPTER_FLD_BASE_ALGORITHM_HPP

#include "4C_config.hpp"

#include "4C_adapter_fld_fluid.hpp"
#include "4C_linalg_map.hpp"
#include "4C_linalg_vector.hpp"
#include "4C_utils_parameter_list.fwd.hpp"

#include <memory>

FOUR_C_NAMESPACE_OPEN

namespace Adapter
{
  /// fluid field solver
  class FluidBaseAlgorithm
  {
   public:
    /// constructor which distinguishes different discretizations for different fluids in
    /// multi-fluid field problems // rauch 09/13
    explicit FluidBaseAlgorithm(const Teuchos::ParameterList& prbdyn,
        const Teuchos::ParameterList& fdyn, const std::string& disname, bool isale,
        bool init = true);  // initialize time-integration scheme immediately
    // remark: parameter init allows for distinguishing an immediate initialization of all members
    // and state vectors and a
    //         later initialization which enables a later modifications of the maps

    /// second constructor (special version for turbulent flows with separate inflow
    /// section for generation of turbulent inflow profiles)
    explicit FluidBaseAlgorithm(const Teuchos::ParameterList& prbdyn,
        const std::shared_ptr<Core::FE::Discretization> discret);

    /// virtual destructor to support polymorph destruction
    virtual ~FluidBaseAlgorithm() = default;

    /// access to fluid field solver
    const std::shared_ptr<Fluid>& fluid_field() { return fluid_; }

    /// set the initial flow field in the fluid
    void set_initial_flow_field(const Teuchos::ParameterList& fdyn);

   private:
    //! setup fluid algorithm (overriding some fluid parameters with
    //! values specified in given problem-dependent ParameterList)
    /**
     * \note In this function the linear solver object is generated. For pure fluid problems or
     * fluid meshtying (no block matrix) the FLUID SOLVER block from the 4C input file is used.
     * For fluid meshtying (block matrix) the MESHTYING SOLVER block is used as main solver object
     * with a block preconditioner (BGS or SIMPLE type). The block preconditioners use the
     * information form the FLUID SOLVER and the FLUID PRESSURE SOLVER block for the velocity and
     * pressure dofs.
     */
    void setup_fluid(const Teuchos::ParameterList& prbdyn, const Teuchos::ParameterList& fdyn,
        const std::string& disname, bool isale, bool init);

    /// set the initial turbulent inflow field in the fluid
    void set_initial_inflow_field(const Teuchos::ParameterList& fdyn);

    /// setup second fluid algorithm (overriding some fluid parameters with
    /// values specified in given problem-dependent Turbulent Inflow ParameterList)
    /// separate discretization for inflow generation
    void setup_inflow_fluid(const Teuchos::ParameterList& prbdyn,
        const std::shared_ptr<Core::FE::Discretization> discret);

    //! set parameters in list required for all schemes
    void set_general_parameters(const std::shared_ptr<Teuchos::ParameterList> fluidtimeparams,
        const Teuchos::ParameterList& prbdyn, const Teuchos::ParameterList& fdyn);

    /// fluid field solver
    std::shared_ptr<Fluid> fluid_;
  };

}  // namespace Adapter

FOUR_C_NAMESPACE_CLOSE

#endif
