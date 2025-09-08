// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_INPAR_FSI_HPP
#define FOUR_C_INPAR_FSI_HPP

#include "4C_config.hpp"

#include "4C_io_input_spec.hpp"

#include <vector>

FOUR_C_NAMESPACE_OPEN

// forward declaration
namespace Core::Conditions
{
  class ConditionDefinition;
}

/*----------------------------------------------------------------------*/
/* The coupling methods for FSI. */
/*----------------------------------------------------------------------*/
enum FsiCoupling
{
  fsi_coupling_freesurface = -1,
  fsi_coupling_undefined = 0,
  fsi_iter_stagg_fixed_rel_param = 4,  /*!< fixed-point solver with fixed relaxation parameter */
  fsi_iter_stagg_AITKEN_rel_param = 5, /*!< fixed-point solver with Aitken relaxation parameter */
  fsi_iter_stagg_steep_desc = 6, /* fixed-point solver with steepest descent relaxation parameter */
  fsi_iter_stagg_CHEB_rel_param = 7,
  fsi_iter_stagg_AITKEN_rel_force = 8,
  fsi_iter_stagg_steep_desc_force = 9,
  fsi_iter_stagg_Newton_FD = 10,
  fsi_iter_stagg_Newton_I = 11,
  fsi_iter_monolithicfluidsplit = 13,
  fsi_iter_monolithicstructuresplit,
  fsi_iter_mortar_monolithicstructuresplit,
  fsi_iter_mortar_monolithicfluidsplit,
  fsi_iter_xfem_monolithic,
  fsi_iter_stagg_NLCG, /*!< nonlinear CG solver (pretty much steepest descent with finite difference
                          Jacobian */
  fsi_iter_stagg_MFNK_FD,  /*!< matrix free Newton Krylov with finite difference Jacobian */
  fsi_iter_stagg_MFNK_FSI, /*!< matrix free Newton Krylov with FSI specific Jacobian */
  fsi_iter_stagg_MPE,      /*!< minimal polynomial extrapolation */
  fsi_iter_fluidfluid_monolithicstructuresplit,
  fsi_iter_fluidfluid_monolithicfluidsplit,
  fsi_iter_fluidfluid_monolithicstructuresplit_nonox,
  fsi_iter_fluidfluid_monolithicfluidsplit_nonox,
  fsi_iter_sliding_monolithicfluidsplit,
  fsi_iter_sliding_monolithicstructuresplit,
  fsi_iter_mortar_monolithicfluidsplit_saddlepoint
};

namespace Inpar
{
  namespace FSI
  {
    /// Type of partitioned coupling for FSI problems
    enum PartitionedCouplingMethod
    {
      DirichletNeumannSlideale,
      DirichletNeumann,
      DirichletNeumannVolCoupl
    };

    /// Linear preconditioning algorithm for monolithic block system
    enum LinearBlockSolver
    {
      PreconditionedKrylov,  ///< BGS(AMG)
      LinalgSolver           ///< use Core::LinAlg::Solver interface
    };

    /// Projection methods for sliding fluid-structure interface
    enum SlideALEProj
    {
      ALEprojection_none,
      ALEprojection_curr,
      ALEprojection_ref,
      ALEprojection_rot_z,
      ALEprojection_rot_zsphere
    };

    /// Preconditioner for constraint fsi problem
    enum PrecConstr
    {
      Simple,
      Simplec
    };

    /// @name Solution technique and related

    /// type of norm to check for convergence of newton loop
    enum ConvNorm
    {
      convnorm_abs,  ///< absolute norm
      convnorm_rel,  ///< relative norm
      convnorm_mix   ///< mixed absolute-relative norm
    };

    /// type of norm to check for convergence
    enum BinaryOp
    {
      bop_and  ///<  and
    };

    /// type of fluid field auxiliary time integrator for time adaptivity
    enum FluidMethod
    {
      timada_fld_none,            ///< no time adaptivity based on fluid field
      timada_fld_expleuler,       ///< use explicit Euler as auxiliary time integrator
      timada_fld_adamsbashforth2  ///< use Adams-Bashforth-2 as auxiliary time integrator
    };

    /// Handling of non-converged nonlinear solver
    enum DivContAct
    {
      divcont_stop,        ///< abort simulation
      divcont_continue,    ///< continue nevertheless
      divcont_halve_step,  ///< halve time step and carry on with simulation
      divcont_revert_dt    ///< revert time step size to previous one
    };

    /// Verbosity of FSI algorithm
    enum Verbosity
    {
      verbosity_subproblem,  ///< output for FSI as a subproblem
      verbosity_low,         ///< write only nonlinear solver status
      verbosity_medium,      ///< write nonlinear and linear solver status
      verbosity_full,        ///< write everything
    };

    /// Coupling variable for Dirichlet-Neumann algorithms
    enum CoupVarPart
    {
      disp,  /*!< displacement coupling*/
      force, /*!< force coupling*/
      vel    /*!<velocity coupling*/
    };

    //@}

    /// set the fsi parameters
    std::vector<Core::IO::InputSpec> set_valid_parameters();

    /// set specific fsi conditions
    void set_valid_conditions(std::vector<Core::Conditions::ConditionDefinition>& condlist);

  }  // namespace FSI

}  // namespace Inpar

/*----------------------------------------------------------------------*/
FOUR_C_NAMESPACE_CLOSE

#endif
