// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_INPAR_LEVELSET_HPP
#define FOUR_C_INPAR_LEVELSET_HPP

#include "4C_config.hpp"

#include "4C_io_input_spec.hpp"

#include <vector>

FOUR_C_NAMESPACE_OPEN

// forward declaration
namespace Core::Conditions
{
  class ConditionDefinition;
}  // namespace Core::Conditions
namespace Inpar
{
  namespace ScaTra
  {
    /// type of reinitialization for level set function
    enum ReInitialAction
    {
      reinitaction_none,
      reinitaction_signeddistancefunction,
      reinitaction_sussman,
      reinitaction_ellipticeq
    };

    /// type of signum function for reinitialization based on solving Sussman's equation
    enum SmoothedSignType
    {
      signtype_nonsmoothed,
      signtype_SussmanFatemi1999,
      signtype_SussmanSmerekaOsher1994,
      signtype_PengEtAl1999
    };

    /// form of linearization for reinitialization based on solving Sussman's equation
    enum LinReinit
    {
      newton,
      fixed_point
    };

    /// type of characteristic element length for reinitialization based on solving Sussman's
    /// equation
    enum CharEleLengthReinit
    {
      root_of_volume_reinit,
      streamlength_reinit
    };

    /// form of velocity for reinitialization based on solving Sussman's equation
    enum VelReinit
    {
      vel_reinit_integration_point_based,
      vel_reinit_node_based
    };

    /// form of artificial diffusion for reinitialization based on solving Sussman's equation
    enum ArtDiff
    {
      artdiff_none,
      artdiff_isotropic,
      artdiff_crosswind
    };

    /// type of reinitialization convergence check
    enum ReInitialStationaryCheck
    {
      reinit_stationarycheck_L1normintegrated,
      reinit_stationarycheck_numsteps
    };

    /// compute error compared to analytical solution
    enum CalcErrorLevelSet
    {
      calcerror_no_ls,
      calcerror_initial_field
    };

    /// problem dimension in case of quasi 2D
    enum LSDim
    {
      ls_3D,
      ls_2Dx,
      ls_2Dy,
      ls_2Dz
    };

    /// diffusivity for elliptic reinitialization
    enum DiffFunc
    {
      hyperbolic,
      hyperbolic_smoothed_positive,
      hyperbolic_clipped_05,
      hyperbolic_clipped_1
    };

  }  // namespace ScaTra

  namespace LevelSet
  {
    /// set the levelset parameters
    std::vector<Core::IO::InputSpec> set_valid_parameters();

    /// set specific level set conditions
    void set_valid_conditions(std::vector<Core::Conditions::ConditionDefinition>& condlist);
  }  // namespace LevelSet
}  // namespace Inpar

/*----------------------------------------------------------------------*/
FOUR_C_NAMESPACE_CLOSE

#endif
