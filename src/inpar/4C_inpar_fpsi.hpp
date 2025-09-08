// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_INPAR_FPSI_HPP
#define FOUR_C_INPAR_FPSI_HPP

#include "4C_config.hpp"

#include "4C_io_input_spec.hpp"

#include <vector>

FOUR_C_NAMESPACE_OPEN

// forward declaration
namespace Core::Conditions
{
  class ConditionDefinition;
}  // namespace Core::Conditions
enum FpsiCouplingType
{
  fpsi_monolithic_plain,
  partitioned
};

namespace Inpar
{
  namespace FPSI
  {
    // type of partitioned coupling for FPSI problems
    enum PartitionedCouplingMethod
    {
      RobinNeumann,
      monolithic,
      nocoupling
    };

    // type of norm to check for convergence
    enum ConvergenceNorm
    {
      absoluteconvergencenorm,            // compare absolute value with single tolerance
      absoluteconvergencenorm_sys_split,  // compare absolute value with correction of systemsize
                                          // with different tolerances for each field
      relativconvergencenorm_sys  // compare relative value with correction of systemsize with
                                  // single tolerance
    };

    // type of norm to check for convergence
    enum BinaryOp
    {
      bop_and,
      bop_or
    };

    enum FluidFieldHierarchy
    {
      fluid,
      porofluid
    };

    /// set the fpsi parameters
    Core::IO::InputSpec set_valid_parameters();

    /// set specific fpsi conditions
    void set_valid_conditions(std::vector<Core::Conditions::ConditionDefinition>& condlist);

  }  // namespace FPSI

}  // namespace Inpar

/*----------------------------------------------------------------------*/
FOUR_C_NAMESPACE_CLOSE

#endif
