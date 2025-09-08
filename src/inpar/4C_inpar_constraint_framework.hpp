// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_INPAR_CONSTRAINT_FRAMEWORK_HPP
#define FOUR_C_INPAR_CONSTRAINT_FRAMEWORK_HPP

#include "4C_config.hpp"

#include "4C_io_input_spec.hpp"


FOUR_C_NAMESPACE_OPEN

namespace Inpar::Constraints
{
  /// type of the submodel for constraintmodels
  enum class SubModelType
  {
    submodel_undefined,    ///< default
    submodel_pbc_rve,      ///< apply periodic displacement bcs on
    submodel_embeddedmesh  ///< apply embedded mesh bcs
  };

  /// type of employed constraint enforcement strategy
  enum class Strategy
  {
    penalty_regularization  ///< penalty method
  };

  enum class EmbeddedMeshCouplingStrategy
  {
    //! Default value
    none,
    //! Mortar method
    mortar
  };

  enum class EmbeddedMeshConstraintEnforcement
  {
    //! Default value
    none,
    //! Penalty method
    penalty
  };

  /**
   * \brief Shape function for the mortar Lagrange-multiplicators for solid to solid embedded
   * coupling
   */
  enum class SolidToSolidMortarShapefunctions
  {
    //! Default value.
    none,
    //! Linear Lagrange elements.
    quad4,
    //! Quadratic Lagrange elements.
    quad9,
    //! Quadratic NURBS elements.
    nurbs9
  };

  /**
  \brief Set constraint parameters
  */
  Core::IO::InputSpec set_valid_parameters();
}  // namespace Inpar::Constraints

FOUR_C_NAMESPACE_CLOSE

#endif
