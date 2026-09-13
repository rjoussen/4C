// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_SOLVER_NONLIN_NOX_FLOATING_POINT_EXCEPTION_HPP
#define FOUR_C_SOLVER_NONLIN_NOX_FLOATING_POINT_EXCEPTION_HPP

#include "4C_config.hpp"

#include <sstream>

FOUR_C_NAMESPACE_OPEN

namespace NOX
{
  namespace Nln
  {
    /** \brief Class which is used to disable and enable the inherent floating
     *  point exception checks.
     *
     *  This class wraps the functionality needed to avoid system abortion calls
     *  because of floating point exceptions. This is mainly used in globalization
     *  techniques where a trial point may lead to a floating point exception but
     *  is rejected afterwards.
     *
     *  */
    class FloatingPointException
    {
     public:
      /** \brief Wrapper function which checks, clears and re-enables the
       *  floating point exception checks if desired
       *
       *  */
      int postcompute(std::ostream& os) const;

      /** \brief Disable the automatic program abortion if a floating point
       *  exception is detected
       *
       *  \note This routine is only executed if the boolean \c shall_be_caught_
       *  is set to TRUE.
       *
       *  */
      void disable() const;

      /** \brief Check for any set floating point exception flags
       *
       *  */
      static int check_and_print(std::ostream& os);

      /** \brief Clear all potentially set floating point flags
       *
       *  \note This routine is only executed if the boolean \c shall_be_caught_
       *  is set to TRUE.
       *
       *  */
      void clear() const;

      /** \brief Enable the automatic program abortion if a floating point
       *  exception is detected
       *
       *  \note This routine is only executed if the boolean \c shall_be_caught_
       *  is set to TRUE.
       *
       *  */
      void enable() const;
    };
  }  // namespace Nln
}  // namespace NOX

FOUR_C_NAMESPACE_CLOSE

#endif
