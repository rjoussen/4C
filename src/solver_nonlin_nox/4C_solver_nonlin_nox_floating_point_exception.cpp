// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_solver_nonlin_nox_floating_point_exception.hpp"

#include <NOX_Utils.H>
#ifdef FOUR_C_ENABLE_FE_TRAPPING
#include <fenv.h>
#endif

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::Nln::FloatingPointException::disable() const
{
#ifdef FOUR_C_ENABLE_FE_TRAPPING
  fedisableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);
#endif
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::Nln::FloatingPointException::clear() const
{
#ifdef FOUR_C_ENABLE_FE_TRAPPING
  feclearexcept(FE_ALL_EXCEPT);
#endif
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
int NOX::Nln::FloatingPointException::check_and_print(std::ostream& os)
{
  int exception_occurred = 0;
#ifdef FOUR_C_ENABLE_FE_TRAPPING
  if (fetestexcept(FE_INVALID))
  {
    os << "Floating Point Exception: INVALID\n";
    ++exception_occurred;
  }
  if (fetestexcept(FE_OVERFLOW))
  {
    os << "Floating Point Exception: OVERFLOW\n";
    ++exception_occurred;
  }
  if (fetestexcept(FE_DIVBYZERO))
  {
    os << "Floating Point Exception: DIVISION BY ZERO\n";
    ++exception_occurred;
  }
#endif
  return exception_occurred;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::Nln::FloatingPointException::enable() const
{
#ifdef FOUR_C_ENABLE_FE_TRAPPING
  feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);
#endif
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
int NOX::Nln::FloatingPointException::postcompute(std::ostream& os) const
{
  const int err = check_and_print(os);
  if (err)
    os << ::NOX::Utils::fill(40, '-') << "\n"
       << "Caught floating point exceptions = " << err << std::endl;

  clear();
  enable();

  return err;
}

FOUR_C_NAMESPACE_CLOSE
