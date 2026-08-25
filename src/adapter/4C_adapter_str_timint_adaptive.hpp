// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_ADAPTER_STR_TIMINT_ADAPTIVE_HPP
#define FOUR_C_ADAPTER_STR_TIMINT_ADAPTIVE_HPP

/*----------------------------------------------------------------------*/
/* headers */
#include "4C_config.hpp"

#include "4C_adapter_str_wrapper.hpp"

FOUR_C_NAMESPACE_OPEN

// forward declarations:
namespace Solid
{
  class TimInt;
  class TimAda;
}  // namespace Solid

/*----------------------------------------------------------------------*/
/* adapting adapter */
namespace Adapter
{
  /*====================================================================*/
  /*!
   * \brief Adapter to adaptive implicit structural time integration
   *

   */
  class StructureTimIntAda : public StructureWrapper
  {
   public:
    /// Constructor
    StructureTimIntAda(std::shared_ptr<Solid::TimAda> sta, std::shared_ptr<Structure> sti);

    /// @name Time step helpers
    //@{

    /// Integrate from \f$t_1\f$ to \f$t_2\f$
    void integrate() override;

    /// prepare output (i.e. calculate stresses, strains, energies)
    void prepare_output(bool force_prepare) override;

    /// output results
    virtual void output();

    //@}

   protected:
    //! Access routines
    //{@

    std::shared_ptr<Solid::TimAda> str_ada() const { return sta_; }

    //@}

   private:
    /// the actual structure algorithm
    std::shared_ptr<Solid::TimAda> sta_;  // Solid::TimAda is the old time integration

  };  // class StructureTimIntAda

}  // namespace Adapter

/*----------------------------------------------------------------------*/
FOUR_C_NAMESPACE_CLOSE

#endif
