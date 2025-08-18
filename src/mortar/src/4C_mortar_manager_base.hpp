// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_MORTAR_MANAGER_BASE_HPP
#define FOUR_C_MORTAR_MANAGER_BASE_HPP

#include "4C_config.hpp"

#include "4C_linalg_map.hpp"
#include "4C_linalg_vector.hpp"
#include "4C_utils_parameter_list.fwd.hpp"

#include <memory>

FOUR_C_NAMESPACE_OPEN

// forward declarations
namespace Core::IO
{
  class DiscretizationWriter;
  class DiscretizationReader;
}  // namespace Core::IO

namespace Core::FE
{
  class Discretization;
}  // namespace Core::FE

namespace Core::Elements
{
  class Element;
}

namespace Mortar
{
  // forward declarations
  class StrategyBase;
  class Node;
  class Element;

  /*!
  \brief Abstract base class to control all mortar coupling

  */
  class ManagerBase
  {
   public:
    //! @name Enums and Friends
    //@}

    /*!
    \brief Standard Constructor

    The base class constructor is empty.

    One needs a derived class for a concrete implementation of the Manager
    class into a given FE code environment (see e.g. contact_manager.H and
    contact_manager.cpp for the 4C mortar contact implementation or
    contact_meshtying_manager.H and meshtying_manager.coo for the 4C mortar
    meshtying implementation).

    This constructor then has to be fed with a discretization that is expected
    to carry at least two mortar boundary conditions (one is only sufficient
    in the case of self contact simulations). It extracts all mortar boundary
    conditions, constructs one or multiple mortar interfaces and stores them.

    It also builds the corresponding strategy solver object and stores a
    reference in the strategy_ member variable.

    */
    ManagerBase();

    /*!
    \brief Destructor

    */
    virtual ~ManagerBase() = default;

    //! @name Access methods

    /*!
    \brief Get communicator

    */
    MPI_Comm get_comm() const { return comm_; }

    /*!
    \brief Return the object for the solving strategy.

    All necessary steps for the computation algorithm
    have to be specialized in subclasses of StrategyBase

    */
    Mortar::StrategyBase& get_strategy() { return *strategy_; }

    //@}

    //! @name Purely virtual functions
    //! @{

    //! Write interface quantities for postprocessing
    virtual void postprocess_quantities(Core::IO::DiscretizationWriter& output) = 0;

    /*!
    \brief Write results for visualization separately for each meshtying/contact interface

    Call each interface, such that each interface can handle its own output of results.

    \param[in] outputParams Parameter list with stuff required by interfaces to write output
    */
    virtual void postprocess_quantities_per_interface(
        std::shared_ptr<Teuchos::ParameterList> outputParams) = 0;

    //! Read restart data from disk
    virtual void read_restart(Core::IO::DiscretizationReader& reader,
        std::shared_ptr<Core::LinAlg::Vector<double>> dis,
        std::shared_ptr<Core::LinAlg::Vector<double>> zero) = 0;

    //! Write restart data to disk
    virtual void write_restart(
        Core::IO::DiscretizationWriter& output, bool forcedrestart = false) = 0;

    //! @}

   protected:
    // don't want cctor (= operator impossible anyway for abstract class)
    ManagerBase(const ManagerBase& old) = delete;

    //! Communicator
    MPI_Comm comm_;

    //! Strategy object
    std::shared_ptr<Mortar::StrategyBase> strategy_;

  };  // class ManagerBase
}  // namespace Mortar


FOUR_C_NAMESPACE_CLOSE

#endif
