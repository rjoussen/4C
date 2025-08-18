// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_CONTACT_MANAGER_HPP
#define FOUR_C_CONTACT_MANAGER_HPP

#include "4C_config.hpp"

#include "4C_mortar_manager_base.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Core::FE
{
  class Discretization;
}  // namespace Core::FE

namespace Core::LinAlg
{
  template <typename T>
  class Vector;
}

namespace CONTACT
{
  // forward declarations
  class Element;

  /*!
  \brief 4C implementation of main class to control all contact

  */
  class Manager : public Mortar::ManagerBase
  {
   public:
    /*!
    \brief Standard Constructor

    The constructor takes a discretization  that is expected to have at least
    two contact boundary conditions. It extracts all contact boundary conditions
    and constructs one or multiple contact interfaces from them and stores them.

    All interfaces will be fill_complete in order to get their internal discretization ready for
    further usage. This step also takes care of extending the interface ghosting depending on the
    user's choice.

    In addition, it creates the necessary solver strategy object which handles
    the whole contact evaluation.

    \param discret (in): A discretization containing contact boundary conditions
    \param alphaf (in): Generalized-alpha parameter (set to 0.0 by default)

    */
    Manager(Core::FE::Discretization& discret, double alphaf = 0.0);

    //! @name Access methods
    //! @{

    /*!
    \brief Get discretization

    */
    const Core::FE::Discretization& discret() const { return discret_; };

    //! @}

    //! @name Evaluation methods
    //! @{

    /*!
    \brief Write restart information for contact

    The additionally necessary restart information in the contact
    case are the current Lagrange multiplier values and the current
    active set status of each slave node.

    \param output (in): IO::discretization writer for restart
    \param forcedrestart (in): Force writing of restart data?
    */
    void write_restart(Core::IO::DiscretizationWriter& output, bool forcedrestart = false) final;

    /*!
    \brief Read restart information for contact

    This method has the inverse functionality of write_restart, as
    it reads the activetoggle / lmold vectors and stores the restart
    status into each slave node. Moreover, all global maps concerning
    the active set and the old mortar matrices D,M are rebuilt based
    on the restart information.

    \param reader (in): IO::discretization reader for restart
    \param dis (in)   : global dof displacement vector
    \param zero (in)  : global dof zero vector

    */
    void read_restart(Core::IO::DiscretizationReader& reader,
        std::shared_ptr<Core::LinAlg::Vector<double>> dis,
        std::shared_ptr<Core::LinAlg::Vector<double>> zero) final;

    /*!
    \brief Write interface quantities for postprocessing

    \param output (in): IO::discretization writer for restart

    */
    void postprocess_quantities(Core::IO::DiscretizationWriter& output) final;

    //! [derived]
    void postprocess_quantities_per_interface(
        std::shared_ptr<Teuchos::ParameterList> outputParams) final;

    /*!
    \brief Reconnect Contact Element -- Parent Element Pointers

    As during the Restart the initial created structural elements are destructed and created again,
    the pointer of these elements changes and therefore needs to be reconnected
    */
    void reconnect_parent_elements();

    /*!
    \brief Set Parent Elements for Poro Face Elements

    \param slavetype type of slave elements --> = (-1; //1 poro, 0 struct, -1 default)
    \param mastertype type of master elements --> = (-1; //1 poro, 0 struct, -1 default)
    \param[out] cele Reference to pointer of contact face element
    \param[out] ele Reference to pointer of contact parent element

    */
    void set_poro_parent_element(int& slavetype, int& mastertype, CONTACT::Element& cele,
        std::shared_ptr<Core::Elements::Element> ele);

    /*!
    \brief Find Physical Type (Poro or Structure) of Poro Interface

    \param poromaster ??
    \param poroslave ??
    \param structmaster ??
    \param structslave ??
    \param slavetype ??
    \param mastertype ??
    */
    void find_poro_interface_types(bool& poromaster, bool& poroslave, bool& structmaster,
        bool& structslave, int& slavetype, int& mastertype) const;

    //! @}

   protected:
    //! the underlying problem discretization
    const Core::FE::Discretization& discret_;

   private:
    /*!
     \brief Read and check contact input parameters

     All specified contact-related input parameters are read from the
     Global::Problem::instance() and stored into a local variable of
     type Teuchos::ParameterList. Invalid parameter combinations are
     sorted out and throw a FOUR_C_THROW.

     */
    bool read_and_check_input(Teuchos::ParameterList& cparams) const;

    //! don't want operator=
    Manager operator=(const Manager& old) = delete;

    //! don't want copy constructor
    Manager(const Manager& old) = delete;

  };  // class Manager
}  // namespace CONTACT

FOUR_C_NAMESPACE_CLOSE

#endif
