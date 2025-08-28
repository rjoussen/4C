// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_FLUID_TIMINT_PORO_HPP
#define FOUR_C_FLUID_TIMINT_PORO_HPP


#include "4C_config.hpp"

#include "4C_fluid_implicit_integration.hpp"

FOUR_C_NAMESPACE_OPEN

namespace FLD
{
  class TimIntPoro : public virtual FluidImplicitTimeInt
  {
   public:
    //! Standard Constructor
    TimIntPoro(const std::shared_ptr<Core::FE::Discretization>& actdis,
        const std::shared_ptr<Core::LinAlg::Solver>& solver,
        const std::shared_ptr<Teuchos::ParameterList>& params,
        const std::shared_ptr<Core::IO::DiscretizationWriter>& output, bool alefluid = false);

    /*!
    \brief initialization

    */
    void init() override;

    /*!
    \brief call elements to calculate system matrix/rhs and assemble

    */
    void assemble_mat_and_rhs() override;

    /*!
    \brief read restart data
    */
    void read_restart(int step) override;

    //! @name Set general parameter in class f3Parameter
    /*!

    \brief parameter (fix over all time step) are set in this method.
           Therefore, these parameter are accessible in the fluid element
           and in the fluid boundary element*/
    virtual void set_element_custom_parameter();

    //! set the initial porosity field
    void set_initial_porosity_field(const PoroElast::InitialField init,  //!< type of initial field
        const int startfuncno  //!< number of spatial function
        ) override;

    /*!
    \brief update iterative increment

    */
    void update_iter_incrementally(
        std::shared_ptr<const Core::LinAlg::Vector<double>> vel  //!< input residual velocities
        ) override;

    /*!
    \brief update configuration and output to file/screen

    */
    void output() override;

    /*!
    \brief Do some poro-specific stuff in assemble_mat_and_rhs

    */
    virtual void poro_int_update();

    /*!
    \brief Set custom parameters in the respective time integration class (Loma, RedModels...)

    */
    void set_custom_ele_params_assemble_mat_and_rhs(Teuchos::ParameterList& eleparams) override;

    /*!

    \brief parameter (fix over all time step) are set in this method.
           Therefore, these parameter are accessible in the fluid element
           and in the fluid boundary element

    */
    void tim_int_calculate_acceleration() override;

    /*!

    \brief parameter (fix over all time step) are set in this method.
           Therefore, these parameter are accessible in the fluid element
           and in the fluid boundary element

    */
    void set_element_general_fluid_parameter() override;

    /*!

    \brief parameter (fix over all time step) are set in this method.
           Therefore, these parameter are accessible in the fluid element
           and in the fluid boundary element

    */
    void set_element_turbulence_parameters() override;


   protected:
    //! initial porosity (poroelasticity)
    std::shared_ptr<Core::LinAlg::Vector<double>> init_porosity_field_;

   private:
  };

}  // namespace FLD


FOUR_C_NAMESPACE_CLOSE

#endif
