// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_FLUID_ELE_BOUNDARY_PARENT_CALC_HPP
#define FOUR_C_FLUID_ELE_BOUNDARY_PARENT_CALC_HPP


#include "4C_config.hpp"

#include "4C_fem_general_utils_integration.hpp"
#include "4C_fem_general_utils_local_connectivity_matrices.hpp"
#include "4C_fluid_ele_parameter_std.hpp"
#include "4C_fluid_ele_parameter_timint.hpp"
#include "4C_linalg_serialdensevector.hpp"
#include "4C_utils_singleton_owner.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Core::FE
{
  class Discretization;
}  // namespace Core::FE

namespace Discret
{
  namespace Elements
  {
    class FluidBoundary;

    /// Interface base class for FluidBoundaryParent
    /*!
      This class exists to provide a common interface for all template
      versions of FluidBoundaryParent. The only function
      this class actually defines is Impl, which returns a pointer to
      the appropriate version of FluidBoundaryImpl.
    */
    class FluidBoundaryParentInterface
    {
     public:
      //! empty constructor
      FluidBoundaryParentInterface() {}

      //! empty destructor
      virtual ~FluidBoundaryParentInterface() = default;
      virtual void flow_dep_pressure_bc(Discret::Elements::FluidBoundary* ele1,
          Teuchos::ParameterList& params, Core::FE::Discretization& discretization,
          std::vector<int>& lm, Core::LinAlg::SerialDenseMatrix& elemat1,
          Core::LinAlg::SerialDenseVector::Base& elevec1) = 0;

      virtual void slip_supp_bc(Discret::Elements::FluidBoundary* ele1,
          Teuchos::ParameterList& params, Core::FE::Discretization& discretization,
          std::vector<int>& lm, Core::LinAlg::SerialDenseMatrix& elemat1,
          Core::LinAlg::SerialDenseVector::Base& elevec1) = 0;

      virtual void navier_slip_bc(Discret::Elements::FluidBoundary* ele1,
          Teuchos::ParameterList& params, Core::FE::Discretization& discretization,
          std::vector<int>& lm, Core::LinAlg::SerialDenseMatrix& elemat1,
          Core::LinAlg::SerialDenseVector::Base& elevec1) = 0;

      virtual void evaluate_weak_dbc(Discret::Elements::FluidBoundary* ele1,
          Teuchos::ParameterList& params, Core::FE::Discretization& discretization,
          std::vector<int>& lm, Core::LinAlg::SerialDenseMatrix& elemat1,
          Core::LinAlg::SerialDenseVector::Base& elevec1) = 0;

      virtual void estimate_nitsche_trace_max_eigenvalue(Core::Elements::FaceElement* ele1,
          Teuchos::ParameterList& params, Core::FE::Discretization& discretization,
          std::vector<int>& lm, Core::LinAlg::SerialDenseMatrix& elemat1,
          Core::LinAlg::SerialDenseMatrix& elemat2) = 0;

      virtual void mix_hyb_dirichlet(Discret::Elements::FluidBoundary* ele1,
          Teuchos::ParameterList& params, Core::FE::Discretization& discretization,
          std::vector<int>& lm, Core::LinAlg::SerialDenseMatrix& elemat1,
          Core::LinAlg::SerialDenseVector::Base& elevec1) = 0;

      /// Internal implementation class for FluidBoundaryParent elements
      static FluidBoundaryParentInterface* impl(Core::Elements::FaceElement* ele);
    };


    /// Internal FluidBoundaryParent element implementation
    /*!
      This internal class keeps all the working arrays needed to
      calculate the FluidBoundaryParent element. Additionally the method Sysmat()
      provides a clean and fast element implementation.

    */
    template <Core::FE::CellType distype>
    class FluidBoundaryParent : public FluidBoundaryParentInterface
    {
      friend class FluidEleParameter;

     public:
      //! Singleton access method
      static FluidBoundaryParentInterface* instance(
          Core::Utils::SingletonAction action = Core::Utils::SingletonAction::create);

      /// Constructor with number of nodes
      FluidBoundaryParent();

      void flow_dep_pressure_bc(Discret::Elements::FluidBoundary* surfele,
          Teuchos::ParameterList& params, Core::FE::Discretization& discretization,
          std::vector<int>& lm, Core::LinAlg::SerialDenseMatrix& elemat,
          Core::LinAlg::SerialDenseVector::Base& elevec) override;

      void slip_supp_bc(Discret::Elements::FluidBoundary* surfele, Teuchos::ParameterList& params,
          Core::FE::Discretization& discretization, std::vector<int>& lm,
          Core::LinAlg::SerialDenseMatrix& elemat,
          Core::LinAlg::SerialDenseVector::Base& elevec) override;

      void navier_slip_bc(Discret::Elements::FluidBoundary* surfele, Teuchos::ParameterList& params,
          Core::FE::Discretization& discretization, std::vector<int>& lm,
          Core::LinAlg::SerialDenseMatrix& elemat,
          Core::LinAlg::SerialDenseVector::Base& elevec) override;

      void evaluate_weak_dbc(Discret::Elements::FluidBoundary* surfele,
          Teuchos::ParameterList& params, Core::FE::Discretization& discretization,
          std::vector<int>& lm, Core::LinAlg::SerialDenseMatrix& elemat,
          Core::LinAlg::SerialDenseVector::Base& elevec) override;

      void estimate_nitsche_trace_max_eigenvalue(Core::Elements::FaceElement* surfele,
          Teuchos::ParameterList& params, Core::FE::Discretization& discretization,
          std::vector<int>& lm, Core::LinAlg::SerialDenseMatrix& elemat1,
          Core::LinAlg::SerialDenseMatrix& elemat2) override;

      void mix_hyb_dirichlet(Discret::Elements::FluidBoundary* surfele,
          Teuchos::ParameterList& params, Core::FE::Discretization& discretization,
          std::vector<int>& lm, Core::LinAlg::SerialDenseMatrix& elemat,
          Core::LinAlg::SerialDenseVector::Base& elevec) override;


     private:
      template <Core::FE::CellType bdistype, Core::FE::CellType pdistype>
      void flow_dep_pressure_bc(Discret::Elements::FluidBoundary* ele,
          Teuchos::ParameterList& params, Core::FE::Discretization& discretization,
          std::vector<int>& lm, Core::LinAlg::SerialDenseMatrix& elemat,
          Core::LinAlg::SerialDenseVector::Base& elevec);

      template <Core::FE::CellType bdistype, Core::FE::CellType pdistype>
      void slip_supp_bc(Discret::Elements::FluidBoundary* ele, Teuchos::ParameterList& params,
          Core::FE::Discretization& discretization, std::vector<int>& lm,
          Core::LinAlg::SerialDenseMatrix& elemat, Core::LinAlg::SerialDenseVector::Base& elevec);

      template <Core::FE::CellType bdistype, Core::FE::CellType pdistype>
      void navier_slip_bc(Discret::Elements::FluidBoundary* ele, Teuchos::ParameterList& params,
          Core::FE::Discretization& discretization, std::vector<int>& lm,
          Core::LinAlg::SerialDenseMatrix& elemat, Core::LinAlg::SerialDenseVector::Base& elevec);

      template <Core::FE::CellType bdistype, Core::FE::CellType pdistype>
      void evaluate_weak_dbc(Discret::Elements::FluidBoundary* ele, Teuchos::ParameterList& params,
          Core::FE::Discretization& discretization, std::vector<int>& lm,
          Core::LinAlg::SerialDenseMatrix& elemat, Core::LinAlg::SerialDenseVector::Base& elevec);

      template <Core::FE::CellType bdistype, Core::FE::CellType pdistype>
      void estimate_nitsche_trace_max_eigenvalue(Core::Elements::FaceElement* ele,
          Teuchos::ParameterList& params, Core::FE::Discretization& discretization,
          std::vector<int>& lm, Core::LinAlg::SerialDenseMatrix& elemat1,
          Core::LinAlg::SerialDenseMatrix& elemat2);

      template <Core::FE::CellType bdistype, Core::FE::CellType pdistype>
      void mix_hyb_dirichlet(Discret::Elements::FluidBoundary* ele, Teuchos::ParameterList& params,
          Core::FE::Discretization& discretization, std::vector<int>& lm,
          Core::LinAlg::SerialDenseMatrix& elemat, Core::LinAlg::SerialDenseVector::Base& elevec);

      void get_density_and_viscosity(std::shared_ptr<const Core::Mat::Material> material,
          const double pscaaf, const double thermpressaf, const double rateofstrain);


      //! pointer to parameter list
      Discret::Elements::FluidEleParameter* fldpara_;
      Discret::Elements::FluidEleParameterTimInt* fldparatimint_;


      //! infinitesimal area element drs
      double drs_;
      //! integration factor
      double fac_;
      //! physical viscosity
      double visc_;
      //! density at t_(n+alpha_F) or t_(n+1)
      double densaf_;

      //-----------------------------------------------------------------
      //-----------------------------------------------------------------
      //
      //                    SPALDINGS LAW OF THE WALL
      //
      //-----------------------------------------------------------------
      //-----------------------------------------------------------------
      //-----------------------------------------------------------------
      //  evaluate the residual of Spaldings law of the wall
      //                                           (private) gammi 11/09
      //-----------------------------------------------------------------
      double spalding_residual(
          const double y, const double visc, const double tau_B, const double normu)
      {
        // get dimensionless velocity
        const double up = sqrt(normu / tau_B);

        // constants
        const double chi = 0.4;
        const double B = 5.5;

        //      +
        // get y , a dimensionless boundary layer thickness
        const double yp = y * sqrt(tau_B * normu) / visc;

        return (yp - (up + exp(-chi * B) *
                               (exp(chi * up) - 1.0 -
                                   chi * up * (1 + chi * up / 2.0 * (1.0 + chi * up / 3.0)))));
      }

      //-----------------------------------------------------------------
      //  evaluate the residual of Spaldings law of the wall
      //                                           (private) gammi 11/09
      //-----------------------------------------------------------------
      double spalding_residual_utau(
          const double y, const double visc, const double utau, const double normu)
      {
        // get dimensionless velocity
        if (abs(utau) < 1.0E-14) FOUR_C_THROW("utau is zero!");
        const double up = normu / utau;

        // constants
        const double chi = 0.4;
        const double B = 5.5;

        //      +
        // get y , a dimensionless boundary layer thickness
        if (visc < 1.0E-14) FOUR_C_THROW("visc is zero or negative!");
        const double yp = y * utau / visc;

        return (yp - (up + exp(-chi * B) *
                               (exp(chi * up) - 1.0 -
                                   chi * up * (1 + chi * up / 2.0 * (1.0 + chi * up / 3.0)))));
      }

      //-----------------------------------------------------------------
      //     evaluate the residual of Spaldings law of the wall
      //                                           (private) gammi 11/09
      //-----------------------------------------------------------------
      double jacobian_spalding_residual(
          const double y, const double visc, const double tau_B, const double normu)
      {
        // constants
        const double chi = 0.4;
        const double B = 5.5;

        // get dimensionless velocity
        const double up = sqrt(normu / tau_B);

        // compute the derivative of the Spalding residual w.r.t. tau_B
        double drdtauB = y / (2.0 * visc * sqrt(tau_B)) * sqrt(normu);

        drdtauB +=
            (1 + chi * exp(-chi * B) * (exp(chi * up) - 1.0 - chi * up * (1.0 + 0.5 * chi * up))) *
            0.5 * sqrt(normu) / (sqrt(tau_B) * sqrt(tau_B) * sqrt(tau_B));

        return (drdtauB);
      }

      //-----------------------------------------------------------------
      //     evaluate the residual of Spaldings law of the wall
      //                                           (private) gammi 11/09
      //-----------------------------------------------------------------
      double jacobian_spalding_residual_utau(
          const double y, const double visc, const double utau, const double normu)
      {
        // get dimensionless velocity
        const double up = normu / utau;

        // constants
        const double chi = 0.4;
        const double B = 5.5;

        //                    +
        // get derivative of y , a dimensionless boundary layer thickness
        const double dyplus_dutau = y / visc;

        //                                   +
        // get derivative of function f wrt u
        const double df_duplus = 1 + exp(-chi * B) * (chi * exp(chi * up) - chi - chi * chi * up -
                                                         chi * chi * chi * up * up / 2.0);

        //                    +
        // get derivative of u  wrt u_tau
        const double duplus_dutau = -normu / (utau * utau);

        return (dyplus_dutau - df_duplus * duplus_dutau);
      }

      //-----------------------------------------------------------------
      //     evaluate the residual of Spaldings law of the wall
      //                                           (private) gammi 11/09
      //-----------------------------------------------------------------
      double jacobian_spalding_residual_u(
          const double y, const double visc, const double utau, const double normu)
      {
        // get dimensionless velocity
        const double up = normu / utau;

        // constants
        const double chi = 0.4;
        const double B = 5.5;

        //                    +
        // get derivative of y , a dimensionless boundary layer thickness
        const double dyplus_du = 0.0;

        //                                   +
        // get derivative of function f wrt u
        const double df_duplus = 1 + exp(-chi * B) * (chi * exp(chi * up) - chi - chi * chi * up -
                                                         chi * chi * chi * up * up / 2.0);

        //                    +
        // get derivative of u  wrt u
        const double duplus_du = 1.0 / utau;

        return (dyplus_du - df_duplus * duplus_du);
      }

      //-----------------------------------------------------------------
      //     evaluate the residual of Spaldings law of the wall
      //                                           (private) gammi 11/09
      //-----------------------------------------------------------------
      double jacobian_spalding_residual_uplus(
          const double y, const double visc, const double utau, const double normu)
      {
        // get dimensionless velocity
        const double up = normu / utau;

        // constants
        const double chi = 0.4;
        const double B = 5.5;

        //                    +
        // get derivative of y , a dimensionless boundary layer thickness
        const double dyplus_duplus = 0.0;
        /*
        // Spaldings law of the wall
        //                     /                                                  \
        //                    |                           /     +\ 2    /     +\ 3 |
        //                    |       +                  | chi*u  |    | chi*u  |  |
        //  +   +    -chi*B   |  chi*u               +    \      /      \      /   |
        // y - u  - e       * | e       - 1.0 - chi*u  - ----------- - ----------- | = 0
        //                    |                              2.0           6.0     |
        //                    |                                                    |
        //                     \                                                  /
        */

        //                                   +
        // get derivative of function f wrt u
        const double df_duplus = 1 + exp(-chi * B) * (chi * exp(chi * up) - chi - chi * chi * up -
                                                         chi * chi * chi * up * up / 2.0);

        return (dyplus_duplus - df_duplus);
      }
    };  // end class fluid3boundaryImpl
  }  // namespace Elements
}  // namespace Discret

FOUR_C_NAMESPACE_CLOSE

#endif
