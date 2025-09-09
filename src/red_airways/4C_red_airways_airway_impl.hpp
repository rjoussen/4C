// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_RED_AIRWAYS_AIRWAY_IMPL_HPP
#define FOUR_C_RED_AIRWAYS_AIRWAY_IMPL_HPP

#include "4C_config.hpp"

#include "4C_fem_general_utils_local_connectivity_matrices.hpp"
#include "4C_linalg_fixedsizematrix.hpp"
#include "4C_red_airways_elementbase.hpp"

FOUR_C_NAMESPACE_OPEN


namespace Discret
{
  namespace Elements
  {
    /// Interface base class for airway_impl
    /*!
      This class exists to provide a common interface for all template
      versions of airway_impl. The only function
      this class actually defines is Impl, which returns a pointer to
      the appropriate version of airway_impl.
     */
    class RedAirwayImplInterface
    {
     public:
      /// Empty constructor
      RedAirwayImplInterface() {}
      /// Empty destructor
      virtual ~RedAirwayImplInterface() = default;  /// Evaluate the element
      virtual int evaluate(RedAirway* ele, Teuchos::ParameterList& params,
          Core::FE::Discretization& discretization, std::vector<int>& lm,
          Core::LinAlg::SerialDenseMatrix& elemat1, Core::LinAlg::SerialDenseMatrix& elemat2,
          Core::LinAlg::SerialDenseVector& elevec1, Core::LinAlg::SerialDenseVector& elevec2,
          Core::LinAlg::SerialDenseVector& elevec3, std::shared_ptr<Core::Mat::Material> mat) = 0;

      virtual void initial(RedAirway* ele, Teuchos::ParameterList& params,
          Core::FE::Discretization& discretization, std::vector<int>& lm,
          Core::LinAlg::SerialDenseVector& radii_in, Core::LinAlg::SerialDenseVector& radii_out,
          std::shared_ptr<const Core::Mat::Material> material) = 0;

      virtual void evaluate_terminal_bc(RedAirway* ele, Teuchos::ParameterList& params,
          Core::FE::Discretization& discretization, std::vector<int>& lm,
          Core::LinAlg::SerialDenseVector& elevec1, std::shared_ptr<Core::Mat::Material> mat) = 0;

      virtual void calc_flow_rates(RedAirway* ele, Teuchos::ParameterList& params,
          Core::FE::Discretization& discretization, std::vector<int>& lm,
          std::shared_ptr<Core::Mat::Material> mat) = 0;

      virtual void calc_elem_volume(RedAirway* ele, Teuchos::ParameterList& params,
          Core::FE::Discretization& discretization, std::vector<int>& lm,
          std::shared_ptr<Core::Mat::Material> mat) = 0;

      virtual void get_coupled_values(RedAirway* ele, Teuchos::ParameterList& params,
          Core::FE::Discretization& discretization, std::vector<int>& lm,
          std::shared_ptr<Core::Mat::Material> material) = 0;

      /// Internal implementation class for airway element
      static RedAirwayImplInterface* impl(Discret::Elements::RedAirway* airway);
    };


    /// Internal airway implementation
    /*!
      This internal class keeps all the working arrays needed to
      calculate the airway element. Additionally the method Sysmat()
      provides a clean and fast element implementation.

      <h3>Purpose</h3>

    */

    template <Core::FE::CellType distype>
    class AirwayImpl : public RedAirwayImplInterface
    {
     public:
      /// Default constructor. Will need to query data from ParameterList.
      AirwayImpl() = default;

      //! number of nodes
      static constexpr int iel = Core::FE::num_nodes(distype);


      /// Evaluate
      /*!
        The evaluate function for the general airway case.
       */
      int evaluate(RedAirway* ele, Teuchos::ParameterList& params,
          Core::FE::Discretization& discretization, std::vector<int>& lm,
          Core::LinAlg::SerialDenseMatrix& elemat1, Core::LinAlg::SerialDenseMatrix& elemat2,
          Core::LinAlg::SerialDenseVector& elevec1, Core::LinAlg::SerialDenseVector& elevec2,
          Core::LinAlg::SerialDenseVector& elevec3,
          std::shared_ptr<Core::Mat::Material> mat) override;

      /*!
        \brief Calculate virtual trajectory xnp and state of airway (open/closed)

        \param ele              (i) the element those matrix is calculated
        \param epn              (i) nodal pressure at n
        \param dt               (i) timestep
        */
      void evaluate_collapse(RedAirway* ele, Core::LinAlg::SerialDenseVector& epn,
          Teuchos::ParameterList& params, double dt);
      /*!
        \brief Calculate Pextn and Pextnp

        \param ele              (i) the element those matrix is calculated
        \param epn              (i) nodal pressure at n
        \param dt               (i) timestep
        */
      void compute_pext(RedAirway* ele, const Core::LinAlg::Vector<double>& on,
          const Core::LinAlg::Vector<double>& pnp, Teuchos::ParameterList& params);


      void evaluate_terminal_bc(RedAirway* ele, Teuchos::ParameterList& params,
          Core::FE::Discretization& discretization, std::vector<int>& lm,
          Core::LinAlg::SerialDenseVector& rhs, std::shared_ptr<Core::Mat::Material> mat) override;

      /*!
        \brief get the initial values of the degrees of freedome at the node

        \param ele              (i) the element those matrix is calculated
        \param eqnp             (i) nodal volumetric flow rate at n+1
        \param evelnp           (i) nodal velocity at n+1
        \param eareanp          (i) nodal cross-sectional area at n+1
        \param eprenp           (i) nodal pressure at n+1
        \param estif            (o) element matrix to calculate
        \param eforce           (o) element rhs to calculate
        \param material         (i) airway material/dimesion
        \param time             (i) current simulation time
        \param dt               (i) timestep
        */
      void initial(RedAirway* ele, Teuchos::ParameterList& params,
          Core::FE::Discretization& discretization, std::vector<int>& lm,
          Core::LinAlg::SerialDenseVector& radii_in, Core::LinAlg::SerialDenseVector& radii_out,
          std::shared_ptr<const Core::Mat::Material> material) override;

      /*!
       \Essential functions to compute the results of essential matrices
      */
      void calc_flow_rates(RedAirway* ele, Teuchos::ParameterList& params,
          Core::FE::Discretization& discretization, std::vector<int>& lm,
          std::shared_ptr<Core::Mat::Material> mat) override;

      /*!
       \Essential functions to compute the volume of elements
      */
      void calc_elem_volume(RedAirway* ele, Teuchos::ParameterList& params,
          Core::FE::Discretization& discretization, std::vector<int>& lm,
          std::shared_ptr<Core::Mat::Material> mat) override;

      /*!
       \Essential functions to evaluate the coupled results
      */
      void get_coupled_values(RedAirway* ele, Teuchos::ParameterList& params,
          Core::FE::Discretization& discretization, std::vector<int>& lm,
          std::shared_ptr<Core::Mat::Material> material) override;

      static const Core::Conditions::Condition* get_first_condition(
          const Core::FE::Discretization& discretization, const Core::Nodes::Node* node,
          const std::string& name);
    };

  }  // namespace Elements
}  // namespace Discret

FOUR_C_NAMESPACE_CLOSE

#endif
