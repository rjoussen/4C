// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_ART_NET_ARTERY_ELE_CALC_LIN_EXP_HPP
#define FOUR_C_ART_NET_ARTERY_ELE_CALC_LIN_EXP_HPP

#include "4C_config.hpp"

#include "4C_art_net_artery.hpp"
#include "4C_art_net_artery_ele_action.hpp"
#include "4C_art_net_artery_ele_calc.hpp"
#include "4C_fem_general_utils_local_connectivity_matrices.hpp"
#include "4C_linalg_fixedsizematrix.hpp"

FOUR_C_NAMESPACE_OPEN


namespace Discret
{
  namespace Elements
  {
    /// Internal artery implementation
    /*!
      This internal class keeps all the working arrays needed to
      calculate the artery element. Additionally the method Sysmat()
      provides a clean and fast element implementation.

      <h3>Purpose</h3>

    */

    template <Core::FE::CellType distype>
    class ArteryEleCalcLinExp : public ArteryEleCalc<distype>
    {
     private:
      using my = ArteryEleCalc<distype>;

      /// private constructor, since we are a Singleton.
      ArteryEleCalcLinExp(const int numdofpernode, const std::string& disname);

     public:
      //! Singleton access method
      static ArteryEleCalcLinExp<distype>* instance(
          const int numdofpernode, const std::string& disname);

      int evaluate(Artery* ele, Teuchos::ParameterList& params,
          Core::FE::Discretization& discretization, Core::Elements::LocationArray& la,
          Core::LinAlg::SerialDenseMatrix& elemat1, Core::LinAlg::SerialDenseMatrix& elemat2,
          Core::LinAlg::SerialDenseVector& elevec1, Core::LinAlg::SerialDenseVector& elevec2,
          Core::LinAlg::SerialDenseVector& elevec3,
          std::shared_ptr<Core::Mat::Material> mat) override;

      int scatra_evaluate(Artery* ele, Teuchos::ParameterList& params,
          Core::FE::Discretization& discretization, std::vector<int>& lm,
          Core::LinAlg::SerialDenseMatrix& elemat1, Core::LinAlg::SerialDenseMatrix& elemat2,
          Core::LinAlg::SerialDenseVector& elevec1, Core::LinAlg::SerialDenseVector& elevec2,
          Core::LinAlg::SerialDenseVector& elevec3,
          std::shared_ptr<Core::Mat::Material> mat) override;

      int evaluate_service(Artery* ele, const Arteries::Action action,
          Teuchos::ParameterList& params, Core::FE::Discretization& discretization,
          Core::Elements::LocationArray& la, Core::LinAlg::SerialDenseMatrix& elemat1,
          Core::LinAlg::SerialDenseMatrix& elemat2, Core::LinAlg::SerialDenseVector& elevec1,
          Core::LinAlg::SerialDenseVector& elevec2, Core::LinAlg::SerialDenseVector& elevec3,
          std::shared_ptr<Core::Mat::Material> mat) override;

      /*!
        \brief calculate element matrix and rhs

        \param ele              (i) the element those matrix is calculated
        \param eqnp             (i) nodal volumetric flow rate at n+1
        \param evelnp           (i) nodal velocity at n+1
        \param eareanp          (i) nodal cross-sectional area at n+1
        \param eprenp           (i) nodal pressure at n+1
        \param estif            (o) element matrix to calculate
        \param eforce           (o) element rhs to calculate
        \param material         (i) artery material/dimesion
        \param time             (i) current simulation time
        \param dt               (i) timestep
        */
      void sysmat(Artery* ele, const Core::LinAlg::Matrix<my::iel_, 1>& eqnp,
          const Core::LinAlg::Matrix<my::iel_, 1>& eareanp,
          Core::LinAlg::Matrix<2 * my::iel_, 2 * my::iel_>& sysmat,
          Core::LinAlg::Matrix<2 * my::iel_, 1>& rhs,
          std::shared_ptr<const Core::Mat::Material> material, double dt);

      void scatra_sysmat(Artery* ele, const Core::LinAlg::Matrix<2 * my::iel_, 1>& escatran,
          const Core::LinAlg::Matrix<my::iel_, 1>& ewfnp,
          const Core::LinAlg::Matrix<my::iel_, 1>& ewbnp,
          const Core::LinAlg::Matrix<my::iel_, 1>& eareanp,
          const Core::LinAlg::Matrix<my::iel_, 1>& earean,
          Core::LinAlg::Matrix<2 * my::iel_, 2 * my::iel_>& sysmat,
          Core::LinAlg::Matrix<2 * my::iel_, 1>& rhs, const Core::Mat::Material& material,
          double dt);

      virtual bool solve_riemann(Artery* ele, Teuchos::ParameterList& params,
          Core::FE::Discretization& discretization, std::vector<int>& lm,
          std::shared_ptr<const Core::Mat::Material> mat);

      virtual void evaluate_terminal_bc(Artery* ele, Teuchos::ParameterList& params,
          Core::FE::Discretization& discretization, std::vector<int>& lm,
          std::shared_ptr<Core::Mat::Material> mat);

      virtual void evaluate_scatra_bc(Artery* ele, Teuchos::ParameterList& params,
          Core::FE::Discretization& discretization, std::vector<int>& lm,
          std::shared_ptr<Core::Mat::Material> material);

      virtual void calc_postprocessing_values(Artery* ele, Teuchos::ParameterList& params,
          Core::FE::Discretization& discretization, std::vector<int>& lm,
          std::shared_ptr<Core::Mat::Material> mat);

      virtual void calc_scatra_from_scatra_fw(Artery* ele, Teuchos::ParameterList& params,
          Core::FE::Discretization& discretization, std::vector<int>& lm,
          std::shared_ptr<Core::Mat::Material> material);

      virtual void evaluate_wf_and_wb(Artery* ele, Teuchos::ParameterList& params,
          Core::FE::Discretization& discretization, std::vector<int>& lm,
          std::shared_ptr<Core::Mat::Material> material);

      virtual void solve_scatra_analytically(Artery* ele, Teuchos::ParameterList& params,
          Core::FE::Discretization& discretization, std::vector<int>& lm,
          std::shared_ptr<Core::Mat::Material> material);

      /*!
        \brief get the initial values of the degrees of freedome at the node

        \param ele              (i) the element those matrix is calculated
        \param eqnp             (i) nodal volumetric flow rate at n+1
        \param evelnp           (i) nodal velocity at n+1
        \param eareanp          (i) nodal cross-sectional area at n+1
        \param eprenp           (i) nodal pressure at n+1
        \param estif            (o) element matrix to calculate
        \param eforce           (o) element rhs to calculate
        \param material         (i) artery material/dimesion
        \param time             (i) current simulation time
        \param dt               (i) timestep
        */
      virtual void initial(Artery* ele, Teuchos::ParameterList& params,
          Core::FE::Discretization& discretization, std::vector<int>& lm,
          std::shared_ptr<const Core::Mat::Material> material);

      /*!
       \Essential functions to compute the results of essential matrices
      */
     private:
      //! nodal volumetric flow rate at time step "n"
      Core::LinAlg::Matrix<my::iel_, 1> qn_;
      //! nodal cross-sectional area at time step "n"
      Core::LinAlg::Matrix<my::iel_, 1> an_;
      //! vector containing the initial cross-sectional area at the element nodes
      Core::LinAlg::Matrix<my::iel_, 1> area0_;
      //! vector containing the initial thickness at the element nodes
      Core::LinAlg::Matrix<my::iel_, 1> th_;
      //! vector containing the initial Youngs modulus at the element nodes
      Core::LinAlg::Matrix<my::iel_, 1> young_;
      //! vector containing the fixed external pressure
      Core::LinAlg::Matrix<my::iel_, 1> pext_;
    };

  }  // namespace Elements
}  // namespace Discret

FOUR_C_NAMESPACE_CLOSE

#endif
