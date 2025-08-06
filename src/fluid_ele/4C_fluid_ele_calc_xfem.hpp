// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_FLUID_ELE_CALC_XFEM_HPP
#define FOUR_C_FLUID_ELE_CALC_XFEM_HPP

#include "4C_config.hpp"

#include "4C_fluid_ele_calc.hpp"
#include "4C_fluid_ele_calc_xfem_coupling.hpp"
#include "4C_utils_singleton_owner.hpp"
#include "4C_xfem_condition_manager.hpp"

FOUR_C_NAMESPACE_OPEN

namespace XFEM
{
  class CouplingBase;
}

namespace Discret
{
  namespace Elements
  {
    class FluidEleParameterXFEM;

    /// Fluid element interface coupling implementation with XFEM
    /*!
      This internal class keeps all the working arrays needed to
      calculate the interface stabilization for XFEM with fluid elements.
      The method element_xfem_interface_hybrid_lm() provides a clean and fast element implementation
      for Mixed/Stress/Hybrid interface coupling, using either Cauchy stress-based or viscous
      stress-based Lagrange multipliers.
      The method element_xfem_interface_nit() provides a clean and fast element implementation
      for interface coupling using Nitsche's method (only the cut element sided mortaring
      for xfluid and xfluidfluid applications).
      The method ElementXfemInterfaceNIT2() provides a clean and fast element implementation
      for interface coupling using Nitsche's method also for non-"cut element sided mortaring"
      (can be used for embedded mortaring in case of Xfluidfluid).

      <h3>Purpose</h3>

      The fluid element will allocate exactly one object of this class for all
      fluid elements with the same number of nodes in the mesh. This
      allows us to use exactly matching working arrays (and keep them
      around.)

      The code is meant to be as clean as possible. This is the only way
      to keep it fast. The number of working arrays has to be reduced to
      a minimum so that the element fits into the cache. (There might be
      room for improvements.)

      <h3>Usability</h3>

      The calculations are done by the EvaluateXfemInterface...() methods.

     */
    template <Core::FE::CellType distype>
    class FluidEleCalcXFEM : public FluidEleCalc<distype>
    {
      /// private Constructor since we are a Singleton.
      FluidEleCalcXFEM();

      /// private copy Constructor since we are a Singleton.
      FluidEleCalcXFEM(FluidEleCalcXFEM const& copy);

      /// private assignment operator since we are a Singleton.
      FluidEleCalcXFEM& operator=(FluidEleCalcXFEM const& copy);

      using my = FluidEleCalc<distype>;
      using my::nen_;
      using my::nsd_;
      using my::numdofpernode_;

     private:
      /// pointer to the cast object, fluid parameter list for XFEM
      Discret::Elements::FluidEleParameterXFEM* fldparaxfem_;


      /// number of stress-dof
      static constexpr int numstressdof_ = Core::FE::DisTypeToNumDeriv2<distype>::numderiv2;

      //! @name useful constants for DOF-index numbering
      //@{
      static constexpr unsigned Velx = 0;
      static constexpr unsigned Vely = 1;
      static constexpr unsigned Velz = 2;
      static constexpr unsigned Pres = 3;

      static constexpr unsigned Sigmaxx = 0;
      static constexpr unsigned Sigmaxy = 1;
      static constexpr unsigned Sigmaxz = 2;
      static constexpr unsigned Sigmayx = 1;
      static constexpr unsigned Sigmayy = 3;
      static constexpr unsigned Sigmayz = 4;
      static constexpr unsigned Sigmazx = 2;
      static constexpr unsigned Sigmazy = 4;
      static constexpr unsigned Sigmazz = 5;
      //@}

      /// get stress dof-index
      unsigned stress_index(unsigned xi, unsigned xj)
      {
        return (xi * xj > 0) ? xi + xj + 1 : xi + xj;
      }

     public:
      /// Singleton access method
      static FluidEleCalcXFEM<distype>* instance(
          Core::Utils::SingletonAction action = Core::Utils::SingletonAction::create);

      /// evaluate the XFEM cut element
      int evaluate_xfem(Discret::Elements::Fluid* ele, Core::FE::Discretization& discretization,
          const std::vector<int>& lm, Teuchos::ParameterList& params,
          std::shared_ptr<Core::Mat::Material>& mat, Core::LinAlg::SerialDenseMatrix& elemat1,
          Core::LinAlg::SerialDenseMatrix& elemat2, Core::LinAlg::SerialDenseVector& elevec1,
          Core::LinAlg::SerialDenseVector& elevec2, Core::LinAlg::SerialDenseVector& elevec3,
          const std::vector<Core::FE::GaussIntegration>& intpoints,
          const Cut::plain_volumecell_set& cells, bool offdiag = false) override;

      /// evaluate the shape functions in the XFEM
      int integrate_shape_function_xfem(Discret::Elements::Fluid* ele,
          Core::FE::Discretization& discretization, const std::vector<int>& lm,
          Core::LinAlg::SerialDenseVector& elevec1,
          const std::vector<Core::FE::GaussIntegration>& intpoints,
          const Cut::plain_volumecell_set& cells) override;

      /// error computation
      int compute_error(Discret::Elements::Fluid* ele, Teuchos::ParameterList& params,
          std::shared_ptr<Core::Mat::Material>& mat, Core::FE::Discretization& discretization,
          std::vector<int>& lm, Core::LinAlg::SerialDenseVector& ele_dom_norms) override;

      int compute_error(Discret::Elements::Fluid* ele, Teuchos::ParameterList& params,
          std::shared_ptr<Core::Mat::Material>& mat, Core::FE::Discretization& discretization,
          std::vector<int>& lm, Core::LinAlg::SerialDenseVector& ele_dom_norms,
          const Core::FE::GaussIntegration& intpoints) override;

      int compute_error_interface(Discret::Elements::Fluid* ele,  ///< fluid element
          Core::FE::Discretization& dis,                          ///< background discretization
          const std::vector<int>& lm,                             ///< element local map
          const std::shared_ptr<XFEM::ConditionManager>& cond_manager,  ///< XFEM condition manager
          std::shared_ptr<Core::Mat::Material>& mat,                    ///< material
          Core::LinAlg::SerialDenseVector& ele_interf_norms,  /// squared element interface norms
          const std::map<int, std::vector<Cut::BoundaryCell*>>& bcells,  ///< boundary cells
          const std::map<int, std::vector<Core::FE::GaussIntegration>>&
              bintpoints,                          ///< boundary integration points
          const Cut::plain_volumecell_set& vcSet,  ///< set of plain volume cells
          Teuchos::ParameterList& params           ///< parameter list
          ) override;

      /// add terms from mixed/hybrid Lagrange multiplier coupling approach to element matrix and
      /// rhs
      void element_xfem_interface_hybrid_lm(Discret::Elements::Fluid* ele,  ///< fluid element
          Core::FE::Discretization& dis,  ///< background discretization
          const std::vector<int>& lm,     ///< element local map
          const std::shared_ptr<XFEM::ConditionManager>& cond_manager,   ///< XFEM condition manager
          const std::vector<Core::FE::GaussIntegration>& intpoints,      ///< element gauss points
          const std::map<int, std::vector<Cut::BoundaryCell*>>& bcells,  ///< boundary cells
          const std::map<int, std::vector<Core::FE::GaussIntegration>>&
              bintpoints,  ///< boundary integration points
          const std::map<int, std::vector<int>>&
              patchcouplm,  ///< lm vectors for coupling elements, key= global coupling side-Id
          std::map<int, std::vector<Core::LinAlg::SerialDenseMatrix>>&
              side_coupling,                          ///< side coupling matrices
          Teuchos::ParameterList& params,             ///< parameter list
          std::shared_ptr<Core::Mat::Material>& mat,  ///< material
          Core::LinAlg::SerialDenseMatrix& elemat1,  ///< local system matrix of intersected element
          Core::LinAlg::SerialDenseVector&
              elevec1,                             ///< local element vector of intersected element
          Core::LinAlg::SerialDenseMatrix& Cuiui,  ///< coupling matrix of a side with itself
          const Cut::plain_volumecell_set& vcSet   ///< set of plain volume cells
          ) override;

      /// add Nitsche (NIT) interface condition to element matrix and rhs
      void element_xfem_interface_nit(Discret::Elements::Fluid* ele,  ///< fluid element
          Core::FE::Discretization& dis,                              ///< background discretization
          const std::vector<int>& lm,                                 ///< element local map
          const std::shared_ptr<XFEM::ConditionManager>& cond_manager,   ///< XFEM condition manager
          const std::map<int, std::vector<Cut::BoundaryCell*>>& bcells,  ///< boundary cells
          const std::map<int, std::vector<Core::FE::GaussIntegration>>&
              bintpoints,  ///< boundary integration points
          const std::map<int, std::vector<int>>& patchcouplm,
          Teuchos::ParameterList& params,                    ///< parameter list
          std::shared_ptr<Core::Mat::Material>& mat_master,  ///< material for the coupled side
          std::shared_ptr<Core::Mat::Material>& mat_slave,   ///< material for the coupled side
          Core::LinAlg::SerialDenseMatrix& elemat1,          ///< element matrix
          Core::LinAlg::SerialDenseVector& elevec1,          ///< element vector
          const Cut::plain_volumecell_set& vcSet,            ///< volumecell sets in this element
          std::map<int, std::vector<Core::LinAlg::SerialDenseMatrix>>&
              side_coupling,                       ///< side coupling matrices
          Core::LinAlg::SerialDenseMatrix& Cuiui,  ///< ui-ui coupling matrix
          bool evaluated_cut  ///< the CUT was updated before this evaluation is called
          ) override;

      ///
      void calculate_continuity_xfem(Discret::Elements::Fluid* ele,  ///< fluid element
          Core::FE::Discretization& dis,                             ///< discretization
          const std::vector<int>& lm,                                ///< local map
          Core::LinAlg::SerialDenseVector& elevec1,                  ///< element vector
          const Core::FE::GaussIntegration& intpoints                ///< integration points
          ) override;

      ///
      void calculate_continuity_xfem(Discret::Elements::Fluid* ele,  ///< fluid element
          Core::FE::Discretization& dis,                             ///< discretization
          const std::vector<int>& lm,                                ///< local map
          Core::LinAlg::SerialDenseVector& elevec1                   ///< element vector
          ) override;

     private:
      //! evaluate analytical reference solution
      void analytical_reference(const int calcerr,   ///< which reference solution
          const int calcerrfunctno,                  ///< error function number
          Core::LinAlg::Matrix<nsd_, 1>& u,          ///< exact jump vector (coupled)
          Core::LinAlg::Matrix<nsd_, nsd_>& grad_u,  ///< exact velocity gradient
          double& p,                                 ///< exact pressure
          Core::LinAlg::Matrix<nsd_, 1>& xyzint,     ///< xyz position of gaussian point
          const double& t,                           ///< time
          std::shared_ptr<Core::Mat::Material> mat = nullptr);

      //! get the interface jump vectors for velocity and traction at the Gaussian point
      void get_interface_jump_vectors(
          const XFEM::EleCoupCond& coupcond,  ///< coupling condition for given interface side
          std::shared_ptr<XFEM::CouplingBase> coupling,  ///< coupling object
          Core::LinAlg::Matrix<nsd_, 1>&
              ivelint_jump,  ///< prescribed interface jump vector for velocity
          Core::LinAlg::Matrix<nsd_, 1>&
              itraction_jump,  ///< prescribed interface jump vector for traction
          Core::LinAlg::Matrix<nsd_, nsd_>& proj_tangential,  ///< tangential projection matrix
          Core::LinAlg::Matrix<nsd_, nsd_>&
              LB_proj_matrix,  ///< prescribed projection matrix for laplace-beltrami problems
          const Core::LinAlg::Matrix<nsd_, 1>& x,       ///< global coordinates of Gaussian point
          const Core::LinAlg::Matrix<nsd_, 1>& normal,  ///< normal vector at Gaussian point
          Discret::Elements::XFLUID::SlaveElementInterface<distype>*
              si,                                       ///< side implementation for cutter element
          Core::LinAlg::Matrix<3, 1>& rst,              ///< local coordinates of GP for bg element
          double& kappa_m,                              ///< fluid sided weighting
          double& visc_m,                               ///< fluid sided weighting
          double& visc_s,                               ///< slave sided dynamic viscosity
          Core::LinAlg::Matrix<3, 1>& rst_slave,        ///< local coord of gp in slave element
          std::vector<double>& eledisp,                 ///< slave element displacement vector
          Core::Elements::Element* coupl_ele = nullptr  ///< slave coupling element
      );

      //! get the interface jump vectors for velocity and traction at the Gaussian point for
      //! previous time step
      void get_interface_jump_vectors_old_state(
          const XFEM::EleCoupCond& coupcond,  ///< coupling condition for given interface side
          XFEM::CouplingBase& coupling,       ///< coupling object
          Core::LinAlg::Matrix<nsd_, 1>&
              ivelintn_jump,  ///< prescribed interface jump vector for velocity
          Core::LinAlg::Matrix<nsd_, 1>&
              itractionn_jump,  ///< prescribed interface jump vector for traction
          const Core::LinAlg::Matrix<nsd_, 1>& x,       ///< global coordinates of Gaussian point
          const Core::LinAlg::Matrix<nsd_, 1>& normal,  ///< normal vector at Gaussian point
          Discret::Elements::XFLUID::SlaveElementInterface<distype>&
              si,                          ///< side implementation for cutter element
          const double& presn_m,           ///< coupling master pressure
          Core::LinAlg::Matrix<3, 1>& rst  ///< local coordinates of GP for bg element
      );

      //! build the patch coupling matrix Cuiui containing Cuiui for all cutting sides
      void nit_build_patch_cuiui(
          Core::LinAlg::SerialDenseMatrix&
              Cuiui,  ///< ui-ui patch coupling matrix containing Cuiui for all cutting sides
          std::map<int, std::vector<Core::LinAlg::SerialDenseMatrix>>&
              Cuiui_coupling  ///< Cuiui matrices for all cutting sides
      );

      //! compute viscous part of Nitsche's penalty term scaling for Nitsche's method
      double nit_compute_visc_penalty_stabfac(
          const Core::FE::CellType ele_distype,  ///< the discretization type of the element w.r.t
                                                 ///< which the stabilization factor is computed
          const double inv_hk,                   ///< the inverse characteristic element length
          const double kappa1,                   ///< Weight parameter (parameter +/master side)
          const double kappa2                    ///< Weight parameter (parameter -/slave  side)
      );

      //! get the constant which satisfies the trace inequality depending on the spatial dimension
      //! and polynomial order of the element
      double nit_get_trace_estimate_constant(const Core::FE::CellType ele_distype);

      //! prepare coupling matrices, that include contributions from convective stabilization
      void hybrid_lm_create_special_contribution_matrices(
          XFEM::ConditionManager& cond_manager,  ///< XFEM condition manager
          std::set<int>& begids,                 ///< ids of intersecting boundary elements
          std::map<int, std::vector<Core::LinAlg::SerialDenseMatrix>>&
              side_coupling_extra  ///< contributions to coupling matrices from convective
                                   ///< stabilizations
      );

      //! evaluate Neumann boundary condition
      void evaluate_neumann(const double& timefacfac,    ///< theta*dt
          const Core::LinAlg::Matrix<nen_, 1>& funct_m,  ///< coupling master shape functions
          const Core::LinAlg::Matrix<nsd_, 1>&
              itraction_jump,  ///< prescribed interface traction, jump height for coupled problems
          Core::LinAlg::SerialDenseMatrix::Base& elevec1  ///< element vector
      );

      //! build traction vector w.r.t fluid domain
      void build_traction_vector(Core::LinAlg::Matrix<nsd_, 1>& traction,  ///< traction vector
          double& press,                         ///< pressure at gaussian point
          Core::LinAlg::Matrix<nsd_, 1>& normal  ///< normal vector
      );

      //! assemble side's interface force
      void assemble_interface_force(
          Core::LinAlg::Vector<double>& iforcecol,  ///< interface force column vector
          Core::FE::Discretization& cutdis,         ///< cut discretization
          std::vector<int>& lm,                     ///< local dof map
          Core::LinAlg::SerialDenseVector& iforce   ///< interface force vector
      );

      //! evaluate shape function and derivative at point with local coordinates rst
      void eval_func_and_deriv(Core::LinAlg::Matrix<3, 1>& rst);

      //! build matrices from volume-based terms for Cauchy & viscous stress-based mixed/hybrid
      //! LM-coupling
      void hybrid_lm_build_vol_based(const std::vector<Core::FE::GaussIntegration>& intpoints,
          const Cut::plain_volumecell_set& cells,
          const Core::LinAlg::Matrix<nsd_, nen_>& evelaf,  ///< element velocity
          const Core::LinAlg::Matrix<nen_, 1>& epreaf,     ///< element pressure
          Core::LinAlg::Matrix<nen_, nen_>& bK_ss,         ///< block K_ss matrix
          Core::LinAlg::Matrix<nen_, nen_>& invbK_ss,      ///< inverse of block K_ss matrix
          Core::LinAlg::BlockMatrix<Core::LinAlg::Matrix<nen_, nen_>, numstressdof_,
              numdofpernode_>& K_su,  ///< K_su matrix
          Core::LinAlg::BlockMatrix<Core::LinAlg::Matrix<nen_, 1>, numstressdof_, 1>&
              rhs_s,  ///< rhs_s vector
          Core::LinAlg::BlockMatrix<Core::LinAlg::Matrix<nen_, nen_>, numdofpernode_,
              numstressdof_>& K_us,  ///< K_us matrix
          Core::LinAlg::BlockMatrix<Core::LinAlg::Matrix<nen_, nen_>, nsd_, nsd_>&
              K_uu,  ///< K_uu matrix
          Core::LinAlg::BlockMatrix<Core::LinAlg::Matrix<nen_, 1>, nsd_, 1>&
              rhs_uu,                    ///< rhs_u(u) vector
          const bool is_MHVS = false,    ///< viscous (true) or Cauchy (false) stress-based LM
          const double mhvs_param = 1.0  ///< stabilizing parameter for viscous stress-based LM
      );

      //! evaluate matrices from volume-based terms for Cauchy stress-based mixed/hybrid LM coupling
      //! at current Gauss-point
      void mhcs_evaluate_vol_based(
          const Core::LinAlg::Matrix<nsd_, nen_>& evelaf,  ///< element velocity
          const Core::LinAlg::Matrix<nen_, 1>& epreaf,     ///< element pressure
          Core::LinAlg::Matrix<nen_, nen_>& bK_ss,         ///< block K_ss matrix
          Core::LinAlg::Matrix<nen_, nen_>& invbK_ss,      ///< inverse of block K_ss matrix
          Core::LinAlg::BlockMatrix<Core::LinAlg::Matrix<nen_, nen_>, numstressdof_,
              numdofpernode_>& K_su,  ///< K_su matrix
          Core::LinAlg::BlockMatrix<Core::LinAlg::Matrix<nen_, 1>, numstressdof_, 1>&
              rhs_s  ///< rhs_s vector
      );

      //! evaluate matrices from volume-based terms for viscous stress-based mixed/hybrid LM
      //! coupling at current Gauss-point
      void mhvs_evaluate_vol_based(
          const Core::LinAlg::Matrix<nsd_, nen_>& evelaf,  ///< element velocity
          Core::LinAlg::Matrix<nen_, nen_>& bK_ss,         ///< block K_ss matrix
          Core::LinAlg::Matrix<nen_, nen_>& invbK_ss,      ///< inverse of block K_ss matrix
          Core::LinAlg::BlockMatrix<Core::LinAlg::Matrix<nen_, nen_>, numstressdof_,
              numdofpernode_>& K_su,  ///< K_su matrix
          Core::LinAlg::BlockMatrix<Core::LinAlg::Matrix<nen_, 1>, numstressdof_, 1>&
              rhs_s,  ///< rhs_s vector
          Core::LinAlg::BlockMatrix<Core::LinAlg::Matrix<nen_, nen_>, numdofpernode_,
              numstressdof_>& K_us,  ///< K_us matrix
          Core::LinAlg::BlockMatrix<Core::LinAlg::Matrix<nen_, nen_>, nsd_, nsd_>&
              K_uu,  ///< K_uu matrix
          Core::LinAlg::BlockMatrix<Core::LinAlg::Matrix<nen_, 1>, nsd_, 1>&
              rhs_uu,               ///< rhs_uu vector
          const double& mhvs_param  ///< stabilizing parameter
      );

      //! evaluate matrices from surface-based terms for Cauchy & viscous stress-based mixed/hybrid
      //! LM coupling at current Gauss-point
      void hybrid_lm_evaluate_surf_based(Discret::Elements::XFLUID::HybridLMInterface<distype>& si,
          const Core::LinAlg::Matrix<nen_, nen_>& bK_ss,
          Core::LinAlg::BlockMatrix<Core::LinAlg::Matrix<nen_, nen_>, numstressdof_,
              numdofpernode_>& K_su,
          Core::LinAlg::BlockMatrix<Core::LinAlg::Matrix<nen_, nen_>, numdofpernode_,
              numstressdof_>& K_us,
          Core::LinAlg::BlockMatrix<Core::LinAlg::Matrix<nen_, 1>, numstressdof_, 1>& rhs_s,
          const Core::LinAlg::Matrix<nen_, 1>& epreaf,
          Core::LinAlg::BlockMatrix<Core::LinAlg::Matrix<nen_, nen_>, nsd_, nsd_>& K_uu,
          Core::LinAlg::BlockMatrix<Core::LinAlg::Matrix<nen_, 1>, nsd_, 1>& rhs_uu,
          Core::LinAlg::BlockMatrix<Core::LinAlg::Matrix<nen_, nen_>, nsd_, 1>& G_up,
          Core::LinAlg::BlockMatrix<Core::LinAlg::Matrix<nen_, nen_>, 1, nsd_>& G_pu,
          Core::LinAlg::BlockMatrix<Core::LinAlg::Matrix<nen_, 1>, nsd_, 1>& rhs_up,
          Core::LinAlg::Matrix<nen_, 1>& rhs_pu, const Core::LinAlg::Matrix<nsd_, 1>& normal,
          const double& timesurffac, const Core::LinAlg::Matrix<nsd_, 1>& ivelint_jump,
          const Core::LinAlg::Matrix<nsd_, 1>& itraction_jump, const bool eval_side_coupling,
          const bool is_MHVS);

      //! Initiates dummy variables and calls FLuidEleCalc get_material_params routine.
      void get_material_parameters_volume_cell(std::shared_ptr<const Core::Mat::Material> material,
          double& densaf,  // done
          double& viscaf,  // done
          double& gamma    // done
      );

      // Density on each side of interface
      double densaf_master_;  ///< density at master side
      double densaf_slave_;   ///< density at slave side

      // Viscosity on each side of interface
      double viscaf_master_;  ///< viscosity at master side
      double viscaf_slave_;   ///< viscosity at slave side

      // Surface tension on each side (If it is not the same, error should be thrown.)
      double gamma_m_;  ///< surface tension coefficient at master side
      double gamma_s_;  ///< surface tension coefficient at master side

      Core::LinAlg::Matrix<nsd_, nen_>
          evelaf_;  ///< element velocity at time t^n+af, implemented also in base class, not
                    ///< accessible via extract_values_from_global_vector
      Core::LinAlg::Matrix<nen_, 1>
          epreaf_;  ///< element pressure at time t^n+af, implemented also in base class, not
                    ///< accessible via extract_values_from_global_vector

      Core::LinAlg::Matrix<nsd_, nen_>
          eveln_;  ///< element velocity at time t^n, implemented also in base class, not accessible
                   ///< via extract_values_from_global_vector
      Core::LinAlg::Matrix<nen_, 1>
          epren_;  ///< element velocity at time t^n, implemented also in base
                   ///< class, not accessible via extract_values_from_global_vector

      Core::LinAlg::Matrix<nsd_, 1> ivelint_jump_;        ///< interface velocity jump at t^n+1
      Core::LinAlg::Matrix<nsd_, 1> itraction_jump_;      ///< interface traction jump at t^n+1
      Core::LinAlg::Matrix<nsd_, nsd_> proj_tangential_;  ///< tangential projection matrix at t^n+1
      Core::LinAlg::Matrix<nsd_, nsd_>
          lb_proj_matrix_;  ///< interface matrix jump (for Laplace-Beltrami) at t^n+1

      std::vector<Core::LinAlg::SerialDenseMatrix>
          solid_stress_;  ///< hold information about solid stress ([0]...traction,
                          ///< [1]...dtraction_dv, [2-4]...d2traction_dv2)

      Core::LinAlg::Matrix<nsd_, 1> ivelintn_jump_;    ///< interface velocity jump at t^n
      Core::LinAlg::Matrix<nsd_, 1> itractionn_jump_;  ///< interface traction jump at t^n

      Core::LinAlg::Matrix<nsd_, 1>
          velint_s_;  ///< velocity of slave side at time t^n+1 at integration point
      Core::LinAlg::Matrix<nsd_, 1>
          velintn_s_;  ///< velocity of slave side at time t^n at integration point

      Core::LinAlg::Matrix<nsd_, 1>
          rst_;  ///< local coordinates of Gaussian point w.r.t background element

      Core::LinAlg::Matrix<nsd_, 1> normal_;  ///< normal vector
      Core::LinAlg::Matrix<nsd_, 1> x_side_;  ///< gauss-point coordinates
      Core::LinAlg::Matrix<nsd_, 1>
          x_gp_lin_;  ///< gauss-point in xyz-system on linearized interface
    };


  }  // namespace Elements
}  // namespace Discret

FOUR_C_NAMESPACE_CLOSE

#endif
