// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_CONTACT_ABSTRACT_DATA_CONTAINER_HPP
#define FOUR_C_CONTACT_ABSTRACT_DATA_CONTAINER_HPP

#include "4C_config.hpp"

#include "4C_mortar_strategy_base.hpp"

FOUR_C_NAMESPACE_OPEN

namespace CONTACT
{
  /*! \brief Data container object for the abstract strategy
   *
   *  This object makes it possible to interchange and share the current state of the
   *  contact simulation between different strategy objects. By using this the
   *  actual strategy stays stateless!
   *
   */
  class AbstractStrategyDataContainer : public Mortar::StrategyDataContainer
  {
   public:
    //! constructor
    AbstractStrategyDataContainer();

    //! @name Accessors
    //!@{

    //! Return parallel unbalance factors (evaluation time) for current time step \f$t_{n+1}\f$
    std::vector<double>& unbalance_time_factors() { return unbalance_evaluation_time_; };
    const std::vector<double>& unbalance_time_factors() const
    {
      return unbalance_evaluation_time_;
    };

    //! Return parallel unbalance factors (number of slave elements) for current time step
    //! \f$t_{n+1}\f$
    std::vector<int>& unbalance_element_factors() { return unbalance_num_slave_elements_; };
    const std::vector<int>& unbalance_element_factors() const
    {
      return unbalance_num_slave_elements_;
    };

    //! return global Lagrange mult. dof row map (of all interfaces)
    std::shared_ptr<Core::LinAlg::Map>& global_lm_dof_row_map_ptr() { return glmdofrowmap_; };
    std::shared_ptr<const Core::LinAlg::Map> global_lm_dof_row_map_ptr() const
    {
      return glmdofrowmap_;
    };

    //! return global reference dof row map for self contact Lagr. multipliers (of all interfaces)
    std::shared_ptr<Core::LinAlg::Map>& global_self_contact_ref_dof_row_map_ptr()
    {
      return gscrefdofrowmap_;
    };
    std::shared_ptr<const Core::LinAlg::Map> global_self_contact_ref_dof_row_map_ptr() const
    {
      return gscrefdofrowmap_;
    };

    //! return global self-contact Lagrange mult. dof row map (of all interfaces)
    std::shared_ptr<Core::LinAlg::Map>& global_self_contact_lm_dof_row_map_ptr()
    {
      return gsclmdofrowmap_;
    };
    std::shared_ptr<const Core::LinAlg::Map> global_self_contact_lm_dof_row_map_ptr() const
    {
      return gsclmdofrowmap_;
    };

    //! return global slave node row map (of all interfaces)
    std::shared_ptr<Core::LinAlg::Map>& global_slave_node_row_map_ptr() { return gsnoderowmap_; };
    std::shared_ptr<const Core::LinAlg::Map> global_slave_node_row_map_ptr() const
    {
      return gsnoderowmap_;
    };

    //! return global master node row map (of all interfaces)
    std::shared_ptr<Core::LinAlg::Map>& global_master_node_row_map_ptr() { return gmnoderowmap_; };
    std::shared_ptr<const Core::LinAlg::Map> global_master_node_row_map_ptr() const
    {
      return gmnoderowmap_;
    };

    //! return global slave dof row map (of all interfaces)
    std::shared_ptr<Core::LinAlg::Map>& global_slave_dof_row_map_ptr() { return gsdofrowmap_; };
    std::shared_ptr<const Core::LinAlg::Map> global_slave_dof_row_map_ptr() const
    {
      return gsdofrowmap_;
    };

    //! return global master dof row map (of all interfaces)
    std::shared_ptr<Core::LinAlg::Map>& global_master_dof_row_map_ptr() { return gmdofrowmap_; };
    std::shared_ptr<const Core::LinAlg::Map> global_master_dof_row_map_ptr() const
    {
      return gmdofrowmap_;
    };

    //! return global internal dof row map
    std::shared_ptr<Core::LinAlg::Map>& global_internal_dof_row_map_ptr() { return gndofrowmap_; };
    std::shared_ptr<const Core::LinAlg::Map> global_internal_dof_row_map_ptr() const
    {
      return gndofrowmap_;
    };

    //! return global slave and master dof row map (s+m map)
    std::shared_ptr<Core::LinAlg::Map>& global_slave_master_dof_row_map_ptr()
    {
      return gsmdofrowmap_;
    };
    std::shared_ptr<const Core::LinAlg::Map> global_slave_master_dof_row_map_ptr() const
    {
      return gsmdofrowmap_;
    };

    //! return global displacement dof row map (s+m+n map)
    std::shared_ptr<Core::LinAlg::Map>& global_disp_dof_row_map_ptr() { return gdisprowmap_; };
    std::shared_ptr<const Core::LinAlg::Map> global_disp_dof_row_map_ptr() const
    {
      return gdisprowmap_;
    };

    //! return global active slave node row map (of all interfaces)
    std::shared_ptr<Core::LinAlg::Map>& global_active_node_row_map_ptr() { return gactivenodes_; };
    std::shared_ptr<const Core::LinAlg::Map> global_active_node_row_map_ptr() const
    {
      return gactivenodes_;
    };
    Core::LinAlg::Map& global_active_node_row_map()
    {
      if (!gactivenodes_) FOUR_C_THROW("The gactivenodes_ is not initialized!");
      return *gactivenodes_;
    }

    //! return global active slave dof row map (of all interfaces)
    std::shared_ptr<Core::LinAlg::Map>& global_active_dof_row_map_ptr() { return gactivedofs_; };
    std::shared_ptr<const Core::LinAlg::Map> global_active_dof_row_map_ptr() const
    {
      return gactivedofs_;
    };
    Core::LinAlg::Map& global_active_dof_row_map()
    {
      if (!gactivedofs_) FOUR_C_THROW("The gAugActiveSlaveDofsPtr_ is not initialized!");
      return *gactivedofs_;
    }


    //! return global inactive slave node row map (of all interfaces)
    std::shared_ptr<Core::LinAlg::Map>& global_inactive_node_row_map_ptr()
    {
      return ginactivenodes_;
    };
    std::shared_ptr<const Core::LinAlg::Map> global_inactive_node_row_map_ptr() const
    {
      return ginactivenodes_;
    };
    Core::LinAlg::Map& global_inactive_node_row_map()
    {
      if (!ginactivenodes_) FOUR_C_THROW("The ginactivenodes_ is not initialized!");
      return *ginactivenodes_;
    }

    //! return global inactive slave dof row map (of all interfaces)
    std::shared_ptr<Core::LinAlg::Map>& global_inactive_dof_row_map_ptr()
    {
      return ginactivedofs_;
    };
    std::shared_ptr<const Core::LinAlg::Map> global_inactive_dof_row_map_ptr() const
    {
      return ginactivedofs_;
    };
    Core::LinAlg::Map& global_inactive_dof_row_map()
    {
      if (!ginactivedofs_) FOUR_C_THROW("The ginactivedofs_ is not initialized!");
      return *ginactivedofs_;
    }


    //! return global active slave dof row map in normal direction (of all interfaces)
    std::shared_ptr<Core::LinAlg::Map>& global_active_n_dof_row_map_ptr() { return gactiven_; };
    std::shared_ptr<const Core::LinAlg::Map> global_active_n_dof_row_map_ptr() const
    {
      return gactiven_;
    };
    Core::LinAlg::Map& global_active_n_dof_row_map()
    {
      if (!gactiven_) FOUR_C_THROW("The gactiven_ is not initialized!");
      return *gactiven_;
    }

    //! return global active slave dof row map in tangential direction (of all interfaces)
    std::shared_ptr<Core::LinAlg::Map>& global_active_t_dof_row_map_ptr() { return gactivet_; };
    std::shared_ptr<const Core::LinAlg::Map> global_active_t_dof_row_map_ptr() const
    {
      return gactivet_;
    };
    Core::LinAlg::Map& global_active_t_dof_row_map()
    {
      if (!gactivet_) FOUR_C_THROW("The gactivet_ is not initialized!");
      return *gactivet_;
    }

    //! return global slip slave node row map (of all interfaces)
    std::shared_ptr<Core::LinAlg::Map>& global_slip_node_row_map_ptr() { return gslipnodes_; };
    std::shared_ptr<const Core::LinAlg::Map> global_slip_node_row_map_ptr() const
    {
      return gslipnodes_;
    };

    //! return global slip slave dof row map (of all interfaces)
    std::shared_ptr<Core::LinAlg::Map>& global_slip_dof_row_map_ptr() { return gslipdofs_; };
    std::shared_ptr<const Core::LinAlg::Map> global_slip_dof_row_map_ptr() const
    {
      return gslipdofs_;
    };

    //! return global slip slave dof row map in tangential direction (of all interfaces)
    std::shared_ptr<Core::LinAlg::Map>& global_slip_t_dof_row_map_ptr() { return gslipt_; };
    std::shared_ptr<const Core::LinAlg::Map> global_slip_t_dof_row_map_ptr() const
    {
      return gslipt_;
    };

    //! return global slave dof row map associated with vertex nodes
    std::shared_ptr<Core::LinAlg::Map>& global_slave_dof_vertex_row_map_ptr()
    {
      return gsdof_vertex_;
    };
    std::shared_ptr<const Core::LinAlg::Map> global_slave_dof_vertex_row_map_ptr() const
    {
      return gsdof_vertex_;
    };

    //! return global slave dof row map associated with edge nodes
    std::shared_ptr<Core::LinAlg::Map>& global_slave_dof_edge_row_map_ptr() { return gsdof_edge_; };
    std::shared_ptr<const Core::LinAlg::Map> global_slave_dof_edge_row_map_ptr() const
    {
      return gsdof_edge_;
    };

    //! return global slave dof row map associated with surface nodes
    std::shared_ptr<Core::LinAlg::Map>& global_slave_dof_surface_row_map_ptr()
    {
      return gsdof_surf_;
    };
    std::shared_ptr<const Core::LinAlg::Map> global_slave_dof_surface_row_map_ptr() const
    {
      return gsdof_surf_;
    };

    //! return global LM dof row map (before parallel redistribution)
    std::shared_ptr<Core::LinAlg::Map>& non_redist_global_lm_dof_row_map_ptr()
    {
      return non_redist_glmdofrowmap_;
    };
    std::shared_ptr<const Core::LinAlg::Map> non_redist_global_lm_dof_row_map_ptr() const
    {
      return non_redist_glmdofrowmap_;
    };

    //! return global slave dof row map (before parallel redistribution)
    std::shared_ptr<Core::LinAlg::Map>& non_redist_global_slave_dof_row_map_ptr()
    {
      return non_redist_gsdofrowmap_;
    };
    std::shared_ptr<const Core::LinAlg::Map> non_redist_global_slave_dof_row_map_ptr() const
    {
      return non_redist_gsdofrowmap_;
    };

    //! return global master dof row map (before parallel redistribution)
    std::shared_ptr<Core::LinAlg::Map>& non_redist_global_master_dof_row_map_ptr()
    {
      return non_redist_gmdofrowmap_;
    };
    std::shared_ptr<const Core::LinAlg::Map> non_redist_global_master_dof_row_map_ptr() const
    {
      return non_redist_gmdofrowmap_;
    };

    //! return global slave and master dof row map (before parallel redistribution)
    std::shared_ptr<Core::LinAlg::Map>& non_redist_global_slave_master_dof_row_map_ptr()
    {
      return non_redist_gsmdofrowmap_;
    };
    std::shared_ptr<const Core::LinAlg::Map> non_redist_global_slave_master_dof_row_map_ptr() const
    {
      return non_redist_gsmdofrowmap_;
    };

    //! return global dirichlet toggle of all slave dofs (before parallel redistribution)
    std::shared_ptr<Core::LinAlg::Vector<double>>&
    non_redist_global_slave_dirich_toggle_dof_row_map_ptr()
    {
      return non_redist_gsdirichtoggle_;
    };
    std::shared_ptr<const Core::LinAlg::Vector<double>>
    non_redist_global_slave_dirich_toggle_dof_row_map_ptr() const
    {
      return non_redist_gsdirichtoggle_;
    };

    //! return initial col ele map for binning strategy (s m)
    std::vector<std::shared_ptr<Core::LinAlg::Map>>& initial_sl_ma_ele_col_map()
    {
      return initial_elecolmap_;
    };
    const std::vector<std::shared_ptr<Core::LinAlg::Map>>& initial_sl_ma_ele_col_map() const
    {
      return initial_elecolmap_;
    };

    //! return global Mortar matrix D
    std::shared_ptr<Core::LinAlg::SparseMatrix>& d_matrix_ptr() { return dmatrix_; };
    std::shared_ptr<const Core::LinAlg::SparseMatrix> d_matrix_ptr() const { return dmatrix_; };
    Core::LinAlg::SparseMatrix& d_matrix()
    {
      if (!dmatrix_) FOUR_C_THROW("The dmatrix_ is not initialized!");
      return *dmatrix_;
    }

    //! return global Mortar matrix M
    std::shared_ptr<Core::LinAlg::SparseMatrix>& m_matrix_ptr() { return mmatrix_; };
    std::shared_ptr<const Core::LinAlg::SparseMatrix> m_matrix_ptr() const { return mmatrix_; };
    Core::LinAlg::SparseMatrix& m_matrix()
    {
      if (!mmatrix_) FOUR_C_THROW("The mmatrix_ is not initialized!");
      return *mmatrix_;
    }

    //! return global weighted gap vector g
    std::shared_ptr<Core::LinAlg::Vector<double>>& w_gap_ptr() { return wgap_; };
    std::shared_ptr<const Core::LinAlg::Vector<double>> w_gap_ptr() const { return wgap_; };
    Core::LinAlg::Vector<double>& w_gap()
    {
      if (!wgap_) FOUR_C_THROW("The wgap_ is not initialized!");
      return *wgap_;
    }

    //! return global tangential rhs vector
    std::shared_ptr<Core::LinAlg::Vector<double>>& tang_rhs_ptr() { return tangrhs_; };
    std::shared_ptr<const Core::LinAlg::Vector<double>> tang_rhs_ptr() const { return tangrhs_; };

    //! return global inactive rhs vector
    std::shared_ptr<Core::LinAlg::Vector<double>>& inactive_rhs_ptr() { return inactiverhs_; };
    std::shared_ptr<const Core::LinAlg::Vector<double>> inactive_rhs_ptr() const
    {
      return inactiverhs_;
    };
    Core::LinAlg::Vector<double>& inactive_rhs()
    {
      if (!inactiverhs_) FOUR_C_THROW("The inactiverhs_ is not initialized!");
      return *inactiverhs_;
    }

    //! Return the structural contact right-hand-side contributions of the current time step
    //! \f$t_{n+1}\f$
    std::shared_ptr<Core::LinAlg::Vector<double>>& str_contact_rhs_ptr()
    {
      return str_contact_rhs_ptr_;
    }
    std::shared_ptr<const Core::LinAlg::Vector<double>> str_contact_rhs_ptr() const
    {
      return str_contact_rhs_ptr_;
    }
    Core::LinAlg::Vector<double>& str_contact_rhs()
    {
      if (!str_contact_rhs_ptr_) FOUR_C_THROW("The strContactRhsPtr_ is not initialized!");
      return *str_contact_rhs_ptr_;
    }

    //! return global constraint rhs vector (only for saddlepoint problems)
    std::shared_ptr<Core::LinAlg::Vector<double>>& constr_rhs_ptr() { return constrrhs_; };
    std::shared_ptr<const Core::LinAlg::Vector<double>> constr_rhs_ptr() const
    {
      return constrrhs_;
    };
    Core::LinAlg::Vector<double>& constr_rhs()
    {
      if (!constrrhs_) FOUR_C_THROW("The constrrhs_ is not initialized!");
      return *constrrhs_;
    }

    //! return global Matrix LinD containing slave fc derivatives
    std::shared_ptr<Core::LinAlg::SparseMatrix>& d_lin_matrix_ptr() { return lindmatrix_; };
    std::shared_ptr<const Core::LinAlg::SparseMatrix> d_lin_matrix_ptr() const
    {
      return lindmatrix_;
    };
    Core::LinAlg::SparseMatrix& d_lin_matrix()
    {
      if (!lindmatrix_) FOUR_C_THROW("The augDnLinMatrixPtr_ is not initialized!");
      return *lindmatrix_;
    }

    //! return global Matrix LinM containing master fc derivatives
    std::shared_ptr<Core::LinAlg::SparseMatrix>& m_lin_matrix_ptr() { return linmmatrix_; };
    std::shared_ptr<const Core::LinAlg::SparseMatrix> m_lin_matrix_ptr() const
    {
      return linmmatrix_;
    };
    Core::LinAlg::SparseMatrix& m_lin_matrix()
    {
      if (!linmmatrix_) FOUR_C_THROW("The augMnLinMatrixPtr_ is not initialized!");
      return *linmmatrix_;
    }

    //! return global Matrix kteffnew containing modified jacobian
    std::shared_ptr<Core::LinAlg::SparseMatrix>& kteffnew_matrix_ptr() { return kteffnew_; };
    std::shared_ptr<const Core::LinAlg::SparseMatrix> kteffnew_matrix_ptr() const
    {
      return kteffnew_;
    };
    Core::LinAlg::SparseMatrix& kteffnew_matrix()
    {
      if (!kteffnew_) FOUR_C_THROW("The kteffnewMatrixPtr is not initialized!");
      return *kteffnew_;
    }

    //! return global Mortar matrix D (last end-point \f$t_{n}\f$)
    std::shared_ptr<Core::LinAlg::SparseMatrix>& old_d_matrix_ptr() { return dold_; };
    std::shared_ptr<const Core::LinAlg::SparseMatrix> old_d_matrix_ptr() const { return dold_; };

    //! return global Mortar matrix M (last end-point \f$t_{n}\f$)
    std::shared_ptr<Core::LinAlg::SparseMatrix>& old_m_matrix_ptr() { return mold_; };
    std::shared_ptr<const Core::LinAlg::SparseMatrix> old_m_matrix_ptr() const { return mold_; };

    //! return current vector of Lagrange multipliers at \f$t_{n+1}\f$
    std::shared_ptr<Core::LinAlg::Vector<double>>& lm_ptr() { return z_; };
    std::shared_ptr<const Core::LinAlg::Vector<double>> lm_ptr() const { return z_; };

    //! return old vector of Lagrange multipliers at \f$t_{n}\f$
    std::shared_ptr<Core::LinAlg::Vector<double>>& old_lm_ptr() { return zold_; };
    std::shared_ptr<const Core::LinAlg::Vector<double>> old_lm_ptr() const { return zold_; };

    /*! \brief Return Lagrange multiplier vector increment
     *
     *  \remark This is NOT the increment of z_ between \f$t_{n+1}\f$ and \f$t_{n}\f$!) */
    std::shared_ptr<Core::LinAlg::Vector<double>>& lm_incr_ptr() { return zincr_; };
    std::shared_ptr<const Core::LinAlg::Vector<double>> lm_incr_ptr() const { return zincr_; };

    //! return vector of Lagrange multipliers from last Uzawa step
    std::shared_ptr<Core::LinAlg::Vector<double>>& lm_uzawa_ptr() { return zuzawa_; };
    std::shared_ptr<const Core::LinAlg::Vector<double>> lm_uzawa_ptr() const { return zuzawa_; };

    //! return vector of normal contact forces at \f$t_{n+1}\f$
    std::shared_ptr<Core::LinAlg::Vector<double>>& stress_normal_ptr() { return stressnormal_; };
    std::shared_ptr<const Core::LinAlg::Vector<double>> stress_normal_ptr() const
    {
      return stressnormal_;
    };

    //! return vector of tangential contact forces at \f$t_{n+1}\f$
    std::shared_ptr<Core::LinAlg::Vector<double>>& stress_tangential_ptr()
    {
      return stresstangential_;
    };
    std::shared_ptr<const Core::LinAlg::Vector<double>> stress_tangential_ptr() const
    {
      return stresstangential_;
    };

    //! return vector of normal contact forces at \f$t_{n+1}\f$
    std::shared_ptr<Core::LinAlg::Vector<double>>& force_normal_ptr() { return forcenormal_; };
    std::shared_ptr<const Core::LinAlg::Vector<double>> force_normal_ptr() const
    {
      return forcenormal_;
    };

    //! return vector of tangential contact forces at \f$t_{n+1}\f$
    std::shared_ptr<Core::LinAlg::Vector<double>>& force_tangential_ptr()
    {
      return forcetangential_;
    };
    std::shared_ptr<const Core::LinAlg::Vector<double>> force_tangential_ptr() const
    {
      return forcetangential_;
    };

    //! return time step index at \f$t_{n+1}\f$
    int& step_np() { return stepnp_; };
    int step_np() const { return stepnp_; };

    //! return non-linear (Newton) iteration index
    int& nln_iter() { return iter_; };
    int nln_iter() const { return iter_; };

    //! return flag indicating global contact status
    bool& is_in_contact() { return isincontact_; };
    bool is_in_contact() const { return isincontact_; };

    //! return flag indicating global contact status of this time step (history)
    bool& was_in_contact() { return wasincontact_; };
    bool was_in_contact() const { return wasincontact_; };

    //! return flag indicating global contact status of last time step
    bool& was_in_contact_last_time_step() { return wasincontactlts_; };
    bool was_in_contact_last_time_step() const { return wasincontactlts_; };

    //! return flag indicating potential self contact
    bool& is_self_contact() { return isselfcontact_; };
    bool is_self_contact() const { return isselfcontact_; };

    //! return flag for frictional contact
    bool& is_friction() { return friction_; };
    bool is_friction() const { return friction_; };

    //! return flag for nonsmooth contact
    bool& is_non_smooth_contact() { return non_smooth_contact_; };
    bool is_non_smooth_contact() const { return non_smooth_contact_; };

    //! return flag indicating whether trafo should be applied
    bool& is_dual_quad_slave_trafo() { return dualquadslavetrafo_; };
    bool is_dual_quad_slave_trafo() const { return dualquadslavetrafo_; };

    //! return transformation matrix T for dual quad 3D case
    std::shared_ptr<Core::LinAlg::SparseMatrix>& trafo_ptr() { return trafo_; };
    std::shared_ptr<const Core::LinAlg::SparseMatrix> trafo_ptr() const { return trafo_; };

    //! return inverse trafo matrix T^(-1) for dual quad 3D case
    std::shared_ptr<Core::LinAlg::SparseMatrix>& inv_trafo_ptr() { return invtrafo_; };
    std::shared_ptr<const Core::LinAlg::SparseMatrix> inv_trafo_ptr() const { return invtrafo_; };

    //! return modified global Mortar matrix D
    std::shared_ptr<Core::LinAlg::SparseMatrix>& modified_d_matrix_ptr() { return dmatrixmod_; };
    std::shared_ptr<const Core::LinAlg::SparseMatrix> modified_d_matrix_ptr() const
    {
      return dmatrixmod_;
    };

    //! return modified global Mortar matrix Dold
    std::shared_ptr<Core::LinAlg::SparseMatrix>& old_modified_d_matrix_ptr() { return doldmod_; };
    std::shared_ptr<const Core::LinAlg::SparseMatrix> old_modified_d_matrix_ptr() const
    {
      return doldmod_;
    };

    //! return integration time
    double& int_time() { return inttime_; };
    double int_time() const { return inttime_; };

    //! return mean interface velocity
    std::vector<double>& mean_interface_vels() { return ivel_; };
    const std::vector<double>& mean_interface_vels() const { return ivel_; };

    //! return current used solving strategy
    CONTACT::SolvingStrategy& sol_type() { return stype_; };
    CONTACT::SolvingStrategy sol_type() const { return stype_; };

    //! return direction in which the contact constraints are formulated
    CONTACT::ConstraintDirection& constr_direction() { return constr_direction_; };
    CONTACT::ConstraintDirection constr_direction() const { return constr_direction_; };

    Inpar::Mortar::ParallelRedist& par_type() { return partype_; };
    Inpar::Mortar::ParallelRedist par_type() const { return partype_; };

    //!@}

   private:
    //! global Lagrange multiplier dof row map (of all interfaces)
    std::shared_ptr<Core::LinAlg::Map> glmdofrowmap_;

    //! global reference dof row map for self contact Lagrange multipliers (of all interfaces)
    std::shared_ptr<Core::LinAlg::Map> gscrefdofrowmap_;

    //! global Lagrange mult. dof row map for self contact (of all interfaces)
    std::shared_ptr<Core::LinAlg::Map> gsclmdofrowmap_;

    //! global slave node row map (of all interfaces)
    std::shared_ptr<Core::LinAlg::Map> gsnoderowmap_;

    //! global master node row map (of all interfaces)
    std::shared_ptr<Core::LinAlg::Map> gmnoderowmap_;

    //! global slave dof row map (of all interfaces)
    std::shared_ptr<Core::LinAlg::Map> gsdofrowmap_;

    //! global master dof row map (of all interfaces)
    std::shared_ptr<Core::LinAlg::Map> gmdofrowmap_;

    //! global internal dof row map
    std::shared_ptr<Core::LinAlg::Map> gndofrowmap_;

    //! global slave and master dof row map (s+m map)
    std::shared_ptr<Core::LinAlg::Map> gsmdofrowmap_;

    //! global displacement dof row map (s+m+n map)
    std::shared_ptr<Core::LinAlg::Map> gdisprowmap_;

    //! @name Active set
    //!@{

    //! global active slave node row map (of all interfaces)
    std::shared_ptr<Core::LinAlg::Map> gactivenodes_;

    //! global active slave dof row map (of all interfaces)
    std::shared_ptr<Core::LinAlg::Map> gactivedofs_;

    //! global inactive slave node row map (of all interfaces)
    std::shared_ptr<Core::LinAlg::Map> ginactivenodes_;

    //! global inactive slave dof row map (of all interfaces)
    std::shared_ptr<Core::LinAlg::Map> ginactivedofs_;

    //! global active slave dof row map in normal direction (of all interfaces)
    std::shared_ptr<Core::LinAlg::Map> gactiven_;

    //! global dof row map of matrix T (of all interfaces)
    std::shared_ptr<Core::LinAlg::Map> gactivet_;

    //! global slip slave node row map (of all interfaces)
    std::shared_ptr<Core::LinAlg::Map> gslipnodes_;

    //! global slip slave dof row map (of all interfaces)
    std::shared_ptr<Core::LinAlg::Map> gslipdofs_;

    //! global slip slave dof row map in tangential direction (of all interfaces)
    std::shared_ptr<Core::LinAlg::Map> gslipt_;

    //!@}

    //! global slave dof row map of vertex nodes
    std::shared_ptr<Core::LinAlg::Map> gsdof_vertex_;

    //! global slave dof row map of edge nodes
    std::shared_ptr<Core::LinAlg::Map> gsdof_edge_;

    //! global slave dof row map of surface nodes
    std::shared_ptr<Core::LinAlg::Map> gsdof_surf_;

    //! @name Parallel redistribution
    //!@{

    /*! Max-to-min ratio of evaluation time across all processes for current time step \f$t_{n+1}\f$
     */
    std::vector<double> unbalance_evaluation_time_;

    /*! Max-to-min ratio of number of row slave elements across all processes for current time step
     * \f$t_{n+1}\f$
     */
    std::vector<int> unbalance_num_slave_elements_;

    //! global LM dof row map (before parallel redistribution)
    std::shared_ptr<Core::LinAlg::Map> non_redist_glmdofrowmap_;

    //! global slave dof row map (before parallel redistribution)
    std::shared_ptr<Core::LinAlg::Map> non_redist_gsdofrowmap_;

    //! global master dof row map (before parallel redistribution)
    std::shared_ptr<Core::LinAlg::Map> non_redist_gmdofrowmap_;

    //! global slave and master dof row map (before parallel redistribution)
    std::shared_ptr<Core::LinAlg::Map> non_redist_gsmdofrowmap_;

    //! global dirichlet toggle of all slave dofs (before parallel redistribution)
    std::shared_ptr<Core::LinAlg::Vector<double>> non_redist_gsdirichtoggle_;

    //! parallel redistribution type
    Inpar::Mortar::ParallelRedist partype_;

    //!@}

    //! @name Binning strategy
    //!@{

    //! initial col ele map for binning strategy (s m)
    std::vector<std::shared_ptr<Core::LinAlg::Map>> initial_elecolmap_;

    //!@}

    //! global Mortar matrix \f$D\f$
    std::shared_ptr<Core::LinAlg::SparseMatrix> dmatrix_;

    //! global Mortar matrix \f$M\f$
    std::shared_ptr<Core::LinAlg::SparseMatrix> mmatrix_;

    //! global weighted gap vector \f$g\f$
    std::shared_ptr<Core::LinAlg::Vector<double>> wgap_;

    //! global tangential right-hand side vector (formulation with incremental #z_)
    std::shared_ptr<Core::LinAlg::Vector<double>> tangrhs_;

    /*! \brief Global inactive right-hand side vector
     *
     * This is used for the formulation with incremental #z_ and saddle point system.
     */
    std::shared_ptr<Core::LinAlg::Vector<double>> inactiverhs_;

    //! structural contact right-hand-side vector at \f$t_{n+1}\f$
    std::shared_ptr<Core::LinAlg::Vector<double>> str_contact_rhs_ptr_;

    //! global constraint right-hand side vector (only for saddlepoint problems)
    std::shared_ptr<Core::LinAlg::Vector<double>> constrrhs_;

    //! global Matrix LinD containing slave fc derivatives
    std::shared_ptr<Core::LinAlg::SparseMatrix> lindmatrix_;

    //! global Matrix LinM containing master fc derivatives
    std::shared_ptr<Core::LinAlg::SparseMatrix> linmmatrix_;

    //! global K matrix
    std::shared_ptr<Core::LinAlg::SparseMatrix> kteffnew_;

    //! global Mortar matrix D (last end-point \f$t_{n}\f$)
    std::shared_ptr<Core::LinAlg::SparseMatrix> dold_;

    //! global Mortar matrix M (last end-point \f$t_{n}\f$)
    std::shared_ptr<Core::LinAlg::SparseMatrix> mold_;

    //! current vector of Lagrange multipliers at \f$t_{n+1}\f$
    std::shared_ptr<Core::LinAlg::Vector<double>> z_;

    //! old vector of Lagrange multipliers at \f$t_{n}\f$
    std::shared_ptr<Core::LinAlg::Vector<double>> zold_;

    /*! \brief Lagrange multiplier vector increment within SaddlePointSolve
     *
     *  \remark This is \em not the increment of #z_ between \f$t_{n+1}\f$ and \f$t_{n}\f$!)
     */
    std::shared_ptr<Core::LinAlg::Vector<double>> zincr_;

    //! vector of Lagrange multipliers from last Uzawa step
    std::shared_ptr<Core::LinAlg::Vector<double>> zuzawa_;

    //! vector of normal contact forces at \f$t_{n+1}\f$
    std::shared_ptr<Core::LinAlg::Vector<double>> stressnormal_;

    //! vector of tangential contact forces at \f$t_{n+1}\f$
    std::shared_ptr<Core::LinAlg::Vector<double>> stresstangential_;

    //! vector of normal contact forces at \f$t_{n+1}\f$
    std::shared_ptr<Core::LinAlg::Vector<double>> forcenormal_;

    //! vector of tangential contact forces at \f$t_{n+1}\f$
    std::shared_ptr<Core::LinAlg::Vector<double>> forcetangential_;

    //! @name Counters and indices
    //!@{

    //! time step index at \f$t_{n+1}\f$
    int stepnp_;

    //! Nonlinear iteration index, e.g. Newton iteration
    int iter_;

    //!@}

    //! @name Status flags
    //!@{

    //! flag indicating global contact status
    bool isincontact_;

    //! flag indicating global contact status of this time step (history)
    bool wasincontact_;

    //! flag indicating global contact status of last time step
    bool wasincontactlts_;

    //! flag indicating potential self contact
    bool isselfcontact_;

    //! flag for frictional contact
    bool friction_;

    //! flag for non-smooth contact
    bool non_smooth_contact_;

    //! flag indicating whether trafo should be applied
    bool dualquadslavetrafo_;

    //!@}

    //! transformation matrix T for dual quad 3D case
    std::shared_ptr<Core::LinAlg::SparseMatrix> trafo_;

    //! inverse trafo matrix T^(-1) for dual quad 3D case
    std::shared_ptr<Core::LinAlg::SparseMatrix> invtrafo_;

    //! modified global Mortar matrix D
    std::shared_ptr<Core::LinAlg::SparseMatrix> dmatrixmod_;

    //! modified global Mortar matrix Dold
    std::shared_ptr<Core::LinAlg::SparseMatrix> doldmod_;

    /*! \brief Integration time
     *
     * \todo Is this the wall clock time required to perform the mortar integration?
     */
    double inttime_;

    //! mean interface velocity
    std::vector<double> ivel_;

    //! current used solving strategy
    CONTACT::SolvingStrategy stype_;

    //! direction in which the contact constraints are formulated
    CONTACT::ConstraintDirection constr_direction_;

  };  // class AbstractStrategyDataContainer


}  // namespace CONTACT

FOUR_C_NAMESPACE_CLOSE

#endif