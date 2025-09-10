// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_BEAMINTERACTION_POTENTIAL_PAIR_BEAM_TO_BEAM_HPP
#define FOUR_C_BEAMINTERACTION_POTENTIAL_PAIR_BEAM_TO_BEAM_HPP

#include "4C_config.hpp"

#include "4C_beaminteraction_potential_pair_base.hpp"
#include "4C_linalg_fixedsizematrix.hpp"

#include <Sacado.hpp>

FOUR_C_NAMESPACE_OPEN


// forward declarations
namespace Core::LinAlg
{
  class SerialDenseVector;
  class SerialDenseMatrix;
}  // namespace Core::LinAlg

namespace Discret
{
  namespace Elements
  {
    class Beam3Base;
  }
}  // namespace Discret

namespace BeamInteraction
{
  /*!
   \brief class for potential-based interaction between two 3D beam elements
   */
  template <unsigned int numnodes, unsigned int numnodalvalues, typename T>
  class BeamToBeamPotentialPair : public BeamPotentialPair
  {
   public:
    //! @name Friends
    // no friend classes defined
    //@}

    //! @name Constructors and destructors and related methods
    /*!
    \brief Standard Constructor
    */
    BeamToBeamPotentialPair();



    //! Setup
    void setup() override;
    //@}

    //! @name Derived methods from base class
    /*!
    \brief Evaluate this contact element pair, return value indicates whether pair is active,
           i.e. non-zero values for force and stiffmat are returned
    */
    bool evaluate(Core::LinAlg::SerialDenseVector* forcevec1,
        Core::LinAlg::SerialDenseVector* forcevec2, Core::LinAlg::SerialDenseMatrix* stiffmat11,
        Core::LinAlg::SerialDenseMatrix* stiffmat12, Core::LinAlg::SerialDenseMatrix* stiffmat21,
        Core::LinAlg::SerialDenseMatrix* stiffmat22,
        const std::vector<const Core::Conditions::Condition*> linechargeconds, const double k,
        const double m) override;

    /*
    \brief Update state of translational nodal DoFs (absolute positions and tangents) of both
    elements
    */
    void reset_state(double time, const std::vector<double>& centerline_dofvec_ele1,
        const std::vector<double>& centerline_dofvec_ele2) override;

    /*!
    \brief Get coordinates of all interacting points on element1 and element2
    */
    void get_all_interacting_point_coords_element1(
        std::vector<Core::LinAlg::Matrix<3, 1, double>>& coords) const override
    {
      coords = centerline_coords_gp_1_;
    }

    void get_all_interacting_point_coords_element2(
        std::vector<Core::LinAlg::Matrix<3, 1, double>>& coords) const override
    {
      coords = centerline_coords_gp_2_;
    }

    /*!
    \brief Get forces at all interacting points on element1 and element2
    */
    void get_forces_at_all_interacting_points_element1(
        std::vector<Core::LinAlg::Matrix<3, 1, double>>& forces) const override
    {
      forces = forces_pot_gp_1_;
    }

    void get_forces_at_all_interacting_points_element2(
        std::vector<Core::LinAlg::Matrix<3, 1, double>>& forces) const override
    {
      forces = forces_pot_gp_2_;
    }

    /*!
    \brief Get moments at all interacting points on element1 and element2
    */
    void get_moments_at_all_interacting_points_element1(
        std::vector<Core::LinAlg::Matrix<3, 1, double>>& moments) const override
    {
      moments = moments_pot_gp_1_;
    }

    void get_moments_at_all_interacting_points_element2(
        std::vector<Core::LinAlg::Matrix<3, 1, double>>& moments) const override
    {
      moments = moments_pot_gp_2_;
    }

    /*!
    \brief Get interaction free energy / potential
    */
    double get_energy() const override { return interaction_potential_; }

    /** \brief print this beam potential-based element pair to screen
     *
     *  */
    void print(std::ostream& out) const override;

    /** \brief print this beam potential element pair to screen
     *
     *  */
    void print_summary_one_line_per_active_segment_pair(std::ostream& out) const override;
    //@}

    //! @name Access methods
    /*!
    \brief Get ptr to first beam element
    */
    inline const Discret::Elements::Beam3Base* beam_element1() const { return beam_element1_; };

    /*!
    \brief Get ptr to second beam element
    */
    inline const Discret::Elements::Beam3Base* beam_element2() const { return beam_element2_; };
    //@}

   private:
    //! @name Private evaluation methods

    /** \brief Evaluate forces and stiffness contribution resulting from potential-based interaction
     *         using double length specific potential with approximation for large separations
     *
     */
    void evaluate_fpotand_stiffpot_double_length_specific_large_sep_approx(
        Core::LinAlg::Matrix<3 * numnodes * numnodalvalues, 1, T>& force_pot1,
        Core::LinAlg::Matrix<3 * numnodes * numnodalvalues, 1, T>& force_pot2,
        Core::LinAlg::SerialDenseMatrix* stiffmat11, Core::LinAlg::SerialDenseMatrix* stiffmat12,
        Core::LinAlg::SerialDenseMatrix* stiffmat21, Core::LinAlg::SerialDenseMatrix* stiffmat22);

    /** \brief compute contributions to analytic linearization (i.e. stiffness matrices) at current
     *         Gauss point pair
     *
     */
    void evaluate_stiffpot_analytic_contributions_large_sep_approx(
        Core::LinAlg::Matrix<3, 1, double> const& dist, double const& norm_dist,
        double const& norm_dist_exp1, double q1q2_JacFac_GaussWeights,
        Core::LinAlg::Matrix<1, numnodes * numnodalvalues, double> const& N1_i_GP1,
        Core::LinAlg::Matrix<1, numnodes * numnodalvalues, double> const& N2_i_GP2,
        Core::LinAlg::SerialDenseMatrix& stiffmat11, Core::LinAlg::SerialDenseMatrix& stiffmat12,
        Core::LinAlg::SerialDenseMatrix& stiffmat21,
        Core::LinAlg::SerialDenseMatrix& stiffmat22) const;

    /** \brief compute contributions to analytic linearization (i.e. stiffness matrices) at current
     *         Gauss point pair
     *
     */
    void evaluate_stiffpot_analytic_contributions_large_sep_approx(
        Core::LinAlg::Matrix<3, 1, Sacado::Fad::DFad<double>> const& dist,
        Sacado::Fad::DFad<double> const& norm_dist, Sacado::Fad::DFad<double> const& norm_dist_exp1,
        double q1q2_JacFac_GaussWeights,
        Core::LinAlg::Matrix<1, numnodes * numnodalvalues, double> const& N1_i_GP1,
        Core::LinAlg::Matrix<1, numnodes * numnodalvalues, double> const& N2_i_GP2,
        Core::LinAlg::SerialDenseMatrix& stiffmat11, Core::LinAlg::SerialDenseMatrix& stiffmat12,
        Core::LinAlg::SerialDenseMatrix& stiffmat21,
        Core::LinAlg::SerialDenseMatrix& stiffmat22) const
    {
      // this is a dummy since for type FAD, no analytic linearization is required
    }

    /** \brief Evaluate forces and stiffness contribution resulting from potential-based interaction
     *         using double length specific potential with approximation for small separations
     *
     */
    void evaluate_fpotand_stiffpot_double_length_specific_small_sep_approx(
        Core::LinAlg::Matrix<3 * numnodes * numnodalvalues, 1, T>& force_pot1,
        Core::LinAlg::Matrix<3 * numnodes * numnodalvalues, 1, T>& force_pot2,
        Core::LinAlg::SerialDenseMatrix* stiffmat11, Core::LinAlg::SerialDenseMatrix* stiffmat12,
        Core::LinAlg::SerialDenseMatrix* stiffmat21, Core::LinAlg::SerialDenseMatrix* stiffmat22);


    /** \brief compute contributions to analytic linearization (i.e. stiffness matrices) at current
     *         Gauss point pair
     *
     */
    void evaluate_stiffpot_analytic_contributions_double_length_specific_small_sep_approx(
        Core::LinAlg::Matrix<3, 1, double> const& dist, double const& norm_dist, double const& gap,
        double const& gap_regularized, double const& gap_exp1, double q1q2_JacFac_GaussWeights,
        Core::LinAlg::Matrix<1, numnodes * numnodalvalues, double> const& N1_i_GP1,
        Core::LinAlg::Matrix<1, numnodes * numnodalvalues, double> const& N2_i_GP2,
        Core::LinAlg::SerialDenseMatrix& stiffmat11, Core::LinAlg::SerialDenseMatrix& stiffmat12,
        Core::LinAlg::SerialDenseMatrix& stiffmat21,
        Core::LinAlg::SerialDenseMatrix& stiffmat22) const;

    /** \brief compute contributions to analytic linearization (i.e. stiffness matrices) at current
     *         Gauss point pair
     *
     */
    void evaluate_stiffpot_analytic_contributions_double_length_specific_small_sep_approx(
        Core::LinAlg::Matrix<3, 1, Sacado::Fad::DFad<double>> const& dist,
        Sacado::Fad::DFad<double> const& norm_dist, Sacado::Fad::DFad<double> const& gap,
        Sacado::Fad::DFad<double> const& gap_regularized, Sacado::Fad::DFad<double> const& gap_exp1,
        double q1q2_JacFac_GaussWeights,
        Core::LinAlg::Matrix<1, numnodes * numnodalvalues, double> const& N1_i_GP1,
        Core::LinAlg::Matrix<1, numnodes * numnodalvalues, double> const& N2_i_GP2,
        Core::LinAlg::SerialDenseMatrix& stiffmat11, Core::LinAlg::SerialDenseMatrix& stiffmat12,
        Core::LinAlg::SerialDenseMatrix& stiffmat21,
        Core::LinAlg::SerialDenseMatrix& stiffmat22) const
    {
      // this is a dummy since for type FAD, no analytic linearization is required
    }

    /** \brief Evaluate forces and stiffness contribution resulting from potential-based interaction
     *         using single length specific potential with approximation for small separations
     *
     */
    void evaluate_fpotand_stiffpot_single_length_specific_small_sep_approx(
        Core::LinAlg::Matrix<3 * numnodes * numnodalvalues, 1, T>& force_pot1,
        Core::LinAlg::Matrix<3 * numnodes * numnodalvalues, 1, T>& force_pot2,
        Core::LinAlg::SerialDenseMatrix* stiffmat11, Core::LinAlg::SerialDenseMatrix* stiffmat12,
        Core::LinAlg::SerialDenseMatrix* stiffmat21, Core::LinAlg::SerialDenseMatrix* stiffmat22);

    /** \brief Evaluate all quantities for the full disk-cylinder potential law
     *
     */
    bool evaluate_full_disk_cylinder_potential(T& interaction_potential_GP,
        Core::LinAlg::Matrix<3 * numnodes * numnodalvalues, 1, T>& force_pot_slave_GP,
        Core::LinAlg::Matrix<3 * numnodes * numnodalvalues, 1, T>& force_pot_master_GP,
        Core::LinAlg::Matrix<3, 1, T> const& r_slave,
        Core::LinAlg::Matrix<3, 1, T> const& r_xi_slave,
        Core::LinAlg::Matrix<3, 1, T> const& t1_slave,
        Core::LinAlg::Matrix<3, 1, T> const& r_master,
        Core::LinAlg::Matrix<3, 1, T> const& r_xi_master,
        Core::LinAlg::Matrix<3, 1, T> const& r_xixi_master,
        Core::LinAlg::Matrix<3, 1, T> const& t1_master, T alpha, T cos_alpha,
        Core::LinAlg::Matrix<3, 1, T> const& dist_ul,
        Core::LinAlg::Matrix<1, 3, T> const& xi_master_partial_r_slave,
        Core::LinAlg::Matrix<1, 3, T> const& xi_master_partial_r_master,
        Core::LinAlg::Matrix<1, 3, T> const& xi_master_partial_r_xi_master,
        double prefactor_visualization_data,
        Core::LinAlg::Matrix<3, 1, double>& vtk_force_pot_slave_GP,
        Core::LinAlg::Matrix<3, 1, double>& vtk_force_pot_master_GP,
        Core::LinAlg::Matrix<3, 1, double>& vtk_moment_pot_slave_GP,
        Core::LinAlg::Matrix<3, 1, double>& vtk_moment_pot_master_GP,
        double rho1rho2_JacFac_GaussWeight,
        Core::LinAlg::Matrix<1, numnodes * numnodalvalues, double> const& N_i_slave,
        Core::LinAlg::Matrix<1, numnodes * numnodalvalues, double> const& N_i_xi_slave,
        Core::LinAlg::Matrix<1, numnodes * numnodalvalues, T> const& N_i_master,
        Core::LinAlg::Matrix<1, numnodes * numnodalvalues, T> const& N_i_xi_master);

    /** \brief Evaluate all quantities for the simple disk-cylinder potential law */
    bool evaluate_simple_disk_cylinder_potential(Core::LinAlg::Matrix<3, 1, T> const& dist_ul,
        T norm_dist_ul, T alpha, T cos_alpha, Core::LinAlg::Matrix<3, 1, T> const& r_slave,
        Core::LinAlg::Matrix<3, 1, T> const& r_xi_slave, T norm_r_xi_slave,
        Core::LinAlg::Matrix<3, 1, T> const& t_slave, Core::LinAlg::Matrix<3, 1, T> const& r_master,
        Core::LinAlg::Matrix<3, 1, T> const& r_xi_master, T norm_r_xi_master,
        Core::LinAlg::Matrix<3, 1, T> const& r_xixi_master,
        Core::LinAlg::Matrix<3, 1, T> const& t_master, T xi_master,
        Core::LinAlg::Matrix<1, 3, T> const& xi_master_partial_r_slave,
        Core::LinAlg::Matrix<1, 3, T> const& xi_master_partial_r_master,
        Core::LinAlg::Matrix<1, 3, T> const& xi_master_partial_r_xi_master,
        std::optional<double> potential_reduction_length, double length_prior_right,
        double length_prior_left, T& interaction_potential_GP,
        Core::LinAlg::Matrix<1, numnodes * numnodalvalues, double> const& N_i_slave,
        Core::LinAlg::Matrix<1, numnodes * numnodalvalues, double> const& N_i_xi_slave,
        Core::LinAlg::Matrix<1, numnodes * numnodalvalues, T> const& N_i_master,
        Core::LinAlg::Matrix<1, numnodes * numnodalvalues, T> const& N_i_xi_master,
        Core::LinAlg::Matrix<1, numnodes * numnodalvalues, T> const& N_i_xixi_master,
        Core::LinAlg::Matrix<3 * numnodes * numnodalvalues, 1, T>& force_pot_slave_GP,
        Core::LinAlg::Matrix<3 * numnodes * numnodalvalues, 1, T>& force_pot_master_GP,
        double prefactor_visualization_data, double rho1rho2_JacFac_GaussWeight,
        Core::LinAlg::Matrix<3, 1, double>& vtk_force_pot_slave_GP,
        Core::LinAlg::Matrix<3, 1, double>& vtk_force_pot_master_GP,
        Core::LinAlg::Matrix<3, 1, double>& vtk_moment_pot_slave_GP,
        Core::LinAlg::Matrix<3, 1, double>& vtk_moment_pot_master_GP,
        Core::LinAlg::SerialDenseMatrix* stiffmat11, Core::LinAlg::SerialDenseMatrix* stiffmat12,
        Core::LinAlg::SerialDenseMatrix* stiffmat21, Core::LinAlg::SerialDenseMatrix* stiffmat22);

    /** \brief scale given stiffness matrices by given scalar factor
     *
     */
    void scale_stiffpot_analytic_contributions_if_required(double const& scalefactor,
        Core::LinAlg::SerialDenseMatrix& stiffmat11, Core::LinAlg::SerialDenseMatrix& stiffmat12,
        Core::LinAlg::SerialDenseMatrix& stiffmat21,
        Core::LinAlg::SerialDenseMatrix& stiffmat22) const;

    /** \brief scale given stiffness matrices by given scalar factor
     *
     */
    void scale_stiffpot_analytic_contributions_if_required(
        Sacado::Fad::DFad<double> const& scalefactor, Core::LinAlg::SerialDenseMatrix& stiffmat11,
        Core::LinAlg::SerialDenseMatrix& stiffmat12, Core::LinAlg::SerialDenseMatrix& stiffmat21,
        Core::LinAlg::SerialDenseMatrix& stiffmat22) const
    {
      // this is a dummy since for type FAD, no analytic linearization is required
    }

    /** \brief compute linearization (i.e. stiffness matrices) of given force vectors
     *         using automatic differentiation based on Sacado::Fad package
     *
     */
    void calc_stiffmat_automatic_differentiation_if_required(
        Core::LinAlg::Matrix<3 * numnodes * numnodalvalues, 1, double> const& force_pot1,
        Core::LinAlg::Matrix<3 * numnodes * numnodalvalues, 1, double> const& force_pot2,
        Core::LinAlg::SerialDenseMatrix& stiffmat11, Core::LinAlg::SerialDenseMatrix& stiffmat12,
        Core::LinAlg::SerialDenseMatrix& stiffmat21,
        Core::LinAlg::SerialDenseMatrix& stiffmat22) const
    {
      // this is a dummy since for type double, no automatic differentiation is required
    }

    /** \brief compute linearization (i.e. stiffness matrices) of given force vectors
     *         using automatic differentiation based on Sacado::Fad package
     *
     */
    void calc_stiffmat_automatic_differentiation_if_required(
        Core::LinAlg::Matrix<3 * numnodes * numnodalvalues, 1, Sacado::Fad::DFad<double>> const&
            force_pot1,
        Core::LinAlg::Matrix<3 * numnodes * numnodalvalues, 1, Sacado::Fad::DFad<double>> const&
            force_pot2,
        Core::LinAlg::SerialDenseMatrix& stiffmat11, Core::LinAlg::SerialDenseMatrix& stiffmat12,
        Core::LinAlg::SerialDenseMatrix& stiffmat21,
        Core::LinAlg::SerialDenseMatrix& stiffmat22) const;

    /** \brief add contributions from linearization of parameter coordinate on master beam xi_master
     *         if determined via point-to-curve projection
     *         using automatic differentiation based on Sacado::Fad package
     *
     */
    void add_stiffmat_contributions_xi_master_automatic_differentiation_if_required(
        Core::LinAlg::Matrix<3 * numnodes * numnodalvalues, 1, double> const& force_pot1,
        Core::LinAlg::Matrix<3 * numnodes * numnodalvalues, 1, double> const& force_pot2,
        Core::LinAlg::Matrix<1, 3 * numnodes * numnodalvalues, double> const&
            lin_xi_master_slaveDofs,
        Core::LinAlg::Matrix<1, 3 * numnodes * numnodalvalues, double> const&
            lin_xi_master_masterDofs,
        Core::LinAlg::SerialDenseMatrix& stiffmat11, Core::LinAlg::SerialDenseMatrix& stiffmat12,
        Core::LinAlg::SerialDenseMatrix& stiffmat21,
        Core::LinAlg::SerialDenseMatrix& stiffmat22) const
    {
      // this is a dummy since for type double, no automatic differentiation is required
    }

    /** \brief add contributions from linearization of parameter coordinate on master beam xi_master
     *         if determined via point-to-curve projection
     *         using automatic differentiation based on Sacado::Fad package
     *
     */
    void add_stiffmat_contributions_xi_master_automatic_differentiation_if_required(
        Core::LinAlg::Matrix<3 * numnodes * numnodalvalues, 1, Sacado::Fad::DFad<double>> const&
            force_pot1,
        Core::LinAlg::Matrix<3 * numnodes * numnodalvalues, 1, Sacado::Fad::DFad<double>> const&
            force_pot2,
        Core::LinAlg::Matrix<1, 3 * numnodes * numnodalvalues, Sacado::Fad::DFad<double>> const&
            lin_xi_master_slaveDofs,
        Core::LinAlg::Matrix<1, 3 * numnodes * numnodalvalues, Sacado::Fad::DFad<double>> const&
            lin_xi_master_masterDofs,
        Core::LinAlg::SerialDenseMatrix& stiffmat11, Core::LinAlg::SerialDenseMatrix& stiffmat12,
        Core::LinAlg::SerialDenseMatrix& stiffmat21,
        Core::LinAlg::SerialDenseMatrix& stiffmat22) const;

    /** \brief compute discrete force vectors using automatic differentiation
     *
     */
    void calc_fpot_gausspoint_automatic_differentiation_if_required(
        Core::LinAlg::Matrix<3 * numnodes * numnodalvalues, 1, double>& force_pot1,
        Core::LinAlg::Matrix<3 * numnodes * numnodalvalues, 1, double>& force_pot2,
        Sacado::Fad::DFad<double> const& interaction_potential,
        Sacado::Fad::DFad<double> const& potential_reduction_factor,
        Core::LinAlg::Matrix<1, 3 * numnodes * numnodalvalues, Sacado::Fad::DFad<double>> const&
            lin_xi_master_slaveDofs,
        Core::LinAlg::Matrix<1, 3 * numnodes * numnodalvalues, Sacado::Fad::DFad<double>> const&
            lin_xi_master_masterDofs) const;

    /** \brief compute discrete force vectors using automatic differentiation
     *
     */
    void calc_fpot_gausspoint_automatic_differentiation_if_required(
        Core::LinAlg::Matrix<3 * numnodes * numnodalvalues, 1, Sacado::Fad::DFad<double>>&
            force_pot1,
        Core::LinAlg::Matrix<3 * numnodes * numnodalvalues, 1, Sacado::Fad::DFad<double>>&
            force_pot2,
        Sacado::Fad::DFad<double> const& interaction_potential,
        Sacado::Fad::DFad<double> const& potential_reduction_factor,
        Core::LinAlg::Matrix<1, 3 * numnodes * numnodalvalues, Sacado::Fad::DFad<double>> const&
            lin_xi_master_slaveDofs,
        Core::LinAlg::Matrix<1, 3 * numnodes * numnodalvalues, Sacado::Fad::DFad<double>> const&
            lin_xi_master_masterDofs) const;

    /** \brief compute discrete force vectors using automatic differentiation
     *
     */
    void calc_fpot_gausspoint_automatic_differentiation_if_required(

        Core::LinAlg::Matrix<3 * numnodes * numnodalvalues, 1, double>& force_pot1,
        Core::LinAlg::Matrix<3 * numnodes * numnodalvalues, 1, double>& force_pot2,
        double const& interaction_potential, double const& potential_reduction_factor,
        Core::LinAlg::Matrix<1, 3 * numnodes * numnodalvalues, double> const&
            lin_xi_master_slaveDofs,
        Core::LinAlg::Matrix<1, 3 * numnodes * numnodalvalues, double> const&
            lin_xi_master_masterDofs) const
    {
      // this is a dummy since for type double, no automatic differentiation is available
    }

    /** \brief evaluate analytic linearization (i.e. stiffness matrices) at current Gauss point
     *
     *  \param stiffmat11 d (Residuum_vec_1) / d (dof_vec_1)
     *
     */
    void evaluate_stiffpot_analytic_contributions_single_length_specific_small_sep_approx_simple(
        Core::LinAlg::Matrix<1, numnodes * numnodalvalues, double> const& N_i_slave,
        Core::LinAlg::Matrix<1, numnodes * numnodalvalues, double> const& N_i_xi_slave,
        Core::LinAlg::Matrix<1, numnodes * numnodalvalues, double> const& N_i_master,
        Core::LinAlg::Matrix<1, numnodes * numnodalvalues, double> const& N_i_xi_master,
        Core::LinAlg::Matrix<1, numnodes * numnodalvalues, double> const& N_i_xixi_master,
        double const& xi_master, Core::LinAlg::Matrix<3, 1, double> const& r_xi_slave,
        Core::LinAlg::Matrix<3, 1, double> const& r_xi_master,
        Core::LinAlg::Matrix<3, 1, double> const& r_xixi_master, double const& norm_dist_ul,
        Core::LinAlg::Matrix<3, 1, double> const& normal_ul, double const& pot_red_fac,
        double const& pot_red_fac_deriv_xi_master, double const& pot_red_fac_2ndderiv_xi_master,
        double const& pot_ia, double const& pot_ia_deriv_gap_ul,
        double const& pot_ia_deriv_cos_alpha, double const& pot_ia_2ndderiv_gap_ul,
        double const& pot_ia_deriv_gap_ul_deriv_cos_alpha, double const& pot_ia_2ndderiv_cos_alpha,
        Core::LinAlg::Matrix<3, 1, double> const& gap_ul_deriv_r_slave,
        Core::LinAlg::Matrix<3, 1, double> const& gap_ul_deriv_r_master,
        Core::LinAlg::Matrix<3, 1, double> const& cos_alpha_deriv_r_slave,
        Core::LinAlg::Matrix<3, 1, double> const& cos_alpha_deriv_r_master,
        Core::LinAlg::Matrix<3, 1, double> const& cos_alpha_deriv_r_xi_slave,
        Core::LinAlg::Matrix<3, 1, double> const& cos_alpha_deriv_r_xi_master,
        Core::LinAlg::Matrix<1, 3, double> const& xi_master_partial_r_slave,
        Core::LinAlg::Matrix<1, 3, double> const& xi_master_partial_r_master,
        Core::LinAlg::Matrix<1, 3, double> const& xi_master_partial_r_xi_master,
        Core::LinAlg::SerialDenseMatrix& stiffmat11, Core::LinAlg::SerialDenseMatrix& stiffmat12,
        Core::LinAlg::SerialDenseMatrix& stiffmat21,
        Core::LinAlg::SerialDenseMatrix& stiffmat22) const;

    /** \brief evaluate analytic linearization (i.e. stiffness matrices) at current Gauss point
     *
     */
    void evaluate_stiffpot_analytic_contributions_single_length_specific_small_sep_approx_simple(
        Core::LinAlg::Matrix<1, numnodes * numnodalvalues, double> const& N_i_slave,
        Core::LinAlg::Matrix<1, numnodes * numnodalvalues, double> const& N_i_xi_slave,
        Core::LinAlg::Matrix<1, numnodes * numnodalvalues, Sacado::Fad::DFad<double>> const&
            N_i_master,
        Core::LinAlg::Matrix<1, numnodes * numnodalvalues, Sacado::Fad::DFad<double>> const&
            N_i_xi_master,
        Core::LinAlg::Matrix<1, numnodes * numnodalvalues, Sacado::Fad::DFad<double>> const&
            N_i_xixi_master,
        Sacado::Fad::DFad<double> const& xi_master,
        Core::LinAlg::Matrix<3, 1, Sacado::Fad::DFad<double>> const& r_xi_slave,
        Core::LinAlg::Matrix<3, 1, Sacado::Fad::DFad<double>> const& r_xi_master,
        Core::LinAlg::Matrix<3, 1, Sacado::Fad::DFad<double>> const& r_xixi_master,
        Sacado::Fad::DFad<double> const& norm_dist_ul,
        Core::LinAlg::Matrix<3, 1, Sacado::Fad::DFad<double>> const& normal_ul,
        Sacado::Fad::DFad<double> const& pot_red_fac,
        Sacado::Fad::DFad<double> const& pot_red_fac_deriv_xi_master,
        Sacado::Fad::DFad<double> const& pot_red_fac_2ndderiv_xi_master,
        Sacado::Fad::DFad<double> const& pot_ia,
        Sacado::Fad::DFad<double> const& pot_ia_deriv_gap_ul,
        Sacado::Fad::DFad<double> const& pot_ia_deriv_cos_alpha,
        Sacado::Fad::DFad<double> const& pot_ia_2ndderiv_gap_ul,
        Sacado::Fad::DFad<double> const& pot_ia_deriv_gap_ul_deriv_cos_alpha,
        Sacado::Fad::DFad<double> const& pot_ia_2ndderiv_cos_alpha,
        Core::LinAlg::Matrix<3, 1, Sacado::Fad::DFad<double>> const& gap_ul_deriv_r_slave,
        Core::LinAlg::Matrix<3, 1, Sacado::Fad::DFad<double>> const& gap_ul_deriv_r_master,
        Core::LinAlg::Matrix<3, 1, Sacado::Fad::DFad<double>> const& cos_alpha_deriv_r_slave,
        Core::LinAlg::Matrix<3, 1, Sacado::Fad::DFad<double>> const& cos_alpha_deriv_r_master,
        Core::LinAlg::Matrix<3, 1, Sacado::Fad::DFad<double>> const& cos_alpha_deriv_r_xi_slave,
        Core::LinAlg::Matrix<3, 1, Sacado::Fad::DFad<double>> const& cos_alpha_deriv_r_xi_master,
        Core::LinAlg::Matrix<1, 3, Sacado::Fad::DFad<double>> const& xi_master_partial_r_slave,
        Core::LinAlg::Matrix<1, 3, Sacado::Fad::DFad<double>> const& xi_master_partial_r_master,
        Core::LinAlg::Matrix<1, 3, Sacado::Fad::DFad<double>> const& xi_master_partial_r_xi_master,
        Core::LinAlg::SerialDenseMatrix& stiffmat11, Core::LinAlg::SerialDenseMatrix& stiffmat12,
        Core::LinAlg::SerialDenseMatrix& stiffmat21,
        Core::LinAlg::SerialDenseMatrix& stiffmat22) const
    {
      // this is a dummy since for type FAD, no analytic linearization is required
    }

    /*!
    \brief Calculate shape function values for given parameter values
    Todo call more general utils method
    */
    void get_shape_functions(
        std::vector<Core::LinAlg::Matrix<1, numnodes * numnodalvalues, double>>& N1_i,
        std::vector<Core::LinAlg::Matrix<1, numnodes * numnodalvalues, double>>& N2_i,
        std::vector<Core::LinAlg::Matrix<1, numnodes * numnodalvalues, double>>& N1_i_xi,
        std::vector<Core::LinAlg::Matrix<1, numnodes * numnodalvalues, double>>& N2_i_xi,
        Core::FE::IntegrationPoints1D& gausspoints) const;

    /*!
    \brief Compute coordinates of centerline point
    */
    template <typename T2>
    void compute_centerline_position(Core::LinAlg::Matrix<3, 1, T>& r,
        const Core::LinAlg::Matrix<1, numnodes * numnodalvalues, T2>& N_i,
        const Core::LinAlg::Matrix<3 * numnodes * numnodalvalues, 1, T> eledofvec) const;

    /*!
    \brief Compute tangent vector at centerline point
    */
    template <typename T2>
    void compute_centerline_tangent(Core::LinAlg::Matrix<3, 1, T>& r_xi,
        const Core::LinAlg::Matrix<1, numnodes * numnodalvalues, T2>& N_i_xi,
        const Core::LinAlg::Matrix<3 * numnodes * numnodalvalues, 1, T> eledofvec) const;

    /** \brief set primary variables for FAD if required
     *
     */
    void set_automatic_differentiation_variables_if_required(
        Core::LinAlg::Matrix<3 * numnodes * numnodalvalues, 1, double>& ele1centerlinedofvec,
        Core::LinAlg::Matrix<3 * numnodes * numnodalvalues, 1, double>& ele2centerlinedofvec)
    {
      // do nothing in case of type double (analytic differentiation)
    }

    /** \brief set primary variables for FAD if required
     *
     */
    void set_automatic_differentiation_variables_if_required(
        Core::LinAlg::Matrix<3 * numnodes * numnodalvalues, 1, Sacado::Fad::DFad<double>>&
            ele1centerlinedofvec,
        Core::LinAlg::Matrix<3 * numnodes * numnodalvalues, 1, Sacado::Fad::DFad<double>>&
            ele2centerlinedofvec);

    /** \brief set primary variables including xi_master for FAD if required
     *
     */
    void set_automatic_differentiation_variables_if_required(
        Core::LinAlg::Matrix<3 * numnodes * numnodalvalues, 1, double>& ele1centerlinedofvec,
        Core::LinAlg::Matrix<3 * numnodes * numnodalvalues, 1, double>& ele2centerlinedofvec,
        double& xi_master)
    {
      // do nothing in case of type double (analytic differentiation)
    }

    /** \brief set primary variables including xi_master for FAD if required
     *
     */
    void set_automatic_differentiation_variables_if_required(
        Core::LinAlg::Matrix<3 * numnodes * numnodalvalues, 1, Sacado::Fad::DFad<double>>&
            ele1centerlinedofvec,
        Core::LinAlg::Matrix<3 * numnodes * numnodalvalues, 1, Sacado::Fad::DFad<double>>&
            ele2centerlinedofvec,
        Sacado::Fad::DFad<double>& xi_master);

    /** \brief estimate whether the elements' separation is much more than the cutoff distance
     *
     *  The idea is to get a conservative lower bound estimate for the minimal centerline separation
     *  of both elements, which doesn't require any interpolation, i.e., information on the
     *  deformation of the elements and is thus cheap to compute before we begin with the actual,
     *  expensive evaluation of the pair on Gauss point level.
     *
     */
    bool are_elements_much_more_separated_than_cutoff_distance();

    //@}

   private:
    //! @name member variables

    //! first element of pair
    Discret::Elements::Beam3Base const* beam_element1_;

    //! second element of pair
    Discret::Elements::Beam3Base const* beam_element2_;

    //! line charge conditions
    std::vector<const Core::Conditions::Condition*> linechargeconds_;

    //! current time
    double time_;

    //! current node coordinates of the two elements
    Core::LinAlg::Matrix<3 * numnodes * numnodalvalues, 1, T> ele1pos_;
    Core::LinAlg::Matrix<3 * numnodes * numnodalvalues, 1, T> ele2pos_;

    //! parameters of the applied (point-point) potential law Phi(r)=k_ * r^(-m_)
    double k_;
    double m_;

    //! initial element lengths
    double ele1length_;
    double ele2length_;

    //! Cross-section radius of beam 1
    double radius1_;

    //! Cross-section radius of beam 2
    double radius2_;
    //@}

    //! @name data storage for visualization output

    //! centerline coordinate vector of interacting points of element 1 and 2
    std::vector<Core::LinAlg::Matrix<3, 1, double>> centerline_coords_gp_1_;
    std::vector<Core::LinAlg::Matrix<3, 1, double>> centerline_coords_gp_2_;

    //! resulting forces at interacting points of element 1 and 2
    std::vector<Core::LinAlg::Matrix<3, 1, double>> forces_pot_gp_1_;
    std::vector<Core::LinAlg::Matrix<3, 1, double>> forces_pot_gp_2_;

    //! resulting moments at interacting points of element 1 and 2
    std::vector<Core::LinAlg::Matrix<3, 1, double>> moments_pot_gp_1_;
    std::vector<Core::LinAlg::Matrix<3, 1, double>> moments_pot_gp_2_;

    //! total interaction potential of this pair
    double interaction_potential_;
    //@}
  };
}  // namespace BeamInteraction

FOUR_C_NAMESPACE_CLOSE

#endif
