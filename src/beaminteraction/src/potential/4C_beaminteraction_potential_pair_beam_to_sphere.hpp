// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_BEAMINTERACTION_POTENTIAL_PAIR_BEAM_TO_SPHERE_HPP
#define FOUR_C_BEAMINTERACTION_POTENTIAL_PAIR_BEAM_TO_SPHERE_HPP

#include "4C_config.hpp"

#include "4C_beaminteraction_beam_to_beam_contact_defines.hpp"
#include "4C_beaminteraction_beam_to_beam_contact_tangentsmoothing.hpp"
#include "4C_beaminteraction_potential_pair_base.hpp"
#include "4C_fem_general_utils_integration.hpp"
#include "4C_fem_general_utils_local_connectivity_matrices.hpp"
#include "4C_linalg_fixedsizematrix.hpp"
#include "4C_linalg_sparsematrix.hpp"

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
    class Rigidsphere;
    class Beam3Base;
  }  // namespace Elements
}  // namespace Discret

namespace BeamInteraction
{
  /*!
   \brief class for potential-based interaction of a 3D beam element and a rigid sphere
   */

  template <unsigned int numnodes, unsigned int numnodalvalues>
  class BeamToSpherePotentialPair : public BeamInteraction::BeamPotentialPair
  {
   public:
    //! @name Friends
    // no friend classes defined
    //@}

    //! @name Constructors and destructors and related methods
    /*!
    \brief Standard Constructor
    */
    BeamToSpherePotentialPair();



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
      FOUR_C_THROW("not implemented yet");
    }

    void get_all_interacting_point_coords_element2(
        std::vector<Core::LinAlg::Matrix<3, 1, double>>& coords) const override
    {
      FOUR_C_THROW("not implemented yet");
    }

    /*!
    \brief Get forces at all interacting points on element1 and element2
    */
    void get_forces_at_all_interacting_points_element1(
        std::vector<Core::LinAlg::Matrix<3, 1, double>>& forces) const override
    {
      FOUR_C_THROW("not implemented yet");
    }

    void get_forces_at_all_interacting_points_element2(
        std::vector<Core::LinAlg::Matrix<3, 1, double>>& forces) const override
    {
      FOUR_C_THROW("not implemented yet");
    }

    /*!
    \brief Get moments at all interacting points on element1 and element2
    */
    void get_moments_at_all_interacting_points_element1(
        std::vector<Core::LinAlg::Matrix<3, 1, double>>& moments) const override
    {
      FOUR_C_THROW("not implemented yet");
    }

    void get_moments_at_all_interacting_points_element2(
        std::vector<Core::LinAlg::Matrix<3, 1, double>>& moments) const override
    {
      FOUR_C_THROW("not implemented yet");
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
    \brief Get first element (beam)
    */
    inline const Discret::Elements::Beam3Base* beam_element() { return beam_element_; };

    /*!
    \brief Get second element (sphere)
    */
    inline const Discret::Elements::Rigidsphere* sphere_element() { return sphere_element_; };
    //@}

   private:
    //! @name member variables

    //! first element of pair
    Discret::Elements::Beam3Base const* beam_element_;

    //! second element of pair
    Discret::Elements::Rigidsphere const* sphere_element_;

    //! line and point charge condition
    std::vector<const Core::Conditions::Condition*> chargeconds_;

    //! current time
    double time_;

    //! current absolute Dof values of the two elements
    Core::LinAlg::Matrix<3 * numnodes * numnodalvalues, 1, TYPE> ele1pos_;
    Core::LinAlg::Matrix<3, 1, TYPE> ele2pos_;

    //! parameters of the applied (point-point) potential law Phi(r)=k_ * r^(-m_)
    double k_;
    double m_;

    //! beam element arc-length in stress-free reference configuration
    double beamele_reflength_;

    //! Cross-section radius of beam
    double radius1_;

    //! Cross-section radius of sphere
    double radius2_;

    //! resulting forces on element 1 and 2
    Core::LinAlg::Matrix<3 * numnodes * numnodalvalues, 1, TYPE> fpot1_;
    Core::LinAlg::Matrix<3, 1, TYPE> fpot2_;

    //! stiffness contributions
    Core::LinAlg::Matrix<3 * numnodes * numnodalvalues, 3 * numnodes * numnodalvalues + 3, TYPE>
        stiffpot1_;
    Core::LinAlg::Matrix<3, 3 * numnodes * numnodalvalues + 3, TYPE> stiffpot2_;

    //! total interaction potential of this pair
    double interaction_potential_;

    //@}

    //! @name Private evaluation methods

    /*!
    \brief Evaluate forces and stiffness contribution resulting from potential-based interaction
    */
    void evaluate_fpotand_stiffpot_large_sep_approx();

    /*!
    \brief Calculate shape function values for given parameter values
        Todo call more general utils method
    */
    void get_shape_functions(std::vector<Core::LinAlg::Matrix<1, numnodes * numnodalvalues>>& N1_i,
        std::vector<Core::LinAlg::Matrix<1, numnodes * numnodalvalues>>& N1_i_xi,
        Core::FE::IntegrationPoints1D& gausspoints);

    /*!
    \brief Compute coordinates of centreline points from the discretization
        Todo call more general utils method
    */
    void compute_coords(Core::LinAlg::Matrix<3, 1, TYPE>& r,
        const Core::LinAlg::Matrix<1, numnodes * numnodalvalues>& N_i,
        const Core::LinAlg::Matrix<3 * numnodes * numnodalvalues, 1, TYPE> elepos);

    //@}
  };
}  // namespace BeamInteraction

FOUR_C_NAMESPACE_CLOSE

#endif
