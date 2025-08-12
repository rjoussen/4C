// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_BEAMINTERACTION_BEAM_TO_SPHERE_CONTACT_PAIR_HPP
#define FOUR_C_BEAMINTERACTION_BEAM_TO_SPHERE_CONTACT_PAIR_HPP

#include "4C_config.hpp"

#include "4C_beaminteraction_beam_to_beam_contact_defines.hpp"
#include "4C_beaminteraction_contact_pair.hpp"
#include "4C_fem_condition.hpp"
#include "4C_fem_general_element.hpp"
#include "4C_fem_general_node.hpp"
#include "4C_fem_general_utils_integration.hpp"
#include "4C_linalg_fixedsizematrix.hpp"
#include "4C_linalg_serialdensematrix.hpp"
#include "4C_linalg_serialdensevector.hpp"
#include "4C_linalg_sparsematrix.hpp"
#include "4C_linalg_vector.hpp"
#include "4C_utils_fad.hpp"

#include <memory>

FOUR_C_NAMESPACE_OPEN

// forward declaration ...
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
  template <unsigned int numnodes, unsigned int numnodalvalues>
  class BeamToSphereContactPair : public BeamContactPair
  {
   public:
    //! @name Friends
    // no friend classes defined
    //@}

    //! @name Constructors and destructors and related methods

    /*!
    \brief Standard Constructor
    */
    BeamToSphereContactPair();



    //! Setup
    void setup() override;

    //@}

    //! @name Derived methods from base class
    /*!
    \brief things that need to be done in a separate loop before the actual evaluation loop
           over all contact pairs
    */
    void pre_evaluate() override;

    /*!
    \brief Evaluate this contact element pair, return value indicates whether pair is active,
           i.e. non-zero values for force and stiffmat are returned
    */
    bool evaluate(Core::LinAlg::SerialDenseVector* forcevec1,
        Core::LinAlg::SerialDenseVector* forcevec2, Core::LinAlg::SerialDenseMatrix* stiffmat11,
        Core::LinAlg::SerialDenseMatrix* stiffmat12, Core::LinAlg::SerialDenseMatrix* stiffmat21,
        Core::LinAlg::SerialDenseMatrix* stiffmat22) override;

    /*
    \brief Update state of translational nodal DoFs (absolute positions (and tangents)) of both
    elements
    */
    void reset_state(const std::vector<double>& centerline_dofvec_ele1,
        const std::vector<double>& centerline_dofvec_ele2) override;

    /** \brief print information about this beam contact element pair to screen
     *
     */
    void print(std::ostream& out) const override;

    /**
     * \brief Returns the type of this beam to sphere contact pair.
     */
    ContactPairType get_type() const override final
    {
      return ContactPairType::beam_to_sphere_contact;
    }

    /** \brief print this beam contact element pair to screen
     *
     */
    void print_summary_one_line_per_active_segment_pair(std::ostream& out) const override;
    //@}

    //! @name Access methods

    /*!
    \brief Get beam element
    */
    inline Discret::Elements::Beam3Base const* beam_element() { return beam_element_; };

    /*!
    \brief Get sphere element
    */
    inline Discret::Elements::Rigidsphere const* sphere_element() { return sphere_element_; };

    /*!
    \brief Get flag indicating whether contact is active (true) or inactive (false)
    */
    inline bool get_contact_flag() const override
    {
      return (contactflag_ or nodalcontactflag_[0] or nodalcontactflag_[1]);
    };

    /*!
    \brief Get number of active contact point pairs on this element pair
    */
    unsigned int get_num_all_active_contact_point_pairs() const override { return 1; }


    /*!
    \brief Get coordinates of all active contact points on element1 and element2
    */
    inline void get_all_active_contact_point_coords_element1(
        std::vector<Core::LinAlg::Matrix<3, 1, double>>& coords) const override
    {
      FOUR_C_THROW("not implemented yet!");
    }

    inline void get_all_active_contact_point_coords_element2(
        std::vector<Core::LinAlg::Matrix<3, 1, double>>& coords) const override
    {
      FOUR_C_THROW("not implemented yet!");
    }

    /*!
    \brief Get all (scalar) contact forces of this contact pair
    */
    inline void get_all_active_contact_forces(std::vector<double>& forces,
        std::vector<double>& angles, std::vector<int>& types) const override
    {
      FOUR_C_THROW("not implemented yet!");
    }

    /*!
    \brief Get all (scalar) gap values of this contact pair
    */
    void get_all_active_contact_gaps(std::vector<double>& gaps) const override
    {
      FOUR_C_THROW("not implemented yet!");
    }

    /*!
    \brief Get energy of penalty contact.
    */
    double get_energy() const override
    {
      FOUR_C_THROW("not implemented yet!");
      return 0.0;
    }
    //@}


   private:
    //! @name Private evaluation methods

    /*!
    \brief Find contact point via closest point projection
    */
    void closest_point_projection();

    /*!
    \brief Utility method for CPP (evaluate nonlinear function f)
    */
    void evaluate_orthogonality_condition(TYPE& f, const Core::LinAlg::Matrix<3, 1, TYPE>& delta_x,
        const double norm_delta_x, const Core::LinAlg::Matrix<3, 1, TYPE>& dx1);

    /*!
    \brief Utility method for CPP (evaluate Jacobian of nonlinear function f)
    */
    void evaluate_lin_orthogonality_condition(TYPE& df, Core::LinAlg::Matrix<3, 1, TYPE>& delta_x,
        const double norm_delta_x, const Core::LinAlg::Matrix<3, 1, TYPE>& dx1,
        const Core::LinAlg::Matrix<3, 1, TYPE>& ddx1);

    /*!
    \brief Evaluate and assemble contact forces
    */
    void evaluate_fc_contact(Core::LinAlg::SerialDenseVector& forcevec1,
        Core::LinAlg::SerialDenseVector& forcevec2, const double& pp, const TYPE& gap,
        const Core::LinAlg::Matrix<3, 1, TYPE>& normal,
        const Core::LinAlg::Matrix<1, numnodes * numnodalvalues, TYPE>& N1_i,
        const bool contactactive);

    /*!
    \brief Evaluate and assemble contact stiffness
    */
    void evaluate_stiffc_contact(Core::LinAlg::SerialDenseMatrix& stiffmat11,
        Core::LinAlg::SerialDenseMatrix& stiffmat12, Core::LinAlg::SerialDenseMatrix& stiffmat21,
        Core::LinAlg::SerialDenseMatrix& stiffmat22, const double& pp, const TYPE& gap,
        const Core::LinAlg::Matrix<3, 1, TYPE>& normal, const TYPE& norm,
        const Core::LinAlg::Matrix<3, 1, TYPE>& x1, const Core::LinAlg::Matrix<3, 1, TYPE>& x2,
        const Core::LinAlg::Matrix<3, 1, TYPE>& dx1, const Core::LinAlg::Matrix<3, 1, TYPE>& ddx1,
        const Core::LinAlg::Matrix<1, numnodes * numnodalvalues, TYPE>& N1_i,
        const Core::LinAlg::Matrix<1, numnodes * numnodalvalues, TYPE>& N1_i_xi,
        const Core::LinAlg::Matrix<1, numnodes * numnodalvalues, TYPE>& N1_i_xixi,
        bool activecontact, bool linxi = true);

    /*!
    \brief Compute normal vector in contact point
    */
    void compute_normal(Core::LinAlg::Matrix<3, 1, TYPE>& normal, TYPE& gap, TYPE& norm,
        const Core::LinAlg::Matrix<3, 1, TYPE>& x1, const Core::LinAlg::Matrix<3, 1, TYPE>& x2);

    /*!
    \brief Evaluate gap function
    */
    void compute_gap(TYPE& gap, const TYPE& norm);

    /*!
    \brief Compute radius of cross section based on moment of inertia
    */
    void compute_ele_radius(double& radius, const double& moi);

    /*!
    \brief Compute coordinates and their derivatives from the discretization
    */
    void compute_coords_and_derivs(Core::LinAlg::Matrix<3, 1, TYPE>& x1,
        Core::LinAlg::Matrix<3, 1, TYPE>& x2, Core::LinAlg::Matrix<3, 1, TYPE>& dx1,
        Core::LinAlg::Matrix<3, 1, TYPE>& ddx1,
        const Core::LinAlg::Matrix<1, numnodes * numnodalvalues, TYPE>& N1_i,
        const Core::LinAlg::Matrix<1, numnodes * numnodalvalues, TYPE>& N1_i_xi,
        const Core::LinAlg::Matrix<1, numnodes * numnodalvalues, TYPE>& N1_i_xixi);

    /*!
    \brief Get shape functions and their derivatives at eta
    */
    void get_shape_functions(Core::LinAlg::Matrix<1, numnodes * numnodalvalues, TYPE>& N1_i,
        Core::LinAlg::Matrix<1, numnodes * numnodalvalues, TYPE>& N1_i_xi,
        Core::LinAlg::Matrix<1, numnodes * numnodalvalues, TYPE>& N1_i_xixi, const TYPE& eta);

    /*!
    \brief Compute linearizations of contact point
    */
    void compute_lin_xi(Core::LinAlg::Matrix<3 * numnodes * numnodalvalues + 3, 1, TYPE>& delta_xi,
        const Core::LinAlg::Matrix<3, 1, TYPE>& x1, const Core::LinAlg::Matrix<3, 1, TYPE>& x2,
        const Core::LinAlg::Matrix<3, 1, TYPE>& dx1, const Core::LinAlg::Matrix<3, 1, TYPE>& ddx1,
        const Core::LinAlg::Matrix<1, numnodes * numnodalvalues, TYPE>& N1_i,
        const Core::LinAlg::Matrix<1, numnodes * numnodalvalues, TYPE>& N1_i_xi);

    /*!
    \brief Compute linearization of gap
    */
    void compute_lin_gap(
        Core::LinAlg::Matrix<3 * numnodes * numnodalvalues + 3, 1, TYPE>& delta_gap,
        Core::LinAlg::Matrix<3 * numnodes * numnodalvalues + 3, 1, TYPE>& delta_xi,
        const Core::LinAlg::Matrix<3, 1, TYPE>& x1, const Core::LinAlg::Matrix<3, 1, TYPE>& x2,
        const Core::LinAlg::Matrix<3, 1, TYPE>& dx1,
        const Core::LinAlg::Matrix<1, numnodes * numnodalvalues, TYPE>& N1_i, const TYPE& normdist,
        const Core::LinAlg::Matrix<3, 1, TYPE>& normal, const TYPE& norm, const TYPE& gap,
        Core::LinAlg::Matrix<3, 3 * numnodes * numnodalvalues + 3, TYPE>& delta_x1_minus_x2);

    /*!
    \brief Compute linearization of normal
    */

    void compute_lin_normal(
        Core::LinAlg::Matrix<3, 3 * numnodes * numnodalvalues + 3, TYPE>& delta_normal,
        const Core::LinAlg::Matrix<3 * numnodes * numnodalvalues + 3, 1, TYPE>& delta_xi,
        const Core::LinAlg::Matrix<3, 1, TYPE>& normal, const TYPE& norm_delta_x,
        const Core::LinAlg::Matrix<3, 1, TYPE>& x1_xi,
        const Core::LinAlg::Matrix<1, numnodes * numnodalvalues, TYPE>& N1_i);

    /*!
    \brief Compute normal contact disctance
    */
    void compute_distance(Core::LinAlg::Matrix<3, 1, TYPE>& distance, TYPE& normdist,
        const Core::LinAlg::Matrix<3, 1, TYPE>& normal, const TYPE& norm);

    /*!
    \brief Check if contact is active and set flag accordingly
    */
    void check_and_set_contact_status();

    //@}

   private:
    //! @name member variables

    //! first element of contact pair
    Discret::Elements::Beam3Base const* beam_element_;

    //! second element of contact pair
    Discret::Elements::Rigidsphere const* sphere_element_;

    //! current node coordinates of the two elements
    Core::LinAlg::Matrix<3 * numnodes * numnodalvalues, 1, TYPE> ele1pos_;
    Core::LinAlg::Matrix<3, 1, TYPE> ele2pos_;

    //! beam element arc-length in stress-free reference configuration
    double beamele_reflength_;

    //! Cross-section radius of beam
    double radius1_;

    //! Cross-section radius of sphere
    double radius2_;

    //! gap function
    TYPE gap_;

    //! flag indicating contact (active/inactive)
    bool contactflag_;

    //! flag indicating contact of beam end points (nodes) (active/inactive)
    std::vector<bool> nodalcontactflag_;

    //! resulting nodal contact forces on ele 1/2
    Core::LinAlg::Matrix<3 * numnodes * numnodalvalues, 1, TYPE> fc1_;
    Core::LinAlg::Matrix<3, 1, TYPE> fc2_;

    //! coordinates of contact point on center lines of beams
    Core::LinAlg::Matrix<3, 1, TYPE> x1_;
    Core::LinAlg::Matrix<3, 1, TYPE> x2_;

    //! parameter value of contact point on beam element
    TYPE xicontact_;

    //! normal vector of current time step
    Core::LinAlg::Matrix<3, 1, TYPE> normal_;

    //@}
  };
}  // namespace BeamInteraction

FOUR_C_NAMESPACE_CLOSE

#endif
