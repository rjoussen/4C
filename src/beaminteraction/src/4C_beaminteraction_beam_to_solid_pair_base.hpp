// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_BEAMINTERACTION_BEAM_TO_SOLID_PAIR_BASE_HPP
#define FOUR_C_BEAMINTERACTION_BEAM_TO_SOLID_PAIR_BASE_HPP


#include "4C_config.hpp"

#include "4C_beaminteraction_contact_pair.hpp"
#include "4C_geometry_pair_element.hpp"
#include "4C_geometry_pair_utility_classes.hpp"

FOUR_C_NAMESPACE_OPEN


namespace BeamInteraction
{
  /**
   * \brief Base class for beam to solid interactions.
   * @tparam scalar_type Scalar FAD type to be used in this pair.
   * @tparam segments_scalar_type Scalar FAD type to be used for the beam-to-solid segments.
   * @tparam beam Type from GeometryPair::ElementDiscretization... representing the beam.
   * @tparam solid Type from GeometryPair::ElementDiscretization... representing the solid.
   */
  template <typename ScalarType, typename SegmentsScalarType, typename Beam, typename Solid>
  class BeamToSolidPairBase : public BeamContactPair
  {
   public:
    /**
     * \brief Standard Constructor
     */
    BeamToSolidPairBase();


    /**
     * \brief Setup the contact pair (derived).
     *
     * This method sets the beam reference positions for this pair.
     */
    void setup() override;

    /**
     * \brief Evaluate this contact element pair.
     * @param forcevec1 (out) Force vector on element 1.
     * @param forcevec2 (out) Force vector on element 2.
     * @param stiffmat11 (out) Stiffness contributions on element 1 - element 1.
     * @param stiffmat12 (out) Stiffness contributions on element 1 - element 2.
     * @param stiffmat21 (out) Stiffness contributions on element 2 - element 1.
     * @param stiffmat22 (out) Stiffness contributions on element 2 - element 2.
     * @return True if pair is in contact.
     */
    bool evaluate(Core::LinAlg::SerialDenseVector* forcevec1,
        Core::LinAlg::SerialDenseVector* forcevec2, Core::LinAlg::SerialDenseMatrix* stiffmat11,
        Core::LinAlg::SerialDenseMatrix* stiffmat12, Core::LinAlg::SerialDenseMatrix* stiffmat21,
        Core::LinAlg::SerialDenseMatrix* stiffmat22) override
    {
      return false;
    };

    /**
     * \brief Update state of translational nodal DoFs (absolute positions and tangents) of the beam
     * element.
     * @param beam_centerline_dofvec
     * @param solid_nodal_dofvec
     */
    void reset_state(const std::vector<double>& beam_centerline_dofvec,
        const std::vector<double>& solid_nodal_dofvec) override;

    /**
     * \brief Set the restart displacement in this pair.
     *
     * If coupling interactions should be evaluated w.r.t the restart state, this method will set
     * them in the pair accordingly.
     *
     * @param centerline_restart_vec_ (in) Vector with the centerline displacements at the restart
     * step, for all contained elements (Vector of vector).
     */
    void set_restart_displacement(
        const std::vector<std::vector<double>>& centerline_restart_vec_) override;

    /**
     * \brief Print information about this beam contact element pair to screen.
     */
    void print(std::ostream& out) const override;

    /**
     * \brief Print this beam contact element pair to screen.
     */
    void print_summary_one_line_per_active_segment_pair(std::ostream& out) const override;

    /**
     * \brief Check if this pair is in contact. The correct value is only returned after
     * pre_evaluate and Evaluate are run on the geometry pair.
     * @return true if it is in contact.
     */
    inline bool get_contact_flag() const override
    {
      // The element pair is assumed to be active when we have at least one active contact point
      if (line_to_3D_segments_.size() > 0)
        return true;
      else
        return false;
    };

    /**
     * \brief Get number of active contact point pairs on this element pair. Not yet implemented.
     */
    unsigned int get_num_all_active_contact_point_pairs() const override
    {
      FOUR_C_THROW("get_num_all_active_contact_point_pairs not yet implemented!");
      return 0;
    };

    /**
     * \brief Get coordinates of all active contact points on element1. Not yet implemented.
     */
    void get_all_active_contact_point_coords_element1(
        std::vector<Core::LinAlg::Matrix<3, 1, double>>& coords) const override
    {
      FOUR_C_THROW("get_all_active_contact_point_coords_element1 not yet implemented!");
    }

    /**
     * \brief Get coordinates of all active contact points on element2. Not yet implemented.
     */
    void get_all_active_contact_point_coords_element2(
        std::vector<Core::LinAlg::Matrix<3, 1, double>>& coords) const override
    {
      FOUR_C_THROW("get_all_active_contact_point_coords_element2 not yet implemented!");
    }

    /**
     * \brief Get all (scalar) contact forces of this contact pair. Not yet implemented.
     */
    void get_all_active_beam_to_beam_visualization_values(std::vector<double>& forces,
        std::vector<double>& gaps, std::vector<double>& angles,
        std::vector<int>& types) const override
    {
      FOUR_C_THROW("get_all_active_contact_forces not yet implemented!");
    }

    /**
     * \brief Get energy of penalty contact. Not yet implemented.
     */
    double get_energy() const override
    {
      FOUR_C_THROW("get_energy not implemented yet!");
      return 0.0;
    }

    /**
     * \brief Returns the type of this beam contact element pair.
     */
    ContactPairType get_type() const override { return ContactPairType::beam_to_solid_base; }

   protected:
    /**
     * \brief This function evaluates the beam position at an integration point for the pairs.
     *
     * This is needed because the cross section pairs have 3 parameter coordinates on the beam and
     * the other pairs have 1. This method is mainly used for visualization.
     *
     * @param integration_point (in) Integration where the position should be evaluated.
     * @param r_beam (out) Position on the beam.
     * @param reference (in) True -> the reference position is calculated, False -> the current
     * position is calculated.
     */
    virtual void evaluate_beam_position_double(
        const GeometryPair::ProjectionPoint1DTo3D<double>& integration_point,
        Core::LinAlg::Matrix<3, 1, double>& r_beam, bool reference) const;

   protected:
    //! Vector with the segments of the line to 3D pair.
    std::vector<GeometryPair::LineSegment<SegmentsScalarType>> line_to_3D_segments_;

    //! Current nodal positions (and tangents) of the beam.
    GeometryPair::ElementData<Beam, ScalarType> ele1pos_;

    //! Reference nodal positions (and tangents) of the beam.
    GeometryPair::ElementData<Beam, double> ele1posref_;
  };
}  // namespace BeamInteraction

FOUR_C_NAMESPACE_CLOSE

#endif
