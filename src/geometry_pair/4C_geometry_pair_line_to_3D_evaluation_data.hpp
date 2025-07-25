// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_GEOMETRY_PAIR_LINE_TO_3D_EVALUATION_DATA_HPP
#define FOUR_C_GEOMETRY_PAIR_LINE_TO_3D_EVALUATION_DATA_HPP


#include "4C_config.hpp"

#include "4C_fem_general_utils_integration.hpp"
#include "4C_geometry_pair_evaluation_data_base.hpp"
#include "4C_geometry_pair_input.hpp"
#include "4C_geometry_pair_utility_classes.hpp"

#include <set>

FOUR_C_NAMESPACE_OPEN

namespace GeometryPair
{
  /**
   * \brief Class to manage inout parameters and evaluation data for line to 3D interactions.
   */
  class LineTo3DEvaluationData : public GeometryEvaluationDataBase
  {
   public:
    /**
     * \brief Constructor (derived).
     */
    LineTo3DEvaluationData(const Teuchos::ParameterList& input_parameter_list);


    /**
     * \brief Clear the evaluation data.
     */
    void clear() override;

    /**
     * \brief Reset the evaluation data tracker objects.
     */
    void reset_tracker();

    /**
     * \brief Get the segmentation strategy.
     * @return flag for segmentation strategy.
     */
    inline GeometryPair::LineTo3DStrategy get_strategy() const { return strategy_; }

    /**
     * \brief Get the number of search points for segmentation search.
     * @return number of search points.
     */
    inline unsigned int get_number_of_search_points() const { return n_search_points_; }

    /**
     * \brief Get the flaf on what to do if not all Gauss points of a segment project valid.
     */
    inline GeometryPair::NotAllGaussPointsProjectValidAction
    get_not_all_gauss_points_project_valid_action() const
    {
      return not_all_gauss_points_project_valid_action_;
    }

    /**
     * \brief Get the Gauss rule to be used for Gauss point projection method.
     * @return Gauss rule
     */
    inline Core::FE::GaussRule1D get_gauss_rule() const { return gauss_rule_; };

    /**
     * \brief Get the Gauss rule to be used for Gauss point projection method.
     * @return Gauss rule
     */
    inline Core::FE::IntegrationPoints1D get_gauss_points() const
    {
      return Core::FE::IntegrationPoints1D(gauss_rule_);
    };

    /**
     * \brief Get the number of Gauss points.
     * @return Gauss rule
     */
    inline int get_number_of_gauss_points() const { return get_gauss_points().nquad; };

    /**
     * \brief Returns the number of integration points along the circumference of the line cross
     * section.
     * @return number of integration points in circumferential direction.
     */
    inline unsigned int get_number_of_integration_points_circumference() const
    {
      return integration_points_circumference_;
    }

    /**
     * \brief Return a reference to the gauss point projection tracker.
     * @return Projection tracker.
     */
    inline std::map<int, std::vector<bool>>& get_gauss_point_projection_tracker()
    {
      return gauss_point_projection_tracker_;
    };

    /**
     * \brief Return a reference to the segment tracker.
     * @return Segment tracker.
     */
    std::map<int, std::set<LineSegment<double>>>& get_segment_tracker()
    {
      return segment_tracker_;
    };

   private:
    //! Strategy to be used for contact search.
    GeometryPair::LineTo3DStrategy strategy_;

    //! Gauss rule for Gauss point projection method.
    Core::FE::GaussRule1D gauss_rule_;

    //! Number of integration points in the circumferential direction of the line cross section.
    unsigned int integration_points_circumference_;

    //! Gauss point projection tracking vector.
    std::map<int, std::vector<bool>> gauss_point_projection_tracker_;

    //! Number of points for segmentation search.
    unsigned int n_search_points_;

    //! What to do if not all Gauss points of a segment project valid
    GeometryPair::NotAllGaussPointsProjectValidAction not_all_gauss_points_project_valid_action_;

    //! Segment tracking vector for segmentation. We use double in this case, because otherwise the
    //! class would have to be templated on the type of this tracker.
    std::map<int, std::set<LineSegment<double>>> segment_tracker_;
  };
}  // namespace GeometryPair

FOUR_C_NAMESPACE_CLOSE

#endif
