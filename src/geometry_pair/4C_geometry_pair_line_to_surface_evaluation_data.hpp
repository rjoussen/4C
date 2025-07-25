// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_GEOMETRY_PAIR_LINE_TO_SURFACE_EVALUATION_DATA_HPP
#define FOUR_C_GEOMETRY_PAIR_LINE_TO_SURFACE_EVALUATION_DATA_HPP


#include "4C_config.hpp"

#include "4C_geometry_pair_line_to_3D_evaluation_data.hpp"

#include <memory>
#include <unordered_map>

// Forward declarations.

FOUR_C_NAMESPACE_OPEN

namespace GeometryPair
{
  class FaceElement;
}

namespace GeometryPair
{
  /**
   * \brief Class to manage input parameters and evaluation data for line to surface interactions.
   */
  class LineToSurfaceEvaluationData : public LineTo3DEvaluationData
  {
   public:
    /**
     * \brief Constructor (derived).
     */
    LineToSurfaceEvaluationData(const Teuchos::ParameterList& input_parameter_list);

    /**
     * \brief Reset the evaluation data (derived).
     */
    void clear() override;

    /**
     * \brief Setup the surface data.
     *
     * \param discret (in) Pointer to the discretization.
     * \param face_elements (in) Map to all face elements in this condition on this rank.
     */
    void setup(const std::shared_ptr<const Core::FE::Discretization>& discret,
        const std::unordered_map<int, std::shared_ptr<GeometryPair::FaceElement>>& face_elements);

    /**
     * \brief Calculate the averaged nodal normals.
     */
    void set_state(const std::shared_ptr<const Core::LinAlg::Vector<double>>& displacement_col_np);

    /**
     * \brief Get a reference to the face element map.
     */
    const std::unordered_map<int, std::shared_ptr<GeometryPair::FaceElement>>& get_face_elements()
        const
    {
      return face_elements_;
    }

    /**
     * \brief Return the strategy to be used for the surface normals.
     */
    GeometryPair::SurfaceNormals get_surface_normal_strategy() const
    {
      return surface_normal_strategy_;
    }

   private:
    //! A map of all face elements needed for this surface.
    std::unordered_map<int, std::shared_ptr<GeometryPair::FaceElement>> face_elements_;

    //! Strategy to be used for surface normals.
    GeometryPair::SurfaceNormals surface_normal_strategy_;
  };
}  // namespace GeometryPair

FOUR_C_NAMESPACE_CLOSE

#endif
