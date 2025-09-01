// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_CUT_FACET_INTEGRATION_HPP
#define FOUR_C_CUT_FACET_INTEGRATION_HPP

#include "4C_config.hpp"

#include "4C_cut_element.hpp"
#include "4C_cut_enum.hpp"
#include "4C_cut_mesh.hpp"

FOUR_C_NAMESPACE_OPEN

// #define DIRECTDIV_EXTENDED_DEBUG_OUTPUT
// #define TRIANGULATE_ALL_FACETS_FOR_DIVERGENCECELLS


namespace Cut
{
  class Element;
  class Facet;

  /*!
  \brief This class performs the integration of base functions over the facet. The points of the
  facet should be arranged in anti-clockwise manner when looking the facet away from the volume
  this ensures outward normal vector when divergence theorem is used
  */
  class FacetIntegration
  {
   public:
    FacetIntegration(Facet* face1, Element* element1, const Cut::Point::PointPosition posi,
        bool bcellInt, bool global)
        : face1_(face1),     // facet under consideration
          elem1_(element1),  // the element for which the facet is a part of
          position_(posi)    // position
    {
    }

    /*!
    \brief Computes the equation of the plane that contains this facet
    */
    std::vector<double> equation_plane(const std::vector<std::vector<double>>& cornersLocal);

    /*!
    \brief Generate Gaussian points over the considered facet by triangulating it. This is used
    when DirectDivergence option is used for Gauss point generation
    */
    void divergence_integration_rule_new(Mesh& mesh, Core::FE::CollectedGaussPoints& cgp);

   private:
    /*!
    \brief Check whether the vertices numbering of the facet is clockwise
    */
    void is_clockwise(
        const std::vector<double>& eqn_plane, const std::vector<std::vector<double>>& cornersLocal);

    /*!
    \brief Generate boundary cells for the considered facet. May need to perform triangulatio
    */
    void generate_divergence_cells_new(bool divergenceRule, Mesh& mesh,
        std::list<std::shared_ptr<BoundaryCell>>& divCells,
        const std::vector<Point*>& cornersGlobal);


    /*!
    \brief Temporarily create Tri3 cell. This is not stored in Mesh
    */
    void temporary_tri3(
        const std::vector<Point*>& corners, std::list<std::shared_ptr<BoundaryCell>>& divCells);

    /*!
    \brief Temporarily create Quad4 cell. This is not stored in Mesh
    */
    void temporary_quad4(
        const std::vector<Point*>& corners, std::list<std::shared_ptr<BoundaryCell>>& divCells);

    //! considered facet
    Facet* face1_;

    //! background element which was cut to produce this facet
    Element* elem1_;

    //! position of the facet
    const Cut::Point::PointPosition position_;

    //! True if nodes of facet are arranged to give inward normal
    bool clockwise_;

    //! equation of plane that contains the facet
    std::vector<double> eqn_plane_;

    std::list<std::shared_ptr<BoundaryCell>> boundarycells_;
  };
}  // namespace Cut


FOUR_C_NAMESPACE_CLOSE

#endif
