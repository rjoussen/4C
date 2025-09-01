// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_CUT_VOLUME_INTEGRATION_HPP
#define FOUR_C_CUT_VOLUME_INTEGRATION_HPP

#include "4C_config.hpp"

#include "4C_cut_enum.hpp"
#include "4C_cut_facet_integration.hpp"
#include "4C_cut_volumecell.hpp"
#include "4C_linalg_serialdensematrix.hpp"
#include "4C_linalg_serialdensevector.hpp"
#include "4C_linalg_utils_sparse_algebra_math.hpp"

#include <vector>

FOUR_C_NAMESPACE_OPEN


namespace Cut
{
  class Element;
  class Facet;
  class VolumeCell;
  // While performing the volume integration, the points of the facet should be arranged in
  // anti-clockwise manner when looking the surface away from the body this ensures outward normal
  // vector
  /*!
  \brief Class to construct integration rule over the volumecell in element local coordinates.
  It performs integration of base functions over the volume, distribution of Gauss points and
  solution of moment fitting matrix to arrive at the Gauss weights
   */
  class VolumeIntegration
  {
   public:
    VolumeIntegration(VolumeCell* volcell, Element* elem, const Cut::Point::PointPosition posi)
        : volcell_(volcell), elem1_(elem), position_(posi)
    {
    }

    /*!
    \brief Returns the location of Gauss points distributed over the volumecell
    */
    std::vector<std::vector<double>> get_gauss_point_location() { return gauss_pts_; }

    /*!
    \brief Check whether the point with this element Local coordinates is inside, outside or
    boundary of this volumecell. The output std::string will be either "outside", "inside" or
    "onBoundary"
    */
    std::string is_point_inside(Core::LinAlg::Matrix<3, 1>& rst);

   private:
    /*
    \brief Check whether the intersection point, which is in the plane containing the facet,
    actually lies with in the facet area
    */
    int pnpoly(const std::vector<std::vector<double>>& xp, const Core::LinAlg::Matrix<3, 1>& pt,
        Cut::ProjectionDirection projType);

    //! considered volumecell
    VolumeCell* volcell_;

    //! background element that contains this volumecell
    Element* elem1_;

    //! position (inside or outside) of volumecell
    const Cut::Point::PointPosition position_;

    //! position of Gauss points
    std::vector<std::vector<double>> gauss_pts_;

    //! equations of plane in which facets of volumecell are contained
    std::vector<std::vector<double>> eqn_facets_;
  };
}  // namespace Cut


FOUR_C_NAMESPACE_CLOSE

#endif
