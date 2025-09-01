// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_cut_volume_integration.hpp"

#include "4C_cut_boundingbox.hpp"
#include "4C_cut_enum.hpp"
#include "4C_cut_options.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>

FOUR_C_NAMESPACE_OPEN

/*-------------------------------------------------------------------------------------------------------------------*
???
*--------------------------------------------------------------------------------------------------------------------*/
int Cut::VolumeIntegration::pnpoly(const std::vector<std::vector<double>>& xp,
    const Core::LinAlg::Matrix<3, 1>& pt, Cut::ProjectionDirection projType)
{
  int npol = xp.size();
  int ind1 = 1, ind2 = 2;
  if (projType == Cut::proj_y)
  {
    ind1 = 2;
    ind2 = 0;
  }
  else if (projType == Cut::proj_z)
  {
    ind1 = 0;
    ind2 = 1;
  }

  // check whether given point is one of the corner points
  for (int i = 0; i < npol; i++)
  {
    if (fabs(xp[i][ind1] - pt(ind1, 0)) < 1e-8 && fabs(xp[i][ind2] - pt(ind2, 0)) < 1e-8) return 1;
  }

  // check whether given point is on the edge
  for (int i = 0; i < npol; i++)
  {
    std::vector<double> end1 = xp[i];
    std::vector<double> end2 = xp[(i + 1) % npol];

    // line with constant ind1
    if (fabs(end1[ind1] - end2[ind1]) < 1e-8 &&
        fabs(end1[ind1] - pt(ind1, 0)) < 1e-8)  // pt is on infinite line check
    {
      // pt is within limits check
      if (pt(ind2, 0) > std::min(end1[ind2], end2[ind2]) &&
          pt(ind2, 0) < std::max(end1[ind2], end2[ind2]))
        return 1;
      /*else
        return 0;*/
    }
    // line with constant ind2
    if (fabs(end1[ind2] - end2[ind2]) < 1e-8 && fabs(end1[ind2] - pt(ind2, 0)) < 1e-8)
    {
      if (pt(ind1, 0) > std::min(end1[ind1], end2[ind1]) &&
          pt(ind1, 0) < std::max(end1[ind1], end2[ind1]))
        return 1;
      /*else
        return 0;*/
    }

    // pt fall on edge if it satisfies eqn of edge
    if (pt(ind2, 0) > std::min(end1[ind2], end2[ind2]) &&  // check within limits
        pt(ind2, 0) < std::max(end1[ind2], end2[ind2]))
    {
      double xpre = end1[ind1] + (end2[ind1] - end1[ind1]) / (end2[ind2] - end1[ind2]) *
                                     (pt(ind2, 0) - end1[ind2]);
      if (fabs(xpre - pt(ind1, 0)) < 1e-8)  // check eqn of line is satisfied with given pt
        return 1;
    }
  }

  // check for inside
  int i, j, c = 0;
  for (i = 0, j = npol - 1; i < npol; j = i++)
  {
    if ((((xp[i][ind2] <= pt(ind2, 0)) && (pt(ind2, 0) < xp[j][ind2])) ||
            ((xp[j][ind2] <= pt(ind2, 0)) && (pt(ind2, 0) < xp[i][ind2]))) &&
        (pt(ind1, 0) < (xp[j][ind1] - xp[i][ind1]) * (pt(ind2, 0) - xp[i][ind2]) /
                               (xp[j][ind2] - xp[i][ind2]) +
                           xp[i][ind1]))
      c = !c;
  }
  return c;
}

/*-------------------------------------------------------------------------------------*
 * Check whether the point with this element Local coordinates is inside,              *
 * outside or on boundary of this volumecell                            sudhakar 07/12 *
 *-------------------------------------------------------------------------------------*/
std::string Cut::VolumeIntegration::is_point_inside(Core::LinAlg::Matrix<3, 1>& rst)
{
  const plain_facet_set& facete = volcell_->facets();

  //--------------------------------------------------------------------------------
  // Step 1: Classify facets into facets having zero normal in x-direction and other
  //--------------------------------------------------------------------------------
  std::vector<plain_facet_set::const_iterator> XFacets;     // facets with non-zero nx
  std::vector<plain_facet_set::const_iterator> NotXFacets;  // facets with zero nx
  std::vector<std::vector<double>> Eqnplane;                // eqn of plane for all facets

  for (plain_facet_set::const_iterator i = facete.begin(); i != facete.end(); i++)
  {
    Facet* fe = *i;
    std::vector<std::vector<double>> cornersLocal;
    fe->corner_points_local(elem1_, cornersLocal);

    FacetIntegration faee1(fe, elem1_, position_, false, false);
    std::vector<double> eqnFacet = faee1.equation_plane(cornersLocal);
    Eqnplane.push_back(eqnFacet);

    if (fabs(eqnFacet[0]) > 1e-10)
      XFacets.push_back(i);
    else
      NotXFacets.push_back(i);
  }

  //-------------------------------------------------------------------------
  // Step 2: Shoot a ray along x-direction and find all intersection points
  //-------------------------------------------------------------------------
  std::map<double, int> Xintersect;
  // ray intersects only XFacets
  for (unsigned i = 0; i < XFacets.size(); i++)
  {
    // check whether the infinite ray thru given (y,z) intersect the facet
    Facet* fe = *XFacets[i];
    std::vector<std::vector<double>> cornersLocal;
    fe->corner_points_local(elem1_, cornersLocal);
    int cutno = pnpoly(cornersLocal, rst, Cut::proj_x);
    if (cutno == 1)
    {
      // find x-value of intersection point, (yInt,zInt) = (y,z) of given pt
      const int idx = std::distance(facete.begin(), XFacets[i]);
      std::vector<double> eqn = Eqnplane[idx];
      double x = (eqn[3] - eqn[1] * rst(1, 0) - eqn[2] * rst(2, 0)) / eqn[0];
      if (fabs((x - rst(0, 0)) / x) < 1e-8)  // pt is on one of the XFacets
        return "onBoundary";
      Xintersect[x] = i;
    }
  }

  //-------------------------------------------------------------------------------------
  // Step 3: based on relative location of intersection points w.r to given point, decide
  //-------------------------------------------------------------------------------------
  int numInter = Xintersect.size();

  // Check to see if point is in x-z or x-y oriented facets
  for (unsigned i = 0; i < NotXFacets.size(); i++)
  {
    const int idx = std::distance(facete.begin(), NotXFacets[i]);
    std::vector<double> eqn = Eqnplane[idx];
    // check pt on x-y plane facet
    if (fabs(eqn[2]) < 1e-10)  // make sure it is x-y facet
    {
      if (fabs((eqn[3] / eqn[1] - rst(1, 0)) / rst(1, 0)) < 1e-10)  // make sure pt is on same plane
      {
        Facet* fe = *NotXFacets[i];
        std::vector<std::vector<double>> cornersLocal;
        fe->corner_points_local(elem1_, cornersLocal);
        int cutno = pnpoly(cornersLocal, rst, Cut::proj_y);
        if (cutno == 1)  // make sure pt is within facet area
        {
          return "onBoundary";
        }
      }
    }
    // check pt on x-z plane facet
    if (fabs(eqn[1]) < 1e-10)  // make sure it is x-z facet
    {
      if (fabs((eqn[3] / eqn[2] - rst(2, 0)) / rst(2, 0)) < 1e-10)
      {
        Facet* fe = *NotXFacets[i];
        std::vector<std::vector<double>> cornersLocal;
        fe->corner_points_local(elem1_, cornersLocal);
        int cutno = pnpoly(cornersLocal, rst, Cut::proj_z);
        if (cutno == 1)
        {
          return "onBoundary";
        }
      }
    }
  }

  // if pt is not on x-y or x-z facet and no intersection --> outside
  if (numInter == 0) return "outside";

  // add given point --> to sort intersection pts along with given pt
  Xintersect[rst(0, 0)] = -1;
  numInter++;


  std::map<double, int>::iterator it = Xintersect.find(rst(0, 0));

  // all intersecting facets are right side to given pt
  if (it == Xintersect.begin()) return "outside";

  std::map<double, int>::iterator itEnd = Xintersect.end();
  --itEnd;
  // all intersecting facets are left side to given pt
  if (it == itEnd) return "outside";

  // if the no of facets right side of given pt is even number --> outside
  int rightFacets = std::distance(it, itEnd);
  if (rightFacets % 2 == 0) return "outside";
  return "inside";
}

FOUR_C_NAMESPACE_CLOSE
