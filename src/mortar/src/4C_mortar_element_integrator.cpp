// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_fem_general_utils_integration.hpp"
#include "4C_mortar_element.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 |  ctor (public)                                             popp 08/08|
 *----------------------------------------------------------------------*/
Mortar::ElementIntegrator::ElementIntegrator(Core::FE::CellType eletype)
{
  //*********************************************************************
  // Create integration points according to eletype!
  // Note that our standard Gauss rules are:
  // 5  points: for integrals on 1D lines                 (1,2,3,4,5)
  // 7  points: for integrals on 2D 1st triangles         (1,3,6,7,12,16,37)
  // 16 points: for integrals on 2D 1st triangles         (1,3,6,7,12,16,37)
  // 9  points: for integrals on 2D 1st order quadrilaterals
  // 25 points: for integrals on 2D 2nd order quadrilaterals
  //**********************************************************************
  std::shared_ptr<Core::FE::IntegrationPoints2D> rule2d;

  switch (eletype)
  {
    case Core::FE::CellType::line2:
    case Core::FE::CellType::line3:
    case Core::FE::CellType::nurbs2:
    case Core::FE::CellType::nurbs3:
    {
      const Core::FE::IntegrationPoints1D intpoints(Core::FE::GaussRule1D::line_5point);
      ngp_ = intpoints.nquad;
      coords_.reshape(n_gp(), 2);
      weights_.resize(n_gp());
      for (int i = 0; i < n_gp(); ++i)
      {
        coords_(i, 0) = intpoints.qxg[i][0];
        coords_(i, 1) = 0.0;
        weights_[i] = intpoints.qwgt[i];
      }
      break;
    }
    case Core::FE::CellType::tri3:
      rule2d = std::make_shared<Core::FE::IntegrationPoints2D>(Core::FE::GaussRule2D::tri_7point);
      break;
    case Core::FE::CellType::tri6:
      rule2d = std::make_shared<Core::FE::IntegrationPoints2D>(Core::FE::GaussRule2D::tri_16point);
      break;
    case Core::FE::CellType::quad4:
      rule2d = std::make_shared<Core::FE::IntegrationPoints2D>(Core::FE::GaussRule2D::quad_9point);
      break;
    case Core::FE::CellType::quad8:
    case Core::FE::CellType::quad9:
    case Core::FE::CellType::nurbs4:
    case Core::FE::CellType::nurbs9:
      rule2d = std::make_shared<Core::FE::IntegrationPoints2D>(Core::FE::GaussRule2D::quad_25point);
      break;
    default:
      FOUR_C_THROW("ElementIntegrator: This contact element type is not implemented!");
  }  // switch(eletype)

  // save Gauss points for all 2D rules
  if (rule2d != nullptr)
  {
    ngp_ = rule2d->nquad;
    coords_.reshape(n_gp(), 2);
    weights_.resize(n_gp());
    for (int i = 0; i < n_gp(); ++i)
    {
      coords_(i, 0) = rule2d->qxg[i][0];
      coords_(i, 1) = rule2d->qxg[i][1];
      weights_[i] = rule2d->qwgt[i];
    }
  }

  return;
}

FOUR_C_NAMESPACE_CLOSE
