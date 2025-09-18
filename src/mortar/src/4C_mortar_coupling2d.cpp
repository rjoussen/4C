// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_mortar_coupling2d.hpp"

#include "4C_fem_discretization.hpp"
#include "4C_linalg_serialdensevector.hpp"
#include "4C_linalg_utils_densematrix_inverse.hpp"
#include "4C_linalg_utils_densematrix_multiply.hpp"
#include "4C_mortar_defines.hpp"
#include "4C_mortar_element.hpp"
#include "4C_mortar_integrator.hpp"
#include "4C_mortar_node.hpp"
#include "4C_mortar_projector.hpp"

#include <utility>

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Mortar::Coupling2d::Coupling2d(Core::FE::Discretization& idiscret, int dim, bool quad,
    Teuchos::ParameterList& params, Mortar::Element& sele, Mortar::Element& mele)
    : idiscret_(idiscret),
      dim_(dim),
      quad_(quad),
      imortar_(params),
      sele_(sele),
      mele_(mele),
      overlap_(false)
{
  // initialize variables
  hasproj_.resize(4);
  xiproj_.resize(4);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
MPI_Comm Mortar::Coupling2d::get_comm() const { return idiscret_.get_comm(); }


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool Mortar::Coupling2d::rough_check_orient()
{
  // we first need the master element center
  std::array<double, 2> loccenter = {0.0, 0.0};

  // compute the unit normal vector at the slave element center
  std::array<double, 3> nsc = {0.0, 0.0, 0.0};
  slave_element().compute_unit_normal_at_xi(loccenter.data(), nsc.data());

  // compute the unit normal vector at the master element center
  std::array<double, 3> nmc = {0.0, 0.0, 0.0};
  master_element().compute_unit_normal_at_xi(loccenter.data(), nmc.data());

  // check orientation of the two normals
  double dot = nsc[0] * nmc[0] + nsc[1] * nmc[1] + nsc[2] * nmc[2];
  return (dot < -0.1);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool Mortar::Coupling2d::project()
{
  // initialize projection status
  hasproj_[0] = false;  // slave 0 end node
  hasproj_[1] = false;  // slave 1 end node
  hasproj_[2] = false;  // master 0 end node
  hasproj_[3] = false;  // master 1 end node

  // rough check of orientation of element centers
  // if slave and master element center normals form an
  // angle > 90 degrees the pair will not be considered further
  bool orient = rough_check_orient();
  if (!orient) return false;

  // get slave and master element nodes
  Core::Nodes::Node** mysnodes = slave_element().nodes();
  if (!mysnodes) FOUR_C_THROW("IntegrateOverlap: Null pointer for mysnodes!");
  Core::Nodes::Node** mymnodes = master_element().nodes();
  if (!mymnodes) FOUR_C_THROW("IntegrateOverlap: Null pointer for mymnodes!");

  // project slave nodes onto master element
  for (int i = 0; i < slave_element().num_node(); ++i)
  {
    auto* snode = dynamic_cast<Mortar::Node*>(mysnodes[i]);
    std::array<double, 2> xi = {0.0, 0.0};

    if (slave_element().shape() == Core::FE::CellType::nurbs3)
    {
      std::array<double, 2> xinode = {0., 0.};
      if (i == 0)
      {
        xinode[0] = -1.;
      }
      if (i == 1)
      {
        xinode[0] = +1.;
      }

      // for nurbs we need to use the Gauss point projector, since the actual spatial coords
      // of the point to be projected is calculated by N*X using shape functions N and CP coords X
      Mortar::Projector::impl(slave_element(), mele_)
          ->project_gauss_point_2d(slave_element(), xinode.data(), mele_, xi.data());
    }
    else
    {
      // TODO random?
      Mortar::Projector::impl(slave_element())
          ->project_nodal_normal(*snode, master_element(), xi.data());
    }

    // save projection if it is feasible
    // we need an expanded feasible domain in order to check pathological
    // cases due to round-off error and iteration tolerances later!
    if ((-1.0 - MORTARPROJTOL <= xi[0]) && (xi[0] <= 1.0 + MORTARPROJTOL))
    {
      // for element overlap only the outer nodes are of interest
      if (i < 2)
      {
        hasproj_[i] = true;
        xiproj_[i] = xi[0];
      }
      // nevertheless we need the inner node projection status later (weighted gap)
      snode->has_proj() = true;
    }
  }

  // project master nodes onto slave element
  for (int i = 0; i < 2; ++i)
  {
    auto* mnode = dynamic_cast<Mortar::Node*>(mymnodes[i]);
    std::array<double, 2> xi = {0.0, 0.0};

    if (master_element().shape() == Core::FE::CellType::nurbs3)
    {
      std::array<double, 2> xinode = {0., 0.};
      if (i == 0)
      {
        xinode[0] = -1.;
      }
      if (i == 1)
      {
        xinode[0] = +1.;
      }

      // for nurbs, we introduce a dummy mortar node at the actual spatial position of the master
      // side element boundary. Hence, we need that location
      std::array<double, 3> xm{};
      Core::LinAlg::SerialDenseVector mval(mele_.num_node());
      Core::LinAlg::SerialDenseMatrix deriv(mele_.num_node(), 1);
      mele_.evaluate_shape(xinode.data(), mval, deriv, mele_.num_node());

      for (int mn = 0; mn < master_element().num_node(); mn++)
      {
        auto* mnode2 = dynamic_cast<Mortar::Node*>(mymnodes[mn]);
        for (int dim = 0; dim < 2; ++dim) xm[dim] += mval(mn) * mnode2->xspatial()[dim];
      }
      std::vector<int> mdofs(2);
      Mortar::Node tmp_node(mnode->id(), xm, mnode->owner(), mdofs, false);
      Mortar::Projector::impl(slave_element())
          ->project_element_normal(tmp_node, slave_element(), xi.data());
    }
    else
    {
      // TODO random?
      Mortar::Projector::impl(slave_element())
          ->project_element_normal(*mnode, slave_element(), xi.data());
    }

    // save projection if it is feasible
    // we need an expanded feasible domain in order to check pathological
    // cases due to round-off error and iteration tolerances later!!!
    if ((-1.0 - MORTARPROJTOL <= xi[0]) && (xi[0] <= 1.0 + MORTARPROJTOL))
    {
      // for element overlap only the outer nodes are of interest
      if (i < 2)
      {
        hasproj_[i + 2] = true;
        xiproj_[i + 2] = xi[0];
      }
    }
  }

  return true;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool Mortar::Coupling2d::detect_overlap()
{
  /**********************************************************************/
  /* OVERLAP CASES                                                      */
  /* Depending on mxi and sxi overlap will be decided!                  */
  /* Even for 3noded elements only the two end nodes matter in 2D!      */
  /* There are several cases how the 2 elements can overlap. Handle all */
  /* of them, including the ones that they don't overlap at all!        */
  /**********************************************************************/

  // For the non-overlapping cases, the possibility of an identical local
  // node numbering direction for both sides is taken into account!!
  // (this can happen, when elements far from each other are projected,
  // which actually should be impossible due to the search radius
  // condition in the potential coupling pair search above!
  // But you never know...)
  // For the overlapping cases, it is a prerequisite that the two local
  // node numbering directions are opposite!!
  // (this is the case, when the elements are sufficiently near each other,
  // which is ensured by only processing nodes that fulfill the
  // search radius condition above!)
  // CAUTION: The bool output variable in this method is a REAL output
  // variable, determining whether there is an overlap or not!
  // initialize local working variables
  bool overlap = false;
  double sxia = 0.0;
  double sxib = 0.0;
  double mxia = 0.0;
  double mxib = 0.0;

  // local working copies of input variables
  bool s0hasproj = hasproj_[0];
  bool s1hasproj = hasproj_[1];
  bool m0hasproj = hasproj_[2];
  bool m1hasproj = hasproj_[3];

  std::vector<double> sprojxi(2);
  sprojxi[0] = xiproj_[0];
  sprojxi[1] = xiproj_[1];

  std::vector<double> mprojxi(2);
  mprojxi[0] = xiproj_[2];
  mprojxi[1] = xiproj_[3];

  /* CASE 1 (NO OVERLAP):
   no feasible projection found for any of the 4 outer element nodes  */

  if (!s0hasproj && !s1hasproj && !m0hasproj && !m1hasproj)
  {
    // do nothing
  }

  /* CASES 2-5 (NO OVERLAP):
   feasible projection found only for 1 of the 4 outer element nodes
   (this can happen due to the necessary projection tolerance!!!)     */

  else if (s0hasproj && !s1hasproj && !m0hasproj && !m1hasproj)
  {
    if ((-1.0 + MORTARPROJTOL <= sprojxi[0]) && (sprojxi[0] <= 1.0 - MORTARPROJTOL))
    {
      std::cout << "SElement Node IDs: " << (slave_element().nodes()[0])->id() << " "
                << (slave_element().nodes()[1])->id() << '\n';
      std::cout << "MElement Node IDs: " << (master_element().nodes()[0])->id() << " "
                << (master_element().nodes()[1])->id() << '\n';
      std::cout << "SPROJXI_0: " << sprojxi[0] << " SPROJXI_1: " << sprojxi[1] << '\n';
      std::cout << "MPROJXI_0: " << mprojxi[0] << " MPROJXI_1: " << mprojxi[1] << '\n';
      FOUR_C_THROW("IntegrateOverlap: Significant overlap ignored S{} M{}!", slave_element().id(),
          master_element().id());
    }
  }

  else if (!s0hasproj && s1hasproj && !m0hasproj && !m1hasproj)
  {
    if ((-1.0 + MORTARPROJTOL <= sprojxi[1]) && (sprojxi[1] <= 1.0 - MORTARPROJTOL))
    {
      std::cout << "SElement Node IDs: " << (slave_element().nodes()[0])->id() << " "
                << (slave_element().nodes()[1])->id() << '\n';
      std::cout << "MElement Node IDs: " << (master_element().nodes()[0])->id() << " "
                << (master_element().nodes()[1])->id() << '\n';
      std::cout << "SPROJXI_0: " << sprojxi[0] << " SPROJXI_1: " << sprojxi[1] << '\n';
      std::cout << "MPROJXI_0: " << mprojxi[0] << " MPROJXI_1: " << mprojxi[1] << '\n';
      FOUR_C_THROW("IntegrateOverlap: Significant overlap ignored S{} M{}!", slave_element().id(),
          master_element().id());
    }
  }

  else if (!s0hasproj && !s1hasproj && m0hasproj && !m1hasproj)
  {
    if ((-1.0 + MORTARPROJTOL <= mprojxi[0]) && (mprojxi[0] <= 1.0 - MORTARPROJTOL))
    {
      std::cout << "SElement Node IDs: " << (slave_element().nodes()[0])->id() << " "
                << (slave_element().nodes()[1])->id() << '\n';
      std::cout << "MElement Node IDs: " << (master_element().nodes()[0])->id() << " "
                << (master_element().nodes()[1])->id() << '\n';
      std::cout << "SPROJXI_0: " << sprojxi[0] << " SPROJXI_1: " << sprojxi[1] << '\n';
      std::cout << "MPROJXI_0: " << mprojxi[0] << " MPROJXI_1: " << mprojxi[1] << '\n';
      FOUR_C_THROW("IntegrateOverlap: Significant overlap ignored S{} M{}!", slave_element().id(),
          master_element().id());
    }
  }

  else if (!s0hasproj && !s1hasproj && !m0hasproj && m1hasproj)
  {
    if ((-1.0 + MORTARPROJTOL <= mprojxi[1]) && (mprojxi[1] <= 1.0 - MORTARPROJTOL))
    {
      std::cout << "SElement Node IDs: " << (slave_element().nodes()[0])->id() << " "
                << (slave_element().nodes()[1])->id() << '\n';
      std::cout << "MElement Node IDs: " << (master_element().nodes()[0])->id() << " "
                << (master_element().nodes()[1])->id() << '\n';
      std::cout << "SPROJXI_0: " << sprojxi[0] << " SPROJXI_1: " << sprojxi[1] << '\n';
      std::cout << "MPROJXI_0: " << mprojxi[0] << " MPROJXI_1: " << mprojxi[1] << '\n';
      FOUR_C_THROW("IntegrateOverlap: Significant overlap ignored S{} M{}!", slave_element().id(),
          master_element().id());
    }
  }

  /* CASE 6 (OVERLAP):
   feasible projection found for all 4 outer element nodes
   (this can happen due to the necessary projection tolerance!!!)     */

  else if (s0hasproj && s1hasproj && m0hasproj && m1hasproj)
  {
    overlap = true;

    // switch, since nurbs might not be ordered anti-clockwise!!!
    if (slave_element().normal_fac() * master_element().normal_fac() > 0.)
    {
      // internal case 1 for global CASE 6
      // (equivalent to global CASE 7, slave fully projects onto master)
      if ((sprojxi[0] < 1.0) && (sprojxi[1] > -1.0))
      {
        sxia = -1.0;
        sxib = 1.0;
        mxia = sprojxi[1];  // local node numbering always anti-clockwise!!!
        mxib = sprojxi[0];
        // std::cout << "Problem solved with internal case 1!" << std::endl;
      }

      // internal case 2 for global CASE 6
      // (equivalent to global CASE 8, master fully projects onto slave)
      else if ((mprojxi[0] < 1.0) && (mprojxi[1] > -1.0))
      {
        mxia = -1.0;
        mxib = 1.0;
        sxia = mprojxi[1];  // local node numbering always anti-clockwise!!!
        sxib = mprojxi[0];
        // std::cout << "Problem solved with internal case 2!" << std::endl;
      }

      // internal case 3 for global CASE 6
      // (equivalent to global CASE 9, both nodes no. 0 project successfully)
      else if ((sprojxi[0] < 1.0 + MORTARPROJLIM) && (mprojxi[0] < 1.0 + MORTARPROJLIM))
      {
        sxia = -1.0;
        sxib = mprojxi[0];  // local node numbering always anti-clockwise!!!
        mxia = -1.0;
        mxib = sprojxi[0];
        // std::cout << "Problem solved with internal case 3!" << std::endl;
      }

      // internal case 4 for global CASE 6
      // (equivalent to global CASE 10, both nodes no. 1 project successfully)
      else if ((sprojxi[1] > -1.0 - MORTARPROJLIM) && (mprojxi[1] > -1.0 - MORTARPROJLIM))
      {
        sxia = mprojxi[1];
        sxib = 1.0;  // local node numbering always anti-clockwise!!!
        mxia = sprojxi[1];
        mxib = 1.0;
        // std::cout << "Problem solved with internal case 4!" << std::endl;
      }

      // unknown internal case for global CASE 6
      else
      {
        std::cout << "Mortar::Coupling2d::DetectOverlap " << '\n'
                  << "has detected '4 projections'-case for Sl./Ma. pair " << slave_element().id()
                  << "/" << master_element().id() << '\n';
        std::cout << "SElement Node IDs: " << (slave_element().nodes()[0])->id() << " "
                  << (slave_element().nodes()[1])->id() << '\n';
        std::cout << "MElement Node IDs: " << (master_element().nodes()[0])->id() << " "
                  << (master_element().nodes()[1])->id() << '\n';
        std::cout << "SPROJXI_0: " << sprojxi[0] << " SPROJXI_1: " << sprojxi[1] << '\n';
        std::cout << "MPROJXI_0: " << mprojxi[0] << " MPROJXI_1: " << mprojxi[1] << '\n';
        FOUR_C_THROW("DetectOverlap: Unknown overlap case found in global case 6!");
      }
    }

    else
    {
      // fully projecting slave element: equivalent to global case 7
      if ((sprojxi[0] > -1.) && (sprojxi[1] < 1.))
      {
        sxia = -1.0;
        sxib = 1.0;
        mxia = sprojxi[0];
        mxib = sprojxi[1];
      }
      // fully projecting master element: equivalent to global case 8
      else if ((mprojxi[0] > -1.) && (mprojxi[1] < 1.))
      {
        mxia = -1.0;
        mxib = 1.0;
        sxia = mprojxi[0];
        sxib = mprojxi[1];
      }
      // equivalent to global case 15
      else if ((sprojxi[1] < 1.0 + MORTARPROJLIM) && (mprojxi[0] < 1.0 + MORTARPROJLIM))
      {
        sxia = mprojxi[0];
        sxib = 1.;
        mxia = -1.;
        mxib = sprojxi[1];
      }
      // equivalent to global case 16
      else if ((sprojxi[0] > -1.0 - MORTARPROJLIM) && (mprojxi[1] > -1.0 - MORTARPROJLIM))
      {
        sxia = -1.;
        sxib = mprojxi[1];
        mxia = sprojxi[0];
        mxib = 1.;
      }
      else
      {
        std::cout << "Mortar::Coupling2d::DetectOverlap " << '\n'
                  << "has detected '4 projections'-case for Sl./Ma. pair " << slave_element().id()
                  << "/" << master_element().id() << '\n';
        std::cout << "SElement Node IDs: " << (slave_element().nodes()[0])->id() << " "
                  << (slave_element().nodes()[1])->id() << '\n';
        std::cout << "MElement Node IDs: " << (master_element().nodes()[0])->id() << " "
                  << (master_element().nodes()[1])->id() << '\n';
        std::cout << "SPROJXI_0: " << sprojxi[0] << " SPROJXI_1: " << sprojxi[1] << '\n';
        std::cout << "MPROJXI_0: " << mprojxi[0] << " MPROJXI_1: " << mprojxi[1] << '\n';
        FOUR_C_THROW("DetectOverlap: Unknown overlap case found in global case 6!");
      }
    }
  }

  /* CASES 7-8 (OVERLAP):
   feasible projections found for both nodes of one element, this
   means one of the two elements is projecting fully onto the other!  */

  else if (s0hasproj && s1hasproj && !m0hasproj && !m1hasproj)
  {
    overlap = true;
    sxia = -1.0;
    sxib = 1.0;
    // nurbs may not be numbered anti-clockwise
    if (slave_element().normal_fac() * master_element().normal_fac() > 0.)
    {
      mxia = sprojxi[1];  // local node numbering always anti-clockwise!!!
      mxib = sprojxi[0];
    }
    else
    {
      mxia = sprojxi[0];
      mxib = sprojxi[1];
    }
  }

  else if (!s0hasproj && !s1hasproj && m0hasproj && m1hasproj)
  {
    overlap = true;
    mxia = -1.0;
    mxib = 1.0;
    // nurbs may not be numbered anti-clockwise
    if (slave_element().normal_fac() * master_element().normal_fac() > 0.)
    {
      sxia = mprojxi[1];  // local node numbering always anti-clockwise!!!
      sxib = mprojxi[0];
    }
    else
    {
      sxia = mprojxi[0];
      sxib = mprojxi[1];
    }
  }

  /* CASES 9-10 (OVERLAP):
   feasible projections found for one node of each element, due to
   node numbering only identical local node ID pairs possible!        */

  else if (s0hasproj && !s1hasproj && m0hasproj && !m1hasproj)
  {
    // do the two elements really have an overlap?
    if ((sprojxi[0] > -1.0 + MORTARPROJLIM) && (mprojxi[0] > -1.0 + MORTARPROJLIM))
    {
      overlap = true;
      sxia = -1.0;
      sxib = mprojxi[0];  // local node numbering always anti-clockwise!!!
      mxia = -1.0;
      mxib = sprojxi[0];
    }
    if (slave_element().normal_fac() * master_element().normal_fac() < 0.)
    {
      std::cout << "SElement: " << slave_element().node_ids()[0] << " "
                << slave_element().node_ids()[1] << '\n';
      std::cout << "MElement: " << master_element().node_ids()[0] << " "
                << master_element().node_ids()[1] << '\n';
      std::cout << "s0: " << s0hasproj << " s1: " << s1hasproj << '\n';
      std::cout << "m0: " << m0hasproj << " m1: " << m1hasproj << '\n';
      FOUR_C_THROW("integrate_overlap: Unknown overlap case found!");
    }
  }

  else if (!s0hasproj && s1hasproj && !m0hasproj && m1hasproj)
  {
    // do the two elements really have an overlap?
    if ((sprojxi[1] < 1.0 - MORTARPROJLIM) && (mprojxi[1] < 1.0 - MORTARPROJLIM))
    {
      overlap = true;
      sxia = mprojxi[1];
      sxib = 1.0;  // local node numbering always anti-clockwise!!!
      mxia = sprojxi[1];
      mxib = 1.0;
    }
    if (slave_element().normal_fac() * master_element().normal_fac() < 0.)
    {
      std::cout << "SElement: " << slave_element().node_ids()[0] << " "
                << slave_element().node_ids()[1] << '\n';
      std::cout << "MElement: " << master_element().node_ids()[0] << " "
                << master_element().node_ids()[1] << '\n';
      std::cout << "s0: " << s0hasproj << " s1: " << s1hasproj << '\n';
      std::cout << "m0: " << m0hasproj << " m1: " << m1hasproj << '\n';
      FOUR_C_THROW("integrate_overlap: Unknown overlap case found!");
    }
  }

  /* CASES 11-14 (OVERLAP):
   feasible projections found for 3 out of the total 4 nodes,
   this can either lead to cases 7/8 or 9/10!                         */
  else if (s0hasproj && s1hasproj && m0hasproj && !m1hasproj)
  {
    overlap = true;

    // switch, since nurbs might not be ordered anti-clockwise!!!
    if (slave_element().normal_fac() * master_element().normal_fac() > 0.)
    {
      // equivalent to global case 7
      if (mprojxi[0] > 1.0)
      {
        sxia = -1.0;
        sxib = 1.0;
        mxia = sprojxi[1];  // local node numbering always anti-clockwise!!!
        mxib = sprojxi[0];
      }
      // equivalent to global case 9
      else
      {
        sxia = -1.0;
        sxib = mprojxi[0];  // local node numbering always anti-clockwise!!!
        mxia = -1.0;
        mxib = sprojxi[0];
      }
    }
    else
    {
      // equivalent to global case 7
      if (mprojxi[0] < -1.)
      {
        sxia = -1.0;
        sxib = 1.0;
        mxia = sprojxi[0];
        mxib = sprojxi[1];
      }
      // equivalent to global case 15
      else
      {
        sxia = mprojxi[0];
        sxib = 1.;
        mxia = -1.;
        mxib = sprojxi[1];
      }
    }
  }

  else if (s0hasproj && s1hasproj && !m0hasproj && m1hasproj)
  {
    overlap = true;

    // switch, since nurbs might not be ordered anti-clockwise!!!
    if (slave_element().normal_fac() * master_element().normal_fac() > 0.)
    {
      // equivalent to global case 7
      if (mprojxi[1] < -1.0)
      {
        sxia = -1.0;
        sxib = 1.0;
        mxia = sprojxi[1];  // local node numbering always anti-clockwise!!!
        mxib = sprojxi[0];
      }
      // equivalent to global case 10
      else
      {
        sxia = mprojxi[1];
        sxib = 1.0;  // local node numbering always anti-clockwise!!!
        mxia = sprojxi[1];
        mxib = 1.0;
      }
    }

    else
    {
      // equivalent to global case 7
      if (mprojxi[1] > 1.)
      {
        sxia = -1.;
        sxib = 1.;
        mxia = sprojxi[0];
        mxib = sprojxi[1];
      }
      // equivalent to global case 16
      else
      {
        sxia = -1.;
        sxib = mprojxi[1];
        mxia = sprojxi[0];
        mxib = 1.;
      }
    }
  }

  else if (s0hasproj && !s1hasproj && m0hasproj && m1hasproj)
  {
    overlap = true;

    // switch, since nurbs might not be ordered anti-clockwise!!!
    if (slave_element().normal_fac() * master_element().normal_fac() > 0.)
    {
      // equivalent to global case 8
      if (sprojxi[0] > 1.0)
      {
        mxia = -1.0;
        mxib = 1.0;
        sxia = mprojxi[1];  // local node numbering always anti-clockwise!!!
        sxib = mprojxi[0];
      }
      // equivalent to global case 9
      else
      {
        sxia = -1.0;
        sxib = mprojxi[0];  // local node numbering always anti-clockwise!!!
        mxia = -1.0;
        mxib = sprojxi[0];
      }
    }
    else
    {
      // equivalent to global case 8
      if (sprojxi[0] < -1.)
      {
        mxia = -1.0;
        mxib = 1.0;
        sxia = mprojxi[0];
        sxib = mprojxi[1];
      }
      // equivalent to global case 16
      else
      {
        sxia = -1.;
        sxib = mprojxi[1];
        mxia = sprojxi[0];
        mxib = 1.;
      }
    }
  }

  else if (!s0hasproj && s1hasproj && m0hasproj && m1hasproj)
  {
    overlap = true;

    // switch, since nurbs might not be ordered anti-clockwise!!!
    if (slave_element().normal_fac() * master_element().normal_fac() > 0.)
    {
      // equivalent to global case 8
      if (sprojxi[1] < -1.0)
      {
        mxia = -1.0;
        mxib = 1.0;
        sxia = mprojxi[1];  // local node numbering always anti-clockwise!!!
        sxib = mprojxi[0];
      }
      // equivalent to global case 10
      else
      {
        sxia = mprojxi[1];
        sxib = 1.0;  // local node numbering always anti-clockwise!!!
        mxia = sprojxi[1];
        mxib = 1.0;
      }
    }
    else
    {
      // equivalent to global case 8
      if (sprojxi[1] > 1.)
      {
        mxia = -1.;
        mxib = 1.;
        sxia = mprojxi[0];
        sxib = mprojxi[1];
      }
      // equivalent to global case 15
      else
      {
        sxia = mprojxi[0];
        sxib = 1.;
        mxia = -1.;
        mxib = sprojxi[1];
      }
    }
  }

  /* CASES 15-16 (OVERLAP):
   feasible projections found for one node of each element, due to
   node numbering opposite local node ID pairs possible only for nurbs! */

  else if (!s0hasproj && s1hasproj && m0hasproj && !m1hasproj)
  {
    // only possible, if slave and master side have opposite normal-fac
    if (slave_element().normal_fac() * master_element().normal_fac() > 0.)
    {
      std::cout << "SElement: " << slave_element().node_ids()[0] << " "
                << slave_element().node_ids()[1] << '\n';
      std::cout << "MElement: " << master_element().node_ids()[0] << " "
                << master_element().node_ids()[1] << '\n';
      std::cout << "s0: " << s0hasproj << " s1: " << s1hasproj << '\n';
      std::cout << "m0: " << m0hasproj << " m1: " << m1hasproj << '\n';
      FOUR_C_THROW("integrate_overlap: Unknown overlap case found!");
    }
    if (sprojxi[0] < -1. || mprojxi[0] > 1.)
      overlap = false;
    else
    {
      overlap = true;
      sxia = mprojxi[0];
      sxib = 1.;
      mxia = -1.;
      mxib = sprojxi[1];
    }
  }

  else if (s0hasproj && !s1hasproj && !m0hasproj && m1hasproj)
  {
    // only possible, if slave and master side have opposite normal-fac
    if (slave_element().normal_fac() * master_element().normal_fac() > 0.)
    {
      std::cout << "SElement: " << slave_element().node_ids()[0] << " "
                << slave_element().node_ids()[1] << '\n';
      std::cout << "MElement: " << master_element().node_ids()[0] << " "
                << master_element().node_ids()[1] << '\n';
      std::cout << "s0: " << s0hasproj << " s1: " << s1hasproj << '\n';
      std::cout << "m0: " << m0hasproj << " m1: " << m1hasproj << '\n';
      FOUR_C_THROW("integrate_overlap: Unknown overlap case found!");
    }
    if (sprojxi[0] > 1.)
      overlap = false;
    else
    {
      overlap = true;
      sxia = -1.;
      sxib = mprojxi[1];
      mxia = sprojxi[0];
      mxib = 1.;
    }
  }

  /* CASE DEFAULT: unknown overlap case                                  */
  else
  {
    std::cout << "SElement: " << slave_element().node_ids()[0] << " "
              << slave_element().node_ids()[1] << '\n';
    std::cout << "MElement: " << master_element().node_ids()[0] << " "
              << master_element().node_ids()[1] << '\n';
    std::cout << "s0: " << s0hasproj << " s1: " << s1hasproj << '\n';
    std::cout << "m0: " << m0hasproj << " m1: " << m1hasproj << '\n';
    FOUR_C_THROW("integrate_overlap: Unknown overlap case found!");
  }

  // check for 1:1 node projections and for infeasible limits
  if ((sxia < -1.0) || (sxib > 1.0) || (mxia < -1.0) || (mxib > 1.0))
  {
    if (abs(sxia + 1.0) < MORTARPROJLIM) sxia = -1.0;
    if (abs(sxib - 1.0) < MORTARPROJLIM) sxib = 1.0;
    if (abs(mxia + 1.0) < MORTARPROJLIM) mxia = -1.0;
    if (abs(mxib - 1.0) < MORTARPROJLIM) mxib = 1.0;

    if ((sxia < -1.0) || (sxib > 1.0) || (mxia < -1.0) || (mxib > 1.0))
    {
      //      std::cout << "Slave: " << sxia << " " << sxib << std::endl;
      //      std::cout << "Master: " << mxia << " " << mxib << std::endl;
      //      FOUR_C_THROW("integrate_overlap: Determined infeasible limits!");
      std::cout << "WARNING: integrate_overlap: Determined infeasible limits!" << '\n';
      overlap = false;
    }
  }

  // update integration limits in xiproj_
  xiproj_[0] = sxia;
  xiproj_[1] = sxib;
  xiproj_[2] = mxia;
  xiproj_[3] = mxib;

  // store overlap information
  overlap_ = overlap;

  return true;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool Mortar::Coupling2d::integrate_overlap(
    const std::shared_ptr<Mortar::ParamsInterface>& mparams_ptr)
{
  // explicitly defined shape function type needed
  if (shape_fcn() == Inpar::Mortar::shape_undefined)
    FOUR_C_THROW("IntegrateOverlap called without specific shape function defined!");

  /**********************************************************************/
  /* INTEGRATION                                                        */
  /* Depending on overlap and the xiproj_ entries integrate the Mortar  */
  /* matrices D and M on the overlap of the current sl / ma pair.       */
  /**********************************************************************/

  // no integration if no overlap
  if (!overlap_) return false;

  // set segmentation status of all slave nodes
  // (hassegment_ of a slave node is true if ANY segment/cell
  // is integrated that contributes to this slave node)
  int nnodes = slave_element().num_node();
  Core::Nodes::Node** mynodes = slave_element().nodes();
  if (!mynodes) FOUR_C_THROW("Null pointer!");
  for (int k = 0; k < nnodes; ++k)
  {
    auto* mycnode = dynamic_cast<Mortar::Node*>(mynodes[k]);
    if (!mycnode) FOUR_C_THROW("Null pointer!");
    mycnode->has_segment() = true;
  }

  // local working copies of input variables
  double sxia = xiproj_[0];
  double sxib = xiproj_[1];
  double mxia = xiproj_[2];
  double mxib = xiproj_[3];

  // *******************************************************************
  // different options for mortar integration
  // *******************************************************************
  // (1) no quadratic element(s) involved -> linear LM interpolation
  // (2) quadratic element(s) involved -> quadratic LM interpolation
  // (3) quadratic element(s) involved -> linear LM interpolation
  // (4) quadratic element(s) involved -> piecew. linear LM interpolation
  // *******************************************************************
  Inpar::Mortar::LagMultQuad lmtype = lag_mult_quad();

  // *******************************************************************
  // cases (1), (2) and (3)
  // *******************************************************************
  if (!quad() || (quad() && lmtype == Inpar::Mortar::lagmult_quad) ||
      (quad() && lmtype == Inpar::Mortar::lagmult_lin) ||
      (quad() && lmtype == Inpar::Mortar::lagmult_const))
  {
    // do the overlap integration (integrate and linearize both M and gap)
    Mortar::Integrator::impl(slave_element(), master_element(), interface_params())
        ->integrate_segment_2d(
            slave_element(), sxia, sxib, master_element(), mxia, mxib, get_comm());
  }

  // *******************************************************************
  // case (4)
  // *******************************************************************
  else if (quad() && lmtype == Inpar::Mortar::lagmult_pwlin)
  {
    FOUR_C_THROW("Piecewise linear LM interpolation not (yet?) implemented in 2D");
  }

  // *******************************************************************
  // undefined case
  // *******************************************************************
  else if (quad() && lmtype == Inpar::Mortar::lagmult_undefined)
  {
    FOUR_C_THROW(
        "Lagrange multiplier interpolation for quadratic elements undefined\n"
        "If you are using 2nd order mortar elements, you need to specify LM_QUAD in MORTAR "
        "COUPLING section");
  }

  // *******************************************************************
  // other cases
  // *******************************************************************
  else
  {
    FOUR_C_THROW("integrate_overlap: Invalid case for 2D mortar coupling LM interpolation");
  }

  return true;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Mortar::Coupling2dManager::Coupling2dManager(Core::FE::Discretization& idiscret, int dim, bool quad,
    Teuchos::ParameterList& params, Mortar::Element* sele, std::vector<Mortar::Element*> mele)
    : idiscret_(idiscret),
      dim_(dim),
      quad_(quad),
      imortar_(params),
      sele_(sele),
      mele_(std::move(mele))
{
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Mortar::Coupling2dManager::integrate_coupling(
    const std::shared_ptr<Mortar::ParamsInterface>& mparams_ptr)
{
  // decide which type of numerical integration scheme

  //**********************************************************************
  // STANDARD INTEGRATION (SEGMENTS)
  //**********************************************************************
  if (int_type() == Inpar::Mortar::inttype_segments)
  {
    // loop over all master elements associated with this slave element
    for (int m = 0; m < (int)master_elements().size(); ++m)
    {
      // create Coupling2d object and push back
      coupling().push_back(std::make_shared<Coupling2d>(
          idiscret_, dim_, quad_, imortar_, slave_element(), master_element(m)));

      // project the element pair
      coupling()[m]->project();

      // check for element overlap
      coupling()[m]->detect_overlap();
    }

    // special treatment of boundary elements
    // calculate consistent dual shape functions for this element
    consistent_dual_shape();

    // do mortar integration
    for (int m = 0; m < (int)master_elements().size(); ++m)
      coupling()[m]->integrate_overlap(mparams_ptr);
  }

  //**********************************************************************
  // FAST INTEGRATION (ELEMENTS)
  //**********************************************************************
  else if (int_type() == Inpar::Mortar::inttype_elements ||
           int_type() == Inpar::Mortar::inttype_elements_BS)
  {
    if ((int)master_elements().size() == 0) return;

    // bool for boundary segmentation
    bool boundary_ele = false;

    // *******************************************************************
    // different options for mortar integration
    // *******************************************************************
    // (1) no quadratic element(s) involved -> linear LM interpolation
    // (2) quadratic element(s) involved -> quadratic LM interpolation
    // (3) quadratic element(s) involved -> linear LM interpolation
    // (4) quadratic element(s) involved -> piecew. linear LM interpolation
    // *******************************************************************
    Inpar::Mortar::LagMultQuad lmtype = lag_mult_quad();

    // *******************************************************************
    // cases (1), (2) and (3)
    // *******************************************************************
    if (!quad() || (quad() && lmtype == Inpar::Mortar::lagmult_quad) ||
        (quad() && lmtype == Inpar::Mortar::lagmult_lin))
    {
      Mortar::Integrator::impl(slave_element(), master_element(0), imortar_)
          ->integrate_ele_based_2d(
              slave_element(), master_elements(), &boundary_ele, idiscret_.get_comm());

      // Perform Boundary Segmentation if required
      if (int_type() == Inpar::Mortar::inttype_elements_BS)
      {
        if (boundary_ele)
        {
          // std::cout << "Boundary segmentation for element: " << SlaveElement().Id() << "\n" ;
          // switch, if consistent boundary modification chosen
          if (Teuchos::getIntegralValue<Inpar::Mortar::ConsistentDualType>(
                  imortar_, "LM_DUAL_CONSISTENT") != Inpar::Mortar::consistent_none &&
              shape_fcn() != Inpar::Mortar::shape_standard  // so for petrov-Galerkin and dual
          )
          {
            // loop over all master elements associated with this slave element
            for (int m = 0; m < (int)master_elements().size(); ++m)
            {
              // create Coupling2d object and push back
              coupling().push_back(std::make_shared<Coupling2d>(
                  idiscret_, dim_, quad_, imortar_, slave_element(), master_element(m)));

              // project the element pair
              coupling()[m]->project();

              // check for element overlap
              coupling()[m]->detect_overlap();
            }

            // calculate consistent dual shape functions for this element
            consistent_dual_shape();

            // do mortar integration
            for (int m = 0; m < (int)master_elements().size(); ++m)
              coupling()[m]->integrate_overlap(mparams_ptr);
          }

          else
          {
            for (int m = 0; m < (int)master_elements().size(); ++m)
            {
              // create Coupling2d object and push back
              coupling().push_back(std::make_shared<Coupling2d>(
                  idiscret_, dim_, quad_, imortar_, slave_element(), master_element(m)));

              // project the element pair
              coupling()[m]->project();

              // check for element overlap
              coupling()[m]->detect_overlap();

              // integrate the element overlap
              coupling()[m]->integrate_overlap(mparams_ptr);
            }
          }
        }
        else
        {
          // nothing
        }
      }
      else
      {
        // nothing
      }
    }
  }
  //**********************************************************************
  // INVALID
  //**********************************************************************
  else
  {
    FOUR_C_THROW("Invalid type of numerical integration");
  }

  // free memory of consistent dual shape function coefficient matrix
  slave_element().mo_data().reset_dual_shape();
  slave_element().mo_data().reset_deriv_dual_shape();
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool Mortar::Coupling2dManager::evaluate_coupling(
    const std::shared_ptr<Mortar::ParamsInterface>& mparams_ptr)
{
  if (master_elements().size() == 0) return false;

  // decide which type of coupling should be evaluated
  auto algo = Teuchos::getIntegralValue<Inpar::Mortar::AlgorithmType>(imortar_, "ALGORITHM");

  //*********************************
  // Mortar Contact
  //*********************************
  if (algo == Inpar::Mortar::algorithm_mortar or algo == Inpar::Mortar::algorithm_gpts)
  {
    integrate_coupling(mparams_ptr);
  }

  //*********************************
  // Error
  //*********************************
  else
    FOUR_C_THROW("chose contact algorithm not supported!");

  return true;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Mortar::Coupling2dManager::consistent_dual_shape()
{
  // For standard shape functions no modification is necessary
  // A switch earlier in the process improves computational efficiency
  auto consistent =
      Teuchos::getIntegralValue<Inpar::Mortar::ConsistentDualType>(imortar_, "LM_DUAL_CONSISTENT");
  if (shape_fcn() == Inpar::Mortar::shape_standard || consistent == Inpar::Mortar::consistent_none)
    return;

  // Consistent modification not yet checked for constant LM interpolation
  if (quad() && lag_mult_quad() == Inpar::Mortar::lagmult_const &&
      consistent != Inpar::Mortar::consistent_none)
    FOUR_C_THROW(
        "ERROR: Consistent dual shape functions not yet checked for constant LM interpolation!");

  // not implemented for nurbs yet
  if (slave_element().shape() == Core::FE::CellType::nurbs3)
    FOUR_C_THROW("Consistent dual shape functions not yet implemented for nurbs");

  // do nothing if there are no coupling pairs
  if (coupling().size() == 0) return;

  // detect entire overlap
  double ximin = 1.0;
  double ximax = -1.0;
  std::map<int, double> dximin;
  std::map<int, double> dximax;
  std::vector<std::map<int, double>> ximaps(4);

  // loop over all master elements associated with this slave element
  for (const auto& coupling : coupling())
  {
    // go on, if this s/m pair has no overlap
    if (not coupling->overlap()) continue;

    double sxia = coupling->xi_proj()[0];
    double sxib = coupling->xi_proj()[1];

    // get element contact integration area
    // and for contact derivatives of beginning and end
    if (sxia < ximin) ximin = sxia;
    if (sxib > ximax) ximax = sxib;
  }

  // no overlap: the applied dual shape functions don't matter, as the integration domain is void
  if ((ximax == -1.0 && ximin == 1.0) || (ximax - ximin < 4. * MORTARINTLIM)) return;

  // fully projecting element: no modification necessary
  if (ximin == -1. && ximax == 1.) return;

  // calculate consistent dual schape functions (see e.g. Cichosz et.al.:
  // Consistent treatment of boundaries with mortar contact formulations, CMAME 2010

  // get number of nodes of present slave element
  int nnodes = slave_element().num_node();

  // compute entries to bi-ortho matrices me/de with Gauss quadrature
  Mortar::ElementIntegrator integrator(slave_element().shape());

  // prepare for calculation of dual shape functions
  Core::LinAlg::SerialDenseMatrix me(nnodes, nnodes, true);
  Core::LinAlg::SerialDenseMatrix de(nnodes, nnodes, true);

  for (int gp = 0; gp < integrator.n_gp(); ++gp)
  {
    Core::LinAlg::SerialDenseVector sval(nnodes);
    Core::LinAlg::SerialDenseMatrix sderiv(nnodes, 1, true);

    // coordinates and weight
    std::array<double, 2> eta = {integrator.coordinate(gp, 0), 0.0};
    double wgt = integrator.weight(gp);

    // coordinate transformation sxi->eta (slave Mortar::Element->Overlap)
    std::array<double, 2> sxi = {0.0, 0.0};
    sxi[0] = 0.5 * (1.0 - eta[0]) * ximin + 0.5 * (1.0 + eta[0]) * ximax;

    // evaluate trace space shape functions
    if (lag_mult_quad() == Inpar::Mortar::lagmult_lin)
    {
      slave_element().evaluate_shape_lag_mult_lin(
          Inpar::Mortar::shape_standard, sxi.data(), sval, sderiv, nnodes);
    }
    else
      slave_element().evaluate_shape(sxi.data(), sval, sderiv, nnodes);

    // evaluate the two slave side Jacobians
    double dxdsxi = slave_element().jacobian(sxi.data());
    double dsxideta = -0.5 * ximin + 0.5 * ximax;

    // integrate dual shape matrices de, me and their linearizations
    for (int j = 0; j < nnodes; ++j)
    {
      // de and linearization
      de(j, j) += wgt * sval[j] * dxdsxi * dsxideta;

      // me and linearization
      for (int k = 0; k < nnodes; ++k)
      {
        me(j, k) += wgt * sval[j] * sval[k] * dxdsxi * dsxideta;
      }
    }
  }

  // declare dual shape functions coefficient matrix
  Core::LinAlg::SerialDenseMatrix ae(nnodes, nnodes, true);

  // compute matrix A_e for linear interpolation of quadratic element
  if (lag_mult_quad() == Inpar::Mortar::lagmult_lin)
  {
    // how many non-zero nodes
    const int nnodeslin = 2;

    // reduce me to non-zero nodes before inverting
    Core::LinAlg::Matrix<nnodeslin, nnodeslin> melin;
    for (int j = 0; j < nnodeslin; ++j)
      for (int k = 0; k < nnodeslin; ++k) melin(j, k) = me(j, k);

    // invert bi-ortho matrix melin
    Core::LinAlg::inverse(melin);

    // re-inflate inverse of melin to full size
    Core::LinAlg::SerialDenseMatrix invme(nnodes, nnodes, true);
    for (int j = 0; j < nnodeslin; ++j)
      for (int k = 0; k < nnodeslin; ++k) invme(j, k) = melin(j, k);

    // get solution matrix with dual parameters
    Core::LinAlg::multiply(ae, de, invme);
  }
  // compute matrix A_e for all other cases
  else
    Core::LinAlg::invert_and_multiply_by_cholesky(me, de, ae);

  // store ae matrix in slave element data container
  slave_element().mo_data().dual_shape() = std::make_shared<Core::LinAlg::SerialDenseMatrix>(ae);
}

FOUR_C_NAMESPACE_CLOSE
