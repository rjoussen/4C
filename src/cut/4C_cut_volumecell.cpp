// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_cut_volumecell.hpp"

#include "4C_cut_boundarycell.hpp"
#include "4C_cut_cycle.hpp"
#include "4C_cut_direct_divergence.hpp"
#include "4C_cut_integrationcell.hpp"
#include "4C_cut_kernel.hpp"
#include "4C_cut_mesh.hpp"
#include "4C_cut_options.hpp"
#include "4C_cut_point.hpp"
#include "4C_cut_tetmesh.hpp"
#include "4C_cut_triangulateFacet.hpp"
#include "4C_cut_utils.hpp"
#include "4C_cut_volume_integration.hpp"

#include <Teuchos_TimeMonitor.hpp>

#include <algorithm>

FOUR_C_NAMESPACE_OPEN

int Cut::VolumeCell::hex8totet4[5][4] = {
    {0, 1, 3, 4}, {1, 2, 3, 6}, {4, 5, 1, 6}, {6, 7, 3, 4}, {1, 6, 3, 4}};

int Cut::VolumeCell::wedge6totet4[3][4] = {{0, 1, 2, 3}, {3, 4, 1, 5}, {1, 5, 2, 3}};


int Cut::VolumeCell::pyramid5totet4[2][4] = {{0, 1, 3, 4}, {1, 2, 3, 4}};


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Cut::VolumeCell::VolumeCell(const plain_facet_set& facets,
    const std::map<std::pair<Point*, Point*>, plain_facet_set>& volume_lines, Element* element)
    : element_(element),
      position_(Point::undecided),
      facets_(facets),
      volume_(0.0),
      is_negligible_small_(false)
{
  for (plain_facet_set::const_iterator i = facets_.begin(); i != facets_.end(); ++i)
  {
    Facet* f = *i;
    f->register_entity(this);
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool Cut::VolumeCell::is_equal(const plain_facet_set& vcell) const
{
  bool isequal = false;
  for (plain_facet_set::const_iterator ci = facets_.begin(); ci != facets_.end(); ++ci)
  {
    const Facet& f = **ci;
    for (plain_facet_set::const_iterator cii = vcell.begin(); cii != vcell.end(); ++cii)
    {
      const Facet& ff = **cii;
      const std::vector<Point*>& ffpoints = ff.points();

      isequal = (f.points().size() == ffpoints.size() and f.contains(ffpoints));
      if (isequal) break;
    }

    if (not isequal) return false;
  }

  // if it reaches this point, the volume cells are identical
  return isequal;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Cut::VolumeCell::neighbors(Point* p, const plain_volumecell_set& cells,
    const plain_volumecell_set& done, plain_volumecell_set& connected, plain_element_set& elements)
{
  if (done.count(this) == 0)
  {
    // this volume is included
    connected.insert(this);
    elements.insert(element_);

    // Do the facets that include the point first. This ensures we choose the
    // right volumes (the ones attached to the point), if there are multiple
    // connections possible (we are faced with a thin structure cut.)

    for (plain_facet_set::const_iterator i = facets_.begin(); i != facets_.end(); ++i)
    {
      Facet* f = *i;
      if (p == nullptr or f->contains(p))
      {
        f->neighbors(p, cells, done, connected, elements);
      }
    }

    if (p != nullptr)
    {
      for (plain_facet_set::const_iterator i = facets_.begin(); i != facets_.end(); ++i)
      {
        Facet* f = *i;
        if (not f->contains(p))
        {
          f->neighbors(p, cells, done, connected, elements);
        }
      }
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
// without check for elements
void Cut::VolumeCell::neighbors(Point* p, const plain_volumecell_set& cells,
    const plain_volumecell_set& done, plain_volumecell_set& connected)
{
  if (done.count(this) == 0)
  {
    // this volume is included
    connected.insert(this);

    // Do the facets that include the point first. This ensures we choose the
    // right volumes (the ones attached to the point), if there are multiple
    // connections possible (we are faced with a thin structure cut.)

    for (plain_facet_set::const_iterator i = facets_.begin(); i != facets_.end(); ++i)
    {
      Facet* f = *i;

      if (p == nullptr or f->contains(p))
      {
        f->neighbors(p, cells, done, connected);
      }
    }

    if (p != nullptr)
    {
      for (plain_facet_set::const_iterator i = facets_.begin(); i != facets_.end(); ++i)
      {
        Facet* f = *i;
        if (not f->contains(p))
        {
          f->neighbors(p, cells, done, connected);
        }
      }
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Cut::VolumeCell::get_all_points(Mesh& mesh, PointSet& cut_points)
{
  for (plain_facet_set::iterator i = facets_.begin(); i != facets_.end(); ++i)
  {
    Facet* f = *i;
    f->get_all_points(mesh, cut_points);
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool Cut::VolumeCell::contains(Point* p)
{
  for (plain_facet_set::const_iterator i = facets_.begin(); i != facets_.end(); ++i)
  {
    Facet* f = *i;
    if (f->contains(p))
    {
      return true;
    }
  }
  return false;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool Cut::VolumeCell::contains(Core::LinAlg::Matrix<3, 1>& x)
{
  if (integrationcells_.size() == 0)
    FOUR_C_THROW(
        "no integrationcells for volumecell stored, implement Contains check without "
        "integrationcells");

  for (Cut::plain_integrationcell_set::iterator it = integrationcells_.begin();
      it != integrationcells_.end(); it++)
  {
    Cut::IntegrationCell* intcell = *it;

    if (intcell->contains(x)) return true;
  }

  return false;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Cut::VolumeCell::create_tet4_integration_cells(Mesh& mesh,
    const std::vector<std::vector<Point*>>& tets,
    const std::map<Facet*, std::vector<Point*>>& sides_xyz)
{
  for (std::vector<std::vector<Point*>>::const_iterator i = tets.begin(); i != tets.end(); ++i)
  {
    const std::vector<Point*>& tet = *i;
    if (tet.size() != 4)
    {
      FOUR_C_THROW("tet expected");
    }
    new_tet4_cell(mesh, tet);
  }

  for (std::map<Facet*, std::vector<Point*>>::const_iterator i = sides_xyz.begin();
      i != sides_xyz.end(); ++i)
  {
    Facet* f = i->first;
    const std::vector<Point*>& points = i->second;

    std::size_t length = points.size();
    if (length % 3 != 0) FOUR_C_THROW("expect list of triangles");

    length /= 3;
    std::vector<Point*> p(3);
    for (std::size_t i = 0; i < length; ++i)  // loop the list of triangles
    {
      std::copy_n(points.data() + 3 * i, 3, p.data());
      new_tri3_cell(mesh, f, p);  // create tri3 cell
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Cut::VolumeCell::get_integration_cells(plain_integrationcell_set& cells)
{
  std::copy(
      integrationcells_.begin(), integrationcells_.end(), std::inserter(cells, cells.begin()));
}

/// get a map of boundary cells for all cutting sides, key= side-Id, value= vector of boundary cells
/// note that the boundary cells of subsides with the same side id are stored now in one key
void Cut::VolumeCell::get_boundary_cells(std::map<int, std::vector<Cut::BoundaryCell*>>& bcells)
{
  for (plain_boundarycell_set::iterator i = bcells_.begin(); i != bcells_.end(); ++i)
  {
    BoundaryCell* bc = *i;
    Facet* f = bc->get_facet();
    if (f->on_cut_side())
    {
      int sid = f->side_id();  // f->OnCutSide => sid>-1
      // usually there are more facets with the same side id as the cutting sides have been
      // subdivided into subsides which produce own facets for each subside
      bcells[sid].push_back(bc);
    }
  }
}


/// get a map of boundary cells for all cutting sides, key= side-Id, value= vector of boundary cells
/// note that the boundary cells of subsides with the same side id are stored now in one key
void Cut::VolumeCell::get_boundary_cells_to_be_integrated(
    std::map<int, std::vector<Cut::BoundaryCell*>>& bcells)
{
  for (plain_boundarycell_set::iterator i = bcells_.begin(); i != bcells_.end(); ++i)
  {
    BoundaryCell* bc = *i;
    Facet* f = bc->get_facet();
    // Get all bc's for cuts from only the outside vc's
    //  as to not integrate twice over the same surface
    if ((f->on_cut_side() and position() == Cut::Point::outside))
    {
      int sid = f->side_id();  // f->OnCutSide => sid>-1
      // usually there are more facets with the same side id as the cutting sides have been
      // subdivided into subsides which produce own facets for each subside
      bcells[sid].push_back(bc);

      //@ Christoph: Here one could add marked sides for the cut-mesh!
    }
    else if (f->on_marked_background_side())
    {
      // Loop over all marked actions and extract bc's for corresponding coupling object.
      for (std::map<Cut::MarkedActions, int>::iterator markit =
               f->parent_side()->get_markedsidemap().begin();
          markit != f->parent_side()->get_markedsidemap().end(); ++markit)
      {
        if (markit->first == Cut::mark_and_create_boundarycells)
          bcells[markit->second].push_back(bc);
      }
    }
  }
}

/*----------------------------------------------------------------------------*
 * SideId() of Facet (used for timeintegration)
 *----------------------------------------------------------------------------*/
void Cut::VolumeCell::collect_cut_sides(plain_int_set& cutside_ids)
{
  for (plain_facet_set::iterator i = facets_.begin(); i != facets_.end(); ++i)
  {
    Facet* f = *i;
    if (f->on_cut_side())
    {
      int sid = f->side_id();  // f->OnCutSide => sid>-1
      // usually there are more facets with the same side id as the cutting sides have been
      // subdivided into subsides which produce own facets for each subside
      cutside_ids.insert(sid);
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
int Cut::VolumeCell::get_parent_element_id() const { return element_->get_parent_id(); }

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Cut::VolumeCell::connect_nodal_dof_sets(bool include_inner)
{
  //   if ( Empty() )
  //     return;
  if (not include_inner and position() != Point::outside) return;

  const std::vector<Node*>& nodes = element_->nodes();
  nodaldofset_.reserve(nodes.size());

  for (std::vector<Node*>::const_iterator i = nodes.begin(); i != nodes.end(); ++i)
  {
    Node* n = *i;
    nodaldofset_.push_back(n->dof_set_number(this));
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Cut::VolumeCell::position(Point::PointPosition position)
{
  if (position_ != position)
  {
    position_ = position;

    for (plain_facet_set::const_iterator i = facets_.begin(); i != facets_.end(); ++i)
    {
      Facet* f = *i;
      Point::PointPosition fp = f->position();
      if (fp == Point::undecided)
      {
        f->position(position);
      }
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Cut::VolumeCell::print(std::ostream& stream) const
{
  stream << "\n==========================================\n";
  stream << "=== VolumeCell ( address: " << std::setw(10) << this << " ) ===\n";
  stream << "# VolumeCell: "
         << " pos: " << Point::point_position_to_string(position_) << " "
         << "#facets: " << facets_.size() << " "
         << "#intcells: " << integrationcells_.size() << " "
         << "#bcells: " << bcells_.size() << "\n";
  unsigned count = 0;
  for (plain_facet_set::const_iterator i = facets_.begin(); i != facets_.end(); ++i)
  {
    Facet* f = *i;
    stream << "\n# Facet " << count++ << " of VolumeCell:\n";
    f->print(stream);
  }

  count = 0;
  for (plain_boundarycell_set::const_iterator i = bcells_.begin(); i != bcells_.end(); ++i)
  {
    BoundaryCell* bcell = *i;
    stream << "\n# BoundaryCell " << count++ << " of VolumeCell:\n";
    bcell->print(stream);
  }

  count = 0;
  for (plain_integrationcell_set::const_iterator i = integrationcells_.begin();
      i != integrationcells_.end(); ++i)
  {
    IntegrationCell* icell = *i;
    stream << "\n# IntegrationCell " << count++ << " of VolumeCell:\n";
    icell->print(stream);
  }

  stream << "\n";
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Cut::VolumeCell::new_boundary_cell(
    Mesh& mesh, Core::FE::CellType shape, Facet* f, const std::vector<Point*>& x)
{
  if (facets_.count(f) == 0)
  {
    FOUR_C_THROW("facet does not belong to volume cell");
  }
  switch (shape)
  {
    case Core::FE::CellType::point1:
      new_point1_cell(mesh, f, x);
      break;
    case Core::FE::CellType::line2:
      new_line2_cell(mesh, f, x);
      break;
    case Core::FE::CellType::tri3:
      new_tri3_cell(mesh, f, x);
      break;
    case Core::FE::CellType::quad4:
      new_quad4_cell(mesh, f, x);
      break;
    default:
      FOUR_C_THROW("Unsupported shape ( shape = {} )", Core::FE::cell_type_to_string(shape));
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Cut::VolumeCell::new_point1_cell(Mesh& mesh, Facet* f, const std::vector<Point*>& x)
{
  f->new_point1_cell(mesh, this, x, bcells_);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Cut::VolumeCell::new_line2_cell(Mesh& mesh, Facet* f, const std::vector<Point*>& x)
{
  f->new_line2_cell(mesh, this, x, bcells_);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Cut::VolumeCell::new_tri3_cell(Mesh& mesh, Facet* f, const std::vector<Point*>& x)
{
  f->new_tri3_cell(mesh, this, x, bcells_);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Cut::VolumeCell::new_quad4_cell(Mesh& mesh, Facet* f, const std::vector<Point*>& x)
{
  f->new_quad4_cell(mesh, this, x, bcells_);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Cut::VolumeCell::new_arbitrary_cell(Mesh& mesh, Facet* f, const std::vector<Point*>& x,
    const Core::FE::GaussIntegration& gp, const Core::LinAlg::Matrix<3, 1>& normal)
{
  f->new_arbitrary_cell(mesh, this, x, bcells_, gp, normal);
}

/*double Cut::VolumeCell::Volume()
{
  double volume = 0;
  for ( plain_integrationcell_set::iterator i=integrationcells_.begin(); i!=integrationcells_.end();
++i )
  {
    IntegrationCell * ic = *i;
    volume += ic->Volume();
  }
  return volume;
}*/

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
int Cut::VolumeCell::num_gauss_points(Core::FE::CellType shape)
{
  int numgp = 0;

  for (plain_integrationcell_set::const_iterator i = integrationcells_.begin();
      i != integrationcells_.end(); ++i)
  {
    IntegrationCell* ic = *i;

    // Create (unmodified) gauss points for integration cell with requested
    // polynomial order. This is supposed to be fast, since there is a cache.
    Core::FE::GaussIntegration gi(ic->shape(), ic->cubature_degree(shape));

    // we just need the number of points per cell
    numgp += gi.num_points();
  }

  return numgp;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Cut::VolumeCell::disconnect()
{
  for (plain_facet_set::iterator i = facets_.begin(); i != facets_.end(); ++i)
  {
    Facet* f = *i;
    f->disconnect_volume(this);
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Cut::VolumeCell::new_integration_cell(
    Mesh& mesh, Core::FE::CellType shape, const std::vector<Point*>& x)
{
  switch (shape)
  {
    // --- 1-D elements ---
    case Core::FE::CellType::line2:
      new_line2_cell(mesh, x);
      break;
    // --- 2-D elements ---
    case Core::FE::CellType::tri3:
      new_tri3_cell(mesh, x);
      break;
    case Core::FE::CellType::quad4:
      new_quad4_cell(mesh, x);
      break;
    // --- 3-D elements ---
    case Core::FE::CellType::hex8:
      new_hex8_cell(mesh, x);
      break;
    case Core::FE::CellType::tet4:
      new_tet4_cell(mesh, x);
      break;
    case Core::FE::CellType::wedge6:
      new_wedge6_cell(mesh, x);
      break;
    case Core::FE::CellType::pyramid5:
      new_pyramid5_cell(mesh, x);
      break;
    default:
      FOUR_C_THROW("Unsupported shape ( shape = {} )", Core::FE::cell_type_to_string(shape));
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Cut::VolumeCell::new_line2_cell(Mesh& mesh, const std::vector<Point*>& points)
{
  Point::PointPosition position = VolumeCell::position();
  integrationcells_.insert(mesh.new_line2_cell(position, points, this));
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Cut::VolumeCell::new_tri3_cell(Mesh& mesh, const std::vector<Point*>& points)
{
  Point::PointPosition position = VolumeCell::position();
  integrationcells_.insert(mesh.new_tri3_cell(position, points, this));
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Cut::VolumeCell::new_quad4_cell(Mesh& mesh, const std::vector<Point*>& points)
{
  Point::PointPosition position = VolumeCell::position();
  integrationcells_.insert(mesh.new_quad4_cell(position, points, this));
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Cut::VolumeCell::new_hex8_cell(Mesh& mesh, const std::vector<Point*>& points)
{
  Point::PointPosition position = VolumeCell::position();
  if (mesh.create_options().gen_hex8())
  {
    integrationcells_.insert(mesh.new_hex8_cell(position, points, this));
  }
  else
  {
    std::vector<Point*> tet4_points(4);
    for (int i = 0; i < 5; ++i)
    {
      set_tet_points(hex8totet4[i], points, tet4_points);
      integrationcells_.insert(mesh.new_tet4_cell(position, tet4_points, this));
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Cut::IntegrationCell* Cut::VolumeCell::new_tet4_cell(Mesh& mesh, const std::vector<Point*>& points)
{
  Point::PointPosition position = VolumeCell::position();
  IntegrationCell* ic = mesh.new_tet4_cell(position, points, this);
  integrationcells_.insert(ic);
  return ic;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Cut::VolumeCell::new_wedge6_cell(Mesh& mesh, const std::vector<Point*>& points)
{
  Point::PointPosition position = VolumeCell::position();
  if (mesh.create_options().gen_wedge6())
  {
    integrationcells_.insert(mesh.new_wedge6_cell(position, points, this));
  }
  else
  {
    std::vector<Point*> tet4_points(4);
    for (int i = 0; i < 3; ++i)
    {
      set_tet_points(wedge6totet4[i], points, tet4_points);
      integrationcells_.insert(mesh.new_tet4_cell(position, tet4_points, this));
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Cut::VolumeCell::new_pyramid5_cell(Mesh& mesh, const std::vector<Point*>& points)
{
  Point::PointPosition position = VolumeCell::position();
  if (mesh.create_options().gen_pyramid5())
  {
    integrationcells_.insert(mesh.new_pyramid5_cell(position, points, this));
  }
  else
  {
    std::vector<Point*> tet4_points(4);
    for (int i = 0; i < 2; ++i)
    {
      set_tet_points(pyramid5totet4[i], points, tet4_points);
      integrationcells_.insert(mesh.new_tet4_cell(position, tet4_points, this));
    }
  }
}

/*--------------------------------------------------------------------*
 * Check whether the point is inside, outside or on the boundary
 * of this volumecelll                                    sudhakar 07/12
 *--------------------------------------------------------------------*/
std::string Cut::VolumeCell::is_this_point_inside(Point* pt)
{
  Core::LinAlg::Matrix<3, 1> xglo;
  pt->coordinates(xglo.data());
  std::string inside = is_this_point_inside(xglo);
  return inside;
}

/*-----------------------------------------------------------------------------------------------*
 * Check whether the point with this global coordinates is inside, outside or on the boundary
 * of this volumecell                                                               sudhakar 07/12
 *-----------------------------------------------------------------------------------------------*/
std::string Cut::VolumeCell::is_this_point_inside(Core::LinAlg::Matrix<3, 1>& xglo)
{
  Core::LinAlg::Matrix<3, 1> xloc;
  element_->local_coordinates(xglo, xloc);

  const Cut::Point::PointPosition posi = position();
  if (posi == 0) FOUR_C_THROW("undefined position for the volumecell");

  VolumeIntegration vc(this, element_, posi, 0);
  std::string inside = vc.is_point_inside(xloc);
  return inside;
}

/*
 * There is a flaw in this test. If there are points added through TetMeshIntersection on the facet
 * boundary the test will fail.
 *
 * In order to fix this test one has to check if the remaining lines are connected. I.e. are the
 * points added from the tesselation connected correctly?
 *
 */
void Cut::VolumeCell::test_surface()
{
  if (empty())
  {
    // This is an artificial cell with zero volume. It should not exist in the
    // first place.
    return;
  }

  // see if all lines are closed
  //
  // This finds all the degenerated cases that where dropped before. Thus the
  // test complains a lot.

  for (plain_facet_set::iterator j = facets_.begin(); j != facets_.end(); ++j)
  {
    Facet* f = *j;

    if (f->on_cut_side())
    {
      if (f->is_triangulated())
      {
        //        std::cout << "f->Triangulation().size(): " << f->Triangulation().size();
        //
      }
      if (f->has_holes())
      {
        //
      }

      point_line_set lines;

      const std::vector<Point*>& points = f->points();
      Cycle cycle(points);
      cycle.add(lines);

      point_line_set facetlines = lines;

      for (plain_boundarycell_set::iterator i = bcells_.begin(); i != bcells_.end(); ++i)
      {
        BoundaryCell* bc = *i;
        if (bc->get_facet() == f)
        {
          //          std::cout << "Printing boundary cell: " << std::endl;
          //          bc->print();
          //          std::cout << std::endl;
          const std::vector<Point*>& points = bc->points();
          Cycle cycle(points);
          cycle.add(lines);
        }
      }

      if (lines.size() != 0)
      {
        //        std::cout << "Problem Facet: " << std::endl;
        //        f->print(std::cout);
        //        std::cout << std::endl;
        //        std::cout << "lines.size(): " << lines.size() << std::endl;
        //        for(unsigned k=0; k< lines.size(); k++)
        //        {
        //          std::cout << "#: " << k <<std::endl;
        //          lines[k].first->print(std::cout);
        //          lines[k].second->print(std::cout);
        //          std::cout << std::endl;
        //        }
        //
        // Q: What line in the facetlines is not connected?
        int numberoflines = 0;
        std::vector<int> facetlineindex;
        for (unsigned k = 0; k < lines.size(); k++)
        {
          for (unsigned l = 0; l < facetlines.size(); l++)
          {
            FOUR_C_THROW(
                "not supported when using std::set in definition of point_line_set. Adapt this. "
                "Anyway this routine is not bugfree! This about before using this function!!!");
          }
        }
        std::cout << "numberoflines: " << numberoflines << std::endl;

        // test that no line is the same.
        for (unsigned k = 0; k < facetlineindex.size(); k++)
        {
          for (unsigned l = 0; l < facetlineindex.size(); l++)
            if (facetlineindex[k] == facetlineindex[l] and k != l)
              FOUR_C_THROW("volume cut facets not closed!!");
        }

        //        //Find the connection.
        //        for(unsigned k=0; k< numberoflines; k++)
        //        {
        //          int index1 = facetlines[facetlineindex[k]].first->Id();
        //          int index2 = facetlines[facetlineindex[k]].second->Id();
        //
        //
        //
        //
        //        }


        // One would have to check if the facetline(s) found here is/are connected in lines. Would
        // probably be best implemented
        // with some sort of tree structure.

        FOUR_C_THROW("volume cut facets not closed");
      }
    }
  }
}

/*-------------------------------------------------------------------------------------*
                Write the bounding lines of volumecell details for visualization
                Gausspoints of moment fitting are not included
*--------------------------------------------------------------------------------------*/
void Cut::VolumeCell::dump_gmsh(std::ofstream& file)
{
  const plain_facet_set& facete = facets();

  file << "View \"Volume Cell \" {\n";
  for (plain_facet_set::const_iterator j = facete.begin(); j != facete.end(); j++)
  {
    Facet* ref = *j;
#ifndef OUTPUT_GLOBAL_DIVERGENCE_CELLS
    std::vector<std::vector<double>> corners;
    ref->CornerPointsLocal(parent_element(), corners);
#else
    const std::vector<std::vector<double>> corners = ref->corner_points_global(parent_element());
#endif

    for (unsigned i = 0; i < corners.size(); i++)
    {
      const std::vector<double> coords1 = corners[i];
      const std::vector<double> coords2 = corners[(i + 1) % corners.size()];
      file << "SL(" << coords1[0] << "," << coords1[1] << "," << coords1[2] << "," << coords2[0]
           << "," << coords2[1] << "," << coords2[2] << ")"
           << "{0,0};\n";
    }
  }
  file << "};\n";
  file << "View[PostProcessing.NbViews-1].ColorTable = { {0,0,255} };\n";  // Changing color to red
  file << "View[PostProcessing.NbViews-1].Light=0;\n";                     // Disable the lighting
  file << "View[PostProcessing.NbViews-1].ShowScale=0;\n";                 // Disable legend
  file << "View[PostProcessing.NbViews-1].LineWidth = 3.0;";               // increase line width
}

/*---------------------------------------------------------------------------------------*
 * Write the geometry of the volumecell based on facet surfaces            sudhakar 07/15
 * Can be used to check if the geometry of vc is correct or not
 *---------------------------------------------------------------------------------------*/
void Cut::VolumeCell::dump_gmsh_solid(std::ofstream& file, Mesh& mesh)
{
  const plain_facet_set& facete = facets();
  file << "View \"Volume Cell \" {\n";
  for (plain_facet_set::const_iterator j = facete.begin(); j != facete.end(); j++)
  {
    Facet* fac = *j;
    std::vector<Point*> corners = fac->corner_points();

    if (corners.size() == 3)
    {
      file << "ST(";
      for (unsigned ipt = 0; ipt < corners.size(); ipt++)
      {
        Point* pt = corners[ipt];
        const double* x = pt->x();
        file << x[0] << "," << x[1] << "," << x[2];
        if (ipt != (corners.size() - 1)) file << ",";
      }
      file << "){0.0,0.0,0.0};\n";
    }
    else
    {
      if (!fac->is_triangulated()) fac->do_triangulation(mesh, corners);
      const std::vector<std::vector<Point*>>& triangulation = fac->triangulation();

      for (std::vector<std::vector<Point*>>::const_iterator j = triangulation.begin();
          j != triangulation.end(); ++j)
      {
        std::vector<Point*> tri = *j;

        if (tri.size() == 3)
          file << "ST(";
        else
          file << "SQ(";

        for (unsigned ipt = 0; ipt < tri.size(); ipt++)
        {
          Point* pt = tri[ipt];
          const double* x = pt->x();
          file << x[0] << "," << x[1] << "," << x[2];
          if (ipt != (tri.size() - 1)) file << ",";
        }
        if (tri.size() == 3)
          file << "){0.0,0.0,0.0};\n";
        else
          file << "){0.0,0.0,0.0,0.0};\n";
      }
    }
  }
  file << "};\n";
}

/*--------------------------------------------------------------------------------------------------------*
        write the boundaries of volumecell and the positions of Gauss points when using
        Moment fitting for visualization a separate file with "Mom_volcell" prefix is generated
        for every volumecell as the gausspoint distribution can be clearly seen
*---------------------------------------------------------------------------------------------------------*/
void Cut::VolumeCell::dump_gmsh_gauss_points_mom_fit(
    const std::vector<std::vector<double>>& gauspts)
{
  static int sideno = 0;
  sideno++;

  std::stringstream str;
  str << "Mom_volcell" << sideno << ".pos";
  std::ofstream file(str.str().c_str());

  dump_gmsh(file);

  file << "View \" Gauss Points \" {\n";
  for (unsigned i = 0; i < gauspts.size(); i++)
  {
    file << "SP(" << gauspts[i][0] << "," << gauspts[i][1] << "," << gauspts[i][2] << ",1){0.0};\n";
  }
  file << "};\n";
  file << "View[PostProcessing.NbViews-1].ColorTable = { {0,0,0} };\n";  // Changing color to black
  file << "View[PostProcessing.NbViews-1].Light=0;\n";                   // Disable the lighting
  file << "View[PostProcessing.NbViews-1].ShowScale=0;\n";               // Disable legend
  file << "View[PostProcessing.NbViews-1].PointSize = 4.0;";             // fix point size
  file.close();
}

/*--------------------------------------------------------------------------------------------------------*
        write the boundaries of volumecell and the positions of Gauss points for visualization
        a separate file when using tessellation
*---------------------------------------------------------------------------------------------------------*/
void Cut::VolumeCell::dump_gmsh_gauss_points_tessellation()
{
  static int sideno = 0;
  sideno++;

  std::stringstream str;
  str << "MTes_volcell" << sideno << ".pos";
  std::ofstream file(str.str().c_str());

  dump_gmsh(file);

  std::shared_ptr<Core::FE::GaussPointsComposite> gpc =
      std::make_shared<Core::FE::GaussPointsComposite>(0);

  const plain_integrationcell_set& cells = integration_cells();
  for (plain_integrationcell_set::const_iterator i = cells.begin(); i != cells.end(); ++i)
  {
    Cut::IntegrationCell* ic = *i;
    switch (ic->shape())
    {
      case Core::FE::CellType::hex8:
      {
        std::shared_ptr<Core::FE::GaussPoints> gp = create_projected<Core::FE::CellType::hex8>(ic);
        gpc->append(gp);
        break;
      }
      case Core::FE::CellType::tet4:
      {
        std::shared_ptr<Core::FE::GaussPoints> gp = create_projected<Core::FE::CellType::tet4>(ic);
        gpc->append(gp);
        break;
      }
      default:
      {
        FOUR_C_THROW("Include this element here");
        break;
      }
    }
  }

  Core::FE::GaussIntegration gpv(gpc);

  file << "View \" Gauss Points \" {\n";
  for (Core::FE::GaussIntegration::iterator iquad = gpv.begin(); iquad != gpv.end(); ++iquad)
  {
    const Core::LinAlg::Matrix<3, 1> eta(iquad.point());
    file << "SP(" << eta(0, 0) << "," << eta(1, 0) << "," << eta(2, 0) << ",1){0.0};\n";
  }

  file << "};\n";
  file << "View[PostProcessing.NbViews-1].ColorTable = { {0,0,0} };\n";  // Changing color to black
  file << "View[PostProcessing.NbViews-1].Light=0;\n";                   // Disable the lighting
  file << "View[PostProcessing.NbViews-1].ShowScale=0;\n";               // Disable legend
  file << "View[PostProcessing.NbViews-1].PointSize = 4.0;";             // fix point size
  file.close();
}

/*----------------------------------------------------------------------------------------------------------*
 * Perform integration of a pre-defined function over this vc using gauss points          Sudhakar
 *01/13 generated from tessellation
 *----------------------------------------------------------------------------------------------------------*/
void Cut::VolumeCell::integrate_specific_functions_tessellation()
{
  std::shared_ptr<Core::FE::GaussPointsComposite> gpc =
      std::make_shared<Core::FE::GaussPointsComposite>(0);

  const plain_integrationcell_set& cells = integration_cells();
  for (plain_integrationcell_set::const_iterator i = cells.begin(); i != cells.end(); ++i)
  {
    Cut::IntegrationCell* ic = *i;
    switch (ic->shape())
    {
      case Core::FE::CellType::hex8:
      {
        std::shared_ptr<Core::FE::GaussPoints> gp = create_projected<Core::FE::CellType::hex8>(ic);
        gpc->append(gp);
        break;
      }
      case Core::FE::CellType::tet4:
      {
        std::shared_ptr<Core::FE::GaussPoints> gp = create_projected<Core::FE::CellType::tet4>(ic);
        gpc->append(gp);
        break;
      }
      default:
      {
        FOUR_C_THROW("Include this element here");
        break;
      }
    }
  }

  Core::FE::GaussIntegration gpv(gpc);

  double intVal = 0.0;
  for (Core::FE::GaussIntegration::iterator iquad = gpv.begin(); iquad != gpv.end(); ++iquad)
  {
    double weight = iquad.weight();

    const Core::LinAlg::Matrix<3, 1> eta(iquad.point());
    double xx = eta(0, 0);
    double yy = eta(1, 0);
    double zz = eta(2, 0);

    intVal +=
        (pow(xx, 6) + xx * pow(yy, 4) * zz + xx * xx * yy * yy * zz * zz + pow(zz, 6)) * weight;
  }
  std::cout << std::setprecision(20) << "TESSELLATION Integration = " << intVal << "\n";
}

template <Core::FE::CellType distype>
std::shared_ptr<Core::FE::GaussPoints> Cut::VolumeCell::create_projected(Cut::IntegrationCell* ic)
{
  const unsigned nen = Core::FE::num_nodes(distype);

  Core::LinAlg::Matrix<3, nen> xie;

  const std::vector<Cut::Point*>& cpoints = ic->points();
  if (cpoints.size() != nen) FOUR_C_THROW("non-matching number of points");

  for (unsigned i = 0; i < nen; ++i)
  {
    Cut::Point* p = cpoints[i];
    Core::LinAlg::Matrix<3, 1> xg, xi;
    p->coordinates(xg.data());
    element_->local_coordinates(xg, xi);
    std::copy(xi.data(), xi.data() + 3, &xie(0, i));
  }

  std::shared_ptr<Core::FE::GaussPoints> gp = Core::FE::GaussIntegration::create_projected<distype>(
      xie, ic->cubature_degree(element_->shape()));
  return gp;
}

/*------------------------------------------------------------------------------------------------------*
    convert the Gaussian points and weights into appropriate Gauss rule as per 4C implementation
*-------------------------------------------------------------------------------------------------------*/
std::shared_ptr<Core::FE::GaussPoints> Cut::VolumeCell::gauss_points_fitting()
{
  std::shared_ptr<Core::FE::CollectedGaussPoints> cgp =
      std::make_shared<Core::FE::CollectedGaussPoints>(0);

  for (unsigned i = 0; i < gauss_pts_.size(); i++)
  {
    Core::LinAlg::Matrix<3, 1> xe, xei;
    xe(0, 0) = gauss_pts_[i][0];
    xe(1, 0) = gauss_pts_[i][1];
    xe(2, 0) = gauss_pts_[i][2];

    cgp->append(xe, weights_(i));
  }

  return cgp;
}

/*--------------------------------------------------------------------------------------------*
                 Generate boundary cells for the cut facets of the volumecell
*---------------------------------------------------------------------------------------------*/
void Cut::VolumeCell::generate_boundary_cells(Mesh& mesh, const Cut::Point::PointPosition posi,
    Element* elem, int BaseNos, Cut::BCellGaussPts BCellgausstype)
{
  // TEUCHOS_FUNC_TIME_MONITOR( "Cut::VolumeCell::generate_boundary_cells" );


  // TODO: we have to restructure the creation of boundary cells.
  // bcs should not be stored for a volumecell but for the respective facet, see comments in
  // f->GetBoundaryCells

  const plain_facet_set& facete = facets();
  for (plain_facet_set::const_iterator i = facete.begin(); i != facete.end(); i++)
  {
    Facet* fac = *i;

    if (fac->on_boundary_cell_side() == false)  // we need boundary cells only for the cut facets
      continue;

    // For LevelSetSides generate Boundary Cells in own way.
    if (fac->belongs_to_level_set_side())
    {
      generate_boundary_cells_level_set_side(mesh, posi, elem, fac, BaseNos, BCellgausstype);
      continue;
    }

    //--------------------------------------------------------------------
    // Normal vector from parent side is used to identify whether normals
    // from facet is in appropriate direction or not
    //--------------------------------------------------------------------
    const Side* parside = fac->parent_side();
    const std::vector<Node*>& par_nodes = parside->nodes();

    std::vector<double> eqnpar(4), eqnfac(4);
    bool reverse = false;

    std::vector<Point*> pts_par(par_nodes.size());
    for (unsigned parnode = 0; parnode < par_nodes.size(); parnode++)
      pts_par[parnode] = par_nodes[parnode]->point();
    eqnpar = Kernel::eqn_plane_of_polygon(pts_par);

    std::vector<Point*> corners = fac->corner_points();
    std::vector<Point*> cornersTemp(corners);

    // when finding eqn of plane for the facet, inline points should not be taken
    Cut::Kernel::delete_inline_pts(cornersTemp);

    if (cornersTemp.size() != 0)
    {
      eqnfac = Kernel::eqn_plane_of_polygon(corners);
      reverse = to_reverse(posi, eqnpar, eqnfac);
    }

    // For Marked sides the boundary-cells on outside vc's need to be the same as for inside.
    if (fac->on_marked_background_side() and posi == Cut::Point::outside) reverse = !reverse;

    if (reverse)  // normal from facet is in wrong direction
    {
      std::reverse(corners.begin(), corners.end());  // change ordering to correct this
      std::reverse(cornersTemp.begin(),
          cornersTemp.end());  // Todo: Ager ... I don't understand the sense of this cornersTemp???
    }

    // if no of corners are 3 or 4, just add them as boundary integrationcells directly
    if (corners.size() == 3)
    {
      double areaCell = Cut::Kernel::get_area_tri(corners);
      if (areaCell <
          REF_AREA_BCELL)  // What is this Ref_AREA_CELL TOLERANCE?! Make it viable for GLOBAL!!!
        continue;
      new_tri3_cell(mesh, fac, corners);
    }
    else
    {
      if (BCellgausstype == Cut::BCellGaussPts_Tessellation)  // generate boundarycell
                                                              // gausspoints by triangulation
      {
        if (!fac->is_triangulated()) fac->do_triangulation(mesh, corners);
        const std::vector<std::vector<Point*>>& triangulation = fac->triangulation();

        for (std::vector<std::vector<Point*>>::const_iterator j = triangulation.begin();
            j != triangulation.end(); ++j)
        {
          std::vector<Point*> tri = *j;

          if (tri.size() == 3)
          {
            double areaCell = Cut::Kernel::get_area_tri(tri);
            if (areaCell < REF_AREA_BCELL) continue;
            new_tri3_cell(mesh, fac, tri);
          }
          else if (tri.size() == 4)
          {
            double areaCell = Cut::Kernel::get_area_convex_quad(tri);
            if (areaCell < REF_AREA_BCELL) continue;
            new_quad4_cell(mesh, fac, tri);
          }
          else
            FOUR_C_THROW("Triangulation created neither tri3 or quad4");
        }
      }
      else if (BCellgausstype ==
               Cut::BCellGaussPts_MomentFitting)  // generate boundarycell gausspoints by
                                                  // solving moment fitting equations
      {
        FOUR_C_THROW("Not supported.");
      }
    }
  }
}

/*--------------------------------------------------------------------------------------------*
                 Generate boundary cells for facet cut by level set side

               COMMENT:  Might need to rethink BC-creation, as it generates a lot of tris now.
                         Could probably be enough with quads some times?
*---------------------------------------------------------------------------------------------*/
void Cut::VolumeCell::generate_boundary_cells_level_set_side(Mesh& mesh,
    const Cut::Point::PointPosition posi, Element* elem, Facet* fac, int BaseNos,
    Cut::BCellGaussPts BCellgausstype)
{
  if (not fac->belongs_to_level_set_side())
    FOUR_C_THROW("Why would you call BC-creation for LS-Side without a LS side?");

  if (BCellgausstype == Cut::BCellGaussPts_MomentFitting)
    FOUR_C_THROW("Not supported for BC-Cell creation for LevelSetSides.");

  // Is the facet split/triangulated and if it consists of 4 corners is it planar.
  //  Then decompose and create Boundary Cells from triangulation/split.
  //  Otherwise, create quad4/tri3 boundary cell.
  bool istriangulated = (fac->is_triangulated() or fac->is_facet_split());
  bool notsimpleshape =
      !(fac->corner_points().size() == 4 and (fac->is_planar(mesh, fac->corner_points())));

  // UNCOMMENT THIS FOR ONLY TRIANGLES ON SURFACE!
  // DO FULL TRIANGULATION OF BC-SURFACE
  //---------------------------------------
  //  if(fac->CornerPoints().size()>3)
  //  {
  ////    std::cout << "fac->CornerPoints().size()>3?!" << std::endl;
  //    if(not fac->IsTriangulated())
  //    {
  ////      std::cout << "not fac->IsTriangulated()" << std::endl;
  //      fac->DoTriangulation(mesh,fac->Points());  //triangulate everything!
  ////      std::cout << "fac->IsTriangulated(): " << fac->IsTriangulated() << std::endl;
  //    }
  //  }
  //  notsimpleshape = true; //Create BC from triangulation solely
  //  istriangulated = true;
  //---------------------------------------

  //  if( (fac->IsTriangulated() or fac->IsFacetSplit()) and (fac->CornerPoints().size() == 4 and
  //  !(fac->is_planar(mesh, fac->CornerPoints())))  )
  if (istriangulated and notsimpleshape)
  {
    std::vector<std::vector<Point*>> facet_triang;
    if (fac->is_triangulated())
      facet_triang = fac->triangulation();
    else
      facet_triang = fac->get_split_cells();

    for (std::vector<std::vector<Point*>>::const_iterator j = facet_triang.begin();
        j != facet_triang.end(); ++j)
    {
      std::vector<Point*> tri = *j;
      std::vector<Point*> tri_temp(tri);  // could be quad?

      //    // when finding eqn of plane for the facet, inline points should not be taken
      //    Cut::Kernel::DeleteInlinePts( cornersTemp );

      std::vector<double> fac_tri_normal(4);
      fac_tri_normal = Kernel::eqn_plane_of_polygon(tri_temp);

      Core::LinAlg::Matrix<3, 1> ls_coord(tri_temp[1]->x());
      const std::vector<double> fac_ls_normal =
          elem->get_level_set_gradient(ls_coord);  // fac->GetLevelSetFacetNormal(elem);
      double dotProduct = fac_tri_normal[0] * fac_ls_normal[0] +
                          fac_tri_normal[1] * fac_ls_normal[1] +
                          fac_tri_normal[2] * fac_ls_normal[2];
      if (posi == Cut::Point::outside)
      {
        if (dotProduct > 0.0)
        {
          std::reverse(tri_temp.begin(), tri_temp.end());
        }
      }
      else if (posi == Cut::Point::inside)
      {
        if (dotProduct < 0.0)  // ( < ) should be correct solution.
          std::reverse(tri_temp.begin(), tri_temp.end());
      }
      if (tri_temp.size() == 3)
      {
        double areaCell = Cut::Kernel::get_area_tri(tri_temp);
        if (areaCell < REF_AREA_BCELL)
        {
          std::cout << "BCell NOT ADDED! areaCell: " << areaCell << std::endl;
          continue;
        }
        new_tri3_cell(mesh, fac, tri_temp);
      }
      else if (tri_temp.size() == 4)
      {
        double areaCell = Cut::Kernel::get_area_convex_quad(tri_temp);
        if (areaCell < REF_AREA_BCELL)
        {
          std::cout << "BCell NOT ADDED! areaCell: " << areaCell << std::endl;
          continue;
        }
        new_quad4_cell(mesh, fac, tri_temp);
      }
      else
        FOUR_C_THROW("Triangulation created neither tri3 or quad4");
    }
  }
  else
  {
    // For non-triangulated facets-> i.e. planar quad4 or tri3. (does quad4 exist here?)

    std::vector<Point*> tri = fac->corner_points();
    std::vector<Point*> tri_temp(tri);  // could be quad?


    //    // when finding eqn of plane for the facet, inline points should not be taken
    //    Cut::Kernel::DeleteInlinePts( cornersTemp );

    // Fix normal direction
    std::vector<double> fac_tri_normal(4);
    fac_tri_normal = Kernel::eqn_plane_of_polygon(tri_temp);

    Core::LinAlg::Matrix<3, 1> ls_coord(tri_temp[1]->x());
    const std::vector<double> fac_ls_normal =
        elem->get_level_set_gradient(ls_coord);  // fac->GetLevelSetFacetNormal(elem);
    double dotProduct = fac_tri_normal[0] * fac_ls_normal[0] +
                        fac_tri_normal[1] * fac_ls_normal[1] + fac_tri_normal[2] * fac_ls_normal[2];
    if (posi == Cut::Point::outside)
    {
      if (dotProduct > 0.0)
      {
        std::reverse(tri_temp.begin(), tri_temp.end());
      }
    }
    else if (posi == Cut::Point::inside)
    {
      if (dotProduct < 0.0)  // ( < ) should be correct solution.
        std::reverse(tri_temp.begin(), tri_temp.end());
    }

    // Add boundary cell
    if (tri_temp.size() == 3)
    {
      double areaCell = Cut::Kernel::get_area_tri(tri_temp);
      if (areaCell < REF_AREA_BCELL)
      {
        std::cout << "BCell NOT ADDED! areaCell: " << areaCell << std::endl;
      }
      new_tri3_cell(mesh, fac, tri_temp);
    }
    else if (tri_temp.size() == 4)
    {
      double areaCell = Cut::Kernel::get_area_convex_quad(tri_temp);
      if (areaCell < REF_AREA_BCELL)
      {
        std::cout << "BCell NOT ADDED! areaCell: " << areaCell << std::endl;
      }
      new_quad4_cell(mesh, fac, tri_temp);
    }
    else
      FOUR_C_THROW("Triangulation created neither tri3 or quad4");
  }
}

/*--------------------------------------------------------------------------------------------------------*
    This is to check whether the corner points of the cut side facet is aligned to give outward
normal
*---------------------------------------------------------------------------------------------------------*/
bool Cut::VolumeCell::to_reverse(const Cut::Point::PointPosition posi,
    const std::vector<double>& parEqn, const std::vector<double>& facetEqn)
{
  bool reverse = false;

  // position is inside
  if (posi == Cut::Point::outside)  //-3 before...
  {
    if (fabs(parEqn[0]) > TOL_EQN_PLANE && parEqn[0] * facetEqn[0] > 0.0)
      reverse = true;
    else if (fabs(parEqn[1]) > TOL_EQN_PLANE && parEqn[1] * facetEqn[1] > 0.0)
      reverse = true;
    else if (fabs(parEqn[2]) > TOL_EQN_PLANE && parEqn[2] * facetEqn[2] > 0.0)
      reverse = true;
    else
      reverse = false;
  }

  // position is outside
  else if (posi == Cut::Point::inside)  //-2 before...
  {
    if (fabs(parEqn[0]) > TOL_EQN_PLANE && parEqn[0] * facetEqn[0] < 0.0)
      reverse = true;
    else if (fabs(parEqn[1]) > TOL_EQN_PLANE && parEqn[1] * facetEqn[1] < 0.0)
      reverse = true;
    else if (fabs(parEqn[2]) > TOL_EQN_PLANE && parEqn[2] * facetEqn[2] < 0.0)
      reverse = true;
    else
      reverse = false;
  }
  return reverse;
}

/*------------------------------------------------------------------------------------------*
   When DirectDivergence method is used for gauss point generation, for every gauss point
   on the facet, an internal gauss rule is to be generated to find the modified integrand
*-------------------------------------------------------------------------------------------*/
std::shared_ptr<Core::FE::GaussPoints> Cut::VolumeCell::generate_internal_gauss_rule(
    std::shared_ptr<Core::FE::GaussPoints>& gp)
{
  // TEUCHOS_FUNC_TIME_MONITOR( "Cut::VolumeCell::generate_internal_gauss_rule" );


  Core::FE::GaussIntegration grule(gp);

  std::shared_ptr<Core::FE::CollectedGaussPoints> cgp =
      std::make_shared<Core::FE::CollectedGaussPoints>(0);

  for (Core::FE::GaussIntegration::iterator quadint = grule.begin(); quadint != grule.end();
      ++quadint)
  {
    const Core::LinAlg::Matrix<3, 1> etaFacet(
        quadint.point());  // coordinates and weight of main gauss point
    Core::LinAlg::Matrix<3, 1> intpt(etaFacet);

    Core::FE::GaussIntegration gi(Core::FE::CellType::line2,
        (DIRECTDIV_GAUSSRULE - 1));  // internal gauss rule for interval (-1,1)

    // x-coordinate of main Gauss point is projected in the reference plane
    double xbegin = (ref_eqn_plane_[3] - ref_eqn_plane_[1] * etaFacet(1, 0) -
                        ref_eqn_plane_[2] * etaFacet(2, 0)) /
                    ref_eqn_plane_[0];

    double jac = fabs(xbegin - etaFacet(0, 0)) * 0.5;  // jacobian for 1D transformation rule

    // -----------------------------------------------------------------------------
    // project internal gauss point from interval (-1,1) to the actual interval
    // -----------------------------------------------------------------------------
    for (Core::FE::GaussIntegration::iterator iqu = gi.begin(); iqu != gi.end(); ++iqu)
    {
      const Core::LinAlg::Matrix<1, 1> eta(iqu.point());
      double weight = iqu.weight();

      double xmid = 0.5 * (xbegin + etaFacet(0, 0));
      intpt(0, 0) = (xmid - xbegin) * eta(0, 0) + xmid;  // location of internal gauss points

      weight = weight * jac;  // weight of internal gauss points
      if (xbegin > etaFacet(0, 0)) weight = -1.0 * weight;

      weight =
          weight * quadint.weight();  // multiply with weight of main Gauss points so that internal
                                      // and main pts can be combined into a single data structure
      cgp->append(intpt, weight);
    }
  }
  return cgp;
}

/*------------------------------------------------------------------------------------------*
   Moment fitting equations are solved at each volume cell to construct integration rules
*-------------------------------------------------------------------------------------------*/
void Cut::VolumeCell::moment_fit_gauss_weights(
    Element* elem, Mesh& mesh, bool include_inner, Cut::BCellGaussPts BCellgausstype)
{
#ifdef LOCAL
  // position is used to decide whether the ordering of points are in clockwise or not
  if (Position() == Point::undecided)
  {
    if (!set_position_cut_side_based())
    {
      FOUR_C_THROW("undefined position for the volumecell");
    }
  }

  // if the volumecell is inside and include_inner is false, no need to compute the Gaussian points
  // as this vc will never be computed in xfem algorithm
  if (Position() == Point::inside && include_inner == false) return;

  int BaseNos = 84;  // number of base functions to be used in the integration
  VolumeIntegration vc_inte(this, elem, Position(), BaseNos);

  weights_ = vc_inte.compute_weights();            // obtain the integration weight at all points
  gaussPts_ = vc_inte.get_gauss_point_location();  // get the coordinates of all the Gauss points

  gp_ = gauss_points_fitting();  // convert the weight and the location to Gauss rule

  // generate boundary cells -- when using tessellation this is automatically done
  generate_boundary_cells(mesh, Position(), elem, BaseNos, BCellgausstype);

  // std::cout<<"MOMENT FITTING ::: Number of points = "<<weights_.Length()<<"\n";
#else

  // std::cout << "DirectDivergence Is Used!!!! MomFitting not functional in Global coordinates. "
  // << std::endl;

  std::cout << "MomFitting not functional in Global coordinates. NOTHING IS DONE!" << std::endl;

  //  direct_divergence_gauss_rule( elem, mesh, include_inner, BCellgausstype );
#endif
}

/*---------------------------------------------------------------------------------------------------------------*
                     The facets that have non-zero x-component normal is triangulated. sudhakar
03/12 The gauss integration rules are generated by applying divergence theorem The reference facet
is identified which will be used to find the modified integral in fluid integration
*----------------------------------------------------------------------------------------------------------------*/
void Cut::VolumeCell::direct_divergence_gauss_rule(
    Element* elem, Mesh& mesh, bool include_inner, Cut::BCellGaussPts BCellgausstype)
{
  if (elem->shape() != Core::FE::CellType::hex8 && elem->shape() != Core::FE::CellType::hex20)
    FOUR_C_THROW("direct_divergence_gauss_rule: Just hex8 and hex20 available yet in DD!");

  if (BCellgausstype != Cut::BCellGaussPts_Tessellation)
    FOUR_C_THROW(
        "direct_divergence_gauss_rule: just Cut::BCellGaussPts_Tessellation supported "
        "at "
        "the "
        "moment!");

  // position is used to decide whether the ordering of points are in clockwise or not
  if (position() == Point::undecided)
  {
    if (!set_position_cut_side_based())
    {
      FOUR_C_THROW("undefined position for the volumecell");
    }
  }

  // if the volumecell is inside and includeinner is false, no need to compute the Gaussian points
  // as this vc will never be computed in xfem algorithm
  if (position() == Point::inside and include_inner == false) return;

  // If the Volume Cell consists of less than 4 facets, it can't span a volume in 3D.
  if (facets().size() < 4)
    FOUR_C_THROW(
        "If the Volume Cell consists of less than 4 facets, it can't span a volume in 3D?");
  // return;

  is_negligible_small_ = false;


  DirectDivergence dd(this, elem, position(), mesh);

  ref_eqn_plane_.reserve(4);  // it has to store a,b,c,d in ax+by+cz=d

  std::shared_ptr<Core::FE::GaussPoints> gp =
      dd.vc_integration_rule(ref_eqn_plane_);  // compute main gauss points
  gp_ =
      generate_internal_gauss_rule(gp);  // compute internal gauss points for every main gauss point

  // compute volume of this cell
  // also check whether generated gauss rule predicts volume accurately
  // also check this vc can be eliminated due to its very small volume
  bool isNegVol = false;
  {
    Core::FE::GaussIntegration gpi(gp_);
    dd.debug_volume(gpi, isNegVol);

    // then this vol is extremely small that we erase the gauss points
    if (isNegVol)
    {
      gp_.reset();
      is_negligible_small_ = true;
    }
  }

  if (not isNegVol)
  {
#ifdef LOCAL

#else
    // we have generated the integration rule in global coordinates of the element
    // Now we map this rule to local coordinates since the weak form evaluation is done on local
    // coord
    project_gauss_points_to_local_coordinates();
#endif
  }
  // generate boundary cells -- when using tessellation this is automatically done
  generate_boundary_cells(mesh, position(), elem, 0, BCellgausstype);
}

/*----------------------------------------------------------------------------------------------------*
 * Project the integration rule generated on global coordinate system of sudhakar 05/15 the
 *background element to its local coordinates
 *----------------------------------------------------------------------------------------------------*/
void Cut::VolumeCell::project_gauss_points_to_local_coordinates()
{
  if (element_->shape() != Core::FE::CellType::hex8)
    FOUR_C_THROW(
        "Currently Direct divergence in global coordinates works only for hex8 elements\n");

  Core::FE::GaussIntegration intpoints(gp_);

  if (element_->is_shadow() && (element_->get_quad_shape() == Core::FE::CellType::hex20 ||
                                   element_->get_quad_shape() == Core::FE::CellType::hex27))
  {
    switch (element_->get_quad_shape())
    {
      case Core::FE::CellType::hex20:
      {
        Core::LinAlg::Matrix<3, 20> xyze;
        element_->coordinates_quad(xyze.data());
        gp_ = Core::FE::GaussIntegration::project_gauss_points_global_to_local<
            Core::FE::CellType::hex20>(xyze, intpoints, false);
        break;
      }
      case Core::FE::CellType::hex27:
      {
        Core::LinAlg::Matrix<3, 27> xyze;
        element_->coordinates_quad(xyze.data());
        gp_ = Core::FE::GaussIntegration::project_gauss_points_global_to_local<
            Core::FE::CellType::hex27>(xyze, intpoints, false);
        break;
      }
      default:
      {
        FOUR_C_THROW(
            "Currently Direct divergence in global coordinates works only for hex8, hex20, hex27 "
            "elements\n");
        break;
      }
    }
  }
  else
  {
    Core::LinAlg::Matrix<3, 8> xyze;
    element_->coordinates(xyze.data());
    gp_ =
        Core::FE::GaussIntegration::project_gauss_points_global_to_local<Core::FE::CellType::hex8>(
            xyze, intpoints, false);
  }
}

/*-------------------------------------------------------------------------------------*
| Return Ids of all the points associated with this volumecell           shahmiri 06/12
*--------------------------------------------------------------------------------------*/
const std::set<int>& Cut::VolumeCell::volume_cell_point_ids()
{
  if (vcpoints_ids_.size() != 0)
  {
    return vcpoints_ids_;
  }
  else
  {
    const plain_facet_set& facete = facets();

    // loop over facets
    for (plain_facet_set::const_iterator i = facete.begin(); i != facete.end(); i++)
    {
      Facet* fe = *i;
      const std::vector<Point*>& corners = fe->corner_points();

      for (std::vector<Point*>::const_iterator c = corners.begin(); c != corners.end(); c++)
      {
        Point* pt = *c;
        vcpoints_ids_.insert(pt->id());
      }
    }
  }

  if (vcpoints_ids_.size() == 0) FOUR_C_THROW("The size of volumecell points is zero!!");

  return vcpoints_ids_;
}

/*-------------------------------------------------------------------------------------*
| Find Position of the Volumecell based on the orientation of the cut_sides   ager 08/15
*--------------------------------------------------------------------------------------*/
bool Cut::VolumeCell::set_position_cut_side_based()
{
  if (position() != Point::undecided)
    FOUR_C_THROW(
        "Do not call FindPositionCutSideBased() if Position for volumecell is already set ({})!",
        (int)position());

  std::map<Facet*, bool> outsidenormal;

  const plain_facet_set& facets = VolumeCell::facets();

  // First find a facet based on a background side!
  Facet* f;
  for (plain_facet_set::const_iterator i = facets.begin(); i != facets.end(); ++i)
  {
    f = *i;
    if (!f->on_cut_side())
    {
      //-----
      // STEP 1: Get centroid of the parent element
      //-----
      Core::LinAlg::Matrix<3, 1> elecen;
      parent_element()->element_center(elecen);

      //-----
      // STEP 2: Get geometric center point of the facet
      // For concave facets, this is not the actual center, but it does not matter
      //-----
      Core::LinAlg::Matrix<3, 1> facecen;
      for (std::vector<Point*>::const_iterator fit = f->corner_points().begin();
          fit != f->corner_points().end(); fit++)
      {
        for (unsigned dim = 0; dim < 3; dim++) facecen(dim, 0) += (*fit)->x()[dim];
      }

      for (unsigned dim = 0; dim < 3; dim++)
        facecen(dim, 0) = facecen(dim, 0) / f->corner_points().size();

      //-----
      // STEP 3: Construct a vector that points FROM element centre TO the facet centre
      // This reference vector points outside the background element! (for LOCAL this is a unit
      // normal vector)
      //-----
      Core::LinAlg::Matrix<3, 1> ref_vec;
      ref_vec.update(1.0, facecen, -1.0, elecen);

      //-----
      // STEP 4: get point from parent side as the normal orientation should be calculated from the
      // side!!!
      // get unit normal vector
      std::vector<double> eqn_plane = Kernel::eqn_plane_of_polygon(f->corner_points());


      //-----
      // STEP 5: Take dot product with the normal of facet
      // If both are in the opposite direction, then the facet nodes are arranged clockwise
      //-----
      Core::LinAlg::Matrix<3, 1> norm_fac;
      for (unsigned dim = 0; dim < 3; dim++) norm_fac(dim, 0) = eqn_plane[dim];

      double dotProduct =
          ref_vec.dot(norm_fac);  // > 0 ...Side Normal is pointing outside of the Element //< 0
                                  // Side Normal points inside the element!
      if (abs(dotProduct) < BASICTOL)
      {
        std::cout << "Reference Vector: (" << ref_vec(0, 0) << "," << ref_vec(1, 0) << ","
                  << ref_vec(2, 0) << ")" << std::endl;
        std::cout << "Facet Vector: (" << norm_fac(0, 0) << "," << norm_fac(1, 0) << ","
                  << norm_fac(2, 0) << ")" << std::endl;
        FOUR_C_THROW("Check this really small dotProduct! {}", dotProduct);
      }

      outsidenormal[f] = (dotProduct > 0);
    }  // if not cutfacet
  }  // for facets

  int iter = 0;
  bool done = false;
  while (outsidenormal.size() != facets.size() && iter < 1000)
  {
    for (plain_facet_set::const_iterator i = facets.begin(); i != facets.end(); ++i)
    {
      Facet* ff = *i;
      if (ff->on_cut_side() && outsidenormal.find(ff) == outsidenormal.end())  // cutside
      {
        for (std::map<Facet*, bool>::iterator on = outsidenormal.begin(); on != outsidenormal.end();
            ++on)
        {
          bool consistent_normal = false;
          if (ff->have_consistent_normal(on->first, consistent_normal))
          {
            if (consistent_normal)
              outsidenormal[ff] = outsidenormal[on->first];
            else
              outsidenormal[ff] = !outsidenormal[on->first];
            done = true;  // If we can a least identify the direction on one CutSide facet this is
                          // fine in principle
            break;
          }
        }
      }
    }
    iter++;
  }

  if (iter == 1000 && !done)
  {
    FOUR_C_THROW(
        "set_position_cut_side_based failed: too many iterations (theoretically a facet with many "
        "points could also lead to this)!");
    return false;
  }

  Point::PointPosition posi = Point::undecided;
  for (std::map<Facet*, bool>::iterator on = outsidenormal.begin(); on != outsidenormal.end(); ++on)
  {
    if (on->first->on_cut_side())  // cutside
    {
      double prod = 0.0;
      std::vector<double> eqn_plane_facet =
          Kernel::eqn_plane_of_polygon(on->first->corner_points());
      std::vector<Point*> spoints(on->first->parent_side()->nodes().size());

      for (unsigned int i = 0; i < on->first->parent_side()->nodes().size(); ++i)
        spoints[i] = on->first->parent_side()->nodes()[i]->point();

      // get unit normal vector
      std::vector<double> eqn_plane_side = Kernel::eqn_plane_of_polygon(spoints);

      for (unsigned int dim = 0; dim < 3; ++dim) prod += eqn_plane_facet[dim] * eqn_plane_side[dim];

      if ((on->second && prod > 0) || (!on->second && prod < 0))  // this means that the
      {
        if (posi != Point::undecided && posi != Point::inside)
          FOUR_C_THROW(
              "set_position_cut_side_based: posi != Point::undecided && posi != Point::inside (Are "
              "all "
              "you Cut Sides oriented correct?)");
        // FOUR_C_THROW("set_position_cut_side_based: posi != Point::undecided && posi !=
        // Point::inside (Are all you Cut Sides oriented correct?)");
        posi = Point::inside;
      }
      else
      {
        if (posi != Point::undecided && posi != Point::outside)
          FOUR_C_THROW(
              "set_position_cut_side_based: posi != Point::undecided && posi != Point::outside "
              "(Are "
              "all you Cut Sides oriented correct?)");
        // FOUR_C_THROW("set_position_cut_side_based: posi != Point::undecided && posi !=
        // Point::outside (Are all you Cut Sides oriented correct?)");
        posi = Point::outside;
      }
    }
  }

  if (posi != Point::undecided)
  {
    position_ = posi;
    return true;
  }
  else
  {
    return false;
  }
}

FOUR_C_NAMESPACE_CLOSE
