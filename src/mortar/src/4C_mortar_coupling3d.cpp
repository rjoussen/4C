// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_mortar_coupling3d.hpp"

#include "4C_comm_mpi_utils.hpp"
#include "4C_fem_discretization.hpp"
#include "4C_linalg_serialdensevector.hpp"
#include "4C_linalg_utils_densematrix_inverse.hpp"
#include "4C_linalg_utils_densematrix_multiply.hpp"
#include "4C_mortar_defines.hpp"
#include "4C_mortar_integrator.hpp"
#include "4C_mortar_node.hpp"
#include "4C_mortar_projector.hpp"
#include "4C_mortar_utils.hpp"

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*
 |  ctor (public)                                             popp 06/09|
 *----------------------------------------------------------------------*/
Mortar::Coupling3d::Coupling3d(Core::FE::Discretization& idiscret, int dim, bool quad,
    Teuchos::ParameterList& params, Mortar::Element& sele, Mortar::Element& mele)
    : idiscret_(idiscret),
      dim_(dim),
      shapefcn_(Teuchos::getIntegralValue<Mortar::ShapeFcn>(params, "LM_SHAPEFCN")),
      quad_(quad),
      lmquadtype_(Teuchos::getIntegralValue<Mortar::LagMultQuad>(params, "LM_QUAD")),
      sele_(sele),
      mele_(mele),
      imortar_(params),
      lauxn_(-1.0)
{
  // empty constructor body
  return;
}


/*----------------------------------------------------------------------*
 |  get communicator  (public)                                popp 06/09|
 *----------------------------------------------------------------------*/
MPI_Comm Mortar::Coupling3d::get_comm() const { return idiscret_.get_comm(); }


/*----------------------------------------------------------------------*
 |  Evaluate coupling (3D)                                    popp 03/09|
 *----------------------------------------------------------------------*/
bool Mortar::Coupling3d::evaluate_coupling()
{
  // rough check whether element centers are "near"
  // whether or not quadratic 3d coupling is performed, we only
  // check the distance of the parent slave and master elements
  bool near = rough_check_centers();
  if (!near) return false;

  // tolerance for polygon clipping
  double tol = 0.0;

  // map to store projection parameter alpha for each master node
  // (currently only necessary for non-AUXPLANE case)
  std::map<int, double> projpar;

  // rough check of orientation of element centers
  // if slave and master element center normals form an
  // angle > 90 degrees the pair will not be considered further
  bool orient = rough_check_orient();
  if (!orient) return false;

  // project slave element nodes onto auxiliary plane
  project_slave();

  // project master element nodes onto auxiliary plane
  project_master();

  // tolerance for polygon clipping
  const double sminedge = source_int_element().min_edge_size();
  const double mminedge = target_int_element().min_edge_size();
  tol = MORTARCLIPTOL * std::min(sminedge, mminedge);

  // do polygon clipping
  clip() = Utils::clip_polygons(source_vertices(), target_vertices(), auxn(), tol);
  int clipsize = (int)(clip().size());

  // proceed only if clipping polygon is at least a triangle
  if (clipsize < 3) return false;

  // proceed only if clipping polygon has non-zero area
  if (polygon_area() < MORTARINTLIM * source_element_area()) return false;

  // check / set  projection status of slave nodes
  has_proj_status();

  // do triangulation (+linearization) of clip polygon
  triangulation(projpar, tol);

  return true;
}


/*----------------------------------------------------------------------*
 |  Rough check if elements are near (with centers)           popp 11/08|
 *----------------------------------------------------------------------*/
bool Mortar::Coupling3d::rough_check_centers()
{
  const double some = source_int_element().max_edge_size();
  const double mme = target_int_element().max_edge_size();
  const double near = 2.0 * std::max(some, mme);

  double loccs[2] = {0.0, 0.0};
  Core::FE::CellType dts = source_int_element().shape();
  if (dts == Core::FE::CellType::tri3 || dts == Core::FE::CellType::tri6)
  {
    loccs[0] = 1.0 / 3.0;
    loccs[1] = 1.0 / 3.0;
  }
  double loccm[2] = {0.0, 0.0};
  Core::FE::CellType dtm = target_int_element().shape();
  if (dtm == Core::FE::CellType::tri3 || dtm == Core::FE::CellType::tri6)
  {
    loccm[0] = 1.0 / 3.0;
    loccm[1] = 1.0 / 3.0;
  }

  double sc[3] = {0.0, 0.0, 0.0};
  double mc[3] = {0.0, 0.0, 0.0};
  source_int_element().local_to_global(loccs, sc, 0);
  target_int_element().local_to_global(loccm, mc, 0);

  const double cdist = sqrt((mc[0] - sc[0]) * (mc[0] - sc[0]) + (mc[1] - sc[1]) * (mc[1] - sc[1]) +
                            (mc[2] - sc[2]) * (mc[2] - sc[2]));
  if (cdist >= near)
    return false;
  else
    return true;
}


/*----------------------------------------------------------------------*
 |  Rough check if elements are near (with master nodes)      popp 11/09|
 *----------------------------------------------------------------------*/
bool Mortar::Coupling3d::rough_check_nodes()
{
  // project master nodes onto auxiliary plane
  int nnodes = target_int_element().num_node();
  Core::Nodes::Node** mynodes = target_int_element().nodes();
  if (!mynodes) FOUR_C_THROW("rough_check_nodes: Null pointer!");

  // prepare check
  bool near = false;
  const double some = source_int_element().max_edge_size();
  const double mme = target_int_element().max_edge_size();
  const double limit = 0.3 * std::max(some, mme);

  for (int i = 0; i < nnodes; ++i)
  {
    Node* mycnode = dynamic_cast<Node*>(mynodes[i]);
    if (!mycnode) FOUR_C_THROW("rough_check_nodes: Null pointer!");

    // first build difference of point and element center
    // and then dot product with unit normal at center
    double dist = (mycnode->xspatial()[0] - auxc()[0]) * auxn()[0] +
                  (mycnode->xspatial()[1] - auxc()[1]) * auxn()[1] +
                  (mycnode->xspatial()[2] - auxc()[2]) * auxn()[2];

    // only near if at least one node is "close" to auxiliary plane
    // (signed distance in order to capture intrusion correctly)
    if (dist < limit) near = true;
  }

  return near;
}


/*----------------------------------------------------------------------*
 |  Rough check if elements are near (with normals)           popp 11/10|
 *----------------------------------------------------------------------*/
bool Mortar::Coupling3d::rough_check_orient()
{
  // compute auxiliary plane for 3D coupling
  auxiliary_plane();

  // we first need the master element center:
  // for quad4, quad8, quad9 elements: xi = eta = 0.0
  // for tri3, tri6 elements: xi = eta = 1/3
  double loccenter[2] = {0.0, 0.0};

  Core::FE::CellType dt = target_int_element().shape();
  if (dt == Core::FE::CellType::tri3 || dt == Core::FE::CellType::tri6)
  {
    loccenter[0] = 1.0 / 3;
    loccenter[1] = 1.0 / 3;
  }
  else if (dt == Core::FE::CellType::quad4 || dt == Core::FE::CellType::quad8 ||
           dt == Core::FE::CellType::quad9)
  {
    loccenter[0] = 0.0;
    loccenter[1] = 0.0;
  }
  else
    FOUR_C_THROW("rough_check_orient called for unknown element type");

  // compute the unit normal vector at the master element center
  double nmc[3] = {0.0, 0.0, 0.0};
  target_int_element().compute_unit_normal_at_xi(loccenter, nmc);

  // check orientation with respect to slave element
  double dot = nmc[0] * auxn()[0] + nmc[1] * auxn()[1] + nmc[2] * auxn()[2];
  if (dot < -1.0e-12)
    return true;
  else
    return false;
}


/*----------------------------------------------------------------------*
 |  Build auxiliary plane from slave element (public)         popp 11/08|
 *----------------------------------------------------------------------*/
bool Mortar::Coupling3d::auxiliary_plane()
{
  // we first need the element center:
  // for quad4, quad8, quad9 elements: xi = eta = 0.0
  // for tri3, tri6 elements: xi = eta = 1/3
  double loccenter[2] = {0.0, 0.0};

  Core::FE::CellType dt = source_int_element().shape();
  if (dt == Core::FE::CellType::tri3 || dt == Core::FE::CellType::tri6)
  {
    loccenter[0] = 1.0 / 3.0;
    loccenter[1] = 1.0 / 3.0;
  }
  else if (dt == Core::FE::CellType::quad4 || dt == Core::FE::CellType::quad8 ||
           dt == Core::FE::CellType::quad9)
  {
    loccenter[0] = 0.0;
    loccenter[1] = 0.0;
  }
  else
    FOUR_C_THROW("auxiliary_plane called for unknown element type");

  // compute element center via shape fct. interpolation
  source_int_element().local_to_global(loccenter, auxc(), 0);

  // we then compute the unit normal vector at the element center
  lauxn() = source_int_element().compute_unit_normal_at_xi(loccenter, auxn());


  // calculate auxplane with cpp normal!
  //  Core::LinAlg::SerialDenseVector val(SlaveIntElement().num_node());
  //  Core::LinAlg::SerialDenseMatrix deriv(SlaveIntElement().num_node(),2,true);
  //  SlaveIntElement().evaluate_shape(loccenter, val, deriv, SlaveIntElement().num_node(),false);
  //  Auxn()[0] = 0.0;
  //  Auxn()[1] = 0.0;
  //  Auxn()[2] = 0.0;
  //
  //  // interpolate between nodal normals
  //  for(int i=0;i<SlaveIntElement().num_node();++i)
  //  {
  //    Node* snode = dynamic_cast<Node*>(SlaveIntElement().Nodes()[i]);
  //    Auxn()[0] += val[i]*snode->MoData().n()[0];
  //    Auxn()[1] += val[i]*snode->MoData().n()[1];
  //    Auxn()[2] += val[i]*snode->MoData().n()[2];
  //  }
  //
  //  // get length (not used in the following)
  //  Lauxn() = sqrt(Auxn()[0]*Auxn()[0] + Auxn()[1]*Auxn()[1] + Auxn()[2]*Auxn()[2]);
  //
  //  // create unit normal
  //  Auxn()[0] /= Lauxn();
  //  Auxn()[1] /= Lauxn();
  //  Auxn()[2] /= Lauxn();

  // std::cout << "Slave Element: " << SlaveIntElement().Id() << std::endl;
  // std::cout << "->Center: " << Auxc()[0] << " " << Auxc()[1] << " " << Auxc()[2] << std::endl;
  //  std::cout << "->Normal cpp: " << Auxn()[0] << " " << Auxn()[1] << " " << Auxn()[2] <<
  //  std::endl;

  return true;
}


/*----------------------------------------------------------------------*
 |  Project slave element onto auxiliary plane (public)       popp 11/08|
 *----------------------------------------------------------------------*/
bool Mortar::Coupling3d::project_slave()
{
  // project slave nodes onto auxiliary plane
  int nnodes = source_int_element().num_node();
  Core::Nodes::Node** mynodes = source_int_element().nodes();
  if (!mynodes) FOUR_C_THROW("project_slave: Null pointer!");

  // initialize storage for slave coords + their ids
  std::vector<double> vertices(3);
  std::vector<int> snodeids(1);

  for (int i = 0; i < nnodes; ++i)
  {
    Node* mycnode = dynamic_cast<Node*>(mynodes[i]);
    if (!mycnode) FOUR_C_THROW("project_slave: Null pointer!");

    // first build difference of point and element center
    // and then dot product with unit normal at center
    const double dist = (mycnode->xspatial()[0] - auxc()[0]) * auxn()[0] +
                        (mycnode->xspatial()[1] - auxc()[1]) * auxn()[1] +
                        (mycnode->xspatial()[2] - auxc()[2]) * auxn()[2];

    // compute projection
    for (int k = 0; k < 3; ++k) vertices[k] = mycnode->xspatial()[k] - dist * auxn()[k];

    // get node id, too
    snodeids[0] = mycnode->id();

    // store into vertex data structure
    source_vertices().emplace_back(vertices, Vertex::source, snodeids, nullptr);

    // std::cout << "->RealNode(S) " << mycnode->Id() << ": " << mycnode->xspatial()[0] << " " <<
    // mycnode->xspatial()[1] << " " << mycnode->xspatial()[2] << std::endl; std::cout <<
    // "->ProjNode(S) " << mycnode->Id() << ": " << vertices[0] << " " << vertices[1] << " " <<
    // vertices[2] << std::endl;
  }

  return true;
}


/*----------------------------------------------------------------------*
 |  Project master element onto auxiliary plane (public)      popp 11/08|
 *----------------------------------------------------------------------*/
bool Mortar::Coupling3d::project_master()
{
  // project master nodes onto auxiliary plane
  int nnodes = target_int_element().num_node();
  Core::Nodes::Node** mynodes = target_int_element().nodes();
  if (!mynodes) FOUR_C_THROW("project_master: Null pointer!");

  // initialize storage for master coords + their ids
  std::vector<double> vertices(3);
  std::vector<int> mnodeids(1);

  for (int i = 0; i < nnodes; ++i)
  {
    Node* mycnode = dynamic_cast<Node*>(mynodes[i]);
    if (!mycnode) FOUR_C_THROW("project_master: Null pointer!");

    // first build difference of point and element center
    // and then dot product with unit normal at center
    const double dist = (mycnode->xspatial()[0] - auxc()[0]) * auxn()[0] +
                        (mycnode->xspatial()[1] - auxc()[1]) * auxn()[1] +
                        (mycnode->xspatial()[2] - auxc()[2]) * auxn()[2];

    // compute projection
    for (int k = 0; k < 3; ++k) vertices[k] = mycnode->xspatial()[k] - dist * auxn()[k];

    // get node id, too
    mnodeids[0] = mycnode->id();

    // store into vertex data structure
    target_vertices().emplace_back(vertices, Vertex::projtarget, mnodeids, nullptr);

    // std::cout << "->RealNode(M) " << mycnode->Id() << ": " << mycnode->xspatial()[0] << " " <<
    // mycnode->xspatial()[1] << " " << mycnode->xspatial()[2] << std::endl; std::cout <<
    // "->ProjNode(M) " << mycnode->Id() << ": " << vertices[0] << " " << vertices[1] << " " <<
    // vertices[2] << std::endl;
  }

  return true;
}


/*----------------------------------------------------------------------*
 |  Compute and return area of clipping polygon (3D)          popp 11/08|
 *----------------------------------------------------------------------*/
double Mortar::Coupling3d::polygon_area()
{
  // initialize and check for trivial case
  double area = 0.0;
  int clipsize = (int)(clip().size());
  if (clipsize < 3) return area;

  // loop over all vertices to compute area
  for (int k = 0; k < clipsize; ++k)
  {
    // current and next vertex vector
    std::vector<double> vk = clip()[k].coord();
    std::vector<double> vkp1;
    if (k == clipsize - 1)
      vkp1 = clip()[0].coord();
    else
      vkp1 = clip()[k + 1].coord();

    // cross product
    std::array<double, 3> cross = {0.0, 0.0, 0.0};
    cross[0] = vk[1] * vkp1[2] - vk[2] * vkp1[1];
    cross[1] = vk[2] * vkp1[0] - vk[0] * vkp1[2];
    cross[2] = vk[0] * vkp1[1] - vk[1] * vkp1[0];

    // add scalar product with negative(!) normal to area
    // (because the clip polygon is oriented clockwise w.r.t. normal)
    area -= 0.5 * (cross[0] * auxn()[0] + cross[1] * auxn()[1] + cross[2] * auxn()[2]);
  }

  // when areas are close to zero, we might get sign problems
  // (as we do not want to check negative / positive orientation
  // here anyway, but only close to zero areas, we simply return
  // the absolute value!)
  return abs(area);
}


/*----------------------------------------------------------------------*
 |  Compute and return area of slave element (3D)             popp 11/08|
 *----------------------------------------------------------------------*/
double Mortar::Coupling3d::source_element_area() const
{
  // initialize
  double selearea = source_int_element().mo_data().area();

  return selearea;
}


/*----------------------------------------------------------------------*
 |  Check /set projection status of slave nodes (3D)          popp 11/08|
 *----------------------------------------------------------------------*/
bool Mortar::Coupling3d::has_proj_status()
{
  // check all nodes
  int nnodes = source_int_element().num_node();
  Core::Nodes::Node** mynodes = source_int_element().nodes();
  if (!mynodes) FOUR_C_THROW("has_proj_status: Null pointer!");

  // loop over all slave nodes
  for (int i = 0; i < nnodes; ++i)
  {
    Node* mycnode = dynamic_cast<Node*>(mynodes[i]);
    if (!mycnode) FOUR_C_THROW("has_proj_status: Null pointer!");

    // loop over all vertices of clip polygon
    for (int j = 0; j < (int)(clip().size()); ++j)
    {
      bool identical = false;

      // check if this clip vertex is slave-type and has the
      // current slave node id
      if ((int)(clip()[j].v_type()) == Vertex::source)
        if (mycnode->id() == clip()[j].nodeids()[0]) identical = true;

      // set hasproj to true, if so
      if (identical) mycnode->has_proj() = true;
    }
  }

  return true;
}


/*----------------------------------------------------------------------*
 |  Triangulation of clip polygon (3D)                        popp 11/08|
 *----------------------------------------------------------------------*/
bool Mortar::Coupling3d::triangulation(std::map<int, double>& projpar, double tol)
{
  // number of nodes
  const int nsrows = source_element().num_node();
  const int nmrows = target_element().num_node();

  // preparations
  int clipsize = (int)(clip().size());
  std::vector<std::vector<Core::Gen::Pairedvector<int, double>>> linvertex(
      clipsize, std::vector<Core::Gen::Pairedvector<int, double>>(3, 3 * nsrows + 3 * nmrows));

  // get integration type
  auto tri_type = Teuchos::getIntegralValue<Mortar::Triangulation>(imortar_, "TRIANGULATION");

  //**********************************************************************
  // (1) Linearization of clip vertex coordinates (only for contact)
  //**********************************************************************
  vertex_linearization(linvertex, projpar);

  switch (tri_type)
  {
    //**********************************************************************
    // (2) Triangulation of clip polygon (DELAUNAY-based (new))
    //**********************************************************************
    case Mortar::triangulation_delaunay:
      if (!delaunay_triangulation(linvertex, tol))
        // (3) Backup triangulation of clip polygon (CENTER-based (old))
        center_triangulation(linvertex, tol);
      break;

    //**********************************************************************
    // (3) Backup triangulation of clip polygon (CENTER-based (old))
    //**********************************************************************
    case Mortar::triangulation_center:
      center_triangulation(linvertex, tol);
      break;

    default:
      FOUR_C_THROW("unknown triangulation type");
      break;
  }
  return true;
}

/*----------------------------------------------------------------------*
 |  Triangulation of clip polygon (3D) - DELAUNAY             popp 08/11|
 *----------------------------------------------------------------------*/
bool Mortar::Coupling3d::delaunay_triangulation(
    std::vector<std::vector<Core::Gen::Pairedvector<int, double>>>& linvertex, double tol)
{
  // preparations
  cells().resize(0);
  int clipsize = (int)(clip().size());
  bool out = false;
  if (out) std::cout << "***\nPolygon with " << clipsize << " vertices\n***" << std::endl;

  //**********************************************************************
  // (1) Trivial clipping polygon -> IntCells
  //**********************************************************************
  // clip polygon = triangle
  // no triangulation necessary -> 1 IntCell
  if (clipsize == 3)
  {
    // IntCell vertices = clip polygon vertices
    Core::LinAlg::Matrix<3, 3> coords;
    for (int i = 0; i < clipsize; ++i)
      for (int k = 0; k < 3; ++k) coords(k, i) = clip()[i].coord()[k];

    // create IntCell object and push back
    cells().push_back(std::make_shared<IntCell>(0, 3, coords, auxn(), Core::FE::CellType::tri3,
        linvertex[0], linvertex[1], linvertex[2], get_deriv_auxn()));

    // get out of here
    return true;
  }

  //**********************************************************************
  // (2) General clipping polygon: Triangulation -> IntCells
  //**********************************************************************
  // store Delaunay triangles here
  std::vector<std::vector<int>> triangles;

  // start with first triangle v0,v1,v2
  std::vector<int> currtriangle(3);
  currtriangle[0] = 0;
  currtriangle[1] = 1;
  currtriangle[2] = 2;
  triangles.push_back(currtriangle);

  // build Delaunay triangulation recursively (starting from a triangle
  // and then adding all remaining nodes of the clipping polygon 1-by-1)
  // loop over clip vertices v3,..,vN
  for (int c = 3; c < clipsize; ++c)
  {
    // current size of triangulated polygon
    int currsize = c + 1;

    // add next triangle v(c-1),vc,v0
    std::vector<int> nexttriangle(3);
    nexttriangle[0] = c - 1;
    nexttriangle[1] = c;
    nexttriangle[2] = 0;
    triangles.push_back(nexttriangle);

    if (out)
    {
      std::cout << "-> checking " << currsize << "th node of polygon" << std::endl;
      std::cout << "-> found " << (int)triangles.size() << " triangles" << std::endl;
      for (int k = 0; k < (int)triangles.size(); ++k)
        std::cout << "      Triangle: " << triangles[k][0] << " " << triangles[k][1] << " "
                  << triangles[k][2] << " " << std::endl;
    }

    // check Delaunay criterion for all triangles and sort
    // triangles accordingly (good / bad)
    int numt = (int)triangles.size();
    std::vector<bool> bad(numt);
    std::vector<double> close(numt);
    for (int t = 0; t < numt; ++t)
    {
      // initialize values indicating a close decision
      // these are needed to later introduce some tolerance into
      // the Delaunay criterion decision (which is needed in the
      // cases where this decision becomes non-unique / arbitrary).
      close[t] = 1.0e12;

      // indices of current triangle
      int idx0 = triangles[t][0];
      int idx1 = triangles[t][1];
      int idx2 = triangles[t][2];

      // coordinates of current triangle
      Core::LinAlg::Matrix<3, 3> coords;
      for (int k = 0; k < 3; ++k)
      {
        coords(k, 0) = clip()[idx0].coord()[k];
        coords(k, 1) = clip()[idx1].coord()[k];
        coords(k, 2) = clip()[idx2].coord()[k];
      }

      // define center of circumcircle of current triangle
      const double x1 = coords(0, 0);
      const double y1 = coords(1, 0);
      const double z1 = coords(2, 0);
      const double x2 = coords(0, 1);
      const double y2 = coords(1, 1);
      const double z2 = coords(2, 1);
      const double x3 = coords(0, 2);
      const double y3 = coords(1, 2);
      const double z3 = coords(2, 2);

      // a=vector P1P2, b=vector P2P3
      const double a1 = x2 - x1;
      const double a2 = y2 - y1;
      const double a3 = z2 - z1;
      const double b1 = x3 - x2;
      const double b2 = y3 - y2;
      const double b3 = z3 - z2;

      // normal vector of plane P1P2P3 via cross product
      const double no1 = a2 * b3 - b2 * a3;
      const double no2 = a3 * b1 - b3 * a1;
      const double no3 = a1 * b2 - b1 * a2;

      // perpendicular bisector of P1P2 via cross product
      const double c1 = a2 * no3 - no2 * a3;
      const double c2 = a3 * no1 - no3 * a1;
      const double c3 = a1 * no2 - no1 * a2;

      // perpendicular bisector of P2P3 via cross product
      const double d1 = b2 * no3 - no2 * b3;
      const double d2 = b3 * no1 - no3 * b1;
      const double d3 = b1 * no2 - no1 * b2;

      // mid-points of P1P2 and P2P3
      const double m1 = (x1 + x2) * 0.5;
      const double m2 = (y1 + y2) * 0.5;
      const double m3 = (z1 + z2) * 0.5;
      const double n1 = (x2 + x3) * 0.5;
      const double n2 = (y2 + y3) * 0.5;
      const double n3 = (z2 + z3) * 0.5;

      // try to minimize error
      int direction = 0;
      if (abs(auxn()[0]) >= abs(auxn()[1]) && abs(auxn()[0]) >= abs(auxn()[2])) direction = 1;
      if (abs(auxn()[1]) >= abs(auxn()[0]) && abs(auxn()[1]) >= abs(auxn()[2])) direction = 2;
      if (abs(auxn()[2]) >= abs(auxn()[0]) && abs(auxn()[2]) >= abs(auxn()[1])) direction = 3;
      if (direction == 0) FOUR_C_THROW("Did not find best direction");

      // intersection of the two perpendicular bisections
      // (solution of m1+s*c1 = n1+t*d1 and m2+s*c2 = n2+t*d2)
      double s = 0.0;
      if (direction == 1)
      {
        // minimize error in yz-plane by solving
        // m2+s*c2 = n2+t*d2 and m3+s*c3 = n3+t*d3
        s = (m3 * d2 - n3 * d2 - d3 * m2 + d3 * n2) / (c2 * d3 - c3 * d2);
      }
      else if (direction == 2)
      {
        // minimize error in xz-plane by solving
        // m1+s*c1 = n1+t*d1 and m3+s*c3 = n3+t*d3
        s = (m3 * d1 - n3 * d1 - d3 * m1 + d3 * n1) / (c1 * d3 - c3 * d1);
      }
      else /* (direction==3)*/
      {
        // minimize error in xy-plane by solving
        // m1+s*c1 = n1+t*d1 and m2+s*c2 = n2+t*d2
        s = (m2 * d1 - n2 * d1 - d2 * m1 + d2 * n1) / (c1 * d2 - c2 * d1);
      }

      // center of the circumcircle
      const double xcenter = m1 + s * c1;
      const double ycenter = m2 + s * c2;
      const double zcenter = m3 + s * c3;

      // radius of the circumcircle
      const double radius1 =
          sqrt((xcenter - x1) * (xcenter - x1) + (ycenter - y1) * (ycenter - y1) +
               (zcenter - z1) * (zcenter - z1));
      const double radius2 =
          sqrt((xcenter - x2) * (xcenter - x2) + (ycenter - y2) * (ycenter - y2) +
               (zcenter - z2) * (zcenter - z2));
      const double radius3 =
          sqrt((xcenter - x3) * (xcenter - x3) + (ycenter - y3) * (ycenter - y3) +
               (zcenter - z3) * (zcenter - z3));

      // check radius computation
      if (abs(radius2 - radius1) > tol || abs(radius3 - radius1) > tol)
      {
        std::cout << "***WARNING*** Delaunay triangulation failed (no well-defined circumcircles)"
                  << " -> using backup" << std::endl;

        // if Delaunay triangulation failed, use old center-based
        // triangulation as backup (therefore return false)
        return false;
      }

      // check Delaunay criterion for all other vertices
      // (of current polygon, NOT the full clipping polygon)
      for (int k = 0; k < currsize; ++k)
      {
        // no check needed for triangle vertices
        if (k == idx0 || k == idx1 || k == idx2) continue;

        // compute distance
        const double dist =
            sqrt((xcenter - clip()[k].coord()[0]) * (xcenter - clip()[k].coord()[0]) +
                 (ycenter - clip()[k].coord()[1]) * (ycenter - clip()[k].coord()[1]) +
                 (zcenter - clip()[k].coord()[2]) * (zcenter - clip()[k].coord()[2]));

        // monitor critical Delaunay criterion decision
        // (necessary to avoid inconsistent good/bad grouping later)
        const double diff = abs(dist - radius1);
        if (diff < close[t]) close[t] = diff;

        // check for bad triangle (without tolerance)
        if (dist < radius1) bad[t] = true;
      }
    }

    // make good/bad decision consistent (with tolerance)
    // (problems might occur if more than 3 vertices on circumcircle)
    for (int t = 0; t < numt; ++t)
    {
      // check if this good decision was really close
      if (!bad[t] && close[t] < tol)
      {
        // check if any bad decision was really close, too
        bool foundpartner = false;
        for (int u = 0; u < numt; ++u)
        {
          if (bad[u] && close[u] < tol) foundpartner = true;
        }

        // set good->bad if partner found
        if (foundpartner) bad[t] = true;
      }
    }

    // now we build vector of all good / bad triangles
    std::vector<std::vector<int>> goodtriangles;
    std::vector<std::vector<int>> badtriangles;
    for (int t = 0; t < numt; ++t)
    {
      if (bad[t])
        badtriangles.push_back(triangles[t]);
      else
        goodtriangles.push_back(triangles[t]);
    }

    if (out)
    {
      std::cout << "   -> found " << (int)goodtriangles.size() << " GOOD triangles" << std::endl;
      for (int k = 0; k < (int)goodtriangles.size(); ++k)
        std::cout << "      Triangle: " << goodtriangles[k][0] << " " << goodtriangles[k][1] << " "
                  << goodtriangles[k][2] << " " << std::endl;
      std::cout << "   -> found " << (int)badtriangles.size() << " BAD triangles" << std::endl;
      for (int k = 0; k < (int)badtriangles.size(); ++k)
        std::cout << "      Triangle: " << badtriangles[k][0] << " " << badtriangles[k][1] << " "
                  << badtriangles[k][2] << " " << std::endl;
    }

    // find vertices in bad triangles: ALL vertices
    // find vertices in bad triangles: NOT connected with current vertex
    std::vector<int> badv;
    std::vector<int> ncv;
    for (int t = 0; t < numt; ++t)
    {
      if (bad[t])
      {
        // indices of current bad triangle
        int idx0 = triangles[t][0];
        int idx1 = triangles[t][1];
        int idx2 = triangles[t][2];

        // collect ALL vertices
        bool foundbefore0 = false;
        for (int k = 0; k < (int)badv.size(); ++k)
        {
          if (badv[k] == idx0) foundbefore0 = true;
        }
        if (!foundbefore0) badv.push_back(idx0);

        bool foundbefore1 = false;
        for (int k = 0; k < (int)badv.size(); ++k)
        {
          if (badv[k] == idx1) foundbefore1 = true;
        }
        if (!foundbefore1) badv.push_back(idx1);

        bool foundbefore2 = false;
        for (int k = 0; k < (int)badv.size(); ++k)
        {
          if (badv[k] == idx2) foundbefore2 = true;
        }
        if (!foundbefore2) badv.push_back(idx2);

        // indices of current vertex neighbors
        int neighbor0 = c - 1;
        int neighbor1 = 0;

        // collect NOT connected vertices
        if (idx0 != c && idx0 != neighbor0 && idx0 != neighbor1)
        {
          bool foundbefore = false;
          for (int k = 0; k < (int)ncv.size(); ++k)
          {
            if (ncv[k] == idx0) foundbefore = true;
          }
          if (!foundbefore) ncv.push_back(idx0);
        }
        if (idx1 != c && idx1 != neighbor0 && idx1 != neighbor1)
        {
          bool foundbefore = false;
          for (int k = 0; k < (int)ncv.size(); ++k)
          {
            if (ncv[k] == idx1) foundbefore = true;
          }
          if (!foundbefore) ncv.push_back(idx1);
        }
        if (idx2 != c && idx2 != neighbor0 && idx2 != neighbor1)
        {
          bool foundbefore = false;
          for (int k = 0; k < (int)ncv.size(); ++k)
          {
            if (ncv[k] == idx2) foundbefore = true;
          }
          if (!foundbefore) ncv.push_back(idx2);
        }
      }
    }

    // build triangles formed by current vertex and ncv vertices
    std::vector<std::vector<int>> addtriangles;
    for (int k = 0; k < (int)ncv.size(); ++k)
    {
      // find ncv vertex neighbor0
      bool validneighbor0 = false;
      int off0 = 0;
      int neighbor0 = 0;
      do
      {
        // set neighbor
        neighbor0 = ncv[k] - 1 - off0;
        if ((ncv[k] - off0) == 0) neighbor0 = currsize - 1 - off0;

        // check if neighbor is in bad vertices
        for (int k = 0; k < (int)badv.size(); ++k)
        {
          if (badv[k] == neighbor0) validneighbor0 = true;
        }

        // increase counter
        ++off0;

      } while (!validneighbor0);

      // find ncv vertex neighbor1
      bool validneighbor1 = false;
      int off1 = 0;
      int neighbor1 = 0;
      do
      {
        // set neighbor
        neighbor1 = ncv[k] + 1 + off1;
        if ((ncv[k] + off1) == currsize - 1) neighbor1 = 0 + off1;

        // check if neighbor is in bad vertices
        for (int k = 0; k < (int)badv.size(); ++k)
        {
          if (badv[k] == neighbor1) validneighbor1 = true;
        }

        // increase counter
        ++off1;

      } while (!validneighbor1);

      // plausibility check
      if (neighbor0 == c || neighbor1 == c) FOUR_C_THROW("Connected nodes not possible here");

      // add triangles
      std::vector<int> add1(3);
      add1[0] = c;
      add1[1] = ncv[k];
      add1[2] = neighbor0;
      addtriangles.push_back(add1);
      std::vector<int> add2(3);
      add2[0] = c;
      add2[1] = ncv[k];
      add2[2] = neighbor1;
      addtriangles.push_back(add2);
    }

    if (out)
    {
      std::cout << "   -> we have " << (int)addtriangles.size() << " potential NEW triangles"
                << std::endl;
      for (int k = 0; k < (int)addtriangles.size(); ++k)
        std::cout << "      Triangle: " << addtriangles[k][0] << " " << addtriangles[k][1] << " "
                  << addtriangles[k][2] << " " << std::endl;
    }

    // collapse addtriangles (remove double entries)
    int nadd = 0;
    for (int k = 0; k < (int)addtriangles.size(); ++k)
    {
      bool addbefore = false;
      int idx0 = addtriangles[k][0];
      int idx1 = addtriangles[k][1];
      int idx2 = addtriangles[k][2];

      // check against all other goodtriangles
      for (int l = 0; l < (int)goodtriangles.size(); ++l)
      {
        // do not check against itself
        int lidx0 = goodtriangles[l][0];
        int lidx1 = goodtriangles[l][1];
        int lidx2 = goodtriangles[l][2];

        if (idx0 == lidx0 && idx1 == lidx1 && idx2 == lidx2) addbefore = true;
        if (idx0 == lidx0 && idx1 == lidx2 && idx2 == lidx1) addbefore = true;
        if (idx0 == lidx1 && idx1 == lidx0 && idx2 == lidx2) addbefore = true;
        if (idx0 == lidx1 && idx1 == lidx2 && idx2 == lidx0) addbefore = true;
        if (idx0 == lidx2 && idx1 == lidx0 && idx2 == lidx1) addbefore = true;
        if (idx0 == lidx2 && idx1 == lidx1 && idx2 == lidx0) addbefore = true;
      }

      // add to good triangles
      if (!addbefore)
      {
        nadd++;
        goodtriangles.push_back(addtriangles[k]);
      }
    }

    // store final triangulation
    triangles.resize(0);
    for (int k = 0; k < (int)goodtriangles.size(); ++k) triangles.push_back(goodtriangles[k]);

    if (out)
    {
      std::cout << "   -> added " << nadd << " NEW triangles" << std::endl;
      std::cout << "   -> overall " << (int)triangles.size() << " triangles" << std::endl;
      for (int k = 0; k < (int)triangles.size(); ++k)
        std::cout << "      Triangle: " << triangles[k][0] << " " << triangles[k][1] << " "
                  << triangles[k][2] << " " << std::endl;
    }
  }

  // create intcells for all triangle
  int numt = (int)triangles.size();
  for (int t = 0; t < numt; ++t)
  {
    // indices of current triangle
    int idx0 = triangles[t][0];
    int idx1 = triangles[t][1];
    int idx2 = triangles[t][2];

    // coordinates of current triangle
    Core::LinAlg::Matrix<3, 3> coords;
    for (int k = 0; k < 3; ++k)
    {
      coords(k, 0) = clip()[idx0].coord()[k];
      coords(k, 1) = clip()[idx1].coord()[k];
      coords(k, 2) = clip()[idx2].coord()[k];
    }

    // create IntCell object and push back
    cells().push_back(std::make_shared<IntCell>(t, 3, coords, auxn(), Core::FE::CellType::tri3,
        linvertex[idx0], linvertex[idx1], linvertex[idx2], get_deriv_auxn()));
  }

  // double check number of triangles
  if (numt != clipsize - 2)
  {
    std::cout << "***WARNING*** Delaunay triangulation failed (" << clipsize << " vertices, "
              << numt << " triangles)"
              << " -> using backup" << std::endl;

    // if Delaunay triangulation failed, use old center-based
    // triangulation as backup (therefore return false)
    return false;
  }

  // triangulation successful
  return true;
}


/*----------------------------------------------------------------------*
 |  Triangulation of clip polygon (3D) - CENTER               popp 08/11|
 *----------------------------------------------------------------------*/
bool Mortar::Coupling3d::center_triangulation(
    std::vector<std::vector<Core::Gen::Pairedvector<int, double>>>& linvertex, double tol)
{
  // preparations
  cells().resize(0);
  int clipsize = (int)(clip().size());
  std::vector<Core::Gen::Pairedvector<int, double>> lincenter(
      3, (target_element().num_node() + source_element().num_node()) * 3);

  //**********************************************************************
  // (1) Trivial clipping polygon -> IntCells
  //**********************************************************************
  // clip polygon = triangle
  // no triangulation necessary -> 1 IntCell
  if (clipsize == 3)
  {
    // IntCell vertices = clip polygon vertices
    Core::LinAlg::Matrix<3, 3> coords;
    for (int i = 0; i < clipsize; ++i)
      for (int k = 0; k < 3; ++k) coords(k, i) = clip()[i].coord()[k];

    // create IntCell object and push back
    cells().push_back(std::make_shared<IntCell>(0, 3, coords, auxn(), Core::FE::CellType::tri3,
        linvertex[0], linvertex[1], linvertex[2], get_deriv_auxn()));

    // get out of here
    return true;
  }

  //**********************************************************************
  // (2) Find center of clipping polygon (centroid formula)
  //**********************************************************************
  std::vector<double> clipcenter(3);
  for (int k = 0; k < 3; ++k) clipcenter[k] = 0.0;
  double fac = 0.0;

  // first we need node averaged center
  std::array<double, 3> nac = {0.0, 0.0, 0.0};
  for (int i = 0; i < clipsize; ++i)
    for (int k = 0; k < 3; ++k) nac[k] += (clip()[i].coord()[k] / clipsize);

  // loop over all triangles of polygon
  for (int i = 0; i < clipsize; ++i)
  {
    std::array<double, 3> xi_i = {0.0, 0.0, 0.0};
    std::array<double, 3> xi_ip1 = {0.0, 0.0, 0.0};

    // standard case
    if (i < clipsize - 1)
    {
      for (int k = 0; k < 3; ++k) xi_i[k] = clip()[i].coord()[k];
      for (int k = 0; k < 3; ++k) xi_ip1[k] = clip()[i + 1].coord()[k];
    }
    // last vertex of clip polygon
    else
    {
      for (int k = 0; k < 3; ++k) xi_i[k] = clip()[clipsize - 1].coord()[k];
      for (int k = 0; k < 3; ++k) xi_ip1[k] = clip()[0].coord()[k];
    }

    // triangle area
    std::array<double, 3> diff1 = {0.0, 0.0, 0.0};
    std::array<double, 3> diff2 = {0.0, 0.0, 0.0};
    for (int k = 0; k < 3; ++k) diff1[k] = xi_ip1[k] - xi_i[k];
    for (int k = 0; k < 3; ++k) diff2[k] = xi_i[k] - nac[k];

    std::array<double, 3> cross = {0.0, 0.0, 0.0};
    cross[0] = diff1[1] * diff2[2] - diff1[2] * diff2[1];
    cross[1] = diff1[2] * diff2[0] - diff1[0] * diff2[2];
    cross[2] = diff1[0] * diff2[1] - diff1[1] * diff2[0];

    double Atri = 0.5 * sqrt(cross[0] * cross[0] + cross[1] * cross[1] + cross[2] * cross[2]);

    // add contributions to clipcenter and fac
    fac += Atri;
    for (int k = 0; k < 3; ++k) clipcenter[k] += 1.0 / 3.0 * (xi_i[k] + xi_ip1[k] + nac[k]) * Atri;
  }

  // final clipcenter
  for (int k = 0; k < 3; ++k) clipcenter[k] /= fac;

  //**********************************************************************
  // (3) Linearization of clip center coordinates (only for contact)
  //**********************************************************************
  center_linearization(linvertex, lincenter);

  //**********************************************************************
  // (4) General clipping polygon: Triangulation -> IntCells
  //**********************************************************************
  // center-based triangulation if clip polygon > quadrilateral
  // No. of IntCells is equal to no. of clip polygon vertices
  for (int num = 0; num < clipsize; ++num)
  {
    // the first vertex is always the clip center
    // the second vertex is always the current clip vertex
    Core::LinAlg::Matrix<3, 3> coords;
    for (int k = 0; k < 3; ++k)
    {
      coords(k, 0) = clipcenter[k];
      coords(k, 1) = clip()[num].coord()[k];
    }

    // the third vertex is the next vertex on clip polygon
    int numplus1 = num + 1;
    if (num == clipsize - 1)
    {
      for (int k = 0; k < 3; ++k) coords(k, 2) = clip()[0].coord()[k];
      numplus1 = 0;
    }
    else
      for (int k = 0; k < 3; ++k) coords(k, 2) = clip()[num + 1].coord()[k];

    // create IntCell object and push back
    cells().push_back(std::make_shared<IntCell>(num, 3, coords, auxn(), Core::FE::CellType::tri3,
        lincenter, linvertex[num], linvertex[numplus1], get_deriv_auxn()));
  }

  // triangulation successful
  return true;
}


/*----------------------------------------------------------------------*
 |  Integration of cells (3D)                                 popp 11/08|
 *----------------------------------------------------------------------*/
bool Mortar::Coupling3d::integrate_cells(
    const std::shared_ptr<Mortar::ParamsInterface>& mparams_ptr)
{
  /**********************************************************************/
  /* INTEGRATION                                                        */
  /* Integrate the Mortar matrices M (and possibly D) on the current    */
  /* integration cell of the slave / master (integration) element pair  */
  /**********************************************************************/

  // nothing to do here, if there are no cells
  if (cells().size() == 0) return false;

  // loop over all integration cells
  for (int i = 0; i < (int)(cells().size()); ++i)
  {
    // integrate cell only if it has a non-zero area
    if (cells()[i]->area() < MORTARINTLIM * source_element_area()) continue;

    // set segmentation status of all slave nodes
    // (hassegment_ of a slave node is true if ANY segment/cell
    // is integrated that contributes to this slave node)
    int nnodes = source_int_element().num_node();
    Core::Nodes::Node** mynodes = source_int_element().nodes();
    if (!mynodes) FOUR_C_THROW("Null pointer!");
    for (int k = 0; k < nnodes; ++k)
    {
      Mortar::Node* mycnode = dynamic_cast<Mortar::Node*>(mynodes[k]);
      if (!mycnode) FOUR_C_THROW("Null pointer!");
      mycnode->has_segment() = true;
    }

    // *******************************************************************
    // different options for mortar integration
    // *******************************************************************
    // (1) no quadratic element(s) involved -> linear LM interpolation
    // (2) quadratic element(s) involved -> quadratic LM interpolation
    // (3) quadratic element(s) involved -> linear LM interpolation
    // (4) quadratic element(s) involved -> piecew. linear LM interpolation
    // *******************************************************************
    Mortar::LagMultQuad lmtype = lag_mult_quad();

    // *******************************************************************
    // case (1)
    // *******************************************************************
    if (!quad())
    {
      // call integrator
      Mortar::Integrator::impl(source_element(), target_element(), interface_params())
          ->integrate_cell_3d_aux_plane(
              source_element(), target_element(), cells()[i], auxn(), get_comm());
    }

    // *******************************************************************
    // cases (2) and (3)
    // *******************************************************************
    else if (quad() && (lmtype == Mortar::lagmult_quad || lmtype == Mortar::lagmult_lin ||
                           lmtype == Mortar::lagmult_const))
    {
      // dynamic_cast to make sure to pass in IntElement&
      Mortar::IntElement& sintref = dynamic_cast<Mortar::IntElement&>(source_int_element());
      Mortar::IntElement& mintref = dynamic_cast<Mortar::IntElement&>(target_int_element());

      Mortar::Integrator::impl(source_element(), target_element(), interface_params())
          ->integrate_cell_3d_aux_plane_quad(
              source_element(), target_element(), sintref, mintref, cells()[i], auxn());
    }

    // *******************************************************************
    // case (4)
    // *******************************************************************
    else if (quad() && lmtype == Mortar::lagmult_pwlin)
    {
      // check for dual shape functions
      if (shape_fcn() == Mortar::shape_dual)
        FOUR_C_THROW(
            "ERROR: Piecewise linear LM interpolation not yet implemented for DUAL 3D quadratic "
            "mortar");

      // dynamic_cast to make sure to pass in IntElement&
      Mortar::IntElement& sintref = dynamic_cast<Mortar::IntElement&>(source_int_element());
      Mortar::IntElement& mintref = dynamic_cast<Mortar::IntElement&>(target_int_element());

      Mortar::Integrator::impl(source_element(), target_element(), interface_params())
          ->integrate_cell_3d_aux_plane_quad(
              source_element(), target_element(), sintref, mintref, cells()[i], auxn());
    }

    // *******************************************************************
    // undefined case
    // *******************************************************************
    else if (quad() && lmtype == Mortar::lagmult_undefined)
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
      FOUR_C_THROW("IntegrateCells: Invalid case for 3D mortar coupling LM interpolation");
    }
    // *******************************************************************
  }

  return true;
}


/*----------------------------------------------------------------------*
 |  Integration of cells (3D)                                 popp 07/11|
 *----------------------------------------------------------------------*/
void Mortar::Coupling3d::gmsh_output_cells(int lid) const
{
  // every processor writes its own cell file
  int proc = Core::Communication::my_mpi_rank(idiscret_.get_comm());
  int nproc = Core::Communication::num_mpi_ranks(idiscret_.get_comm());

  // write each integration cell only once
  // (no overlap, only owner of slave element writes output)
  if (proc != source_element().owner()) return;

  // construct unique filename for gmsh output
  // first index = time step index
  std::ostringstream filename;
  filename << "o/gmsh_output/cells_" << proc << ".pos";

  // do output to file in c-style
  FILE* fp = nullptr;

  // static variable
  static int count = 0;

  // open file
  if (count == 0)
    fp = fopen(filename.str().c_str(), "w");
  else
    fp = fopen(filename.str().c_str(), "a");

  // plot current integration cell
  const Core::LinAlg::Matrix<3, 3>& coord = cells()[lid]->coords();

  // write output to temporary std::stringstream
  std::stringstream gmshfilecontent;

  // header and dummy elements
  if (count == 0)
  {
    // header
    gmshfilecontent << "View \"Integration Cells Proc " << proc << "\" {" << std::endl;

    // dummy element 1
    gmshfilecontent << "ST(" << std::scientific << 0.0 << "," << 0.0 << "," << 0.0 << "," << 0.0
                    << "," << 0.0 << "," << 0.0 << "," << 0.0 << "," << 0.0 << "," << 0.0 << ")";
    gmshfilecontent << "{" << std::scientific << 0 << "," << 0 << "," << 0 << "};" << std::endl;

    // dummy element 1
    gmshfilecontent << "ST(" << std::scientific << 0.0 << "," << 0.0 << "," << 0.0 << "," << 0.0
                    << "," << 0.0 << "," << 0.0 << "," << 0.0 << "," << 0.0 << "," << 0.0 << ")";
    gmshfilecontent << "{" << std::scientific << nproc - 1 << "," << nproc - 1 << "," << nproc - 1
                    << "};" << std::endl;
  }

  // plot cell itself
  gmshfilecontent << "ST(" << std::scientific << coord(0, 0) << "," << coord(1, 0) << ","
                  << coord(2, 0) << "," << coord(0, 1) << "," << coord(1, 1) << "," << coord(2, 1)
                  << "," << coord(0, 2) << "," << coord(1, 2) << "," << coord(2, 2) << ")";
  gmshfilecontent << "{" << std::scientific << proc << "," << proc << "," << proc << "};"
                  << std::endl;

  // move everything to gmsh post-processing files and close them
  fprintf(fp, "%s", gmshfilecontent.str().c_str());
  fclose(fp);

  // increase static variable
  count += 1;

  return;
}


/*----------------------------------------------------------------------*
 | Split Mortar::Elements->IntElements for 3D quad. coupling    popp 03/09|
 *----------------------------------------------------------------------*/
bool Mortar::Coupling3dQuadManager::split_int_elements(
    Mortar::Element& ele, std::vector<std::shared_ptr<Mortar::IntElement>>& auxele)
{
  // *********************************************************************
  // do splitting for given element
  // *********************************************************** quad9 ***
  if (ele.shape() == Core::FE::CellType::quad9)
  {
    // split into for quad4 elements
    int numnode = 4;
    Core::FE::CellType dt = Core::FE::CellType::quad4;

    // first integration element
    // containing parent nodes 0,4,8,7
    int nodeids[4] = {0, 0, 0, 0};
    nodeids[0] = ele.node_ids()[0];
    nodeids[1] = ele.node_ids()[4];
    nodeids[2] = ele.node_ids()[8];
    nodeids[3] = ele.node_ids()[7];

    std::vector<Core::Nodes::Node*> nodes(4);
    nodes[0] = ele.nodes()[0];
    nodes[1] = ele.nodes()[4];
    nodes[2] = ele.nodes()[8];
    nodes[3] = ele.nodes()[7];

    auxele.push_back(std::make_shared<IntElement>(
        0, ele.id(), ele.owner(), &ele, dt, numnode, nodeids, nodes, ele.is_source(), false));

    // second integration element
    // containing parent nodes 4,1,5,8
    nodeids[0] = ele.node_ids()[4];
    nodeids[1] = ele.node_ids()[1];
    nodeids[2] = ele.node_ids()[5];
    nodeids[3] = ele.node_ids()[8];

    nodes[0] = ele.nodes()[4];
    nodes[1] = ele.nodes()[1];
    nodes[2] = ele.nodes()[5];
    nodes[3] = ele.nodes()[8];

    auxele.push_back(std::make_shared<IntElement>(
        1, ele.id(), ele.owner(), &ele, dt, numnode, nodeids, nodes, ele.is_source(), false));

    // third integration element
    // containing parent nodes 8,5,2,6
    nodeids[0] = ele.node_ids()[8];
    nodeids[1] = ele.node_ids()[5];
    nodeids[2] = ele.node_ids()[2];
    nodeids[3] = ele.node_ids()[6];

    nodes[0] = ele.nodes()[8];
    nodes[1] = ele.nodes()[5];
    nodes[2] = ele.nodes()[2];
    nodes[3] = ele.nodes()[6];

    auxele.push_back(std::make_shared<IntElement>(
        2, ele.id(), ele.owner(), &ele, dt, numnode, nodeids, nodes, ele.is_source(), false));

    // fourth integration element
    // containing parent nodes 7,8,6,3
    nodeids[0] = ele.node_ids()[7];
    nodeids[1] = ele.node_ids()[8];
    nodeids[2] = ele.node_ids()[6];
    nodeids[3] = ele.node_ids()[3];

    nodes[0] = ele.nodes()[7];
    nodes[1] = ele.nodes()[8];
    nodes[2] = ele.nodes()[6];
    nodes[3] = ele.nodes()[3];

    auxele.push_back(std::make_shared<IntElement>(
        3, ele.id(), ele.owner(), &ele, dt, numnode, nodeids, nodes, ele.is_source(), false));
  }

  // *********************************************************** quad8 ***
  else if (ele.shape() == Core::FE::CellType::quad8)
  {
    // split into four tri3 elements and one quad4 element
    int numnodetri = 3;
    int numnodequad = 4;
    Core::FE::CellType dttri = Core::FE::CellType::tri3;
    Core::FE::CellType dtquad = Core::FE::CellType::quad4;

    // first integration element
    // containing parent nodes 0,4,7
    int nodeids[3] = {0, 0, 0};
    nodeids[0] = ele.node_ids()[0];
    nodeids[1] = ele.node_ids()[4];
    nodeids[2] = ele.node_ids()[7];

    std::vector<Core::Nodes::Node*> nodes(3);
    nodes[0] = ele.nodes()[0];
    nodes[1] = ele.nodes()[4];
    nodes[2] = ele.nodes()[7];

    auxele.push_back(std::make_shared<IntElement>(
        0, ele.id(), ele.owner(), &ele, dttri, numnodetri, nodeids, nodes, ele.is_source(), false));

    // second integration element
    // containing parent nodes 1,5,4
    nodeids[0] = ele.node_ids()[1];
    nodeids[1] = ele.node_ids()[5];
    nodeids[2] = ele.node_ids()[4];

    nodes[0] = ele.nodes()[1];
    nodes[1] = ele.nodes()[5];
    nodes[2] = ele.nodes()[4];

    auxele.push_back(std::make_shared<IntElement>(
        1, ele.id(), ele.owner(), &ele, dttri, numnodetri, nodeids, nodes, ele.is_source(), false));

    // third integration element
    // containing parent nodes 2,6,5
    nodeids[0] = ele.node_ids()[2];
    nodeids[1] = ele.node_ids()[6];
    nodeids[2] = ele.node_ids()[5];

    nodes[0] = ele.nodes()[2];
    nodes[1] = ele.nodes()[6];
    nodes[2] = ele.nodes()[5];

    auxele.push_back(std::make_shared<IntElement>(
        2, ele.id(), ele.owner(), &ele, dttri, numnodetri, nodeids, nodes, ele.is_source(), false));

    // fourth integration element
    // containing parent nodes 3,7,6
    nodeids[0] = ele.node_ids()[3];
    nodeids[1] = ele.node_ids()[7];
    nodeids[2] = ele.node_ids()[6];

    nodes[0] = ele.nodes()[3];
    nodes[1] = ele.nodes()[7];
    nodes[2] = ele.nodes()[6];

    auxele.push_back(std::make_shared<IntElement>(
        3, ele.id(), ele.owner(), &ele, dttri, numnodetri, nodeids, nodes, ele.is_source(), false));

    // fifth integration element
    // containing parent nodes 4,5,6,7
    int nodeidsquad[4] = {0, 0, 0, 0};
    nodeidsquad[0] = ele.node_ids()[4];
    nodeidsquad[1] = ele.node_ids()[5];
    nodeidsquad[2] = ele.node_ids()[6];
    nodeidsquad[3] = ele.node_ids()[7];

    std::vector<Core::Nodes::Node*> nodesquad(4);
    nodesquad[0] = ele.nodes()[4];
    nodesquad[1] = ele.nodes()[5];
    nodesquad[2] = ele.nodes()[6];
    nodesquad[3] = ele.nodes()[7];

    auxele.push_back(std::make_shared<IntElement>(4, ele.id(), ele.owner(), &ele, dtquad,
        numnodequad, nodeidsquad, nodesquad, ele.is_source(), false));
  }

  // ************************************************************ tri6 ***
  else if (ele.shape() == Core::FE::CellType::tri6)
  {
    // split into four tri3 elements
    int numnode = 3;
    Core::FE::CellType dt = Core::FE::CellType::tri3;

    // first integration element
    // containing parent nodes 0,3,5
    int nodeids[3] = {0, 0, 0};
    nodeids[0] = ele.node_ids()[0];
    nodeids[1] = ele.node_ids()[3];
    nodeids[2] = ele.node_ids()[5];

    std::vector<Core::Nodes::Node*> nodes(3);
    nodes[0] = ele.nodes()[0];
    nodes[1] = ele.nodes()[3];
    nodes[2] = ele.nodes()[5];

    auxele.push_back(std::make_shared<IntElement>(
        0, ele.id(), ele.owner(), &ele, dt, numnode, nodeids, nodes, ele.is_source(), false));

    // second integration element
    // containing parent nodes 3,1,4
    nodeids[0] = ele.node_ids()[3];
    nodeids[1] = ele.node_ids()[1];
    nodeids[2] = ele.node_ids()[4];

    nodes[0] = ele.nodes()[3];
    nodes[1] = ele.nodes()[1];
    nodes[2] = ele.nodes()[4];

    auxele.push_back(std::make_shared<IntElement>(
        1, ele.id(), ele.owner(), &ele, dt, numnode, nodeids, nodes, ele.is_source(), false));

    // third integration element
    // containing parent nodes 5,4,2
    nodeids[0] = ele.node_ids()[5];
    nodeids[1] = ele.node_ids()[4];
    nodeids[2] = ele.node_ids()[2];

    nodes[0] = ele.nodes()[5];
    nodes[1] = ele.nodes()[4];
    nodes[2] = ele.nodes()[2];

    auxele.push_back(std::make_shared<IntElement>(
        2, ele.id(), ele.owner(), &ele, dt, numnode, nodeids, nodes, ele.is_source(), false));

    // fourth integration element
    // containing parent nodes 4,5,3
    nodeids[0] = ele.node_ids()[4];
    nodeids[1] = ele.node_ids()[5];
    nodeids[2] = ele.node_ids()[3];

    nodes[0] = ele.nodes()[4];
    nodes[1] = ele.nodes()[5];
    nodes[2] = ele.nodes()[3];

    auxele.push_back(std::make_shared<IntElement>(
        3, ele.id(), ele.owner(), &ele, dt, numnode, nodeids, nodes, ele.is_source(), false));
  }

  // *********************************************************** quad4 ***
  else if (ele.shape() == Core::FE::CellType::quad4)
  {
    // 1:1 conversion to IntElement
    std::vector<Core::Nodes::Node*> nodes(4);
    nodes[0] = ele.nodes()[0];
    nodes[1] = ele.nodes()[1];
    nodes[2] = ele.nodes()[2];
    nodes[3] = ele.nodes()[3];

    auxele.push_back(std::make_shared<IntElement>(0, ele.id(), ele.owner(), &ele, ele.shape(),
        ele.num_node(), ele.node_ids(), nodes, ele.is_source(), false));
  }

  // ************************************************************ tri3 ***
  else if (ele.shape() == Core::FE::CellType::tri3)
  {
    // 1:1 conversion to IntElement
    std::vector<Core::Nodes::Node*> nodes(3);
    nodes[0] = ele.nodes()[0];
    nodes[1] = ele.nodes()[1];
    nodes[2] = ele.nodes()[2];

    auxele.push_back(std::make_shared<IntElement>(0, ele.id(), ele.owner(), &ele, ele.shape(),
        ele.num_node(), ele.node_ids(), nodes, ele.is_source(), false));
  }

  // ************************************************************ nurbs9 ***
  else if (ele.shape() == Core::FE::CellType::nurbs9)
  {
    // create one IntElement from one nurbs9 element
    // new nodes are created as the images of the corners of the parameter
    // space on the actual surface.
    // nodes are ordered anti-clockwise using the info of the normal-fac
    // of the nurbs element

    // create pseudo-nodes
    std::vector<Core::Nodes::Node> pseudo_nodes;
    pseudo_nodes.clear();
    std::vector<Core::Nodes::Node*> pseudo_nodes_ptr;

    // parameter space coords of pseudo nodes
    double pseudo_nodes_param_coords[4][2];
    std::array<int, 4> id = {0, 0, 0, 0};
    bool rewind = false;

    if (ele.normal_fac() == 1.)
    {
      rewind = false;
      pseudo_nodes_param_coords[0][0] = -1.;
      pseudo_nodes_param_coords[0][1] = -1.;
      id[0] = ele.node_ids()[0];
      pseudo_nodes_param_coords[1][0] = +1.;
      pseudo_nodes_param_coords[1][1] = -1.;
      id[1] = ele.node_ids()[2];
      pseudo_nodes_param_coords[2][0] = +1.;
      pseudo_nodes_param_coords[2][1] = +1.;
      id[2] = ele.node_ids()[8];
      pseudo_nodes_param_coords[3][0] = -1.;
      pseudo_nodes_param_coords[3][1] = +1.;
      id[3] = ele.node_ids()[6];
    }
    else if (ele.normal_fac() == -1.)
    {
      rewind = true;
      pseudo_nodes_param_coords[0][0] = -1.;
      pseudo_nodes_param_coords[0][1] = -1.;
      id[0] = ele.node_ids()[0];
      pseudo_nodes_param_coords[1][0] = -1.;
      pseudo_nodes_param_coords[1][1] = +1.;
      id[1] = ele.node_ids()[6];
      pseudo_nodes_param_coords[2][0] = +1.;
      pseudo_nodes_param_coords[2][1] = +1.;
      id[2] = ele.node_ids()[8];
      pseudo_nodes_param_coords[3][0] = +1.;
      pseudo_nodes_param_coords[3][1] = -1.;
      id[3] = ele.node_ids()[2];
    }
    else
      FOUR_C_THROW("don't know what to do with this ele.NormalFac()");


    Core::LinAlg::SerialDenseVector sval(9);
    Core::LinAlg::SerialDenseMatrix sderiv(9, 2);
    std::vector<int> empty_dofs(3, -1);

    for (int i = 0; i < 4; ++i)
    {
      double xi[2] = {pseudo_nodes_param_coords[i][0], pseudo_nodes_param_coords[i][1]};
      std::vector<double> xspatial(dim_, 0.0);
      ele.evaluate_shape(xi, sval, sderiv, 9, true);
      for (int dim = 0; dim < dim_; ++dim)
        for (int n = 0; n < ele.num_node(); ++n)
          xspatial[dim] += sval(n) * dynamic_cast<Mortar::Node*>(ele.nodes()[n])->xspatial()[dim];
      pseudo_nodes.push_back(Mortar::Node(-1, xspatial, ele.owner(), empty_dofs, ele.is_source()));
    }

    for (int i = 0; i < 4; ++i) pseudo_nodes_ptr.push_back(&(pseudo_nodes[i]));

    auxele.push_back(std::make_shared<IntElement>(0, ele.id(), ele.owner(), &ele,
        Core::FE::CellType::quad4, 4, &(id[0]), pseudo_nodes_ptr, ele.is_source(), rewind));
  }

  // ********************************************************* invalid ***
  else
    FOUR_C_THROW("split_int_elements called for unknown element shape!");

  // *********************************************************************

  return true;
}


/*----------------------------------------------------------------------*
 |  ctor (public)                                             popp 06/09|
 *----------------------------------------------------------------------*/
Mortar::Coupling3dQuad::Coupling3dQuad(Core::FE::Discretization& idiscret, int dim, bool quad,
    Teuchos::ParameterList& params, Mortar::Element& sele, Mortar::Element& mele,
    Mortar::IntElement& sintele, Mortar::IntElement& mintele)
    : Mortar::Coupling3d(idiscret, dim, quad, params, sele, mele),
      sintele_(sintele),
      mintele_(mintele)
{
  //  3D quadratic coupling only for quadratic ansatz type
  if (!Coupling3dQuad::quad()) FOUR_C_THROW("Coupling3dQuad called for non-quadratic ansatz!");

  return;
}


/*----------------------------------------------------------------------*
 |  get communicator  (public)                               farah 01/13|
 *----------------------------------------------------------------------*/
MPI_Comm Mortar::Coupling3dManager::get_comm() const { return idiscret_.get_comm(); }


/*----------------------------------------------------------------------*
 |  ctor (public)                                             popp 06/09|
 *----------------------------------------------------------------------*/
Mortar::Coupling3dManager::Coupling3dManager(Core::FE::Discretization& idiscret, int dim, bool quad,
    Teuchos::ParameterList& params, Mortar::Element* sele, std::vector<Mortar::Element*> mele)
    : idiscret_(idiscret),
      dim_(dim),
      integrationtype_(Teuchos::getIntegralValue<Mortar::IntType>(params, "INTTYPE")),
      shapefcn_(Teuchos::getIntegralValue<Mortar::ShapeFcn>(params, "LM_SHAPEFCN")),
      lmdualconsistent_(
          Teuchos::getIntegralValue<Mortar::ConsistentDualType>(params, "LM_DUAL_CONSISTENT")),
      quad_(quad),
      imortar_(params),
      sele_(sele),
      mele_(mele)
{
  return;
}


/*----------------------------------------------------------------------*
 |  ctor (public) -- empty                                   farah 01/13|
 *----------------------------------------------------------------------*/
Mortar::Coupling3dQuadManager::Coupling3dQuadManager(Core::FE::Discretization& idiscret, int dim,
    bool quad, Teuchos::ParameterList& params, Mortar::Element* sele,
    std::vector<Mortar::Element*> mele)
    : Coupling3dManager(idiscret, dim, quad, params, sele, mele)
{
  return;
}


/*----------------------------------------------------------------------*
 |  Evaluate coupling pairs                                  farah 10/14|
 *----------------------------------------------------------------------*/
bool Mortar::Coupling3dManager::evaluate_coupling(
    std::shared_ptr<Mortar::ParamsInterface> mparams_ptr)
{
  // check of we need to start the real coupling
  if (target_elements().size() == 0) return false;

  // decide which type of coupling should be evaluated
  auto algo = Teuchos::getIntegralValue<Mortar::AlgorithmType>(imortar_, "ALGORITHM");

  //*********************************
  // Mortar Contact
  //*********************************
  if (algo == Mortar::algorithm_mortar or algo == Mortar::algorithm_gpts)
    integrate_coupling(mparams_ptr);

  //*********************************
  // Error
  //*********************************
  else
    FOUR_C_THROW("chosen contact algorithm not supported!");

  return true;
}


/*----------------------------------------------------------------------*
 |  Evaluate mortar-coupling pairs                            popp 03/09|
 *----------------------------------------------------------------------*/
void Mortar::Coupling3dManager::integrate_coupling(
    const std::shared_ptr<Mortar::ParamsInterface>& mparams_ptr)
{
  // decide which type of numerical integration scheme

  //**********************************************************************
  // STANDARD INTEGRATION (SEGMENTS)
  //**********************************************************************
  if (int_type() == Mortar::inttype_segments)
  {
    // loop over all master elements associated with this slave element
    for (int m = 0; m < (int)target_elements().size(); ++m)
    {
      // create Coupling3d object and push back
      coupling().push_back(std::make_shared<Coupling3d>(
          idiscret_, dim_, false, imortar_, source_element(), target_element(m)));
      // do coupling
      coupling()[m]->evaluate_coupling();
    }

    // special treatment of boundary elements
    // calculate consistent dual shape functions for this element
    consist_dual_shape();

    // integrate cells
    for (int m = 0; m < (int)target_elements().size(); ++m)
      coupling()[m]->integrate_cells(mparams_ptr);
  }
  //**********************************************************************
  // FAST INTEGRATION (ELEMENTS)
  //**********************************************************************
  else if (int_type() == Mortar::inttype_elements || int_type() == Mortar::inttype_elements_BS)
  {
    if ((int)target_elements().size() == 0) return;

    if (!quad())
    {
      bool boundary_ele = false;

      // integrate D and M -- 2 Cells
      Mortar::Integrator::impl(source_element(), target_element(0), imortar_)
          ->integrate_ele_based_3d(
              source_element(), target_elements(), &boundary_ele, idiscret_.get_comm());

      if (int_type() == Mortar::inttype_elements_BS)
      {
        if (boundary_ele == true)
        {
          if (lmdualconsistent_ != Mortar::consistent_none)
          {
            // loop over all master elements associated with this slave element
            for (int m = 0; m < (int)target_elements().size(); ++m)
            {
              // create Coupling3d object and push back
              coupling().push_back(std::make_shared<Coupling3d>(
                  idiscret_, dim_, false, imortar_, source_element(), target_element(m)));
              // do coupling
              coupling()[m]->evaluate_coupling();

              // integrate cells
              coupling()[m]->integrate_cells(mparams_ptr);
            }
          }

          // consistent boundary treatment
          else
          {
            // loop over all master elements associated with this slave element
            for (int m = 0; m < (int)target_elements().size(); ++m)
            {
              // create Coupling3d object and push back
              coupling().push_back(std::make_shared<Coupling3d>(
                  idiscret_, dim_, false, imortar_, source_element(), target_element(m)));
              // do coupling
              coupling()[m]->evaluate_coupling();
            }

            // consistent dual shape functions
            consist_dual_shape();

            // integrate cells
            for (int m = 0; m < (int)target_elements().size(); ++m)
              coupling()[m]->integrate_cells(mparams_ptr);
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
    else
    {
      FOUR_C_THROW(
          "ERROR: You should not be here! This coupling manager is not able to perform mortar "
          "coupling for high-order elements.");
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
  source_element().mo_data().reset_dual_shape();
  source_element().mo_data().reset_deriv_dual_shape();

  return;
}


/*----------------------------------------------------------------------*
 |  Evaluate coupling pairs for Quad-coupling                farah 01/13|
 *----------------------------------------------------------------------*/
void Mortar::Coupling3dQuadManager::integrate_coupling(
    const std::shared_ptr<Mortar::ParamsInterface>& mparams_ptr)
{
  // decide which type of numerical integration scheme

  //**********************************************************************
  // STANDARD INTEGRATION (SEGMENTS)
  //**********************************************************************
  if (int_type() == Mortar::inttype_segments)
  {
    // build linear integration elements from quadratic Mortar::Elements
    std::vector<std::shared_ptr<Mortar::IntElement>> sauxelements;
    std::vector<std::vector<std::shared_ptr<Mortar::IntElement>>> mauxelements(
        target_elements().size());
    split_int_elements(source_element(), sauxelements);

    // loop over all master elements associated with this slave element
    for (int m = 0; m < (int)target_elements().size(); ++m)
    {
      // build linear integration elements from quadratic Mortar::Elements
      mauxelements[m].resize(0);
      split_int_elements(*target_elements()[m], mauxelements[m]);

      // loop over all IntElement pairs for coupling
      for (int i = 0; i < (int)sauxelements.size(); ++i)
      {
        for (int j = 0; j < (int)mauxelements[m].size(); ++j)
        {
          // create instance of coupling class
          coupling().push_back(std::make_shared<Coupling3dQuad>(idiscret_, dim_, true, imortar_,
              source_element(), *target_elements()[m], *sauxelements[i], *mauxelements[m][j]));

          // do coupling
          coupling()[coupling().size() - 1]->evaluate_coupling();
        }  // for maux
      }  // for saux
    }  // for m

    consist_dual_shape();

    // do integration
    for (int i = 0; i < (int)coupling().size(); ++i) coupling()[i]->integrate_cells(mparams_ptr);
  }
  //**********************************************************************
  // FAST INTEGRATION (ELEMENTS)
  //**********************************************************************
  else if (int_type() == Mortar::inttype_elements || int_type() == Mortar::inttype_elements_BS)
  {
    if ((int)target_elements().size() == 0) return;

    bool boundary_ele = false;

    // integrate D and M -- 2 Cells
    Mortar::Integrator::impl(source_element(), target_element(0), imortar_)
        ->integrate_ele_based_3d(
            source_element(), target_elements(), &boundary_ele, idiscret_.get_comm());

    if (int_type() == Mortar::inttype_elements_BS)
    {
      if (boundary_ele == true)
      {
        // loop over all master elements associated with this slave element
        for (int m = 0; m < (int)target_elements().size(); ++m)
        {
          // build linear integration elements from quadratic Mortar::Elements
          std::vector<std::shared_ptr<Mortar::IntElement>> sauxelements;
          std::vector<std::shared_ptr<Mortar::IntElement>> mauxelements;
          split_int_elements(source_element(), sauxelements);
          split_int_elements(*target_elements()[m], mauxelements);

          // loop over all IntElement pairs for coupling
          for (int i = 0; i < (int)sauxelements.size(); ++i)
          {
            for (int j = 0; j < (int)mauxelements.size(); ++j)
            {
              // create instance of coupling class
              Mortar::Coupling3dQuad coup(idiscret_, dim_, true, imortar_, source_element(),
                  *target_elements()[m], *sauxelements[i], *mauxelements[j]);
              // do coupling
              coup.evaluate_coupling();

              // Integrate Cells
              coup.integrate_cells(mparams_ptr);
            }  // for maux
          }  // for saux
        }  // for m
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
  //**********************************************************************
  // INVALID
  //**********************************************************************
  else
  {
    FOUR_C_THROW("Invalid type of numerical integration");
  }

  // free memory of consistent dual shape function coefficient matrix
  source_element().mo_data().reset_dual_shape();
  source_element().mo_data().reset_deriv_dual_shape();

  return;
}

/*----------------------------------------------------------------------*
 |  Calculate dual shape functions                           seitz 07/13|
 *----------------------------------------------------------------------*/
void Mortar::Coupling3dManager::consist_dual_shape()
{
  // For standard shape functions no modification is necessary
  // A switch earlier in the process improves computational efficiency
  if (shape_fcn() == Mortar::shape_standard || lmdualconsistent_ == Mortar::consistent_none) return;

  // Consistent modification not yet checked for constant LM interpolation
  if (quad() == true && lag_mult_quad() == Mortar::lagmult_const &&
      lmdualconsistent_ != Mortar::consistent_none)
    FOUR_C_THROW(
        "ERROR: Consistent dual shape functions not yet checked for constant LM interpolation!");

  if (coupling().size() == 0) return;

  if (int_type() == Mortar::inttype_segments && lmdualconsistent_ == Mortar::consistent_boundary)
  {
    // check if fully projecting
    bool boundary_ele = false;

    Mortar::ElementIntegrator integrator(source_element().shape());
    for (int gp = 0; gp < integrator.n_gp(); ++gp)
    {
      // coordinates and weight
      std::array<double, 2> eta = {integrator.coordinate(gp, 0), integrator.coordinate(gp, 1)};

      // note that the third component of sxi is necessary!
      // (although it will always be 0.0 of course)
      // double tempsxi[3] = {0.0, 0.0, 0.0};
      double sxi[2] = {0.0, 0.0};
      double mxi[2] = {0.0, 0.0};
      double projalpha = 0.0;

      sxi[0] = eta[0];
      sxi[1] = eta[1];

      // check for Boundary Segmentation
      bool projectable_gp = false;

      // discretization type of master element
      Core::FE::CellType dt = coupling()[0]->target_element().shape();

      //*******************************************************************
      // loop over meles
      //*******************************************************************
      for (int nummaster = 0; nummaster < (int)coupling().size(); ++nummaster)
      {
        // project Gauss point onto master element
        Mortar::Projector::impl(source_element(), coupling()[nummaster]->target_element())
            ->project_gauss_point_3d(
                source_element(), sxi, coupling()[nummaster]->target_element(), mxi, projalpha);

        bool is_on_mele = true;

        // check GP projection
        double tol = 0.00;
        if (dt == Core::FE::CellType::quad4 || dt == Core::FE::CellType::quad8 ||
            dt == Core::FE::CellType::quad9)
        {
          if (mxi[0] < -1.0 - tol || mxi[1] < -1.0 - tol || mxi[0] > 1.0 + tol ||
              mxi[1] > 1.0 + tol)
            is_on_mele = false;
        }
        else
        {
          if (mxi[0] < -tol || mxi[1] < -tol || mxi[0] > 1.0 + tol || mxi[1] > 1.0 + tol ||
              mxi[0] + mxi[1] > 1.0 + 2 * tol)
            is_on_mele = false;
        }
        if (is_on_mele == true)
        {
          projectable_gp = true;
          break;
        }
      }  // loop over meles

      if (projectable_gp == false)
      {
        boundary_ele = true;
        break;
      }
    }  // gp-loop
    if (boundary_ele == false) return;
  }

  // get number of nodes of present slave element
  int nnodes = source_element().num_node();

  // initialize Jacobian determinant
  double detg = 0.0;

  // initialize matrices de and me
  Core::LinAlg::SerialDenseMatrix me(nnodes, nnodes, true);
  Core::LinAlg::SerialDenseMatrix de(nnodes, nnodes, true);

  // loop over all master elements associated with this slave element
  for (int m = 0; m < (int)target_elements().size(); ++m)
  {
    // loop over all integration cells
    for (int c = 0; c < (int)coupling()[m]->cells().size(); ++c)
    {
      std::shared_ptr<Mortar::IntCell> currcell = coupling()[m]->cells()[c];

      // create an integrator for this cell
      for (int gp = 0;
          gp < Mortar::Integrator::impl(source_element(), target_element(m), imortar_)->n_gp();
          ++gp)
      {
        // coordinates and weight
        double eta[2] = {Mortar::Integrator::impl(source_element(), target_element(m), imortar_)
                             ->coordinate(gp, 0),
            Mortar::Integrator::impl(source_element(), target_element(m), imortar_)
                ->coordinate(gp, 1)};
        double wgt =
            Mortar::Integrator::impl(source_element(), target_element(m), imortar_)->weight(gp);

        // get global Gauss point coordinates
        double globgp[3] = {0.0, 0.0, 0.0};
        currcell->local_to_global(eta, globgp, 0);

        // project Gauss point onto slave integration element
        double sxi[2] = {0.0, 0.0};
        double sprojalpha = 0.0;

        // TODO random?
        Mortar::Projector::impl(source_element())
            ->project_gauss_point_auxn_3d(
                globgp, coupling()[m]->auxn(), source_element(), sxi, sprojalpha);

        // create vector for shape function evaluation
        Core::LinAlg::SerialDenseVector sval(nnodes);
        Core::LinAlg::SerialDenseMatrix sderiv(nnodes, 2, true);

        // evaluate trace space shape functions at Gauss point
        if (lag_mult_quad() == Mortar::lagmult_lin)
          source_element().evaluate_shape_lag_mult_lin(
              Mortar::shape_standard, sxi, sval, sderiv, nnodes);
        else
          source_element().evaluate_shape(sxi, sval, sderiv, nnodes);

        detg = currcell->jacobian();

        // computing de, me
        for (int j = 0; j < nnodes; ++j)
        {
          double fac;
          fac = sval[j] * wgt;

          // computing de
          de(j, j) += fac * detg;

          for (int k = 0; k < nnodes; ++k)
          {
            // computing me
            fac = wgt * sval[j] * sval[k];
            me(j, k) += fac * detg;
          }
        }
      }
    }  // cells
  }  // mele

  // in case of no overlap just return, as there is no integration area
  // and therefore the consistent dual shape functions are not defined.
  // This doesn't matter, as there is no associated integration domain anyway
  if (Core::LinAlg::det_long(me) == 0) return;

  // declare dual shape functions coefficient matrix
  Core::LinAlg::SerialDenseMatrix ae(nnodes, nnodes, true);

  // compute matrix A_e for linear interpolation of quadratic element
  if (lag_mult_quad() == Mortar::lagmult_lin)
  {
    // declare and initialize to zero inverse of Matrix M_e
    Core::LinAlg::SerialDenseMatrix meinv(nnodes, nnodes, true);

    if (source_element().shape() == Core::FE::CellType::tri6)
    {
      // reduce me to non-zero nodes before inverting
      Core::LinAlg::Matrix<3, 3> melin;
      for (int j = 0; j < 3; ++j)
        for (int k = 0; k < 3; ++k) melin(j, k) = me(j, k);

      // invert bi-ortho matrix melin
      Core::LinAlg::inverse(melin);

      // re-inflate inverse of melin to full size
      for (int j = 0; j < 3; ++j)
        for (int k = 0; k < 3; ++k) meinv(j, k) = melin(j, k);
    }
    else if (source_element().shape() == Core::FE::CellType::quad8 ||
             source_element().shape() == Core::FE::CellType::quad9)
    {
      // reduce me to non-zero nodes before inverting
      Core::LinAlg::Matrix<4, 4> melin;
      for (int j = 0; j < 4; ++j)
        for (int k = 0; k < 4; ++k) melin(j, k) = me(j, k);

      // invert bi-ortho matrix melin
      Core::LinAlg::inverse(melin);

      // re-inflate inverse of melin to full size
      for (int j = 0; j < 4; ++j)
        for (int k = 0; k < 4; ++k) meinv(j, k) = melin(j, k);
    }
    else
      FOUR_C_THROW("incorrect element shape for linear interpolation of quadratic element!");

    // get solution matrix with dual parameters
    Core::LinAlg::multiply(ae, de, meinv);
  }
  // compute matrix A_e for all other cases
  else
    Core::LinAlg::invert_and_multiply_by_cholesky(me, de, ae);

  // store ae matrix in slave element data container
  source_element().mo_data().dual_shape() = std::make_shared<Core::LinAlg::SerialDenseMatrix>(ae);

  return;
}

FOUR_C_NAMESPACE_CLOSE
