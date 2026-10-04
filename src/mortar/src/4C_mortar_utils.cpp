// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_mortar_utils.hpp"

#include "4C_comm_exporter.hpp"
#include "4C_fem_nurbs_discretization.hpp"
#include "4C_fem_nurbs_discretization_control_point.hpp"
#include "4C_fem_nurbs_discretization_knotvector.hpp"
#include "4C_linalg_sparsematrix.hpp"
#include "4C_linalg_utils_densematrix_communication.hpp"
#include "4C_linalg_utils_sparse_algebra_create.hpp"
#include "4C_linalg_utils_sparse_algebra_manipulation.hpp"
#include "4C_linalg_utils_sparse_algebra_math.hpp"
#include "4C_structure_new_timint_base.hpp"
#include "4C_utils_exceptions.hpp"

#include <algorithm>
#include <numeric>

FOUR_C_NAMESPACE_OPEN

int Mortar::sort_convex_hull_points(bool out, Core::LinAlg::SerialDenseMatrix& transformed,
    std::vector<Vertex>& collconvexhull, std::vector<Vertex>& respoly, double& tol)
{
  // The tolerance is not used while building the hull, since points could then be removed based
  // on neighbors that are not part of the final polygon. Removing the straightest point first
  // keeps real corners of tiny polygons, in which several points are below the tolerance.
  const int np = static_cast<int>(collconvexhull.size());

  // sort points lexicographically w.r.t. their x- and y-values
  auto is_smaller = [&](int a, int b)
  {
    if (transformed(0, a) != transformed(0, b)) return transformed(0, a) < transformed(0, b);
    return transformed(1, a) < transformed(1, b);
  };
  std::vector<int> sorted(np);
  std::iota(sorted.begin(), sorted.end(), 0);
  std::sort(sorted.begin(), sorted.end(), is_smaller);

  // a-b-c rotates clockwise
  auto is_clockwise = [&](int a, int b, int c)
  {
    const double abx = transformed(0, b) - transformed(0, a);
    const double aby = transformed(1, b) - transformed(1, a);
    const double acx = transformed(0, c) - transformed(0, a);
    const double acy = transformed(1, c) - transformed(1, a);
    return abx * acy - aby * acx < 0.0;
  };

  // (1) exact convex hull: upper chain from left to right, then lower chain from right to left
  std::vector<int> hull;
  for (int i = 0; i < np; ++i)
  {
    while (hull.size() >= 2 and not is_clockwise(hull[hull.size() - 2], hull.back(), sorted[i]))
      hull.pop_back();
    hull.push_back(sorted[i]);
  }
  const std::size_t upper_size = hull.size();
  for (int i = np - 2; i >= 0; --i)
  {
    while (hull.size() > upper_size and
           not is_clockwise(hull[hull.size() - 2], hull.back(), sorted[i]))
      hull.pop_back();
    hull.push_back(sorted[i]);
  }
  // the lower chain ends at the starting point again
  if (hull.size() > 1) hull.pop_back();

  // (2) remove points that form an almost straight line with their neighbors, starting with the
  // straightest one (if less than three points remain, there is no clip polygon)
  // (the cross product of the two adjacent edges is negative for clockwise corners)
  // Note: the cross product is an area, but tol is a length (scaled with the element size), as in
  // the inside/outside checks of the polygon clipping. Thus, on small polygons, real corners can
  // be removed. Dividing the cross product by the distance between the two neighbors would turn it
  // into the distance of the point from the edge between them, consistent with the clipping. This
  // is not done yet, since it changes the results of many tests slightly.
  auto clockwise_rotation = [&](std::size_t k)
  {
    const int a = hull[(k + hull.size() - 1) % hull.size()];
    const int b = hull[k];
    const int c = hull[(k + 1) % hull.size()];
    const double abx = transformed(0, b) - transformed(0, a);
    const double aby = transformed(1, b) - transformed(1, a);
    const double bcx = transformed(0, c) - transformed(0, b);
    const double bcy = transformed(1, c) - transformed(1, b);
    return -(abx * bcy - aby * bcx);
  };
  while (hull.size() >= 3)
  {
    std::size_t straightest = 0;
    for (std::size_t k = 1; k < hull.size(); ++k)
      if (clockwise_rotation(k) < clockwise_rotation(straightest)) straightest = k;
    if (clockwise_rotation(straightest) >= tol) break;
    hull.erase(hull.begin() + straightest);
  }

  // start again at the point with the smallest x-value, if the starting point was removed
  std::rotate(hull.begin(), std::min_element(hull.begin(), hull.end(), is_smaller), hull.end());

  for (int i : hull)
  {
    const Vertex& current = collconvexhull[i];
    respoly.push_back(Vertex(current.coord(), current.v_type(), current.nodeids(), nullptr, nullptr,
        false, false, nullptr, -1.0));

    if (out)
      std::cout << "Clip polygon point " << i << "\t" << transformed(0, i) << "\t"
                << transformed(1, i) << std::endl;
  }

  // number of points removed from convex hull
  return np - static_cast<int>(respoly.size());
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Mortar::Utils::create_volume_ghosting(const Core::FE::Discretization& dis_src,
    const std::vector<std::shared_ptr<Core::FE::Discretization>>& voldis,
    std::vector<std::pair<int, int>> material_links, bool check_on_in, bool check_on_exit)
{
  if (voldis.size() == 0) return;

  if (check_on_in)
    for (int c = 1; c < (int)voldis.size(); ++c)
      if (voldis.at(c)->element_row_map()->same_as(*voldis.at(0)->element_row_map()) == false)
        FOUR_C_THROW("row maps on input do not coincide");

  const Core::LinAlg::Map* ielecolmap = dis_src.element_col_map();

  // 1 Ghost all Volume Element + Nodes,for all col elements in dis_src
  for (unsigned disidx = 0; disidx < voldis.size(); ++disidx)
  {
    std::vector<int> rdata;

    // Fill rdata with existing colmap
    const Core::LinAlg::Map* elecolmap = voldis[disidx]->element_col_map();
    const std::shared_ptr<Core::LinAlg::Map> allredelecolmap =
        Core::LinAlg::allreduce_e_map(*voldis[disidx]->element_row_map());

    for (int i = 0; i < elecolmap->num_my_elements(); ++i)
    {
      int gid = elecolmap->gid(i);
      rdata.push_back(gid);
    }

    // Find elements, which are ghosted on the interface but not in the volume discretization
    for (int i = 0; i < ielecolmap->num_my_elements(); ++i)
    {
      int gid = ielecolmap->gid(i);

      Core::Elements::Element* ele = dis_src.g_element(gid);
      if (!ele) FOUR_C_THROW("Cannot find element with gid %", gid);
      Core::Elements::FaceElement* faceele = dynamic_cast<Core::Elements::FaceElement*>(ele);
      if (!faceele) FOUR_C_THROW("source element is not a face element");
      int volgid = faceele->parent_element_id();
      // Ghost the parent element additionally
      if (elecolmap->lid(volgid) == -1 &&
          allredelecolmap->lid(volgid) !=
              -1)  // Volume discretization has not Element on this proc but on another
        rdata.push_back(volgid);
    }

    // re-build element column map
    Core::LinAlg::Map newelecolmap(
        -1, (int)rdata.size(), rdata.data(), 0, voldis[disidx]->get_comm());
    rdata.clear();

    // redistribute the volume discretization according to the
    // new (=old) element column layout & and ghost also nodes!
    voldis[disidx]->extended_ghosting(newelecolmap, true, true, true, false);  // no check!!!
  }

  // 2 Reconnect Face Element -- Parent Element Pointers to first dis in dis_tar
  {
    const Core::LinAlg::Map* elecolmap = voldis[0]->element_col_map();

    for (int i = 0; i < ielecolmap->num_my_elements(); ++i)
    {
      int gid = ielecolmap->gid(i);

      Core::Elements::Element* ele = dis_src.g_element(gid);
      if (!ele) FOUR_C_THROW("Cannot find element with gid %", gid);
      Core::Elements::FaceElement* faceele = dynamic_cast<Core::Elements::FaceElement*>(ele);
      if (!faceele) FOUR_C_THROW("source element is not a face element");
      int volgid = faceele->parent_element_id();

      if (elecolmap->lid(volgid) == -1)  // Volume discretization has not Element
        FOUR_C_THROW("create_volume_ghosting: Element {} does not exist on this Proc!", volgid);

      Core::Elements::Element* vele = voldis[0]->g_element(volgid);
      if (!vele) FOUR_C_THROW("Cannot find element with gid %", volgid);

      faceele->set_parent_target_element(vele, faceele->face_parent_number());

      if (voldis.size() == 2)
      {
        const auto* elecolmap2 = voldis[1]->element_col_map();
        if (elecolmap2->lid(volgid) == -1)
          faceele->set_parent_source_element(nullptr, -1);
        else
        {
          auto* volele = voldis[1]->g_element(volgid);
          if (volele == nullptr) FOUR_C_THROW("Cannot find element with gid %", volgid);
          faceele->set_parent_source_element(volele, faceele->face_parent_number());
        }
      }
    }
  }

  if (check_on_exit)
    for (int c = 1; c < (int)voldis.size(); ++c)
    {
      if (voldis.at(c)->element_row_map()->same_as(*voldis.at(0)->element_row_map()) == false)
        FOUR_C_THROW("row maps on exit do not coincide");
      if (voldis.at(c)->element_col_map()->same_as(*voldis.at(0)->element_col_map()) == false)
        FOUR_C_THROW("col maps on exit do not coincide");
    }

  // 3 setup material pointers between newly ghosted elements
  for (std::vector<std::pair<int, int>>::const_iterator m = material_links.begin();
      m != material_links.end(); ++m)
  {
    std::shared_ptr<Core::FE::Discretization> dis_src_mat = voldis.at(m->first);
    std::shared_ptr<Core::FE::Discretization> dis_tar_mat = voldis.at(m->second);

    for (int i = 0; i < dis_tar_mat->num_my_col_elements(); ++i)
    {
      Core::Elements::Element* targetele = dis_tar_mat->l_col_element(i);
      const int gid = targetele->id();

      Core::Elements::Element* sourceele = dis_src_mat->g_element(gid);

      targetele->add_material(sourceele->material());
    }
  }
}



/*----------------------------------------------------------------------*
 |  Prepare mortar element for nurbs-case                    farah 11/14|
 *----------------------------------------------------------------------*/
void Mortar::Utils::prepare_nurbs_element(Core::FE::Discretization& discret,
    std::shared_ptr<Core::Elements::Element> ele, Mortar::Element& cele, int dim)
{
  Core::FE::Nurbs::NurbsDiscretization* nurbsdis =
      dynamic_cast<Core::FE::Nurbs::NurbsDiscretization*>(&(discret));

  std::shared_ptr<Core::FE::Nurbs::Knotvector> knots = (*nurbsdis).get_knot_vector();
  std::vector<Core::LinAlg::SerialDenseVector> parentknots(dim);
  std::vector<Core::LinAlg::SerialDenseVector> mortarknots(dim - 1);

  double normalfac = 0.0;
  std::shared_ptr<Core::Elements::FaceElement> faceele =
      std::dynamic_pointer_cast<Core::Elements::FaceElement>(ele);
  bool zero_size = knots->get_boundary_ele_and_parent_knots(parentknots, mortarknots, normalfac,
      faceele->parent_target_element()->id(), faceele->face_target_number());

  // store nurbs specific data to node
  cele.zero_sized() = zero_size;
  cele.knots() = mortarknots;
  cele.normal_fac() = normalfac;

  return;
}


/*----------------------------------------------------------------------*
 |  Prepare mortar node for nurbs-case                       farah 11/14|
 *----------------------------------------------------------------------*/
void Mortar::Utils::prepare_nurbs_node(Core::Nodes::Node* node, Mortar::Node& mnode)
{
  Core::FE::Nurbs::ControlPoint* cp = dynamic_cast<Core::FE::Nurbs::ControlPoint*>(node);

  mnode.nurbs_w() = cp->w();

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Mortar::Utils::mortar_matrix_condensation(std::shared_ptr<Core::LinAlg::SparseMatrix>& k,
    const std::shared_ptr<const Core::LinAlg::SparseMatrix>& p_row,
    const std::shared_ptr<const Core::LinAlg::SparseMatrix>& p_col)
{
  // prepare maps by making a deep copy of the map and wrap it in a shared_ptr
  auto gsrow = std::make_shared<Core::LinAlg::Map>(p_row->range_map());
  auto gmrow = std::make_shared<Core::LinAlg::Map>(p_row->domain_map());
  auto gscol = std::make_shared<Core::LinAlg::Map>(p_col->range_map());
  auto gmcol = std::make_shared<Core::LinAlg::Map>(p_col->domain_map());

  std::shared_ptr<Core::LinAlg::Map> gsmrow = Core::LinAlg::merge_map(gsrow, gmrow, false);
  std::shared_ptr<Core::LinAlg::Map> gnrow = Core::LinAlg::split_map(k->range_map(), *gsmrow);

  std::shared_ptr<Core::LinAlg::Map> gsmcol = Core::LinAlg::merge_map(gscol, gmcol, false);
  std::shared_ptr<Core::LinAlg::Map> gncol = Core::LinAlg::split_map(k->domain_map(), *gsmcol);

  /*--------------------------------------------------------------------*/
  /* Split kteff into 3x3 block matrix                                  */
  /*--------------------------------------------------------------------*/
  // we want to split k into 3 groups s,m,n = 9 blocks
  std::shared_ptr<Core::LinAlg::SparseMatrix> kss = nullptr;
  std::shared_ptr<Core::LinAlg::SparseMatrix> ksm = nullptr;
  std::shared_ptr<Core::LinAlg::SparseMatrix> ksn = nullptr;
  std::shared_ptr<Core::LinAlg::SparseMatrix> kms = nullptr;
  std::shared_ptr<Core::LinAlg::SparseMatrix> kmm = nullptr;
  std::shared_ptr<Core::LinAlg::SparseMatrix> kmn = nullptr;
  std::shared_ptr<Core::LinAlg::SparseMatrix> kns = nullptr;
  std::shared_ptr<Core::LinAlg::SparseMatrix> knm = nullptr;
  std::shared_ptr<Core::LinAlg::SparseMatrix> knn = nullptr;

  // temporarily we need the blocks ksmsm, ksmn, knsm
  // (FIXME: because a direct SplitMatrix3x3 is still missing!)
  std::shared_ptr<Core::LinAlg::SparseMatrix> ksmsm = nullptr;
  std::shared_ptr<Core::LinAlg::SparseMatrix> ksmn = nullptr;
  std::shared_ptr<Core::LinAlg::SparseMatrix> knsm = nullptr;

  // some temporary std::shared_ptrs
  std::shared_ptr<Core::LinAlg::Map> tempmap;
  std::shared_ptr<Core::LinAlg::SparseMatrix> tempmtx1 = nullptr;
  std::shared_ptr<Core::LinAlg::SparseMatrix> tempmtx2 = nullptr;

  // split
  Core::LinAlg::split_matrix2x2(k, gsmrow, gnrow, gsmcol, gncol, ksmsm, ksmn, knsm, knn);
  Core::LinAlg::split_matrix2x2(ksmsm, gsrow, gmrow, gscol, gmcol, kss, ksm, kms, kmm);
  Core::LinAlg::split_matrix2x2(ksmn, gsrow, gmrow, gncol, tempmap, ksn, tempmtx1, kmn, tempmtx2);
  Core::LinAlg::split_matrix2x2(knsm, gnrow, tempmap, gscol, gmcol, kns, knm, tempmtx1, tempmtx2);

  std::shared_ptr<Core::LinAlg::SparseMatrix> kteffnew =
      std::make_shared<Core::LinAlg::SparseMatrix>(
          k->row_map(), 81, true, false, k->get_matrixtype());

  // build new stiffness matrix
  Core::LinAlg::matrix_add(*knn, false, 1.0, *kteffnew, 1.0);
  Core::LinAlg::matrix_add(*knm, false, 1.0, *kteffnew, 1.0);
  Core::LinAlg::matrix_add(*kmn, false, 1.0, *kteffnew, 1.0);
  Core::LinAlg::matrix_add(*kmm, false, 1.0, *kteffnew, 1.0);
  Core::LinAlg::matrix_add(
      *Core::LinAlg::matrix_multiply(*kns, false, *p_col, false, true, false, true), false, 1.,
      *kteffnew, 1.);
  Core::LinAlg::matrix_add(
      *Core::LinAlg::matrix_multiply(*p_row, true, *ksn, false, true, false, true), false, 1.,
      *kteffnew, 1.);
  Core::LinAlg::matrix_add(
      *Core::LinAlg::matrix_multiply(*kms, false, *p_col, false, true, false, true), false, 1.,
      *kteffnew, 1.);
  Core::LinAlg::matrix_add(
      *Core::LinAlg::matrix_multiply(*p_row, true, *ksm, false, true, false, true), false, 1.,
      *kteffnew, 1.);
  Core::LinAlg::matrix_add(
      *Core::LinAlg::matrix_multiply(*p_row, true,
          *Core::LinAlg::matrix_multiply(*kss, false, *p_col, false, true, false, true), false,
          true, false, true),
      false, 1., *kteffnew, 1.);

  if (p_row == p_col)
    Core::LinAlg::matrix_add(
        *Core::LinAlg::create_identity_matrix(*gsrow), false, 1., *kteffnew, 1.);

  kteffnew->complete(k->domain_map(), k->range_map());

  // return new matrix
  k = kteffnew;

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Mortar::Utils::mortar_rhs_condensation(
    Core::LinAlg::Vector<double>& rhs, Core::LinAlg::SparseMatrix& p)
{
  // prepare maps
  std::shared_ptr<Core::LinAlg::Map> gsdofrowmap = std::const_pointer_cast<Core::LinAlg::Map>(
      Core::Utils::shared_ptr_from_ref<const Core::LinAlg::Map>(p.range_map()));
  std::shared_ptr<Core::LinAlg::Map> gmdofrowmap = std::const_pointer_cast<Core::LinAlg::Map>(
      Core::Utils::shared_ptr_from_ref<const Core::LinAlg::Map>(p.domain_map()));

  Core::LinAlg::Vector<double> fs(*gsdofrowmap);
  Core::LinAlg::Vector<double> fm_cond(*gmdofrowmap);
  Core::LinAlg::export_to(rhs, fs);
  Core::LinAlg::Vector<double> fs_full(rhs.get_map());
  Core::LinAlg::export_to(fs, fs_full);
  rhs.update(-1., fs_full, 1.);

  p.multiply(true, fs, fm_cond);

  Core::LinAlg::Vector<double> fm_cond_full(rhs.get_map());
  Core::LinAlg::export_to(fm_cond, fm_cond_full);
  rhs.update(1., fm_cond_full, 1.);

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Mortar::Utils::mortar_recover(Core::LinAlg::Vector<double>& inc, Core::LinAlg::SparseMatrix& p)
{
  // prepare maps
  std::shared_ptr<Core::LinAlg::Map> gsdofrowmap = std::const_pointer_cast<Core::LinAlg::Map>(
      Core::Utils::shared_ptr_from_ref<const Core::LinAlg::Map>(p.range_map()));
  std::shared_ptr<Core::LinAlg::Map> gmdofrowmap = std::const_pointer_cast<Core::LinAlg::Map>(
      Core::Utils::shared_ptr_from_ref<const Core::LinAlg::Map>(p.domain_map()));

  Core::LinAlg::Vector<double> m_inc(*gmdofrowmap);
  Core::LinAlg::export_to(inc, m_inc);

  Core::LinAlg::Vector<double> s_inc(*gsdofrowmap);
  p.multiply(false, m_inc, s_inc);
  Core::LinAlg::Vector<double> s_inc_full(inc.get_map());
  Core::LinAlg::export_to(s_inc, s_inc_full);
  inc.update(1., s_inc_full, 1.);

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Mortar::Utils::mortar_matrix_condensation(
    std::shared_ptr<Core::LinAlg::BlockSparseMatrixBase>& k,
    const std::vector<std::shared_ptr<Core::LinAlg::SparseMatrix>>& p)
{
  std::shared_ptr<Core::LinAlg::BlockSparseMatrixBase> cond_mat =
      std::make_shared<Core::LinAlg::BlockSparseMatrix<Core::LinAlg::DefaultBlockMatrixStrategy>>(
          k->domain_extractor(), k->range_extractor(), 81, false, true);

  for (int row = 0; row < k->rows(); ++row)
    for (int col = 0; col < k->cols(); ++col)
    {
      std::shared_ptr<Core::LinAlg::SparseMatrix> new_matrix =
          std::make_shared<Core::LinAlg::SparseMatrix>(k->matrix(row, col));
      mortar_matrix_condensation(new_matrix, p.at(row), p.at(col) /*,row!=col*/);
      cond_mat->assign(row, col, Core::LinAlg::DataAccess::Copy, *new_matrix);
    }

  cond_mat->complete();

  k = cond_mat;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Mortar::Utils::mortar_rhs_condensation(Core::LinAlg::Vector<double>& rhs,
    const std::vector<std::shared_ptr<Core::LinAlg::SparseMatrix>>& p)
{
  for (unsigned i = 0; i < p.size(); mortar_rhs_condensation(rhs, *p[i++]));
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Mortar::Utils::mortar_recover(Core::LinAlg::Vector<double>& inc,
    const std::vector<std::shared_ptr<Core::LinAlg::SparseMatrix>>& p)
{
  for (unsigned i = 0; i < p.size(); mortar_recover(inc, *p[i++]));
}

FOUR_C_NAMESPACE_CLOSE
