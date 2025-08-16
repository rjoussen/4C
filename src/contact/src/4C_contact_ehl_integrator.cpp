// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_contact_ehl_integrator.hpp"

#include "4C_contact_element.hpp"
#include "4C_contact_nitsche_integrator.hpp"  // for CONTACT::Utils:: functions
#include "4C_contact_node.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::IntegratorEhl::integrate_gp_3d(Mortar::Element& sele, Mortar::Element& mele,
    Core::LinAlg::SerialDenseVector& sval, Core::LinAlg::SerialDenseVector& lmval,
    Core::LinAlg::SerialDenseVector& mval, Core::LinAlg::SerialDenseMatrix& sderiv,
    Core::LinAlg::SerialDenseMatrix& mderiv, Core::LinAlg::SerialDenseMatrix& lmderiv,
    Core::Gen::Pairedvector<int, Core::LinAlg::SerialDenseMatrix>& dualmap, double& wgt,
    double& jac, Core::Gen::Pairedvector<int, double>& derivjac, double* normal,
    std::vector<Core::Gen::Pairedvector<int, double>>& dnmap_unit, double& gap,
    Core::Gen::Pairedvector<int, double>& deriv_gap, double* sxi, double* mxi,
    std::vector<Core::Gen::Pairedvector<int, double>>& derivsxi,
    std::vector<Core::Gen::Pairedvector<int, double>>& derivmxi)
{
  // check bound
  bool bound = false;
  for (int i = 0; i < sele.num_node(); ++i)
    if (dynamic_cast<CONTACT::Node*>(sele.nodes()[i])->is_on_boundor_ce())
      FOUR_C_THROW("no boundary modification for EHL implemented");

  // is quadratic case?
  bool quad = sele.is_quad();

  // weighted gap
  gp_3d_w_gap(sele, sval, lmval, &gap, jac, wgt, quad);
  for (int j = 0; j < sele.num_node(); ++j)
    gp_g_lin(j, sele, mele, sval, mval, lmval, sderiv, lmderiv, gap, normal, jac, wgt, deriv_gap,
        derivjac, derivsxi, dualmap);

  // integrate D and M matrix
  gp_dm(sele, mele, lmval, sval, mval, jac, wgt, bound);
  gp_3d_dm_lin(sele, mele, sval, mval, lmval, sderiv, mderiv, lmderiv, wgt, jac, derivsxi, derivmxi,
      derivjac, dualmap);

  // get second derivative of shape function
  Core::LinAlg::SerialDenseMatrix ssecderiv(sele.num_node(), 3);
  sele.evaluate2nd_deriv_shape(sxi, ssecderiv, sele.num_node());

  // weighted surface gradient
  gp_weighted_surf_grad_and_deriv(
      sele, sxi, derivsxi, lmval, lmderiv, dualmap, sval, sderiv, ssecderiv, wgt, jac, derivjac);

  //  // weighted tangential velocity (average and relative)
  gp_weighted_av_rel_vel(sele, mele, sval, lmval, mval, sderiv, mderiv, lmderiv, dualmap, wgt, jac,
      derivjac, normal, dnmap_unit, gap, deriv_gap, sxi, mxi, derivsxi, derivmxi);
}



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::IntegratorEhl::integrate_gp_2d(Mortar::Element& sele, Mortar::Element& mele,
    Core::LinAlg::SerialDenseVector& sval, Core::LinAlg::SerialDenseVector& lmval,
    Core::LinAlg::SerialDenseVector& mval, Core::LinAlg::SerialDenseMatrix& sderiv,
    Core::LinAlg::SerialDenseMatrix& mderiv, Core::LinAlg::SerialDenseMatrix& lmderiv,
    Core::Gen::Pairedvector<int, Core::LinAlg::SerialDenseMatrix>& dualmap, double& wgt,
    double& jac, Core::Gen::Pairedvector<int, double>& derivjac, double* normal,
    std::vector<Core::Gen::Pairedvector<int, double>>& dnmap_unit, double& gap,
    Core::Gen::Pairedvector<int, double>& deriv_gap, double* sxi, double* mxi,
    std::vector<Core::Gen::Pairedvector<int, double>>& derivsxi,
    std::vector<Core::Gen::Pairedvector<int, double>>& derivmxi)
{
  FOUR_C_THROW("2D EHL integration not supported");
}


void CONTACT::IntegratorEhl::gp_weighted_surf_grad_and_deriv(Mortar::Element& sele,
    const double* xi, const std::vector<Core::Gen::Pairedvector<int, double>>& dsxigp,
    const Core::LinAlg::SerialDenseVector& lmval, const Core::LinAlg::SerialDenseMatrix& lmderiv,
    const Core::Gen::Pairedvector<int, Core::LinAlg::SerialDenseMatrix>& dualmap,
    const Core::LinAlg::SerialDenseVector& sval, const Core::LinAlg::SerialDenseMatrix& sderiv,
    const Core::LinAlg::SerialDenseMatrix& sderiv2, const double& wgt, const double& jac,
    const Core::Gen::Pairedvector<int, double>& jacintcellmap)
{
  // empty local basis vectors
  std::vector<std::vector<double>> gxi(2, std::vector<double>(3, 0));

  // metrics routine gives local basis vectors
  sele.metrics(xi, gxi.at(0).data(), gxi.at(1).data());

  if (n_dim() == 2)
  {
    gxi.at(1).at(0) = gxi.at(1).at(1) = 0.;
    gxi.at(1).at(2) = 1.;
  }

  Core::LinAlg::Matrix<3, 1> t1(gxi.at(0).data(), true);
  Core::LinAlg::Matrix<3, 1> t2(gxi.at(1).data(), true);
  Core::LinAlg::Matrix<3, 1> n;
  n.cross_product(t1, t2);
  n.scale(1. / n.norm2());
  Core::LinAlg::Matrix<3, 3> covariant_metric;
  for (int i = 0; i < 3; ++i)
  {
    covariant_metric(i, 0) = t1(i);
    covariant_metric(i, 1) = t2(i);
    covariant_metric(i, 2) = n(i);
  }
  Core::LinAlg::Matrix<3, 3> contravariant_metric;
  contravariant_metric.invert(covariant_metric);

  std::vector<std::vector<double>> gxi_contra(2, std::vector<double>(3, 0));
  for (int i = 0; i < n_dim() - 1; ++i)
    for (int d = 0; d < n_dim(); ++d) gxi_contra.at(i).at(d) = contravariant_metric(i, d);

  for (int a = 0; a < sele.num_node(); ++a)
  {
    Core::Nodes::Node* node = sele.nodes()[a];
    CONTACT::Node* cnode = dynamic_cast<CONTACT::Node*>(node);
    if (!cnode) FOUR_C_THROW("this is not a contact node");

    for (int c = 0; c < sele.num_node(); ++c)
    {
      Core::LinAlg::Matrix<3, 1>& tmp =
          cnode->ehl_data()
              .get_surf_grad()[dynamic_cast<CONTACT::Node*>(sele.nodes()[c])->dofs()[0]];
      for (int d = 0; d < n_dim(); ++d)
        for (int al = 0; al < n_dim() - 1; ++al)
          tmp(d) += wgt * jac * lmval(a) * sderiv(c, al) * gxi_contra.at(al).at(d);
    }

    for (int c = 0; c < sele.num_node(); ++c)
    {
      for (auto p = jacintcellmap.begin(); p != jacintcellmap.end(); ++p)
      {
        Core::LinAlg::Matrix<3, 1>& tmp = cnode->ehl_data().get_surf_grad_deriv()[p
                ->first][dynamic_cast<CONTACT::Node*>(sele.nodes()[c])->dofs()[0]];
        for (int d = 0; d < n_dim(); ++d)
          for (int al = 0; al < n_dim() - 1; ++al)
            tmp(d) += wgt * p->second * lmval(a) * sderiv(c, al) * gxi_contra.at(al).at(d);
      }
      for (int e = 0; e < n_dim() - 1; ++e)
        for (auto p = dsxigp.at(e).begin(); p != dsxigp.at(e).end(); ++p)
        {
          Core::LinAlg::Matrix<3, 1>& tmp = cnode->ehl_data().get_surf_grad_deriv()[p
                  ->first][dynamic_cast<CONTACT::Node*>(sele.nodes()[c])->dofs()[0]];
          for (int d = 0; d < n_dim(); ++d)
            for (int al = 0; al < n_dim() - 1; ++al)
              tmp(d) +=
                  wgt * jac * lmderiv(a, e) * p->second * sderiv(c, al) * gxi_contra.at(al).at(d);
        }

      for (auto p = dualmap.begin(); p != dualmap.end(); ++p)
      {
        Core::LinAlg::Matrix<3, 1>& tmp = cnode->ehl_data().get_surf_grad_deriv()[p
                ->first][dynamic_cast<CONTACT::Node*>(sele.nodes()[c])->dofs()[0]];
        for (int d = 0; d < n_dim(); ++d)
          for (int al = 0; al < n_dim() - 1; ++al)
            for (int m = 0; m < sele.num_node(); ++m)
              tmp(d) +=
                  wgt * jac * p->second(a, m) * sval(m) * sderiv(c, al) * gxi_contra.at(al).at(d);
      }
    }
  }
}

void CONTACT::IntegratorEhl::gp_weighted_av_rel_vel(Mortar::Element& sele, Mortar::Element& mele,
    const Core::LinAlg::SerialDenseVector& sval, const Core::LinAlg::SerialDenseVector& lmval,
    const Core::LinAlg::SerialDenseVector& mval, const Core::LinAlg::SerialDenseMatrix& sderiv,
    const Core::LinAlg::SerialDenseMatrix& mderiv, const Core::LinAlg::SerialDenseMatrix& lmderiv,
    const Core::Gen::Pairedvector<int, Core::LinAlg::SerialDenseMatrix>& dualmap, const double& wgt,
    const double& jac, const Core::Gen::Pairedvector<int, double>& derivjac, const double* normal,
    const std::vector<Core::Gen::Pairedvector<int, double>>& dnmap_unit, const double& gap,
    const Core::Gen::Pairedvector<int, double>& deriv_gap, const double* sxi, const double* mxi,
    const std::vector<Core::Gen::Pairedvector<int, double>>& derivsxi,
    const std::vector<Core::Gen::Pairedvector<int, double>>& derivmxi)
{
  constexpr int dim = 3;
  if (IntegratorEhl::n_dim() != 3)
    FOUR_C_THROW("dimension inconsistency, or is this not implemented for all spatial dimensions?");

  Core::LinAlg::Matrix<dim, 1> t1, t2;
  std::vector<Core::Gen::Pairedvector<int, double>> dt1, dt2;
  Core::LinAlg::Matrix<dim, 1> relVel;
  std::vector<Core::Gen::Pairedvector<int, double>> relVel_deriv(
      dim, sele.num_node() * dim + mele.num_node() * dim + derivsxi[0].size() + derivmxi[0].size());
  double vt1, vt2;
  Core::Gen::Pairedvector<int, double> dvt1(0);
  Core::Gen::Pairedvector<int, double> dvt2(0);

  CONTACT::Utils::build_tangent_vectors<dim>(normal, dnmap_unit, t1.data(), dt1, t2.data(), dt2);
  CONTACT::Utils::rel_vel_invariant<dim>(sele, sxi, derivsxi, sval, sderiv, mele, mxi, derivmxi,
      mval, mderiv, gap, deriv_gap, relVel, relVel_deriv, -.5);

  CONTACT::Utils::vector_scalar_product<dim>(t1, dt1, relVel, relVel_deriv, vt1, dvt1);
  CONTACT::Utils::vector_scalar_product<dim>(t2, dt2, relVel, relVel_deriv, vt2, dvt2);

  for (int i = 0; i < sele.num_node(); ++i)
  {
    CONTACT::Node* cnode = dynamic_cast<CONTACT::Node*>(sele.nodes()[i]);
    for (int d = 0; d < dim; ++d)
      cnode->ehl_data().get_weighted_rel_tang_vel()(d) +=
          jac * wgt * lmval(i) * (vt1 * t1(d) + vt2 * t2(d));

    for (auto p = derivjac.begin(); p != derivjac.end(); ++p)
    {
      Core::LinAlg::Matrix<3, 1>& tmp =
          cnode->ehl_data().get_weighted_rel_tang_vel_deriv()[p->first];
      for (int d = 0; d < dim; ++d)
        tmp(d) += p->second * wgt * lmval(i) * (vt1 * t1(d) + vt2 * t2(d));
    }

    for (int e = 0; e < dim - 1; ++e)
      for (auto p = derivsxi.at(e).begin(); p != derivsxi.at(e).end(); ++p)
      {
        Core::LinAlg::Matrix<3, 1>& tmp =
            cnode->ehl_data().get_weighted_rel_tang_vel_deriv()[p->first];
        for (int d = 0; d < dim; ++d)
          tmp(d) += jac * wgt * lmderiv(i, e) * p->second * (vt1 * t1(d) + vt2 * t2(d));
      }

    for (auto p = dualmap.begin(); p != dualmap.end(); ++p)
    {
      Core::LinAlg::Matrix<3, 1>& tmp =
          cnode->ehl_data().get_weighted_rel_tang_vel_deriv()[p->first];
      for (int d = 0; d < dim; ++d)
        for (int m = 0; m < sele.num_node(); ++m)
          tmp(d) += jac * wgt * p->second(i, m) * sval(m) * (vt1 * t1(d) + vt2 * t2(d));
    }

    for (auto p = dvt1.begin(); p != dvt1.end(); ++p)
    {
      Core::LinAlg::Matrix<3, 1>& tmp =
          cnode->ehl_data().get_weighted_rel_tang_vel_deriv()[p->first];
      for (int d = 0; d < dim; ++d) tmp(d) += jac * wgt * lmval(i) * p->second * t1(d);
    }

    for (int d = 0; d < dim; ++d)
      for (auto p = dt1.at(d).begin(); p != dt1.at(d).end(); ++p)
        cnode->ehl_data().get_weighted_rel_tang_vel_deriv()[p->first](d) +=
            jac * wgt * lmval(i) * vt1 * p->second;

    for (auto p = dvt2.begin(); p != dvt2.end(); ++p)
    {
      Core::LinAlg::Matrix<3, 1>& tmp =
          cnode->ehl_data().get_weighted_rel_tang_vel_deriv()[p->first];
      for (int d = 0; d < dim; ++d) tmp(d) += jac * wgt * lmval(i) * p->second * t2(d);
    }

    for (int d = 0; d < dim; ++d)
      for (auto p = dt2.at(d).begin(); p != dt2.at(d).end(); ++p)
        cnode->ehl_data().get_weighted_rel_tang_vel_deriv()[p->first](d) +=
            jac * wgt * lmval(i) * vt2 * p->second;
  }

  CONTACT::Utils::rel_vel<dim>(sele, sval, sderiv, derivsxi, -1., relVel, relVel_deriv);
  CONTACT::Utils::vector_scalar_product<dim>(t1, dt1, relVel, relVel_deriv, vt1, dvt1);
  CONTACT::Utils::vector_scalar_product<dim>(t2, dt2, relVel, relVel_deriv, vt2, dvt2);

  for (int i = 0; i < sele.num_node(); ++i)
  {
    CONTACT::Node* cnode = dynamic_cast<CONTACT::Node*>(sele.nodes()[i]);
    for (int d = 0; d < dim; ++d)
      cnode->ehl_data().get_weighted_av_tang_vel()(d) -=
          jac * wgt * lmval(i) * (vt1 * t1(d) + vt2 * t2(d));

    for (auto p = derivjac.begin(); p != derivjac.end(); ++p)
    {
      Core::LinAlg::Matrix<3, 1>& tmp =
          cnode->ehl_data().get_weighted_av_tang_vel_deriv()[p->first];
      for (int d = 0; d < dim; ++d)
        tmp(d) -= p->second * wgt * lmval(i) * (vt1 * t1(d) + vt2 * t2(d));
    }

    for (int e = 0; e < dim - 1; ++e)
      for (auto p = derivsxi.at(e).begin(); p != derivsxi.at(e).end(); ++p)
      {
        Core::LinAlg::Matrix<3, 1>& tmp =
            cnode->ehl_data().get_weighted_av_tang_vel_deriv()[p->first];
        for (int d = 0; d < dim; ++d)
          tmp(d) -= jac * wgt * lmderiv(i, e) * p->second * (vt1 * t1(d) + vt2 * t2(d));
      }

    for (auto p = dualmap.begin(); p != dualmap.end(); ++p)
    {
      Core::LinAlg::Matrix<3, 1>& tmp =
          cnode->ehl_data().get_weighted_av_tang_vel_deriv()[p->first];
      for (int d = 0; d < dim; ++d)
        for (int m = 0; m < sele.num_node(); ++m)
          tmp(d) -= jac * wgt * p->second(i, m) * sval(m) * (vt1 * t1(d) + vt2 * t2(d));
    }

    for (auto p = dvt1.begin(); p != dvt1.end(); ++p)
    {
      Core::LinAlg::Matrix<3, 1>& tmp =
          cnode->ehl_data().get_weighted_av_tang_vel_deriv()[p->first];
      for (int d = 0; d < dim; ++d) tmp(d) -= jac * wgt * lmval(i) * p->second * t1(d);
    }

    for (int d = 0; d < dim; ++d)
      for (auto p = dt1.at(d).begin(); p != dt1.at(d).end(); ++p)
        cnode->ehl_data().get_weighted_av_tang_vel_deriv()[p->first](d) -=
            jac * wgt * lmval(i) * vt1 * p->second;

    for (auto p = dvt2.begin(); p != dvt2.end(); ++p)
    {
      Core::LinAlg::Matrix<3, 1>& tmp =
          cnode->ehl_data().get_weighted_av_tang_vel_deriv()[p->first];
      for (int d = 0; d < dim; ++d) tmp(d) -= jac * wgt * lmval(i) * p->second * t2(d);
    }

    for (int d = 0; d < dim; ++d)
      for (auto p = dt2.at(d).begin(); p != dt2.at(d).end(); ++p)
        cnode->ehl_data().get_weighted_av_tang_vel_deriv()[p->first](d) -=
            jac * wgt * lmval(i) * vt2 * p->second;
  }
}

FOUR_C_NAMESPACE_CLOSE
