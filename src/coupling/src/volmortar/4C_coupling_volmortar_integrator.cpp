// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_coupling_volmortar_integrator.hpp"

#include "4C_comm_mpi_utils.hpp"
#include "4C_coupling_volmortar_cell.hpp"
#include "4C_coupling_volmortar_defines.hpp"
#include "4C_coupling_volmortar_shape.hpp"
#include "4C_fem_discretization.hpp"
#include "4C_fem_general_utils_integration.hpp"
#include "4C_linalg_serialdensematrix.hpp"
#include "4C_linalg_serialdensevector.hpp"
#include "4C_linalg_sparsematrix.hpp"
#include "4C_mortar_calc_utils.hpp"
#include "4C_mortar_coupling3d_classes.hpp"

#include <Teuchos_StandardParameterEntryValidators.hpp>

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*
 |  ctor (public)                                            farah 02/15|
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype_s>
Coupling::VolMortar::VolMortarIntegratorEleBased<distype_s>::VolMortarIntegratorEleBased(
    Teuchos::ParameterList& params)
{
  // get type of quadratic modification
  dualquad_ = Teuchos::getIntegralValue<DualQuad>(params, "DUALQUAD");

  // get type of quadratic modification
  shape_ = Teuchos::getIntegralValue<Shapefcn>(params, "SHAPEFCN");
}

/*----------------------------------------------------------------------*
 |  Initialize gauss points for ele-based integration        farah 02/15|
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype_s>
void Coupling::VolMortar::VolMortarIntegratorEleBased<distype_s>::initialize_gp()
{
  // init shape of integration domain
  Core::FE::CellType intshape = distype_s;

  //*******************************
  // choose Gauss rule accordingly
  //*******************************
  switch (intshape)
  {
    //*******************************
    //               2D
    //*******************************
    case Core::FE::CellType::tri3:
    {
      Core::FE::GaussRule2D mygaussrule = Core::FE::GaussRule2D::tri_7point;

      const Core::FE::IntegrationPoints2D intpoints(mygaussrule);
      ngp_ = intpoints.nquad;
      coords_.reshape(ngp_, 2);
      weights_.resize(ngp_);
      for (int i = 0; i < ngp_; ++i)
      {
        coords_(i, 0) = intpoints.qxg[i][0];
        coords_(i, 1) = intpoints.qxg[i][1];
        weights_[i] = intpoints.qwgt[i];
      }
      break;
    }
    case Core::FE::CellType::tri6:
    {
      Core::FE::GaussRule2D mygaussrule = Core::FE::GaussRule2D::tri_12point;

      const Core::FE::IntegrationPoints2D intpoints(mygaussrule);
      ngp_ = intpoints.nquad;
      coords_.reshape(ngp_, 2);
      weights_.resize(ngp_);
      for (int i = 0; i < ngp_; ++i)
      {
        coords_(i, 0) = intpoints.qxg[i][0];
        coords_(i, 1) = intpoints.qxg[i][1];
        weights_[i] = intpoints.qwgt[i];
      }
      break;
    }
    case Core::FE::CellType::quad4:
    {
      Core::FE::GaussRule2D mygaussrule = Core::FE::GaussRule2D::quad_64point;

      const Core::FE::IntegrationPoints2D intpoints(mygaussrule);
      ngp_ = intpoints.nquad;
      coords_.reshape(ngp_, 2);
      weights_.resize(ngp_);
      for (int i = 0; i < ngp_; ++i)
      {
        coords_(i, 0) = intpoints.qxg[i][0];
        coords_(i, 1) = intpoints.qxg[i][1];
        weights_[i] = intpoints.qwgt[i];
      }
      break;
    }
    case Core::FE::CellType::quad8:
    {
      Core::FE::GaussRule2D mygaussrule = Core::FE::GaussRule2D::quad_64point;

      const Core::FE::IntegrationPoints2D intpoints(mygaussrule);
      ngp_ = intpoints.nquad;
      coords_.reshape(ngp_, 2);
      weights_.resize(ngp_);
      for (int i = 0; i < ngp_; ++i)
      {
        coords_(i, 0) = intpoints.qxg[i][0];
        coords_(i, 1) = intpoints.qxg[i][1];
        weights_[i] = intpoints.qwgt[i];
      }
      break;
    }
    case Core::FE::CellType::quad9:
    {
      Core::FE::GaussRule2D mygaussrule = Core::FE::GaussRule2D::quad_64point;

      const Core::FE::IntegrationPoints2D intpoints(mygaussrule);
      ngp_ = intpoints.nquad;
      coords_.reshape(ngp_, 2);
      weights_.resize(ngp_);
      for (int i = 0; i < ngp_; ++i)
      {
        coords_(i, 0) = intpoints.qxg[i][0];
        coords_(i, 1) = intpoints.qxg[i][1];
        weights_[i] = intpoints.qwgt[i];
      }
      break;
    }
    //*******************************
    //               3D
    //*******************************
    case Core::FE::CellType::tet4:
    {
      Core::FE::GaussRule3D mygaussrule = Core::FE::GaussRule3D::tet_45point;

      const Core::FE::IntegrationPoints3D intpoints(mygaussrule);
      ngp_ = intpoints.nquad;
      coords_.reshape(ngp_, 3);
      weights_.resize(ngp_);
      for (int i = 0; i < ngp_; ++i)
      {
        coords_(i, 0) = intpoints.qxg[i][0];
        coords_(i, 1) = intpoints.qxg[i][1];
        coords_(i, 2) = intpoints.qxg[i][2];
        weights_[i] = intpoints.qwgt[i];
      }
      break;
    }
    case Core::FE::CellType::tet10:
    {
      Core::FE::GaussRule3D mygaussrule = Core::FE::GaussRule3D::tet_45point;

      const Core::FE::IntegrationPoints3D intpoints(mygaussrule);
      ngp_ = intpoints.nquad;
      coords_.reshape(ngp_, 3);
      weights_.resize(ngp_);
      for (int i = 0; i < ngp_; ++i)
      {
        coords_(i, 0) = intpoints.qxg[i][0];
        coords_(i, 1) = intpoints.qxg[i][1];
        coords_(i, 2) = intpoints.qxg[i][2];
        weights_[i] = intpoints.qwgt[i];
      }
      break;
    }
    case Core::FE::CellType::hex8:
    {
      Core::FE::GaussRule3D mygaussrule = Core::FE::GaussRule3D::hex_27point;

      const Core::FE::IntegrationPoints3D intpoints(mygaussrule);
      ngp_ = intpoints.nquad;
      coords_.reshape(ngp_, 3);
      weights_.resize(ngp_);
      for (int i = 0; i < ngp_; ++i)
      {
        coords_(i, 0) = intpoints.qxg[i][0];
        coords_(i, 1) = intpoints.qxg[i][1];
        coords_(i, 2) = intpoints.qxg[i][2];
        weights_[i] = intpoints.qwgt[i];
      }
      break;
    }
    case Core::FE::CellType::hex20:
    {
      Core::FE::GaussRule3D mygaussrule = Core::FE::GaussRule3D::hex_125point;

      const Core::FE::IntegrationPoints3D intpoints(mygaussrule);
      ngp_ = intpoints.nquad;
      coords_.reshape(ngp_, 3);
      weights_.resize(ngp_);
      for (int i = 0; i < ngp_; ++i)
      {
        coords_(i, 0) = intpoints.qxg[i][0];
        coords_(i, 1) = intpoints.qxg[i][1];
        coords_(i, 2) = intpoints.qxg[i][2];
        weights_[i] = intpoints.qwgt[i];
      }
      break;
    }
    case Core::FE::CellType::hex27:
    {
      Core::FE::GaussRule3D mygaussrule = Core::FE::GaussRule3D::hex_125point;

      const Core::FE::IntegrationPoints3D intpoints(mygaussrule);
      ngp_ = intpoints.nquad;
      coords_.reshape(ngp_, 3);
      weights_.resize(ngp_);
      for (int i = 0; i < ngp_; ++i)
      {
        coords_(i, 0) = intpoints.qxg[i][0];
        coords_(i, 1) = intpoints.qxg[i][1];
        coords_(i, 2) = intpoints.qxg[i][2];
        weights_[i] = intpoints.qwgt[i];
      }
      break;
    }
    case Core::FE::CellType::pyramid5:
    {
      Core::FE::GaussRule3D mygaussrule = Core::FE::GaussRule3D::pyramid_8point;

      const Core::FE::IntegrationPoints3D intpoints(mygaussrule);
      ngp_ = intpoints.nquad;
      coords_.reshape(ngp_, 3);
      weights_.resize(ngp_);
      for (int i = 0; i < ngp_; ++i)
      {
        coords_(i, 0) = intpoints.qxg[i][0];
        coords_(i, 1) = intpoints.qxg[i][1];
        coords_(i, 2) = intpoints.qxg[i][2];
        weights_[i] = intpoints.qwgt[i];
      }
      break;
    }
    //*******************************
    //            Default
    //*******************************
    default:
    {
      FOUR_C_THROW("ERROR: VolMortarIntegrator: This element type is not implemented!");
      break;
    }
  }  // switch(eletype)

  return;
}

/*----------------------------------------------------------------------*
 |  Initialize gauss points for ele-based integration        farah 02/15|
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype_s>
void Coupling::VolMortar::VolMortarIntegratorEleBased<distype_s>::integrate_ele_based_3d(
    Core::Elements::Element& sele, std::vector<int>& foundeles, Core::LinAlg::SparseMatrix& D,
    Core::LinAlg::SparseMatrix& M, const Core::FE::Discretization& Adis,
    const Core::FE::Discretization& Bdis, int dofseta, int dofsetb,
    const Core::LinAlg::Map& PAB_dofrowmap, const Core::LinAlg::Map& PAB_dofcolmap)
{
  //**********************************************************************
  // loop over all Gauss points for integration
  //**********************************************************************
  for (int gp = 0; gp < ngp_; ++gp)
  {
    // coordinates and weight
    double eta[3] = {coords_(gp, 0), coords_(gp, 1), 0.0};

    if (ndim_ == 3) eta[2] = coords_(gp, 2);

    double wgt = weights_[gp];
    double jac = 0.0;
    double globgp[3] = {0.0, 0.0, 0.0};

    // quantities for eval. outside gp
    double gpdist = 1.0e12;
    int gpid = 0;
    double AuxXi[3] = {0.0, 0.0, 0.0};

    // evaluate the integration cell Jacobian
    jac = Utils::jacobian<distype_s>(eta, sele);

    // get global Gauss point coordinates
    Utils::local_to_global<distype_s>(sele, eta, globgp);

    // map gp into A and B para space
    double Axi[3] = {0.0, 0.0, 0.0};
    Mortar::Utils::global_to_local<distype_s>(sele, globgp, Axi);

    // loop over beles
    for (int found = 0; found < (int)foundeles.size(); ++found)
    {
      // get master element
      Core::Elements::Element* Bele = Bdis.g_element(foundeles[found]);
      Core::FE::CellType shape = Bele->shape();

      bool proj = false;

      switch (shape)
      {
        //************************************************
        //                    2D
        //************************************************
        case Core::FE::CellType::tri3:
        {
          proj = vol_mortar_ele_based_gp<distype_s, Core::FE::CellType::tri3>(sele, Bele, foundeles,
              found, gpid, jac, wgt, gpdist, Axi, AuxXi, globgp, dualquad_, shape_, D, M, Adis,
              Bdis, dofseta, dofsetb, PAB_dofrowmap, PAB_dofcolmap);

          break;
        }
        case Core::FE::CellType::tri6:
        {
          proj = vol_mortar_ele_based_gp<distype_s, Core::FE::CellType::tri6>(sele, Bele, foundeles,
              found, gpid, jac, wgt, gpdist, Axi, AuxXi, globgp, dualquad_, shape_, D, M, Adis,
              Bdis, dofseta, dofsetb, PAB_dofrowmap, PAB_dofcolmap);

          break;
        }
        case Core::FE::CellType::quad4:
        {
          proj = vol_mortar_ele_based_gp<distype_s, Core::FE::CellType::quad4>(sele, Bele,
              foundeles, found, gpid, jac, wgt, gpdist, Axi, AuxXi, globgp, dualquad_, shape_, D, M,
              Adis, Bdis, dofseta, dofsetb, PAB_dofrowmap, PAB_dofcolmap);

          break;
        }
        case Core::FE::CellType::quad8:
        {
          proj = vol_mortar_ele_based_gp<distype_s, Core::FE::CellType::quad8>(sele, Bele,
              foundeles, found, gpid, jac, wgt, gpdist, Axi, AuxXi, globgp, dualquad_, shape_, D, M,
              Adis, Bdis, dofseta, dofsetb, PAB_dofrowmap, PAB_dofcolmap);

          break;
        }
        case Core::FE::CellType::quad9:
        {
          proj = vol_mortar_ele_based_gp<distype_s, Core::FE::CellType::quad9>(sele, Bele,
              foundeles, found, gpid, jac, wgt, gpdist, Axi, AuxXi, globgp, dualquad_, shape_, D, M,
              Adis, Bdis, dofseta, dofsetb, PAB_dofrowmap, PAB_dofcolmap);
          break;
        }
        //************************************************
        //                    3D
        //************************************************
        case Core::FE::CellType::hex8:
        {
          proj = vol_mortar_ele_based_gp<distype_s, Core::FE::CellType::hex8>(sele, Bele, foundeles,
              found, gpid, jac, wgt, gpdist, Axi, AuxXi, globgp, dualquad_, shape_, D, M, Adis,
              Bdis, dofseta, dofsetb, PAB_dofrowmap, PAB_dofcolmap);

          break;
        }
        case Core::FE::CellType::hex20:
        {
          proj = vol_mortar_ele_based_gp<distype_s, Core::FE::CellType::hex20>(sele, Bele,
              foundeles, found, gpid, jac, wgt, gpdist, Axi, AuxXi, globgp, dualquad_, shape_, D, M,
              Adis, Bdis, dofseta, dofsetb, PAB_dofrowmap, PAB_dofcolmap);

          break;
        }
        case Core::FE::CellType::hex27:
        {
          proj = vol_mortar_ele_based_gp<distype_s, Core::FE::CellType::hex27>(sele, Bele,
              foundeles, found, gpid, jac, wgt, gpdist, Axi, AuxXi, globgp, dualquad_, shape_, D, M,
              Adis, Bdis, dofseta, dofsetb, PAB_dofrowmap, PAB_dofcolmap);

          break;
        }
        case Core::FE::CellType::tet4:
        {
          proj = vol_mortar_ele_based_gp<distype_s, Core::FE::CellType::tet4>(sele, Bele, foundeles,
              found, gpid, jac, wgt, gpdist, Axi, AuxXi, globgp, dualquad_, shape_, D, M, Adis,
              Bdis, dofseta, dofsetb, PAB_dofrowmap, PAB_dofcolmap);

          break;
        }
        case Core::FE::CellType::tet10:
        {
          proj = vol_mortar_ele_based_gp<distype_s, Core::FE::CellType::tet10>(sele, Bele,
              foundeles, found, gpid, jac, wgt, gpdist, Axi, AuxXi, globgp, dualquad_, shape_, D, M,
              Adis, Bdis, dofseta, dofsetb, PAB_dofrowmap, PAB_dofcolmap);

          break;
        }
        case Core::FE::CellType::pyramid5:
        {
          proj = vol_mortar_ele_based_gp<distype_s, Core::FE::CellType::pyramid5>(sele, Bele,
              foundeles, found, gpid, jac, wgt, gpdist, Axi, AuxXi, globgp, dualquad_, shape_, D, M,
              Adis, Bdis, dofseta, dofsetb, PAB_dofrowmap, PAB_dofcolmap);

          break;
        }
        default:
        {
          FOUR_C_THROW("ERROR: unknown shape!");
          break;
        }
      }

      // if gp evaluated break ele loop
      if (proj == true)
        break;
      else
        continue;
    }  // beles
  }  // end gp loop
  return;
}


/*----------------------------------------------------------------------*
 |  possible elements for ele-based integration              farah 02/15|
 *----------------------------------------------------------------------*/
template class Coupling::VolMortar::VolMortarIntegratorEleBased<Core::FE::CellType::quad4>;
template class Coupling::VolMortar::VolMortarIntegratorEleBased<Core::FE::CellType::quad8>;
template class Coupling::VolMortar::VolMortarIntegratorEleBased<Core::FE::CellType::quad9>;

template class Coupling::VolMortar::VolMortarIntegratorEleBased<Core::FE::CellType::tri3>;
template class Coupling::VolMortar::VolMortarIntegratorEleBased<Core::FE::CellType::tri6>;

template class Coupling::VolMortar::VolMortarIntegratorEleBased<Core::FE::CellType::hex8>;
template class Coupling::VolMortar::VolMortarIntegratorEleBased<Core::FE::CellType::hex20>;
template class Coupling::VolMortar::VolMortarIntegratorEleBased<Core::FE::CellType::hex27>;

template class Coupling::VolMortar::VolMortarIntegratorEleBased<Core::FE::CellType::tet4>;
template class Coupling::VolMortar::VolMortarIntegratorEleBased<Core::FE::CellType::tet10>;

template class Coupling::VolMortar::VolMortarIntegratorEleBased<Core::FE::CellType::pyramid5>;


/*----------------------------------------------------------------------*
 |  gp evaluation                                            farah 02/15|
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype_s, Core::FE::CellType distype_m>
bool Coupling::VolMortar::vol_mortar_ele_based_gp(Core::Elements::Element& sele,
    Core::Elements::Element* mele, std::vector<int>& foundeles, int& found, int& gpid, double& jac,
    double& wgt, double& gpdist, double* Axi, double* AuxXi, double* globgp, DualQuad& dq,
    Shapefcn& shape, Core::LinAlg::SparseMatrix& D, Core::LinAlg::SparseMatrix& M,
    const Core::FE::Discretization& Adis, const Core::FE::Discretization& Bdis, int dofseta,
    int dofsetb, const Core::LinAlg::Map& PAB_dofrowmap, const Core::LinAlg::Map& PAB_dofcolmap)
{
  //! ns_: number of slave element nodes
  static const int ns_ = Core::FE::num_nodes(distype_s);

  //! nm_: number of master element nodes
  static const int nm_ = Core::FE::num_nodes(distype_m);

  // create empty vectors for shape fct. evaluation
  Core::LinAlg::Matrix<ns_, 1> sval_A;
  Core::LinAlg::Matrix<nm_, 1> mval_A;
  Core::LinAlg::Matrix<ns_, 1> lmval_A;

  double Bxi[3] = {0.0, 0.0, 0.0};

  bool converged = true;
  Mortar::Utils::global_to_local<distype_m>(*mele, globgp, Bxi, converged);
  if (!converged and found != ((int)foundeles.size() - 1)) return false;

  // save distance of gp
  double l = sqrt(Bxi[0] * Bxi[0] + Bxi[1] * Bxi[1] + Bxi[2] * Bxi[2]);
  if (l < gpdist)
  {
    gpdist = l;
    gpid = foundeles[found];
    AuxXi[0] = Bxi[0];
    AuxXi[1] = Bxi[1];
    AuxXi[2] = Bxi[2];
  }

  // Check parameter space mapping
  bool proj = check_mapping<distype_s, distype_m>(sele, *mele, Axi, Bxi);

  // if gp outside continue or eval nearest gp
  if (!proj and (found != ((int)foundeles.size() - 1)))
    return false;
  else if (!proj and found == ((int)foundeles.size() - 1))
  {
    Bxi[0] = AuxXi[0];
    Bxi[1] = AuxXi[1];
    Bxi[2] = AuxXi[2];
    mele = Bdis.g_element(gpid);
  }

  // for "master" side
  Utils::shape_function<distype_s>(sval_A, Axi, dq);
  Utils::shape_function<distype_m>(mval_A, Bxi);

  // evaluate Lagrange multiplier shape functions (on slave element)
  Utils::dual_shape_function<distype_s>(lmval_A, Axi, sele, dq);

  // compute cell D/M matrix ****************************************
  // dual shape functions
  for (int j = 0; j < ns_; ++j)
  {
    Core::Nodes::Node* cnode = sele.nodes()[j];
    if (cnode->owner() != Core::Communication::my_mpi_rank(Adis.get_comm())) continue;

    const int nsdof = Adis.num_dof(dofseta, cnode);

    if (shape == shape_std)
    {
      for (int j = 0; j < ns_; ++j)
      {
        Core::Nodes::Node* cnode = sele.nodes()[j];
        int nsdof = Adis.num_dof(dofseta, cnode);

        // loop over slave dofs
        for (int jdof = 0; jdof < nsdof; ++jdof)
        {
          int row = Adis.dof(dofseta, cnode, jdof);

          // integrate M
          for (int k = 0; k < nm_; ++k)
          {
            Core::Nodes::Node* mnode = mele->nodes()[k];
            int nmdof = Bdis.num_dof(dofsetb, mnode);

            for (int kdof = 0; kdof < nmdof; ++kdof)
            {
              int col = Bdis.dof(dofsetb, mnode, kdof);

              // multiply the two shape functions
              double prod = sval_A(j) * mval_A(k) * jac * wgt;

              // dof to dof
              if (jdof == kdof)
              {
                if (abs(prod) > VOLMORTARINTTOL) M.assemble(prod, row, col);
              }
            }
          }

          // integrate D
          for (int k = 0; k < ns_; ++k)
          {
            Core::Nodes::Node* snode = sele.nodes()[k];
            int nddof = Adis.num_dof(dofseta, snode);

            for (int kdof = 0; kdof < nddof; ++kdof)
            {
              // multiply the two shape functions
              double prod = sval_A(j) * sval_A(k) * jac * wgt;

              // dof to dof
              if (jdof == kdof)
              {
                if (abs(prod) > VOLMORTARINTTOL) D.assemble(prod, row, row);
              }
            }
          }
        }
      }
    }
    else if (shape == shape_dual)
    {
      // loop over slave dofs
      for (int jdof = 0; jdof < nsdof; ++jdof)
      {
        const int row = Adis.dof(dofseta, cnode, jdof);

        if (not PAB_dofrowmap.my_gid(row)) continue;

        // integrate D
        const double prod2 = lmval_A(j) * sval_A(j) * jac * wgt;
        if (abs(prod2) > VOLMORTARINTTOL) D.assemble(prod2, row, row);

        // integrate M
        for (int k = 0; k < nm_; ++k)
        {
          Core::Nodes::Node* mnode = mele->nodes()[k];
          const int col = Bdis.dof(dofsetb, mnode, jdof);

          if (not PAB_dofcolmap.my_gid(col)) continue;

          // multiply the two shape functions
          const double prod = lmval_A(j) * mval_A(k) * jac * wgt;

          if (abs(prod) > VOLMORTARINTTOL) M.assemble(prod, row, col);
        }
      }
    }
    else
    {
      FOUR_C_THROW("ERROR: Unknown shape!");
    }
  }

  return true;
}


/*----------------------------------------------------------------------*
 |  ctor (public)                                            farah 01/14|
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype_s, Core::FE::CellType distype_m>
Coupling::VolMortar::VolMortarIntegrator<distype_s, distype_m>::VolMortarIntegrator(
    Teuchos::ParameterList& params)
{
  // get type of quadratic modification
  dualquad_ = Teuchos::getIntegralValue<DualQuad>(params, "DUALQUAD");

  // get type of quadratic modification
  shape_ = Teuchos::getIntegralValue<Shapefcn>(params, "SHAPEFCN");

  // define gp rule
  initialize_gp();
}


/*----------------------------------------------------------------------*
 |  Initialize gauss points                                  farah 01/14|
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype_s, Core::FE::CellType distype_m>
void Coupling::VolMortar::VolMortarIntegrator<distype_s, distype_m>::initialize_gp(
    bool integrateele, int domain, Core::FE::CellType shape)
{
  // init shape of integration domain
  Core::FE::CellType intshape = Core::FE::CellType::dis_none;

  if (integrateele)
  {
    if (domain == 0)
      intshape = distype_s;
    else if (domain == 1)
      intshape = distype_m;
    else
      FOUR_C_THROW("integration domain not specified!");
  }
  else
  {
    if (ndim_ == 2)
      intshape = Core::FE::CellType::tri3;
    else if (ndim_ == 3)
      intshape = shape;
    else
      FOUR_C_THROW("wrong dimension!");
  }

  //*******************************
  // choose Gauss rule accordingly
  //*******************************
  switch (intshape)
  {
    case Core::FE::CellType::tri3:
    {
      Core::FE::GaussRule2D mygaussrule = Core::FE::GaussRule2D::tri_7point;

      const Core::FE::IntegrationPoints2D intpoints(mygaussrule);
      ngp_ = intpoints.nquad;
      coords_.reshape(ngp_, 2);
      weights_.resize(ngp_);
      for (int i = 0; i < ngp_; ++i)
      {
        coords_(i, 0) = intpoints.qxg[i][0];
        coords_(i, 1) = intpoints.qxg[i][1];
        weights_[i] = intpoints.qwgt[i];
      }
      break;
    }
    case Core::FE::CellType::tet4:
    {
      Core::FE::GaussRule3D mygaussrule = Core::FE::GaussRule3D::tet_45point;

      const Core::FE::IntegrationPoints3D intpoints(mygaussrule);
      ngp_ = intpoints.nquad;
      coords_.reshape(ngp_, 3);
      weights_.resize(ngp_);
      for (int i = 0; i < ngp_; ++i)
      {
        coords_(i, 0) = intpoints.qxg[i][0];
        coords_(i, 1) = intpoints.qxg[i][1];
        coords_(i, 2) = intpoints.qxg[i][2];
        weights_[i] = intpoints.qwgt[i];
      }
      break;
    }
    case Core::FE::CellType::tet10:
    {
      Core::FE::GaussRule3D mygaussrule = Core::FE::GaussRule3D::tet_45point;

      const Core::FE::IntegrationPoints3D intpoints(mygaussrule);
      ngp_ = intpoints.nquad;
      coords_.reshape(ngp_, 3);
      weights_.resize(ngp_);
      for (int i = 0; i < ngp_; ++i)
      {
        coords_(i, 0) = intpoints.qxg[i][0];
        coords_(i, 1) = intpoints.qxg[i][1];
        coords_(i, 2) = intpoints.qxg[i][2];
        weights_[i] = intpoints.qwgt[i];
      }
      break;
    }
    case Core::FE::CellType::hex8:
    {
      Core::FE::GaussRule3D mygaussrule = Core::FE::GaussRule3D::hex_27point;

      const Core::FE::IntegrationPoints3D intpoints(mygaussrule);
      ngp_ = intpoints.nquad;
      coords_.reshape(ngp_, 3);
      weights_.resize(ngp_);
      for (int i = 0; i < ngp_; ++i)
      {
        coords_(i, 0) = intpoints.qxg[i][0];
        coords_(i, 1) = intpoints.qxg[i][1];
        coords_(i, 2) = intpoints.qxg[i][2];
        weights_[i] = intpoints.qwgt[i];
      }
      break;
    }
    case Core::FE::CellType::hex20:
    {
      Core::FE::GaussRule3D mygaussrule = Core::FE::GaussRule3D::hex_125point;

      const Core::FE::IntegrationPoints3D intpoints(mygaussrule);
      ngp_ = intpoints.nquad;
      coords_.reshape(ngp_, 3);
      weights_.resize(ngp_);
      for (int i = 0; i < ngp_; ++i)
      {
        coords_(i, 0) = intpoints.qxg[i][0];
        coords_(i, 1) = intpoints.qxg[i][1];
        coords_(i, 2) = intpoints.qxg[i][2];
        weights_[i] = intpoints.qwgt[i];
      }
      break;
    }
    case Core::FE::CellType::hex27:
    {
      Core::FE::GaussRule3D mygaussrule = Core::FE::GaussRule3D::hex_125point;

      const Core::FE::IntegrationPoints3D intpoints(mygaussrule);
      ngp_ = intpoints.nquad;
      coords_.reshape(ngp_, 3);
      weights_.resize(ngp_);
      for (int i = 0; i < ngp_; ++i)
      {
        coords_(i, 0) = intpoints.qxg[i][0];
        coords_(i, 1) = intpoints.qxg[i][1];
        coords_(i, 2) = intpoints.qxg[i][2];
        weights_[i] = intpoints.qwgt[i];
      }
      break;
    }
    case Core::FE::CellType::pyramid5:
    {
      Core::FE::GaussRule3D mygaussrule = Core::FE::GaussRule3D::pyramid_8point;

      const Core::FE::IntegrationPoints3D intpoints(mygaussrule);
      ngp_ = intpoints.nquad;
      coords_.reshape(ngp_, 3);
      weights_.resize(ngp_);
      for (int i = 0; i < ngp_; ++i)
      {
        coords_(i, 0) = intpoints.qxg[i][0];
        coords_(i, 1) = intpoints.qxg[i][1];
        coords_(i, 2) = intpoints.qxg[i][2];
        weights_[i] = intpoints.qwgt[i];
      }
      break;
    }
    default:
    {
      FOUR_C_THROW("ERROR: VolMortarIntegrator: This element type is not implemented!");
      break;
    }
  }  // switch(eletype)

  return;
}


/*----------------------------------------------------------------------*
 |  Compute D/M entries for Volumetric Mortar                farah 01/14|
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype_s, Core::FE::CellType distype_m>
void Coupling::VolMortar::VolMortarIntegrator<distype_s, distype_m>::integrate_cells_2d(
    Core::Elements::Element& sele, Core::Elements::Element& mele, Mortar::IntCell& cell,
    Core::LinAlg::SparseMatrix& dmatrix, Core::LinAlg::SparseMatrix& mmatrix,
    const Core::FE::Discretization& slavedis, const Core::FE::Discretization& masterdis,
    int sdofset, int mdofset)
{
  // create empty vectors for shape fct. evaluation
  Core::LinAlg::Matrix<ns_, 1> sval;
  Core::LinAlg::Matrix<nm_, 1> mval;
  Core::LinAlg::Matrix<ns_, 1> lmval;

  //**********************************************************************
  // loop over all Gauss points for integration
  //**********************************************************************
  for (int gp = 0; gp < ngp_; ++gp)
  {
    //    // coordinates and weight
    double eta[2] = {coords_(gp, 0), coords_(gp, 1)};
    double wgt = weights_[gp];

    // get global Gauss point coordinates
    double globgp[3] = {0.0, 0.0, 0.0};
    cell.local_to_global(eta, globgp, 0);

    // map gp into slave and master para space
    double sxi[3] = {0.0, 0.0, 0.0};
    double mxi[3] = {0.0, 0.0, 0.0};
    Mortar::Utils::global_to_local<distype_s>(sele, globgp, sxi);
    Mortar::Utils::global_to_local<distype_m>(mele, globgp, mxi);

    // Check parameter space mapping
    bool proj = check_mapping_2d(sele, mele, sxi, mxi);
    if (proj == false) FOUR_C_THROW("ERROR: Mapping failed!");

    // evaluate trace space shape functions (on both elements)
    Utils::shape_function<distype_s>(sval, sxi);
    Utils::shape_function<distype_m>(mval, mxi);

    // evaluate Lagrange multiplier shape functions (on slave element)
    Utils::dual_shape_function<distype_s>(lmval, sxi, sele);

    // evaluate the integration cell Jacobian
    double jac = cell.jacobian();

    // compute segment D/M matrix ****************************************
    // standard shape functions
    if (shape_ == shape_std)
    {
      for (int j = 0; j < ns_; ++j)
      {
        Core::Nodes::Node* cnode = sele.nodes()[j];
        int nsdof = slavedis.num_dof(sdofset, cnode);

        if (cnode->owner() != Core::Communication::my_mpi_rank(slavedis.get_comm())) continue;

        // loop over slave dofs
        for (int jdof = 0; jdof < nsdof; ++jdof)
        {
          int row = slavedis.dof(sdofset, cnode, jdof);

          ////////////////////////////////////////
          // integrate M
          for (int k = 0; k < nm_; ++k)
          {
            Core::Nodes::Node* mnode = mele.nodes()[k];
            int nmdof = masterdis.num_dof(mdofset, mnode);

            for (int kdof = 0; kdof < nmdof; ++kdof)
            {
              int col = masterdis.dof(mdofset, mnode, kdof);

              // multiply the two shape functions
              double prod = sval(j) * mval(k) * jac * wgt;

              // dof to dof
              if (jdof == kdof)
              {
                if (abs(prod) > VOLMORTARINTTOL) mmatrix.assemble(prod, row, col);
              }
            }
          }

          ////////////////////////////////////////
          // integrate D
          for (int k = 0; k < ns_; ++k)
          {
            Core::Nodes::Node* snode = sele.nodes()[k];
            int nddof = slavedis.num_dof(sdofset, snode);

            for (int kdof = 0; kdof < nddof; ++kdof)
            {
              // int col = slavedis->Dof(sdofset,snode,kdof);

              // multiply the two shape functions
              double prod = sval(j) * sval(k) * jac * wgt;

              // dof to dof
              if (jdof == kdof)
              {
                if (abs(prod) > VOLMORTARINTTOL) dmatrix.assemble(prod, row, row);
              }
            }
          }
        }
      }
    }
    else if (shape_ == shape_dual)
    {
      for (int j = 0; j < ns_; ++j)
      {
        Core::Nodes::Node* cnode = sele.nodes()[j];

        if (cnode->owner() != Core::Communication::my_mpi_rank(slavedis.get_comm())) continue;

        int nsdof = slavedis.num_dof(sdofset, cnode);

        // loop over slave dofs
        for (int jdof = 0; jdof < nsdof; ++jdof)
        {
          int row = slavedis.dof(sdofset, cnode, jdof);

          ////////////////////////////////////////////////////////////////
          // integrate M and D
          for (int k = 0; k < nm_; ++k)
          {
            Core::Nodes::Node* mnode = mele.nodes()[k];
            int nmdof = masterdis.num_dof(mdofset, mnode);

            for (int kdof = 0; kdof < nmdof; ++kdof)
            {
              int col = masterdis.dof(mdofset, mnode, kdof);

              // multiply the two shape functions
              double prod = lmval(j) * mval(k) * jac * wgt;

              // dof to dof
              if (jdof == kdof)
              {
                if (abs(prod) > VOLMORTARINTTOL) mmatrix.assemble(prod, row, col);
                if (abs(prod) > VOLMORTARINTTOL) dmatrix.assemble(prod, row, row);
              }
            }
          }
          ////////////////////////////////////////////////////////////////
        }
      }
    }
    else
    {
      FOUR_C_THROW("ERROR: Unknown shape function!");
    }
  }  // end gp loop

  return;
}


/*----------------------------------------------------------------------*
 |  Compute D/M entries for Volumetric Mortar                farah 01/14|
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype_s, Core::FE::CellType distype_m>
void Coupling::VolMortar::VolMortarIntegrator<distype_s, distype_m>::integrate_cells_3d(
    Core::Elements::Element& Aele, Core::Elements::Element& Bele, Coupling::VolMortar::Cell& cell,
    Core::LinAlg::SparseMatrix& dmatrix_A, Core::LinAlg::SparseMatrix& mmatrix_A,
    Core::LinAlg::SparseMatrix& dmatrix_B, Core::LinAlg::SparseMatrix& mmatrix_B,
    const Core::FE::Discretization& Adis, const Core::FE::Discretization& Bdis, int sdofset_A,
    int mdofset_A, int sdofset_B, int mdofset_B)
{
  if (shape_ == shape_std) FOUR_C_THROW("ERROR: std. shape functions not supported");

  // create empty vectors for shape fct. evaluation
  Core::LinAlg::Matrix<ns_, 1> sval_A;
  Core::LinAlg::Matrix<nm_, 1> mval_A;
  Core::LinAlg::Matrix<ns_, 1> lmval_A;
  Core::LinAlg::Matrix<nm_, 1> lmval_B;

  //**********************************************************************
  // loop over all Gauss points for integration
  //**********************************************************************
  for (int gp = 0; gp < ngp_; ++gp)
  {
    // coordinates and weight
    double eta[3] = {coords_(gp, 0), coords_(gp, 1), coords_(gp, 2)};
    double wgt = weights_[gp];

    // get global Gauss point coordinates
    double globgp[3] = {0.0, 0.0, 0.0};
    cell.local_to_global(eta, globgp);

    // map gp into A and B para space
    double Axi[3] = {0.0, 0.0, 0.0};
    double Bxi[3] = {0.0, 0.0, 0.0};
    Mortar::Utils::global_to_local<distype_s>(Aele, globgp, Axi);
    Mortar::Utils::global_to_local<distype_m>(Bele, globgp, Bxi);

    // evaluate the integration cell Jacobian
    double jac = 0.0;
    if (cell.shape() == Core::FE::CellType::tet4)
      jac = cell.vol();
    else if (cell.shape() == Core::FE::CellType::hex8)
      jac = cell.calc_jac(eta);
    else
      FOUR_C_THROW("used shape not supported in volmortar integrator!");

    // Check parameter space mapping
    // std::cout << "globgp " << globgp[0] <<"  "<< globgp[1] <<"  "<< globgp[2] <<std::endl;

    bool check = check_mapping_3d(Aele, Bele, Axi, Bxi);
    if (!check) continue;

    // evaluate trace space shape functions (on both elements)
    Utils::shape_function<distype_s>(sval_A, Axi);
    Utils::shape_function<distype_m>(mval_A, Bxi);

    // evaluate Lagrange multiplier shape functions (on slave element)
    Utils::dual_shape_function<distype_s>(lmval_A, Axi, Aele, dualquad_);
    Utils::dual_shape_function<distype_m>(lmval_B, Bxi, Bele, dualquad_);

    // compute cell D/M matrix ****************************************
    // dual shape functions
    for (int j = 0; j < ns_; ++j)
    {
      Core::Nodes::Node* cnode = Aele.nodes()[j];
      int nsdof = Adis.num_dof(sdofset_A, cnode);

      // loop over slave dofs
      for (int jdof = 0; jdof < nsdof; ++jdof)
      {
        int row = Adis.dof(sdofset_A, cnode, jdof);

        // integrate M and D
        for (int k = 0; k < nm_; ++k)
        {
          Core::Nodes::Node* mnode = Bele.nodes()[k];
          int nmdof = Bdis.num_dof(mdofset_A, mnode);

          for (int kdof = 0; kdof < nmdof; ++kdof)
          {
            int col = Bdis.dof(mdofset_A, mnode, kdof);

            // multiply the two shape functions
            double prod = lmval_A(j) * mval_A(k) * jac * wgt;

            // dof to dof
            if (jdof == kdof)
            {
              if (abs(prod) > VOLMORTARINTTOL) mmatrix_A.assemble(prod, row, col);
              if (abs(prod) > VOLMORTARINTTOL) dmatrix_A.assemble(prod, row, row);
            }
          }
        }
      }
    }

    // compute cell D/M matrix ****************************************
    // dual shape functions
    for (int j = 0; j < nm_; ++j)
    {
      Core::Nodes::Node* cnode = Bele.nodes()[j];
      int nsdof = Bdis.num_dof(sdofset_B, cnode);

      // loop over slave dofs
      for (int jdof = 0; jdof < nsdof; ++jdof)
      {
        int row = Bdis.dof(sdofset_B, cnode, jdof);

        // integrate M and D
        for (int k = 0; k < ns_; ++k)
        {
          Core::Nodes::Node* mnode = Aele.nodes()[k];
          int nmdof = Adis.num_dof(mdofset_B, mnode);

          for (int kdof = 0; kdof < nmdof; ++kdof)
          {
            int col = Adis.dof(mdofset_B, mnode, kdof);

            // multiply the two shape functions
            double prod = lmval_B(j) * sval_A(k) * jac * wgt;

            // dof to dof
            if (jdof == kdof)
            {
              if (abs(prod) > VOLMORTARINTTOL) mmatrix_B.assemble(prod, row, col);
              if (abs(prod) > VOLMORTARINTTOL) dmatrix_B.assemble(prod, row, row);
            }
          }
        }
      }
    }
  }  // end gp loop

  return;
}


/*----------------------------------------------------------------------*
 |  Compute D/M entries for Volumetric Mortar                farah 04/14|
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype_s, Core::FE::CellType distype_m>
void Coupling::VolMortar::VolMortarIntegrator<distype_s,
    distype_m>::integrate_cells_3d_direct_diveregence(Core::Elements::Element& Aele,
    Core::Elements::Element& Bele, Cut::VolumeCell& vc,
    std::shared_ptr<Core::FE::GaussPoints> intpoints, bool switched_conf,
    Core::LinAlg::SparseMatrix& dmatrix_A, Core::LinAlg::SparseMatrix& mmatrix_A,
    Core::LinAlg::SparseMatrix& dmatrix_B, Core::LinAlg::SparseMatrix& mmatrix_B,
    const Core::FE::Discretization& Adis, const Core::FE::Discretization& Bdis, int sdofset_A,
    int mdofset_A, int sdofset_B, int mdofset_B)
{
  if (shape_ == shape_std) FOUR_C_THROW("ERROR: std. shape functions not supported");

  // create empty vectors for shape fct. evaluation
  Core::LinAlg::Matrix<ns_, 1> sval_A;
  Core::LinAlg::Matrix<nm_, 1> mval_A;
  Core::LinAlg::Matrix<ns_, 1> lmval_A;
  Core::LinAlg::Matrix<nm_, 1> lmval_B;

  //**********************************************************************
  // loop over all Gauss points for integration
  //**********************************************************************
  for (int gp = 0; gp < intpoints->num_points(); ++gp)
  {
    double weight_out = intpoints->weight(gp);
    // coordinates and weight
    double eta[3] = {intpoints->point(gp)[0], intpoints->point(gp)[1], intpoints->point(gp)[2]};

    double globgp[3] = {0.0, 0.0, 0.0};

    if (switched_conf)
      Utils::local_to_global<distype_s>(Aele, eta, globgp);
    else
      Utils::local_to_global<distype_m>(Bele, eta, globgp);

    // map gp into A and B para space
    double Axi[3] = {0.0, 0.0, 0.0};
    double Bxi[3] = {0.0, 0.0, 0.0};
    Mortar::Utils::global_to_local<distype_s>(Aele, globgp, Axi);
    Mortar::Utils::global_to_local<distype_m>(Bele, globgp, Bxi);

    //      std::cout << "-------------------------------------" << std::endl;
    //      std::cout << "globgp= " << globgp[0] << "  " << globgp[1] << "  " << globgp[2] <<
    //      std::endl; std::cout << "eta= " << eta[0] << "  " << eta[1] << "  " << eta[2] <<
    //      std::endl; std::cout << "Axi= " << Axi[0] << "  " << Axi[1] << "  " << Axi[2] <<
    //      std::endl; std::cout << "Bxi= " << Bxi[0] << "  " << Bxi[1] << "  " << Bxi[2] <<
    //      std::endl;

    // evaluate the integration cell Jacobian
    double jac = 0.0;

    if (switched_conf)
      jac = Utils::jacobian<distype_s>(Axi, Aele);
    else
      jac = Utils::jacobian<distype_m>(Bxi, Bele);

    // Check parameter space mapping
    // check_mapping_3d(Aele,Bele,Axi,Bxi);

    // evaluate trace space shape functions (on both elements)
    Utils::shape_function<distype_s>(sval_A, Axi);
    Utils::shape_function<distype_m>(mval_A, Bxi);

    // evaluate Lagrange multiplier shape functions (on slave element)
    Utils::dual_shape_function<distype_s>(lmval_A, Axi, Aele, dualquad_);
    Utils::dual_shape_function<distype_m>(lmval_B, Bxi, Bele, dualquad_);

    // compute cell D/M matrix ****************************************
    // dual shape functions
    for (int j = 0; j < ns_; ++j)
    {
      Core::Nodes::Node* cnode = Aele.nodes()[j];
      int nsdof = Adis.num_dof(sdofset_A, cnode);

      // loop over slave dofs
      for (int jdof = 0; jdof < nsdof; ++jdof)
      {
        int row = Adis.dof(sdofset_A, cnode, jdof);

        // integrate M and D
        for (int k = 0; k < nm_; ++k)
        {
          Core::Nodes::Node* mnode = Bele.nodes()[k];
          int nmdof = Bdis.num_dof(mdofset_A, mnode);

          for (int kdof = 0; kdof < nmdof; ++kdof)
          {
            int col = Bdis.dof(mdofset_A, mnode, kdof);

            // multiply the two shape functions
            double prod = lmval_A(j) * mval_A(k) * jac * weight_out;
            //              std::cout << "PROD1 = " << prod  << " row= " << row << "  col= " << col
            //              << "  j= " << j<<  "  nsdof= " << nsdof<< std::endl; cnode->print(cout);
            // dof to dof
            if (jdof == kdof)
            {
              if (abs(prod) > VOLMORTARINTTOL) mmatrix_A.assemble(prod, row, col);
              if (abs(prod) > VOLMORTARINTTOL) dmatrix_A.assemble(prod, row, row);
            }
          }
        }
      }
    }

    // compute cell D/M matrix ****************************************
    // dual shape functions
    for (int j = 0; j < nm_; ++j)
    {
      Core::Nodes::Node* cnode = Bele.nodes()[j];
      int nsdof = Bdis.num_dof(sdofset_B, cnode);

      // loop over slave dofs
      for (int jdof = 0; jdof < nsdof; ++jdof)
      {
        int row = Bdis.dof(sdofset_B, cnode, jdof);

        // integrate M and D
        for (int k = 0; k < ns_; ++k)
        {
          Core::Nodes::Node* mnode = Aele.nodes()[k];
          int nmdof = Adis.num_dof(mdofset_B, mnode);

          for (int kdof = 0; kdof < nmdof; ++kdof)
          {
            int col = Adis.dof(sdofset_B, mnode, kdof);

            // multiply the two shape functions
            double prod = lmval_B(j) * sval_A(k) * jac * weight_out;
            //              std::cout << "PROD2 = " << prod  << " row= " << row << "  col= " <<
            //              col<< std::endl; cnode->print(cout);

            // dof to dof
            if (jdof == kdof)
            {
              if (abs(prod) > VOLMORTARINTTOL) mmatrix_B.assemble(prod, row, col);
              if (abs(prod) > VOLMORTARINTTOL) dmatrix_B.assemble(prod, row, row);
            }
          }
        }
      }
    }
  }  // end gp loop

  return;
}


/*----------------------------------------------------------------------*
 |  Compute D/M entries for Volumetric Mortar                farah 04/14|
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype_s, Core::FE::CellType distype_m>
void Coupling::VolMortar::VolMortarIntegrator<distype_s, distype_m>::integrate_ele_based_3d_a_dis(
    Core::Elements::Element& Aele, std::vector<int>& foundeles,
    Core::LinAlg::SparseMatrix& dmatrix_A, Core::LinAlg::SparseMatrix& mmatrix_A,
    const Core::FE::Discretization& Adis, const Core::FE::Discretization& Bdis, int dofsetA,
    int dofsetB)
{
  if (shape_ == shape_std) FOUR_C_THROW("ERROR: std. shape functions not supported");

  // create empty vectors for shape fct. evaluation
  Core::LinAlg::Matrix<ns_, 1> sval_A;
  Core::LinAlg::Matrix<nm_, 1> mval_A;
  Core::LinAlg::Matrix<ns_, 1> lmval_A;

  //**********************************************************************
  // loop over all Gauss points for integration
  //**********************************************************************
  for (int gp = 0; gp < ngp_; ++gp)
  {
    // coordinates and weight
    double eta[3] = {coords_(gp, 0), coords_(gp, 1), coords_(gp, 2)};
    double wgt = weights_[gp];
    double jac = 0.0;
    double globgp[3] = {0.0, 0.0, 0.0};

    // quantities for eval. outside gp
    double gpdist = 1.0e12;
    int gpid = 0;
    std::array<double, 3> AuxXi = {0.0, 0.0, 0.0};

    // evaluate the integration cell Jacobian
    jac = Utils::jacobian<distype_s>(eta, Aele);

    // get global Gauss point coordinates
    Utils::local_to_global<distype_s>(Aele, eta, globgp);

    // map gp into A and B para space
    double Axi[3] = {0.0, 0.0, 0.0};
    Mortar::Utils::global_to_local<distype_s>(Aele, globgp, Axi);

    // loop over beles
    for (int found = 0; found < (int)foundeles.size(); ++found)
    {
      // get master element
      Core::Elements::Element* Bele = Bdis.g_element(foundeles[found]);
      double Bxi[3] = {0.0, 0.0, 0.0};

      bool converged = true;
      Mortar::Utils::global_to_local<distype_m>(*Bele, globgp, Bxi, converged);
      if (!converged and found != ((int)foundeles.size() - 1)) continue;

      // save distance of gp
      double l = sqrt(Bxi[0] * Bxi[0] + Bxi[1] * Bxi[1] + Bxi[2] * Bxi[2]);
      if (l < gpdist)
      {
        gpdist = l;
        gpid = foundeles[found];
        AuxXi[0] = Bxi[0];
        AuxXi[1] = Bxi[1];
        AuxXi[2] = Bxi[2];
      }

      // Check parameter space mapping
      bool proj = check_mapping_3d(Aele, *Bele, Axi, Bxi);

      // if gp outside continue or eval nearest gp
      if (!proj and (found != ((int)foundeles.size() - 1)))
        continue;
      else if (!proj and found == ((int)foundeles.size() - 1))
      {
        Bxi[0] = AuxXi[0];
        Bxi[1] = AuxXi[1];
        Bxi[2] = AuxXi[2];
        Bele = Bdis.g_element(gpid);
      }

      // for "master" side
      Utils::shape_function<distype_s>(sval_A, Axi, dualquad_);
      Utils::shape_function<distype_m>(mval_A, Bxi);

      // evaluate Lagrange multiplier shape functions (on slave element)
      Utils::dual_shape_function<distype_s>(lmval_A, Axi, Aele, dualquad_);

      // compute cell D/M matrix ****************************************
      // dual shape functions
      for (int j = 0; j < ns_; ++j)
      {
        Core::Nodes::Node* cnode = Aele.nodes()[j];
        if (cnode->owner() != Core::Communication::my_mpi_rank(Adis.get_comm())) continue;

        int nsdof = Adis.num_dof(dofsetA, cnode);

        // loop over slave dofs
        for (int jdof = 0; jdof < nsdof; ++jdof)
        {
          int row = Adis.dof(dofsetA, cnode, jdof);

          // integrate D
          double prod2 = lmval_A(j) * sval_A(j) * jac * wgt;
          if (abs(prod2) > VOLMORTARINTTOL) dmatrix_A.assemble(prod2, row, row);

          // integrate M
          for (int k = 0; k < nm_; ++k)
          {
            Core::Nodes::Node* mnode = Bele->nodes()[k];
            int col = Bdis.dof(dofsetB, mnode, jdof);

            // multiply the two shape functions
            double prod = lmval_A(j) * mval_A(k) * jac * wgt;

            if (abs(prod) > VOLMORTARINTTOL) mmatrix_A.assemble(prod, row, col);
          }
        }
      }

      break;
    }  // beles
  }  // end gp loop

  return;
}


/*----------------------------------------------------------------------*
 |  Compute D/M entries for Volumetric Mortar                farah 04/14|
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype_s, Core::FE::CellType distype_m>
void Coupling::VolMortar::VolMortarIntegrator<distype_s, distype_m>::integrate_ele_based_3d_b_dis(
    Core::Elements::Element& Bele, std::vector<int>& foundeles,
    Core::LinAlg::SparseMatrix& dmatrix_B, Core::LinAlg::SparseMatrix& mmatrix_B,
    const Core::FE::Discretization& Adis, const Core::FE::Discretization& Bdis, int dofsetA,
    int dofsetB)
{
  if (shape_ == shape_std) FOUR_C_THROW("ERROR: std. shape functions not supported");

  // create empty vectors for shape fct. evaluation
  Core::LinAlg::Matrix<ns_, 1> mval_A;
  Core::LinAlg::Matrix<nm_, 1> sval_B;
  Core::LinAlg::Matrix<nm_, 1> lmval_B;

  //**********************************************************************
  // loop over all Gauss points for integration
  //**********************************************************************
  for (int gp = 0; gp < ngp_; ++gp)
  {
    //    // coordinates and weight
    double eta[3] = {coords_(gp, 0), coords_(gp, 1), coords_(gp, 2)};
    double wgt = weights_[gp];
    double jac = 0.0;
    double globgp[3] = {0.0, 0.0, 0.0};

    // quantities for eval. outside gp
    double gpdist = 1.0e12;
    int gpid = 0;
    double AuxXi[3] = {0.0, 0.0, 0.0};

    // evaluate the integration cell Jacobian
    jac = Utils::jacobian<distype_m>(eta, Bele);

    // get global Gauss point coordinates
    Utils::local_to_global<distype_m>(Bele, eta, globgp);

    // map gp into A and B para space
    double Bxi[3] = {0.0, 0.0, 0.0};
    Mortar::Utils::global_to_local<distype_m>(Bele, globgp, Bxi);

    // loop over beles
    for (int found = 0; found < (int)foundeles.size(); ++found)
    {
      // get master element
      Core::Elements::Element* Aele = Adis.g_element(foundeles[found]);
      double Axi[3] = {0.0, 0.0, 0.0};

      bool converged = true;
      Mortar::Utils::global_to_local<distype_s>(*Aele, globgp, Axi, converged);
      if (!converged and found != ((int)foundeles.size() - 1)) continue;

      // save distance of gp
      double l = sqrt(Axi[0] * Axi[0] + Axi[1] * Axi[1] + Axi[2] * Axi[2]);
      if (l < gpdist)
      {
        gpdist = l;
        gpid = foundeles[found];
        AuxXi[0] = Axi[0];
        AuxXi[1] = Axi[1];
        AuxXi[2] = Axi[2];
      }

      // Check parameter space mapping
      bool proj = check_mapping_3d(*Aele, Bele, Axi, Bxi);

      // if gp outside continue or eval nearest gp
      if (!proj and (found != ((int)foundeles.size() - 1)))
        continue;
      else if (!proj and found == ((int)foundeles.size() - 1))
      {
        Axi[0] = AuxXi[0];
        Axi[1] = AuxXi[1];
        Axi[2] = AuxXi[2];
        Aele = Adis.g_element(gpid);
      }

      // evaluate trace space shape functions (on both elements)
      Utils::shape_function<distype_m>(sval_B, Bxi, dualquad_);
      Utils::shape_function<distype_s>(mval_A, Axi);

      // evaluate Lagrange multiplier shape functions (on slave element)
      Utils::dual_shape_function<distype_m>(lmval_B, Bxi, Bele, dualquad_);
      // compute cell D/M matrix ****************************************
      // dual shape functions
      for (int j = 0; j < nm_; ++j)
      {
        Core::Nodes::Node* cnode = Bele.nodes()[j];
        if (cnode->owner() != Core::Communication::my_mpi_rank(Bdis.get_comm())) continue;

        int nsdof = Bdis.num_dof(dofsetB, cnode);

        // loop over slave dofs
        for (int jdof = 0; jdof < nsdof; ++jdof)
        {
          int row = Bdis.dof(dofsetB, cnode, jdof);

          // integrate D
          double prod2 = lmval_B(j) * sval_B(j) * jac * wgt;
          if (abs(prod2) > VOLMORTARINTTOL) dmatrix_B.assemble(prod2, row, row);

          // integrate M
          for (int k = 0; k < ns_; ++k)
          {
            Core::Nodes::Node* mnode = Aele->nodes()[k];
            int col = Adis.dof(dofsetA, mnode, jdof);

            // multiply the two shape functions
            double prod = lmval_B(j) * mval_A(k) * jac * wgt;

            if (abs(prod) > VOLMORTARINTTOL) mmatrix_B.assemble(prod, row, col);
          }
        }
      }

      break;
    }  // beles
  }  // end gp loop

  return;
}


/*----------------------------------------------------------------------*
 |  Compute D/M entries for Volumetric Mortar                farah 01/14|
 |  This function is for element-wise integration when an               |
 |  element is completely located within an other element               |
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype_s, Core::FE::CellType distype_m>
void Coupling::VolMortar::VolMortarIntegrator<distype_s, distype_m>::integrate_ele_3d(int domain,
    Core::Elements::Element& Aele, Core::Elements::Element& Bele,
    Core::LinAlg::SparseMatrix& dmatrix_A, Core::LinAlg::SparseMatrix& mmatrix_A,
    Core::LinAlg::SparseMatrix& dmatrix_B, Core::LinAlg::SparseMatrix& mmatrix_B,
    const Core::FE::Discretization& Adis, const Core::FE::Discretization& Bdis, int sdofset_A,
    int mdofset_A, int sdofset_B, int mdofset_B)
{
  if (shape_ == shape_std) FOUR_C_THROW("ERROR: std. shape functions not supported");

  // create empty vectors for shape fct. evaluation
  Core::LinAlg::Matrix<ns_, 1> sval_A;
  Core::LinAlg::Matrix<nm_, 1> mval_A;
  Core::LinAlg::Matrix<ns_, 1> lmval_A;
  Core::LinAlg::Matrix<nm_, 1> lmval_B;

  //**********************************************************************
  // loop over all Gauss points for integration
  //**********************************************************************
  for (int gp = 0; gp < ngp_; ++gp)
  {
    //    // coordinates and weight
    double eta[3] = {coords_(gp, 0), coords_(gp, 1), coords_(gp, 2)};
    double wgt = weights_[gp];
    double jac = 0.0;
    double globgp[3] = {0.0, 0.0, 0.0};


    if (domain == 0)
    {
      // evaluate the integration cell Jacobian
      jac = Utils::jacobian<distype_s>(eta, Aele);

      // get global Gauss point coordinates
      Utils::local_to_global<distype_s>(Aele, eta, globgp);
    }
    else if (domain == 1)
    {
      // evaluate the integration cell Jacobian
      jac = Utils::jacobian<distype_m>(eta, Bele);

      // get global Gauss point coordinates
      Utils::local_to_global<distype_m>(Bele, eta, globgp);
    }
    else
      FOUR_C_THROW("wrong domain for integration!");


    // map gp into A and B para space
    double Axi[3] = {0.0, 0.0, 0.0};
    double Bxi[3] = {0.0, 0.0, 0.0};
    Mortar::Utils::global_to_local<distype_s>(Aele, globgp, Axi);
    Mortar::Utils::global_to_local<distype_m>(Bele, globgp, Bxi);

    // Check parameter space mapping
    check_mapping_3d(Aele, Bele, Axi, Bxi);

    // evaluate trace space shape functions (on both elements)
    Utils::shape_function<distype_s>(sval_A, Axi);
    Utils::shape_function<distype_m>(mval_A, Bxi);

    // evaluate Lagrange multiplier shape functions (on slave element)
    Utils::dual_shape_function<distype_s>(lmval_A, Axi, Aele, dualquad_);
    Utils::dual_shape_function<distype_m>(lmval_B, Bxi, Bele, dualquad_);

    // compute cell D/M matrix ****************************************
    // dual shape functions
    for (int j = 0; j < ns_; ++j)
    {
      Core::Nodes::Node* cnode = Aele.nodes()[j];
      int nsdof = Adis.num_dof(sdofset_A, cnode);

      // loop over slave dofs
      for (int jdof = 0; jdof < nsdof; ++jdof)
      {
        int row = Adis.dof(sdofset_A, cnode, jdof);

        // integrate M and D
        for (int k = 0; k < nm_; ++k)
        {
          Core::Nodes::Node* mnode = Bele.nodes()[k];
          int nmdof = Bdis.num_dof(mdofset_A, mnode);

          for (int kdof = 0; kdof < nmdof; ++kdof)
          {
            int col = Bdis.dof(mdofset_A, mnode, kdof);

            // multiply the two shape functions
            double prod = lmval_A(j) * mval_A(k) * jac * wgt;

            // dof to dof
            if (jdof == kdof)
            {
              if (abs(prod) > VOLMORTARINTTOL) mmatrix_A.assemble(prod, row, col);
              if (abs(prod) > VOLMORTARINTTOL) dmatrix_A.assemble(prod, row, row);
            }
          }
        }
      }
    }

    // compute cell D/M matrix ****************************************
    // dual shape functions
    for (int j = 0; j < nm_; ++j)
    {
      Core::Nodes::Node* cnode = Bele.nodes()[j];
      int nsdof = Bdis.num_dof(sdofset_B, cnode);

      // loop over slave dofs
      for (int jdof = 0; jdof < nsdof; ++jdof)
      {
        int row = Bdis.dof(sdofset_B, cnode, jdof);

        // integrate M and D
        for (int k = 0; k < ns_; ++k)
        {
          Core::Nodes::Node* mnode = Aele.nodes()[k];
          int nmdof = Adis.num_dof(mdofset_B, mnode);

          for (int kdof = 0; kdof < nmdof; ++kdof)
          {
            int col = Adis.dof(mdofset_B, mnode, kdof);

            // multiply the two shape functions
            double prod = lmval_B(j) * sval_A(k) * jac * wgt;

            // dof to dof
            if (jdof == kdof)
            {
              if (abs(prod) > VOLMORTARINTTOL) mmatrix_B.assemble(prod, row, col);
              if (abs(prod) > VOLMORTARINTTOL) dmatrix_B.assemble(prod, row, row);
            }
          }
        }
      }
    }

  }  // end gp loop

  return;
}


/*----------------------------------------------------------------------*
 |  Compute D/M entries for Volumetric Mortar                farah 01/14|
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype_s, Core::FE::CellType distype_m>
bool Coupling::VolMortar::VolMortarIntegrator<distype_s, distype_m>::check_mapping_2d(
    Core::Elements::Element& sele, Core::Elements::Element& mele, double* sxi, double* mxi)
{
  // check GP projection (SLAVE)
  const double tol = 1e-10;
  if (distype_s == Core::FE::CellType::quad4 || distype_s == Core::FE::CellType::quad8 ||
      distype_s == Core::FE::CellType::quad9)
  {
    if (sxi[0] < -1.0 - tol || sxi[1] < -1.0 - tol || sxi[0] > 1.0 + tol || sxi[1] > 1.0 + tol)
    {
      std::cout << "\n***Warning: Gauss point projection outside!";
      std::cout << "Slave ID: " << sele.id() << " Master ID: " << mele.id() << std::endl;
      std::cout << "Slave GP projection: " << sxi[0] << " " << sxi[1] << std::endl;
      return false;
    }
  }
  else if (distype_s == Core::FE::CellType::tri3 || distype_s == Core::FE::CellType::tri6)
  {
    if (sxi[0] < -tol || sxi[1] < -tol || sxi[0] > 1.0 + tol || sxi[1] > 1.0 + tol ||
        sxi[0] + sxi[1] > 1.0 + 2 * tol)
    {
      std::cout << "\n***Warning: Gauss point projection outside!";
      std::cout << "Slave ID: " << sele.id() << " Master ID: " << mele.id() << std::endl;
      std::cout << "Slave GP projection: " << sxi[0] << " " << sxi[1] << std::endl;
      return false;
    }
  }
  else
    FOUR_C_THROW("Wrong element type!");

  // check GP projection (MASTER)
  if (distype_m == Core::FE::CellType::quad4 || distype_m == Core::FE::CellType::quad8 ||
      distype_m == Core::FE::CellType::quad9)
  {
    if (mxi[0] < -1.0 - tol || mxi[1] < -1.0 - tol || mxi[0] > 1.0 + tol || mxi[1] > 1.0 + tol)
    {
      std::cout << "\n***Warning: Gauss point projection outside!";
      std::cout << "Slave ID: " << sele.id() << " Master ID: " << mele.id() << std::endl;
      std::cout << "Master GP projection: " << mxi[0] << " " << mxi[1] << std::endl;
      return false;
    }
  }
  else if (distype_s == Core::FE::CellType::tri3 || distype_s == Core::FE::CellType::tri6)
  {
    if (mxi[0] < -tol || mxi[1] < -tol || mxi[0] > 1.0 + tol || mxi[1] > 1.0 + tol ||
        mxi[0] + mxi[1] > 1.0 + 2 * tol)
    {
      std::cout << "\n***Warning: Gauss point projection outside!";
      std::cout << "Slave ID: " << sele.id() << " Master ID: " << mele.id() << std::endl;
      std::cout << "Master GP projection: " << mxi[0] << " " << mxi[1] << std::endl;
      return false;
    }
  }
  else
    FOUR_C_THROW("Wrong element type!");

  return true;
}


/*----------------------------------------------------------------------*
 |  Compute D/M entries for Volumetric Mortar                farah 01/14|
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype_s, Core::FE::CellType distype_m>
bool Coupling::VolMortar::VolMortarIntegrator<distype_s, distype_m>::check_mapping_3d(
    Core::Elements::Element& sele, Core::Elements::Element& mele, double* sxi, double* mxi)
{
  // check GP projection (SLAVE)
  double tol = 1e-5;
  if (distype_s == Core::FE::CellType::hex8 || distype_s == Core::FE::CellType::hex20 ||
      distype_s == Core::FE::CellType::hex27)
  {
    if (sxi[0] < -1.0 - tol || sxi[1] < -1.0 - tol || sxi[2] < -1.0 - tol || sxi[0] > 1.0 + tol ||
        sxi[1] > 1.0 + tol || sxi[2] > 1.0 + tol)
    {
      //      std::cout << "\n***Warning: Gauss point projection outside!";
      //      std::cout << "Slave ID: " << sele.Id() << " Master ID: " << mele.Id() << std::endl;
      //      std::cout << "Slave GP projection: " << sxi[0] << " " << sxi[1] << " " << sxi[2] <<
      //      std::endl;
      //
      //      for(int i=0;i<sele.num_node();++i)
      //      {
      //        std::cout << "create vertex " << sele.Nodes()[i]->X()[0] <<"  "<<
      //        sele.Nodes()[i]->X()[1] <<"  "<< sele.Nodes()[i]->X()[2] <<std::endl;
      //      }
      //      std::cout << "------------" << std::endl;
      //      for(int i=0;i<mele.num_node();++i)
      //      {
      //        std::cout << "create vertex " << mele.Nodes()[i]->X()[0] <<"  "<<
      //        mele.Nodes()[i]->X()[1] <<"  "<< mele.Nodes()[i]->X()[2] <<std::endl;
      //      }

      return false;
    }
  }
  else if (distype_s == Core::FE::CellType::tet4 || distype_s == Core::FE::CellType::tet10)
  {
    if (sxi[0] < 0.0 - tol || sxi[1] < 0.0 - tol || sxi[2] < 0.0 - tol ||
        (sxi[0] + sxi[1] + sxi[2]) > 1.0 + tol)
    {
      //      std::cout << "\n***Warning: Gauss point projection outside!";
      //      std::cout << "Slave ID: " << sele.Id() << " Master ID: " << mele.Id() << std::endl;
      //      std::cout << "Slave GP projection: " << sxi[0] << " " << sxi[1] << " " << sxi[2] <<
      //      std::endl; for(int i=0;i<sele.num_node();++i)
      //      {
      //        std::cout << "create vertex " << sele.Nodes()[i]->X()[0] <<"  "<<
      //        sele.Nodes()[i]->X()[1] <<"  "<< sele.Nodes()[i]->X()[2] <<std::endl;
      //      }
      //      std::cout << "------------" << std::endl;
      //      for(int i=0;i<mele.num_node();++i)
      //      {
      //        std::cout << "create vertex " << mele.Nodes()[i]->X()[0] <<"  "<<
      //        mele.Nodes()[i]->X()[1] <<"  "<< mele.Nodes()[i]->X()[2] <<std::endl;
      //      }
      return false;
    }
  }
  else if (distype_s == Core::FE::CellType::pyramid5)
  {
    if (sxi[2] < 0.0 - tol || -sxi[0] + sxi[2] > 1.0 + tol || sxi[0] + sxi[2] > 1.0 + tol ||
        -sxi[1] + sxi[2] > 1.0 + tol || sxi[1] + sxi[2] > 1.0 + tol)
    {
      //      std::cout << "\n***Warning: Gauss point projection outside!";
      //      std::cout << "Slave ID: " << sele.Id() << " Master ID: " << mele.Id() << std::endl;
      //      std::cout << "Slave GP projection: " << sxi[0] << " " << sxi[1] << " " << sxi[2] <<
      //      std::endl; for(int i=0;i<sele.num_node();++i)
      //      {
      //        std::cout << "create vertex " << sele.Nodes()[i]->X()[0] <<"  "<<
      //        sele.Nodes()[i]->X()[1] <<"  "<< sele.Nodes()[i]->X()[2] <<std::endl;
      //      }
      //      std::cout << "------------" << std::endl;
      //      for(int i=0;i<mele.num_node();++i)
      //      {
      //        std::cout << "create vertex " << mele.Nodes()[i]->X()[0] <<"  "<<
      //        mele.Nodes()[i]->X()[1] <<"  "<< mele.Nodes()[i]->X()[2] <<std::endl;
      //      }
      return false;
    }
  }
  else
    FOUR_C_THROW("Wrong element type!");

  // check GP projection (MASTER)
  if (distype_m == Core::FE::CellType::hex8 || distype_m == Core::FE::CellType::hex20 ||
      distype_m == Core::FE::CellType::hex27)
  {
    if (mxi[0] < -1.0 - tol || mxi[1] < -1.0 - tol || mxi[2] < -1.0 - tol || mxi[0] > 1.0 + tol ||
        mxi[1] > 1.0 + tol || mxi[2] > 1.0 + tol)
    {
      //      std::cout << "\n***Warning: Gauss point projection outside!";
      //      std::cout << "Slave ID: " << sele.Id() << " Master ID: " << mele.Id() << std::endl;
      //      std::cout << "Master GP projection: " << mxi[0] << " " << mxi[1] << " " << mxi[2] <<
      //      std::endl; for(int i=0;i<sele.num_node();++i)
      //      {
      //        std::cout << "create vertex " << sele.Nodes()[i]->X()[0] <<"  "<<
      //        sele.Nodes()[i]->X()[1] <<"  "<< sele.Nodes()[i]->X()[2] <<std::endl;
      //      }
      //      std::cout << "------------" << std::endl;
      //      for(int i=0;i<mele.num_node();++i)
      //      {
      //        std::cout << "create vertex " << mele.Nodes()[i]->X()[0] <<"  "<<
      //        mele.Nodes()[i]->X()[1] <<"  "<< mele.Nodes()[i]->X()[2] <<std::endl;
      //      }
      return false;
    }
  }
  else if (distype_m == Core::FE::CellType::tet4 || distype_m == Core::FE::CellType::tet10)
  {
    if (mxi[0] < 0.0 - tol || mxi[1] < 0.0 - tol || mxi[2] < 0.0 - tol ||
        (mxi[0] + mxi[1] + mxi[2]) > 1.0 + tol)
    {
      //      std::cout << "\n***Warning: Gauss point projection outside!";
      //      std::cout << "Slave ID: " << sele.Id() << " Master ID: " << mele.Id() << std::endl;
      //      std::cout << "Master GP projection: " << mxi[0] << " " << mxi[1] << " " << mxi[2] <<
      //      std::endl; for(int i=0;i<sele.num_node();++i)
      //      {
      //        std::cout << "create vertex " << sele.Nodes()[i]->X()[0] <<"  "<<
      //        sele.Nodes()[i]->X()[1] <<"  "<< sele.Nodes()[i]->X()[2] <<std::endl;
      //      }
      //      std::cout << "------------" << std::endl;
      //      for(int i=0;i<mele.num_node();++i)
      //      {
      //        std::cout << "create vertex " << mele.Nodes()[i]->X()[0] <<"  "<<
      //        mele.Nodes()[i]->X()[1] <<"  "<< mele.Nodes()[i]->X()[2] <<std::endl;
      //      }
      return false;
    }
  }
  else if (distype_m == Core::FE::CellType::pyramid5)
  {
    if (mxi[2] < 0.0 - tol || -mxi[0] + mxi[2] > 1.0 + tol || mxi[0] + mxi[2] > 1.0 + tol ||
        -mxi[1] + mxi[2] > 1.0 + tol || mxi[1] + mxi[2] > 1.0 + tol)
    {
      //      std::cout << "\n***Warning: Gauss point projection outside!";
      //      std::cout << "Slave ID: " << sele.Id() << " Master ID: " << mele.Id() << std::endl;
      //      std::cout << "Master GP projection: " << mxi[0] << " " << mxi[1] << " " << mxi[2] <<
      //      std::endl; for(int i=0;i<sele.num_node();++i)
      //      {
      //        std::cout << "create vertex " << sele.Nodes()[i]->X()[0] <<"  "<<
      //        sele.Nodes()[i]->X()[1] <<"  "<< sele.Nodes()[i]->X()[2] <<std::endl;
      //      }
      //      std::cout << "------------" << std::endl;
      //      for(int i=0;i<mele.num_node();++i)
      //      {
      //        std::cout << "create vertex " << mele.Nodes()[i]->X()[0] <<"  "<<
      //        mele.Nodes()[i]->X()[1] <<"  "<< mele.Nodes()[i]->X()[2] <<std::endl;
      //      }
      return false;
    }
  }
  else
    FOUR_C_THROW("Wrong element type!");

  return true;
}


/*----------------------------------------------------------------------*
 |  possible slave/master element pairs                       farah 01/14|
 *----------------------------------------------------------------------*/
// slave quad4
template class Coupling::VolMortar::VolMortarIntegrator<Core::FE::CellType::quad4,
    Core::FE::CellType::quad4>;
template class Coupling::VolMortar::VolMortarIntegrator<Core::FE::CellType::quad4,
    Core::FE::CellType::tri3>;

// slave tri3
template class Coupling::VolMortar::VolMortarIntegrator<Core::FE::CellType::tri3,
    Core::FE::CellType::quad4>;
template class Coupling::VolMortar::VolMortarIntegrator<Core::FE::CellType::tri3,
    Core::FE::CellType::tri3>;

// slave hex8
template class Coupling::VolMortar::VolMortarIntegrator<Core::FE::CellType::hex8,
    Core::FE::CellType::tet4>;
template class Coupling::VolMortar::VolMortarIntegrator<Core::FE::CellType::hex8,
    Core::FE::CellType::tet10>;
template class Coupling::VolMortar::VolMortarIntegrator<Core::FE::CellType::hex8,
    Core::FE::CellType::hex8>;
template class Coupling::VolMortar::VolMortarIntegrator<Core::FE::CellType::hex8,
    Core::FE::CellType::hex27>;
template class Coupling::VolMortar::VolMortarIntegrator<Core::FE::CellType::hex8,
    Core::FE::CellType::hex20>;
template class Coupling::VolMortar::VolMortarIntegrator<Core::FE::CellType::hex8,
    Core::FE::CellType::pyramid5>;

// slave hex20
template class Coupling::VolMortar::VolMortarIntegrator<Core::FE::CellType::hex20,
    Core::FE::CellType::tet4>;
template class Coupling::VolMortar::VolMortarIntegrator<Core::FE::CellType::hex20,
    Core::FE::CellType::tet10>;
template class Coupling::VolMortar::VolMortarIntegrator<Core::FE::CellType::hex20,
    Core::FE::CellType::hex8>;
template class Coupling::VolMortar::VolMortarIntegrator<Core::FE::CellType::hex20,
    Core::FE::CellType::hex27>;
template class Coupling::VolMortar::VolMortarIntegrator<Core::FE::CellType::hex20,
    Core::FE::CellType::hex20>;
template class Coupling::VolMortar::VolMortarIntegrator<Core::FE::CellType::hex20,
    Core::FE::CellType::pyramid5>;

// slave hex27
template class Coupling::VolMortar::VolMortarIntegrator<Core::FE::CellType::hex27,
    Core::FE::CellType::tet4>;
template class Coupling::VolMortar::VolMortarIntegrator<Core::FE::CellType::hex27,
    Core::FE::CellType::tet10>;
template class Coupling::VolMortar::VolMortarIntegrator<Core::FE::CellType::hex27,
    Core::FE::CellType::hex8>;
template class Coupling::VolMortar::VolMortarIntegrator<Core::FE::CellType::hex27,
    Core::FE::CellType::hex27>;
template class Coupling::VolMortar::VolMortarIntegrator<Core::FE::CellType::hex27,
    Core::FE::CellType::hex20>;
template class Coupling::VolMortar::VolMortarIntegrator<Core::FE::CellType::hex27,
    Core::FE::CellType::pyramid5>;

// slave tet4
template class Coupling::VolMortar::VolMortarIntegrator<Core::FE::CellType::tet4,
    Core::FE::CellType::tet4>;
template class Coupling::VolMortar::VolMortarIntegrator<Core::FE::CellType::tet4,
    Core::FE::CellType::tet10>;
template class Coupling::VolMortar::VolMortarIntegrator<Core::FE::CellType::tet4,
    Core::FE::CellType::hex8>;
template class Coupling::VolMortar::VolMortarIntegrator<Core::FE::CellType::tet4,
    Core::FE::CellType::hex27>;
template class Coupling::VolMortar::VolMortarIntegrator<Core::FE::CellType::tet4,
    Core::FE::CellType::hex20>;
template class Coupling::VolMortar::VolMortarIntegrator<Core::FE::CellType::tet4,
    Core::FE::CellType::pyramid5>;

// slave tet10
template class Coupling::VolMortar::VolMortarIntegrator<Core::FE::CellType::tet10,
    Core::FE::CellType::tet4>;
template class Coupling::VolMortar::VolMortarIntegrator<Core::FE::CellType::tet10,
    Core::FE::CellType::tet10>;
template class Coupling::VolMortar::VolMortarIntegrator<Core::FE::CellType::tet10,
    Core::FE::CellType::hex8>;
template class Coupling::VolMortar::VolMortarIntegrator<Core::FE::CellType::tet10,
    Core::FE::CellType::hex27>;
template class Coupling::VolMortar::VolMortarIntegrator<Core::FE::CellType::tet10,
    Core::FE::CellType::hex20>;
template class Coupling::VolMortar::VolMortarIntegrator<Core::FE::CellType::tet10,
    Core::FE::CellType::pyramid5>;

// slave pyramid 5
template class Coupling::VolMortar::VolMortarIntegrator<Core::FE::CellType::pyramid5,
    Core::FE::CellType::tet4>;
template class Coupling::VolMortar::VolMortarIntegrator<Core::FE::CellType::pyramid5,
    Core::FE::CellType::tet10>;
template class Coupling::VolMortar::VolMortarIntegrator<Core::FE::CellType::pyramid5,
    Core::FE::CellType::hex8>;
template class Coupling::VolMortar::VolMortarIntegrator<Core::FE::CellType::pyramid5,
    Core::FE::CellType::hex27>;
template class Coupling::VolMortar::VolMortarIntegrator<Core::FE::CellType::pyramid5,
    Core::FE::CellType::hex20>;
template class Coupling::VolMortar::VolMortarIntegrator<Core::FE::CellType::pyramid5,
    Core::FE::CellType::pyramid5>;

/*----------------------------------------------------------------------*
 |  ctor (public)                                            farah 06/14|
 *----------------------------------------------------------------------*/
Coupling::VolMortar::ConsInterpolator::ConsInterpolator()
{
  // empty
}


/*----------------------------------------------------------------------*
 |  interpolate (public)                                     farah 06/14|
 *----------------------------------------------------------------------*/
void Coupling::VolMortar::ConsInterpolator::interpolate(Core::Nodes::Node* node,
    Core::LinAlg::SparseMatrix& pmatrix, const Core::FE::Discretization& nodediscret,
    const Core::FE::Discretization& elediscret, std::vector<int>& foundeles,
    std::pair<int, int>& dofset, const Core::LinAlg::Map& P_dofrowmap,
    const Core::LinAlg::Map& P_dofcolmap)
{
  // check ownership
  if (node->owner() != Core::Communication::my_mpi_rank(nodediscret.get_comm())) return;

  // map gp into A and B para space
  double nodepos[3];

  const auto node_x = node->x();
  nodepos[0] = node_x[0];
  nodepos[1] = node_x[1];
  if (node_x.size() > 2)
    nodepos[2] = node_x[2];
  else
    nodepos[2] = 0.0;

  double dist = 1.0e12;
  int eleid = 0;
  double AuxXi[3] = {0.0, 0.0, 0.0};

  // element loop (brute force)
  for (int found = 0; found < (int)foundeles.size(); ++found)
  {
    bool proj = false;

    // get master element
    Core::Elements::Element* ele = elediscret.g_element(foundeles[found]);

    switch (ele->shape())
    {
      // 2D --------------------------------------------
      case Core::FE::CellType::tri3:
      {
        proj = cons_interpolator_eval<Core::FE::CellType::tri3>(node, ele, pmatrix, nodediscret,
            elediscret, foundeles, found, eleid, dist, AuxXi, nodepos, dofset, P_dofrowmap,
            P_dofcolmap);
        break;
      }
      case Core::FE::CellType::tri6:
      {
        proj = cons_interpolator_eval<Core::FE::CellType::tri6>(node, ele, pmatrix, nodediscret,
            elediscret, foundeles, found, eleid, dist, AuxXi, nodepos, dofset, P_dofrowmap,
            P_dofcolmap);
        break;
      }
      case Core::FE::CellType::quad4:
      {
        proj = cons_interpolator_eval<Core::FE::CellType::quad4>(node, ele, pmatrix, nodediscret,
            elediscret, foundeles, found, eleid, dist, AuxXi, nodepos, dofset, P_dofrowmap,
            P_dofcolmap);
        break;
      }
      case Core::FE::CellType::quad8:
      {
        proj = cons_interpolator_eval<Core::FE::CellType::quad8>(node, ele, pmatrix, nodediscret,
            elediscret, foundeles, found, eleid, dist, AuxXi, nodepos, dofset, P_dofrowmap,
            P_dofcolmap);
        break;
      }
      case Core::FE::CellType::quad9:
      {
        proj = cons_interpolator_eval<Core::FE::CellType::quad9>(node, ele, pmatrix, nodediscret,
            elediscret, foundeles, found, eleid, dist, AuxXi, nodepos, dofset, P_dofrowmap,
            P_dofcolmap);
        break;
      }

      // 3D --------------------------------------------
      case Core::FE::CellType::hex8:
      {
        proj = cons_interpolator_eval<Core::FE::CellType::hex8>(node, ele, pmatrix, nodediscret,
            elediscret, foundeles, found, eleid, dist, AuxXi, nodepos, dofset, P_dofrowmap,
            P_dofcolmap);
        break;
      }
      case Core::FE::CellType::hex20:
      {
        proj = cons_interpolator_eval<Core::FE::CellType::hex20>(node, ele, pmatrix, nodediscret,
            elediscret, foundeles, found, eleid, dist, AuxXi, nodepos, dofset, P_dofrowmap,
            P_dofcolmap);
        break;
      }
      case Core::FE::CellType::hex27:
      {
        proj = cons_interpolator_eval<Core::FE::CellType::hex27>(node, ele, pmatrix, nodediscret,
            elediscret, foundeles, found, eleid, dist, AuxXi, nodepos, dofset, P_dofrowmap,
            P_dofcolmap);
        break;
      }
      case Core::FE::CellType::tet4:
      {
        proj = cons_interpolator_eval<Core::FE::CellType::tet4>(node, ele, pmatrix, nodediscret,
            elediscret, foundeles, found, eleid, dist, AuxXi, nodepos, dofset, P_dofrowmap,
            P_dofcolmap);
        break;
      }
      case Core::FE::CellType::tet10:
      {
        proj = cons_interpolator_eval<Core::FE::CellType::tet10>(node, ele, pmatrix, nodediscret,
            elediscret, foundeles, found, eleid, dist, AuxXi, nodepos, dofset, P_dofrowmap,
            P_dofcolmap);
        break;
      }
      case Core::FE::CellType::pyramid5:
      {
        proj = cons_interpolator_eval<Core::FE::CellType::pyramid5>(node, ele, pmatrix, nodediscret,
            elediscret, foundeles, found, eleid, dist, AuxXi, nodepos, dofset, P_dofrowmap,
            P_dofcolmap);
        break;
      }
      default:
      {
        FOUR_C_THROW("ERROR: unknown shape!");
        break;
      }
    }  // end switch

    // if node evaluated break ele loop
    if (proj == true)
      break;
    else
      continue;

    break;
  }

  return;
}

/*----------------------------------------------------------------------*
 |  node evaluation                                          farah 02/15|
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
bool Coupling::VolMortar::cons_interpolator_eval(Core::Nodes::Node* node,
    Core::Elements::Element* ele, Core::LinAlg::SparseMatrix& pmatrix,
    const Core::FE::Discretization& nodediscret, const Core::FE::Discretization& elediscret,
    std::vector<int>& foundeles, int& found, int& eleid, double& dist, double* AuxXi,
    double* nodepos, std::pair<int, int>& dofset, const Core::LinAlg::Map& P_dofrowmap,
    const Core::LinAlg::Map& P_dofcolmap)
{
  //! ns_: number of slave element nodes
  static const int n_ = Core::FE::num_nodes(distype);

  double xi[3] = {0.0, 0.0, 0.0};
  bool converged = true;

  Mortar::Utils::global_to_local<distype>(*ele, nodepos, xi, converged);

  // no convergence of local newton?
  if (!converged and found != ((int)foundeles.size() - 1)) return false;

  // save distance of gp
  double l = sqrt(xi[0] * xi[0] + xi[1] * xi[1] + xi[2] * xi[2]);
  if (l < dist)
  {
    dist = l;
    eleid = foundeles[found];
    AuxXi[0] = xi[0];
    AuxXi[1] = xi[1];
    AuxXi[2] = xi[2];
  }

  // Check parameter space mapping
  bool proj = check_mapping<distype>(*ele, xi);

  // if node outside --> continue or eval nearest gp
  if (!proj and found != ((int)foundeles.size() - 1))
  {
    return false;
  }
  else if (!proj and found == ((int)foundeles.size() - 1))
  {
    xi[0] = AuxXi[0];
    xi[1] = AuxXi[1];
    xi[2] = AuxXi[2];
    ele = elediscret.g_element(eleid);
  }

  // get values
  Core::LinAlg::Matrix<n_, 1> val;
  Utils::shape_function<distype>(val, xi);

  int nsdof = nodediscret.num_dof(dofset.first, node);

  // loop over slave dofs
  for (int jdof = 0; jdof < nsdof; ++jdof)
  {
    const int row = nodediscret.dof(dofset.first, node, jdof);

    if (not P_dofrowmap.my_gid(row)) continue;

    for (int k = 0; k < ele->num_node(); ++k)
    {
      Core::Nodes::Node* bnode = ele->nodes()[k];
      const int col = elediscret.dof(dofset.second, bnode, jdof);

      if (not P_dofcolmap.my_gid(col)) continue;

      const double val2 = val(k);
      // if (abs(val2)>VOLMORTARINTTOL)
      pmatrix.assemble(val2, row, col);
    }
  }

  return true;
}


/*----------------------------------------------------------------------*
 |  possible elements for interpolation                      farah 06/14|
 *----------------------------------------------------------------------*/
// template class Coupling::VolMortar::ConsInterpolator<Core::FE::CellType::quad4>;
// template class Coupling::VolMortar::ConsInterpolator<Core::FE::CellType::quad8>;
// template class Coupling::VolMortar::ConsInterpolator<Core::FE::CellType::quad9>;
//
// template class Coupling::VolMortar::ConsInterpolator<Core::FE::CellType::tri3>;
// template class Coupling::VolMortar::ConsInterpolator<Core::FE::CellType::tri6>;
//
// template class Coupling::VolMortar::ConsInterpolator<Core::FE::CellType::hex8>;
// template class Coupling::VolMortar::ConsInterpolator<Core::FE::CellType::hex20>;
// template class Coupling::VolMortar::ConsInterpolator<Core::FE::CellType::hex27>;
//
// template class Coupling::VolMortar::ConsInterpolator<Core::FE::CellType::tet4>;
// template class Coupling::VolMortar::ConsInterpolator<Core::FE::CellType::tet10>;

FOUR_C_NAMESPACE_CLOSE
