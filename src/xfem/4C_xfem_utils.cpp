// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_xfem_utils.hpp"

#include "4C_fem_discretization_faces.hpp"
#include "4C_fem_general_element.hpp"
#include "4C_fem_general_extract_values.hpp"
#include "4C_linalg_utils_sparse_algebra_math.hpp"
#include "4C_mat_list.hpp"
#include "4C_mat_material_factory.hpp"
#include "4C_mat_newtonianfluid.hpp"
#include "4C_material_base.hpp"
#include "4C_rebalance_binning_based.hpp"

FOUR_C_NAMESPACE_OPEN

void XFEM::Utils::extract_node_vectors(Core::FE::Discretization& dis,
    std::map<int, Core::LinAlg::Matrix<3, 1>>& nodevecmap,
    std::shared_ptr<Core::LinAlg::Vector<double>> idispnp)
{
  std::shared_ptr<const Core::LinAlg::Vector<double>> dispcol =
      Core::Rebalance::get_col_version_of_row_vector(dis, idispnp);
  nodevecmap.clear();

  for (int lid = 0; lid < dis.num_my_col_nodes(); ++lid)
  {
    const Core::Nodes::Node* node = dis.l_col_node(lid);
    std::vector<int> lm;
    dis.dof(node, lm);
    std::vector<double> mydisp = Core::FE::extract_values(*dispcol, lm);
    if (mydisp.size() < 3) FOUR_C_THROW("we need at least 3 dofs here");

    Core::LinAlg::Matrix<3, 1> currpos;
    currpos(0) = node->x()[0] + mydisp[0];
    currpos(1) = node->x()[1] + mydisp[1];
    currpos(2) = node->x()[2] + mydisp[2];
    nodevecmap.insert(std::make_pair(node->id(), currpos));
  }
}

// -------------------------------------------------------------------
// set master and slave parameters (winter 01/2015)
// -------------------------------------------------------------------
void XFEM::Utils::get_volume_cell_material(Core::Elements::Element* actele,
    std::shared_ptr<Core::Mat::Material>& mat, Cut::Point::PointPosition position)
{
  int position_id = 0;
  if (position == Cut::Point::inside)  // minus domain, Omega^i with i<j
    position_id = 1;
  else if (position != Cut::Point::outside)  // plus domain, \Omega^j with j>i
    FOUR_C_THROW("Volume cell is either undecided or on surface. That can't be good....");

  std::shared_ptr<Core::Mat::Material> material = actele->material();

  if (material->material_type() == Core::Materials::m_matlist)
  {
    // get material list for this element
    const Mat::MatList* matlist = static_cast<const Mat::MatList*>(material.get());
    int numofmaterials = matlist->num_mat();

    // Error messages
    if (numofmaterials > 2)
    {
      FOUR_C_THROW("More than two materials is currently not supported.");
    }

    // set default id in list of materials
    int matid = -1;
    matid = matlist->mat_id(position_id);
    mat = matlist->material_by_id(matid);
  }
  else
  {
    mat = material;
  }

  return;
}

/*----------------------------------------------------------------------*
 | Checks if Materials in parent and neighbor element are identical     |
 |                                                         winter 01/15 |
 *----------------------------------------------------------------------*/
void XFEM::Utils::safety_check_materials(
    std::shared_ptr<Core::Mat::Material>& pmat, std::shared_ptr<Core::Mat::Material>& nmat)
{
  //------------------------------ see whether materials in patch are equal

  if (pmat->material_type() != nmat->material_type())
    FOUR_C_THROW(" not the same material for master and slave parent element");

  if (pmat->material_type() == Core::Materials::m_matlist)
    FOUR_C_THROW(
        "A matlist has been found in edge based stabilization! If you are running XTPF, check "
        "calls as this should NOT happen!!!");

  if (pmat->material_type() != Core::Materials::m_carreauyasuda &&
      pmat->material_type() != Core::Materials::m_modpowerlaw &&
      pmat->material_type() != Core::Materials::m_herschelbulkley &&
      pmat->material_type() != Core::Materials::m_fluid)
    FOUR_C_THROW("Material law for parent element is not a fluid");

  if (pmat->material_type() == Core::Materials::m_fluid)
  {
    {
      const Mat::NewtonianFluid* actmat_p = static_cast<const Mat::NewtonianFluid*>(pmat.get());
      const double pvisc = actmat_p->viscosity();
      const double pdens = actmat_p->density();

      const Mat::NewtonianFluid* actmat_m = static_cast<const Mat::NewtonianFluid*>(nmat.get());
      const double nvisc = actmat_m->viscosity();
      const double ndens = actmat_m->density();

      if (std::abs(nvisc - pvisc) > 1e-14)
      {
        std::cout << "Parent element viscosity: " << pvisc
                  << " ,neighbor element viscosity: " << nvisc << std::endl;
        FOUR_C_THROW("parent and neighbor element do not have the same viscosity!");
      }
      if (std::abs(ndens - pdens) > 1e-14)
      {
        std::cout << "Parent element density: " << pdens << " ,neighbor element density: " << ndens
                  << std::endl;
        FOUR_C_THROW("parent and neighbor element do not have the same density!");
      }
    }
  }
  else
  {
    FOUR_C_THROW("up to now I expect a FLUID (m_fluid) material for edge stabilization\n");
  }

  return;
}

//! Extract a quantity for an element
void XFEM::Utils::extract_quantity_at_element(Core::LinAlg::SerialDenseMatrix& element_vector,
    const Core::Elements::Element* element,
    const Core::LinAlg::MultiVector<double>& global_col_vector, Core::FE::Discretization& dis,
    const int nds_vector, const int nsd)
{
  // get the other nds-set which is connected to the current one via this boundary-cell
  Core::Elements::LocationArray la(dis.num_dof_sets());
  element->location_vector(dis, la);

  const size_t numnode = element->num_node();

  if (la[nds_vector].lm_.size() != numnode)
  {
    std::cout << "la[nds_vector].lm_.size(): " << la[nds_vector].lm_.size() << std::endl;
    FOUR_C_THROW("assume a unique level-set dof in cutterdis-Dofset per node");
  }

  std::vector<double> local_vector =
      Core::FE::extract_values(global_col_vector, la[nds_vector].lm_);

  if (local_vector.size() != nsd * numnode)
    FOUR_C_THROW("wrong size of (potentially resized) local matrix!");

  // copy local to normal....
  Core::LinAlg::copy(local_vector.data(), element_vector);
}


//! Extract a quantity for a node
void XFEM::Utils::extract_quantity_at_node(Core::LinAlg::SerialDenseMatrix& element_vector,
    const Core::Nodes::Node* node, const Core::LinAlg::MultiVector<double>& global_col_vector,
    Core::FE::Discretization& dis, const int nds_vector, const unsigned int nsd)
{
  const std::vector<int> lm = dis.dof(nds_vector, node);
  if (lm.size() != 1) FOUR_C_THROW("assume a unique level-set dof in cutterdis-Dofset");

  std::vector<double> local_vector = Core::FE::extract_values(global_col_vector, lm);

  if (local_vector.size() != nsd) FOUR_C_THROW("wrong size of (potentially resized) local matrix!");

  // copy local to nvec....
  Core::LinAlg::copy(local_vector.data(), element_vector);
}

FOUR_C_NAMESPACE_CLOSE
