// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_fluid_ele_calc_xfem_coupling.hpp"

#include "4C_bele_bele3.hpp"
#include "4C_fluid_ele_calc_xfem_coupling_impl.hpp"
#include "4C_utils_enum.hpp"

FOUR_C_NAMESPACE_OPEN

using namespace Discret::Elements::XFLUID;

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
void SlaveElementInterface<distype>::define_state_names(
    Core::FE::CellType slave_distype,  ///< coupling slave discretization type
    std::string& disp_statename,       ///< name of displacement state at current step
    std::string& vel_statename,        ///< name of velocity state at current step
    std::string& veln_statename        ///< name of velocity state at previous step
)
{
  switch (slave_distype)
  {
      //  case Core::FE::CellType::tri3:
      //  case Core::FE::CellType::tri6:
    case Core::FE::CellType::quad4:
    case Core::FE::CellType::quad8:
    case Core::FE::CellType::quad9:
    {
      disp_statename = std::string("idispnp");
      vel_statename = std::string("ivelnp");
      veln_statename = std::string("iveln");
      break;
    }
    case Core::FE::CellType::hex8:
    case Core::FE::CellType::hex20:
    case Core::FE::CellType::hex27:
    {
      disp_statename = std::string("dispnp");
      vel_statename = std::string("velaf");
      veln_statename = std::string("veln");
      break;
    }
    default:
      FOUR_C_THROW("Unsupported element shape {}", slave_distype);
      break;
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
std::shared_ptr<SlaveElementInterface<distype>>
SlaveElementInterface<distype>::create_slave_element_representation(
    Core::Elements::Element* slave_ele,  ///< coupling slave element
    Core::LinAlg::SerialDenseMatrix&
        slave_xyz  ///< global node coordinates of coupling slave element
)
{
  SlaveElementInterface* sla = nullptr;

  // get number of dofs for this slave element
  const unsigned numdofpernode = slave_ele->num_dof_per_node(*slave_ele->nodes()[0]);

  if (numdofpernode == 3)
  {
    switch (slave_ele->shape())
    {
      case Core::FE::CellType::quad4:
      {
        using SlaveEleType = SlaveElementRepresentation<distype, Core::FE::CellType::quad4, 3>;
        sla = new SlaveEleType(slave_xyz);
        break;
      }
      case Core::FE::CellType::quad8:
      {
        using SlaveEleType = SlaveElementRepresentation<distype, Core::FE::CellType::quad8, 3>;
        sla = new SlaveEleType(slave_xyz);
        break;
      }
      case Core::FE::CellType::quad9:
      {
        using SlaveEleType = SlaveElementRepresentation<distype, Core::FE::CellType::quad9, 3>;
        sla = new SlaveEleType(slave_xyz);
        break;
      }
      default:
        FOUR_C_THROW("Unsupported boundary element shape {}", slave_ele->shape());
        break;
    }
  }
  else if (numdofpernode ==
           4)  // volumetric coupling partners only required for fluid-fluid coupling
  {
    switch (slave_ele->shape())
    {
      case Core::FE::CellType::quad4:
      {
        using SlaveEleType = SlaveElementRepresentation<distype, Core::FE::CellType::quad4, 4>;
        sla = new SlaveEleType(slave_xyz);
        break;
      }
      case Core::FE::CellType::quad8:
      {
        using SlaveEleType = SlaveElementRepresentation<distype, Core::FE::CellType::quad8, 4>;
        sla = new SlaveEleType(slave_xyz);
        break;
      }
      case Core::FE::CellType::quad9:
      {
        using SlaveEleType = SlaveElementRepresentation<distype, Core::FE::CellType::quad9, 4>;
        sla = new SlaveEleType(slave_xyz);
        break;
      }
      case Core::FE::CellType::hex8:
      {
        using SlaveEleType = SlaveElementRepresentation<distype, Core::FE::CellType::hex8, 4>;
        sla = new SlaveEleType(slave_xyz);
        break;
      }
      case Core::FE::CellType::hex20:
      {
        using SlaveEleType = SlaveElementRepresentation<distype, Core::FE::CellType::hex20, 4>;
        sla = new SlaveEleType(slave_xyz);
        break;
      }
      case Core::FE::CellType::hex27:
      {
        using SlaveEleType = SlaveElementRepresentation<distype, Core::FE::CellType::hex27, 4>;
        sla = new SlaveEleType(slave_xyz);
        break;
      }
      default:
        FOUR_C_THROW("unsupported boundary element shape {}", slave_ele->shape());
        break;
    }
  }

  return std::shared_ptr<SlaveElementInterface<distype>>(sla);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
std::shared_ptr<NitscheInterface<distype>>
NitscheInterface<distype>::create_nitsche_coupling_x_fluid_wdbc(
    Core::LinAlg::SerialDenseMatrix& C_umum, Core::LinAlg::SerialDenseMatrix& rhC_um,
    const Discret::Elements::FluidEleParameterXFEM& fldparaxfem)
{
  NitscheInterface* nit = nullptr;
  using NitscheCouplType = NitscheCoupling<distype, Core::FE::CellType::dis_none, 3>;
  nit = new NitscheCouplType(C_umum, rhC_um, fldparaxfem);

  return std::shared_ptr<NitscheInterface<distype>>(nit);
}

template <Core::FE::CellType distype>
std::shared_ptr<NitscheInterface<distype>>
NitscheInterface<distype>::create_nitsche_coupling_x_fluid_wdbc(
    Core::LinAlg::SerialDenseMatrix& C_umum, Core::LinAlg::SerialDenseVector& rhC_um,
    const Discret::Elements::FluidEleParameterXFEM& fldparaxfem)
{
  NitscheInterface* nit = nullptr;
  using NitscheCouplType = NitscheCoupling<distype, Core::FE::CellType::dis_none, 3>;
  nit = new NitscheCouplType(C_umum, rhC_um, fldparaxfem);

  return std::shared_ptr<NitscheInterface<distype>>(nit);
}
/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
std::shared_ptr<NitscheInterface<distype>>
NitscheInterface<distype>::create_nitsche_coupling_x_fluid_wdbc(Core::Elements::Element* bele,
    Core::LinAlg::SerialDenseMatrix& bele_xyz, Core::LinAlg::SerialDenseMatrix& C_umum,
    Core::LinAlg::SerialDenseMatrix& rhC_um,
    const Discret::Elements::FluidEleParameterXFEM& fldparaxfem)
{
  NitscheInterface* nit = nullptr;

  // get number of dofs for this boundary element
  const unsigned numdofpernode = bele->num_dof_per_node(*bele->nodes()[0]);

  // three dofs per node, for standard Dirichlet coupling, four dofs per node for background
  // geometry coupling
  if (numdofpernode == 3)
  {
    switch (bele->shape())
    {
      case Core::FE::CellType::quad4:
      {
        using NitscheCouplType = NitscheCoupling<distype, Core::FE::CellType::quad4, 3>;
        nit = new NitscheCouplType(bele_xyz, C_umum, rhC_um, fldparaxfem);
        break;
      }
      case Core::FE::CellType::quad8:
      {
        using NitscheCouplType = NitscheCoupling<distype, Core::FE::CellType::quad8, 3>;
        nit = new NitscheCouplType(bele_xyz, C_umum, rhC_um, fldparaxfem);
        break;
      }
      case Core::FE::CellType::quad9:
      {
        using NitscheCouplType = NitscheCoupling<distype, Core::FE::CellType::quad9, 3>;
        nit = new NitscheCouplType(bele_xyz, C_umum, rhC_um, fldparaxfem);
        break;
      }
      default:
        FOUR_C_THROW("Unsupported boundary element shape {}", bele->shape());
        break;
    }
  }
  else if (numdofpernode == 4)
  {
    switch (bele->shape())
    {
      case Core::FE::CellType::quad4:
      {
        using NitscheCouplType = NitscheCoupling<distype, Core::FE::CellType::quad4, 4>;
        nit = new NitscheCouplType(bele_xyz, C_umum, rhC_um, fldparaxfem);
        break;
      }
      case Core::FE::CellType::quad8:
      {
        using NitscheCouplType = NitscheCoupling<distype, Core::FE::CellType::quad8, 4>;
        nit = new NitscheCouplType(bele_xyz, C_umum, rhC_um, fldparaxfem);
        break;
      }
      case Core::FE::CellType::quad9:
      {
        using NitscheCouplType = NitscheCoupling<distype, Core::FE::CellType::quad9, 4>;
        nit = new NitscheCouplType(bele_xyz, C_umum, rhC_um, fldparaxfem);
        break;
      }
      default:
        FOUR_C_THROW("Unsupported boundary element shape {}", bele->shape());
        break;
    }
  }
  else
    FOUR_C_THROW("Unsupported number of {} nodes for coupling slave element.", numdofpernode);

  return std::shared_ptr<NitscheInterface<distype>>(nit);
}

template <Core::FE::CellType distype>
std::shared_ptr<NitscheInterface<distype>>
NitscheInterface<distype>::create_nitsche_coupling_x_fluid_wdbc(Core::Elements::Element* bele,
    Core::LinAlg::SerialDenseMatrix& bele_xyz, Core::LinAlg::SerialDenseMatrix& C_umum,
    Core::LinAlg::SerialDenseVector& rhC_um,
    const Discret::Elements::FluidEleParameterXFEM& fldparaxfem)
{
  NitscheInterface* nit = nullptr;

  // get number of dofs for this boundary element
  const unsigned numdofpernode = bele->num_dof_per_node(*bele->nodes()[0]);

  // three dofs per node, for standard Dirichlet coupling, four dofs per node for background
  // geometry coupling
  if (numdofpernode == 3)
  {
    switch (bele->shape())
    {
      case Core::FE::CellType::quad4:
      {
        using NitscheCouplType = NitscheCoupling<distype, Core::FE::CellType::quad4, 3>;
        nit = new NitscheCouplType(bele_xyz, C_umum, rhC_um, fldparaxfem);
        break;
      }
      case Core::FE::CellType::quad8:
      {
        using NitscheCouplType = NitscheCoupling<distype, Core::FE::CellType::quad8, 3>;
        nit = new NitscheCouplType(bele_xyz, C_umum, rhC_um, fldparaxfem);
        break;
      }
      case Core::FE::CellType::quad9:
      {
        using NitscheCouplType = NitscheCoupling<distype, Core::FE::CellType::quad9, 3>;
        nit = new NitscheCouplType(bele_xyz, C_umum, rhC_um, fldparaxfem);
        break;
      }
      default:
        FOUR_C_THROW("Unsupported boundary element shape {}", bele->shape());
        break;
    }
  }
  else if (numdofpernode == 4)
  {
    switch (bele->shape())
    {
      case Core::FE::CellType::quad4:
      {
        using NitscheCouplType = NitscheCoupling<distype, Core::FE::CellType::quad4, 4>;
        nit = new NitscheCouplType(bele_xyz, C_umum, rhC_um, fldparaxfem);
        break;
      }
      case Core::FE::CellType::quad8:
      {
        using NitscheCouplType = NitscheCoupling<distype, Core::FE::CellType::quad8, 4>;
        nit = new NitscheCouplType(bele_xyz, C_umum, rhC_um, fldparaxfem);
        break;
      }
      case Core::FE::CellType::quad9:
      {
        using NitscheCouplType = NitscheCoupling<distype, Core::FE::CellType::quad9, 4>;
        nit = new NitscheCouplType(bele_xyz, C_umum, rhC_um, fldparaxfem);
        break;
      }
      default:
        FOUR_C_THROW("Unsupported boundary element shape {}", bele->shape());
        break;
    }
  }
  else
    FOUR_C_THROW("Unsupported number of {} nodes for coupling slave element.", numdofpernode);

  return std::shared_ptr<NitscheInterface<distype>>(nit);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
std::shared_ptr<NitscheInterface<distype>>
NitscheInterface<distype>::create_nitsche_coupling_x_fluid_sided(Core::Elements::Element* bele,
    Core::LinAlg::SerialDenseMatrix& bele_xyz, Core::LinAlg::SerialDenseMatrix& C_umum,
    Core::LinAlg::SerialDenseMatrix& C_usum, Core::LinAlg::SerialDenseMatrix& C_umus,
    Core::LinAlg::SerialDenseMatrix& C_usus, Core::LinAlg::SerialDenseMatrix& rhC_um,
    Core::LinAlg::SerialDenseMatrix& rhC_us,
    const Discret::Elements::FluidEleParameterXFEM& fldparaxfem)
{
  NitscheInterface* nit = nullptr;

  // get number of dofs for this boundary element
  const unsigned numdofpernode = bele->num_dof_per_node(*bele->nodes()[0]);

  // three dofs per node, for monolithic XFSI, four dofs per node for background-sided fluid-fluid
  // coupling
  if (numdofpernode == 3)
  {
    switch (bele->shape())
    {
      case Core::FE::CellType::quad4:
      {
        using NitscheCouplType = NitscheCoupling<distype, Core::FE::CellType::quad4, 3>;
        nit = new NitscheCouplType(
            bele_xyz, C_umum, C_usum, C_umus, C_usus, rhC_um, rhC_us, fldparaxfem);
        break;
      }
      case Core::FE::CellType::quad8:
      {
        using NitscheCouplType = NitscheCoupling<distype, Core::FE::CellType::quad8, 3>;
        nit = new NitscheCouplType(
            bele_xyz, C_umum, C_usum, C_umus, C_usus, rhC_um, rhC_us, fldparaxfem);
        break;
      }
      case Core::FE::CellType::quad9:
      {
        using NitscheCouplType = NitscheCoupling<distype, Core::FE::CellType::quad9, 3>;
        nit = new NitscheCouplType(
            bele_xyz, C_umum, C_usum, C_umus, C_usus, rhC_um, rhC_us, fldparaxfem);
        break;
      }
      default:
        FOUR_C_THROW("Unsupported boundary element shape {}", bele->shape());
        break;
    }
  }
  else if (numdofpernode == 4)
  {
    switch (bele->shape())
    {
      case Core::FE::CellType::quad4:
      {
        using NitscheCouplType = NitscheCoupling<distype, Core::FE::CellType::quad4, 4>;
        nit = new NitscheCouplType(
            bele_xyz, C_umum, C_usum, C_umus, C_usus, rhC_um, rhC_us, fldparaxfem);
        break;
      }
      case Core::FE::CellType::quad8:
      {
        using NitscheCouplType = NitscheCoupling<distype, Core::FE::CellType::quad8, 4>;
        nit = new NitscheCouplType(
            bele_xyz, C_umum, C_usum, C_umus, C_usus, rhC_um, rhC_us, fldparaxfem);
        break;
      }
      case Core::FE::CellType::quad9:
      {
        using NitscheCouplType = NitscheCoupling<distype, Core::FE::CellType::quad9, 4>;
        nit = new NitscheCouplType(
            bele_xyz, C_umum, C_usum, C_umus, C_usus, rhC_um, rhC_us, fldparaxfem);
        break;
      }
      default:
        FOUR_C_THROW("Unsupported boundary element shape {}", bele->shape());
        break;
    }
  }
  else
    FOUR_C_THROW("Unsupported number of {} nodes for coupling slave element.", numdofpernode);

  return std::shared_ptr<NitscheInterface<distype>>(nit);
}

template <Core::FE::CellType distype>
std::shared_ptr<NitscheInterface<distype>>
NitscheInterface<distype>::create_nitsche_coupling_x_fluid_sided(Core::Elements::Element* bele,
    Core::LinAlg::SerialDenseMatrix& bele_xyz, Core::LinAlg::SerialDenseMatrix& C_umum,
    Core::LinAlg::SerialDenseMatrix& C_usum, Core::LinAlg::SerialDenseMatrix& C_umus,
    Core::LinAlg::SerialDenseMatrix& C_usus, Core::LinAlg::SerialDenseVector& rhC_um,
    Core::LinAlg::SerialDenseMatrix& rhC_us,
    const Discret::Elements::FluidEleParameterXFEM& fldparaxfem)
{
  NitscheInterface* nit = nullptr;

  // get number of dofs for this boundary element
  const unsigned numdofpernode = bele->num_dof_per_node(*bele->nodes()[0]);

  // three dofs per node, for monolithic XFSI, four dofs per node for background-sided fluid-fluid
  // coupling
  if (numdofpernode == 3)
  {
    switch (bele->shape())
    {
      case Core::FE::CellType::quad4:
      {
        using NitscheCouplType = NitscheCoupling<distype, Core::FE::CellType::quad4, 3>;
        nit = new NitscheCouplType(
            bele_xyz, C_umum, C_usum, C_umus, C_usus, rhC_um, rhC_us, fldparaxfem);
        break;
      }
      case Core::FE::CellType::quad8:
      {
        using NitscheCouplType = NitscheCoupling<distype, Core::FE::CellType::quad8, 3>;
        nit = new NitscheCouplType(
            bele_xyz, C_umum, C_usum, C_umus, C_usus, rhC_um, rhC_us, fldparaxfem);
        break;
      }
      case Core::FE::CellType::quad9:
      {
        using NitscheCouplType = NitscheCoupling<distype, Core::FE::CellType::quad9, 3>;
        nit = new NitscheCouplType(
            bele_xyz, C_umum, C_usum, C_umus, C_usus, rhC_um, rhC_us, fldparaxfem);
        break;
      }
      default:
        FOUR_C_THROW("Unsupported boundary element shape {}", bele->shape());
        break;
    }
  }
  else if (numdofpernode == 4)
  {
    switch (bele->shape())
    {
      case Core::FE::CellType::quad4:
      {
        using NitscheCouplType = NitscheCoupling<distype, Core::FE::CellType::quad4, 4>;
        nit = new NitscheCouplType(
            bele_xyz, C_umum, C_usum, C_umus, C_usus, rhC_um, rhC_us, fldparaxfem);
        break;
      }
      case Core::FE::CellType::quad8:
      {
        using NitscheCouplType = NitscheCoupling<distype, Core::FE::CellType::quad8, 4>;
        nit = new NitscheCouplType(
            bele_xyz, C_umum, C_usum, C_umus, C_usus, rhC_um, rhC_us, fldparaxfem);
        break;
      }
      case Core::FE::CellType::quad9:
      {
        using NitscheCouplType = NitscheCoupling<distype, Core::FE::CellType::quad9, 4>;
        nit = new NitscheCouplType(
            bele_xyz, C_umum, C_usum, C_umus, C_usus, rhC_um, rhC_us, fldparaxfem);
        break;
      }
      default:
        FOUR_C_THROW("Unsupported boundary element shape {}", bele->shape());
        break;
    }
  }
  else
    FOUR_C_THROW("Unsupported number of {} nodes for coupling slave element.", numdofpernode);

  return std::shared_ptr<NitscheInterface<distype>>(nit);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
std::shared_ptr<NitscheInterface<distype>>
NitscheInterface<distype>::create_nitsche_coupling_two_sided(Core::Elements::Element* vele,
    Core::LinAlg::SerialDenseMatrix& vele_xyz, Core::LinAlg::SerialDenseMatrix& C_umum,
    Core::LinAlg::SerialDenseMatrix& C_usum, Core::LinAlg::SerialDenseMatrix& C_umus,
    Core::LinAlg::SerialDenseMatrix& C_usus, Core::LinAlg::SerialDenseMatrix& rhC_um,
    Core::LinAlg::SerialDenseMatrix& rhC_us,
    const Discret::Elements::FluidEleParameterXFEM& fldparaxfem)
{
  NitscheInterface* nit = nullptr;

  // get number of dofs for the embedded element
  const unsigned numdofpernode = vele->num_dof_per_node(*vele->nodes()[0]);

  if (numdofpernode == 4)
  {
    switch (vele->shape())
    {
      case Core::FE::CellType::hex8:
      {
        using NitscheCouplType = NitscheCoupling<distype, Core::FE::CellType::hex8, 4>;
        nit = new NitscheCouplType(
            vele_xyz, C_umum, C_usum, C_umus, C_usus, rhC_um, rhC_us, fldparaxfem);
        break;
      }
      case Core::FE::CellType::hex20:
      {
        using NitscheCouplType = NitscheCoupling<distype, Core::FE::CellType::hex20, 4>;
        nit = new NitscheCouplType(
            vele_xyz, C_umum, C_usum, C_umus, C_usus, rhC_um, rhC_us, fldparaxfem);
        break;
      }
      case Core::FE::CellType::hex27:
      {
        using NitscheCouplType = NitscheCoupling<distype, Core::FE::CellType::hex27, 4>;
        nit = new NitscheCouplType(
            vele_xyz, C_umum, C_usum, C_umus, C_usus, rhC_um, rhC_us, fldparaxfem);
        break;
      }
      default:
        FOUR_C_THROW("Unsupported volume element shape {}.", vele->shape());
        break;
    }
  }
  else if (numdofpernode == 3)
  {
    switch (vele->shape())
    {
      case Core::FE::CellType::hex8:
      {
        using NitscheCouplType = NitscheCoupling<distype, Core::FE::CellType::hex8, 3>;
        nit = new NitscheCouplType(
            vele_xyz, C_umum, C_usum, C_umus, C_usus, rhC_um, rhC_us, fldparaxfem);
        break;
      }
      default:
        // expecting 3 dofs per slave element node as this is fluid-solid coupling
        FOUR_C_THROW("Unsupported volume element shape {}.", vele->shape());
        break;
    }
  }

  return std::shared_ptr<NitscheInterface<distype>>(nit);
}

template <Core::FE::CellType distype>
std::shared_ptr<NitscheInterface<distype>>
NitscheInterface<distype>::create_nitsche_coupling_two_sided(Core::Elements::Element* vele,
    Core::LinAlg::SerialDenseMatrix& vele_xyz, Core::LinAlg::SerialDenseMatrix& C_umum,
    Core::LinAlg::SerialDenseMatrix& C_usum, Core::LinAlg::SerialDenseMatrix& C_umus,
    Core::LinAlg::SerialDenseMatrix& C_usus, Core::LinAlg::SerialDenseVector& rhC_um,
    Core::LinAlg::SerialDenseMatrix& rhC_us,
    const Discret::Elements::FluidEleParameterXFEM& fldparaxfem)
{
  NitscheInterface* nit = nullptr;

  // get number of dofs for the embedded element
  const unsigned numdofpernode = vele->num_dof_per_node(*vele->nodes()[0]);

  if (numdofpernode == 4)
  {
    switch (vele->shape())
    {
      case Core::FE::CellType::hex8:
      {
        using NitscheCouplType = NitscheCoupling<distype, Core::FE::CellType::hex8, 4>;
        nit = new NitscheCouplType(
            vele_xyz, C_umum, C_usum, C_umus, C_usus, rhC_um, rhC_us, fldparaxfem);
        break;
      }
      case Core::FE::CellType::hex20:
      {
        using NitscheCouplType = NitscheCoupling<distype, Core::FE::CellType::hex20, 4>;
        nit = new NitscheCouplType(
            vele_xyz, C_umum, C_usum, C_umus, C_usus, rhC_um, rhC_us, fldparaxfem);
        break;
      }
      case Core::FE::CellType::hex27:
      {
        using NitscheCouplType = NitscheCoupling<distype, Core::FE::CellType::hex27, 4>;
        nit = new NitscheCouplType(
            vele_xyz, C_umum, C_usum, C_umus, C_usus, rhC_um, rhC_us, fldparaxfem);
        break;
      }
      default:
        FOUR_C_THROW("Unsupported volume element shape {}.", vele->shape());
        break;
    }
  }
  else if (numdofpernode == 3)
  {
    switch (vele->shape())
    {
      case Core::FE::CellType::hex8:
      {
        using NitscheCouplType = NitscheCoupling<distype, Core::FE::CellType::hex8, 3>;
        nit = new NitscheCouplType(
            vele_xyz, C_umum, C_usum, C_umus, C_usus, rhC_um, rhC_us, fldparaxfem);
        break;
      }
      default:
        // expecting 3 dofs per slave element node as this is fluid-solid coupling
        FOUR_C_THROW("Unsupported volume element shape {}.", vele->shape());
        break;
    }
  }

  return std::shared_ptr<NitscheInterface<distype>>(nit);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
std::shared_ptr<HybridLMInterface<distype>>
HybridLMInterface<distype>::create_hybrid_lm_coupling_x_fluid_wdbc(
    bool is_viscAdjointSymmetric  ///< flag that indicates equal signs of Nitsche's standard &
                                  ///< adjoint viscous term
)
{
  HybridLMInterface* hybridlm = nullptr;
  using HybridLMCouplType = HybridLMCoupling<distype, Core::FE::CellType::dis_none, 3>;
  hybridlm = new HybridLMCouplType(is_viscAdjointSymmetric);

  return std::shared_ptr<HybridLMInterface<distype>>(hybridlm);
}
/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
std::shared_ptr<HybridLMInterface<distype>>
HybridLMInterface<distype>::create_hybrid_lm_coupling_x_fluid_wdbc(
    Core::Elements::Element* bele,              ///< boundary element
    Core::LinAlg::SerialDenseMatrix& bele_xyz,  ///< global node coordinates of boundary element
    bool is_viscAdjointSymmetric  ///< flag that indicates equal signs of Nitsche's standard &
                                  ///< adjoint viscous term
)
{
  // get number of dofs for this boundary element
  const unsigned numdofpernode = bele->num_dof_per_node(*bele->nodes()[0]);

  // three dofs per node, for standard Dirichlet coupling, four dofs per node for fluid-fluid
  // coupling
  if (numdofpernode != 3)
  {
    FOUR_C_THROW("Unsupported number of {} nodes for standard Dirichlet coupling.", numdofpernode);
  }

  switch (bele->shape())
  {
    case Core::FE::CellType::quad4:
    {
      using HybridLMCouplType = HybridLMCoupling<distype, Core::FE::CellType::quad4, 3>;
      return std::make_shared<HybridLMCouplType>(bele_xyz, is_viscAdjointSymmetric);
      break;
    }
    case Core::FE::CellType::quad8:
    {
      using HybridLMCouplType = HybridLMCoupling<distype, Core::FE::CellType::quad8, 3>;
      return std::make_shared<HybridLMCouplType>(bele_xyz, is_viscAdjointSymmetric);
      break;
    }
    case Core::FE::CellType::quad9:
    {
      using HybridLMCouplType = HybridLMCoupling<distype, Core::FE::CellType::quad9, 3>;
      return std::make_shared<HybridLMCouplType>(bele_xyz, is_viscAdjointSymmetric);
      break;
    }
    default:
      FOUR_C_THROW("Unsupported boundary element shape {}", bele->shape());
      break;
  }

  return nullptr;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
std::shared_ptr<HybridLMInterface<distype>>
HybridLMInterface<distype>::create_hybrid_lm_coupling_x_fluid_sided(
    Core::Elements::Element* bele,              ///< boundary element
    Core::LinAlg::SerialDenseMatrix& bele_xyz,  ///< global node coordinates of boundary element
    Core::LinAlg::SerialDenseMatrix& C_usum,    ///< C_usum coupling matrix
    Core::LinAlg::SerialDenseMatrix& C_umus,    ///< C_umus coupling matrix
    Core::LinAlg::SerialDenseMatrix& rhC_us,    ///< C_us coupling rhs
    Core::LinAlg::SerialDenseMatrix& G_s_us,    ///< \f$G_{u^s \sigma}\f$ coupling matrix
    Core::LinAlg::SerialDenseMatrix& G_us_s,    ///< \f$G_{\sigma u^s}\f$ coupling matrix
    bool is_viscAdjointSymmetric  ///< flag that indicates equal signs of Nitsche's standard &
                                  ///< adjoint viscous term
)
{
  HybridLMInterface* hlm = nullptr;

  // get number of dofs for this boundary element
  const unsigned numdofpernode = bele->num_dof_per_node(*bele->nodes()[0]);

  // three dofs per node, for monolithic XFSI, four dofs per node for backgrounds-sided fluid-fluid
  // coupling
  if (numdofpernode == 3)
  {
    switch (bele->shape())
    {
      case Core::FE::CellType::quad4:
      {
        using HybridLMCouplType = HybridLMCoupling<distype, Core::FE::CellType::quad4, 3>;
        hlm = new HybridLMCouplType(
            bele_xyz, C_usum, C_umus, rhC_us, G_s_us, G_us_s, is_viscAdjointSymmetric);
        break;
      }
      case Core::FE::CellType::quad8:
      {
        using HybridLMCouplType = HybridLMCoupling<distype, Core::FE::CellType::quad8, 3>;
        hlm = new HybridLMCouplType(
            bele_xyz, C_usum, C_umus, rhC_us, G_s_us, G_us_s, is_viscAdjointSymmetric);
        break;
      }
      case Core::FE::CellType::quad9:
      {
        using HybridLMCouplType = HybridLMCoupling<distype, Core::FE::CellType::quad9, 3>;
        hlm = new HybridLMCouplType(
            bele_xyz, C_usum, C_umus, rhC_us, G_s_us, G_us_s, is_viscAdjointSymmetric);
        break;
      }
      default:
        FOUR_C_THROW("Unsupported boundary element shape {}", bele->shape());
        break;
    }
  }
  else if (numdofpernode == 4)
  {
    switch (bele->shape())
    {
      case Core::FE::CellType::quad4:
      {
        using HybridLMCouplType = HybridLMCoupling<distype, Core::FE::CellType::quad4, 4>;
        hlm = new HybridLMCouplType(
            bele_xyz, C_usum, C_umus, rhC_us, G_s_us, G_us_s, is_viscAdjointSymmetric);
        break;
      }
      case Core::FE::CellType::quad8:
      {
        using HybridLMCouplType = HybridLMCoupling<distype, Core::FE::CellType::quad8, 4>;
        hlm = new HybridLMCouplType(
            bele_xyz, C_usum, C_umus, rhC_us, G_s_us, G_us_s, is_viscAdjointSymmetric);
        break;
      }
      case Core::FE::CellType::quad9:
      {
        using HybridLMCouplType = HybridLMCoupling<distype, Core::FE::CellType::quad9, 4>;
        hlm = new HybridLMCouplType(
            bele_xyz, C_usum, C_umus, rhC_us, G_s_us, G_us_s, is_viscAdjointSymmetric);
        break;
      }
      default:
        FOUR_C_THROW("Unsupported boundary element shape {}", bele->shape());
        break;
    }
  }
  else
    FOUR_C_THROW("Unsupported number of {} nodes for coupling slave element.", numdofpernode);

  return std::shared_ptr<HybridLMInterface<distype>>(hlm);
}

template class Discret::Elements::XFLUID::SlaveElementInterface<Core::FE::CellType::hex8>;
template class Discret::Elements::XFLUID::SlaveElementInterface<Core::FE::CellType::hex20>;
template class Discret::Elements::XFLUID::SlaveElementInterface<Core::FE::CellType::hex27>;
template class Discret::Elements::XFLUID::SlaveElementInterface<Core::FE::CellType::tet4>;
template class Discret::Elements::XFLUID::SlaveElementInterface<Core::FE::CellType::tet10>;
template class Discret::Elements::XFLUID::SlaveElementInterface<Core::FE::CellType::wedge6>;
template class Discret::Elements::XFLUID::SlaveElementInterface<Core::FE::CellType::wedge15>;
// template class
// Discret::Elements::XFLUID::SlaveElementInterface<Core::FE::CellType::pyramid5>;

template class Discret::Elements::XFLUID::NitscheInterface<Core::FE::CellType::hex8>;
template class Discret::Elements::XFLUID::NitscheInterface<Core::FE::CellType::hex20>;
template class Discret::Elements::XFLUID::NitscheInterface<Core::FE::CellType::hex27>;
template class Discret::Elements::XFLUID::NitscheInterface<Core::FE::CellType::tet4>;
template class Discret::Elements::XFLUID::NitscheInterface<Core::FE::CellType::tet10>;
template class Discret::Elements::XFLUID::NitscheInterface<Core::FE::CellType::wedge6>;
template class Discret::Elements::XFLUID::NitscheInterface<Core::FE::CellType::wedge15>;
// template class
// Discret::Elements::XFLUID::NitscheInterface<Core::FE::CellType::pyramid5>;

template class Discret::Elements::XFLUID::HybridLMInterface<Core::FE::CellType::hex8>;
template class Discret::Elements::XFLUID::HybridLMInterface<Core::FE::CellType::hex20>;
template class Discret::Elements::XFLUID::HybridLMInterface<Core::FE::CellType::hex27>;
template class Discret::Elements::XFLUID::HybridLMInterface<Core::FE::CellType::tet4>;
template class Discret::Elements::XFLUID::HybridLMInterface<Core::FE::CellType::tet10>;
template class Discret::Elements::XFLUID::HybridLMInterface<Core::FE::CellType::wedge6>;
template class Discret::Elements::XFLUID::HybridLMInterface<Core::FE::CellType::wedge15>;
// template class
// Discret::Elements::XFLUID::HybridLMInterface<Core::FE::CellType::pyramid5>;

FOUR_C_NAMESPACE_CLOSE
