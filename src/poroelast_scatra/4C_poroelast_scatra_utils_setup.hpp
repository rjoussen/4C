// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_POROELAST_SCATRA_UTILS_SETUP_HPP
#define FOUR_C_POROELAST_SCATRA_UTILS_SETUP_HPP

#include "4C_config.hpp"

#include "4C_fem_dofset_gidbased_wrapper.hpp"
#include "4C_fem_dofset_predefineddofnumber.hpp"
#include "4C_fem_general_utils_createdis.hpp"
#include "4C_global_data.hpp"
#include "4C_poroelast_scatra_utils.hpp"
#include "4C_poroelast_scatra_utils_clonestrategy.hpp"
#include "4C_poroelast_utils_setup.hpp"
#include "4C_rebalance_binning_based.hpp"
#include "4C_solid_3D_ele.hpp"
#include "4C_utils_parameter_list.hpp"

FOUR_C_NAMESPACE_OPEN


namespace PoroElastScaTra
{
  namespace Utils
  {
    //! setup discretization, includes cloning the structure discretization
    template <class PoroCloneStrategy, class PoroScatraCloneStrategy>
    void setup_poro_scatra_discretizations()
    {
      // Scheme    : the structure discretization is received from the input. Then, an ale-fluid
      // disc.is cloned from the struct. one.
      //  After that, an ale-scatra disc. is cloned from the structure discretization.

      Global::Problem* problem = Global::Problem::instance();

      // 1.-Initialization.
      std::shared_ptr<Core::FE::Discretization> structdis = problem->get_dis("structure");
      std::shared_ptr<Core::FE::Discretization> fluiddis = problem->get_dis("porofluid");
      std::shared_ptr<Core::FE::Discretization> scatradis = problem->get_dis("scatra");

      // setup of the discretizations, including clone strategy (do not set material pointers, this
      // will be done here)
      PoroElast::Utils::setup_poro<PoroCloneStrategy>(false);

      // 3.-Access the scatra discretization, make sure it's empty, and fill it by cloning the
      // structural one.
      if (fluiddis->num_global_nodes() == 0) FOUR_C_THROW("Fluid discretization is empty!");

      if (!scatradis->filled()) scatradis->fill_complete();

      if (scatradis->num_global_nodes() == 0)
      {
        // fill scatra discretization by cloning structure discretization
        Core::FE::clone_discretization<PoroScatraCloneStrategy>(
            *structdis, *scatradis, Global::Problem::instance()->cloning_material_map());
        scatradis->fill_complete();

        // assign materials. Order is important here!
        PoroElast::Utils::set_material_pointers_matching_grid(*structdis, *fluiddis);
        PoroElast::Utils::set_material_pointers_matching_grid(*structdis, *scatradis);
        PoroElast::Utils::set_material_pointers_matching_grid(*fluiddis, *scatradis);

        // the problem is two way coupled, thus each discretization must know the other
        // discretization

        // build a proxy of the structure discretization for the scatra field
        std::shared_ptr<Core::DOFSets::DofSetInterface> structdofset =
            structdis->get_dof_set_proxy();
        // build a proxy of the fluid discretization for the scatra field
        std::shared_ptr<Core::DOFSets::DofSetInterface> fluiddofset = fluiddis->get_dof_set_proxy();
        // build a proxy of the fluid discretization for the structure/fluid field
        std::shared_ptr<Core::DOFSets::DofSetInterface> scatradofset =
            scatradis->get_dof_set_proxy();

        // check if ScatraField has 2 discretizations, so that coupling is possible
        if (scatradis->add_dof_set(structdofset) != 1)
          FOUR_C_THROW("unexpected dof sets in scatra field");
        if (scatradis->add_dof_set(fluiddofset) != 2)
          FOUR_C_THROW("unexpected dof sets in scatra field");
        if (structdis->add_dof_set(scatradofset) != 2)
          FOUR_C_THROW("unexpected dof sets in structure field");
        if (fluiddis->add_dof_set(scatradofset) != 2)
          FOUR_C_THROW("unexpected dof sets in fluid field");

        structdis->fill_complete();
        fluiddis->fill_complete();
        scatradis->fill_complete();
      }
      else
      {
        // create vector of discr.
        std::vector<std::shared_ptr<Core::FE::Discretization>> dis;
        dis.push_back(structdis);
        dis.push_back(fluiddis);
        dis.push_back(scatradis);

        Teuchos::ParameterList binning_params =
            Global::Problem::instance()->binning_strategy_params();
        Core::Utils::add_enum_class_to_parameter_list<Core::FE::ShapeFunctionType>(
            "spatial_approximation_type", Global::Problem::instance()->spatial_approximation_type(),
            binning_params);
        Core::Rebalance::rebalance_discretizations_by_binning(binning_params,
            Global::Problem::instance()->output_control_file(), dis, nullptr, nullptr, false);

        // set material pointers
        PoroElast::Utils::set_material_pointers_matching_grid(*structdis, *fluiddis);

        // first call fill_complete for single discretizations.
        // This way the physical dofs are numbered successively
        structdis->fill_complete();
        fluiddis->fill_complete();
        scatradis->fill_complete();

        // build auxiliary dofsets, i.e. pseudo dofs on each discretization
        const int ndofpernode_fluid = fluiddis->num_dof(0, fluiddis->l_row_node(0));
        const int ndofperelement_fluid = 0;
        const int ndofpernode_struct = structdis->num_dof(0, structdis->l_row_node(0));
        const int ndofperelement_struct = 0;
        const int ndofpernode_scatra = scatradis->num_dof(0, scatradis->l_row_node(0));
        const int ndofperelement_scatra = 0;

        std::shared_ptr<Core::DOFSets::DofSetInterface> dofsetaux;
        dofsetaux = std::make_shared<Core::DOFSets::DofSetPredefinedDoFNumber>(
            ndofpernode_scatra, ndofperelement_scatra, 0, true);
        if (structdis->add_dof_set(dofsetaux) != 2)
          FOUR_C_THROW("unexpected dof sets in structure field");
        dofsetaux = std::make_shared<Core::DOFSets::DofSetPredefinedDoFNumber>(
            ndofpernode_scatra, ndofperelement_scatra, 0, true);
        if (fluiddis->add_dof_set(dofsetaux) != 2)
          FOUR_C_THROW("unexpected dof sets in fluid field");
        dofsetaux = std::make_shared<Core::DOFSets::DofSetPredefinedDoFNumber>(
            ndofpernode_struct, ndofperelement_struct, 0, true);
        if (scatradis->add_dof_set(dofsetaux) != 1)
          FOUR_C_THROW("unexpected dof sets in scatra field");
        dofsetaux = std::make_shared<Core::DOFSets::DofSetPredefinedDoFNumber>(
            ndofpernode_fluid, ndofperelement_fluid, 0, true);
        if (scatradis->add_dof_set(dofsetaux) != 2)
          FOUR_C_THROW("unexpected dof sets in scatra field");

        // call assign_degrees_of_freedom also for auxiliary dofsets
        // note: the order of fill_complete() calls determines the gid numbering!
        // 1. structure dofs
        // 2. fluiddis dofs
        // 3. scatradis dofs
        // 4. auxiliary dofs
        structdis->fill_complete({
            .assign_degrees_of_freedom = true,
            .init_elements = false,
            .do_boundary_conditions = false,
        });
        fluiddis->fill_complete({
            .assign_degrees_of_freedom = true,
            .init_elements = false,
            .do_boundary_conditions = false,
        });
        scatradis->fill_complete({
            .assign_degrees_of_freedom = true,
            .init_elements = false,
            .do_boundary_conditions = false,
        });
      }
    }
  }  // namespace Utils
}  // namespace PoroElastScaTra

FOUR_C_NAMESPACE_CLOSE

#endif
