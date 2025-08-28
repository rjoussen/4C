// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_particle_wall_discretization_runtime_vtu_writer.hpp"

#include "4C_fem_discretization.hpp"
#include "4C_fem_general_element.hpp"
#include "4C_fem_general_node.hpp"
#include "4C_global_data.hpp"
#include "4C_io.hpp"
#include "4C_io_control.hpp"
#include "4C_io_discretization_visualization_writer_mesh.hpp"
#include "4C_io_visualization_parameters.hpp"
#include "4C_particle_wall_datastate.hpp"

#include <memory>

FOUR_C_NAMESPACE_OPEN


/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
PARTICLEWALL::WallDiscretizationRuntimeVtuWriter::WallDiscretizationRuntimeVtuWriter(
    const std::shared_ptr<Core::FE::Discretization> walldiscretization,
    const std::shared_ptr<PARTICLEWALL::WallDataState> walldatastate, const double restart_time)
    : walldiscretization_(walldiscretization), walldatastate_(walldatastate)
{
  // construct the writer object
  runtime_vtuwriter_ =
      std::make_unique<Core::IO::DiscretizationVisualizationWriterMesh>(walldiscretization,
          Core::IO::visualization_parameters_factory(
              Global::Problem::instance()->io_params().sublist("RUNTIME VTK OUTPUT"),
              *Global::Problem::instance()->output_control_file(), restart_time));
}

void PARTICLEWALL::WallDiscretizationRuntimeVtuWriter::write_wall_discretization_runtime_output(
    const int step, const double time) const
{
  // reset the writer object
  runtime_vtuwriter_->reset();

  // node displacements
  {
    if (walldatastate_->get_disp_col() != nullptr)
    {
      std::vector<std::optional<std::string>> context(3, "disp");
      runtime_vtuwriter_->append_result_data_vector_with_context(
          *walldatastate_->get_ref_disp_col(), Core::IO::OutputEntity::dof, context);
    }
  }

  // node owner
  {
    auto nodeowner = Core::LinAlg::Vector<double>(*walldiscretization_->node_col_map(), true);
    for (int inode = 0; inode < walldiscretization_->num_my_col_nodes(); ++inode)
    {
      const Core::Nodes::Node* node = walldiscretization_->l_col_node(inode);
      (nodeowner).get_values()[inode] = node->owner();
    }
    runtime_vtuwriter_->append_result_data_vector_with_context(
        nodeowner, Core::IO::OutputEntity::node, {"owner"});
  }

  // element owner
  {
    runtime_vtuwriter_->append_element_owner("owner");
  }

  // element id
  {
    auto eleid = Core::LinAlg::Vector<double>(*walldiscretization_->element_row_map(), true);
    for (int iele = 0; iele < walldiscretization_->num_my_row_elements(); ++iele)
    {
      const Core::Elements::Element* ele = walldiscretization_->l_row_element(iele);
      (eleid).get_values()[iele] = ele->id();
    }
    runtime_vtuwriter_->append_result_data_vector_with_context(
        eleid, Core::IO::OutputEntity::element, {"id"});
  }

  // finalize everything and write all required files to filesystem
  runtime_vtuwriter_->write_to_disk(time, step);
}

FOUR_C_NAMESPACE_CLOSE
