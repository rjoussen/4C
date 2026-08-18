// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_config.hpp"

#include "4C_reduced_lung_1d_pipe_flow_resulttest.hpp"

#include "4C_comm_mpi_utils.hpp"
#include "4C_fem_general_node.hpp"
#include "4C_io_input_parameter_container.hpp"
#include "4C_reduced_lung_1d_pipe_flow_terminal_unit.hpp"
#include "4C_utils_result_test.hpp"

#include <cmath>
#include <numbers>

FOUR_C_NAMESPACE_OPEN
// Constructor
ReducedLung1dPipeFlow::ResultTest::ResultTest(std::shared_ptr<Core::FE::Discretization> dis,
    std::shared_ptr<const Core::LinAlg::Vector<double>> sol,
    std::vector<ReducedLung1DPipe::TerminalUnit::TerminalUnitModel> terminal_units,
    std::unordered_map<int, std::size_t> tu_index_map)
    : Core::Utils::ResultTest("ARTNET"),
      dis_(std::move(dis)),
      sol_(std::move(sol)),
      terminal_units_(std::move(terminal_units)),
      tu_index_map_(std::move(tu_index_map))
{
}

void ReducedLung1dPipeFlow::ResultTest::test_node(
    const Core::IO::InputParameterContainer& container, int& nerr, int& test_count)
{
  std::string dis = container.get<std::string>("DIS");
  if (dis != dis_->name())
  {
    return;
  }

  int gid = container.get<int>("NODE") - 1;
  const int is_local_node = dis_->have_global_node(gid);

  if (const int is_node_of_any_rank = Core::Communication::sum_all(is_local_node, dis_->get_comm());
      is_node_of_any_rank == 0)
  {
    FOUR_C_ASSERT(is_node_of_any_rank != 0, "Node {} does not belong to discretization {}", gid + 1,
        dis_->name());
  }

  if (!dis_->have_global_node(gid)) return;

  Core::Nodes::Node* node_to_check = dis_->g_node(gid);
  if (node_to_check->owner() != Core::Communication::my_mpi_rank(dis_->get_comm())) return;

  const Core::LinAlg::Map& dofmap = sol_->get_map();

  std::string quantity = container.get<std::string>("QUANTITY");

  double result = 0.0;

  if (quantity == "area")
  {
    const int dof = dis_->dof(node_to_check, 0);  // A
    result = sol_->local_values_as_span()[dofmap.lid(dof)];
  }
  else if (quantity == "velocity")
  {
    const int dof = dis_->dof(node_to_check, 1);  // u
    result = sol_->local_values_as_span()[dofmap.lid(dof)];
  }
  else if (quantity == "flow")
  {
    const int dof_A = dis_->dof(node_to_check, 0);
    const int dof_u = dis_->dof(node_to_check, 1);  // u
    result = sol_->local_values_as_span()[dofmap.lid(dof_A)] *
             sol_->local_values_as_span()[dofmap.lid(dof_u)];
  }
  else if (quantity == "impedance" || quantity == "impedance_f_min")
  {
    const auto& index_it = tu_index_map_.find(gid);

    FOUR_C_ASSERT(index_it != tu_index_map_.end(),
        "Node {} has no terminal unit; impedance quantities require a structured-tree "
        "terminal unit",
        gid + 1);

    const auto* structured_tree =
        std::get_if<ReducedLung1DPipe::TerminalUnit::StructuredTreeTerminalUnit>(
            &terminal_units_[index_it->second].model);

    FOUR_C_ASSERT(structured_tree != nullptr,
        "Terminal unit at node {} is not a structured-tree model", gid + 1);

    const auto& spectrum = structured_tree->impedance_spectrum;
    FOUR_C_ASSERT(spectrum.size() > 2,
        "Structured tree at node {} has no impedance spectrum (size {})", gid + 1, spectrum.size());

    // The impedance at zero frequency is a positive real resistance.
    const double z_dc = std::abs(spectrum.front());

    if (quantity == "impedance")
    {
      result = z_dc;
    }
    else  // impedance minimum
    {
      double min_abs = std::abs(spectrum[1]);
      int min_k = 1;
      const std::size_t n_half = spectrum.size() / 2;
      for (std::size_t k = 2; k <= n_half; ++k)
      {
        const double abs_z = std::abs(spectrum[k]);
        if (abs_z < min_abs)
        {
          min_abs = abs_z;
          min_k = static_cast<int>(k);
        }
      }
      result = 2.0 * std::numbers::pi * min_k / (structured_tree->n_cycle * structured_tree->dt);
    }
  }
  else
  {
    FOUR_C_THROW("Unsupported QUANTITY '{}'", quantity);
  }

  nerr += compare_values(result, "NODE", container);
  test_count++;
}

FOUR_C_NAMESPACE_CLOSE