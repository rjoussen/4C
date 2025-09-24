// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_contact_rough_node.hpp"

#include "4C_comm_pack_helpers.hpp"
#include "4C_contact_defines.hpp"
#include "4C_contact_element.hpp"
#include "4C_global_data.hpp"
#include "4C_utils_function.hpp"

#ifdef FOUR_C_WITH_MIRCO

#include <mirco_topology.h>
#include <mirco_topologyutilities.h>

#endif

FOUR_C_NAMESPACE_OPEN

CONTACT::RoughNodeType CONTACT::RoughNodeType::instance_;

Core::Communication::ParObject* CONTACT::RoughNodeType::create(
    Core::Communication::UnpackBuffer& buffer)
{
  std::vector<double> x(3, 0.0);
  std::vector<int> dofs;

  CONTACT::RoughNode* node = new CONTACT::RoughNode(0, x, 0, dofs, false, false, 0, 0, 0, 0, 0, 0);
  node->unpack(buffer);

  return node;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
CONTACT::RoughNode::RoughNode(int id, std::span<const double> coords, const int owner,
    const std::vector<int>& dofs, const bool isslave, const bool initactive,
    const int hurstexponentfunction, int initialtopologystddeviationfunction, int resolution,
    bool randomtopologyflag, bool randomseedflag, int randomgeneratorseed)
    : CONTACT::Node(id, coords, owner, dofs, isslave, initactive),
      hurstexponentfunction_(hurstexponentfunction),
      initialtopologystddeviationfunction_(initialtopologystddeviationfunction),
      resolution_(resolution),
      randomtopologyflag_(randomtopologyflag),
      randomseedflag_(randomseedflag),
      randomgeneratorseed_(randomgeneratorseed)
{
#ifdef FOUR_C_WITH_MIRCO
  if (isslave)
  {
    hurstExponent_ = Global::Problem::instance()
                         ->function_by_id<Core::Utils::FunctionOfSpaceTime>(hurstexponentfunction_)
                         .evaluate(this->x(), 1, this->n_dim());
    initialTopologyStdDeviation_ =
        Global::Problem::instance()
            ->function_by_id<Core::Utils::FunctionOfSpaceTime>(initialtopologystddeviationfunction_)
            .evaluate(this->x(), 1, this->n_dim());

    const int N = pow(2, resolution_);
    topology_.shape(N + 1, N + 1);

    std::string topologyFilePath = "";
    Teuchos::RCP<MIRCO::TopologyGeneration> surfacegenerator;
    // creating the correct surface object
    MIRCO::CreateSurfaceObject(resolution_, initialTopologyStdDeviation_, hurstExponent_,
        randomseedflag_, topologyFilePath, randomtopologyflag_, randomgeneratorseed_,
        surfacegenerator);
    surfacegenerator->GetSurface(topology_.base());

    auto max_and_mean = MIRCO::ComputeMaxAndMean(topology_.base());
    maxTopologyHeight_ = max_and_mean.max_;
  }
#else
  FOUR_C_THROW(
      "You are trying to create a RoughNode with FOUR_C_WITH_MIRCO flag turned off. Please enable "
      "this flag and build 4C again");
#endif
}


/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 *----------------------------------------------------------------------*/
void CONTACT::RoughNode::pack(Core::Communication::PackBuffer& data) const
{
  // pack type of this instance of ParObject
  int type = unique_par_object_id();
  add_to_pack(data, type);

  // add base class Mortar::Node
  CONTACT::Node::pack(data);

  add_to_pack(data, hurstexponentfunction_);
  add_to_pack(data, initialtopologystddeviationfunction_);
  add_to_pack(data, resolution_);
  add_to_pack(data, randomtopologyflag_);
  add_to_pack(data, randomseedflag_);
  add_to_pack(data, randomgeneratorseed_);

  add_to_pack(data, hurstExponent_);
  add_to_pack(data, initialTopologyStdDeviation_);
  add_to_pack(data, topology_);
  add_to_pack(data, maxTopologyHeight_);

  return;
}

/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 *----------------------------------------------------------------------*/
void CONTACT::RoughNode::unpack(Core::Communication::UnpackBuffer& buffer)
{
  Core::Communication::extract_and_assert_id(buffer, unique_par_object_id());

  // extract base class CONTACT::Node
  CONTACT::Node::unpack(buffer);

  extract_from_pack(buffer, hurstexponentfunction_);
  extract_from_pack(buffer, initialtopologystddeviationfunction_);
  extract_from_pack(buffer, resolution_);
  extract_from_pack(buffer, randomtopologyflag_);
  extract_from_pack(buffer, randomseedflag_);
  extract_from_pack(buffer, randomgeneratorseed_);

  extract_from_pack(buffer, hurstExponent_);
  extract_from_pack(buffer, initialTopologyStdDeviation_);
  extract_from_pack(buffer, topology_);
  extract_from_pack(buffer, maxTopologyHeight_);

  // Check

  return;
}

FOUR_C_NAMESPACE_CLOSE
