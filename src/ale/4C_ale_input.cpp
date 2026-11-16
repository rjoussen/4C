// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_ale_input.hpp"

#include "4C_fem_condition_definition.hpp"
#include "4C_io_input_spec_builders.hpp"

FOUR_C_NAMESPACE_OPEN



Core::IO::InputSpec ALE::valid_parameters()
{
  using namespace Core::IO::InputSpecBuilders;
  Core::IO::InputSpec spec = group("ALE DYNAMIC",
      {

          parameter<double>("TIMESTEP", {.description = "time step size", .default_value = 0.1}),

          parameter<int>(
              "NUMSTEP", {.description = "max number of time steps", .default_value = 41}),

          parameter<double>(
              "MAXTIME", {.description = "max simulation time", .default_value = 4.0}),

          parameter<ALE::AleDynamic>(
              "ALE_TYPE", {.description = "ale mesh movement algorithm", .default_value = solid}),

          parameter<double>("MESH_QUALITY_SIZE_WEIGHT",
              {.description =
                      "Weight w of the element size term relative to the element shape term in the "
                      "mesh quality potential (ALE_TYPE mesh_quality). For small deviations from "
                      "the "
                      "target, the potential acts like a linear elastic material with the Poisson "
                      "ratio (9w - 2) / (2 (9w + 1)), e.g., 0 for w = 2/9 and 0.35 for w = 1; "
                      "larger "
                      "values keep the element sizes more firmly. The size term or the motion "
                      "term (MESH_QUALITY_MOTION_WEIGHT) has to be active: the shape term alone "
                      "does not determine the sizes of the elements and the positions of nodes on "
                      "sliding boundaries.",
                  .default_value = 1.0}),

          parameter<ALE::MeshQualitySizeTarget>("MESH_QUALITY_SIZE_TARGET",
              {.description =
                      "Target of the element size in the mesh quality potential (ALE_TYPE "
                      "mesh_quality): initial_mesh, each element keeps its size in the initial "
                      "mesh; function, the size h(x, t) given by MESH_QUALITY_SIZE_FUNCTION at the "
                      "current position of each node, i.e., the target volume of an element is "
                      "h^3 and only its shape is that of the initial mesh (only supported by the "
                      "mesh regularization of solids). E.g., for plane strain with one element "
                      "layer of thickness t and an in-plane element size s, h = (s^2 t)^(1/3).",
                  .default_value = ALE::MeshQualitySizeTarget::initial_mesh}),

          parameter<int>("MESH_QUALITY_SIZE_FUNCTION",
              {.description = "Function h(x, t) of the target element size for "
                              "MESH_QUALITY_SIZE_TARGET function: the edge length of a cube with "
                              "the target volume (positive and continuously differentiable in "
                              "space).",
                  .default_value = -1}),

          parameter<double>("MESH_QUALITY_MOTION_WEIGHT",
              {.description =
                      "Weight of the motion term relative to the element shape term in the mesh "
                      "quality potential (ALE_TYPE mesh_quality): the squared displacement of the "
                      "nodes relative to the element size. It keeps the nodes close to the "
                      "reference configuration of the mesh movement, which is the current "
                      "configuration for the mesh regularization of solids, so that the mesh "
                      "stays close to the material and is only moved where its quality requires "
                      "it. A small weight (e.g., 0.001) without the size term makes each mesh "
                      "movement unique while the element sizes adapt to the shapes, e.g., become "
                      "small around tool edges.",
                  .default_value = 0.0}),

          parameter<bool>("ASSESSMESHQUALITY",
              {.description = "Evaluate element quality measure according to [Oddy et al. 1988]",
                  .default_value = false}),

          parameter<bool>(
              "UPDATEMATRIX", {.description = "Update stiffness matrix in every time step "
                                              "(only for linear/material strategies)",
                                  .default_value = false}),

          parameter<int>("MAXITER",
              {.description = "Maximum number of newton iterations.", .default_value = 1}),
          parameter<double>(
              "TOLRES", {.description = "Absolute tolerance for length scaled L2 residual norm ",
                            .default_value = 1.0e-06}),
          parameter<double>(
              "TOLDISP", {.description = "Absolute tolerance for length scaled L2 increment norm ",
                             .default_value = 1.0e-06}),

          parameter<int>("NUM_INITSTEP", {.description = "", .default_value = 0}),
          parameter<int>("RESTARTEVERY",
              {.description = "write restart data every RESTARTEVERY steps", .default_value = 1}),
          parameter<int>("RESULTSEVERY",
              {.description = "write results every RESULTSTEVERY steps", .default_value = 0}),
          deprecated_selection<ALE::DivContAct>("DIVERCONT",
              {
                  {"stop", divcont_stop},
                  {"continue", divcont_continue},
              },
              {.description = "What to do if nonlinear solver does not converge?",
                  .default_value = divcont_continue}),

          deprecated_selection<ALE::MeshTying>("MESHTYING",
              {
                  {"no", no_meshtying},
                  {"meshtying", meshtying},
                  {"meshsliding", meshsliding},
              },
              {.description = "Flag to (de)activate mesh tying and mesh sliding algorithm",
                  .default_value = no_meshtying}),

          // Initial displacement
          deprecated_selection<ALE::InitialDisp>("INITIALDISP",
              {
                  {"zero_displacement", initdisp_zero_disp},
                  {"displacement_by_function", initdisp_disp_by_function},
              },
              {.description = "Initial displacement for structure problem",
                  .default_value = initdisp_zero_disp}),

          // Function to evaluate initial displacement
          parameter<int>("STARTFUNCNO",
              {.description = "Function for Initial displacement", .default_value = -1}),

          // linear solver id used for scalar ale problems
          parameter<int>(
              "LINEAR_SOLVER", {.description = "number of linear solver used for ale problems...",
                                   .default_value = -1})},
      {.required = false});
  return spec;
}



void ALE::set_valid_conditions(std::vector<Core::Conditions::ConditionDefinition>& condlist)
{
  using namespace Core::IO::InputSpecBuilders;

  /*--------------------------------------------------------------------*/
  // Ale update boundary condition

  Core::Conditions::ConditionDefinition linealeupdate("DESIGN ALE UPDATE LINE CONDITIONS",
      "ALEUPDATECoupling", "ALEUPDATE Coupling", Core::Conditions::ALEUPDATECoupling, true,
      Core::Conditions::geometry_type_line);
  Core::Conditions::ConditionDefinition surfaleupdate("DESIGN ALE UPDATE SURF CONDITIONS",
      "ALEUPDATECoupling", "ALEUPDATE Coupling", Core::Conditions::ALEUPDATECoupling, true,
      Core::Conditions::geometry_type_surface);

  const auto make_ale_update = [&condlist](Core::Conditions::ConditionDefinition& cond)
  {
    cond.add_component(deprecated_selection<std::string>("COUPLING",
        {"lagrange", "heightfunction", "sphereHeightFunction", "meantangentialvelocity",
            "meantangentialvelocityscaled"},
        {.description = "", .default_value = "lagrange"}));
    cond.add_component(parameter<double>("VAL"));
    cond.add_component(parameter<int>("NODENORMALFUNCT"));

    condlist.emplace_back(cond);
  };

  make_ale_update(linealeupdate);
  make_ale_update(surfaleupdate);
}

FOUR_C_NAMESPACE_CLOSE