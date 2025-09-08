// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_inpar_volmortar.hpp"

#include "4C_coupling_volmortar.hpp"
#include "4C_io_input_spec_builders.hpp"
FOUR_C_NAMESPACE_OPEN



Core::IO::InputSpec Inpar::VolMortar::set_valid_parameters()
{
  using namespace Core::IO::InputSpecBuilders;

  /* parameters for volmortar */
  Core::IO::InputSpec spec = group("VOLMORTAR COUPLING",
      {

          deprecated_selection<Coupling::VolMortar::IntType>("INTTYPE",
              {
                  {"Elements", Coupling::VolMortar::inttype_elements},
                  {"elements", Coupling::VolMortar::inttype_elements},
                  {"Segments", Coupling::VolMortar::inttype_segments},
                  {"segments", Coupling::VolMortar::inttype_segments},
              },
              {.description = "Type of numerical integration scheme",
                  .default_value = Coupling::VolMortar::inttype_elements}),


          deprecated_selection<Coupling::VolMortar::CouplingType>("COUPLINGTYPE",
              {
                  {"Volmortar", Coupling::VolMortar::couplingtype_volmortar},
                  {"volmortar", Coupling::VolMortar::couplingtype_volmortar},
                  {"consistentinterpolation", Coupling::VolMortar::couplingtype_coninter},
                  {"consint", Coupling::VolMortar::couplingtype_coninter},
              },
              {.description = "Type of coupling",
                  .default_value = Coupling::VolMortar::couplingtype_volmortar}),

          deprecated_selection<Coupling::VolMortar::Shapefcn>("SHAPEFCN",
              {
                  {"Dual", Coupling::VolMortar::shape_dual},
                  {"dual", Coupling::VolMortar::shape_dual},
                  {"Standard", Coupling::VolMortar::shape_std},
                  {"standard", Coupling::VolMortar::shape_std},
                  {"std", Coupling::VolMortar::shape_std},
              },
              {.description = "Type of employed set of shape functions",
                  .default_value = Coupling::VolMortar::shape_dual}),

          deprecated_selection<Coupling::VolMortar::CutType>("CUTTYPE",
              {
                  {"dd", Coupling::VolMortar::cuttype_directdivergence},
                  {"directdivergence", Coupling::VolMortar::cuttype_directdivergence},
                  {"DirectDivergence", Coupling::VolMortar::cuttype_directdivergence},
                  {"tessellation", Coupling::VolMortar::cuttype_tessellation},
                  {"t", Coupling::VolMortar::cuttype_tessellation},
                  {"Tessellation", Coupling::VolMortar::cuttype_tessellation},
              },
              {.description = "Type of cut procedure/ integration point calculation",
                  .default_value = Coupling::VolMortar::cuttype_directdivergence}),

          deprecated_selection<Coupling::VolMortar::DualQuad>("DUALQUAD",
              {
                  {"nm", Coupling::VolMortar::dualquad_no_mod},
                  {"nomod", Coupling::VolMortar::dualquad_no_mod},
                  {"lm", Coupling::VolMortar::dualquad_lin_mod},
                  {"lin_mod", Coupling::VolMortar::dualquad_lin_mod},
                  {"qm", Coupling::VolMortar::dualquad_quad_mod},
                  {"quad_mod", Coupling::VolMortar::dualquad_quad_mod},
              },
              {.description =
                      "Type of dual shape function for weighting function for quadr. problems",
                  .default_value = Coupling::VolMortar::dualquad_no_mod}),

          parameter<bool>(
              "MESH_INIT", {.description = "If chosen, mesh initialization procedure is performed",
                               .default_value = false}),

          parameter<bool>("KEEP_EXTENDEDGHOSTING",
              {.description = "If chosen, extended ghosting is kept for simulation",
                  .default_value = true})},
      {.required = false});
  return spec;
}

FOUR_C_NAMESPACE_CLOSE