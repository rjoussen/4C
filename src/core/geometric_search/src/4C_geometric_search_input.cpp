// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_geometric_search_input.hpp"

#include "4C_io_input_spec_builders.hpp"

FOUR_C_NAMESPACE_OPEN

Core::IO::InputSpec Core::GeometricSearch::valid_parameters()
{
  using namespace Core::IO::InputSpecBuilders;
  Core::IO::InputSpec spec = group("BOUNDINGVOLUME STRATEGY",
      {

          parameter<double>("BEAM_RADIUS_EXTENSION_FACTOR",
              {.description =
                      "Beams radius is multiplied with the factor and then the bounding volume "
                      "only depending on the beam centerline is extended in all directions (+ and "
                      "-) by that value.",
                  .default_value = 2.0}),

          parameter<double>("SPHERE_RADIUS_EXTENSION_FACTOR",
              {.description =
                      "Bounding volume of the sphere is the sphere center extended by this factor "
                      "times the sphere radius in all directions (+ and -).",
                  .default_value = 2.0}),

          parameter<bool>("WRITE_GEOMETRIC_SEARCH_VISUALIZATION",
              {.description = "If visualization output for the geometric search should be written",
                  .default_value = false})},
      {.required = false});
  return spec;
}

FOUR_C_NAMESPACE_CLOSE