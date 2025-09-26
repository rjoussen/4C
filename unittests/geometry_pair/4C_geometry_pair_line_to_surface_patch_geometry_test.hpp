// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_GEOMETRY_PAIR_LINE_TO_SURFACE_PATCH_GEOMETRY_TEST_HPP
#define FOUR_C_GEOMETRY_PAIR_LINE_TO_SURFACE_PATCH_GEOMETRY_TEST_HPP


#include "4C_fem_discretization.hpp"
#include "4C_geometry_pair_element_faces.hpp"
#include "4C_solid_3D_ele.hpp"
#include "4C_solid_3D_ele_surface.hpp"

#include <memory>
#include <unordered_map>
#include <vector>


namespace
{
  using namespace FourC;

  /**
   * \brief Fill a discretization with the geometry for the unit test.
   */
  template <typename FaceElementType>
  void xtest_surface_patch_quad4(Core::FE::Discretization& discret,
      std::unordered_map<int, std::shared_ptr<GeometryPair::FaceElement>>& face_elements_map)
  {
    using namespace FourC;

    // The order here is the following: First all x-values, then y-values and then z-values.
    std::vector<double> points = {1, 0.914608710081328, 0.613445400072412, 0.666666666666667, 1,
        0.914988288005194, 0.608832472629939, 0.666666666666667, 0.915722552019025,
        0.617885698155344, 0.915328376050839, 0.608008652205731, 1, 0.677669199032765, 1,
        0.663714344169117, 0.305406479825865, 0.333333333333333, 0.297268211988688,
        0.333333333333333, 0.30631403457946, 0.28954114121979, 0.334933600787086, 0.309727499398062,
        0, 3.33066907387547e-16, 0, 2.77555756156289e-16, 0, 0, 0, 0, 1, 0.915372004085148,
        0.603654517990169, 0.666666666666667, 0.914938458270761, 0.597131389631667, 1,
        0.649808694939741, 0.289012515582511, 0.333333333333333, 0.272900260165232,
        0.284835646644582, 0, 0, 0, 0, -0.1, 0.297300857123347, 0.274412402416529,
        -0.0666666666666667, -0.1, 0.294582283217281, 0.270552077872148, -0.0666666666666667,
        0.702426303321499, 0.62077477585621, 0.699629322271778, 0.615075928451673, 1.1,
        0.96203472565012, 1.1, 0.957096003499481, 0.277315184045015, -0.0333333333333334,
        0.277363130592267, -0.0333333333333334, 0.590616526107514, 0.59142030007428,
        0.901354836931972, 0.903848072185622, 0.333333333333333, -4.16333634234434e-17,
        0.333333333333333, -3.46944695195361e-17, 0.666666666666667, 0.666666666666667, 1, 1, -0.1,
        0.291863980063899, 0.265748734619595, -0.0666666666666667, 0.696832619828237,
        0.608208150759378, 1.1, 0.952313306202292, 0.277628539006058, -0.0333333333333333,
        0.592477798188276, 0.907337750184763, 0.333333333333333, 0, 0.666666666666667, 1, -0.3,
        -0.3, -0.3, -0.3, -0.2, -0.131830286364796, -0.112029085619941, -0.183333333333333, -0.3,
        -0.3, -0.0728034375083871, -0.0466025463193502, -0.3, -0.3, -0.025, 0.0110181722478759,
        -0.3, -0.3, -0.11555759043839, -0.166666666666667, -0.3, -0.0676543302987808, -0.3,
        -0.0204807213222538, -0.3, -0.3, -0.15, -0.15, -0.3, -0.15, -0.3, -0.15, -0.1,
        0.0354506598306837, 0.0759203836218857, -0.0666666666666667, 0.153674899548879,
        0.207259346030087, 0.25, 0.323405835580196, 0.0668194593876439, -0.0333333333333333,
        0.159712762721663, 0.244274360917852, 0, 0, 0, 0};

    // Displacements are ordered classically: node1_x, node1_y, node1_z, node2_x, ...
    std::vector<double> displacements(points.size());
    for (unsigned int i = 0; i < displacements.size(); i++) displacements[i] = 0.01 * i;

    // Connectivity of the volume elements.
    std::vector<int> connectivity_volumes = {1, 2, 3, 4, 5, 6, 7, 8, 2, 9, 10, 3, 6, 11, 12, 7, 9,
        13, 14, 10, 11, 15, 16, 12, 4, 3, 17, 18, 8, 7, 19, 20, 3, 10, 21, 17, 7, 12, 22, 19, 10,
        14, 23, 21, 12, 16, 24, 22, 18, 17, 25, 26, 20, 19, 27, 28, 17, 21, 29, 25, 19, 22, 30, 27,
        21, 23, 31, 29, 22, 24, 32, 30, 5, 6, 7, 8, 33, 34, 35, 36, 6, 11, 12, 7, 34, 37, 38, 35,
        11, 15, 16, 12, 37, 39, 40, 38, 8, 7, 19, 20, 36, 35, 41, 42, 7, 12, 22, 19, 35, 38, 43, 41,
        12, 16, 24, 22, 38, 40, 44, 43, 20, 19, 27, 28, 42, 41, 45, 46, 19, 22, 30, 27, 41, 43, 47,
        45, 22, 24, 32, 30, 43, 44, 48, 47};

    // Connectivity of the face elements.
    std::vector<int> connectivity_faces = {42, 36, 35, 41, 36, 33, 34, 35, 41, 35, 38, 43, 35, 34,
        37, 38, 43, 38, 40, 44, 38, 37, 39, 40};

    // Define some shortcuts.
    const unsigned int n_nodes_volume = 8;
    const unsigned int n_nodes_face = 4;
    const unsigned int n_nodes_problem = points.size() / 3;
    const unsigned int n_el_volume = connectivity_volumes.size() / n_nodes_volume;
    const unsigned int n_el_face = connectivity_faces.size() / n_nodes_face;

    // Add all nodes to the discretization.
    {
      std::array<double, 3> coordinates;
      std::shared_ptr<Core::Nodes::Node> new_node;
      for (unsigned int i_node = 0; i_node < n_nodes_problem; i_node++)
      {
        for (unsigned int i_dim = 0; i_dim < 3; i_dim++)
          coordinates[i_dim] = points[i_node + n_nodes_problem * i_dim];

        discret.add_node(coordinates, i_node, nullptr);
      }
    }

    // Add all volume elements to the discretization.
    {
      std::vector<int> node_ids(n_nodes_volume);
      std::shared_ptr<Core::Elements::Element> new_element;
      for (unsigned int i_el = 0; i_el < n_el_volume; i_el++)
      {
        for (unsigned int i_node = 0; i_node < n_nodes_volume; i_node++)
          node_ids[i_node] = connectivity_volumes[i_node + n_nodes_volume * i_el] - 1;
        new_element =
            std::shared_ptr<Core::Elements::Element>(new Discret::Elements::Solid(i_el, 0));
        new_element->set_node_ids(n_nodes_volume, node_ids.data());
        discret.add_element(new_element);
      }
      discret.fill_complete({
          .assign_degrees_of_freedom = true,
          .init_elements = true,
          .do_boundary_conditions = false,
      });
    }

    // Create the face elements.
    {
      std::vector<int> node_ids(n_nodes_face);
      std::vector<Core::Nodes::Node*> element_nodes(n_nodes_face);
      int node_index;
      for (unsigned int i_el = 0; i_el < n_el_face; i_el++)
      {
        for (unsigned int i_node = 0; i_node < n_nodes_face; i_node++)
        {
          node_index = connectivity_faces[i_node + n_nodes_face * i_el] - 1;
          node_ids[i_node] = node_index;
          element_nodes[i_node] = discret.g_node(node_index);
        }

        // Find the volume element for this face element.
        int parent_id = -1;
        for (unsigned int i_el_volume = 0; i_el_volume < n_el_volume; i_el_volume++)
        {
          int n_common_nodes = 0;
          for (unsigned int i_node_volume = 0; i_node_volume < n_nodes_volume; i_node_volume++)
          {
            int volume_node_id = discret.g_element(i_el_volume)->node_ids()[i_node_volume];
            if (std::find(node_ids.begin(), node_ids.end(), volume_node_id) != node_ids.end())
              n_common_nodes += 1;
          }
          if (n_common_nodes == n_nodes_face)
          {
            parent_id = i_el_volume;
            break;
          }
        }

        // Create the Core::Elements::FaceElement.
        auto face_element = std::make_shared<Discret::Elements::SolidSurface>(i_el, 0, n_nodes_face,
            node_ids.data(), element_nodes.data(), discret.g_element(parent_id), 0);

        // Create the geometry pair face element.
        face_elements_map[parent_id] = std::make_shared<FaceElementType>(face_element, true);
      }
    }
  }
}  // namespace



#endif
