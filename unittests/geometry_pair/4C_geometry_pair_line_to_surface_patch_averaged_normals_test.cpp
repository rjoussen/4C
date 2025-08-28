// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <gtest/gtest.h>

#include "4C_comm_mpi_utils.hpp"
#include "4C_geometry_pair_element.hpp"
#include "4C_geometry_pair_element_evaluation_functions.hpp"
#include "4C_geometry_pair_element_faces.hpp"
#include "4C_geometry_pair_line_to_surface_patch_geometry_test.hpp"
#include "4C_geometry_pair_line_to_surface_patch_results_test.hpp"
#include "4C_geometry_pair_scalar_types.hpp"
#include "4C_linalg_vector.hpp"



namespace
{
  using namespace FourC;

  /**
   * \brief Class to test the surface patch functionality of the geometry pairs.
   */
  class GeometryPairLineToSurfacePatchTest : public ::testing::Test
  {
   protected:
    /**
     * \brief Set up the testing environment, is called before each test.
     */
    GeometryPairLineToSurfacePatchTest()
    {
      auto comm = MPI_COMM_WORLD;
      discret_ = std::make_shared<Core::FE::Discretization>("unit_test", comm, 3);
    }

    /**
     * \brief Return a reference to the connected faces of a face element.
     */
    template <typename A>
    std::map<int, GeometryPair::ConnectedFace>& get_connected_faces(A& face_element)
    {
      return face_element.connected_faces_;
    }

    /**
     * \brief Get the number of dofs for the beam.
     */
    template <typename A>
    unsigned int get_n_other_dof(A& face_element)
    {
      return face_element.n_dof_other_element_;
    }

    //! Pointer to the discretization object that holds the geometry for the tests.
    std::shared_ptr<Core::FE::Discretization> discret_;
  };

  /**
   * \brief Test the evaluation of averaged normals on a patch of hex8/quad4 elements.
   */
  TEST_F(GeometryPairLineToSurfacePatchTest, TestSurfacePatchAveragedNormalsQuad4)
  {
    using namespace FourC;

    // Define the type of the face elements.
    using surface = GeometryPair::t_quad4;
    using scalar_type = GeometryPair::line_to_surface_patch_scalar_type;
    using face_element_type = GeometryPair::FaceElementPatchTemplate<surface, scalar_type>;

    // Tolerance for the result tests.
    const double eps = 1e-12;

    // Fill the discretization object with the geometry.
    std::unordered_map<int, std::shared_ptr<GeometryPair::FaceElement>> face_elements_map;
    xtest_surface_patch_quad4<face_element_type>(*discret_, face_elements_map);

    // Load the result vectors.
    std::vector<double> reference_normals, current_normals, position;
    std::vector<std::vector<double>> current_normals_derivative, position_derivative;
    std::vector<std::vector<std::vector<double>>> current_normals_derivative_2,
        position_derivative_2;
    xtest_surface_patch_quad4_results(reference_normals, current_normals,
        current_normals_derivative, current_normals_derivative_2, position, position_derivative,
        position_derivative_2);

    // Face element that will be analyzed.
    const unsigned int investigated_face_element_volume_id = 14;
    std::shared_ptr<face_element_type> face_element = std::dynamic_pointer_cast<face_element_type>(
        face_elements_map[investigated_face_element_volume_id]);

    // Offset in the derivatives for the beam dof.
    const unsigned int dof_offset = get_n_other_dof(*face_element);

    // Setup all face elements and get the patch information.
    for (auto& face_element_map_iterator : face_elements_map)
      face_element_map_iterator.second->setup(discret_, face_elements_map);

    {
      // Check if the GID are correct.
      std::vector<int> patch_dof_gid_reference = {126, 127, 128, 111, 112, 113, 117, 118, 119, 129,
          130, 131, 120, 121, 122, 102, 103, 104, 99, 100, 101, 108, 109, 110, 114, 115, 116};
      EXPECT_EQ(face_element->get_patch_gid().size(), patch_dof_gid_reference.size());
      for (unsigned int i = 0; i < face_element->get_patch_gid().size(); i++)
        EXPECT_EQ(face_element->get_patch_gid()[i], patch_dof_gid_reference[i]);

      // Check if the local node ID map of the connected faces to the main face could be found.
      EXPECT_EQ(get_connected_faces(*face_element).size(), 3);

      EXPECT_EQ(get_connected_faces(*face_element)[10].node_lid_map_.size(), 1);
      EXPECT_EQ(get_connected_faces(*face_element)[10].node_lid_map_[3], 1);
      EXPECT_EQ(get_connected_faces(*face_element)[10].my_node_patch_lid_.size(), 4);
      EXPECT_EQ(get_connected_faces(*face_element)[10].my_node_patch_lid_[0], 5);
      EXPECT_EQ(get_connected_faces(*face_element)[10].my_node_patch_lid_[1], 6);
      EXPECT_EQ(get_connected_faces(*face_element)[10].my_node_patch_lid_[2], 7);
      EXPECT_EQ(get_connected_faces(*face_element)[10].my_node_patch_lid_[3], 1);

      EXPECT_EQ(get_connected_faces(*face_element)[11].node_lid_map_.size(), 2);
      EXPECT_EQ(get_connected_faces(*face_element)[11].node_lid_map_[0], 1);
      EXPECT_EQ(get_connected_faces(*face_element)[11].node_lid_map_[3], 2);
      EXPECT_EQ(get_connected_faces(*face_element)[11].my_node_patch_lid_.size(), 4);
      EXPECT_EQ(get_connected_faces(*face_element)[11].my_node_patch_lid_[0], 1);
      EXPECT_EQ(get_connected_faces(*face_element)[11].my_node_patch_lid_[1], 7);
      EXPECT_EQ(get_connected_faces(*face_element)[11].my_node_patch_lid_[2], 8);
      EXPECT_EQ(get_connected_faces(*face_element)[11].my_node_patch_lid_[3], 2);

      EXPECT_EQ(get_connected_faces(*face_element)[13].node_lid_map_.size(), 2);
      EXPECT_EQ(get_connected_faces(*face_element)[13].node_lid_map_[2], 1);
      EXPECT_EQ(get_connected_faces(*face_element)[13].node_lid_map_[3], 0);
      EXPECT_EQ(get_connected_faces(*face_element)[13].my_node_patch_lid_.size(), 4);
      EXPECT_EQ(get_connected_faces(*face_element)[13].my_node_patch_lid_[0], 4);
      EXPECT_EQ(get_connected_faces(*face_element)[13].my_node_patch_lid_[1], 5);
      EXPECT_EQ(get_connected_faces(*face_element)[13].my_node_patch_lid_[2], 1);
      EXPECT_EQ(get_connected_faces(*face_element)[13].my_node_patch_lid_[3], 0);
    }

    // Calculate the averaged reference normals on the face.
    face_element->calculate_averaged_reference_normals(face_elements_map);
    {
      for (unsigned int i = 0; i < reference_normals.size(); i++)
        EXPECT_NEAR(face_element->get_face_reference_element_data().nodal_normals_(i),
            reference_normals[i], eps);
    }

    // Set the state in the face element, here also the FAD variables for each patch are set.
    Core::LinAlg::Map gid_map(discret_->num_global_nodes() * 3, discret_->num_global_nodes() * 3, 0,
        discret_->get_comm());
    auto displacement_vector = std::make_shared<Core::LinAlg::Vector<double>>(gid_map);
    for (int i = 0; i < displacement_vector->global_length(); i++)
      (*displacement_vector).get_values()[i] = i * 0.01;
    face_element->set_state(displacement_vector, face_elements_map);
    {
      // Check the values of the averaged normals.
      for (unsigned int i_dof = 0; i_dof < 3 * surface::n_nodes_; i_dof++)
      {
        EXPECT_NEAR(Core::FADUtils::cast_to_double(
                        face_element->get_face_element_data().nodal_normals_(i_dof)),
            current_normals[i_dof], eps);
        for (unsigned int i_der = 0; i_der < face_element->get_patch_gid().size(); i_der++)
        {
          EXPECT_NEAR(Core::FADUtils::cast_to_double(
                          face_element->get_face_element_data().nodal_normals_(i_dof).dx(
                              dof_offset + i_der)),
              current_normals_derivative[i_dof][i_der], eps);
          for (unsigned int i_der_2 = 0; i_der_2 < face_element->get_patch_gid().size(); i_der_2++)
          {
            EXPECT_NEAR(Core::FADUtils::cast_to_double(face_element->get_face_element_data()
                                .nodal_normals_(i_dof)
                                .dx(dof_offset + i_der)
                                .dx(dof_offset + i_der_2)),
                current_normals_derivative_2[i_dof][i_der][i_der_2], eps);
          }
        }
      }

      // Check an surface position on the element.
      Core::LinAlg::Matrix<3, 1, double> xi;
      xi(0) = 0.2;
      xi(1) = -0.8;
      xi(2) = 0.69;
      Core::LinAlg::Matrix<3, 1, scalar_type> r;
      GeometryPair::evaluate_surface_position<surface>(
          xi, face_element->get_face_element_data(), r);
      for (unsigned int i_dim = 0; i_dim < 3; i_dim++)
      {
        EXPECT_NEAR(Core::FADUtils::cast_to_double(r(i_dim)), position[i_dim], eps);
        for (unsigned int i_der = 0; i_der < face_element->get_patch_gid().size(); i_der++)
        {
          EXPECT_NEAR(Core::FADUtils::cast_to_double(r(i_dim).dx(dof_offset + i_der)),
              position_derivative[i_dim][i_der], eps);
          for (unsigned int i_der_2 = 0; i_der_2 < face_element->get_patch_gid().size(); i_der_2++)
            EXPECT_NEAR(Core::FADUtils::cast_to_double(
                            r(i_dim).dx(dof_offset + i_der).dx(dof_offset + i_der_2)),
                position_derivative_2[i_dim][i_der][i_der_2], eps);
        }
      }
    }
  }

}  // namespace
