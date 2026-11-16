// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <gtest/gtest.h>

#include "4C_fem_geometry_point_location.hpp"

#include "4C_fem_discretization.hpp"
#include "4C_fem_general_element.hpp"
#include "4C_fem_general_node.hpp"
#include "4C_fem_general_utils_fem_shapefunctions.hpp"
#include "4C_io_pstream.hpp"
#include "4C_unittest_utils_create_discretization_helper_test.hpp"

namespace
{
  using namespace FourC;

  class PointLocationTest : public testing::Test
  {
   public:
    PointLocationTest()
    {
      comm_ = MPI_COMM_WORLD;
      discretization_ = std::make_shared<Core::FE::Discretization>("dummy", comm_, 3);
      TESTING::fill_discretization_hyper_cube(*discretization_, 4, comm_);
      Core::IO::cout.setup(false, false, false, Core::IO::standard, comm_, 0, 0, "dummyFilePrefix");
      discretization_->fill_complete(Core::FE::OptionsFillComplete::none());
    }

    void TearDown() override { Core::IO::cout.close(); }

    //! sheared configuration: x = X + 0.3 * Z * e_x
    static Core::LinAlg::Matrix<3, 1> sheared_position(const Core::Nodes::Node& node)
    {
      Core::LinAlg::Matrix<3, 1> x;
      for (int d = 0; d < 3; ++d) x(d) = node.x()[d];
      x(0) += 0.3 * node.x()[2];
      return x;
    }

    //! map parameter space coordinates of an element to the given configuration
    static Core::LinAlg::Matrix<3, 1> map_to_configuration(const Core::Elements::Element& element,
        const Core::LinAlg::Matrix<3, 1>& xi, const Core::Geo::NodePositionFunction& position)
    {
      Core::LinAlg::Matrix<8, 1> shape_functions;
      Core::FE::shape_function<Core::FE::CellType::hex8>(xi, shape_functions);
      Core::LinAlg::Matrix<3, 1> x(Core::LinAlg::Initialization::zero);
      for (int i = 0; i < 8; ++i) x.update(shape_functions(i), position(*element.nodes()[i]), 1.0);
      return x;
    }

    static Core::LinAlg::Matrix<3, 1> make_point(double x, double y, double z)
    {
      Core::LinAlg::Matrix<3, 1> point;
      point(0) = x;
      point(1) = y;
      point(2) = z;
      return point;
    }

    std::vector<const Core::Elements::Element*> column_elements() const
    {
      std::vector<const Core::Elements::Element*> elements;
      for (int lid = 0; lid < discretization_->num_my_col_elements(); ++lid)
        elements.push_back(discretization_->l_col_element(lid));
      return elements;
    }

    //! locate a single point in all column elements
    std::optional<Core::Geo::PointInElement> locate(
        const Core::LinAlg::Matrix<3, 1>& point, const double tolerance) const
    {
      return Core::Geo::locate_points_in_elements(
          column_elements(), std::span(&point, 1), sheared_position, tolerance)[0];
    }

   protected:
    std::shared_ptr<Core::FE::Discretization> discretization_;
    MPI_Comm comm_;
  };

  TEST_F(PointLocationTest, FindsPointsInEveryElement)
  {
    if (Core::Communication::num_mpi_ranks(comm_) > 1) GTEST_SKIP();

    // a point in every element of the sheared mesh is found and mapped back exactly
    const Core::LinAlg::Matrix<3, 1> xi_target = make_point(0.2, -0.7, 0.4);
    const std::vector<const Core::Elements::Element*> elements = column_elements();
    std::vector<Core::LinAlg::Matrix<3, 1>> points;
    for (const auto* element : elements)
      points.push_back(map_to_configuration(*element, xi_target, sheared_position));

    const auto locations =
        Core::Geo::locate_points_in_elements(elements, points, sheared_position);

    ASSERT_EQ(locations.size(), elements.size());
    for (std::size_t p = 0; p < points.size(); ++p)
    {
      ASSERT_TRUE(locations[p].has_value());
      EXPECT_EQ(locations[p]->element->id(), elements[p]->id());
      for (int d = 0; d < 3; ++d) EXPECT_NEAR(locations[p]->xi(d), xi_target(d), 1.0e-10);
    }
  }

  TEST_F(PointLocationTest, DoesNotFindPointsOutsideOfTheMesh)
  {
    if (Core::Communication::num_mpi_ranks(comm_) > 1) GTEST_SKIP();

    const Core::Elements::Element& element = *discretization_->g_element(0);
    const auto point = map_to_configuration(element, make_point(-3.0, 0.0, 0.0), sheared_position);

    EXPECT_FALSE(locate(point, 1.0e-6).has_value());
  }

  TEST_F(PointLocationTest, AcceptsPointsSlightlyOutsideWithinTolerance)
  {
    if (Core::Communication::num_mpi_ranks(comm_) > 1) GTEST_SKIP();

    const Core::Elements::Element& element = *discretization_->g_element(0);
    const Core::LinAlg::Matrix<3, 1> xi_outside = make_point(-1.0 - 1.0e-9, 0.0, 0.0);
    const auto point = map_to_configuration(element, xi_outside, sheared_position);

    const auto location = locate(point, 1.0e-8);
    ASSERT_TRUE(location.has_value());
    EXPECT_EQ(location->element->id(), 0);

    EXPECT_FALSE(locate(point, 1.0e-12).has_value());
  }

  TEST_F(PointLocationTest, AcceptsPointsOutsideOfTheBoundingBoxWithinTolerance)
  {
    if (Core::Communication::num_mpi_ranks(comm_) > 1) GTEST_SKIP();

    // a point outside of the mesh and of the bounding box of the element (the mesh is not sheared
    // in y)
    const Core::Elements::Element& element = *discretization_->g_element(0);
    const Core::LinAlg::Matrix<3, 1> xi_outside = make_point(0.0, -1.2, 0.0);
    const auto point = map_to_configuration(element, xi_outside, sheared_position);

    const auto location = locate(point, 0.25);
    ASSERT_TRUE(location.has_value());
    EXPECT_EQ(location->element->id(), 0);
    EXPECT_NEAR(location->xi(1), -1.2, 1.0e-10);

    EXPECT_FALSE(locate(point, 0.1).has_value());
  }
}  // namespace
