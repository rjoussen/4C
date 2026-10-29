// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <gtest/gtest.h>

#include "4C_solid_scatra_ele_calc.hpp"

#include "4C_fem_discretization.hpp"
#include "4C_fem_dofset_predefineddofnumber.hpp"
#include "4C_fem_general_element.hpp"
#include "4C_linalg_serialdensematrix.hpp"
#include "4C_linalg_serialdensevector.hpp"
#include "4C_linalg_vector.hpp"
#include "4C_mat_material_factory.hpp"
#include "4C_mat_so3_material.hpp"
#include "4C_mat_thermostvenantkirchhoff.hpp"
#include "4C_solid_ele_calc_displacement_based.hpp"
#include "4C_solid_ele_calc_displacement_based_linear_kinematics.hpp"
#include "4C_solid_ele_calc_fbar.hpp"
#include "4C_solid_scatra_ele.hpp"
#include "4C_unittest_utils_assertions_test.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <memory>

namespace
{
  using namespace FourC;

  //! thermo-elastic material with the temperature-dependent Young's modulus E(T) = 2.1e5 - 50 T
  std::shared_ptr<Mat::ThermoStVenantKirchhoff> create_material()
  {
    Core::IO::InputParameterContainer material_data;
    material_data.add("YOUNG", std::vector<double>{2.1e5, -50.0});
    material_data.add("NUE", 0.3);
    material_data.add("DENS", 1.0);
    material_data.add("THEXPANS", 1.0e-3);
    material_data.add("INITTEMP", 290.0);
    auto parameters = std::dynamic_pointer_cast<Mat::PAR::ThermoStVenantKirchhoff>(std::shared_ptr(
        Mat::make_parameter(1, Core::Materials::MaterialType::m_thermostvenant, material_data)));

    return {new Mat::ThermoStVenantKirchhoff(parameters.get()),
        [parameters](Mat::ThermoStVenantKirchhoff* material) { delete material; }};
  }

  //! nodal data of a slightly distorted hex8 element under finite deformation with varying
  //! velocities and temperatures
  struct Hex8Data
  {
    static constexpr Core::FE::CellType celltype = Core::FE::CellType::hex8;
    static constexpr std::array<std::array<double, 3>, 8> reference_coordinates = {{
        {-1.02, -0.97, -1.04},
        {1.03, -1.01, -0.96},
        {0.98, 1.04, -1.02},
        {-0.96, 0.99, -0.98},
        {-1.01, -1.03, 1.02},
        {0.97, -0.98, 1.03},
        {1.04, 1.02, 0.97},
        {-0.99, 0.96, 1.01},
    }};
    // finite deformation, such that the volume change varies within the element
    static constexpr std::array<std::array<double, 3>, 8> displacements = {{
        {0.03, -0.02, 0.04},
        {-0.01, 0.05, -0.03},
        {0.04, 0.01, 0.02},
        {-0.05, -0.04, 0.01},
        {0.02, 0.03, -0.05},
        {-0.03, 0.02, 0.04},
        {0.01, -0.05, -0.02},
        {0.05, 0.04, 0.03},
    }};
    static constexpr std::array<std::array<double, 3>, 8> velocities = {{
        {0.2, -0.1, 0.3},
        {-0.4, 0.2, 0.1},
        {0.1, 0.3, -0.2},
        {0.3, -0.2, 0.4},
        {-0.2, 0.4, 0.1},
        {0.4, 0.1, -0.3},
        {-0.1, -0.3, 0.2},
        {0.2, 0.2, -0.4},
    }};
    static constexpr std::array<double, 8> temperatures = {
        305.0, 296.0, 309.0, 292.0, 301.0, 307.0, 294.0, 299.0};
  };

  //! nodal data of a slightly distorted quad4 element, the in-plane part of the bottom face of the
  //! hex8 element
  struct Quad4Data
  {
    static constexpr Core::FE::CellType celltype = Core::FE::CellType::quad4;
    static constexpr std::array<std::array<double, 2>, 4> reference_coordinates = {{
        {-1.02, -0.97},
        {1.03, -1.01},
        {0.98, 1.04},
        {-0.96, 0.99},
    }};
    static constexpr std::array<std::array<double, 2>, 4> displacements = {{
        {0.03, -0.02},
        {-0.01, 0.05},
        {0.04, 0.01},
        {-0.05, -0.04},
    }};
    static constexpr std::array<std::array<double, 2>, 4> velocities = {{
        {0.2, -0.1},
        {-0.4, 0.2},
        {0.1, 0.3},
        {0.3, -0.2},
    }};
    static constexpr std::array<double, 4> temperatures = {305.0, 296.0, 309.0, 292.0};
  };

  /*!
   * @brief A single, slightly distorted solid-scatra element with a thermo-elastic material
   *
   * The element is part of a discretization with the displacements and velocities in the first and
   * the temperatures in the second dofset, i.e., the setup of the structural discretization in TSI.
   */
  template <typename ElementData>
  class SolidScatraElementTest : public ::testing::Test
  {
   protected:
    using Data = ElementData;
    static constexpr Core::FE::CellType celltype = Data::celltype;
    static constexpr int num_nodes = Core::FE::num_nodes(celltype);
    static constexpr int num_dim = Core::FE::dim<celltype>;
    static constexpr int num_dof = num_nodes * num_dim;

    using NodalVectors = std::array<std::array<double, num_dim>, num_nodes>;
    using NodalScalars = std::array<double, num_nodes>;

    std::shared_ptr<Mat::So3Material> material_ = create_material();

    // empty solid discretization
    Core::FE::Discretization discretization_{"structure", MPI_COMM_WORLD, num_dim};

    // the solid-scatra element
    std::shared_ptr<Discret::Elements::SolidScatra<num_dim>> element_ =
        std::make_shared<Discret::Elements::SolidScatra<num_dim>>(0, 0);
    Core::Elements::LocationArray location_array_{2};

    // set the nodal displacements, velocities and temperatures in the discretization
    void set_states(const NodalVectors& displacements, const NodalVectors& velocities,
        const NodalScalars& temperatures)
    {
      const auto set_nodal_vectors = [&](const std::string& name, const NodalVectors& values)
      {
        Core::LinAlg::Vector<double> state(*discretization_.dof_row_map(0), true);
        for (int dof = 0; dof < num_dof; ++dof)
        {
          // sort the nodal vectors into the state vector
          state.replace_global_value(
              location_array_[0].lm_[dof], values[dof / num_dim][dof % num_dim]);
        }
        discretization_.set_state(0, name, state);
      };
      set_nodal_vectors("displacement", displacements);
      set_nodal_vectors("velocity", velocities);

      Core::LinAlg::Vector<double> temperature_state(*discretization_.dof_row_map(1), true);
      for (int node = 0; node < num_nodes; ++node)
        temperature_state.replace_global_value(location_array_[1].lm_[node], temperatures[node]);
      discretization_.set_state(1, "temperature", temperature_state);
    }

    void SetUp() override
    {
      // add the element with its distorted reference configuration to the discretization
      std::array<int, num_nodes> node_ids{};
      for (int node = 0; node < num_nodes; ++node)
      {
        discretization_.add_node(Data::reference_coordinates[node], node, nullptr);
        node_ids[node] = node;
      }
      element_->set_node_ids(num_nodes, node_ids.data());
      discretization_.add_element(element_);

      // add the dofsets for the displacements and temperatures to the discretization and fill it
      discretization_.add_dof_set(
          std::make_shared<Core::DOFSets::DofSetPredefinedDoFNumber>(1, 0, 0, true));
      discretization_.fill_complete();
      element_->location_vector(discretization_, location_array_);
    }
  };

  template <typename ElementData, template <Core::FE::CellType> typename Formulation>
  struct TestCase
  {
    using Data = ElementData;
    using SolidFormulation = Formulation<ElementData::celltype>;
  };

  //! evaluation of the element with the given solid formulation, 2D elements in plane strain
  template <typename TestCase>
  class SolidScatraEleCalcTest : public SolidScatraElementTest<typename TestCase::Data>
  {
   protected:
    using Base = SolidScatraElementTest<typename TestCase::Data>;
    using typename Base::Data;
    using typename Base::NodalScalars;
    using typename Base::NodalVectors;
    using ElementCalc =
        Discret::Elements::SolidScatraEleCalc<Base::celltype, typename TestCase::SolidFormulation>;

    //! derivative of the velocities w.r.t. the displacements due to the time integration
    static constexpr double d_velocity_d_displacement = 10.0;

    struct HeatSource
    {
      Core::LinAlg::SerialDenseVector value{Base::num_nodes};
      Core::LinAlg::SerialDenseMatrix d_temperatures{Base::num_nodes, Base::num_nodes};
      Core::LinAlg::SerialDenseMatrix d_displacements{Base::num_nodes, Base::num_dof};
    };

    const double total_time_ = 1.0;
    const double time_step_size_ = 0.1;
    const Mat::EvaluationContext<Base::num_dim> context_{.total_time = &total_time_,
        .time_step_size = &time_step_size_,
        .xi = nullptr,
        .ref_coords = nullptr};

    //! calculation class of the element, 2D elements with the given reference thickness
    static ElementCalc create_element_calc(const double reference_thickness)
    {
      if constexpr (Base::num_dim == 3)
        return ElementCalc{};
      else
        return ElementCalc{reference_thickness, Discret::Elements::PlaneAssumption::plane_strain};
    }

    // the calculation class for the solid-scatra element with the given solid formulation
    ElementCalc element_calc_ = create_element_calc(2.0);

    void SetUp() override
    {
      Base::SetUp();
      element_calc_.setup(*this->material_, Core::IO::InputParameterContainer{});
    }

    //! internal force vector of the element for the given nodal displacements and temperatures
    Core::LinAlg::SerialDenseVector internal_force(
        const NodalVectors& displacements, const NodalScalars& temperatures)
    {
      this->set_states(displacements, Data::velocities, temperatures);
      Teuchos::ParameterList params;
      Core::LinAlg::SerialDenseVector force(Base::num_dof);
      element_calc_.evaluate_nonlinear_force_stiffness_mass(*this->element_, *this->material_,
          this->discretization_, this->location_array_, params, &force, nullptr, nullptr);
      return force;
    }

    //! analytic linearization of the internal force vector wrt the nodal temperatures
    Core::LinAlg::SerialDenseMatrix d_internal_force_d_temperatures(
        const NodalVectors& displacements, const NodalScalars& temperatures)
    {
      this->set_states(displacements, Data::velocities, temperatures);
      Teuchos::ParameterList params;
      Core::LinAlg::SerialDenseMatrix d_force_d_temperatures(Base::num_dof, Base::num_nodes);
      element_calc_.evaluate_d_stress_d_scalar(*this->element_, *this->material_,
          this->discretization_, this->location_array_, params, d_force_d_temperatures);
      return d_force_d_temperatures;
    }

    //! mechanical heat source and its derivatives for the given nodal states
    HeatSource heat_source(ElementCalc& element_calc, const NodalVectors& displacements,
        const NodalVectors& velocities, const NodalScalars& temperatures)
    {
      this->set_states(displacements, velocities, temperatures);
      HeatSource heat_source;
      element_calc.evaluate_mechanical_heat_source(*this->element_, *this->material_,
          this->discretization_, this->location_array_, context_, d_velocity_d_displacement,
          heat_source.value, &heat_source.d_temperatures, &heat_source.d_displacements);
      return heat_source;
    }

    HeatSource heat_source(const NodalVectors& displacements, const NodalVectors& velocities,
        const NodalScalars& temperatures)
    {
      return heat_source(element_calc_, displacements, velocities, temperatures);
    }

    static double max_abs(const Core::LinAlg::SerialDenseMatrix& matrix)
    {
      double result = 0.0;
      for (int i = 0; i < matrix.numRows(); ++i)
        for (int j = 0; j < matrix.numCols(); ++j)
          result = std::max(result, std::abs(matrix(i, j)));
      return result;
    }
  };

  using TestCases =
      ::testing::Types<TestCase<Hex8Data, Discret::Elements::DisplacementBasedFormulation>,
          TestCase<Hex8Data, Discret::Elements::DisplacementBasedLinearKinematicsFormulation>,
          TestCase<Hex8Data, Discret::Elements::FBarFormulation>,
          TestCase<Quad4Data, Discret::Elements::DisplacementBasedFormulation>,
          TestCase<Quad4Data, Discret::Elements::DisplacementBasedLinearKinematicsFormulation>>;
  TYPED_TEST_SUITE(SolidScatraEleCalcTest, TestCases);

  TYPED_TEST(SolidScatraEleCalcTest, InternalForceWrtTemperatureMatchesFiniteDifferences)
  {
    using Data = typename TestFixture::Data;
    constexpr int num_nodes = TestFixture::num_nodes;
    constexpr int num_dof = TestFixture::num_dof;

    const auto d_force_d_temperatures =
        this->d_internal_force_d_temperatures(Data::displacements, Data::temperatures);

    // Here, The internal force is quadratic in the temperature since the Young's modulus is linear
    // in it. Thus, the central difference is exact and the tolerance only accounts for round-off.
    const double perturbation_size = 1.0e-4;
    Core::LinAlg::SerialDenseMatrix d_force_d_temperatures_FD(num_dof, num_nodes);
    for (int node = 0; node < num_nodes; ++node)
    {
      auto temperatures_plus = Data::temperatures;
      auto temperatures_minus = Data::temperatures;
      temperatures_plus[node] += perturbation_size;
      temperatures_minus[node] -= perturbation_size;
      const auto force_plus = this->internal_force(Data::displacements, temperatures_plus);
      const auto force_minus = this->internal_force(Data::displacements, temperatures_minus);
      for (int dof = 0; dof < num_dof; ++dof)
        d_force_d_temperatures_FD(dof, node) =
            (force_plus[dof] - force_minus[dof]) / (2.0 * perturbation_size);
    }

    FOUR_C_EXPECT_NEAR(d_force_d_temperatures, d_force_d_temperatures_FD, 1.0e-6);
  }

  TYPED_TEST(SolidScatraEleCalcTest, HeatSourceWrtTemperatureMatchesFiniteDifferences)
  {
    using Data = typename TestFixture::Data;
    constexpr int num_nodes = TestFixture::num_nodes;

    const auto heat_source =
        this->heat_source(Data::displacements, Data::velocities, Data::temperatures);
    ASSERT_GT(heat_source.value.normInf(), 0.0);

    // The stress-temperature modulus is linear in the temperature, hence the heat source is
    // quadratic in it and the central difference is exact up to round-off.
    const double perturbation_size = 1.0e-3;
    Core::LinAlg::SerialDenseMatrix d_temperatures_FD(num_nodes, num_nodes);
    for (int node = 0; node < num_nodes; ++node)
    {
      auto temperatures_plus = Data::temperatures;
      auto temperatures_minus = Data::temperatures;
      temperatures_plus[node] += perturbation_size;
      temperatures_minus[node] -= perturbation_size;
      const auto plus = this->heat_source(Data::displacements, Data::velocities, temperatures_plus);
      const auto minus =
          this->heat_source(Data::displacements, Data::velocities, temperatures_minus);
      for (int row = 0; row < num_nodes; ++row)
        d_temperatures_FD(row, node) =
            (plus.value[row] - minus.value[row]) / (2.0 * perturbation_size);
    }

    FOUR_C_EXPECT_NEAR(heat_source.d_temperatures, d_temperatures_FD,
        1.0e-8 * TestFixture::max_abs(d_temperatures_FD));
  }

  TYPED_TEST(SolidScatraEleCalcTest, HeatSourceWrtDisplacementsMatchesFiniteDifferences)
  {
    using Data = typename TestFixture::Data;
    constexpr int num_nodes = TestFixture::num_nodes;
    constexpr int num_dim = TestFixture::num_dim;
    constexpr int num_dof = TestFixture::num_dof;

    const auto heat_source =
        this->heat_source(Data::displacements, Data::velocities, Data::temperatures);

    // the velocities depend on the displacements through the time integration
    const double perturbation_size = 1.0e-6;
    Core::LinAlg::SerialDenseMatrix d_displacements_FD(num_nodes, num_dof);
    for (int dof = 0; dof < num_dof; ++dof)
    {
      auto displacements_plus = Data::displacements;
      auto displacements_minus = Data::displacements;
      auto velocities_plus = Data::velocities;
      auto velocities_minus = Data::velocities;
      displacements_plus[dof / num_dim][dof % num_dim] += perturbation_size;
      displacements_minus[dof / num_dim][dof % num_dim] -= perturbation_size;
      velocities_plus[dof / num_dim][dof % num_dim] +=
          TestFixture::d_velocity_d_displacement * perturbation_size;
      velocities_minus[dof / num_dim][dof % num_dim] -=
          TestFixture::d_velocity_d_displacement * perturbation_size;
      const auto plus = this->heat_source(displacements_plus, velocities_plus, Data::temperatures);
      const auto minus =
          this->heat_source(displacements_minus, velocities_minus, Data::temperatures);
      for (int row = 0; row < num_nodes; ++row)
        d_displacements_FD(row, dof) =
            (plus.value[row] - minus.value[row]) / (2.0 * perturbation_size);
    }

    FOUR_C_EXPECT_NEAR(heat_source.d_displacements, d_displacements_FD,
        1.0e-6 * TestFixture::max_abs(d_displacements_FD));
  }

  using PlaneStrainEleCalcTest =
      SolidScatraEleCalcTest<TestCase<Quad4Data, Discret::Elements::DisplacementBasedFormulation>>;

  TEST_F(PlaneStrainEleCalcTest, HeatSourceIsIndependentOfReferenceThickness)
  {
    // the thermal elements integrate per unit thickness, hence also the heat source
    auto element_calc_unit_thickness = create_element_calc(1.0);
    element_calc_unit_thickness.setup(*material_, Core::IO::InputParameterContainer{});

    const auto heat_source =
        this->heat_source(element_calc_, Data::displacements, Data::velocities, Data::temperatures);
    const auto heat_source_unit_thickness = this->heat_source(
        element_calc_unit_thickness, Data::displacements, Data::velocities, Data::temperatures);

    FOUR_C_EXPECT_NEAR(heat_source.value, heat_source_unit_thickness.value,
        1.0e-14 * heat_source_unit_thickness.value.normInf());
  }
}  // namespace
