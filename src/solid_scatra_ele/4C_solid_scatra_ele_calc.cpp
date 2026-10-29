// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_solid_scatra_ele_calc.hpp"

#include "4C_fem_discretization.hpp"
#include "4C_fem_general_cell_type.hpp"
#include "4C_fem_general_cell_type_traits.hpp"
#include "4C_fem_general_extract_values.hpp"
#include "4C_fem_general_utils_interpolation.hpp"
#include "4C_linalg_tensor.hpp"
#include "4C_linalg_tensor_generators.hpp"
#include "4C_mat_monolithic_solid_scalar_material.hpp"
#include "4C_mat_so3_material.hpp"
#include "4C_mat_trait_thermo_solid.hpp"
#include "4C_solid_ele_calc_displacement_based.hpp"
#include "4C_solid_ele_calc_displacement_based_linear_kinematics.hpp"
#include "4C_solid_ele_calc_eas.hpp"
#include "4C_solid_ele_calc_fbar.hpp"
#include "4C_solid_ele_calc_lib.hpp"
#include "4C_solid_ele_calc_lib_formulation.hpp"
#include "4C_solid_ele_calc_lib_integration.hpp"
#include "4C_solid_ele_calc_lib_io.hpp"
#include "4C_solid_ele_calc_lib_plane.hpp"
#include "4C_solid_ele_formulation.hpp"
#include "4C_solid_ele_interface_serializable.hpp"
#include "4C_utils_exceptions.hpp"

#include <Teuchos_ParameterList.hpp>

#include <memory>
#include <optional>

FOUR_C_NAMESPACE_OPEN

namespace
{
  template <typename T>
  T* get_ptr(std::optional<T>& opt)
  {
    return opt.has_value() ? &opt.value() : nullptr;
  }
  template <typename T>
  const T* get_data(const std::optional<std::vector<T>>& opt)
  {
    return opt.has_value() ? opt.value().data() : nullptr;
  }

  template <Core::FE::CellType celltype>
  std::vector<
      Core::LinAlg::SymmetricTensor<double, Core::FE::dim<celltype>, Core::FE::dim<celltype>>>
  evaluate_d_material_stress_d_scalars(Mat::So3Material& solid_material,
      const Discret::Elements::ElementProperties<celltype>& element_properties,
      const Core::LinAlg::Tensor<double, Core::FE::dim<celltype>, Core::FE::dim<celltype>>&
          deformation_gradient,
      const Core::LinAlg::SymmetricTensor<double, Core::FE::dim<celltype>, Core::FE::dim<celltype>>&
          gl_strain,
      Teuchos::ParameterList& params,
      const Mat::EvaluationContext<Core::FE::dim<celltype>>& context, const int gp,
      const int eleGID, const int num_scalars)
  {
    auto* monolithic_material = dynamic_cast<Mat::MonolithicSolidScalarMaterial*>(&solid_material);

    FOUR_C_ASSERT_ALWAYS(
        monolithic_material, "Your material does not allow to evaluate a monolithic ssi material!");

    if constexpr (Core::FE::dim<celltype> == 3)
    {
      // The derivative of the solid stress w.r.t. the scalar is implemented in the normal
      // material Evaluate call by not passing the linearization matrix.
      return monolithic_material->evaluate_d_stress_d_scalars(
          deformation_gradient, gl_strain, params, context, num_scalars, gp, eleGID);
    }
    else
    {
      std::vector<Core::LinAlg::SymmetricTensor<double, 3, 3>> d_stress_d_scalars_3d;
      Discret::Elements::transform_to_3d(solid_material, element_properties, deformation_gradient,
          gl_strain, params, context, gp, eleGID,
          [&](const Core::LinAlg::Tensor<double, 3, 3>& defgrd_3d,
              const Core::LinAlg::SymmetricTensor<double, 3, 3>& gl_strain_3d,
              const Mat::EvaluationContext<3>& context_3d)
          {
            d_stress_d_scalars_3d = monolithic_material->evaluate_d_stress_d_scalars(
                defgrd_3d, gl_strain_3d, params, context_3d, num_scalars, gp, eleGID);
          });

      // only return the 2D part of the tensor
      std::vector<
          Core::LinAlg::SymmetricTensor<double, Core::FE::dim<celltype>, Core::FE::dim<celltype>>>
          d_stress_d_scalars_2d(num_scalars);
      for (int k = 0; k < num_scalars; ++k)
      {
        d_stress_d_scalars_2d[k] = Core::LinAlg::assume_symmetry(Core::LinAlg::Tensor<double, 2, 2>{
            {{d_stress_d_scalars_3d[k](0, 0), d_stress_d_scalars_3d[k](0, 1)},
                {d_stress_d_scalars_3d[k](1, 0), d_stress_d_scalars_3d[k](1, 1)}}});
      }
      return d_stress_d_scalars_2d;
    }
  }

  //! mechanical heat source and its derivatives w.r.t. the kinematics of the element
  template <Core::FE::CellType celltype>
  struct MechanicalHeatSource
  {
    double value = 0.0;
    double derivative_wrt_temperature = 0.0;
    //! derivative w.r.t. the (Green-Lagrange) strain
    Core::LinAlg::SymmetricTensor<double, Core::FE::dim<celltype>, Core::FE::dim<celltype>>
        derivative_wrt_strain{};
    //! derivative w.r.t. the (Green-Lagrange) strain rate
    Core::LinAlg::SymmetricTensor<double, Core::FE::dim<celltype>, Core::FE::dim<celltype>>
        derivative_wrt_strain_rate{};
  };

  /*!
   * @brief Evaluate the mechanical heat source of the material at the kinematics of the element
   *
   * 2D elements are only supported for plane strain, i.e., the out-of-plane strain and its rate
   * vanish. For plane stress, the out-of-plane strain rate would follow from the plane stress
   * condition and depend on the temperature rate.
   */
  template <Core::FE::CellType celltype>
  MechanicalHeatSource<celltype> evaluate_material_mechanical_heat_source(
      Mat::So3Material& solid_material, Mat::Trait::ThermoSolid& thermo_solid,
      const Discret::Elements::ElementProperties<celltype>& element_properties,
      const double temperature,
      const Core::LinAlg::Tensor<double, Core::FE::dim<celltype>, Core::FE::dim<celltype>>&
          deformation_gradient,
      const Core::LinAlg::SymmetricTensor<double, Core::FE::dim<celltype>, Core::FE::dim<celltype>>&
          gl_strain,
      const Core::LinAlg::SymmetricTensor<double, Core::FE::dim<celltype>, Core::FE::dim<celltype>>&
          strain_rate,
      const Mat::EvaluationContext<Core::FE::dim<celltype>>& context, const int gp,
      const int eleGID)
  {
    if constexpr (Core::FE::dim<celltype> == 3)
    {
      const Mat::HeatSource source = thermo_solid.evaluate_mechanical_heat_source(temperature,
          {.strain = gl_strain, .strain_rate = strain_rate, .defgrad = deformation_gradient},
          context, gp, eleGID);
      return {.value = source.value,
          .derivative_wrt_temperature = source.derivative_wrt_temperature,
          .derivative_wrt_strain = source.derivative_wrt_strain,
          .derivative_wrt_strain_rate = source.derivative_wrt_strain_rate};
    }
    else
    {
      FOUR_C_ASSERT_ALWAYS(
          element_properties.plane_assumption == Discret::Elements::PlaneAssumption::plane_strain,
          "The mechanical heat source of 2D elements is only implemented for plane strain.");

      // plane strain does not need any parameters for the transformation to 3D
      Teuchos::ParameterList params;
      Mat::HeatSource source_3d;
      Discret::Elements::transform_to_3d(solid_material, element_properties, deformation_gradient,
          gl_strain, params, context, gp, eleGID,
          [&](const Core::LinAlg::Tensor<double, 3, 3>& defgrd_3d,
              const Core::LinAlg::SymmetricTensor<double, 3, 3>& gl_strain_3d,
              const Mat::EvaluationContext<3>& context_3d)
          {
            source_3d = thermo_solid.evaluate_mechanical_heat_source(temperature,
                {.strain = gl_strain_3d,
                    .strain_rate = Discret::Elements::make_3d_tensor(strain_rate),
                    .defgrad = defgrd_3d},
                context_3d, gp, eleGID);
          });

      // only return the 2D part of the derivatives
      return {.value = source_3d.value,
          .derivative_wrt_temperature = source_3d.derivative_wrt_temperature,
          .derivative_wrt_strain =
              Discret::Elements::extract_2d_part(source_3d.derivative_wrt_strain),
          .derivative_wrt_strain_rate =
              Discret::Elements::extract_2d_part(source_3d.derivative_wrt_strain_rate)};
    }
  }

  template <Core::FE::CellType celltype>
  auto interpolate_quantity_to_point(
      const Discret::Elements::ShapeFunctionsAndDerivatives<celltype>& shape_functions,
      const std::vector<Core::LinAlg::Matrix<Core::FE::num_nodes(celltype), 1>>& nodal_quantities)
  {
    std::vector<double> quantities_at_gp(nodal_quantities.size(), 0.0);

    for (std::size_t k = 0; k < nodal_quantities.size(); ++k)
    {
      quantities_at_gp[k] = shape_functions.shapefunctions_.dot(nodal_quantities[k]);
    }
    return quantities_at_gp;
  }

  template <Core::FE::CellType celltype>
  auto interpolate_quantity_to_point(
      const Discret::Elements::ShapeFunctionsAndDerivatives<celltype>& shape_functions,
      const Core::LinAlg::Matrix<Core::FE::num_nodes(celltype), 1>& nodal_quantity)
  {
    return shape_functions.shapefunctions_.dot(nodal_quantity);
  }


  template <Core::FE::CellType celltype, bool is_scalar>
  auto get_element_quantities(const int num_scalars, const std::vector<double>& quantities_at_dofs)
  {
    if constexpr (is_scalar)
    {
      Core::LinAlg::Matrix<Core::FE::num_nodes(celltype), 1> nodal_quantities(
          Core::LinAlg::Initialization::zero);
      for (int i = 0; i < Core::FE::num_nodes(celltype); ++i)
        nodal_quantities(i, 0) = quantities_at_dofs.at(i);

      return nodal_quantities;
    }
    else
    {
      std::vector<Core::LinAlg::Matrix<Core::FE::num_nodes(celltype), 1>> nodal_quantities(
          num_scalars);

      for (int k = 0; k < num_scalars; ++k)
        for (int i = 0; i < Core::FE::num_nodes(celltype); ++i)
          (nodal_quantities[k])(i, 0) = quantities_at_dofs.at(num_scalars * i + k);

      return nodal_quantities;
    }
  }

  std::optional<int> detect_field_index(const Core::FE::Discretization& discretization,
      const Core::Elements::LocationArray& la, const std::string& field_name)
  {
    std::optional<int> detected_field_index = {};
    for (int field_index = 0; field_index < la.size(); ++field_index)
    {
      if (discretization.has_state(field_index, field_name))
      {
        FOUR_C_ASSERT_ALWAYS(!detected_field_index.has_value(),
            "There are multiple dofsets with the field name {} in the discretization. Found at "
            "least in dofset {} and {}.",
            field_name, *detected_field_index, field_index);

        detected_field_index = field_index;
      }
    }

    return detected_field_index;
  }

  template <Core::FE::CellType celltype, bool is_scalar>
  auto extract_my_nodal_scalars(const Core::Elements::Element& element,
      const Core::FE::Discretization& discretization, const Core::Elements::LocationArray& la,
      const std::string& field_name)
      -> std::optional<
          std::conditional_t<is_scalar, Core::LinAlg::Matrix<Core::FE::num_nodes(celltype), 1>,
              std::vector<Core::LinAlg::Matrix<Core::FE::num_nodes(celltype), 1>>>>
  {
    std::optional<int> field_index = detect_field_index(discretization, la, field_name);
    if (!field_index.has_value())
    {
      return std::nullopt;
    }

    const int num_scalars = discretization.num_dof(*field_index, element.nodes()[0]);

    FOUR_C_ASSERT(
        !is_scalar || num_scalars == 1, "numscalars must be 1 if result type is not a vector!");

    // get quantity from discretization
    std::shared_ptr<const Core::LinAlg::Vector<double>> quantities_np =
        discretization.get_state(*field_index, field_name);

    if (quantities_np == nullptr) FOUR_C_THROW("Cannot get state vector '{}' ", field_name);

    const auto my_quantities = Core::FE::extract_values(*quantities_np, la[*field_index].lm_);

    return get_element_quantities<celltype, is_scalar>(num_scalars, my_quantities);
  }

  template <Core::FE::CellType celltype>
  void prepare_scalar_in_parameter_list(Teuchos::ParameterList& params, const std::string& name,
      const Discret::Elements::ShapeFunctionsAndDerivatives<celltype>& shape_functions,
      const std::optional<Core::LinAlg::Matrix<Core::FE::num_nodes(celltype), 1>>& nodal_quantities)
  {
    if (!nodal_quantities.has_value()) return;

    auto gp_quantities = interpolate_quantity_to_point(shape_functions, *nodal_quantities);

    params.set(name, gp_quantities);
  }

  template <Core::FE::CellType celltype>
  void prepare_scalar_in_parameter_list(Teuchos::ParameterList& params, const std::string& name,
      const Discret::Elements::ShapeFunctionsAndDerivatives<celltype>& shape_functions,
      const std::optional<std::vector<Core::LinAlg::Matrix<Core::FE::num_nodes(celltype), 1>>>&
          nodal_quantities)
  {
    if (!nodal_quantities) return;

    // the value of a Teuchos::ParameterList needs to be printable. Until we get rid of the
    // parameter list here, we wrap it into a std::shared_ptr<> :(
    auto gp_quantities = std::make_shared<std::vector<double>>();
    *gp_quantities = interpolate_quantity_to_point(shape_functions, *nodal_quantities);

    params.set(name, gp_quantities);
  }



  template <Core::FE::CellType celltype, typename SolidFormulation>
  double evaluate_cauchy_n_dir_at_xi(Mat::So3Material& mat,
      Discret::Elements::ShapeFunctionsAndDerivatives<celltype> shape_functions,
      const Core::LinAlg::Tensor<double, 3>& xi,
      const Core::LinAlg::Tensor<double, Core::FE::dim<celltype>, Core::FE::dim<celltype>>&
          deformation_gradient,
      const std::vector<double>& scalars_at_xi, const Core::LinAlg::Tensor<double, 3>& n,
      const Core::LinAlg::Tensor<double, 3>& dir, int eleGID,
      const Discret::Elements::ElementFormulationDerivativeEvaluator<celltype, SolidFormulation>&
          evaluator,
      Discret::Elements::SolidScatraCauchyNDirLinearizations<3>& linearizations)
  {
    Discret::Elements::CauchyNDirLinearizationDependencies<celltype> linearization_dependencies =
        Discret::Elements::get_initialized_cauchy_n_dir_linearization_dependencies(
            evaluator, linearizations);

    Mat::EvaluationContext<3> context{.total_time = nullptr,
        .time_step_size = nullptr,
        .xi = &xi,
        .ref_coords = nullptr};  // maybe compute ref_coords?
    double cauchy_n_dir = mat.evaluate_cauchy_n_dir_and_derivatives(deformation_gradient, n, dir,
        linearizations.solid.d_cauchyndir_dn, linearizations.solid.d_cauchyndir_ddir,
        get_ptr(linearization_dependencies.d_cauchyndir_dF),
        get_ptr(linearization_dependencies.d2_cauchyndir_dF2),
        get_ptr(linearization_dependencies.d2_cauchyndir_dF_dn),
        get_ptr(linearization_dependencies.d2_cauchyndir_dF_ddir), context, eleGID,
        scalars_at_xi.data(), nullptr, nullptr, nullptr);

    // Evaluate pure solid linearizations
    Discret::Elements::evaluate_cauchy_n_dir_linearizations<celltype>(
        linearization_dependencies, linearizations.solid);

    // Evaluate ssi-linearizations
    if (linearizations.d_cauchyndir_ds)
    {
      FOUR_C_ASSERT(linearization_dependencies.d_cauchyndir_dF, "Not all tensors are computed!");
      linearizations.d_cauchyndir_ds->shape(Core::FE::num_nodes(celltype), 1);

      static Core::LinAlg::Matrix<9, 1> d_F_dc{};
      mat.evaluate_linearization_od(deformation_gradient, (scalars_at_xi)[0], d_F_dc);

      double d_cauchyndir_ds_gp = (*linearization_dependencies.d_cauchyndir_dF).dot(d_F_dc);

      Core::LinAlg::Matrix<Core::FE::num_nodes(celltype), 1>(
          linearizations.d_cauchyndir_ds->values(), true)
          .update(d_cauchyndir_ds_gp, shape_functions.shapefunctions_, 1.0);
    }
    return cauchy_n_dir;
  }
}  // namespace

template <Core::FE::CellType celltype, typename SolidFormulation>
Discret::Elements::SolidScatraEleCalc<celltype, SolidFormulation>::SolidScatraEleCalc()
  requires(Core::FE::dim<celltype> == 3)
    : stiffness_matrix_integration_(Core::FE::create_gauss_integration<celltype>(
          get_gauss_rule_stiffness_matrix<celltype>())),
      mass_matrix_integration_(
          Core::FE::create_gauss_integration<celltype>(get_gauss_rule_mass_matrix<celltype>()))
{
}

template <Core::FE::CellType celltype, typename SolidFormulation>
Discret::Elements::SolidScatraEleCalc<celltype, SolidFormulation>::SolidScatraEleCalc(
    const double reference_thickness, const Discret::Elements::PlaneAssumption plane_assumption)
  requires(Core::FE::dim<celltype> == 2)
    : stiffness_matrix_integration_(Core::FE::create_gauss_integration<celltype>(
          get_gauss_rule_stiffness_matrix<celltype>())),
      mass_matrix_integration_(
          Core::FE::create_gauss_integration<celltype>(get_gauss_rule_mass_matrix<celltype>())),
      element_properties_(
          {.reference_thickness = reference_thickness, .plane_assumption = plane_assumption})
{
}

template <Core::FE::CellType celltype, typename SolidFormulation>
void Discret::Elements::SolidScatraEleCalc<celltype, SolidFormulation>::pack(
    Core::Communication::PackBuffer& data) const
{
  Discret::Elements::pack(data, history_data_);
}

template <Core::FE::CellType celltype, typename SolidFormulation>
void Discret::Elements::SolidScatraEleCalc<celltype, SolidFormulation>::unpack(
    Core::Communication::UnpackBuffer& buffer)
{
  Discret::Elements::unpack(buffer, history_data_);
}

template <Core::FE::CellType celltype, typename SolidFormulation>
void Discret::Elements::SolidScatraEleCalc<celltype,
    SolidFormulation>::evaluate_nonlinear_force_stiffness_mass(const Core::Elements::Element& ele,
    Mat::So3Material& solid_material, const Core::FE::Discretization& discretization,
    const Core::Elements::LocationArray& la, Teuchos::ParameterList& params,
    Core::LinAlg::SerialDenseVector* force_vector,
    Core::LinAlg::SerialDenseMatrix* stiffness_matrix, Core::LinAlg::SerialDenseMatrix* mass_matrix)
{
  // Create views to SerialDenseMatrices
  std::optional<Core::LinAlg::Matrix<num_dof_per_ele_, num_dof_per_ele_>> stiff{};
  std::optional<Core::LinAlg::Matrix<num_dof_per_ele_, num_dof_per_ele_>> mass{};
  std::optional<Core::LinAlg::Matrix<num_dof_per_ele_, 1>> force{};
  if (stiffness_matrix != nullptr) stiff.emplace(*stiffness_matrix, true);
  if (mass_matrix != nullptr) mass.emplace(*mass_matrix, true);
  if (force_vector != nullptr) force.emplace(*force_vector, true);

  const ElementNodes<celltype> nodal_coordinates =
      evaluate_element_nodes<celltype>(ele, discretization, la[0].lm_);

  constexpr bool scalars_are_scalar = false;
  std::optional<std::vector<Core::LinAlg::Matrix<Core::FE::num_nodes(celltype), 1>>> nodal_scalars =
      extract_my_nodal_scalars<celltype, scalars_are_scalar>(
          ele, discretization, la, "scalarfield");

  constexpr bool temperature_is_scalar = true;
  std::optional<Core::LinAlg::Matrix<Core::FE::num_nodes(celltype), 1>> nodal_temperatures =
      extract_my_nodal_scalars<celltype, temperature_is_scalar>(
          ele, discretization, la, "temperature");


  bool equal_integration_mass_stiffness =
      compare_gauss_integration(mass_matrix_integration_, stiffness_matrix_integration_);

  evaluate_centroid_coordinates_and_add_to_parameter_list(nodal_coordinates, params);

  const PreparationData<SolidFormulation> preparation_data =
      prepare(ele, nodal_coordinates, history_data_);

  if constexpr (has_condensed_contribution<SolidFormulation>)
  {
    reset_condensed_variable_integration(ele, nodal_coordinates, preparation_data, history_data_);
  }

  double element_mass = 0.0;
  double element_volume = 0.0;
  const double* total_time =
      params.isParameter("total time") ? &params.get<double>("total time") : nullptr;
  const double* time_step_size =
      params.isParameter("delta time") ? &params.get<double>("delta time") : nullptr;
  for_each_gauss_point(nodal_coordinates, element_properties_, stiffness_matrix_integration_,
      [&](const Core::LinAlg::Tensor<double, Core::FE::dim<celltype>>& xi,
          const ShapeFunctionsAndDerivatives<celltype>& shape_functions,
          const JacobianMapping<celltype>& jacobian_mapping, double integration_factor, int gp)
      {
        prepare_scalar_in_parameter_list(params, "scalars", shape_functions, nodal_scalars);
        prepare_scalar_in_parameter_list(
            params, "temperature", shape_functions, nodal_temperatures);

        evaluate(ele, nodal_coordinates, xi, shape_functions, jacobian_mapping, preparation_data,
            history_data_, gp,
            [&](const Core::LinAlg::Tensor<double, Core::FE::dim<celltype>,
                    Core::FE::dim<celltype>>& deformation_gradient,
                const Core::LinAlg::SymmetricTensor<double, Core::FE::dim<celltype>,
                    Core::FE::dim<celltype>>& gl_strain,
                const auto& linearization)
            {
              auto gp_ref_coord = evaluate_reference_coordinate<celltype>(
                  nodal_coordinates.reference_coordinates, shape_functions.shapefunctions_);

              Mat::EvaluationContext<Core::FE::dim<celltype>> context{.total_time = total_time,
                  .time_step_size = time_step_size,
                  .xi = &xi,
                  .ref_coords = &gp_ref_coord};
              const Stress<celltype> stress =
                  evaluate_material_stress<celltype>(solid_material, element_properties_,
                      deformation_gradient, gl_strain, params, context, gp, ele.id());

              if constexpr (has_condensed_contribution<SolidFormulation>)
              {
                integrate_condensed_contribution(
                    linearization, stress, integration_factor, preparation_data, history_data_, gp);
              }

              if (force.has_value())
              {
                add_internal_force_vector(jacobian_mapping, deformation_gradient, linearization,
                    stress, integration_factor, preparation_data, history_data_, gp, *force);
              }

              if (stiff.has_value())
              {
                add_stiffness_matrix(jacobian_mapping, deformation_gradient, xi, shape_functions,
                    linearization, stress, integration_factor, preparation_data, history_data_, gp,
                    *stiff);
              }

              if (mass.has_value())
              {
                if (equal_integration_mass_stiffness)
                {
                  add_mass_matrix(
                      shape_functions, integration_factor, solid_material.density(gp), *mass);
                }
                else
                {
                  element_mass += solid_material.density(gp) * integration_factor;
                  element_volume += integration_factor;
                }
              }
            });
      });

  if constexpr (has_condensed_contribution<SolidFormulation>)
  {
    const auto condensed_contribution_data =
        prepare_condensed_contribution(preparation_data, history_data_);

    if (force.has_value())
    {
      add_condensed_contribution_to_force_vector<celltype>(
          condensed_contribution_data, preparation_data, history_data_, *force);
    }

    if (stiff.has_value())
    {
      add_condensed_contribution_to_stiffness_matrix<celltype>(
          condensed_contribution_data, preparation_data, history_data_, *stiff);
    }
  }

  if (mass.has_value() && !equal_integration_mass_stiffness)
  {
    // integrate mass matrix
    FOUR_C_ASSERT(element_mass > 0, "It looks like the element mass is 0.0");
    for_each_gauss_point<celltype>(nodal_coordinates, element_properties_, mass_matrix_integration_,
        [&](const Core::LinAlg::Tensor<double, Core::FE::dim<celltype>>& xi,
            const ShapeFunctionsAndDerivatives<celltype>& shape_functions,
            const JacobianMapping<celltype>& jacobian_mapping, double integration_factor, int gp)
        {
          add_mass_matrix(
              shape_functions, integration_factor, element_mass / element_volume, *mass);
        });
  }
}

template <Core::FE::CellType celltype, typename SolidFormulation>
void Discret::Elements::SolidScatraEleCalc<celltype, SolidFormulation>::evaluate_d_stress_d_scalar(
    const Core::Elements::Element& ele, Mat::So3Material& solid_material,
    const Core::FE::Discretization& discretization, const Core::Elements::LocationArray& la,
    Teuchos::ParameterList& params, Core::LinAlg::SerialDenseMatrix& stiffness_matrix_dScalar)
{
  const int scatra_column_stride = std::invoke(
      [&]()
      {
        if (params.isParameter("numscatradofspernode"))
        {
          return params.get<int>("numscatradofspernode");
        }
        return 1;
      });


  const ElementNodes<celltype> nodal_coordinates =
      evaluate_element_nodes<celltype>(ele, discretization, la[0].lm_);

  constexpr bool scalars_are_scalar = false;
  std::optional<std::vector<Core::LinAlg::Matrix<Core::FE::num_nodes(celltype), 1>>> nodal_scalars =
      extract_my_nodal_scalars<celltype, scalars_are_scalar>(
          ele, discretization, la, "scalarfield");

  constexpr bool temperature_is_scalar = true;
  std::optional<Core::LinAlg::Matrix<Core::FE::num_nodes(celltype), 1>> nodal_temperatures =
      extract_my_nodal_scalars<celltype, temperature_is_scalar>(
          ele, discretization, la, "temperature");

  evaluate_centroid_coordinates_and_add_to_parameter_list(nodal_coordinates, params);

  const PreparationData<SolidFormulation> preparation_data =
      prepare(ele, nodal_coordinates, history_data_);

  // Check for negative Jacobian determinants
  ensure_positive_jacobian_determinant_at_element_nodes(nodal_coordinates);

  const double* total_time =
      params.isParameter("total time") ? &params.get<double>("total time") : nullptr;
  const double* time_step_size =
      params.isParameter("delta time") ? &params.get<double>("delta time") : nullptr;
  for_each_gauss_point(nodal_coordinates, element_properties_, stiffness_matrix_integration_,
      [&](const Core::LinAlg::Tensor<double, Core::FE::dim<celltype>>& xi,
          const ShapeFunctionsAndDerivatives<celltype>& shape_functions,
          const JacobianMapping<celltype>& jacobian_mapping, double integration_factor, int gp)
      {
        prepare_scalar_in_parameter_list(params, "scalars", shape_functions, nodal_scalars);
        prepare_scalar_in_parameter_list(
            params, "temperature", shape_functions, nodal_temperatures);

        evaluate(ele, nodal_coordinates, xi, shape_functions, jacobian_mapping, preparation_data,
            history_data_, gp,
            [&](const Core::LinAlg::Tensor<double, Core::FE::dim<celltype>,
                    Core::FE::dim<celltype>>& deformation_gradient,
                const Core::LinAlg::SymmetricTensor<double, Core::FE::dim<celltype>,
                    Core::FE::dim<celltype>>& gl_strain,
                const auto& linearization)
            {
              auto gp_ref_coord = evaluate_reference_coordinate<celltype>(
                  nodal_coordinates.reference_coordinates, shape_functions.shapefunctions_);
              Mat::EvaluationContext<Core::FE::dim<celltype>> context{.total_time = total_time,
                  .time_step_size = time_step_size,
                  .xi = &xi,
                  .ref_coords = &gp_ref_coord};

              const auto dSdc = evaluate_d_material_stress_d_scalars<celltype>(solid_material,
                  element_properties_, deformation_gradient, gl_strain, params, context, gp,
                  ele.id(), scatra_column_stride);

              constexpr int num_dof_per_ele =
                  Core::FE::dim<celltype> * Core::FE::num_nodes(celltype);

              // Assemble matrix
              // k_dS = d(internal force)/dc . N, i.e., the internal force vector of the solid
              // formulation evaluated with dS/dc instead of the stress
              Core::LinAlg::Matrix<num_dof_per_ele, 1> BdSdc(Core::LinAlg::Initialization::zero);
              Stress<celltype> d_stress_d_scalar{};
              for (int k = 0; k < scatra_column_stride; ++k)
              {
                BdSdc.put_scalar(0.0);
                d_stress_d_scalar.pk2_ = dSdc[k];
                add_internal_force_vector(jacobian_mapping, deformation_gradient, linearization,
                    d_stress_d_scalar, integration_factor, preparation_data, history_data_, gp,
                    BdSdc);

                for (int row = 0; row < num_dof_per_ele; ++row)
                {
                  const double BdSdc_row = BdSdc(row, 0);
                  for (int col = 0; col < Core::FE::num_nodes(celltype); ++col)
                  {
                    stiffness_matrix_dScalar(row, col * scatra_column_stride + k) +=
                        BdSdc_row * shape_functions.shapefunctions_(col, 0);
                  }
                }
              }
            });
      });
}

template <Core::FE::CellType celltype, typename SolidFormulation>
void Discret::Elements::SolidScatraEleCalc<celltype,
    SolidFormulation>::evaluate_mechanical_heat_source(const Core::Elements::Element& ele,
    Mat::So3Material& solid_material, const Core::FE::Discretization& discretization,
    const Core::Elements::LocationArray& la,
    const Mat::EvaluationContext<Core::FE::dim<celltype>>& context,
    std::optional<double> d_velocity_d_displacement, Core::LinAlg::SerialDenseVector& heat_source,
    Core::LinAlg::SerialDenseMatrix* d_heat_source_d_temperature,
    Core::LinAlg::SerialDenseMatrix* d_heat_source_d_displacement)
{
  auto* thermo_solid = dynamic_cast<Mat::Trait::ThermoSolid*>(&solid_material);
  if (thermo_solid == nullptr) return;

  // The heat source requires the rate of the strain passed to the material and its
  // linearization, which the solid formulation has to provide.
  if constexpr (!has_strain_rate<SolidFormulation>)
  {
    FOUR_C_THROW(
        "The mechanical heat source is not implemented for the element technology of element {}.",
        ele.id());
  }
  else
  {
    // views on the element vectors and matrices (rows: temperature dofs)
    Core::LinAlg::Matrix<num_nodes_, 1> force(heat_source, true);
    std::optional<Core::LinAlg::Matrix<num_nodes_, num_nodes_>> d_force_d_temperature{};
    std::optional<Core::LinAlg::Matrix<num_nodes_, num_dof_per_ele_>> d_force_d_displacement{};
    if (d_heat_source_d_temperature != nullptr)
      d_force_d_temperature.emplace(*d_heat_source_d_temperature, true);
    if (d_heat_source_d_displacement != nullptr)
    {
      FOUR_C_ASSERT_ALWAYS(d_velocity_d_displacement.has_value(),
          "The derivative of the heat source w.r.t. the displacements requires the derivative of "
          "the velocities w.r.t. the displacements.");
      d_force_d_displacement.emplace(*d_heat_source_d_displacement, true);
    }

    const ElementNodes<celltype> nodal_coordinates =
        evaluate_element_nodes<celltype>(ele, discretization, la[0].lm_);

    FOUR_C_ASSERT_ALWAYS(discretization.has_state(0, "velocity"),
        "The mechanical heat source requires the velocities of element {}.", ele.id());
    const Core::LinAlg::Matrix<num_dof_per_ele_, 1> nodal_velocities(
        Core::FE::extract_values_as_array<num_dof_per_ele_>(
            *discretization.get_state(0, "velocity"), la[0].lm_)
            .data());

    constexpr bool temperature_is_scalar = true;
    const std::optional<Core::LinAlg::Matrix<num_nodes_, 1>> nodal_temperatures =
        extract_my_nodal_scalars<celltype, temperature_is_scalar>(
            ele, discretization, la, "temperature");
    FOUR_C_ASSERT_ALWAYS(nodal_temperatures.has_value(),
        "The mechanical heat source requires the temperatures of element {}.", ele.id());

    const PreparationData<SolidFormulation> preparation_data =
        prepare(ele, nodal_coordinates, history_data_);

    // The thermal elements integrate per unit thickness, hence also the heat source of 2D
    // elements.
    ElementProperties<celltype> unit_thickness_properties = element_properties_;
    if constexpr (num_dim_ == 2) unit_thickness_properties.reference_thickness = 1.0;

    for_each_gauss_point(nodal_coordinates, unit_thickness_properties,
        stiffness_matrix_integration_,
        [&](const Core::LinAlg::Tensor<double, num_dim_>& xi,
            const ShapeFunctionsAndDerivatives<celltype>& shape_functions,
            const JacobianMapping<celltype>& jacobian_mapping, double integration_factor, int gp)
        {
          const double temperature =
              interpolate_quantity_to_point(shape_functions, *nodal_temperatures);

          evaluate(ele, nodal_coordinates, xi, shape_functions, jacobian_mapping, preparation_data,
              history_data_, gp,
              [&](const Core::LinAlg::Tensor<double, num_dim_, num_dim_>& deformation_gradient,
                  const Core::LinAlg::SymmetricTensor<double, num_dim_, num_dim_>& gl_strain,
                  const auto& linearization)
              {
                const auto strain_rate =
                    evaluate_strain_rate<SolidFormulation>(jacobian_mapping, deformation_gradient,
                        linearization, nodal_velocities, preparation_data, history_data_, gp);

                const auto gp_ref_coord = evaluate_reference_coordinate<celltype>(
                    nodal_coordinates.reference_coordinates, shape_functions.shapefunctions_);
                const Mat::EvaluationContext<num_dim_> gp_context{.total_time = context.total_time,
                    .time_step_size = context.time_step_size,
                    .xi = &xi,
                    .ref_coords = &gp_ref_coord};

                const MechanicalHeatSource<celltype> source =
                    evaluate_material_mechanical_heat_source<celltype>(solid_material,
                        *thermo_solid, element_properties_, temperature, deformation_gradient,
                        gl_strain, strain_rate.strain_rate, gp_context, gp, ele.id());

                const auto& N = shape_functions.shapefunctions_;

                // the heat source enters the thermal internal force with a negative sign
                force.update(-integration_factor * source.value, N, 1.0);

                if (d_force_d_temperature.has_value())
                {
                  d_force_d_temperature->multiply_nt(
                      -integration_factor * source.derivative_wrt_temperature, N, N, 1.0);
                }

                if (d_force_d_displacement.has_value())
                {
                  Core::LinAlg::Matrix<num_dof_per_ele_, 1> d_source_d_displacement(
                      Core::LinAlg::Initialization::zero);
                  add_d_strain_d_displacements_contraction<SolidFormulation>(jacobian_mapping,
                      deformation_gradient, linearization, source.derivative_wrt_strain, 1.0,
                      preparation_data, history_data_, gp, d_source_d_displacement);
                  add_d_strain_rate_d_displacements_contraction<SolidFormulation>(jacobian_mapping,
                      deformation_gradient, linearization, strain_rate,
                      source.derivative_wrt_strain_rate, *d_velocity_d_displacement, 1.0,
                      preparation_data, history_data_, gp, d_source_d_displacement);

                  d_force_d_displacement->multiply_nt(
                      -integration_factor, N, d_source_d_displacement, 1.0);
                }
              });
        });
  }
}

template <Core::FE::CellType celltype, typename SolidFormulation>
void Discret::Elements::SolidScatraEleCalc<celltype, SolidFormulation>::recover(
    Core::Elements::Element& ele, const Core::FE::Discretization& discretization,
    const Core::Elements::LocationArray& la, Teuchos::ParameterList& params)
{
  if constexpr (has_condensed_contribution<SolidFormulation>)
  {
    Solid::Elements::ParamsInterface& params_interface =
        *std::dynamic_pointer_cast<Solid::Elements::ParamsInterface>(ele.params_interface_ptr());

    const double step_length = params_interface.get_step_length();

    const ElementNodes<celltype> element_nodes =
        evaluate_element_nodes<celltype>(ele, discretization, la[0].lm_);

    const PreparationData<SolidFormulation> preparation_data =
        prepare(ele, element_nodes, history_data_);

    if (params_interface.is_default_step())
    {
      update_condensed_variables(ele, &params_interface, element_nodes,
          get_displacement_increment<celltype>(discretization, la[0].lm_), step_length,
          preparation_data, history_data_);
    }
    else
    {
      correct_condensed_variables_for_linesearch(
          ele, &params_interface, step_length, preparation_data, history_data_);
    }
  }
}

template <Core::FE::CellType celltype, typename SolidFormulation>
void Discret::Elements::SolidScatraEleCalc<celltype, SolidFormulation>::update(
    const Core::Elements::Element& ele, Mat::So3Material& solid_material,
    const Core::FE::Discretization& discretization, const Core::Elements::LocationArray& la,
    Teuchos::ParameterList& params)
{
  const ElementNodes<celltype> nodal_coordinates =
      evaluate_element_nodes<celltype>(ele, discretization, la[0].lm_);

  constexpr bool scalars_are_scalar = false;
  std::optional<std::vector<Core::LinAlg::Matrix<Core::FE::num_nodes(celltype), 1>>> nodal_scalars =
      extract_my_nodal_scalars<celltype, scalars_are_scalar>(
          ele, discretization, la, "scalarfield");

  constexpr bool temperature_is_scalar = true;
  std::optional<Core::LinAlg::Matrix<Core::FE::num_nodes(celltype), 1>> nodal_temperatures =
      extract_my_nodal_scalars<celltype, temperature_is_scalar>(
          ele, discretization, la, "temperature");

  evaluate_centroid_coordinates_and_add_to_parameter_list(nodal_coordinates, params);

  const PreparationData<SolidFormulation> preparation_data =
      prepare(ele, nodal_coordinates, history_data_);

  const double* total_time =
      params.isParameter("total time") ? &params.get<double>("total time") : nullptr;
  const double* time_step_size =
      params.isParameter("delta time") ? &params.get<double>("delta time") : nullptr;
  Discret::Elements::for_each_gauss_point(nodal_coordinates, element_properties_,
      stiffness_matrix_integration_,
      [&](const Core::LinAlg::Tensor<double, Core::FE::dim<celltype>>& xi,
          const ShapeFunctionsAndDerivatives<celltype>& shape_functions,
          const JacobianMapping<celltype>& jacobian_mapping, double integration_factor, int gp)
      {
        prepare_scalar_in_parameter_list(params, "scalars", shape_functions, nodal_scalars);
        prepare_scalar_in_parameter_list(
            params, "temperature", shape_functions, nodal_temperatures);


        auto gp_ref_coord = evaluate_reference_coordinate<celltype>(
            nodal_coordinates.reference_coordinates, shape_functions.shapefunctions_);

        Mat::EvaluationContext<Core::FE::dim<celltype>> context{.total_time = total_time,
            .time_step_size = time_step_size,
            .xi = &xi,
            .ref_coords = &gp_ref_coord};
        evaluate(ele, nodal_coordinates, xi, shape_functions, jacobian_mapping, preparation_data,
            history_data_, gp,
            [&](const Core::LinAlg::Tensor<double, Core::FE::dim<celltype>,
                    Core::FE::dim<celltype>>& deformation_gradient,
                const Core::LinAlg::SymmetricTensor<double, Core::FE::dim<celltype>,
                    Core::FE::dim<celltype>>& gl_strain,
                const auto& linearization)
            {
              update_material(solid_material, element_properties_, deformation_gradient, params,
                  context, gp, ele.id());
            });
      });

  solid_material.update();
}

template <Core::FE::CellType celltype, typename SolidFormulation>
double Discret::Elements::SolidScatraEleCalc<celltype, SolidFormulation>::calculate_internal_energy(
    const Core::Elements::Element& ele, Mat::So3Material& solid_material,
    const Core::FE::Discretization& discretization, const Core::Elements::LocationArray& la,
    Teuchos::ParameterList& params)
{
  const ElementNodes<celltype> nodal_coordinates =
      evaluate_element_nodes<celltype>(ele, discretization, la[0].lm_);

  constexpr bool scalars_are_scalar = false;
  std::optional<std::vector<Core::LinAlg::Matrix<Core::FE::num_nodes(celltype), 1>>> nodal_scalars =
      extract_my_nodal_scalars<celltype, scalars_are_scalar>(
          ele, discretization, la, "scalarfield");

  constexpr bool temperature_is_scalar = true;
  std::optional<Core::LinAlg::Matrix<Core::FE::num_nodes(celltype), 1>> nodal_temperatures =
      extract_my_nodal_scalars<celltype, temperature_is_scalar>(
          ele, discretization, la, "temperature");

  evaluate_centroid_coordinates_and_add_to_parameter_list(nodal_coordinates, params);

  const PreparationData<SolidFormulation> preparation_data =
      prepare(ele, nodal_coordinates, history_data_);

  double intenergy = 0;
  const double* total_time =
      params.isParameter("total time") ? &params.get<double>("total time") : nullptr;
  const double* time_step_size =
      params.isParameter("delta time") ? &params.get<double>("delta time") : nullptr;
  Discret::Elements::for_each_gauss_point(nodal_coordinates, element_properties_,
      stiffness_matrix_integration_,
      [&](const Core::LinAlg::Tensor<double, Core::FE::dim<celltype>>& xi,
          const ShapeFunctionsAndDerivatives<celltype>& shape_functions,
          const JacobianMapping<celltype>& jacobian_mapping, double integration_factor, int gp)
      {
        prepare_scalar_in_parameter_list(params, "scalars", shape_functions, nodal_scalars);
        prepare_scalar_in_parameter_list(
            params, "temperature", shape_functions, nodal_temperatures);

        evaluate(ele, nodal_coordinates, xi, shape_functions, jacobian_mapping, preparation_data,
            history_data_, gp,
            [&](const Core::LinAlg::Tensor<double, Core::FE::dim<celltype>,
                    Core::FE::dim<celltype>>& deformation_gradient,
                const Core::LinAlg::SymmetricTensor<double, Core::FE::dim<celltype>,
                    Core::FE::dim<celltype>>& gl_strain,
                const auto& linearization)
            {
              auto gp_ref_coord = evaluate_reference_coordinate<celltype>(
                  nodal_coordinates.reference_coordinates, shape_functions.shapefunctions_);

              Mat::EvaluationContext<Core::FE::dim<celltype>> context{.total_time = total_time,
                  .time_step_size = time_step_size,
                  .xi = &xi,
                  .ref_coords = &gp_ref_coord};
              double psi = evaluate_material_strain_energy<celltype>(
                  solid_material, element_properties_, gl_strain, params, context, gp, ele.id());
              intenergy += psi * integration_factor;
            });
      });

  return intenergy;
}

template <Core::FE::CellType celltype, typename SolidFormulation>
void Discret::Elements::SolidScatraEleCalc<celltype, SolidFormulation>::calculate_stress(
    const Core::Elements::Element& ele, Mat::So3Material& solid_material, const StressIO& stressIO,
    const StrainIO& strainIO, const Core::FE::Discretization& discretization,
    const Core::Elements::LocationArray& la, Teuchos::ParameterList& params)
{
  std::vector<char>& serialized_stress_data = stressIO.mutable_data;
  std::vector<char>& serialized_strain_data = strainIO.mutable_data;
  constexpr std::size_t num_str_for_output = 6;
  Core::LinAlg::SerialDenseMatrix stress_data(
      stiffness_matrix_integration_.num_points(), num_str_for_output);
  Core::LinAlg::SerialDenseMatrix strain_data(
      stiffness_matrix_integration_.num_points(), num_str_for_output);

  const ElementNodes<celltype> nodal_coordinates =
      evaluate_element_nodes<celltype>(ele, discretization, la[0].lm_);

  constexpr bool scalars_are_scalar = false;
  std::optional<std::vector<Core::LinAlg::Matrix<Core::FE::num_nodes(celltype), 1>>> nodal_scalars =
      extract_my_nodal_scalars<celltype, scalars_are_scalar>(
          ele, discretization, la, "scalarfield");

  constexpr bool temperature_is_scalar = true;
  std::optional<Core::LinAlg::Matrix<Core::FE::num_nodes(celltype), 1>> nodal_temperatures =
      extract_my_nodal_scalars<celltype, temperature_is_scalar>(
          ele, discretization, la, "temperature");

  evaluate_centroid_coordinates_and_add_to_parameter_list(nodal_coordinates, params);

  const PreparationData<SolidFormulation> preparation_data =
      prepare(ele, nodal_coordinates, history_data_);

  const double* total_time =
      params.isParameter("total time") ? &params.get<double>("total time") : nullptr;
  const double* time_step_size =
      params.isParameter("delta time") ? &params.get<double>("delta time") : nullptr;
  Discret::Elements::for_each_gauss_point(nodal_coordinates, element_properties_,
      stiffness_matrix_integration_,
      [&](const Core::LinAlg::Tensor<double, Core::FE::dim<celltype>>& xi,
          const ShapeFunctionsAndDerivatives<celltype>& shape_functions,
          const JacobianMapping<celltype>& jacobian_mapping, double integration_factor, int gp)
      {
        prepare_scalar_in_parameter_list(params, "scalars", shape_functions, nodal_scalars);
        prepare_scalar_in_parameter_list(
            params, "temperature", shape_functions, nodal_temperatures);

        evaluate(ele, nodal_coordinates, xi, shape_functions, jacobian_mapping, preparation_data,
            history_data_, gp,
            [&](const Core::LinAlg::Tensor<double, Core::FE::dim<celltype>,
                    Core::FE::dim<celltype>>& deformation_gradient,
                const Core::LinAlg::SymmetricTensor<double, Core::FE::dim<celltype>,
                    Core::FE::dim<celltype>>& gl_strain,
                const auto& linearization)
            {
              auto gp_ref_coord = evaluate_reference_coordinate<celltype>(
                  nodal_coordinates.reference_coordinates, shape_functions.shapefunctions_);

              Mat::EvaluationContext<Core::FE::dim<celltype>> context{.total_time = total_time,
                  .time_step_size = time_step_size,
                  .xi = &xi,
                  .ref_coords = &gp_ref_coord};
              const Stress<celltype> stress =
                  evaluate_material_stress<celltype>(solid_material, element_properties_,
                      deformation_gradient, gl_strain, params, context, gp, ele.id());

              assemble_strain_type_to_matrix_row<celltype>(element_properties_, gl_strain, stress,
                  deformation_gradient, strainIO.type, strain_data, gp);
              assemble_stress_type_to_matrix_row(element_properties_, deformation_gradient, stress,
                  stressIO.type, stress_data, gp);
            });
      });

  serialize(stress_data, serialized_stress_data);
  serialize(strain_data, serialized_strain_data);
}

template <Core::FE::CellType celltype, typename SolidFormulation>
double
Discret::Elements::SolidScatraEleCalc<celltype, SolidFormulation>::get_normal_cauchy_stress_at_xi(
    const Core::Elements::Element& ele, Mat::So3Material& solid_material,
    const std::vector<double>& disp, const std::vector<double>& scalars,
    const Core::LinAlg::Tensor<double, 3>& xi, const Core::LinAlg::Tensor<double, 3>& n,
    const Core::LinAlg::Tensor<double, 3>& dir,
    SolidScatraCauchyNDirLinearizations<3>& linearizations)
{
  if constexpr (has_gauss_point_history<SolidFormulation>)
  {
    FOUR_C_THROW(
        "Cannot evaluate the Cauchy stress at xi with an element formulation with Gauss point "
        "history. The element formulation is {}.",
        Core::Utils::get_type_name<SolidFormulation>().c_str());
  }
  else if constexpr (Core::FE::dim<celltype> != 3)
  {
    FOUR_C_THROW(
        "Cannot evaluate the Cauchy stress at xi for an element formulation with spatial "
        "dimension different than 3. The element formulation is {}.",
        Core::Utils::get_type_name<SolidFormulation>().c_str());
  }
  else if constexpr (Core::FE::is_nurbs<celltype>)
  {
    FOUR_C_THROW("Cannot evaluate the Cauchy stress at xi for NURBS elements.");
  }
  else
  {
    // project scalar values to xi
    const auto scalar_values_at_xi = Core::FE::interpolate_to_xi<celltype>(
        Core::LinAlg::make_matrix_view<Core::FE::dim<celltype>, 1>(xi), scalars);


    ElementNodes<celltype> element_nodes = evaluate_element_nodes<celltype>(ele, disp);

    const ShapeFunctionsAndDerivatives<celltype> shape_functions =
        evaluate_shape_functions_and_derivs<celltype>(xi, element_nodes);

    const JacobianMapping<celltype> jacobian_mapping =
        evaluate_jacobian_mapping(shape_functions, element_nodes);

    const PreparationData<SolidFormulation> preparation_data =
        prepare(ele, element_nodes, history_data_);

    return evaluate(ele, element_nodes, xi, shape_functions, jacobian_mapping, preparation_data,
        history_data_,
        [&](const Core::LinAlg::Tensor<double, Core::FE::dim<celltype>, Core::FE::dim<celltype>>&
                deformation_gradient,
            const Core::LinAlg::SymmetricTensor<double, Core::FE::dim<celltype>,
                Core::FE::dim<celltype>>& gl_strain,
            const auto& linearization)
        {
          const ElementFormulationDerivativeEvaluator<celltype, SolidFormulation> evaluator(ele,
              element_nodes, xi, shape_functions, jacobian_mapping, deformation_gradient,
              preparation_data, history_data_);

          return evaluate_cauchy_n_dir_at_xi<celltype>(solid_material, shape_functions, xi,
              deformation_gradient, scalar_values_at_xi, n, dir, ele.id(), evaluator,
              linearizations);
        });
  }
}

template <Core::FE::CellType celltype, typename SolidFormulation>
void Discret::Elements::SolidScatraEleCalc<celltype, SolidFormulation>::setup(
    Mat::So3Material& solid_material, const Core::IO::InputParameterContainer& container)
{
  solid_material.setup(stiffness_matrix_integration_.num_points(), read_fibers(container),
      read_coordinate_system(container));
}

template <Core::FE::CellType celltype, typename SolidFormulation>
void Discret::Elements::SolidScatraEleCalc<celltype, SolidFormulation>::material_post_setup(
    const Core::Elements::Element& ele, Mat::So3Material& solid_material)
{
  Teuchos::ParameterList params{};

  // Check if element has fiber nodes, if so interpolate fibers to Gauss Points and add to params
  interpolate_fibers_to_gauss_points_and_add_to_parameter_list<celltype>(
      stiffness_matrix_integration_, ele, params);

  // Call post_setup of material
  solid_material.post_setup(params, ele.id());
}

template <Core::FE::CellType celltype, typename SolidFormulation>
void Discret::Elements::SolidScatraEleCalc<celltype,
    SolidFormulation>::initialize_gauss_point_data_output(const Core::Elements::Element& ele,
    const Mat::So3Material& solid_material,
    Solid::ModelEvaluator::GaussPointDataOutputManager& gp_data_output_manager) const
{
  FOUR_C_ASSERT(ele.is_params_interface(),
      "This action type should only be called from the new time integration framework!");

  ask_and_add_quantities_to_gauss_point_data_output(
      stiffness_matrix_integration_.num_points(), solid_material, gp_data_output_manager);
}

template <Core::FE::CellType celltype, typename SolidFormulation>
void Discret::Elements::SolidScatraEleCalc<celltype,
    SolidFormulation>::evaluate_gauss_point_data_output(const Core::Elements::Element& ele,
    const Mat::So3Material& solid_material,
    Solid::ModelEvaluator::GaussPointDataOutputManager& gp_data_output_manager) const
{
  FOUR_C_ASSERT(ele.is_params_interface(),
      "This action type should only be called from the new time integration framework!");

  collect_and_assemble_gauss_point_data_output<celltype>(
      stiffness_matrix_integration_, solid_material, ele, gp_data_output_manager);
}

template <Core::FE::CellType celltype, typename SolidFormulation>
void Discret::Elements::SolidScatraEleCalc<celltype, SolidFormulation>::reset_to_last_converged(
    const Core::Elements::Element& ele, Mat::So3Material& solid_material)
{
  solid_material.reset_step();
}

template <Core::FE::CellType... celltypes>
struct VerifyPackable
{
  static constexpr bool are_all_packable =
      (Core::Communication::Packable<Discret::Elements::SolidScatraEleCalc<celltypes,
              Discret::Elements::DisplacementBasedFormulation<celltypes>>> &&
          ...);

  static constexpr bool are_all_unpackable =
      (Core::Communication::Unpackable<Discret::Elements::SolidScatraEleCalc<celltypes,
              Discret::Elements::DisplacementBasedFormulation<celltypes>>> &&
          ...);

  void static_asserts() const
  {
    static_assert(are_all_packable);
    static_assert(are_all_unpackable);
  }
};

template struct VerifyPackable<Core::FE::CellType::hex8, Core::FE::CellType::hex27,
    Core::FE::CellType::tet4, Core::FE::CellType::tet10>;

// explicit instantiations of template classes
// for displacement based formulation
template class Discret::Elements::SolidScatraEleCalc<Core::FE::CellType::hex8,
    Discret::Elements::DisplacementBasedFormulation<Core::FE::CellType::hex8>>;
template class Discret::Elements::SolidScatraEleCalc<Core::FE::CellType::hex27,
    Discret::Elements::DisplacementBasedFormulation<Core::FE::CellType::hex27>>;
template class Discret::Elements::SolidScatraEleCalc<Core::FE::CellType::tet4,
    Discret::Elements::DisplacementBasedFormulation<Core::FE::CellType::tet4>>;
template class Discret::Elements::SolidScatraEleCalc<Core::FE::CellType::tet10,
    Discret::Elements::DisplacementBasedFormulation<Core::FE::CellType::tet10>>;
template class Discret::Elements::SolidScatraEleCalc<Core::FE::CellType::nurbs27,
    Discret::Elements::DisplacementBasedFormulation<Core::FE::CellType::nurbs27>>;

template class Discret::Elements::SolidScatraEleCalc<Core::FE::CellType::quad4,
    Discret::Elements::DisplacementBasedFormulation<Core::FE::CellType::quad4>>;
template class Discret::Elements::SolidScatraEleCalc<Core::FE::CellType::quad9,
    Discret::Elements::DisplacementBasedFormulation<Core::FE::CellType::quad9>>;
template class Discret::Elements::SolidScatraEleCalc<Core::FE::CellType::tri3,
    Discret::Elements::DisplacementBasedFormulation<Core::FE::CellType::tri3>>;
template class Discret::Elements::SolidScatraEleCalc<Core::FE::CellType::tri6,
    Discret::Elements::DisplacementBasedFormulation<Core::FE::CellType::tri6>>;

// for displacement based formulation with linear kinematics
template class Discret::Elements::SolidScatraEleCalc<Core::FE::CellType::hex8,
    Discret::Elements::DisplacementBasedLinearKinematicsFormulation<Core::FE::CellType::hex8>>;
template class Discret::Elements::SolidScatraEleCalc<Core::FE::CellType::hex27,
    Discret::Elements::DisplacementBasedLinearKinematicsFormulation<Core::FE::CellType::hex27>>;
template class Discret::Elements::SolidScatraEleCalc<Core::FE::CellType::tet4,
    Discret::Elements::DisplacementBasedLinearKinematicsFormulation<Core::FE::CellType::tet4>>;
template class Discret::Elements::SolidScatraEleCalc<Core::FE::CellType::tet10,
    Discret::Elements::DisplacementBasedLinearKinematicsFormulation<Core::FE::CellType::tet10>>;
template class Discret::Elements::SolidScatraEleCalc<Core::FE::CellType::nurbs27,
    Discret::Elements::DisplacementBasedLinearKinematicsFormulation<Core::FE::CellType::nurbs27>>;

template class Discret::Elements::SolidScatraEleCalc<Core::FE::CellType::quad4,
    Discret::Elements::DisplacementBasedLinearKinematicsFormulation<Core::FE::CellType::quad4>>;
template class Discret::Elements::SolidScatraEleCalc<Core::FE::CellType::quad9,
    Discret::Elements::DisplacementBasedLinearKinematicsFormulation<Core::FE::CellType::quad9>>;
template class Discret::Elements::SolidScatraEleCalc<Core::FE::CellType::tri3,
    Discret::Elements::DisplacementBasedLinearKinematicsFormulation<Core::FE::CellType::tri3>>;
template class Discret::Elements::SolidScatraEleCalc<Core::FE::CellType::tri6,
    Discret::Elements::DisplacementBasedLinearKinematicsFormulation<Core::FE::CellType::tri6>>;


// FBar based formulation
template class Discret::Elements::SolidScatraEleCalc<Core::FE::CellType::hex8,
    Discret::Elements::FBarFormulation<Core::FE::CellType::hex8>>;

// explicit instantiations for hex8 with EAS
template class Discret::Elements::SolidScatraEleCalc<Core::FE::CellType::hex8,
    Discret::Elements::EASFormulation<Core::FE::CellType::hex8,
        Discret::Elements::EasType::eastype_h8_9, Solid::KinemType::nonlinearTotLag>>;
template class Discret::Elements::SolidScatraEleCalc<Core::FE::CellType::hex8,
    Discret::Elements::EASFormulation<Core::FE::CellType::hex8,
        Discret::Elements::EasType::eastype_h8_21, Solid::KinemType::nonlinearTotLag>>;

// explicit instantiations for quad4 with EAS
template class Discret::Elements::SolidScatraEleCalc<Core::FE::CellType::quad4,
    Discret::Elements::EASFormulation<Core::FE::CellType::quad4,
        Discret::Elements::EasType::eastype_q4_4, Solid::KinemType::nonlinearTotLag>>;


FOUR_C_NAMESPACE_CLOSE