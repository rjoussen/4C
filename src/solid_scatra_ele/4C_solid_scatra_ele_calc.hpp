// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_SOLID_SCATRA_ELE_CALC_HPP
#define FOUR_C_SOLID_SCATRA_ELE_CALC_HPP

#include "4C_config.hpp"

#include "4C_fem_general_element.hpp"
#include "4C_fem_general_utils_gausspoints.hpp"
#include "4C_solid_ele_calc_interface.hpp"
#include "4C_solid_ele_formulation.hpp"
#include "4C_solid_scatra_ele_calc_lib_nitsche.hpp"
#include "4C_structure_new_gauss_point_data_output_manager.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Mat
{
  class So3Material;
}  // namespace Mat
namespace Solid::Elements
{
  class ParamsInterface;
}  // namespace Solid::Elements



namespace Discret::Elements
{
  template <Core::FE::CellType celltype, typename SolidFormulation>
  class SolidScatraEleCalc
  {
   public:
    SolidScatraEleCalc()
      requires(Core::FE::dim<celltype> == 3);
    SolidScatraEleCalc(
        const double reference_thickness, const Discret::Elements::PlaneAssumption plane_assumption)
      requires(Core::FE::dim<celltype> == 2);

    void pack(Core::Communication::PackBuffer& data) const;

    void unpack(Core::Communication::UnpackBuffer& buffer);

    void setup(
        Mat::So3Material& solid_material, const Core::IO::InputParameterContainer& container);

    void material_post_setup(const Core::Elements::Element& ele, Mat::So3Material& solid_material);

    void recover(Core::Elements::Element& ele, const Core::FE::Discretization& discretization,
        const Core::Elements::LocationArray& la, Teuchos::ParameterList& params);

    void update(const Core::Elements::Element& ele, Mat::So3Material& solid_material,
        const Core::FE::Discretization& discretization, const Core::Elements::LocationArray& la,
        Teuchos::ParameterList& params);

    void calculate_stress(const Core::Elements::Element& ele, Mat::So3Material& solid_material,
        const StressIO& stressIO, const StrainIO& strainIO,
        const Core::FE::Discretization& discretization, const Core::Elements::LocationArray& la,
        Teuchos::ParameterList& params);

    double calculate_internal_energy(const Core::Elements::Element& ele,
        Mat::So3Material& solid_material, const Core::FE::Discretization& discretization,
        const Core::Elements::LocationArray& la, Teuchos::ParameterList& params);

    void initialize_gauss_point_data_output(const Core::Elements::Element& ele,
        const Mat::So3Material& solid_material,
        FourC::Solid::ModelEvaluator::GaussPointDataOutputManager& gp_data_output_manager) const;

    void evaluate_gauss_point_data_output(const Core::Elements::Element& ele,
        const Mat::So3Material& solid_material,
        FourC::Solid::ModelEvaluator::GaussPointDataOutputManager& gp_data_output_manager) const;

    void reset_to_last_converged(
        const Core::Elements::Element& ele, Mat::So3Material& solid_material);

    void evaluate_nonlinear_force_stiffness_mass(const Core::Elements::Element& ele,
        Mat::So3Material& solid_material, const Core::FE::Discretization& discretization,
        const Core::Elements::LocationArray& la, Teuchos::ParameterList& params,
        Core::LinAlg::SerialDenseVector* force_vector,
        Core::LinAlg::SerialDenseMatrix* stiffness_matrix,
        Core::LinAlg::SerialDenseMatrix* mass_matrix);

    void evaluate_d_stress_d_scalar(const Core::Elements::Element& ele,
        Mat::So3Material& solid_material, const Core::FE::Discretization& discretization,
        const Core::Elements::LocationArray& la, Teuchos::ParameterList& params,
        Core::LinAlg::SerialDenseMatrix& stiffness_matrix_dScalar);

    /*!
     * @brief Evaluate the mechanical heat source of a thermo-solid material for the thermal
     * equation
     *
     * The heat source is evaluated with the kinematics that the solid formulation passes to the
     * material. Materials with internal variables evaluate it at the state of their last
     * evaluation, i.e., this has to be called after the evaluation of the stress at the same
     * state. The rows of all outputs refer to the temperature dofs (second dofset).
     *
     * @param context (in) : evaluation context, the point-wise entries are set by the element
     * @param d_velocity_d_displacement (in) : derivative of the velocities w.r.t. the
     * displacements due to the time integration (required for @p d_heat_source_d_displacement)
     * @param heat_source (out) : contribution to the thermal internal force
     * @param d_heat_source_d_temperature (out) : its derivative w.r.t. the temperatures
     * (optional)
     * @param d_heat_source_d_displacement (out) : its derivative w.r.t. the displacements
     * (optional)
     */
    void evaluate_mechanical_heat_source(const Core::Elements::Element& ele,
        Mat::So3Material& solid_material, const Core::FE::Discretization& discretization,
        const Core::Elements::LocationArray& la,
        const Mat::EvaluationContext<Core::FE::dim<celltype>>& context,
        std::optional<double> d_velocity_d_displacement,
        Core::LinAlg::SerialDenseVector& heat_source,
        Core::LinAlg::SerialDenseMatrix* d_heat_source_d_temperature,
        Core::LinAlg::SerialDenseMatrix* d_heat_source_d_displacement);

    double get_normal_cauchy_stress_at_xi(const Core::Elements::Element& ele,
        Mat::So3Material& solid_material, const std::vector<double>& disp,
        const std::vector<double>& scalars, const Core::LinAlg::Tensor<double, 3>& xi,
        const Core::LinAlg::Tensor<double, 3>& n, const Core::LinAlg::Tensor<double, 3>& dir,
        SolidScatraCauchyNDirLinearizations<3>& linearizations);

   private:
    /// static values for matrix sizes
    static constexpr int num_nodes_ = Core::FE::num_nodes(celltype);
    static constexpr int num_dim_ = Core::FE::dim<celltype>;
    static constexpr int num_dof_per_ele_ = num_nodes_ * num_dim_;
    static constexpr int num_str_ = num_dim_ * (num_dim_ + 1) / 2;

    Core::FE::GaussIntegration stiffness_matrix_integration_;
    Core::FE::GaussIntegration mass_matrix_integration_;

    ElementProperties<celltype> element_properties_;

    SolidFormulationHistory<SolidFormulation> history_data_{};
  };
}  // namespace Discret::Elements


FOUR_C_NAMESPACE_CLOSE

#endif