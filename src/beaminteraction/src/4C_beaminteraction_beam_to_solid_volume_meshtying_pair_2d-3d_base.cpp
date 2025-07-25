// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_beaminteraction_beam_to_solid_volume_meshtying_pair_2d-3d_base.hpp"

#include "4C_geometry_pair_element_evaluation_functions.hpp"
#include "4C_geometry_pair_line_to_3D_evaluation_data.hpp"
#include "4C_geometry_pair_line_to_volume_gauss_point_projection_cross_section.hpp"
#include "4C_geometry_pair_utility_classes.hpp"

FOUR_C_NAMESPACE_OPEN


/**
 *
 */
template <typename ScalarType, typename Beam, typename Solid>
void BeamInteraction::BeamToSolidVolumeMeshtyingPair2D3DBase<ScalarType, Beam,
    Solid>::create_geometry_pair(const Core::Elements::Element* element1,
    const Core::Elements::Element* element2,
    const std::shared_ptr<GeometryPair::GeometryEvaluationDataBase>& geometry_evaluation_data_ptr)
{
  // Cast the geometry evaluation data to the correct format.
  auto line_to_3d_evaluation_data =
      std::dynamic_pointer_cast<GeometryPair::LineTo3DEvaluationData>(geometry_evaluation_data_ptr);

  // Explicitly create the cross section projection geometry pair here and check that the correct
  // parameter is set in the input file.
  GeometryPair::LineTo3DStrategy strategy = line_to_3d_evaluation_data->get_strategy();
  if (strategy != GeometryPair::LineTo3DStrategy::gauss_point_projection_cross_section)
    FOUR_C_THROW(
        "The 2D-3D beam-to-volume mesh tying pair only works with the cross section projection "
        "geometry pair. This has to be specified in the input file.");
  this->geometry_pair_ = std::make_shared<
      GeometryPair::GeometryPairLineToVolumeGaussPointProjectionCrossSection<double, Beam, Solid>>(
      element1, element2, line_to_3d_evaluation_data);
}

/**
 *
 */
template <typename ScalarType, typename Beam, typename Solid>
void BeamInteraction::BeamToSolidVolumeMeshtyingPair2D3DBase<ScalarType, Beam,
    Solid>::evaluate_beam_position_double(const GeometryPair::ProjectionPoint1DTo3D<double>&
                                              integration_point,
    Core::LinAlg::Matrix<3, 1, double>& r_beam, bool reference) const
{
  auto evaluate_position = [&](const auto& q, auto& r_beam)
  {
    const auto eta = integration_point.get_eta();
    Core::LinAlg::Matrix<3, 3, double> triad;
    get_triad_at_xi_double(eta, triad, reference);
    Core::LinAlg::Matrix<3, 1, double> r_cross_section_ref, r_cross_section_cur;
    r_cross_section_ref(0) = 0.0;
    r_cross_section_ref(1) = integration_point.get_eta_cross_section()(0);
    r_cross_section_ref(2) = integration_point.get_eta_cross_section()(1);
    r_cross_section_cur.multiply(triad, r_cross_section_ref);
    GeometryPair::evaluate_position<Beam>(eta, q, r_beam);
    r_beam += r_cross_section_cur;
  };

  if (reference)
  {
    evaluate_position(this->ele1posref_, r_beam);
  }
  else
  {
    evaluate_position(GeometryPair::ElementDataToDouble<Beam>::to_double(this->ele1pos_), r_beam);
  }
}


/**
 * Explicit template initialization of template class.
 */
namespace BeamInteraction
{
  using namespace GeometryPair;

  template class BeamToSolidVolumeMeshtyingPair2D3DBase<double, t_hermite, t_hex8>;
  template class BeamToSolidVolumeMeshtyingPair2D3DBase<double, t_hermite, t_hex20>;
  template class BeamToSolidVolumeMeshtyingPair2D3DBase<double, t_hermite, t_hex27>;
  template class BeamToSolidVolumeMeshtyingPair2D3DBase<double, t_hermite, t_tet4>;
  template class BeamToSolidVolumeMeshtyingPair2D3DBase<double, t_hermite, t_tet10>;

  template class BeamToSolidVolumeMeshtyingPair2D3DBase<
      line_to_volume_scalar_type<t_hermite, t_hex8>, t_hermite, t_hex8>;
  template class BeamToSolidVolumeMeshtyingPair2D3DBase<
      line_to_volume_scalar_type<t_hermite, t_hex20>, t_hermite, t_hex20>;
  template class BeamToSolidVolumeMeshtyingPair2D3DBase<
      line_to_volume_scalar_type<t_hermite, t_hex27>, t_hermite, t_hex27>;
  template class BeamToSolidVolumeMeshtyingPair2D3DBase<
      line_to_volume_scalar_type<t_hermite, t_tet4>, t_hermite, t_tet4>;
  template class BeamToSolidVolumeMeshtyingPair2D3DBase<
      line_to_volume_scalar_type<t_hermite, t_tet10>, t_hermite, t_tet10>;
}  // namespace BeamInteraction

FOUR_C_NAMESPACE_CLOSE
