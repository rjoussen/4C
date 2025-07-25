// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_geometry_pair_factory.hpp"

#include "4C_geometry_pair_element.hpp"
#include "4C_geometry_pair_line_to_3D_evaluation_data.hpp"
#include "4C_geometry_pair_line_to_surface_evaluation_data.hpp"
#include "4C_geometry_pair_line_to_surface_gauss_point_projection.hpp"
#include "4C_geometry_pair_line_to_surface_segmentation.hpp"
#include "4C_geometry_pair_line_to_volume_gauss_point_projection.hpp"
#include "4C_geometry_pair_line_to_volume_segmentation.hpp"
#include "4C_geometry_pair_scalar_types.hpp"

FOUR_C_NAMESPACE_OPEN

/**
 *
 */
template <typename ScalarType, typename Line, typename Volume>
std::shared_ptr<GeometryPair::GeometryPairBase> GeometryPair::geometry_pair_line_to_volume_factory(
    const Core::Elements::Element* element1, const Core::Elements::Element* element2,
    const std::shared_ptr<GeometryEvaluationDataBase>& geometry_evaluation_data)
{
  // Cast the geometry evaluation data to the correct format.
  auto line_to_3d_evaluation_data =
      std::dynamic_pointer_cast<LineTo3DEvaluationData>(geometry_evaluation_data);

  // Get the strategy for line to volume interaction.
  GeometryPair::LineTo3DStrategy strategy = line_to_3d_evaluation_data->get_strategy();

  // Create the class depending on the strategy.
  switch (strategy)
  {
    case GeometryPair::LineTo3DStrategy::gauss_point_projection_without_boundary_segmentation:
    case GeometryPair::LineTo3DStrategy::gauss_point_projection_boundary_segmentation:
      return std::make_shared<
          GeometryPairLineToVolumeGaussPointProjection<ScalarType, Line, Volume>>(
          element1, element2, line_to_3d_evaluation_data);
    case GeometryPair::LineTo3DStrategy::segmentation:
      return std::make_shared<GeometryPairLineToVolumeSegmentation<ScalarType, Line, Volume>>(
          element1, element2, line_to_3d_evaluation_data);
    default:
    {
      FOUR_C_THROW(
          "The given geometry pair strategy is not universally valid. You might want to create "
          "your pair directly if you need certain features (for example cross section "
          "projection)!");
      return nullptr;
    }
  }
}


/**
 * Explicit template initialization of factory function.
 */
template std::shared_ptr<GeometryPair::GeometryPairBase>
GeometryPair::geometry_pair_line_to_volume_factory<double, GeometryPair::t_hermite,
    GeometryPair::t_hex8>(const Core::Elements::Element*, const Core::Elements::Element*,
    const std::shared_ptr<GeometryEvaluationDataBase>&);
template std::shared_ptr<GeometryPair::GeometryPairBase>
GeometryPair::geometry_pair_line_to_volume_factory<double, GeometryPair::t_hermite,
    GeometryPair::t_hex20>(const Core::Elements::Element*, const Core::Elements::Element*,
    const std::shared_ptr<GeometryEvaluationDataBase>&);
template std::shared_ptr<GeometryPair::GeometryPairBase>
GeometryPair::geometry_pair_line_to_volume_factory<double, GeometryPair::t_hermite,
    GeometryPair::t_hex27>(const Core::Elements::Element*, const Core::Elements::Element*,
    const std::shared_ptr<GeometryEvaluationDataBase>&);
template std::shared_ptr<GeometryPair::GeometryPairBase>
GeometryPair::geometry_pair_line_to_volume_factory<double, GeometryPair::t_hermite,
    GeometryPair::t_tet4>(const Core::Elements::Element*, const Core::Elements::Element*,
    const std::shared_ptr<GeometryEvaluationDataBase>&);
template std::shared_ptr<GeometryPair::GeometryPairBase>
GeometryPair::geometry_pair_line_to_volume_factory<double, GeometryPair::t_hermite,
    GeometryPair::t_tet10>(const Core::Elements::Element*, const Core::Elements::Element*,
    const std::shared_ptr<GeometryEvaluationDataBase>&);
template std::shared_ptr<GeometryPair::GeometryPairBase>
GeometryPair::geometry_pair_line_to_volume_factory<double, GeometryPair::t_hermite,
    GeometryPair::t_nurbs27>(const Core::Elements::Element*, const Core::Elements::Element*,
    const std::shared_ptr<GeometryEvaluationDataBase>&);
template std::shared_ptr<GeometryPair::GeometryPairBase>
GeometryPair::geometry_pair_line_to_volume_factory<double, GeometryPair::t_hermite,
    GeometryPair::t_wedge6>(const Core::Elements::Element*, const Core::Elements::Element*,
    const std::shared_ptr<GeometryEvaluationDataBase>&);


/**
 *
 */
template <typename ScalarType, typename Line, typename Surface>
std::shared_ptr<GeometryPair::GeometryPairBase> GeometryPair::geometry_pair_line_to_surface_factory(
    const Core::Elements::Element* element1, const Core::Elements::Element* element2,
    const std::shared_ptr<GeometryEvaluationDataBase>& geometry_evaluation_data)
{
  // Cast the geometry evaluation data to the correct format.
  auto line_to_surface_evaluation_data =
      std::dynamic_pointer_cast<LineToSurfaceEvaluationData>(geometry_evaluation_data);

  // Get the strategy for line to volume interaction.
  GeometryPair::LineTo3DStrategy strategy = line_to_surface_evaluation_data->get_strategy();

  // Create the class depending on the strategy.
  switch (strategy)
  {
    case GeometryPair::LineTo3DStrategy::gauss_point_projection_without_boundary_segmentation:
    case GeometryPair::LineTo3DStrategy::gauss_point_projection_boundary_segmentation:
      return std::make_shared<
          GeometryPairLineToSurfaceGaussPointProjection<ScalarType, Line, Surface>>(
          element1, element2, line_to_surface_evaluation_data);
    case GeometryPair::LineTo3DStrategy::segmentation:
      return std::make_shared<GeometryPairLineToSurfaceSegmentation<ScalarType, Line, Surface>>(
          element1, element2, line_to_surface_evaluation_data);
    default:
    {
      FOUR_C_THROW("The given geometry pair strategy is not valid.");
      return nullptr;
    }
  }
}

/**
 *
 */
template <typename ScalarType, typename Line, typename Surface>
std::shared_ptr<GeometryPair::GeometryPairBase>
GeometryPair::geometry_pair_line_to_surface_factory_fad(const Core::Elements::Element* element1,
    const Core::Elements::Element* element2,
    const std::shared_ptr<GeometryEvaluationDataBase>& geometry_evaluation_data)
{
  // Create the internal pair as double.
  auto internal_geometry_pair_double =
      std::dynamic_pointer_cast<GeometryPairLineToSurface<double, Line, Surface>>(
          geometry_pair_line_to_surface_factory<double, Line, Surface>(
              element1, element2, geometry_evaluation_data));

  // Create the wrapper.
  return std::make_shared<GeometryPairLineToSurfaceFADWrapper<ScalarType, Line, Surface>>(
      element1, element2, internal_geometry_pair_double);
}


/**
 * Explicit template initialization of factory function.
 */
template std::shared_ptr<GeometryPair::GeometryPairBase>
GeometryPair::geometry_pair_line_to_surface_factory<double, GeometryPair::t_line2,
    GeometryPair::t_tri3>(const Core::Elements::Element*, const Core::Elements::Element*,
    const std::shared_ptr<GeometryEvaluationDataBase>&);
template std::shared_ptr<GeometryPair::GeometryPairBase>
GeometryPair::geometry_pair_line_to_surface_factory<double, GeometryPair::t_line2,
    GeometryPair::t_tri6>(const Core::Elements::Element*, const Core::Elements::Element*,
    const std::shared_ptr<GeometryEvaluationDataBase>&);
template std::shared_ptr<GeometryPair::GeometryPairBase>
GeometryPair::geometry_pair_line_to_surface_factory<double, GeometryPair::t_line2,
    GeometryPair::t_quad4>(const Core::Elements::Element*, const Core::Elements::Element*,
    const std::shared_ptr<GeometryEvaluationDataBase>&);
template std::shared_ptr<GeometryPair::GeometryPairBase>
GeometryPair::geometry_pair_line_to_surface_factory<double, GeometryPair::t_line2,
    GeometryPair::t_quad8>(const Core::Elements::Element*, const Core::Elements::Element*,
    const std::shared_ptr<GeometryEvaluationDataBase>&);
template std::shared_ptr<GeometryPair::GeometryPairBase>
GeometryPair::geometry_pair_line_to_surface_factory<double, GeometryPair::t_line2,
    GeometryPair::t_quad9>(const Core::Elements::Element*, const Core::Elements::Element*,
    const std::shared_ptr<GeometryEvaluationDataBase>&);
template std::shared_ptr<GeometryPair::GeometryPairBase>
GeometryPair::geometry_pair_line_to_surface_factory<double, GeometryPair::t_line2,
    GeometryPair::t_nurbs9>(const Core::Elements::Element*, const Core::Elements::Element*,
    const std::shared_ptr<GeometryEvaluationDataBase>&);

template std::shared_ptr<GeometryPair::GeometryPairBase>
GeometryPair::geometry_pair_line_to_surface_factory_fad<
    GeometryPair::line_to_surface_patch_scalar_type_1st_order, GeometryPair::t_line2,
    GeometryPair::t_tri3>(const Core::Elements::Element*, const Core::Elements::Element*,
    const std::shared_ptr<GeometryEvaluationDataBase>&);
template std::shared_ptr<GeometryPair::GeometryPairBase>
GeometryPair::geometry_pair_line_to_surface_factory_fad<
    GeometryPair::line_to_surface_patch_scalar_type_1st_order, GeometryPair::t_line2,
    GeometryPair::t_tri6>(const Core::Elements::Element*, const Core::Elements::Element*,
    const std::shared_ptr<GeometryEvaluationDataBase>&);
template std::shared_ptr<GeometryPair::GeometryPairBase>
GeometryPair::geometry_pair_line_to_surface_factory_fad<
    GeometryPair::line_to_surface_patch_scalar_type_1st_order, GeometryPair::t_line2,
    GeometryPair::t_quad4>(const Core::Elements::Element*, const Core::Elements::Element*,
    const std::shared_ptr<GeometryEvaluationDataBase>&);
template std::shared_ptr<GeometryPair::GeometryPairBase>
GeometryPair::geometry_pair_line_to_surface_factory_fad<
    GeometryPair::line_to_surface_patch_scalar_type_1st_order, GeometryPair::t_line2,
    GeometryPair::t_quad8>(const Core::Elements::Element*, const Core::Elements::Element*,
    const std::shared_ptr<GeometryEvaluationDataBase>&);
template std::shared_ptr<GeometryPair::GeometryPairBase>
GeometryPair::geometry_pair_line_to_surface_factory_fad<
    GeometryPair::line_to_surface_patch_scalar_type_1st_order, GeometryPair::t_line2,
    GeometryPair::t_quad9>(const Core::Elements::Element*, const Core::Elements::Element*,
    const std::shared_ptr<GeometryEvaluationDataBase>&);
template std::shared_ptr<GeometryPair::GeometryPairBase>
GeometryPair::geometry_pair_line_to_surface_factory_fad<
    GeometryPair::line_to_surface_patch_scalar_type_fixed_size_1st_order<GeometryPair::t_line2,
        GeometryPair::t_nurbs9>,
    GeometryPair::t_line2, GeometryPair::t_nurbs9>(const Core::Elements::Element*,
    const Core::Elements::Element*, const std::shared_ptr<GeometryEvaluationDataBase>&);

template std::shared_ptr<GeometryPair::GeometryPairBase>
GeometryPair::geometry_pair_line_to_surface_factory_fad<
    GeometryPair::line_to_surface_patch_scalar_type, GeometryPair::t_line2, GeometryPair::t_tri3>(
    const Core::Elements::Element*, const Core::Elements::Element*,
    const std::shared_ptr<GeometryEvaluationDataBase>&);
template std::shared_ptr<GeometryPair::GeometryPairBase>
GeometryPair::geometry_pair_line_to_surface_factory_fad<
    GeometryPair::line_to_surface_patch_scalar_type, GeometryPair::t_line2, GeometryPair::t_tri6>(
    const Core::Elements::Element*, const Core::Elements::Element*,
    const std::shared_ptr<GeometryEvaluationDataBase>&);
template std::shared_ptr<GeometryPair::GeometryPairBase>
GeometryPair::geometry_pair_line_to_surface_factory_fad<
    GeometryPair::line_to_surface_patch_scalar_type, GeometryPair::t_line2, GeometryPair::t_quad4>(
    const Core::Elements::Element*, const Core::Elements::Element*,
    const std::shared_ptr<GeometryEvaluationDataBase>&);
template std::shared_ptr<GeometryPair::GeometryPairBase>
GeometryPair::geometry_pair_line_to_surface_factory_fad<
    GeometryPair::line_to_surface_patch_scalar_type, GeometryPair::t_line2, GeometryPair::t_quad8>(
    const Core::Elements::Element*, const Core::Elements::Element*,
    const std::shared_ptr<GeometryEvaluationDataBase>&);
template std::shared_ptr<GeometryPair::GeometryPairBase>
GeometryPair::geometry_pair_line_to_surface_factory_fad<
    GeometryPair::line_to_surface_patch_scalar_type, GeometryPair::t_line2, GeometryPair::t_quad9>(
    const Core::Elements::Element*, const Core::Elements::Element*,
    const std::shared_ptr<GeometryEvaluationDataBase>&);
template std::shared_ptr<GeometryPair::GeometryPairBase>
GeometryPair::geometry_pair_line_to_surface_factory_fad<
    GeometryPair::line_to_surface_patch_scalar_type_fixed_size<GeometryPair::t_line2,
        GeometryPair::t_nurbs9>,
    GeometryPair::t_line2, GeometryPair::t_nurbs9>(const Core::Elements::Element*,
    const Core::Elements::Element*, const std::shared_ptr<GeometryEvaluationDataBase>&);

template std::shared_ptr<GeometryPair::GeometryPairBase>
GeometryPair::geometry_pair_line_to_surface_factory<double, GeometryPair::t_hermite,
    GeometryPair::t_tri3>(const Core::Elements::Element*, const Core::Elements::Element*,
    const std::shared_ptr<GeometryEvaluationDataBase>&);
template std::shared_ptr<GeometryPair::GeometryPairBase>
GeometryPair::geometry_pair_line_to_surface_factory<double, GeometryPair::t_hermite,
    GeometryPair::t_tri6>(const Core::Elements::Element*, const Core::Elements::Element*,
    const std::shared_ptr<GeometryEvaluationDataBase>&);
template std::shared_ptr<GeometryPair::GeometryPairBase>
GeometryPair::geometry_pair_line_to_surface_factory<double, GeometryPair::t_hermite,
    GeometryPair::t_quad4>(const Core::Elements::Element*, const Core::Elements::Element*,
    const std::shared_ptr<GeometryEvaluationDataBase>&);
template std::shared_ptr<GeometryPair::GeometryPairBase>
GeometryPair::geometry_pair_line_to_surface_factory<double, GeometryPair::t_hermite,
    GeometryPair::t_quad8>(const Core::Elements::Element*, const Core::Elements::Element*,
    const std::shared_ptr<GeometryEvaluationDataBase>&);
template std::shared_ptr<GeometryPair::GeometryPairBase>
GeometryPair::geometry_pair_line_to_surface_factory<double, GeometryPair::t_hermite,
    GeometryPair::t_quad9>(const Core::Elements::Element*, const Core::Elements::Element*,
    const std::shared_ptr<GeometryEvaluationDataBase>&);
template std::shared_ptr<GeometryPair::GeometryPairBase>
GeometryPair::geometry_pair_line_to_surface_factory<double, GeometryPair::t_hermite,
    GeometryPair::t_nurbs9>(const Core::Elements::Element*, const Core::Elements::Element*,
    const std::shared_ptr<GeometryEvaluationDataBase>&);

template std::shared_ptr<GeometryPair::GeometryPairBase>
GeometryPair::geometry_pair_line_to_surface_factory_fad<
    GeometryPair::line_to_surface_patch_scalar_type_1st_order, GeometryPair::t_hermite,
    GeometryPair::t_tri3>(const Core::Elements::Element*, const Core::Elements::Element*,
    const std::shared_ptr<GeometryEvaluationDataBase>&);
template std::shared_ptr<GeometryPair::GeometryPairBase>
GeometryPair::geometry_pair_line_to_surface_factory_fad<
    GeometryPair::line_to_surface_patch_scalar_type_1st_order, GeometryPair::t_hermite,
    GeometryPair::t_tri6>(const Core::Elements::Element*, const Core::Elements::Element*,
    const std::shared_ptr<GeometryEvaluationDataBase>&);
template std::shared_ptr<GeometryPair::GeometryPairBase>
GeometryPair::geometry_pair_line_to_surface_factory_fad<
    GeometryPair::line_to_surface_patch_scalar_type_1st_order, GeometryPair::t_hermite,
    GeometryPair::t_quad4>(const Core::Elements::Element*, const Core::Elements::Element*,
    const std::shared_ptr<GeometryEvaluationDataBase>&);
template std::shared_ptr<GeometryPair::GeometryPairBase>
GeometryPair::geometry_pair_line_to_surface_factory_fad<
    GeometryPair::line_to_surface_patch_scalar_type_1st_order, GeometryPair::t_hermite,
    GeometryPair::t_quad8>(const Core::Elements::Element*, const Core::Elements::Element*,
    const std::shared_ptr<GeometryEvaluationDataBase>&);
template std::shared_ptr<GeometryPair::GeometryPairBase>
GeometryPair::geometry_pair_line_to_surface_factory_fad<
    GeometryPair::line_to_surface_patch_scalar_type_1st_order, GeometryPair::t_hermite,
    GeometryPair::t_quad9>(const Core::Elements::Element*, const Core::Elements::Element*,
    const std::shared_ptr<GeometryEvaluationDataBase>&);
template std::shared_ptr<GeometryPair::GeometryPairBase>
GeometryPair::geometry_pair_line_to_surface_factory_fad<
    GeometryPair::line_to_surface_patch_scalar_type_fixed_size_1st_order<GeometryPair::t_hermite,
        GeometryPair::t_nurbs9>,
    GeometryPair::t_hermite, GeometryPair::t_nurbs9>(const Core::Elements::Element*,
    const Core::Elements::Element*, const std::shared_ptr<GeometryEvaluationDataBase>&);

template std::shared_ptr<GeometryPair::GeometryPairBase>
GeometryPair::geometry_pair_line_to_surface_factory_fad<
    GeometryPair::line_to_surface_patch_scalar_type, GeometryPair::t_hermite, GeometryPair::t_tri3>(
    const Core::Elements::Element*, const Core::Elements::Element*,
    const std::shared_ptr<GeometryEvaluationDataBase>&);
template std::shared_ptr<GeometryPair::GeometryPairBase>
GeometryPair::geometry_pair_line_to_surface_factory_fad<
    GeometryPair::line_to_surface_patch_scalar_type, GeometryPair::t_hermite, GeometryPair::t_tri6>(
    const Core::Elements::Element*, const Core::Elements::Element*,
    const std::shared_ptr<GeometryEvaluationDataBase>&);
template std::shared_ptr<GeometryPair::GeometryPairBase> GeometryPair::
    geometry_pair_line_to_surface_factory_fad<GeometryPair::line_to_surface_patch_scalar_type,
        GeometryPair::t_hermite, GeometryPair::t_quad4>(const Core::Elements::Element*,
        const Core::Elements::Element*, const std::shared_ptr<GeometryEvaluationDataBase>&);
template std::shared_ptr<GeometryPair::GeometryPairBase> GeometryPair::
    geometry_pair_line_to_surface_factory_fad<GeometryPair::line_to_surface_patch_scalar_type,
        GeometryPair::t_hermite, GeometryPair::t_quad8>(const Core::Elements::Element*,
        const Core::Elements::Element*, const std::shared_ptr<GeometryEvaluationDataBase>&);
template std::shared_ptr<GeometryPair::GeometryPairBase> GeometryPair::
    geometry_pair_line_to_surface_factory_fad<GeometryPair::line_to_surface_patch_scalar_type,
        GeometryPair::t_hermite, GeometryPair::t_quad9>(const Core::Elements::Element*,
        const Core::Elements::Element*, const std::shared_ptr<GeometryEvaluationDataBase>&);
template std::shared_ptr<GeometryPair::GeometryPairBase>
GeometryPair::geometry_pair_line_to_surface_factory_fad<
    GeometryPair::line_to_surface_patch_scalar_type_fixed_size<GeometryPair::t_hermite,
        GeometryPair::t_nurbs9>,
    GeometryPair::t_hermite, GeometryPair::t_nurbs9>(const Core::Elements::Element*,
    const Core::Elements::Element*, const std::shared_ptr<GeometryEvaluationDataBase>&);

FOUR_C_NAMESPACE_CLOSE
