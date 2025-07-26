// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_geometry_pair_line_to_surface_evaluation_data.hpp"

#include "4C_geometry_pair_element_faces.hpp"
#include "4C_geometry_pair_utility_classes.hpp"
#include "4C_utils_parameter_list.hpp"

FOUR_C_NAMESPACE_OPEN

/**
 *
 */
GeometryPair::LineToSurfaceEvaluationData::LineToSurfaceEvaluationData(
    const Teuchos::ParameterList& input_parameter_list)
    : LineTo3DEvaluationData(input_parameter_list),
      face_elements_(),
      surface_normal_strategy_(GeometryPair::SurfaceNormals::standard)
{
  surface_normal_strategy_ = Teuchos::getIntegralValue<GeometryPair::SurfaceNormals>(
      input_parameter_list, "GEOMETRY_PAIR_SURFACE_NORMALS");
}

/**
 *
 */
void GeometryPair::LineToSurfaceEvaluationData::clear()
{
  // Call reset on the base method.
  LineTo3DEvaluationData::clear();
  face_elements_.clear();
}

/**
 *
 */
void GeometryPair::LineToSurfaceEvaluationData::setup(
    const std::shared_ptr<const Core::FE::Discretization>& discret,
    const std::unordered_map<int, std::shared_ptr<GeometryPair::FaceElement>>& face_elements)
{
  face_elements_ = face_elements;

  for (const auto& face_element_iterator : face_elements_)
    face_element_iterator.second->setup(discret, face_elements_);

  // The averaged reference normals have to be calculated after each face element is set up.
  for (const auto& face_element_iterator : face_elements_)
    if (face_element_iterator.second->is_part_of_pair())
      face_element_iterator.second->calculate_averaged_reference_normals(face_elements_);
}

/**
 *
 */
void GeometryPair::LineToSurfaceEvaluationData::set_state(
    const std::shared_ptr<const Core::LinAlg::Vector<double>>& displacement_col_np)
{
  for (const auto& [id, face_element] : face_elements_)
    if (face_element->is_part_of_pair())
      face_element->set_state(displacement_col_np, face_elements_);
}

FOUR_C_NAMESPACE_CLOSE
