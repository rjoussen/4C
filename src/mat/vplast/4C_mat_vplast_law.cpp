// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_mat_vplast_law.hpp"

#include "4C_global_data.hpp"
#include "4C_legacy_enum_definitions_materials.hpp"
#include "4C_mat_par_bundle.hpp"
#include "4C_mat_vplast_anand.hpp"
#include "4C_mat_vplast_reform_johnsoncook.hpp"

#include <utility>

FOUR_C_NAMESPACE_OPEN


/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
Mat::Viscoplastic::Law::Law(Core::Mat::PAR::Parameter* params) : params_(params) {}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
Mat::Viscoplastic::Law::Law() : params_(nullptr) {}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
std::shared_ptr<Mat::Viscoplastic::Law> Mat::Viscoplastic::Law::factory(int matnum)
{
  // for the sake of safety
  if (Global::Problem::instance()->materials() == nullptr)
    FOUR_C_THROW("List of materials cannot be accessed in the global problem instance.");

  // another safety check
  if (Global::Problem::instance()->materials()->num() == 0)
    FOUR_C_THROW("List of materials in the global problem instance is empty.");

  // retrieve problem instance to read from
  const int probinst = Global::Problem::instance()->materials()->get_read_from_problem();
  // retrieve validated input line of material ID in question
  auto* curmat = Global::Problem::instance(probinst)->materials()->parameter_by_id(matnum);


  // get material type and call corresponding constructors
  const Core::Materials::MaterialType currentMaterialType = curmat->type();
  switch (currentMaterialType)
  {
    case Core::Materials::mvl_reformulated_Johnson_Cook:
    {
      // get pointer to parameter class
      auto* params = dynamic_cast<Mat::Viscoplastic::PAR::ReformulatedJohnsonCook*>(curmat);

      // return pointer to material
      return std::make_shared<Mat::Viscoplastic::ReformulatedJohnsonCook>(params);
    }

    case Core::Materials::mvl_Anand:
    {
      // get pointer to parameter class
      auto* params = dynamic_cast<Mat::Viscoplastic::PAR::Anand*>(curmat);

      // return pointer to material
      return std::make_shared<Mat::Viscoplastic::Anand>(params);
    }


    default:
      FOUR_C_THROW("cannot deal with type %d", curmat->type());
  }
}

FOUR_C_NAMESPACE_CLOSE
