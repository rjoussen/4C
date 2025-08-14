// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_mat_par_bundle.hpp"

#include "4C_comm_pack_helpers.hpp"
#include "4C_mat_elast_summand.hpp"
#include "4C_mat_material_factory.hpp"
#include "4C_material_base.hpp"
#include "4C_utils_enum.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
void Mat::PAR::Bundle::insert(int matid, Core::Utils::LazyPtr<Core::Mat::PAR::Parameter> mat)
{
  if (matmap_.contains(matid)) FOUR_C_THROW("Material with ID {} already exists", matid);
  matmap_.emplace(matid, mat);
}

/*----------------------------------------------------------------------*/
bool Mat::PAR::Bundle::id_exists(int id) const { return !(matmap_.find(id) == matmap_.end()); }

/*----------------------------------------------------------------------*/
Core::Mat::PAR::Parameter* Mat::PAR::Bundle::parameter_by_id(const int num) const
{
  if (matmap_.size() == 0) FOUR_C_THROW("No materials available, num={}", num);

  if (auto it = matmap_.find(num); it != matmap_.end())
    return it->second.get();
  else
    FOUR_C_THROW("Material 'MAT {}' could not be found", num);
}

/*----------------------------------------------------------------------*/
int Mat::PAR::Bundle::first_id_by_type(const Core::Materials::MaterialType type) const
{
  if (auto it = std::find_if(matmap_.begin(), matmap_.end(),
          [&type](const auto& m) { return m.second->type() == type; });
      it != matmap_.end())
    return it->first;
  else
    return -1;
}

FOUR_C_NAMESPACE_CLOSE
