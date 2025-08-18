// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_contact_abstract_data_container.hpp"

FOUR_C_NAMESPACE_OPEN

CONTACT::AbstractStrategyDataContainer::AbstractStrategyDataContainer()
    : glmdofrowmap_(nullptr),
      gsnoderowmap_(nullptr),
      gmnoderowmap_(nullptr),
      gsdofrowmap_(nullptr),
      gmdofrowmap_(nullptr),
      gndofrowmap_(nullptr),
      gsmdofrowmap_(nullptr),
      gdisprowmap_(nullptr),
      gactivenodes_(nullptr),
      gactivedofs_(nullptr),
      ginactivenodes_(nullptr),
      ginactivedofs_(nullptr),
      gactiven_(nullptr),
      gactivet_(nullptr),
      gslipnodes_(nullptr),
      gslipdofs_(nullptr),
      gslipt_(nullptr),
      gsdof_vertex_(nullptr),
      gsdof_edge_(nullptr),
      gsdof_surf_(nullptr),
      unbalance_evaluation_time_(0),
      unbalance_num_slave_elements_(0),
      non_redist_glmdofrowmap_(nullptr),
      non_redist_gsdofrowmap_(nullptr),
      non_redist_gmdofrowmap_(nullptr),
      non_redist_gsmdofrowmap_(nullptr),
      non_redist_gsdirichtoggle_(nullptr),
      partype_(Inpar::Mortar::ParallelRedist::redist_none),
      dmatrix_(nullptr),
      mmatrix_(nullptr),
      wgap_(nullptr),
      tangrhs_(nullptr),
      inactiverhs_(nullptr),
      str_contact_rhs_ptr_(nullptr),
      constrrhs_(nullptr),
      lindmatrix_(nullptr),
      linmmatrix_(nullptr),
      kteffnew_(nullptr),
      dold_(nullptr),
      mold_(nullptr),
      z_(nullptr),
      zold_(nullptr),
      zincr_(nullptr),
      zuzawa_(nullptr),
      stressnormal_(nullptr),
      stresstangential_(nullptr),
      forcenormal_(nullptr),
      forcetangential_(nullptr),
      stepnp_(-1),
      iter_(-1),
      isincontact_(false),
      wasincontact_(false),
      wasincontactlts_(false),
      isselfcontact_(false),
      friction_(false),
      non_smooth_contact_(false),
      dualquadslavetrafo_(false),
      trafo_(nullptr),
      invtrafo_(nullptr),
      dmatrixmod_(nullptr),
      doldmod_(nullptr),
      inttime_(0.0),
      ivel_(0),
      stype_(CONTACT::SolvingStrategy::vague),
      constr_direction_(CONTACT::ConstraintDirection::vague)
{
  return;
}

FOUR_C_NAMESPACE_CLOSE
