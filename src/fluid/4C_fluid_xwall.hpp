// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_FLUID_XWALL_HPP
#define FOUR_C_FLUID_XWALL_HPP

#include "4C_config.hpp"

#include "4C_fem_discretization.hpp"
#include "4C_inpar_fluid.hpp"
#include "4C_linalg_utils_sparse_algebra_manipulation.hpp"
#include "4C_linalg_vector.hpp"
#include "4C_utils_parameter_list.fwd.hpp"

#include <memory>

FOUR_C_NAMESPACE_OPEN

// forward declarations
namespace Core::FE
{
  class Discretization;
}  // namespace Core::FE

namespace Core::LinAlg
{
  class SparseMatrix;
  class Solver;
}  // namespace Core::LinAlg

namespace Core::IO
{
  class DiscretizationReader;
}

namespace FLD
{
  class FluidImplicitTimeInt;
  class TransferTurbulentInflowConditionNodal;
  namespace Utils
  {
    class StressManager;
  }

  class XWall
  {
   public:
    /// Standard Constructor
    XWall(std::shared_ptr<Core::FE::Discretization> dis, int nsd,
        std::shared_ptr<Teuchos::ParameterList>& params,
        std::shared_ptr<Core::LinAlg::MapExtractor> dbcmaps,
        std::shared_ptr<FLD::Utils::StressManager> wssmanager);

    /// Destructor
    virtual ~XWall() = default;

    // set element params for xwall EnrichmentType
    virtual void set_x_wall_params(Teuchos::ParameterList& eleparams);

    // adapt ml nullspace for aggregation for scale separation (mfs)
    void adapt_ml_nullspace(Core::LinAlg::Solver& solver);

    // get output vector of enriched dofs
    std::shared_ptr<Core::LinAlg::Vector<double>> get_output_vector(
        Core::LinAlg::Vector<double>& vel);

    // returns, if Properties for GenAlpha have to be updated
    virtual void update_tau_w(int step, std::shared_ptr<Core::LinAlg::Vector<double>> trueresidual,
        int itnum, std::shared_ptr<Core::LinAlg::Vector<double>> accn,
        std::shared_ptr<Core::LinAlg::Vector<double>> velnp,
        std::shared_ptr<Core::LinAlg::Vector<double>> veln);

    // returns tauw of discret_
    std::shared_ptr<Core::LinAlg::Vector<double>> get_tauw();

    std::shared_ptr<Core::LinAlg::Vector<double>> get_tauw_vector()
    {
      std::shared_ptr<Core::LinAlg::Vector<double>> tauw =
          std::make_shared<Core::LinAlg::Vector<double>>(*(discret_->node_row_map()), true);
      Core::LinAlg::export_to(*tauw_, *tauw);
      return tauw;
    }

    // read restart including wall stresses
    void read_restart(Core::IO::DiscretizationReader& reader);

    // fix residual at Dirichlet-inflow nodes such that the wss can be calculated
    std::shared_ptr<Core::LinAlg::Vector<double>> fix_dirichlet_inflow(
        Core::LinAlg::Vector<double>& trueresidual);

    // set current non-linear iteration number
    void set_iter(int iter)
    {
      iter_ = iter;
      return;
    }

   protected:
    // set element params for xwall EnrichmentType, distributed for xwall discretization
    virtual void set_x_wall_params_xw_dis(Teuchos::ParameterList& eleparams);

    // setup XWall
    void setup();

    // initialize some maps
    void init_x_wall_maps();

    // initialize a toggle vector for the element
    void init_toggle_vector();

    // initialize wall distance
    void init_wall_dist();

    // setup xwall discretization
    void setup_x_wall_dis();

    // setup l2 projection
    void setup_l2_projection();

    // calculate wall shear stress
    void calc_tau_w(
        int step, Core::LinAlg::Vector<double>& velnp, Core::LinAlg::Vector<double>* wss);

    // l2 project vectors
    void l2_project_vector(Core::LinAlg::Vector<double>& veln,
        std::shared_ptr<Core::LinAlg::Vector<double>> velnp,
        std::shared_ptr<Core::LinAlg::Vector<double>> accn);

    // calculate parameter for stabilization parameter mk
    void calc_mk();

    // for inflowchannel
    void transfer_and_save_tauw();

    // for inflowchannel
    void overwrite_transferred_values();

    //! discretisation
    std::shared_ptr<Core::FE::Discretization> discret_;

    //! fluid params
    std::shared_ptr<Teuchos::ParameterList> params_;

    //! manager for wall shear stress
    std::shared_ptr<FLD::Utils::StressManager> mystressmanager_;

    //! the processor ID from the communicator
    int myrank_;

    //! map including all wall nodes, redundant map
    std::shared_ptr<Core::LinAlg::Map> dircolnodemap_;

    //! xwall node row map
    std::shared_ptr<Core::LinAlg::Map> xwallrownodemap_;

    //! xwall node col map
    std::shared_ptr<Core::LinAlg::Map> xwallcolnodemap_;

    //! map including the enriched dofs, row map
    std::shared_ptr<Core::LinAlg::Map> enrdofrowmap_;

    //! map including the unused pressure dofs, row map
    std::shared_ptr<Core::LinAlg::Map> lagrdofrowmap_;

    //! map including all enriched dofs plus unused pressure dofs
    std::shared_ptr<Core::LinAlg::Map> mergedmap_;

    //! wall distance, local vector
    std::shared_ptr<Core::LinAlg::Vector<double>> walldist_;

    //! wall distance, standard node row map
    std::shared_ptr<Core::LinAlg::Vector<double>> wdist_;

    //! wall distance, row map of redistributed discretization
    std::shared_ptr<Core::LinAlg::Vector<double>> wdistxwdis_;

    //! vector on same dofs for tauw
    std::shared_ptr<Core::LinAlg::Vector<double>> tauw_;

    //! vector on same dofs for tauw
    std::shared_ptr<Core::LinAlg::Vector<double>> tauwxwdis_;

    //! vector on same dofs for deltatauw (increment)
    std::shared_ptr<Core::LinAlg::Vector<double>> inctauw_;

    //! vector on same dofs for deltatauw (increment)
    std::shared_ptr<Core::LinAlg::Vector<double>> inctauwxwdis_;

    //! vector on same dofs for old tauw
    std::shared_ptr<Core::LinAlg::Vector<double>> oldtauw_;

    //! vector on same dofs for old tauw
    std::shared_ptr<Core::LinAlg::Vector<double>> oldinctauw_;

    //! matrix projecting the wall shear stress to off-wall nodes
    std::shared_ptr<Core::LinAlg::SparseMatrix> tauwcouplingmattrans_;

    //! toggle vector, standard node row map
    std::shared_ptr<Core::LinAlg::Vector<double>> xwalltoggle_;

    //! toggle vector, xwall discretization
    std::shared_ptr<Core::LinAlg::Vector<double>> xwalltogglexwdis_;

    //! toggle vector, local vector
    std::shared_ptr<Core::LinAlg::Vector<double>> xtoggleloc_;

    //! redistributed xwall discretization
    std::shared_ptr<Core::FE::Discretization> xwdiscret_;

    //! mass matrix for projection
    std::shared_ptr<Core::LinAlg::SparseMatrix> massmatrix_;

    //! solver for projection
    std::shared_ptr<Core::LinAlg::Solver> solver_;

    //! increment of veln during projection
    std::shared_ptr<Core::LinAlg::Vector<double>> incveln_;

    //! increment of velnp during projection
    std::shared_ptr<Core::LinAlg::Vector<double>> incvelnp_;

    //! increment of accn during projection
    std::shared_ptr<Core::LinAlg::Vector<double>> incaccn_;

    //! veln for state of xwall discretization during projection
    std::shared_ptr<Core::LinAlg::Vector<double>> stateveln_;

    //! velnp for state of xwall discretization during projection
    std::shared_ptr<Core::LinAlg::Vector<double>> statevelnp_;

    //! accn for state of xwall discretization during projection
    std::shared_ptr<Core::LinAlg::Vector<double>> stateaccn_;

    //! MK for standard discretization
    std::shared_ptr<Core::LinAlg::Vector<double>> mkstate_;

    //! MK for xwall discretization
    std::shared_ptr<Core::LinAlg::Vector<double>> mkxwstate_;

    std::shared_ptr<Core::LinAlg::Vector<double>> restart_wss_;

    std::shared_ptr<TransferTurbulentInflowConditionNodal> turbulent_inflow_condition_;

    //! increment factor of tauw
    double fac_;

    //! increment of tauw
    double inctauwnorm_;

    //! constant tauw from input file
    double constant_tauw_;

    //! minimum tauw from input file
    double min_tauw_;

    //! number of gauss points in wall-parallel and normal direction
    int gp_norm_;
    int gp_norm_ow_;
    int gp_par_;

    //! viscosity
    double visc_;

    //! density
    double dens_;

    //! when and how to update tauw
    enum Inpar::FLUID::XWallTauwType tauwtype_;

    //! how to calculate tauw
    enum Inpar::FLUID::XWallTauwCalcType tauwcalctype_;

    //! how to blend
    enum Inpar::FLUID::XWallBlendingType blendingtype_;

    //! projection
    bool proj_;

    //! smoothing through aggregation of residual
    bool smooth_res_aggregation_;

    //! mfs on fine scales
    bool fix_residual_on_inflow_;

    //! switch from gradient-based to residual-based calculation of tauw
    int switch_step_;

    int iter_;

  };  // class XWall

  class XWallAleFSI : public XWall
  {
   public:
    /// Standard Constructor
    XWallAleFSI(std::shared_ptr<Core::FE::Discretization> dis, int nsd,
        std::shared_ptr<Teuchos::ParameterList>& params,
        std::shared_ptr<Core::LinAlg::MapExtractor> dbcmaps,
        std::shared_ptr<FLD::Utils::StressManager> wssmanager,
        std::shared_ptr<Core::LinAlg::Vector<double>> dispnp,
        std::shared_ptr<Core::LinAlg::Vector<double>> gridv);

    void update_w_dist_wale();

    // set element params for xwall EnrichmentType
    void set_x_wall_params(Teuchos::ParameterList& eleparams) override;

    // returns, if Properties for GenAlpha have to be updated
    void update_tau_w(int step, std::shared_ptr<Core::LinAlg::Vector<double>> trueresidual,
        int itnum, std::shared_ptr<Core::LinAlg::Vector<double>> accn,
        std::shared_ptr<Core::LinAlg::Vector<double>> velnp,
        std::shared_ptr<Core::LinAlg::Vector<double>> veln) override;

   private:
    // set element params for xwall EnrichmentType, distributed for xwall discretization
    void set_x_wall_params_xw_dis(Teuchos::ParameterList& eleparams) override;
    std::shared_ptr<Core::LinAlg::Vector<double>> mydispnp_;
    std::shared_ptr<Core::LinAlg::Vector<double>> mygridv_;

    //! wall distance, row map of redistributed discretization
    std::shared_ptr<Core::LinAlg::Vector<double>> incwdistxwdis_;
  };

}  // namespace FLD


FOUR_C_NAMESPACE_CLOSE

#endif
