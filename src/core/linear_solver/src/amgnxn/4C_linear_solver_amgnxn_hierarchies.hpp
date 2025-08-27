// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_LINEAR_SOLVER_AMGNXN_HIERARCHIES_HPP
#define FOUR_C_LINEAR_SOLVER_AMGNXN_HIERARCHIES_HPP

#include "4C_config.hpp"

#include "4C_linalg_blocksparsematrix.hpp"
#include "4C_linear_solver_amgnxn_smoothers.hpp"
#include "4C_linear_solver_method_linalg.hpp"
#include "4C_linear_solver_preconditioner_type.hpp"

#include <Epetra_Operator.h>
#include <MueLu.hpp>
#include <MueLu_BaseClass.hpp>
#include <MueLu_Level.hpp>
#include <MueLu_UseDefaultTypes.hpp>
#include <MueLu_Utilities.hpp>
#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN

namespace Core::LinearSolver::AMGNxN
{
  class Hierarchies
  {
   public:
    Hierarchies(Teuchos::RCP<AMGNxN::BlockedMatrix> A,
        std::vector<Teuchos::ParameterList> muelu_params, std::vector<int> num_pdes,
        std::vector<int> null_spaces_dim,
        std::vector<std::shared_ptr<std::vector<double>>> null_spaces_data, int NumLevelAMG,
        std::string verbosity = "off");

    int get_num_level_min();
    int get_num_blocks();
    int get_num_levels(int block);

    Teuchos::RCP<AMGNxN::BlockedMatrix> get_block_matrix();
    Teuchos::RCP<MueLu::Hierarchy<Scalar, LocalOrdinal, GlobalOrdinal, Node>> get_h(int block);

    Teuchos::RCP<Core::LinAlg::SparseMatrix> get_a(int block, int level);
    Teuchos::RCP<Core::LinAlg::SparseMatrix> get_p(int block, int level);
    Teuchos::RCP<Core::LinAlg::SparseMatrix> get_r(int block, int level);
    Teuchos::RCP<AMGNxN::MueluSmootherWrapper> get_s_pre(int block, int level);
    Teuchos::RCP<AMGNxN::MueluSmootherWrapper> get_s_pos(int block, int level);

    std::vector<Teuchos::RCP<Core::LinAlg::SparseMatrix>> get_a(int block);
    std::vector<Teuchos::RCP<Core::LinAlg::SparseMatrix>> get_p(int block);
    std::vector<Teuchos::RCP<Core::LinAlg::SparseMatrix>> get_r(int block);
    std::vector<Teuchos::RCP<AMGNxN::MueluSmootherWrapper>> get_s_pre(int block);
    std::vector<Teuchos::RCP<AMGNxN::MueluSmootherWrapper>> get_s_pos(int block);

    int get_num_pd_es(int block);
    int get_null_space_dim(int block);
    std::shared_ptr<std::vector<double>> get_null_space_data(int block);

   private:
    Teuchos::RCP<AMGNxN::BlockedMatrix> a_;
    std::vector<Teuchos::ParameterList> muelu_params_;
    std::vector<int> num_pdes_;
    std::vector<int> null_spaces_dim_;
    std::vector<std::shared_ptr<std::vector<double>>> null_spaces_data_;
    int num_blocks_;
    int num_level_max_;
    int num_level_min_;
    int num_level_amg_;
    std::vector<Teuchos::RCP<MueLu::Hierarchy<Scalar, LocalOrdinal, GlobalOrdinal, Node>>> h_block_;
    std::vector<std::vector<Teuchos::RCP<Core::LinAlg::SparseMatrix>>> a_block_level_;
    std::vector<std::vector<Teuchos::RCP<Core::LinAlg::SparseMatrix>>> p_block_level_;
    std::vector<std::vector<Teuchos::RCP<Core::LinAlg::SparseMatrix>>> r_block_level_;
    std::vector<std::vector<Teuchos::RCP<AMGNxN::MueluSmootherWrapper>>> s_pre_block_level_;
    std::vector<std::vector<Teuchos::RCP<AMGNxN::MueluSmootherWrapper>>> s_pos_block_level_;
    std::string verbosity_;

    void setup();

    Teuchos::RCP<MueLu::Hierarchy<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
    build_mue_lu_hierarchy(Teuchos::ParameterList paramListFromXml, int numdf, int dimns,
        std::shared_ptr<std::vector<double>> nsdata, Teuchos::RCP<Epetra_CrsMatrix> A_eop,
        int block, int NumBlocks, std::vector<int>& offsets, int offsetFineLevel);

    std::string convert_int(int number)
    {
      std::stringstream ss;
      ss << number;
      return ss.str();
    }
  };

  class Vcycle;

  class MonolithicHierarchy
  {
   public:
    MonolithicHierarchy(Teuchos::RCP<AMGNxN::Hierarchies> H, const Teuchos::ParameterList& params,
        const Teuchos::ParameterList& params_smoothers);

    int get_num_levels();

    Teuchos::RCP<AMGNxN::Hierarchies> get_hierarchies();

    Teuchos::RCP<AMGNxN::BlockedMatrix> get_a(int level);

    Teuchos::RCP<Vcycle> build_v_cycle();

   private:
    Teuchos::RCP<AMGNxN::Hierarchies> h_;
    int num_levels_;
    int num_blocks_;
    std::vector<Teuchos::RCP<AMGNxN::BlockedMatrix>> a_;
    std::vector<Teuchos::RCP<AMGNxN::BlockedMatrix>> p_;
    std::vector<Teuchos::RCP<AMGNxN::BlockedMatrix>> r_;
    std::vector<Teuchos::RCP<AMGNxN::GenericSmoother>> spre_;
    std::vector<Teuchos::RCP<AMGNxN::GenericSmoother>> spos_;
    Teuchos::ParameterList params_;
    Teuchos::ParameterList params_smoothers_;

    void setup();
    Teuchos::RCP<AMGNxN::GenericSmoother> build_smoother(int level);
  };
}  // namespace Core::LinearSolver::AMGNxN

FOUR_C_NAMESPACE_CLOSE

#endif
