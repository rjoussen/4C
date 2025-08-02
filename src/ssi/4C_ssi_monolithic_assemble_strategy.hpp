// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_SSI_MONOLITHIC_ASSEMBLE_STRATEGY_HPP
#define FOUR_C_SSI_MONOLITHIC_ASSEMBLE_STRATEGY_HPP

#include "4C_config.hpp"

#include "4C_coupling_adapter.hpp"
#include "4C_linalg_sparseoperator.hpp"
#include "4C_linalg_vector.hpp"
#include "4C_ssi_utils.hpp"
#include "4C_utils_parameter_list.fwd.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Core::LinAlg
{
  class MultiMapExtractor;
  class Solver;
  class SparseMatrix;
}  // namespace Core::LinAlg

namespace SSI
{
  namespace Utils
  {
    class SSISlaveSideConverter;
  }
  class SsiMono;

  /*!
  We have three options how the global system matrix and the sub matrices are arranged:
  1) System matrix: sparse
    ->Scatra matrix sparse
    ->Structure matrix sparse
  2) System matrix: block
    2a) Scatra matrix block
    ->Structure matrix sparse
    2b) Scatra matrix sparse
    ->Structure matrix sparse

  The inheritance hierarchy is appropriate*/
  class AssembleStrategyBase
  {
   public:
    /**
     * Virtual destructor.
     */
    virtual ~AssembleStrategyBase() = default;

    //! constructor
    explicit AssembleStrategyBase(
        std::shared_ptr<const SSI::Utils::SSIMaps> ssi_maps, const bool is_scatra_manifold);

    //! assemble RHS
    void assemble_rhs(std::shared_ptr<Core::LinAlg::Vector<double>> rhs,
        std::shared_ptr<const Core::LinAlg::Vector<double>> rhs_scatra,
        std::shared_ptr<const Core::LinAlg::Vector<double>> rhs_structure,
        const Core::LinAlg::Vector<double>* rhs_manifold);

    //! assemble ScaTra-ScaTra-Block into system matrix
    virtual void assemble_scatra_scatra(std::shared_ptr<Core::LinAlg::SparseOperator> systemmatrix,
        std::shared_ptr<const Core::LinAlg::SparseOperator> scatra_scatra_matrix) = 0;

    //! assemble ScaTra-Structure-Block into system matrix
    virtual void assemble_scatra_structure(
        std::shared_ptr<Core::LinAlg::SparseOperator> systemmatrix,
        std::shared_ptr<const Core::LinAlg::SparseOperator> scatra_structure_matrix) = 0;

    //! assemble Structure-Structure-Block into system matrix
    virtual void assemble_structure_structure(
        std::shared_ptr<Core::LinAlg::SparseOperator> systemmatrix,
        std::shared_ptr<const Core::LinAlg::SparseMatrix> structure_structure_matrix) = 0;

    //! assemble Structure-ScaTra-Block into system matrix
    virtual void assemble_structure_scatra(
        std::shared_ptr<Core::LinAlg::SparseOperator> systemmatrix,
        std::shared_ptr<const Core::LinAlg::SparseOperator> structure_scatra_matrix) = 0;

    //! assemble ScaTra Manifold-ScaTra Manifold-Block into system matrix
    virtual void assemble_scatramanifold_scatramanifold(
        std::shared_ptr<Core::LinAlg::SparseOperator> systemmatrix,
        std::shared_ptr<const Core::LinAlg::SparseOperator>
            scatramanifold_scatramanifold_matrix) = 0;

    //! assemble ScaTra Manifold-Structure-Block into system matrix
    virtual void assemble_scatramanifold_structure(
        std::shared_ptr<Core::LinAlg::SparseOperator> systemmatrix,
        std::shared_ptr<const Core::LinAlg::SparseOperator> scatramanifold_structure_matrix) = 0;

    //! assemble ScaTra Manifold-ScaTra-Block into system matrix
    virtual void assemble_scatramanifold_scatra(
        std::shared_ptr<Core::LinAlg::SparseOperator> systemmatrix,
        std::shared_ptr<const Core::LinAlg::SparseOperator> scatramanifold_scatra_matrix) = 0;

    //! assemble ScaTra-ScaTra Manifold-Block into system matrix
    virtual void assemble_scatra_scatramanifold(
        std::shared_ptr<Core::LinAlg::SparseOperator> systemmatrix,
        std::shared_ptr<const Core::LinAlg::SparseOperator> scatra_scatramanifold_matrix) = 0;

   protected:
    //! solve additional scatra field on manifolds
    bool is_scatra_manifold() const { return is_scatra_manifold_; }

    //! this object holds all maps relevant to monolithic scalar transport - structure interaction
    std::shared_ptr<const SSI::Utils::SSIMaps> ssi_maps() const { return ssi_maps_; }

   private:
    //! solve additional scatra field on manifolds
    const bool is_scatra_manifold_;

    //! this object holds all maps relevant to monolithic scalar transport - structure interaction
    std::shared_ptr<const SSI::Utils::SSIMaps> ssi_maps_;
  };

  //======================================================================================================
  // SSI problem is organized in sub matrices
  class AssembleStrategyBlock : public AssembleStrategyBase
  {
   public:
    explicit AssembleStrategyBlock(
        std::shared_ptr<const SSI::Utils::SSIMaps> ssi_maps, const bool is_scatra_manifold);

    void assemble_scatra_scatra(std::shared_ptr<Core::LinAlg::SparseOperator> systemmatrix,
        std::shared_ptr<const Core::LinAlg::SparseOperator> scatra_scatra_matrix) override = 0;

    void assemble_scatra_structure(std::shared_ptr<Core::LinAlg::SparseOperator> systemmatrix,
        std::shared_ptr<const Core::LinAlg::SparseOperator> scatra_structure_matrix) override = 0;

    void assemble_structure_structure(std::shared_ptr<Core::LinAlg::SparseOperator> systemmatrix,
        std::shared_ptr<const Core::LinAlg::SparseMatrix> structure_structure_matrix) override = 0;

    void assemble_structure_scatra(std::shared_ptr<Core::LinAlg::SparseOperator> systemmatrix,
        std::shared_ptr<const Core::LinAlg::SparseOperator> structure_scatra_matrix) override = 0;

    void assemble_scatramanifold_scatramanifold(
        std::shared_ptr<Core::LinAlg::SparseOperator> systemmatrix,
        std::shared_ptr<const Core::LinAlg::SparseOperator> scatramanifold_scatramanifold_matrix)
        override = 0;

    void assemble_scatramanifold_structure(
        std::shared_ptr<Core::LinAlg::SparseOperator> systemmatrix,
        std::shared_ptr<const Core::LinAlg::SparseOperator> scatramanifold_structure_matrix)
        override = 0;

    void assemble_scatramanifold_scatra(std::shared_ptr<Core::LinAlg::SparseOperator> systemmatrix,
        std::shared_ptr<const Core::LinAlg::SparseOperator> scatramanifold_scatra_matrix) override =
        0;

    void assemble_scatra_scatramanifold(std::shared_ptr<Core::LinAlg::SparseOperator> systemmatrix,
        std::shared_ptr<const Core::LinAlg::SparseOperator> scatra_scatramanifold_matrix) override =
        0;

   protected:
    //! position of scatra blocks in system matrix
    const std::vector<int>& block_position_scatra() const { return block_position_scatra_; }

    //! position of scatra manifold blocks in system matrix
    const std::vector<int>& block_position_scatra_manifold() const
    {
      return block_position_scatra_manifold_;
    }

    //! position of structure block in system matrix
    int position_structure() const { return position_structure_; };

   private:
    //! position of scatra blocks in system matrix
    const std::vector<int> block_position_scatra_;

    //! position of scatra manifold blocks in system matrix
    std::vector<int> block_position_scatra_manifold_;

    //! position of structure block in system matrix
    const int position_structure_;
  };

  // *********************************************************************************************
  // SSI problem is organized in sparse structure sub matrix and block scatra sub matrix
  class AssembleStrategyBlockBlock : public AssembleStrategyBlock
  {
   public:
    explicit AssembleStrategyBlockBlock(
        std::shared_ptr<const SSI::Utils::SSIMaps> ssi_maps, const bool is_scatra_manifold);

    void assemble_scatra_scatra(std::shared_ptr<Core::LinAlg::SparseOperator> systemmatrix,
        std::shared_ptr<const Core::LinAlg::SparseOperator> scatra_scatra_matrix) override;

    void assemble_scatra_structure(std::shared_ptr<Core::LinAlg::SparseOperator> systemmatrix,
        std::shared_ptr<const Core::LinAlg::SparseOperator> scatra_structure_matrix) override;

    void assemble_structure_structure(std::shared_ptr<Core::LinAlg::SparseOperator> systemmatrix,
        std::shared_ptr<const Core::LinAlg::SparseMatrix> structure_structure_matrix) override;

    void assemble_structure_scatra(std::shared_ptr<Core::LinAlg::SparseOperator> systemmatrix,
        std::shared_ptr<const Core::LinAlg::SparseOperator> structure_scatra_matrix) override;

    void assemble_scatramanifold_scatramanifold(
        std::shared_ptr<Core::LinAlg::SparseOperator> systemmatrix,
        std::shared_ptr<const Core::LinAlg::SparseOperator> scatramanifold_scatramanifold_matrix)
        override;

    void assemble_scatramanifold_structure(
        std::shared_ptr<Core::LinAlg::SparseOperator> systemmatrix,
        std::shared_ptr<const Core::LinAlg::SparseOperator> scatramanifold_structure_matrix)
        override;

    void assemble_scatramanifold_scatra(std::shared_ptr<Core::LinAlg::SparseOperator> systemmatrix,
        std::shared_ptr<const Core::LinAlg::SparseOperator> scatramanifold_scatra_matrix) override;

    void assemble_scatra_scatramanifold(std::shared_ptr<Core::LinAlg::SparseOperator> systemmatrix,
        std::shared_ptr<const Core::LinAlg::SparseOperator> scatra_scatramanifold_matrix) override;
  };

  // *********************************************************************************************
  // SSI problem is organized in sparse sub matrices
  class AssembleStrategyBlockSparse : public AssembleStrategyBlock
  {
   public:
    explicit AssembleStrategyBlockSparse(
        std::shared_ptr<const SSI::Utils::SSIMaps> ssi_maps, const bool is_scatra_manifold);

    void assemble_scatra_scatra(std::shared_ptr<Core::LinAlg::SparseOperator> systemmatrix,
        std::shared_ptr<const Core::LinAlg::SparseOperator> scatra_scatra_matrix) override;

    void assemble_scatra_structure(std::shared_ptr<Core::LinAlg::SparseOperator> systemmatrix,
        std::shared_ptr<const Core::LinAlg::SparseOperator> scatra_structure_matrix) override;

    void assemble_structure_structure(std::shared_ptr<Core::LinAlg::SparseOperator> systemmatrix,
        std::shared_ptr<const Core::LinAlg::SparseMatrix> structure_structure_matrix) override;

    void assemble_structure_scatra(std::shared_ptr<Core::LinAlg::SparseOperator> systemmatrix,
        std::shared_ptr<const Core::LinAlg::SparseOperator> structure_scatra_matrix) override;

    void assemble_scatramanifold_scatramanifold(
        std::shared_ptr<Core::LinAlg::SparseOperator> systemmatrix,
        std::shared_ptr<const Core::LinAlg::SparseOperator> scatramanifold_scatramanifold_matrix)
        override;

    void assemble_scatramanifold_structure(
        std::shared_ptr<Core::LinAlg::SparseOperator> systemmatrix,
        std::shared_ptr<const Core::LinAlg::SparseOperator> scatramanifold_structure_matrix)
        override;

    void assemble_scatramanifold_scatra(std::shared_ptr<Core::LinAlg::SparseOperator> systemmatrix,
        std::shared_ptr<const Core::LinAlg::SparseOperator> scatramanifold_scatra_matrix) override;

    void assemble_scatra_scatramanifold(std::shared_ptr<Core::LinAlg::SparseOperator> systemmatrix,
        std::shared_ptr<const Core::LinAlg::SparseOperator> scatra_scatramanifold_matrix) override;
  };

  //======================================================================================================
  // SSI problem is organized in one sparse matrix
  class AssembleStrategySparse : public AssembleStrategyBase
  {
   public:
    explicit AssembleStrategySparse(
        std::shared_ptr<const SSI::Utils::SSIMaps> ssi_maps, const bool is_scatra_manifold);

    void assemble_scatra_scatra(std::shared_ptr<Core::LinAlg::SparseOperator> systemmatrix,
        std::shared_ptr<const Core::LinAlg::SparseOperator> scatra_scatra_matrix) override;

    void assemble_scatra_structure(std::shared_ptr<Core::LinAlg::SparseOperator> systemmatrix,
        std::shared_ptr<const Core::LinAlg::SparseOperator> scatra_structure_matrix) override;

    void assemble_structure_structure(std::shared_ptr<Core::LinAlg::SparseOperator> systemmatrix,
        std::shared_ptr<const Core::LinAlg::SparseMatrix> structure_structure_matrix) override;

    void assemble_structure_scatra(std::shared_ptr<Core::LinAlg::SparseOperator> systemmatrix,
        std::shared_ptr<const Core::LinAlg::SparseOperator> structure_scatra_matrix) override;

    void assemble_scatramanifold_scatramanifold(
        std::shared_ptr<Core::LinAlg::SparseOperator> systemmatrix,
        std::shared_ptr<const Core::LinAlg::SparseOperator> scatramanifold_scatramanifold_matrix)
        override;

    void assemble_scatramanifold_structure(
        std::shared_ptr<Core::LinAlg::SparseOperator> systemmatrix,
        std::shared_ptr<const Core::LinAlg::SparseOperator> scatramanifold_structure_matrix)
        override;

    void assemble_scatramanifold_scatra(std::shared_ptr<Core::LinAlg::SparseOperator> systemmatrix,
        std::shared_ptr<const Core::LinAlg::SparseOperator> scatramanifold_scatra_matrix) override;

    void assemble_scatra_scatramanifold(std::shared_ptr<Core::LinAlg::SparseOperator> systemmatrix,
        std::shared_ptr<const Core::LinAlg::SparseOperator> scatra_scatramanifold_matrix) override;
  };

  //! build specific assemble strategy
  std::shared_ptr<SSI::AssembleStrategyBase> build_assemble_strategy(
      std::shared_ptr<const SSI::Utils::SSIMaps> ssi_maps, const bool is_scatra_manifold,
      Core::LinAlg::MatrixType matrixtype_ssi, Core::LinAlg::MatrixType matrixtype_scatra);

}  // namespace SSI
FOUR_C_NAMESPACE_CLOSE

#endif
