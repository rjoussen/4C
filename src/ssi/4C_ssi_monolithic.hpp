// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_SSI_MONOLITHIC_HPP
#define FOUR_C_SSI_MONOLITHIC_HPP

#include "4C_config.hpp"

#include "4C_ssi_base.hpp"

#include <Teuchos_Time.hpp>

FOUR_C_NAMESPACE_OPEN

// forward declarations
namespace SSI
{
  enum class ScaTraTimIntType;
}  // namespace SSI

namespace Core::LinAlg
{
  class Solver;
  class Equilibration;
  enum class EquilibrationMethod;
  enum class MatrixType;
}  // namespace Core::LinAlg

namespace SSI
{
  namespace Utils
  {
    class SSIMaps;
    class SSIMatrices;
    class SSIVectors;
  }  // namespace Utils

  class AssembleStrategyBase;
  class ContactStrategyBase;
  class DBCHandlerBase;
  class ManifoldMeshTyingStrategyBase;
  class MeshtyingStrategyBase;
  class ScatraStructureOffDiagCoupling;
  class ScaTraManifoldScaTraFluxEvaluator;

  //! equilibration methods applied to system matrix
  struct SSIMonoEquilibrationMethod
  {
    const Core::LinAlg::EquilibrationMethod global;     //! unique equilibration
    const Core::LinAlg::EquilibrationMethod scatra;     //! equilibration for scatra block
    const Core::LinAlg::EquilibrationMethod structure;  //! equilibration for structure block
  };

  enum class Subproblem : int
  {
    scalar_transport,
    structure,
    manifold
  };

  class SsiMono : public SSIBase
  {
   public:
    //! constructor
    explicit SsiMono(MPI_Comm comm,                     //!< communicator
        const Teuchos::ParameterList& globaltimeparams  //!< parameter list for time integration
    );

    //! return global map of degrees of freedom
    [[nodiscard]] const std::shared_ptr<const Core::LinAlg::Map>& dof_row_map() const;

    void init(MPI_Comm comm, const Teuchos::ParameterList& globaltimeparams,
        const Teuchos::ParameterList& scatraparams, const Teuchos::ParameterList& structparams,
        const std::string& struct_disname, const std::string& scatra_disname, bool isAle) override;

    //! return global map extractor (0: scalar transport, 1: structure, [2: scatra manifold])
    [[nodiscard]] std::shared_ptr<const Core::LinAlg::MultiMapExtractor> maps_sub_problems() const;

    //! return map extractor associated with all degrees of freedom inside the scatra field
    [[nodiscard]] std::shared_ptr<const Core::LinAlg::MultiMapExtractor> block_map_scatra() const;

    //! return map extractor associated with all degrees of freedom inside the scatra manifold field
    [[nodiscard]] std::shared_ptr<const Core::LinAlg::MultiMapExtractor> block_map_scatra_manifold()
        const;

    //! return map extractor associated with all degrees of freedom inside the structural field
    [[nodiscard]] std::shared_ptr<const Core::LinAlg::MultiMapExtractor> block_map_structure()
        const;

    //! return map extractor associated with blocks of global system matrix
    [[nodiscard]] std::shared_ptr<const Core::LinAlg::MultiMapExtractor> block_map_system_matrix()
        const;

    //! Return matrix type of global system matrix
    [[nodiscard]] Core::LinAlg::MatrixType matrix_type() const { return matrixtype_; };

    void read_restart(int restart) override;

    void setup() override;

    void setup_system() override;

    /*!
     * @brief solves the linear system
     *
     * @note in case an equilibration method (scaling of rows and columns) is defined, this is also
     * performed within this call
     */
    void solve_linear_system() const;

    //! this object holds all maps relevant to monolithic scalar transport - structure interaction
    [[nodiscard]] std::shared_ptr<SSI::Utils::SSIMaps> ssi_maps() const { return ssi_maps_; }

    //! return algebraic solver for the global system of equations
    [[nodiscard]] const Core::LinAlg::Solver& solver() const { return *solver_; };

    void timeloop() override;

   private:
    //! strategies for Newton-Raphson convergence check
    class ConvCheckStrategyBase;
    class ConvCheckStrategyElch;
    class ConvCheckStrategyElchScaTraManifold;
    class ConvCheckStrategyStd;

    //! apply the contact contributions to matrices and residuals of the subproblems
    void apply_contact_to_sub_problems() const;

    //! apply the Dirichlet boundary conditions to the ssi system, i.e., matrices and residuals
    void apply_dbc_to_system() const;

    //! apply mesh tying between manifold domains on matrices and residuals
    void apply_manifold_meshtying() const;

    //! perform mesh tying on matrices and residuals as obtained from subproblems
    void apply_meshtying_to_sub_problems() const;

    //! assemble the global system of equations
    void assemble_mat_and_rhs() const;

    //! assemble linearization of scatra residuals to system matrix
    void assemble_mat_scatra() const;

    //! assemble linearization of scatra on manifold residuals to system matrix
    void assemble_mat_scatra_manifold() const;

    //! assemble linearization of structural residuals to system matrix
    void assemble_mat_structure() const;

    //! build null spaces associated with blocks of global system matrix
    void build_null_spaces() const;

    //! calc initial potential field for the monolithic SSI problem including scatra and scatra
    //! manifold fields
    void calc_initial_potential_field();

    //! calc initial time derivative of transported scalars for the monolithic SSI problem including
    //! scatra and scatra manifold fields
    void calc_initial_time_derivative();

    //! call complete on the sub problem matrices
    void complete_subproblem_matrices() const;

    //! distribute the solution to all other fields
    //! \param restore_velocity   restore velocity when structure_field()->set_state() is called
    void distribute_solution_all_fields(bool restore_velocity = false);

    //! evaluate all off-diagonal matrix contributions
    void evaluate_off_diag_contributions() const;

    //! Evaluate ScaTra including copy to corresponding ssi matrix
    void evaluate_scatra() const;

    //! Evaluate ScaTra on manifold incl. coupling with scatra
    void evaluate_scatra_manifold() const;

    //! get matrix and right-hand-side for all subproblems incl. coupling
    void evaluate_subproblems();

    //! build and return vector of equilibration methods for each block of system matrix
    [[nodiscard]] std::vector<Core::LinAlg::EquilibrationMethod> get_block_equilibration() const;

    /*!
     * @note This is only necessary in the first iteration of the simulation, since only there the
     * graph of the matrix changes
     *
     * @return flag indicating if we need to uncomplete the matrices before adding the mesh tying
     * contributions.
     */
    [[nodiscard]] bool is_uncomplete_of_matrices_necessary_for_mesh_tying() const;

    void output() override;

    //! do everything that has to be done once before the first time step
    void prepare_time_loop();

    void prepare_time_step() override;

    //! prepare output for subproblems if needed
    void prepare_output();

    //! print system matrix, rhs, and map of system matrix to file
    void print_system_matrix_rhs_to_mat_lab_format() const;

    //! print time step size, time, and number of the time step
    void print_time_step_info() const;

    //! set scatra manifold solution on scatra field
    void set_scatra_manifold_solution(const Core::LinAlg::Vector<double>& phi) const;

    //! evaluate the current time step using Newton-Raphson iteration
    void newton_loop();

    void update() override;

    //! update ScaTra state within Newton iteration
    void update_iter_scatra() const;

    //! update structure state within the Newton iteration
    void update_iter_structure() const;

    //! Dirichlet boundary condition handler
    std::shared_ptr<SSI::DBCHandlerBase> dbc_handler_;

    //! time for element evaluation and assembly of global system of equations
    double dt_eval_ = 0.0;

    //! time for solution of global system of equations
    double dt_solve_ = 0.0;

    //! equilibration method applied to system matrix
    const struct SSIMonoEquilibrationMethod equilibration_method_;

    //! Evaluation of coupling flux between scatra and manifold on scatra
    std::shared_ptr<SSI::ScaTraManifoldScaTraFluxEvaluator> manifoldscatraflux_;

    //! type of global system matrix in global system of equations
    const Core::LinAlg::MatrixType matrixtype_;

    //! print system matrix, rhs, and map of system matrix to file
    const bool print_matlab_;

    //! relax the tolerance of the linear solver in case it is an iterative solver by scaling the
    //! convergence tolerance with factor @p relax_lin_solver_tolerance_
    const double relax_lin_solver_tolerance_;

    //! relax the tolerance of the linear solver within the first @p relax_lin_solver_step_ steps
    const int relax_lin_solver_iter_step_;

    //! all OD evaluation is in here
    std::shared_ptr<SSI::ScatraStructureOffDiagCoupling> scatrastructure_off_diagcoupling_;

    //! algebraic solver for global system of equations
    std::shared_ptr<Core::LinAlg::Solver> solver_;

    //! this object holds all maps relevant to monolithic scalar transport - structure interaction
    std::shared_ptr<SSI::Utils::SSIMaps> ssi_maps_;

    //! this object holds the system matrix and all subblocks
    std::shared_ptr<SSI::Utils::SSIMatrices> ssi_matrices_;

    //! this object holds the system residuals and increment
    std::shared_ptr<SSI::Utils::SSIVectors> ssi_vectors_;

    //! strategy how to assembly system matrix and rhs
    std::shared_ptr<SSI::AssembleStrategyBase> strategy_assemble_;

    //! strategy how to apply contact contributions to sub matrices and rhs
    std::shared_ptr<SSI::ContactStrategyBase> strategy_contact_;

    //! strategy for Newton-Raphson convergence check
    std::shared_ptr<SSI::SsiMono::ConvCheckStrategyBase> strategy_convcheck_;

    //! all equilibration of global system matrix and RHS is done in here
    std::shared_ptr<Core::LinAlg::Equilibration> strategy_equilibration_;

    //! strategy how to apply mesh tying on manifold domains
    std::shared_ptr<SSI::ManifoldMeshTyingStrategyBase> strategy_manifold_meshtying_;

    //! strategy how to apply mesh tying to system matrix and rhs
    std::unique_ptr<SSI::MeshtyingStrategyBase> strategy_meshtying_;

    //! timer for Newton-Raphson iteration
    std::shared_ptr<Teuchos::Time> timer_;
  };
}  // namespace SSI
FOUR_C_NAMESPACE_CLOSE

#endif
