// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_REDUCED_LUNG_1D_PIPE_FLOW_MAIN_HPP
#define FOUR_C_REDUCED_LUNG_1D_PIPE_FLOW_MAIN_HPP

#include "4C_config.hpp"

#include "4C_fem_discretization.hpp"
#include "4C_linalg_vector.hpp"
#include "4C_reduced_lung_1d_pipe_flow_input.hpp"

FOUR_C_NAMESPACE_OPEN

namespace ReducedLung1dPipeFlow
{
  /**
   * Who owns nodes in a junction.
   */
  struct JunctionInfo
  {
    std::vector<int> node_ids;
    std::vector<int> node_owners;
  };

  /**
   * ID, owner, coordinates of a node
   */
  struct NodeInfo
  {
    int id;
    int owner;
    double x, y, z;
  };

  /**
   * Fills vectors created from nodal maps with data from input file.
   *
   * Input values are listed in @param parameters and accessed to fill the passed vectors.
   * The following vectors are filled with properties for each node, given that they are described
   * elementwise in the inputfile:
   * @param reference_area @param thickness @param Young @param beta @param radius
   * e.g. X = [X_0 X_1 ...]^T The solution
   * vector is set to the initial state with the two DOF/ node:
   *    @param solution = [A_0 u_0 A_1 u_1 ...]^T
   *  The @param discretization is required, as writing to the vectors is
   * performed across MPI ranks.
   */
  void fill_parameters(Parameters& parameters, Core::LinAlg::Vector<double>& solution,
      Core::LinAlg::Vector<double>& reference_area, Core::LinAlg::Vector<double>& thickness,
      Core::LinAlg::Vector<double>& Young, Core::LinAlg::Vector<double>& beta,
      Core::LinAlg::Vector<double>& radius, const Core::FE::Discretization& discretization);

  /**
   * Computes and returns element length.
   */
  double compute_length(const Core::Elements::Element& element);

  /**
   * Computes the matrix of shape functions according to SUPG.
   * The SUPG shape functions are defined as Psi = N + delta * H^T * dNdxi * dxidx.
   * Since the problem has two unknowns, A and u, a matrix of shape functions needs to be defined.
   *
   * The result is stored in @param Psi_matrix .
   * The normal shape functions N are passed in matrix form in @param N_matrix .
   * The associated derivative matrix in @param dNdxi_matrix .
   * The @param flux_jacobian is defined by H = [u A, c^2/ A u] and @param delta is defined by L /
   * (2*lambda_max).
   * The derivative dxidx is defined by 1 / @param L
   */
  void compute_psi_matrix(Core::LinAlg::Matrix<2, 4>& Psi_matrix,
      const Core::LinAlg::Matrix<2, 4>& N_matrix, const Core::LinAlg::Matrix<2, 4>& dNdxi_matrix,
      const Core::LinAlg::Matrix<2, 2>& flux_jacobian, const double L, const double delta);

  /**
   * Get conditions for A and u when flow Q is prescribed. Computed through Newton-Raphson with
   * f = (W_in - W_out)^4 /1024 (rho/beta)^2 /2 * (W_in + W_out) .
   * W_out is known from the domain, W_in needs to be determined from prescribed Q.
   * Constants needed for computation are passed in @param input and the prescribed flow in
   * @param Q_condition .
   * Known parameters at the node are @param boundary_A0 @param characteristic_W_outgoing
   * @param beta
   * The computed conditions for A and u are written to @param A_condition and @param u_condition
   */
  void conditions_from_newton_raphson(const Parameters& input, const double Q_condition,
      const double boundary_A0, const double characteristic_W_outgoing, const double beta,
      double& A_condition, double& u_condition);

  /**
   * Calculate residual vector for Newton-Raphson iterations.
   * @param f residual vector to be computed
   * @param x solution vector
   * Further information needed for the computation of characteristics are the number of connected
   * nodes at a junction @param N_connected_nodes and the geometry or material properties at the
   * respective nodes.
   * The first N rows in f refer to characteristic speed information, N+1 to mass conservation, and
   * the remaining to pressure conservation.
   */
  void compute_residual(Core::LinAlg::SerialDenseVector& f, const int N_connected_nodes,
      const Core::LinAlg::SerialDenseVector& x, const std::vector<double>& junction_normal,
      const std::vector<double>& junction_ref_area_A0,
      const std::vector<double>& junction_characteristic_out,
      const std::vector<double>& junction_beta, const double density_rho);

  /**
   * Jacobian computation for Newton-Raphson iterations.
   * The Jacobian is defined as follows: J = (par f)/(par x), with residual f and @param x the
   * solution vector.
   * The first N (@param N_connected_nodes) rows in f refer to characteristic speed information, N+1
   * to mass conservation, and the remaining to pressure conservation.
   */
  void compute_jacobian(Core::LinAlg::SerialDenseMatrix& jacobian, const int N_connected_nodes,
      const Core::LinAlg::SerialDenseVector& x, const std::vector<double>& junction_normal,
      const std::vector<double>& junction_beta, const double density_rho);

  /**
   * Compute the values for A_np and u_np for all nodes at the junctions.
   * Information of all nodes being part of the junction need to be stored locally
   * (vectors_junction) so that conditions A and u can be computed for all nodes on the process they
   * belong to.
   * @param all_junctions contains the information of nodes of all junctions on the rank
   * @param dof_update is the vector where the conditions for the primary variables are stored so
   * that they can be applied to the rhs vector
   */
  void get_conditions_at_junctions(const double density_rho,
      const Core::LinAlg::Vector<double>& characteristics_for_junction,
      Core::LinAlg::Vector<double>& solution_for_junction,
      const Core::LinAlg::Vector<double>& normals_for_junction,
      const Core::LinAlg::Vector<double>& area0_for_junction,
      const Core::LinAlg::Vector<double>& beta_for_junction,
      const std::vector<JunctionInfo>& all_junctions, Core::LinAlg::Vector<double>& dof_update);

  /**
   * Compute the contribution of the junctions to the rhs vector with +=(F_h - F).
   */
  void update_rhs_with_junction_properties(const double density_rho,
      Core::LinAlg::Vector<double>& rhs_junction,
      const Core::LinAlg::Vector<double>& solution_update_junction,
      const Core::LinAlg::Vector<double>& beta, const Core::LinAlg::Vector<double>& reference_area,
      const Core::LinAlg::Vector<double>& normals, const Core::LinAlg::Vector<double>& solution,
      const std::vector<JunctionInfo>& all_junctions);

  /**
   * Calculate flow velocity and vessel area as primary variables in a pipe over time. Pressure
   * and flux can be determined from these two variables. For calculation, the governing
   * equations are solved on each element and then assembled for the whole geometry. (I) : (par
   * A / par t) + (par (A * u) / par x) = 0 (II) : (par u / par t) + 1/rho (par p / par x) +
   * viscosity * u / A = 0
   */
  void main();


}  // namespace ReducedLung1dPipeFlow

FOUR_C_NAMESPACE_CLOSE

#endif
