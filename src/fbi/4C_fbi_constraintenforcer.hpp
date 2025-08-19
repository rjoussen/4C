// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_FBI_CONSTRAINTENFORCER_HPP
#define FOUR_C_FBI_CONSTRAINTENFORCER_HPP

#include "4C_config.hpp"

#include "4C_linalg_vector.hpp"
#include "4C_utils_exceptions.hpp"

#include <map>
#include <memory>
#include <vector>

FOUR_C_NAMESPACE_OPEN

namespace Adapter
{
  class FSIStructureWrapper;
  class FluidMovingBoundary;

}  // namespace Adapter

namespace BeamInteraction
{
  class BeamToFluidMeshtyingVtkOutputWriter;
}

namespace Core::Binstrategy
{
  class BinningStrategy;
}
namespace Core::FE
{
  class Discretization;
}  // namespace Core::FE
namespace Core::Elements
{
  class Element;
}

namespace FBI
{
  class FBIGeometryCoupler;
}
namespace Core::LinAlg
{
  class SparseMatrix;
  class SparseOperator;
  class MapExtractor;
}  // namespace Core::LinAlg
namespace Adapter
{
  class ConstraintEnforcerFactory;
  class FBIConstraintBridge;

  /**
   * \brief Abstract class to be overloaded by different constraint enforcement techniques for
   * fluid-beam interaction
   *
   * Depending on the constraint enforcement technique used to couple embedded meshes, e.g. through
   * the penalty method, Lagrange multiplier method, Nitsche method, etc., very different
   * information have to be passed to the participating fields.
   * This class is designed to decouple the decision of which information to pass from the actual
   * (partitioned) algorithm.
   *
   * The interface to the outside world is just a Setup and Evaluate routine, as well as the
   * fluid_to_struct() and struct_to_fluid() routines, which only return one vector, as customary
   * for the FSI::DirichletNeumann algorithm. Everything else, like contributions to the stiffness
   * matrix of the fluid, it provides by calling internal functions. To make this possible, the
   * constraint enforcer need information on the two field.
   *
   */
  class FBIConstraintenforcer
  {
    friend ConstraintEnforcerFactory;
    friend class BeamInteraction::BeamToFluidMeshtyingVtkOutputWriter;

   public:
    /// empty destructor
    virtual ~FBIConstraintenforcer() = default;

    /**
     * \brief Sets up the constraint enforcer
     *
     *\param[in] structure wrapper for the structure solver
     *\param[in] fluid moving boundary wrapper for the fluid solver
     */
    virtual void setup(std::shared_ptr<Adapter::FSIStructureWrapper> structure,
        std::shared_ptr<Adapter::FluidMovingBoundary> fluid);

    /** \brief Hand the binning strategy used for the distribution of the fluid mesh
     *  to the object responsible for the element pair search in the FBI framework
     *
     *  \param[in] binning binning strategy object
     */
    void set_binning(std::shared_ptr<Core::Binstrategy::BinningStrategy> binning);

    /**
     * \brief Computes the coupling matrices
     *
     * This is where the magic happens. The stiffness contributions are integrated using information
     * of the structure elements, the fluid elements and their position relative to each other
     */

    virtual void evaluate();

    /**
     * \brief Recomputes all coupling related quantities without performing a search
     */
    virtual void recompute_coupling_without_pair_creation();

    /**
     * \brief Abstractly, we do everything we have to, to introduce the coupling condition into the
     * structure field.
     *
     * Depending on the constraint enforcement strategy, either only an interface force is returned
     * (mortar Lagrange multiplier partitioned, linearized penalty force partitioned), or a force
     * vector as well as a stiffness matrix with additional
     * information is returned (monolithic formulation, full penalty partitioned).
     *
     *
     * \returns structure force vector
     */

    virtual std::shared_ptr<Core::LinAlg::Vector<double>> fluid_to_structure();

    /**
     * \brief Abstractly, we do everything we have to, to introduce the coupling condition into the
     * slave field.
     *
     * Depending on the constraint enforcement strategy, either only an interface force is returned
     * (Mortar-Lagrange multiplier partitioned, linearized penalty force partitioned), or force
     * vector with additional contributions as well as a stiffness matrix with additional
     * information is returned (monolithic formulation, full penalty partitioned, weak Dirichlet).
     *
     *
     * \returns fluid velocity on the whole domain
     */

    virtual std::shared_ptr<Core::LinAlg::Vector<double>> structure_to_fluid(int step);

    /// Interface to do preparations to solve the fluid
    virtual void prepare_fluid_solve() = 0;

    /// Get function for the structure field #structure_
    std::shared_ptr<const Adapter::FSIStructureWrapper> get_structure() const
    {
      return structure_;
    };

    /// Get function for the bridge object #bridge_
    std::shared_ptr<const Adapter::FBIConstraintBridge> get_bridge() const { return bridge_; };

    /// Handle fbi specific output
    virtual void output(double time, int step) = 0;

   protected:
    /** \brief You will have to use the Adapter::ConstraintEnforcerFactory
     *
     * \param[in] bridge an object managing the pair contributins
     * \param[in] geometrycoupler an object managing the search, parallel communication, etc.
     */
    FBIConstraintenforcer(std::shared_ptr<Adapter::FBIConstraintBridge> bridge,
        std::shared_ptr<FBI::FBIGeometryCoupler> geometrycoupler);

    /**
     * \brief Creates all possible interaction pairs
     *
     * \param[in] pairids a map relating all beam element ids to a set of fluid
     * elements ids which they potentially cut
     */
    void create_pairs(std::shared_ptr<std::map<int, std::vector<int>>> pairids);

    /**
     * \brief Resets the state, i.e. the velocity of all interaction pairs
     */
    void reset_all_pair_states();

    /**
     * \brief Extracts current element dofs that are needed for the computations on pair level
     *
     *\param[in] elements elements belonging to the pair
     *\param[out] beam_dofvec current positions and velocities of the beam element
     *\param[out] fluid_dofvec current positions and velocities of the fluid element
     */
    virtual void extract_current_element_dofs(std::vector<Core::Elements::Element const*> elements,
        std::vector<double>& beam_dofvec, std::vector<double>& fluid_dofvec) const;

    /**
     * \brief Computes the contributions to the stiffness matrix of the fluid field.
     *
     * This has to be implemented differently depending on the concrete constraint enforcement
     * strategy.
     *
     * \returns coupling contributions to the fluid system matrix
     */
    virtual std::shared_ptr<const Core::LinAlg::SparseOperator> assemble_fluid_coupling_matrix()
        const
    {
      FOUR_C_THROW("Not yet implemented! This has to be overloaded by a derived class.\n");
      return nullptr;
    };

    /**
     * \brief Computes the contributions to the stiffness matrix of the structure field.
     *
     * This has to be implemented differently depending on the concrete constraint enforcement
     * strategy.
     *
     * \returns coupling contributions to the structure system matrix
     */
    virtual std::shared_ptr<const Core::LinAlg::SparseMatrix> assemble_structure_coupling_matrix()
        const
    {
      FOUR_C_THROW("Not yet implemented! This has to be overloaded by a derived class.\n");
      return nullptr;
    };

    /**
     * \brief Computes the contributions to the rhs of the structure field.
     *
     * This has to be implemented differently depending on the concrete constraint enforcement
     * strategy.
     *
     * \returns coupling contributions to the structure residual
     */
    virtual std::shared_ptr<Core::LinAlg::Vector<double>> assemble_structure_coupling_residual()
        const
    {
      FOUR_C_THROW("Not yet implemented! This has to be overloaded by a derived class.\n");
      return nullptr;
    };

    /**
     * \brief Computes the contributions to the residuum of the fluid field.
     *
     * This has to be implemented differently depending on the concrete constraint enforcement
     * strategy.
     *
     * \returns coupling contributions to the fluid residual
     */
    virtual std::shared_ptr<Core::LinAlg::Vector<double>> assemble_fluid_coupling_residual() const
    {
      FOUR_C_THROW("Not yet implemented! This has to be overloaded by a derived class.\n");
      return nullptr;
    };

    /// Get function for the fluid field #fluid_
    std::shared_ptr<Adapter::FluidMovingBoundary> get_fluid() const { return fluid_; };

    /// Get function for the structure and the fluid discretization in the vector #discretizations_
    std::vector<std::shared_ptr<Core::FE::Discretization>> get_discretizations() const
    {
      return discretizations_;
    }

    /// Get function for the bridge object #bridge_
    std::shared_ptr<Adapter::FBIConstraintBridge> bridge() const { return bridge_; };

    /// Get map extractor to split fluid velocity and pressure values
    std::shared_ptr<const Core::LinAlg::MapExtractor> get_velocity_pressure_splitter() const
    {
      return velocity_pressure_splitter_;
    }

   private:
    FBIConstraintenforcer() = delete;

    /// underlying fluid of the FSI problem
    std::shared_ptr<Adapter::FluidMovingBoundary> fluid_;

    /// underlying structure of the FSI problem
    std::shared_ptr<Adapter::FSIStructureWrapper> structure_;

    /// Vector containing both (fluid and structure) field discretizations
    std::vector<std::shared_ptr<Core::FE::Discretization>> discretizations_;

    /**
     * \brief Object bridging the gap between the specific implementation of the constraint
     * enforcement technique and the specific implementation of the meshtying discretization
     * approach
     */
    std::shared_ptr<Adapter::FBIConstraintBridge> bridge_;

    /**
     * \brief Object handling geometric operations like the search of embedded pairs as well as the
     * parallel communication
     */
    std::shared_ptr<FBI::FBIGeometryCoupler> geometrycoupler_;

    /// Displacement of the structural column nodes on the current proc
    std::shared_ptr<const Core::LinAlg::Vector<double>> column_structure_displacement_;
    /// Velocity of the structural column nodes on the current proc
    std::shared_ptr<const Core::LinAlg::Vector<double>> column_structure_velocity_;
    /// Velocity of the fluid column nodes on the current proc
    std::shared_ptr<const Core::LinAlg::Vector<double>> column_fluid_velocity_;
    /**
     * \brief Extractor to split fluid values into velocities and pressure DOFs
     *
     * velocities  = OtherVector
     * pressure    = CondVector
     */
    std::shared_ptr<Core::LinAlg::MapExtractor> velocity_pressure_splitter_;
  };
}  // namespace Adapter
FOUR_C_NAMESPACE_CLOSE

#endif
