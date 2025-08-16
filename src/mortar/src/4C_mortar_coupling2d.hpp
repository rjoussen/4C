// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_MORTAR_COUPLING2D_HPP
#define FOUR_C_MORTAR_COUPLING2D_HPP

#include "4C_config.hpp"

#include "4C_inpar_mortar.hpp"

#include <mpi.h>
#include <Teuchos_StandardParameterEntryValidators.hpp>

#include <vector>

FOUR_C_NAMESPACE_OPEN

namespace Core::FE
{
  class Discretization;
}  // namespace Core::FE

namespace Mortar
{
  // forward declarations
  class Element;
  class ParamsInterface;

  /*!
   \brief A class representing the framework for mortar coupling of ONE
   slave element and ONE master element of a mortar interface in
   2D. Concretely, this class controls projection, overlap detection
   and finally integration of the mortar coupling matrices D and M
   and possibly of the weighted gap vector g~.

   */

  class Coupling2d
  {
   public:
    /*!
     \brief Constructor with shape function specification

     Constructs an instance of this class and enables custom shape function types.<br>
     Note that this is \b not a collective call as coupling is
     performed in parallel by individual processes.

     */
    Coupling2d(Core::FE::Discretization& idiscret, int dim, bool quad,
        Teuchos::ParameterList& params, Mortar::Element& sele, Mortar::Element& mele);

    /*!
     \brief Destructor

     */
    virtual ~Coupling2d() = default;
    //! @name Access methods

    /*!
     \brief Get interface discretization

     */
    virtual Core::FE::Discretization& discret() const { return idiscret_; }

    /*!
     \brief Get coupling slave element

     */
    virtual Mortar::Element& slave_element() const { return sele_; }

    /*!
     \brief Get coupling master element

     */
    virtual Mortar::Element& master_element() const { return mele_; }

    /*!
     \brief Get problem dimension (here: 2D)

     */
    virtual const int& n_dim() { return dim_; };

    /*!
     \brief Get coupling / FE ansatz type (true = quadratic)

     */
    virtual const bool& quad() { return quad_; };

    /*!
     \brief Return the LM interpolation / testing type for quadratic FE

     */
    Inpar::Mortar::LagMultQuad lag_mult_quad()
    {
      return Teuchos::getIntegralValue<Inpar::Mortar::LagMultQuad>(imortar_, "LM_QUAD");
    }

    /*!
     \brief Return the LM shape fcn type

     */
    Inpar::Mortar::ShapeFcn shape_fcn()
    {
      return Teuchos::getIntegralValue<Inpar::Mortar::ShapeFcn>(imortar_, "LM_SHAPEFCN");
    }

    /*!
     \brief Get interface contact parameter list

     */
    virtual Teuchos::ParameterList& interface_params() { return imortar_; };

    /*!
     \brief Get projection status of the four end nodes

     */
    virtual const std::vector<bool>& has_proj() { return hasproj_; }

    /*!
     \brief Get overlap regions in parameter spaces

     */
    virtual const std::vector<double>& xi_proj() { return xiproj_; }

    /*!
     \brief Get overlap status

     */
    virtual const bool& overlap() { return overlap_; };

    //@}

    //! @name Evaluation methods

    /*!
     \brief Projection of slave / master pair

     This method projects the nodes of the slave Mortar::Element sele_ onto
     the master Mortar::Element mele_ and vice versa. The variable hasproj
     stores a boolean variable for each of the 4 end nodes, indicating
     whether a feasible projection was found or not. The local element
     coordinates of the 4 projection points are stored in xiproj.

     */
    virtual bool project();

    /*!
     \brief Detect overlap of slave / master pair

     This method evaluates the overlap of the current Mortar::Element pair
     sele_ / mele_ based on the projection status of the 4 end nodes
     (hasproj) and the coordinates of the projection points (xiproj).
     According to the detected overlap case, the integration limits
     are determined and written into xiproj and the overlap status
     is stored in the boolean variable overlap.

     */
    virtual bool detect_overlap();

    /*!
     \brief Integrate overlap of slave / master pair

     This method integrates the overlap of the current Mortar::Element
     pair sele_ / mele_ based on the integration limits (xiproj). The
     integration always includes the Mortar matrices D/M and the gap g.

     */
    virtual bool integrate_overlap(const std::shared_ptr<Mortar::ParamsInterface>& mparams_ptr);

    //@}
   private:
    /*!
     \brief Checks roughly whether the two elements are near

     This methods checks the orientation of slave and master element.
     Projection and overlap detection only make sense if the scalar
     product of the two center normals is negative (i.e. if the two
     elements form an angle smaller than 90 degrees).

     */
    virtual bool rough_check_orient();

   protected:
    /*!
     \brief Get communicator

     */
    virtual MPI_Comm get_comm() const;

    // don't want = operator and cctor
    Coupling2d operator=(const Coupling2d& old);
    Coupling2d(const Coupling2d& old);

    Core::FE::Discretization& idiscret_;  // discretization of the contact interface
    int dim_;                             // problem dimension (here: 2D)
    bool quad_;                           // flag indicating coupling type (true = quadratic)
    Teuchos::ParameterList& imortar_;     // containing contact input parameters
    Mortar::Element& sele_;               // slave element to perform coupling for
    Mortar::Element& mele_;               // master element to perform coupling for

    std::vector<bool> hasproj_;   // projection status of the four end nodes
    std::vector<double> xiproj_;  // overlap regions in parameter spaces
    bool overlap_;                // overlap status
  };
  // class Coupling2d

  /*!
   \brief A class representing the framework for mortar coupling of ONE
   slave element and SEVERAL master elements of a mortar interface in
   2D. This class simply stores and manages several Coupling2d objects.

   */

  class Coupling2dManager
  {
   public:
    /*!
     \brief Constructor with shape function specification

     Constructs an instance of this class and enables custom shape function types.<br>
     Note that this is \b not a collective call as coupling is
     performed in parallel by individual processes.

     */
    Coupling2dManager(Core::FE::Discretization& idiscret, int dim, bool quad,
        Teuchos::ParameterList& params, Mortar::Element* sele, std::vector<Mortar::Element*> mele);

    /*!
     \brief Destructor

     */
    virtual ~Coupling2dManager() = default;
    /*!
     \brief Get coupling slave element

     */
    virtual Mortar::Element& slave_element() const { return *sele_; }

    /*!
     \brief Get one specific coupling master element

     */
    virtual Mortar::Element& master_element(int k) const { return *(mele_[k]); }

    /*!
     \brief Get all coupling master elements

     */
    virtual std::vector<Mortar::Element*> master_elements() const { return mele_; }

    /*!
     \brief Get coupling pairs

     */
    virtual std::vector<std::shared_ptr<Mortar::Coupling2d>>& coupling() { return coup_; }

    /*!
     \brief Get type of quadratic LM interpolation

     */
    virtual Inpar::Mortar::LagMultQuad lag_mult_quad()
    {
      return Teuchos::getIntegralValue<Inpar::Mortar::LagMultQuad>(imortar_, "LM_QUAD");
    }

    /*!
     \brief Get integration type

     */
    virtual Inpar::Mortar::IntType int_type()
    {
      return Teuchos::getIntegralValue<Inpar::Mortar::IntType>(imortar_, "INTTYPE");
    }

    /*!
     \brief Evaluate coupling pairs

     */
    virtual bool evaluate_coupling(const std::shared_ptr<Mortar::ParamsInterface>& mparams_ptr);

    /*!
     \brief Get coupling type

     */
    virtual const bool& quad() { return quad_; };

    /*!
     \brief Return the LM shape fcn type

     */
    Inpar::Mortar::ShapeFcn shape_fcn()
    {
      return Teuchos::getIntegralValue<Inpar::Mortar::ShapeFcn>(imortar_, "LM_SHAPEFCN");
    }

    //@}
   protected:
    /*!
     \brief Evaluate mortar coupling pairs

     */
    virtual void integrate_coupling(const std::shared_ptr<Mortar::ParamsInterface>& mparams_ptr);

    /*!
     \brief Calculate consistent dual shape functions in boundary elements

     */
    virtual void consistent_dual_shape();

   protected:
    // don't want = operator and cctor
    Coupling2dManager operator=(const Coupling2dManager& old);
    Coupling2dManager(const Coupling2dManager& old);

    Core::FE::Discretization& idiscret_;  // discretization of the contact interface
    int dim_;                             // problem dimension (here: 2D)
    bool quad_;                           // flag indicating coupling type (true = quadratic)
    Teuchos::ParameterList& imortar_;     // containing contact input parameters
    Mortar::Element* sele_;               // slave element
    std::vector<Mortar::Element*> mele_;  // master elements
    std::vector<std::shared_ptr<Coupling2d>> coup_;  // coupling pairs
  };
  // class Coupling2dManager
}  // namespace Mortar

FOUR_C_NAMESPACE_CLOSE

#endif
