// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_POROELAST_UTILS_HPP
#define FOUR_C_POROELAST_UTILS_HPP

#include "4C_config.hpp"

#include "4C_coupling_volmortar_utils.hpp"
#include "4C_fem_general_element.hpp"
#include "4C_linalg_vector.hpp"
#include "4C_poroelast_input.hpp"
#include "4C_utils_parameter_list.fwd.hpp"

#include <memory>

FOUR_C_NAMESPACE_OPEN

namespace Core::FE
{
  class Discretization;
}  // namespace Core::FE

namespace Core::Elements
{
  class Element;
}

namespace Core::LinAlg
{
  class MapExtractor;
}  // namespace Core::LinAlg

namespace Mat
{
  class Material;
}  // namespace Mat

namespace PoroElast
{
  class PoroBase;

  //! coupling types of porous media problem
  enum Coupltype
  {
    fluidfluid,
    fluidstructure,
    undefined
  };

  //! PoroElast::UTILS: Random stuff that might be helpful when dealing with poroelasticity problems
  namespace Utils
  {
    //! check if element is a poro-element
    bool is_poro_element(const Core::Elements::Element* actele);

    //! check if element is a poro-p1-element
    bool is_poro_p1_element(const Core::Elements::Element* actele);

    std::shared_ptr<Core::LinAlg::MapExtractor> build_poro_splitter(Core::FE::Discretization& dis);

    //! create solution algorithm depending on input file
    std::shared_ptr<PoroElast::PoroBase> create_poro_algorithm(
        const Teuchos::ParameterList& timeparams,  //!< problem parameters (i)
        MPI_Comm comm,                             //!< communicator(i)
        bool setup_solve = true,  //!< setup linear solver for Poroelastic problem (only required if
                                  //!< Solve() is called) (i)
        std::shared_ptr<Core::LinAlg::MapExtractor> porosity_splitter =
            nullptr  //! select porosity splitter
    );

    //! reset Material pointers after redistribution
    void set_material_pointers_matching_grid(
        const Core::FE::Discretization& sourcedis, const Core::FE::Discretization& targetdis);

    /*!
     Create volume ghosting:
     To ease the contact search algorithms we'll afford the luxury to ghost all nodes
     on all processors in the general mortar coupling framework.
     Here we additionally ghost volume elements and their nodes
     (required if an evaluation of gradients is required)!
     */
    void create_volume_ghosting(
        Core::FE::Discretization& idiscret);  // redistributed interface discretization of contact!

    /*! Reconnect Face Element - Parent Element Pointers!
     Parent Element need to be ghosted on the processors where Face Elements
     exist already.
     */
    void reconnect_parent_pointers(Core::FE::Discretization& idiscret,
        Core::FE::Discretization& voldiscret, Core::FE::Discretization* voldiscret2 = nullptr);

    //! Determine norm of vector
    double calculate_vector_norm(const enum PoroElast::VectorNorm norm,  //!< norm to use
        const Core::LinAlg::Vector<double>& vect                         //!< the vector of interest
    );

    //! Set the slave and master elements of the face element
    void set_slave_and_master(const Core::FE::Discretization& voldiscret,
        const Core::FE::Discretization* voldiscret2, const Core::LinAlg::Map* elecolmap,
        Core::Elements::FaceElement* faceele);

    //! strategy for material assignment for non matching meshes with poro

    //! Helper class for assigning materials for volumetric coupling of non conforming meshes (poro)
    /*!
     When coupling two overlapping discretizations, most often one discretization needs access
     to the corresponding element/material on the other side. For conforming meshes this is straight
     forward as there is one unique element on the other side and therefore one unique material,
     which can be accessed. However, for non conforming meshes there are potentially several
     elements overlapping. Therefore, some rule for assigning materials is needed. This class is
     meant to do that. It gets the element to which it shall assign a material and a vector of IDs
     of the overlapping elements of the other discretization.

     In case of poro, we also need the initial porosity type of the structural element to be known
     in the fluid element, which is why there is a special strategy for poro. Note that this is not
     yet working for inhomogeneous material properties.
     */
    class PoroMaterialStrategy : public Coupling::VolMortar::Utils::DefaultMaterialStrategy
    {
     public:
      //! constructor
      PoroMaterialStrategy() = default;

      //! assignment of fluid material to structure material
      void assign_material2_to1(const Coupling::VolMortar::VolMortarCoupl* volmortar,
          Core::Elements::Element* ele1, const std::vector<int>& ids_2,
          std::shared_ptr<Core::FE::Discretization> dis1,
          std::shared_ptr<Core::FE::Discretization> dis2) override;

      //! assignment of structure material to fluid material
      void assign_material1_to2(const Coupling::VolMortar::VolMortarCoupl* volmortar,
          Core::Elements::Element* ele2, const std::vector<int>& ids_1,
          std::shared_ptr<Core::FE::Discretization> dis1,
          std::shared_ptr<Core::FE::Discretization> dis2) override;
    };
  }  // namespace Utils

  void print_logo();

}  // namespace PoroElast

FOUR_C_NAMESPACE_CLOSE

#endif
