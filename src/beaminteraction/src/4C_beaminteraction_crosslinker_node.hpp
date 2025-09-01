// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_BEAMINTERACTION_CROSSLINKER_NODE_HPP
#define FOUR_C_BEAMINTERACTION_CROSSLINKER_NODE_HPP

#include "4C_config.hpp"

#include "4C_comm_parobjectfactory.hpp"
#include "4C_fem_general_node.hpp"
#include "4C_utils_exceptions.hpp"

FOUR_C_NAMESPACE_OPEN

// forward declaration ...
namespace Core::Mat
{
  class Material;
}  // namespace Core::Mat

namespace Mat
{
  class CrosslinkerMat;
}

namespace CrossLinking
{
  /*!
  \brief A class for a crosslinker derived from Core::Nodes::Node
  */
  class CrosslinkerNodeType : public Core::Communication::ParObjectType
  {
   public:
    std::string name() const override { return "CrosslinkerNodeType"; };

    static CrosslinkerNodeType& instance() { return instance_; };

    Core::Communication::ParObject* create(Core::Communication::UnpackBuffer& buffer) override;

   private:
    static CrosslinkerNodeType instance_;
  };

  /*!
   \brief A class containing additional data from crosslinker nodes

   This class contains additional information from crosslinker nodes which are
   needed for correct crosslinking in a biopolymer network simulation. Note they are only
   available on the node's processor (ColMap). The class CrosslinkerNodeDataContainer
   must be declared before the Mortar::Node itself.
   */
  class CrosslinkerNodeDataContainer
  {
   public:
    //! @name Constructors and destructors and related methods

    /*!
     \brief Standard Constructor

     */
    CrosslinkerNodeDataContainer();

    /*!
     \brief Destructor

     */
    virtual ~CrosslinkerNodeDataContainer() = default;
    /*!
     \brief Pack this class so that it can be communicated

     This function packs the datacontainer. This is only called
     when the class has been initialized and the pointer to this
     class exists.

     */
    void pack(Core::Communication::PackBuffer& data) const;

    /*!
     \brief Unpack data from a vector into this class

     This function unpacks the datacontainer. This is only called
     when the class has been initialized and the pointer to this
     class exists.

     */
    void unpack(Core::Communication::UnpackBuffer& buffer);

    //@}

    //! @name Access methods

    /*!
     \brief Get current binding spot status of linker
     */
    const std::vector<std::pair<int, int>>& get_cl_b_spot_status()
    {
#ifdef FOUR_C_ENABLE_ASSERTIONS
      // safety check
      if ((int)clbspots_.size() != 2) FOUR_C_THROW("crosslinker has wrong bspot size");
#endif
      return clbspots_;
    }

    /*!
     \brief Get
     */
    void set_cl_b_spot_status(std::vector<std::pair<int, int>> clbspots)
    {
      clbspots_ = clbspots;
      return;
    }

    /*!
    \brief Get current number of bonds of crosslinker
    */
    const int& get_number_of_bonds() { return numbond_; }

    /*!
    \brief Set current number of bonds of crosslinker
    */
    void set_number_of_bonds(int numbond)
    {
      numbond_ = numbond;
      return;
    }

    //@}

   protected:
    // don't want = operator and cctor
    CrosslinkerNodeDataContainer operator=(const CrosslinkerNodeDataContainer& old);
    CrosslinkerNodeDataContainer(const CrosslinkerNodeDataContainer& old);

    /// gid of element to local number of bspot, [0] and [1] first and second bspot of cl
    std::vector<std::pair<int, int>> clbspots_;
    /// number of active bonds
    int numbond_;
  };
  // class CrosslinkerNodeDataContainer

  /*!
   \brief A class for a crosslinker node derived from Core::Nodes::Node

  This class represents a single crosslinker involved in a biopolymer network simulation.
  * note:
  *
  *
  *
   */

  class CrosslinkerNode : public Core::Nodes::Node
  {
   public:
    //! @name Enums and Friends

    /*!
     \brief The discretization is a friend of Mortar::Node
     */
    friend class Core::FE::Discretization;

    //@}

    //! @name Constructors and destructors and related methods

    /*!
    \brief Standard Constructor

    \param id     (in): A globally unique node id
    \param coords (in): span of nodal coordinates, length 3
    \param owner  (in): Owner of this node.

    */
    CrosslinkerNode(int id, std::span<const double> coords, const int owner);

    /*!
    \brief Copy Constructor

    Makes a deep copy of a CrosslinkerNode

    */
    CrosslinkerNode(const CrossLinking::CrosslinkerNode& old);

    /*!
     \brief Deep copy the derived class and return pointer to it

     */
    CrossLinking::CrosslinkerNode* clone() const override;



    /*!
     \brief Return unique ParObject id

     every class implementing ParObject needs a unique id defined at the
     top of lib/parobject.H.

     */
    int unique_par_object_id() const override
    {
      return CrosslinkerNodeType::instance().unique_par_object_id();
    }

    /*!
     \brief Pack this class so it can be communicated

     \ref pack and \ref unpack are used to communicate this node

     */
    void pack(Core::Communication::PackBuffer& data) const override;

    /*!
     \brief Unpack data from a char vector into this class

     \ref pack and \ref unpack are used to communicate this node

     */
    void unpack(Core::Communication::UnpackBuffer& buffer) override;

    //@}

    //! @name Access methods

    /*!
     \brief Print this Mortar::Node
     */
    void print(std::ostream& os) const override;


    /*!
    \brief Return material of this node

    This method returns the material associated with this crosslinker node

    */
    inline std::shared_ptr<Mat::CrosslinkerMat> get_material() const
    {
      if (mat_ == nullptr) FOUR_C_THROW("No crosslinker material attached.");
      return mat_;
    }

    //  /*!
    //   \brief Initializes the data container of the node
    //
    //   With this function, the container with crosslinker binding specific quantities/information
    //   is initialized.
    //
    //   */
    //  virtual void initialize_data_container();

    /*!
     \brief Set material for crosslinker node

     Matnum needs to be assigned to a crosslinker type in the input file

     */
    virtual void set_material(int const matnum);


    /*!
     \brief Set material for crosslinker node
     */
    virtual void set_material(std::shared_ptr<Core::Mat::Material> material);

    //  /*!
    //   \brief Resets the data container of the node
    //
    //   With this function, the container with crosslinker binding specific quantities/information
    //   is deleted / reset to nullptr pointer
    //
    //   */
    //  virtual void ResetDataContainer();

    //@}


   protected:
    //  /// information of crosslinker binding status, this is different for each crosslinker
    //  //  and may change each time step
    //  std::shared_ptr<CrossLinking::CrosslinkerNodeDataContainer> cldata_;

    /// this contains information that does not change during the simulation time and is
    //  the same for a subset of crosslinker, we only need one object for each subset
    std::shared_ptr<Mat::CrosslinkerMat> mat_;


  };  // class CrosslinkerNode
}  // namespace CrossLinking

// << operator
std::ostream& operator<<(std::ostream& os, const CrossLinking::CrosslinkerNode& crosslinker_node);

FOUR_C_NAMESPACE_CLOSE

#endif
