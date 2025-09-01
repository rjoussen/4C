// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_FEM_NURBS_DISCRETIZATION_CONTROL_POINT_HPP
#define FOUR_C_FEM_NURBS_DISCRETIZATION_CONTROL_POINT_HPP

#include "4C_config.hpp"

#include "4C_comm_parobjectfactory.hpp"
#include "4C_fem_general_node.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Core::FE
{
  namespace Nurbs
  {
    class ControlPointType : public Core::Communication::ParObjectType
    {
     public:
      std::string name() const override { return "ControlPointType"; }

      static ControlPointType& instance() { return instance_; };

      Core::Communication::ParObject* create(Core::Communication::UnpackBuffer& buffer) override;

     private:
      static ControlPointType instance_;
    };

    /*!
    \brief Control points for nurbs surfaces/volumes

    Control points are derived from nodes. For isogeometric
    analysis, they replace the nodes in the (nurbs)
    discretisation.
    They are connected to elements, have degrees of freedom
    etc just like normal nodes (and hence the discretisation
    can handle them correctly), the only difference is that
    their coordinate is just a control point coordinate and
    not a mesh point coordinate.

    */
    class ControlPoint : public Core::Nodes::Node
    {
     public:
      //! @name Enums and Friends

      /*!
      \brief The discretization is a friend of the control point
      */
      friend class Core::FE::Discretization;

      //@}

      //! @name Constructors and destructors and related methods

      /*!
      \brief Standard Constructor

      \param id     (in): A globally unique control point id
      \param coords (in): span of nodal coordinates, length 3
      \param weight (in): nurbs weight
      \param owner  (in): Owner of this node.
      */
      ControlPoint(int id, std::span<const double> coords, const double weight, const int owner);

      /*!
      \brief Copy Constructor

      Makes a deep copy of a control point

      \param old (in): The control point to copy

      */
      ControlPoint(const Core::FE::Nurbs::ControlPoint& old);

      /*!
      \brief Deep copy the derived class and return
             pointer to it

      */
      Core::FE::Nurbs::ControlPoint* clone() const override;


      /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique
      id defined at the top of parobject.H.

      \return the parobject id
      */
      int unique_par_object_id() const override
      {
        return ControlPointType::instance().unique_par_object_id();
      }

      /*!
      \brief Pack this class so it can be communicated

      \ref pack and \ref unpack are used to communicate this CP

      \param data (in/out): a char vector to pack the data into

      */
      void pack(Core::Communication::PackBuffer& data) const override;

      /*!
      \brief Unpack data from a char vector into this class

      \ref pack and \ref unpack are used to communicate this CP

      \param data (in): a char vector to unpack the data from

      */
      void unpack(Core::Communication::UnpackBuffer& buffer) override;

      //@}

      //! @name Additional access methods

      /*!
      \brief Return weight

      \return weight

      */
      virtual inline double w() const { return w_; }

      /*!
      \brief Print this node

      \param os ofstrem
      */
      void print(std::ostream& os) const override;

      //@}

     protected:
      //! nurbs weight
      double w_;

    };  // class ControlPoint

  }  // namespace Nurbs

}  // namespace Core::FE

FOUR_C_NAMESPACE_CLOSE

#endif
