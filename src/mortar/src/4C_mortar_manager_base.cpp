// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_mortar_manager_base.hpp"



FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 |  ctor (public)                                             popp 01/10|
 *----------------------------------------------------------------------*/
Mortar::ManagerBase::ManagerBase()
{
  //**********************************************************************
  // empty constructor (this is an abstract base class)
  //**********************************************************************
  // Setup of the mortar contact library is done by a derived class. This
  // derived class is specific to the FEM code into which the mortar contact
  // library is meant to be integrated. For 4C this is realized via the
  // CONTACT::ContactManager class! There the following actions are performed:
  //**********************************************************************
  // 1) get problem dimension (2D or 3D)
  // 2) read and check contact input parameters
  // 3) read and check contact boundary conditions
  // 4) build contact interfaces
  //**********************************************************************
  // A similar process also applies to mortar meshtying libraries. Again
  // a specific derived class is needed. For 4C this is realized via the
  // CONTACT::MeshtyingManager class!
  //**********************************************************************

  comm_ = MPI_COMM_SELF;
}

FOUR_C_NAMESPACE_CLOSE
