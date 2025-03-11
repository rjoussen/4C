// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_FEM_GENERAL_UTILS_INTEGRATION_HPP
#define FOUR_C_FEM_GENERAL_UTILS_INTEGRATION_HPP

#include "4C_config.hpp"

#include "4C_fem_general_cell_type.hpp"
#include "4C_fem_general_cell_type_traits.hpp"
#include "4C_fem_general_element.hpp"
#include "4C_utils_exceptions.hpp"

#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string>

FOUR_C_NAMESPACE_OPEN

namespace Core::FE
{
  // forward declarations
  struct IntegrationPoints1D;
  struct IntegrationPoints2D;
  struct IntegrationPoints3D;


  //! supported 3d integration rules
  //!
  //! integration rules can be distinguished into open and closed rules
  //! hereby, open rules are rules where all integration points are inside the domain
  //! and none on the boundary of the integrated domain
  enum class GaussRule3D
  {
    undefined,           ///< use this to initialize a gaussrule, if you don't know the value,
                         ///< yet!
    hex_1point,          ///< GAUSS INTEGRATION, DEG.OF PRECISION 1 (open)
    hex_8point,          ///< GAUSS INTEGRATION, DEG.OF PRECISION 3 (open)
    hex_lobatto8point,   ///< DEG.OF PRECISION 1 (closed)
    hex_18point,         ///< GAUSS INTEGRATION, DEG.OF PRECISION 5 in xy-plane, 3 in z-direction
                         ///< (open)
    hex_27point,         ///< GAUSS INTEGRATION, DEG.OF PRECISION 5 (open)
    hex_lobatto27point,  ///< GAUSS-LOBATTO (closed)
    hex_64point,
    hex_lobatto64point,  ///< GAUSS-LOBATTO (closed)
    hex_125point,
    hex_lobatto125point,  ///< GAUSS-LOBATTO (closed)
    hex_216point,
    hex_lobatto216point,  ///< GAUSS-LOBATTO (closed)
    hex_343point,
    hex_lobatto343point,  ///< GAUSS-LOBATTO (closed)
    hex_512point,
    hex_lobatto512point,  ///< GAUSS-LOBATTO(closed)
    hex_729point,
    hex_lobatto729point,  ///< GAUSS-LOBATTO (closed)
    hex_1000point,
    hex_lobatto1000point,    ///< GAUSS-LOBATTO (closed)
    tet_1point,              ///< GAUSS INTEGRATION, DEG.OF PRECISION 1 (open)
    tet_4point,              ///< GAUSS INTEGRATION, DEG.OF PRECISION 2 (open)
    tet_4point_gauss_radau,  ///< GAUSS INTEGRATION, DEG.OF PRECISION 2 (closed)
    tet_5point,              ///< GAUSS INTEGRATION, DEG.OF PRECISION 3 (open)
    tet_10point,             ///< KEAST3           , DEG.OF PRECISION 3 (closed)
    tet_11point,             ///< KEAST5           , DEG.OF PRECISION 4 (open)
    tet_15point,             ///< KEAST6           , DEG.OF PRECISION 5 (closed)
    tet_24point,             ///< KEAST7           , DEG.OF PRECISION 6 (open)
    tet_64point_peano,       ///< GAUSS-LOBATTO    , DEG.OF PRECISION 6 (open)   (Peano, 1982)
    tet_45point,             ///< KEAST9           , DEG.OF PRECISION 8 (open)
    tet_125point_peano,      ///< GAUSS-LOBATTO    , DEG.OF PRECISION 9 (open)   (Peano, 1982)
    tet_343point_peano,      ///< GAUSS-LOBATTO    , DEG.OF PRECISION 13 (open)  (Peano, 1982)
    tet_729point_peano,      ///< GAUSS-LOBATTO    , DEG.OF PRECISION 16 (open)  (Peano, 1982)
    wedge_1point,            ///< GAUSS INTEGRATION, DEG.OF PRECISION ?
    wedge_6point,            ///< GAUSS INTEGRATION, DEG.OF PRECISION ?
    wedge_9point,            ///< GAUSS INTEGRATION, DEG.OF PRECISION ?
    pyramid_1point,          ///< GAUSS INTEGRATION, DEG.OF PRECISION ?
    pyramid_8point           ///< GAUSS INTEGRATION, DEG.OF PRECISION ?
  };

  /*!
   * @brief Evaluates the Gauss rule given the number of Gauss points
   *
   * @note Only Gauss-Legendre integration rules are returned.
   *
   * @tparam distype (in) : discretization type
   * @param numgp (in): number of Gauss points
   * @return GaussRule1D / GaussRule2D / GaussRule3D depending on the dimension of the
   * discretization
   */
  template <Core::FE::CellType distype>
    requires(FE::is_hex<distype>)
  GaussRule3D num_gauss_points_to_gauss_rule(unsigned numgp)
  {
    switch (numgp)
    {
      case 1:
        return GaussRule3D::hex_1point;
      case 8:
        return GaussRule3D::hex_8point;
      case 18:
        return GaussRule3D::hex_18point;
      case 27:
        return GaussRule3D::hex_27point;
      case 64:
        return GaussRule3D::hex_64point;
      case 125:
        return GaussRule3D::hex_125point;
      case 216:
        return GaussRule3D::hex_216point;
      case 343:
        return GaussRule3D::hex_343point;
      case 512:
        return GaussRule3D::hex_512point;
      case 729:
        return GaussRule3D::hex_729point;
      case 1000:
        return GaussRule3D::hex_1000point;
      default:
        FOUR_C_THROW("I don't know the GaussRule3D for hex elements with %d Gauss points", numgp);
    }
  }

  template <Core::FE::CellType distype>
    requires(FE::is_nurbs<distype>)
  GaussRule3D num_gauss_points_to_gauss_rule(unsigned numgp)
  {
    switch (numgp)
    {
      case 27:
        return GaussRule3D::hex_27point;
      default:
        FOUR_C_THROW("I don't know the GaussRule3D for nurbs elements with %d Gauss points", numgp);
    }
  }

  template <Core::FE::CellType distype>
    requires(FE::is_tet<distype>)
  GaussRule3D num_gauss_points_to_gauss_rule(unsigned numgp)
  {
    switch (numgp)
    {
      case 1:
        return GaussRule3D::tet_1point;
      case 4:
        return GaussRule3D::tet_4point;
      case 11:
        return GaussRule3D::tet_11point;
      case 24:
        return GaussRule3D::tet_24point;
      case 45:
        return GaussRule3D::tet_45point;
      default:
        FOUR_C_THROW(
            "I don't know an open GaussRule3D for tet elements with %d Gauss points", numgp);
    }
  }

  template <Core::FE::CellType distype>
    requires(FE::is_wedge<distype>)
  GaussRule3D num_gauss_points_to_gauss_rule(unsigned numgp)
  {
    switch (numgp)
    {
      case 1:
        return GaussRule3D::wedge_1point;
      case 6:
        return GaussRule3D::wedge_6point;
      case 9:
        return GaussRule3D::wedge_9point;
      default:
        FOUR_C_THROW("I don't know the GaussRule3D for wedge elements with %d Gauss points", numgp);
    }
  }

  template <Core::FE::CellType distype>
    requires(FE::is_pyramid<distype>)
  GaussRule3D num_gauss_points_to_gauss_rule(unsigned numgp)
  {
    switch (numgp)
    {
      case 1:
        return GaussRule3D::pyramid_1point;
      case 8:
        return GaussRule3D::pyramid_8point;
      default:
        FOUR_C_THROW(
            "I don't know the GaussRule3D for pyramid elements with %d Gauss points", numgp);
    }
  }

  /// get total number of GPs of the 3D Gauss rule at compile time
  template <GaussRule3D rule>
  struct GaussRule3DToNumGaussPoints
  {
  };
  template <>
  struct GaussRule3DToNumGaussPoints<GaussRule3D::undefined>
  {
    static constexpr unsigned value = 0;
  };
  template <>
  struct GaussRule3DToNumGaussPoints<GaussRule3D::hex_1point>
  {
    static constexpr unsigned value = 1;
  };
  template <>
  struct GaussRule3DToNumGaussPoints<GaussRule3D::hex_8point>
  {
    static constexpr unsigned value = 8;
  };
  template <>
  struct GaussRule3DToNumGaussPoints<GaussRule3D::hex_lobatto8point>
  {
    static constexpr unsigned value = 8;
  };
  template <>
  struct GaussRule3DToNumGaussPoints<GaussRule3D::hex_18point>
  {
    static constexpr unsigned value = 18;
  };
  template <>
  struct GaussRule3DToNumGaussPoints<GaussRule3D::hex_27point>
  {
    static constexpr unsigned value = 27;
  };
  template <>
  struct GaussRule3DToNumGaussPoints<GaussRule3D::hex_lobatto27point>
  {
    static constexpr unsigned value = 27;
  };
  template <>
  struct GaussRule3DToNumGaussPoints<GaussRule3D::hex_64point>
  {
    static constexpr unsigned value = 64;
  };
  template <>
  struct GaussRule3DToNumGaussPoints<GaussRule3D::hex_lobatto64point>
  {
    static constexpr unsigned value = 64;
  };
  template <>
  struct GaussRule3DToNumGaussPoints<GaussRule3D::hex_125point>
  {
    static constexpr unsigned value = 125;
  };
  template <>
  struct GaussRule3DToNumGaussPoints<GaussRule3D::hex_lobatto125point>
  {
    static constexpr unsigned value = 125;
  };
  template <>
  struct GaussRule3DToNumGaussPoints<GaussRule3D::hex_216point>
  {
    static constexpr unsigned value = 216;
  };
  template <>
  struct GaussRule3DToNumGaussPoints<GaussRule3D::hex_lobatto216point>
  {
    static constexpr unsigned value = 216;
  };
  template <>
  struct GaussRule3DToNumGaussPoints<GaussRule3D::hex_343point>
  {
    static constexpr unsigned value = 343;
  };
  template <>
  struct GaussRule3DToNumGaussPoints<GaussRule3D::hex_lobatto343point>
  {
    static constexpr unsigned value = 343;
  };
  template <>
  struct GaussRule3DToNumGaussPoints<GaussRule3D::hex_512point>
  {
    static constexpr unsigned value = 512;
  };
  template <>
  struct GaussRule3DToNumGaussPoints<GaussRule3D::hex_lobatto512point>
  {
    static constexpr unsigned value = 512;
  };
  template <>
  struct GaussRule3DToNumGaussPoints<GaussRule3D::hex_729point>
  {
    static constexpr unsigned value = 729;
  };
  template <>
  struct GaussRule3DToNumGaussPoints<GaussRule3D::hex_lobatto729point>
  {
    static constexpr unsigned value = 729;
  };
  template <>
  struct GaussRule3DToNumGaussPoints<GaussRule3D::hex_1000point>
  {
    static constexpr unsigned value = 1000;
  };
  template <>
  struct GaussRule3DToNumGaussPoints<GaussRule3D::hex_lobatto1000point>
  {
    static constexpr unsigned value = 1000;
  };
  template <>
  struct GaussRule3DToNumGaussPoints<GaussRule3D::tet_1point>
  {
    static constexpr unsigned value = 1;
  };
  template <>
  struct GaussRule3DToNumGaussPoints<GaussRule3D::tet_4point>
  {
    static constexpr unsigned value = 4;
  };
  template <>
  struct GaussRule3DToNumGaussPoints<GaussRule3D::tet_4point_gauss_radau>
  {
    static constexpr unsigned value = 4;
  };
  template <>
  struct GaussRule3DToNumGaussPoints<GaussRule3D::tet_5point>
  {
    static constexpr unsigned value = 5;
  };
  template <>
  struct GaussRule3DToNumGaussPoints<GaussRule3D::tet_10point>
  {
    static constexpr unsigned value = 10;
  };
  template <>
  struct GaussRule3DToNumGaussPoints<GaussRule3D::tet_11point>
  {
    static constexpr unsigned value = 11;
  };
  template <>
  struct GaussRule3DToNumGaussPoints<GaussRule3D::tet_15point>
  {
    static constexpr unsigned value = 15;
  };
  template <>
  struct GaussRule3DToNumGaussPoints<GaussRule3D::tet_24point>
  {
    static constexpr unsigned value = 24;
  };
  template <>
  struct GaussRule3DToNumGaussPoints<GaussRule3D::tet_64point_peano>
  {
    static constexpr unsigned value = 64;
  };
  template <>
  struct GaussRule3DToNumGaussPoints<GaussRule3D::tet_45point>
  {
    static constexpr unsigned value = 45;
  };
  template <>
  struct GaussRule3DToNumGaussPoints<GaussRule3D::tet_125point_peano>
  {
    static constexpr unsigned value = 125;
  };
  template <>
  struct GaussRule3DToNumGaussPoints<GaussRule3D::tet_343point_peano>
  {
    static constexpr unsigned value = 343;
  };
  template <>
  struct GaussRule3DToNumGaussPoints<GaussRule3D::tet_729point_peano>
  {
    static constexpr unsigned value = 729;
  };
  template <>
  struct GaussRule3DToNumGaussPoints<GaussRule3D::wedge_1point>
  {
    static constexpr unsigned value = 1;
  };
  template <>
  struct GaussRule3DToNumGaussPoints<GaussRule3D::wedge_6point>
  {
    static constexpr unsigned value = 6;
  };
  template <>
  struct GaussRule3DToNumGaussPoints<GaussRule3D::wedge_9point>
  {
    static constexpr unsigned value = 9;
  };
  template <>
  struct GaussRule3DToNumGaussPoints<GaussRule3D::pyramid_1point>
  {
    static constexpr unsigned value = 1;
  };
  template <>
  struct GaussRule3DToNumGaussPoints<GaussRule3D::pyramid_8point>
  {
    static constexpr unsigned value = 8;
  };

  //! supported 2d integration rules
  //!
  //! integration rules can be distinguished into open and closed rules
  //! hereby, open rules are rules where all integration points are inside the domain
  //! and none on the boundary of the integrated domain
  enum class GaussRule2D
  {
    undefined,            ///< use this to initialize a gaussrule, if you don't know the value,
                          ///< yet!
    quad_1point,          ///< degree of precision: 1  (open)
    quad_4point,          ///< degree of precision: 3  (open)
    quad_6point,          ///< degree of precision: 3/5  (open) (2 by 3 gauss points for quad6
                          ///< element
    quad_9point,          ///< degree of precision: 5  (open)
    quad_16point,         ///< degree of precision: 7  (open)
    quad_25point,         ///< degree of precision: 9  (open)
    quad_lobatto25point,  ///< degree of precision: 7  (closed)
    quad_36point,         ///< degree of precision: 11 (open)
    quad_49point,         ///< degree of precision: 13 (open)
    quad_64point,         ///< degree of precision: 15 (open)
    quad_81point,         ///< degree of precision: 17 (open)
    quad_100point,        ///< degree of precision: 19 (open)
    quad_256point,        ///< degree of precision: 31 (open)
    quad_400point,        ///< degree of precision: 39 (open)
    quad_1024point,       ///< degree of precision: 63 (open)
    tri_1point,           ///< degree of precision: 1  (open)
    tri_3point,           ///< degree of precision: 2  (open)   (Hughes, The Finite Element Method)
    tri_3point_gauss_radau,  ///< degree of precision: 2  (closed) (Hughes, The Finite
                             ///< Element Method)
    tri_4point,   ///< degree of precision: 3  (open)   (Hughes, The Finite Element Method)
    tri_6point,   ///< degree of precision: 4  (open)   (Hughes, The Finite Element Method)
    tri_7point,   ///< degree of precision: 5  (open)   (Gilbert Strange, George Fix, An
                  ///< Analysis of the Finite Element Method, Cambridge, 1973,ISBN:
                  ///< 096140888X,LC: TA335.S77.)
    tri_12point,  ///< degree of precision: 6  (open)   (Hughes, The Finite Element
                  ///< Method)
    tri_16point,  ///< degree of precision: 8  (open)   (Linbo Zhang, Tao Cui and Hui Liu:
                  ///< A SET OF SYMMETRIC QUADRATURE RULES ON TRIANGLES
                  ///<                                   AND TETRAHEDRA, Journal of
                  ///<                                   Computational Mathematics,
                  ///<                                   Vol.27, No.1, 2009, 89-96.)
    tri_37point,  ///< degree of precision: 13 (open)   (rule from ACM TOMS algorithm
                  ///< #706)
    tri_64point   ///< degree of precision: 15 (open)   (essentially a product of two 8
                  ///< point 1D Gauss-Legendre rules)
  };


  template <Core::FE::CellType distype>
    requires(FE::is_quad<distype>)
  GaussRule2D num_gauss_points_to_gauss_rule(unsigned numgp)
  {
    switch (numgp)
    {
      case 1:
        return GaussRule2D::quad_1point;
      case 4:
        return GaussRule2D::quad_4point;
      case 6:
        return GaussRule2D::quad_6point;
      case 9:
        return GaussRule2D::quad_9point;
      case 16:
        return GaussRule2D::quad_16point;
      case 25:
        return GaussRule2D::quad_25point;
      case 49:
        return GaussRule2D::quad_49point;
      case 64:
        return GaussRule2D::quad_64point;
      case 81:
        return GaussRule2D::quad_81point;
      case 100:
        return GaussRule2D::quad_100point;
      case 256:
        return GaussRule2D::quad_256point;
      case 400:
        return GaussRule2D::quad_400point;
      case 1024:
        return GaussRule2D::quad_1024point;
      default:
        FOUR_C_THROW("I don't know the GaussRule2D for quad elements with %d Gauss points", numgp);
    }
  }

  template <Core::FE::CellType distype>
    requires(FE::is_tri<distype>)
  GaussRule2D num_gauss_points_to_gauss_rule(unsigned numgp)
  {
    switch (numgp)
    {
      case 1:
        return GaussRule2D::tri_1point;
      case 3:
        return GaussRule2D::tri_3point;
      case 4:
        return GaussRule2D::tri_4point;
      case 6:
        return GaussRule2D::tri_6point;
      case 7:
        return GaussRule2D::tri_7point;
      case 12:
        return GaussRule2D::tri_12point;
      case 16:
        return GaussRule2D::tri_16point;
      case 37:
        return GaussRule2D::tri_37point;
      case 64:
        return GaussRule2D::tri_64point;
      default:
        FOUR_C_THROW("I don't know the GaussRule2D for tri elements with %d Gauss points", numgp);
    }
  }

  //! supported 1d integration rules
  enum class GaussRule1D
  {
    undefined,            ///< use this to initialize a gaussrule, if you don't know the value,
                          ///< yet!
    line_1point,          ///< degree of precision: 1  (open)
    line_2point,          ///< degree of precision: 3  (open)
    line_3point,          ///< degree of precision: 5  (open)
    line_4point,          ///< degree of precision: 7  (open)
    line_5point,          ///< degree of precision: 9  (open)
    line_6point,          ///< degree of precision: 11 (open)
    line_7point,          ///< degree of precision: 13 (open)
    line_8point,          ///< degree of precision: 15 (open)
    line_9point,          ///< degree of precision: 17 (open)
    line_10point,         ///< degree of precision: 19 (open)
    line_11point,         ///< degree of precision: 21 (open)
    line_12point,         ///< degree of precision: 23 (open)
    line_13point,         ///< degree of precision: 25 (open)
    line_14point,         ///< degree of precision: 27 (open)
    line_15point,         ///< degree of precision: 29 (open)
    line_16point,         ///< degree of precision: 31 (open)
    line_20point,         ///< degree of precision: 39 (open)
    line_32point,         ///< degree of precision: 63 (open)
    line_50point,         ///< degree of precision: 99 (open)
    line_lobatto2point,   ///< degree of precision: 1  (closed)
    line_lobatto3point,   ///< degree of precision: 3  (closed)
    line_lobatto4point,   ///< degree of precision: 5  (closed)
    line_lobatto5point,   ///< degree of precision: 7  (closed)
    line_lobatto6point,   ///< degree of precision: 9  (closed)
    line_lobatto7point,   ///< degree of precision: 11 (closed)
    line_lobatto8point,   ///< degree of precision: 13 (closed)
    line_lobatto9point,   ///< degree of precision: 15 (closed)
    line_lobatto10point,  ///< degree of precision: 17 (closed)
    line_lobatto11point,  ///< degree of precision: 19 (closed)
    line_lobatto12point,  ///< degree of precision: 21 (closed)
    line_lobatto13point,  ///< degree of precision: 23 (closed)
    line_lobatto14point   ///< degree of precision: 25 (closed)
  };

  template <Core::FE::CellType distype>
    requires(FE::is_line<distype>)
  GaussRule1D num_gauss_points_to_gauss_rule(unsigned numgp)
  {
    switch (numgp)
    {
      case 1:
        return GaussRule1D::line_1point;
      case 2:
        return GaussRule1D::line_2point;
      case 3:
        return GaussRule1D::line_3point;
      case 4:
        return GaussRule1D::line_4point;
      case 5:
        return GaussRule1D::line_5point;
      case 6:
        return GaussRule1D::line_6point;
      case 7:
        return GaussRule1D::line_7point;
      case 8:
        return GaussRule1D::line_8point;
      case 9:
        return GaussRule1D::line_9point;
      case 10:
        return GaussRule1D::line_10point;
      case 11:
        return GaussRule1D::line_11point;
      case 12:
        return GaussRule1D::line_12point;
      case 13:
        return GaussRule1D::line_13point;
      case 14:
        return GaussRule1D::line_14point;
      case 15:
        return GaussRule1D::line_15point;
      case 16:
        return GaussRule1D::line_16point;
      case 20:
        return GaussRule1D::line_20point;
      case 32:
        return GaussRule1D::line_32point;
      case 50:
        return GaussRule1D::line_50point;
      default:
        FOUR_C_THROW("I don't know the GaussRule1D for line elements with %d Gauss points", numgp);
    }
  }

  /*!
   \brief a point with a position and an associated weight
   */
  template <int dim, class Vec>
  class IntegrationPoint
  {
   public:
    //! Standard constructor
    explicit IntegrationPoint(const Vec& pos, const double weight) : position_(pos), weight_(weight)
    {
      return;
    }

    //! copy constructor
    explicit IntegrationPoint(const IntegrationPoint& other)
        : position_(other.position_), weight_(other.weight_)
    {
      return;
    }

    //! return string representation of this class
    std::string to_string() const
    {
      std::stringstream s;
      s << "IntPoint, Position: " << position_ << ", weight: " << weight_;
      return s.str();
    }

    //! get integration point position
    Vec position() const { return position_; }

    //! get integration point weight
    double weight() const { return weight_; }

   private:
    //! 1d, 2d or 3d position of integration point in element coordinates
    const Vec position_;

    //! integration weight
    const double weight_;

    //! Empty constructor
    explicit IntegrationPoint() : position_(), weight_(0.0) { return; }
  };  // class IntegrationPoint


  /*!
  \brief integration parameters

  In this structure the coordinates and weights used by gauss integration are stored.

  \note tetrahedral integration rule weights are stored such that after integration over the
  parameter space, we get the real volume of 0.1666..., NOT the unit volume of 1.0. Integration
  over the hex rules gives an volume of 8.0.
  */
  struct IntegrationPoints3D
  {
    /// constructor
    explicit IntegrationPoints3D(const GaussRule3D intrule  ///< 3d rule
    );

   private:
    static constexpr int max_nquad = 1000;  ///< size of c array
    Core::FE::GaussRule3D intrule_;         ///< associated integration rule name
   public:
    int nquad;                 ///< number of gausspoints
    double qxg[max_nquad][3];  ///< coordinates
    double qwgt[max_nquad];    ///< weights

    /// gauss point coordinates
    const double* point(int point) const
    {
      FOUR_C_ASSERT(point < max_nquad, "Index out of range");
      return qxg[point];
    };

    [[nodiscard]] int num_points() const { return nquad; }

    /// get used integration rule
    Core::FE::GaussRule3D get_int_rule() const { return intrule_; };
  };

  /*!
  \brief integration parameters

  In this structure the coordinates and weights used by gauss integration are stored.

  \note triangular integration rule weights are stored such that after integration over the
  parameter space, we get the real area of 0.5, NOT the unit area of 1.0. Integration over the
  quad rules gives an area of 4.0.
  */
  struct IntegrationPoints2D
  {
    /// constructor
    explicit IntegrationPoints2D(const GaussRule2D intrule  ///< 2d rule
    );

   private:
    static constexpr int max_nquad = 1024;  ///< size of c array
    Core::FE::GaussRule2D intrule_;         ///< associated integration rule name
   public:
    int nquad;                 ///< number of gausspoints
    double qxg[max_nquad][2];  ///< coordinates
    double qwgt[max_nquad];    ///< weights

    /// gauss point coordinates
    const double* point(int point) const
    {
      FOUR_C_ASSERT(point < max_nquad, "Index out of range");
      return qxg[point];
    };

    [[nodiscard]] int num_points() const { return nquad; }

    /// get used integration rule
    Core::FE::GaussRule2D get_int_rule() const { return intrule_; };
  };

  /*!
  \brief integration parameters

  In this structure the coordinates and weights used by numerical integration
  are stored.

      */
  struct IntegrationPoints1D
  {
    /// constructor
    explicit IntegrationPoints1D(const GaussRule1D intrule  ///< 1d rule
    );

   private:
    static constexpr int max_nquad = 50;  ///< size of c array
    Core::FE::GaussRule1D intrule_;       ///< associated integration rule name
   public:
    int nquad;                 ///< number of integration points
    double qxg[max_nquad][1];  ///< coordinates
    double qwgt[max_nquad];    ///< weights

    /// gauss point coordinates
    const double* point(int point) const
    {
      FOUR_C_ASSERT(point < max_nquad, "Index out of range");
      return qxg[point];
    };

    [[nodiscard]] int num_points() const { return nquad; }

    /// get used integration rule
    Core::FE::GaussRule1D get_int_rule() const { return intrule_; };
  };


  namespace Internal
  {
    template <unsigned dim>
    struct DimToIntegrationPoints
    {
    };

    template <>
    struct DimToIntegrationPoints<1>
    {
      using type = Core::FE::IntegrationPoints1D;
    };

    template <>
    struct DimToIntegrationPoints<2>
    {
      using type = Core::FE::IntegrationPoints2D;
    };

    template <>
    struct DimToIntegrationPoints<3>
    {
      using type = Core::FE::IntegrationPoints3D;
    };
  }  // namespace Internal

  template <unsigned dim>
  using IntegrationPoints = typename Internal::DimToIntegrationPoints<dim>::type;



  /*!
  \brief integration points and weights


  This class unifies the access of integration points and corresponding weights
  w.r.t. to the number of space dimensions.
  For each number of space dimensions (i.e., 1,2 and 3) a specialization of this class
  is defined. This avoids performing a switch over the number of space dimensions
  during runtime!
  The already existing structs IntegrationPoints3D,... are reused in the specializations
  in form of member variables of the class.

  */
  template <const int nsd>
  class IntPointsAndWeights
  {
   public:
    /// Constructor
    template <class G>
    explicit IntPointsAndWeights(const G intrule);
  };


  //! specialization for 3D
  template <>
  class IntPointsAndWeights<3>
  {
   public:
    /// Constructor for 3D specialization
    explicit IntPointsAndWeights(const GaussRule3D intrule)
        : intpoints_(IntegrationPoints3D(intrule))
    {
      return;
    };

    //! get reference to the integration points and weights
    const IntegrationPoints3D& ip() const { return intpoints_; };

    inline unsigned num_points() const { return ip().nquad; };
    inline unsigned num_dimension() const { return 3; }
    inline const double* point(int point) const { return ip().qxg[point]; };
    inline double weight(int point) const { return ip().qwgt[point]; };

   private:
    //! integration points and weights for 3D
    IntegrationPoints3D intpoints_;
  };


  //! specialization for 2D
  template <>
  class IntPointsAndWeights<2>
  {
   public:
    /// Constructor for 2D specialization
    explicit IntPointsAndWeights(const GaussRule2D intrule)
        : intpoints_(IntegrationPoints2D(intrule))
    {
      return;
    };

    //! get reference to the integration points and weights
    const IntegrationPoints2D& ip() const { return intpoints_; };

    inline unsigned num_points() const { return ip().nquad; };
    inline unsigned num_dimension() const { return 2; }
    inline const double* point(int point) const { return ip().qxg[point]; };
    inline double weight(int point) const { return ip().qwgt[point]; };

   private:
    //! integration points and weights for 3D
    IntegrationPoints2D intpoints_;
  };


  //! specialization for 1D
  template <>
  class IntPointsAndWeights<1>
  {
   public:
    /// Constructor for 1D specialization
    explicit IntPointsAndWeights(const GaussRule1D intrule)
        : intpoints_(IntegrationPoints1D(intrule))
    {
      return;
    };

    //! get reference to the integration points and weights
    const IntegrationPoints1D& ip() const { return intpoints_; };

    inline unsigned num_points() const { return ip().nquad; };
    inline unsigned num_dimension() const { return 1; }
    inline const double* point(int point) const { return ip().qxg[point]; };
    inline double weight(int point) const { return ip().qwgt[point]; };

   private:
    //! integration points and weights for 3D
    IntegrationPoints1D intpoints_;
  };


}  // namespace Core::FE

FOUR_C_NAMESPACE_CLOSE

#endif
