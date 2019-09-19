// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 2 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.

  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the list of
  copyright holders.
*/
/*!
 * \file
 *
 * \copydoc Opm::VcfvStencil
 */
#ifndef EWOMS_VCFV_STENCIL_HH
#define EWOMS_VCFV_STENCIL_HH

#include <opm/models/utils/quadraturegeometries.hh>

#include <opm/material/common/Unused.hpp>
#include <opm/material/common/Exceptions.hpp>

#include <dune/grid/common/intersectioniterator.hh>
#include <dune/grid/common/mcmgmapper.hh>
#include <dune/geometry/referenceelements.hh>

#if HAVE_DUNE_LOCALFUNCTIONS
#include <dune/localfunctions/lagrange/pqkfactory.hh>
#endif // HAVE_DUNE_LOCALFUNCTIONS

#include <dune/common/version.hh>

#include <vector>

namespace Opm {

/*!
 * \brief The types of reference elements available.
 */
enum ElementType
{
    none,
    simplex,
    cube,
};

/*!
 * \cond SKIP_THIS
 */
template <class Scalar, unsigned dim, unsigned basicGeomType>
class VcfvScvGeometries;

////////////////////
// local geometries for 1D elements
////////////////////
template <class Scalar>
class VcfvScvGeometries<Scalar, /*dim=*/1, ElementType::cube>
{
    enum { dim = 1 };
    enum { numScv = 2 };

public:
    typedef Opm::QuadrialteralQuadratureGeometry<Scalar, dim> ScvLocalGeometry;

    static void init()
    {
        // 1D LINE SEGMENTS
        Scalar scvCorners[numScv][ScvLocalGeometry::numCorners][dim] = {
            { // corners of the first sub control volume
                { 0.0 },
                { 0.5 }
            },

            { // corners of the second sub control volume
                { 0.5 },
                { 1.0 }
            }
        };
        for (unsigned scvIdx = 0; scvIdx < numScv; ++scvIdx)
            scvGeoms_[scvIdx].setCorners(scvCorners[scvIdx], ScvLocalGeometry::numCorners);
    }

    static const ScvLocalGeometry& get(unsigned scvIdx)
    { return scvGeoms_[scvIdx]; }

private:
    static ScvLocalGeometry scvGeoms_[numScv];
};

template <class Scalar>
typename VcfvScvGeometries<Scalar, /*dim=*/1, ElementType::cube>::ScvLocalGeometry
VcfvScvGeometries<Scalar, /*dim=*/1, ElementType::cube>::scvGeoms_[
    VcfvScvGeometries<Scalar, /*dim=*/1, ElementType::cube>::numScv];

template <class Scalar>
class VcfvScvGeometries<Scalar, /*dim=*/1, ElementType::simplex>
{
    enum { dim = 1 };
    enum { numScv = 2 };

public:
    typedef Opm::QuadrialteralQuadratureGeometry<Scalar, dim> ScvLocalGeometry;

    static const ScvLocalGeometry& get(unsigned scvIdx OPM_UNUSED)
    {
        throw std::logic_error("Not implemented: VcfvScvGeometries<Scalar, 1, ElementType::simplex>");
    }
};

////////////////////
// local geometries for 2D elements
////////////////////
template <class Scalar>
class VcfvScvGeometries<Scalar, /*dim=*/2, ElementType::simplex>
{
    enum { dim = 2 };
    enum { numScv = 3 };

public:
    typedef Opm::QuadrialteralQuadratureGeometry<Scalar, dim> ScvLocalGeometry;

    static const ScvLocalGeometry& get(unsigned scvIdx)
    { return scvGeoms_[scvIdx]; }

    static void init()
    {
        // 2D SIMPLEX
        Scalar scvCorners[numScv][ScvLocalGeometry::numCorners][dim] =
        {
            { // SCV 0 corners
                { 0.0,   0.0 },
                { 1.0/2, 0.0 },
                { 1.0/3, 1.0/3 },
                { 0.0,   1.0/2 },
            },

            { // SCV 1 corners
                { 1.0/2, 0.0 },
                { 1.0,   0.0 },
                { 1.0/3, 1.0/3 },
                { 1.0/2, 1.0/2 },
            },

            { // SCV 2 corners
                { 0.0,   1.0/2 },
                { 1.0/3, 1.0/3 },
                { 0.0,   1.0 },
                { 1.0/2, 1.0/2 },
            }
        };

        for (unsigned scvIdx = 0; scvIdx < numScv; ++scvIdx)
            scvGeoms_[scvIdx].setCorners(scvCorners[scvIdx], ScvLocalGeometry::numCorners);
    }

private:
    static ScvLocalGeometry scvGeoms_[numScv];
};

template <class Scalar>
typename VcfvScvGeometries<Scalar, /*dim=*/2, ElementType::simplex>::ScvLocalGeometry
VcfvScvGeometries<Scalar, /*dim=*/2, ElementType::simplex>::scvGeoms_[
    VcfvScvGeometries<Scalar, /*dim=*/2, ElementType::simplex>::numScv];

template <class Scalar>
class VcfvScvGeometries<Scalar, /*dim=*/2, ElementType::cube>
{
    enum { dim = 2 };
    enum { numScv = 4 };

public:
    typedef Opm::QuadrialteralQuadratureGeometry<Scalar, dim> ScvLocalGeometry;

    static const ScvLocalGeometry& get(unsigned scvIdx)
    { return scvGeoms_[scvIdx]; }

    static void init()
    {
        // 2D CUBE
        Scalar scvCorners[numScv][ScvLocalGeometry::numCorners][dim] =
        {
            { // SCV 0 corners
                { 0.0, 0.0 },
                { 0.5, 0.0 },
                { 0.0, 0.5 },
                { 0.5, 0.5 }
            },

            { // SCV 1 corners
                { 0.5, 0.0 },
                { 1.0, 0.0 },
                { 0.5, 0.5 },
                { 1.0, 0.5 }
            },

            { // SCV 2 corners
                { 0.0,  0.5 },
                { 0.5,  0.5 },
                { 0.0,  1.0 },
                { 0.5,  1.0 }
            },

            { // SCV 3 corners
                { 0.5,   0.5 },
                { 1.0,   0.5 },
                { 0.5,   1.0 },
                { 1.0,   1.0 }
            }
        };

        for (unsigned scvIdx = 0; scvIdx < numScv; ++scvIdx)
            scvGeoms_[scvIdx].setCorners(scvCorners[scvIdx], ScvLocalGeometry::numCorners);
    }

    static ScvLocalGeometry scvGeoms_[numScv];
};

template <class Scalar>
typename VcfvScvGeometries<Scalar, /*dim=*/2, ElementType::cube>::ScvLocalGeometry
VcfvScvGeometries<Scalar, /*dim=*/2, ElementType::cube>::scvGeoms_[
    VcfvScvGeometries<Scalar, /*dim=*/2, ElementType::cube>::numScv];

////////////////////
// local geometries for 3D elements
////////////////////
template <class Scalar>
class VcfvScvGeometries<Scalar, /*dim=*/3, ElementType::simplex>
{
    enum { dim = 3 };
    enum { numScv = 4 };

public:
    typedef Opm::QuadrialteralQuadratureGeometry<Scalar, dim> ScvLocalGeometry;

    static const ScvLocalGeometry& get(unsigned scvIdx)
    { return scvGeoms_[scvIdx]; }

    static void init()
    {
        // 3D SIMPLEX
        Scalar scvCorners[numScv][ScvLocalGeometry::numCorners][dim] =
        {
            { // SCV 0 corners
                { 0.0,   0.0,   0.0 },
                { 1.0/2, 0.0,   0.0 },
                { 0.0,   1.0/2, 0.0 },
                { 1.0/3, 1.0/3, 0.0 },

                { 0.0,   0.0,   0.5 },
                { 1.0/3, 0.0,   1.0/3 },
                { 0.0,   1.0/3, 1.0/3 },
                { 1.0/4, 1.0/4, 1.0/4 },
            },

            { // SCV 1 corners
                { 1.0/2, 0.0,   0.0 },
                { 1.0,   0.0,   0.0 },
                { 1.0/3, 1.0/3, 0.0 },
                { 1.0/2, 1.0/2, 0.0 },

                { 1.0/3, 0.0,   1.0/3 },
                { 1.0/2, 0.0,   1.0/2 },
                { 1.0/4, 1.0/4, 1.0/4 },
                { 1.0/3, 1.0/3, 1.0/3 },
            },

            { // SCV 2 corners
                { 0.0,   1.0/2, 0.0 },
                { 1.0/3, 1.0/3, 0.0 },
                { 0.0,   1.0,   0.0 },
                { 1.0/2, 1.0/2, 0.0 },

                { 0.0,   1.0/3, 1.0/3 },
                { 1.0/4, 1.0/4, 1.0/4 },
                { 0.0,   1.0/2, 1.0/2 },
                { 1.0/3, 1.0/3, 1.0/3 },
            },

            { // SCV 3 corners
                { 0.0,   0.0,   1.0/2 },
                { 1.0/3, 0.0,   1.0/3 },
                { 0.0,   1.0/3, 1.0/3 },
                { 1.0/4, 1.0/4, 1.0/4 },

                { 0.0,   0.0,   1.0   },
                { 1.0/2, 0.0,   1.0/2 },
                { 0.0,   1.0/2, 1.0/2 },
                { 1.0/3, 1.0/3, 1.0/3 },
            }
        };

        for (unsigned scvIdx = 0; scvIdx < numScv; ++scvIdx)
            scvGeoms_[scvIdx].setCorners(scvCorners[scvIdx],  ScvLocalGeometry::numCorners);
    }

private:
    static ScvLocalGeometry scvGeoms_[numScv];
};

template <class Scalar>
typename VcfvScvGeometries<Scalar, /*dim=*/3, ElementType::simplex>::ScvLocalGeometry
VcfvScvGeometries<Scalar, /*dim=*/3, ElementType::simplex>::scvGeoms_[
    VcfvScvGeometries<Scalar, /*dim=*/3, ElementType::simplex>::numScv];

template <class Scalar>
class VcfvScvGeometries<Scalar, /*dim=*/3, ElementType::cube>
{
    enum { dim = 3 };
    enum { numScv = 8 };

public:
    typedef Opm::QuadrialteralQuadratureGeometry<Scalar, dim> ScvLocalGeometry;

    static const ScvLocalGeometry& get(unsigned scvIdx)
    { return scvGeoms_[scvIdx]; }

    static void init()
    {
        // 3D CUBE
        Scalar scvCorners[numScv][ScvLocalGeometry::numCorners][dim] =
        {
            { // SCV 0 corners
                { 0.0,   0.0,   0.0 },
                { 1.0/2, 0.0,   0.0 },
                { 0.0,   1.0/2, 0.0 },
                { 1.0/2, 1.0/2, 0.0 },

                { 0.0,   0.0,   1.0/2 },
                { 1.0/2, 0.0,   1.0/2 },
                { 0.0,   1.0/2, 1.0/2 },
                { 1.0/2, 1.0/2, 1.0/2 },
            },

            { // SCV 1 corners
                { 1.0/2, 0.0,   0.0 },
                { 1.0,   0.0,   0.0 },
                { 1.0/2, 1.0/2, 0.0 },
                { 1.0,   1.0/2, 0.0 },

                { 1.0/2, 0.0,   1.0/2 },
                { 1.0,   0.0,   1.0/2 },
                { 1.0/2, 1.0/2, 1.0/2 },
                { 1.0,   1.0/2, 1.0/2 },
            },

            { // SCV 2 corners
                { 0.0,   1.0/2, 0.0 },
                { 1.0/2, 1.0/2, 0.0 },
                { 0.0,   1.0,   0.0 },
                { 1.0/2, 1.0,   0.0 },

                { 0.0,   1.0/2, 1.0/2 },
                { 1.0/2, 1.0/2, 1.0/2 },
                { 0.0,   1.0,   1.0/2 },
                { 1.0/2, 1.0,   1.0/2 },
            },

            { // SCV 3 corners
                { 1.0/2, 1.0/2, 0.0 },
                { 1.0,   1.0/2, 0.0 },
                { 1.0/2, 1.0,   0.0 },
                { 1.0,   1.0,   0.0 },

                { 1.0/2, 1.0/2, 1.0/2 },
                { 1.0,   1.0/2, 1.0/2 },
                { 1.0/2, 1.0,   1.0/2 },
                { 1.0,   1.0,   1.0/2 },
            },

            { // SCV 4 corners
                { 0.0,   0.0,   1.0/2 },
                { 1.0/2, 0.0,   1.0/2 },
                { 0.0,   1.0/2, 1.0/2 },
                { 1.0/2, 1.0/2, 1.0/2 },

                { 0.0,   0.0,   1.0 },
                { 1.0/2, 0.0,   1.0 },
                { 0.0,   1.0/2, 1.0 },
                { 1.0/2, 1.0/2, 1.0 },
            },

            { // SCV 5 corners
                { 1.0/2, 0.0,   1.0/2 },
                { 1.0,   0.0,   1.0/2 },
                { 1.0/2, 1.0/2, 1.0/2 },
                { 1.0,   1.0/2, 1.0/2 },

                { 1.0/2, 0.0,   1.0 },
                { 1.0,   0.0,   1.0 },
                { 1.0/2, 1.0/2, 1.0 },
                { 1.0,   1.0/2, 1.0 },
            },

            { // SCV 6 corners
                { 0.0,   1.0/2, 1.0/2 },
                { 1.0/2, 1.0/2, 1.0/2 },
                { 0.0,   1.0,   1.0/2 },
                { 1.0/2, 1.0,   1.0/2 },

                { 0.0,   1.0/2, 1.0 },
                { 1.0/2, 1.0/2, 1.0 },
                { 0.0,   1.0,   1.0 },
                { 1.0/2, 1.0,   1.0 },
            },

            { // SCV 7 corners
                { 1.0/2, 1.0/2, 1.0/2 },
                { 1.0,   1.0/2, 1.0/2 },
                { 1.0/2, 1.0,   1.0/2 },
                { 1.0,   1.0,   1.0/2 },

                { 1.0/2, 1.0/2, 1.0 },
                { 1.0,   1.0/2, 1.0 },
                { 1.0/2, 1.0,   1.0 },
                { 1.0,   1.0,   1.0 },
            },
        };

        for (unsigned scvIdx = 0; scvIdx < numScv; ++scvIdx)
            scvGeoms_[scvIdx].setCorners(scvCorners[scvIdx], ScvLocalGeometry::numCorners);
    }
private:
    static ScvLocalGeometry scvGeoms_[numScv];
};

template <class Scalar>
typename VcfvScvGeometries<Scalar, /*dim=*/3, ElementType::cube>::ScvLocalGeometry
VcfvScvGeometries<Scalar, /*dim=*/3, ElementType::cube>::scvGeoms_[
    VcfvScvGeometries<Scalar, /*dim=*/3, ElementType::cube>::numScv];

/*!
 * \endcond
 */

/*!
 * \ingroup VcfvDiscretization
 *
 * \brief Represents the finite volume geometry of a single element in
 *        the VCFV discretization.
 *
 * The VCFV discretization is a vertex centered finite volume approach. This
 * means that each vertex corresponds to a control volume which
 * intersects each of the vertex' neighboring elements. If only
 * looking at a single element of the primary grid (which is what this
 * class does), the element is subdivided into multiple fragments of
 * control volumes called sub-control volumes. Each of the element's
 * vertices corresponds to exactly one sub-control volume in this
 * scenario.
 *
 * For the vertex-cented finite volume method the sub-control volumes
 * are constructed by connecting the element's center with each edge
 * of the element.
 */
template <class Scalar, class GridView>
class VcfvStencil
{
    enum{dim = GridView::dimension};
    enum{dimWorld = GridView::dimensionworld};
    enum{maxNC = (dim < 3 ? 4 : 8)};
    enum{maxNE = (dim < 3 ? 4 : 12)};
    enum{maxNF = (dim < 3 ? 1 : 6)};
    enum{maxBF = (dim < 3 ? 8 : 24)};
    typedef typename GridView::ctype CoordScalar;
    typedef typename GridView::Traits::template Codim<0>::Entity            Element;
public:
    typedef typename GridView::Traits::template Codim<dim>::Entity          Entity;
private:
    typedef typename Element::Geometry Geometry;
    typedef Dune::FieldVector<Scalar,dimWorld> DimVector;
    typedef Dune::FieldVector<CoordScalar,dimWorld> GlobalPosition;
    typedef Dune::FieldVector<CoordScalar,dim> LocalPosition;
    typedef typename GridView::IntersectionIterator IntersectionIterator;

    typedef Opm::QuadrialteralQuadratureGeometry<Scalar, dim> ScvLocalGeometry;

#if HAVE_DUNE_LOCALFUNCTIONS
    typedef Dune::PQkLocalFiniteElementCache<CoordScalar, Scalar, dim, 1> LocalFiniteElementCache;
    typedef typename LocalFiniteElementCache::FiniteElementType LocalFiniteElement;
    typedef typename LocalFiniteElement::Traits::LocalBasisType::Traits LocalBasisTraits;
    typedef typename LocalBasisTraits::JacobianType ShapeJacobian;
#endif // HAVE_DUNE_LOCALFUNCTIONS

    Scalar quadrilateralArea(const GlobalPosition& p0,
                             const GlobalPosition& p1,
                             const GlobalPosition& p2,
                             const GlobalPosition& p3)
    {
        return 0.5*std::abs((p3[0] - p1[0])*(p2[1] - p0[1]) - (p3[1] - p1[1])*(p2[0] - p0[0]));
    }

    void crossProduct(DimVector& c, const DimVector& a, const DimVector& b)
    {
        c[0] = a[1]*b[2] - a[2]*b[1];
        c[1] = a[2]*b[0] - a[0]*b[2];
        c[2] = a[0]*b[1] - a[1]*b[0];
    }

    Scalar pyramidVolume(const GlobalPosition& p0,
                         const GlobalPosition& p1,
                         const GlobalPosition& p2,
                         const GlobalPosition& p3,
                         const GlobalPosition& p4)
    {
        DimVector a(p2); a -= p0;
        DimVector b(p3); b -= p1;

        DimVector n;
        crossProduct(n, a, b);

        a = p4; a -= p0;

        return 1.0/6.0*(n*a);
    }

    Scalar prismVolume(const GlobalPosition& p0,
                       const GlobalPosition& p1,
                       const GlobalPosition& p2,
                       const GlobalPosition& p3,
                       const GlobalPosition& p4,
                       const GlobalPosition& p5)
    {
        DimVector a(p4);
        for (unsigned k = 0; k < dimWorld; ++k)
            a[k] -= p0[k];
        DimVector b(p1);
        for (unsigned k = 0; k < dimWorld; ++k)
            b[k] -= p3[k];
        DimVector m;
        crossProduct(m, a, b);

        for (unsigned k = 0; k < dimWorld; ++k)
            a[k] = p1[k] - p0[k];
        for (unsigned k = 0; k < dimWorld; ++k)
            b[k] = p2[k] - p0[k];
        DimVector n;
        crossProduct(n, a, b);
        n += m;

        for (unsigned k = 0; k < dimWorld; ++k)
            a[k] = p5[k] - p0[k];

        return std::abs(1.0/6.0*(n*a));
    }

    Scalar hexahedronVolume(const GlobalPosition& p0,
                            const GlobalPosition& p1,
                            const GlobalPosition& p2,
                            const GlobalPosition& p3,
                            const GlobalPosition& p4,
                            const GlobalPosition& p5,
                            const GlobalPosition& p6,
                            const GlobalPosition& p7)
    {
        return
            prismVolume(p0,p1,p2,p4,p5,p6)
            + prismVolume(p0,p2,p3,p4,p6,p7);
    }

    void normalOfQuadrilateral3D(DimVector& normal,
                                 const GlobalPosition& p0,
                                 const GlobalPosition& p1,
                                 const GlobalPosition& p2,
                                 const GlobalPosition& p3)
    {
        DimVector a(p2);
        for (unsigned k = 0; k < dimWorld; ++k)
            a[k] -= p0[k];
        DimVector b(p3);
        for (unsigned k = 0; k < dimWorld; ++k)
            b[k] -= p1[k];

        crossProduct(normal, a, b);
        normal *= 0.5;
    }

    Scalar quadrilateralArea3D(const GlobalPosition& p0,
                               const GlobalPosition& p1,
                               const GlobalPosition& p2,
                               const GlobalPosition& p3)
    {
        DimVector normal;
        normalOfQuadrilateral3D(normal, p0, p1, p2, p3);
        return normal.two_norm();
    }

    void getFaceIndices(unsigned numElemVertices, unsigned k, unsigned& leftFace, unsigned& rightFace)
    {
        static const unsigned edgeToFaceTet[2][6] = {
            {1, 0, 3, 2, 1, 3},
            {0, 2, 0, 1, 3, 2}
        };
        static const unsigned edgeToFacePyramid[2][8] = {
            {1, 2, 3, 4, 1, 3, 4, 2},
            {0, 0, 0, 0, 3, 2, 1, 4}
        };
        static const unsigned edgeToFacePrism[2][9] = {
            {1, 0, 2, 0, 1, 2, 4, 4, 4},
            {0, 2, 1, 3, 3, 3, 0, 1, 2}
        };
        static const unsigned edgeToFaceHex[2][12] = {
            {0, 2, 3, 1, 4, 1, 2, 4, 0, 5, 5, 3},
            {2, 1, 0, 3, 0, 4, 4, 3, 5, 1, 2, 5}
        };

        switch (numElemVertices) {
        case 4:
            leftFace = edgeToFaceTet[0][k];
            rightFace = edgeToFaceTet[1][k];
            break;
        case 5:
            leftFace = edgeToFacePyramid[0][k];
            rightFace = edgeToFacePyramid[1][k];
            break;
        case 6:
            leftFace = edgeToFacePrism[0][k];
            rightFace = edgeToFacePrism[1][k];
            break;
        case 8:
            leftFace = edgeToFaceHex[0][k];
            rightFace = edgeToFaceHex[1][k];
            break;
        default:
            throw std::logic_error("Not implemented: VcfvStencil::getFaceIndices for "
                                   +std::to_string(numElemVertices)+" vertices");
            break;
        }
    }

    void getEdgeIndices(unsigned numElemVertices, unsigned face, unsigned vert, unsigned& leftEdge, unsigned& rightEdge)
    {
        static const int faceAndVertexToLeftEdgeTet[4][4] = {
                { 0, 0, 2, -1},
                { 0, 0, -1, 3},
                { 1, -1, 1, 3},
                {-1, 2, 2, 4}
        };
        static const int faceAndVertexToRightEdgeTet[4][4] = {
                { 1, 2, 1, -1},
                { 3, 4, -1, 4},
                { 3, -1, 5, 5},
                {-1, 4, 5, 5}
        };
        static const int faceAndVertexToLeftEdgePyramid[5][5] = {
            { 0, 2, 3, 1, -1},
            { 0, -1, 0, -1, 4},
            {-1, 1, -1, 1, 5},
            { 2, 2, -1, -1, 4},
            {-1, -1, 3, 3, 7}
        };
        static const int faceAndVertexToRightEdgePyramid[5][5] = {
            { 2, 1, 0, 3, -1},
            { 4, -1, 6, -1, 6},
            {-1, 5, -1, 7, 7},
            { 4, 5, -1, -1, 5},
            {-1, -1, 6, 7, 6}
        };
        static const int faceAndVertexToLeftEdgePrism[5][6] = {
            { 3, 3, -1, 0, 1, -1},
            { 4, -1, 4, 0, -1, 2},
            {-1, 5, 5, -1, 1, 2},
            { 3, 3, 5, -1, -1, -1},
            {-1, -1, -1, 6, 6, 8}
        };
        static const int faceAndVertexToRightEdgePrism[5][6] = {
            { 0, 1, -1, 6, 6, -1},
            { 0, -1, 2, 7, -1, 7},
            {-1, 1, 2, -1, 8, 8},
            { 4, 5, 4, -1, -1, -1},
            {-1, -1, -1, 7, 8, 7}
        };
        static const int faceAndVertexToLeftEdgeHex[6][8] = {
            { 0, -1, 4, -1, 8, -1, 2, -1},
            {-1, 5, -1, 3, -1, 1, -1, 9},
            { 6, 1, -1, -1, 0, 10, -1, -1},
            {-1, -1, 2, 7, -1, -1, 11, 3},
            { 4, 6, 7, 5, -1, -1, -1, -1},
            {-1, -1, -1, -1, 10, 9, 8, 11}
        };
        static const int faceAndVertexToRightEdgeHex[6][8] = {
            { 4, -1, 2, -1, 0, -1, 8, -1},
            {-1, 1, -1, 5, -1, 9, -1, 3},
            { 0, 6, -1, -1, 10, 1, -1, -1},
            {-1, -1, 7, 3, -1, -1, 2, 11},
            { 6, 5, 4, 7, -1, -1, -1, -1},
            {-1, -1, -1, -1, 8, 10, 11, 9}
        };

        switch (numElemVertices) {
        case 4:
            leftEdge = static_cast<unsigned>(faceAndVertexToLeftEdgeTet[face][vert]);
            rightEdge = static_cast<unsigned>(faceAndVertexToRightEdgeTet[face][vert]);
            break;
        case 5:
            leftEdge = static_cast<unsigned>(faceAndVertexToLeftEdgePyramid[face][vert]);
            rightEdge = static_cast<unsigned>(faceAndVertexToRightEdgePyramid[face][vert]);
            break;
        case 6:
            leftEdge = static_cast<unsigned>(faceAndVertexToLeftEdgePrism[face][vert]);
            rightEdge = static_cast<unsigned>(faceAndVertexToRightEdgePrism[face][vert]);
            break;
        case 8:
            leftEdge = static_cast<unsigned>(faceAndVertexToLeftEdgeHex[face][vert]);
            rightEdge = static_cast<unsigned>(faceAndVertexToRightEdgeHex[face][vert]);
            break;
        default:
            throw std::logic_error("Not implemented: VcfvStencil::getFaceIndices for "
                                   +std::to_string(numElemVertices)+" vertices");
            break;
        }
    }

public:
    //! exported Mapper type
#if DUNE_VERSION_NEWER(DUNE_GRID, 2,6)
    typedef Dune::MultipleCodimMultipleGeomTypeMapper<GridView> Mapper;
#else
    typedef Dune::MultipleCodimMultipleGeomTypeMapper<GridView, Dune::MCMGVertexLayout> Mapper;
#endif

    class ScvGeometry
    {
    public:
        const GlobalPosition center() const
        { return global(localGeometry_->center()); }

        const GlobalPosition corner(unsigned cornerIdx) const
        { return global(localGeometry_->corner(cornerIdx)); }

        const GlobalPosition global(const LocalPosition& localPos) const
        { return element_->geometry().global(localPos); }

        const ScvLocalGeometry& localGeometry() const
        { return *localGeometry_; }

        const ScvLocalGeometry *localGeometry_;
        const Element *element_;
    };

    struct SubControlVolume //! finite volume intersected with element
    {
        const GlobalPosition& globalPos() const
        { return global; }

        const GlobalPosition center() const
        { return geometry_.center(); }

        Scalar volume() const
        { return volume_; }

        const ScvLocalGeometry& localGeometry() const
        { return geometry_.localGeometry(); }

        const ScvGeometry& geometry() const
        { return geometry_; }

        //! local vertex position
        LocalPosition local;
        //! global vertex position
        GlobalPosition global;
        //! space occupied by the sub-control volume
        Scalar volume_;
        //! The geometry of the sub-control volume in local coordinates.
        ScvGeometry geometry_;

        //! derivative of shape function at the center of the sub control volume
        Dune::FieldVector<DimVector, maxNC> gradCenter;
    };

    struct SubControlVolumeFace //! interior face of a sub control volume
    {
        const DimVector& normal() const
        { return normal_; }

        unsigned short interiorIndex() const
        { return i; }

        unsigned short exteriorIndex() const
        { return j; }

        Scalar area() const
        { return area_; }

        const LocalPosition& localPos() const
        { return ipLocal_; }

        const GlobalPosition& integrationPos() const
        { return ipGlobal_; }

        //! scvf seperates corner i and j of elem
        unsigned short i,j;
        //! integration point in local coords
        LocalPosition ipLocal_;
        //! integration point in global coords
        GlobalPosition ipGlobal_;
        //! normal on face pointing to CV j or outward of the domain with length equal to |scvf|
        DimVector normal_;
        //! area of face
        Scalar area_;
    };

    //! compatibility typedef
    typedef SubControlVolumeFace BoundaryFace;

    VcfvStencil(const GridView& gridView, const Mapper& mapper)
        : gridView_(gridView)
        , vertexMapper_(mapper )
        , element_(*gridView.template begin</*codim=*/0>())
    {
        // try to check if the mapper really maps the vertices
        assert(static_cast<int>(gridView.size(/*codim=*/dimWorld)) == static_cast<int>(mapper.size()));

        static bool localGeometriesInitialized = false;
        if (!localGeometriesInitialized) {
            localGeometriesInitialized = true;

            VcfvScvGeometries<Scalar, /*dim=*/1, ElementType::cube>::init();
            VcfvScvGeometries<Scalar, /*dim=*/2, ElementType::cube>::init();
            VcfvScvGeometries<Scalar, /*dim=*/2, ElementType::simplex>::init();
            VcfvScvGeometries<Scalar, /*dim=*/3, ElementType::cube>::init();
            VcfvScvGeometries<Scalar, /*dim=*/3, ElementType::simplex>::init();
        }
    }

    /*!
     * \brief Update the non-geometric part of the stencil.
     *
     * I.e., indices and neighboring information, but nothing else...
     */
    void updateTopology(const Element& e)
    {
        element_ = e;

        numVertices = e.subEntities(/*codim=*/dim);
        numEdges = e.subEntities(/*codim=*/dim-1);
        numFaces = (dim<3)?0:e.subEntities(/*codim=*/1);

        numBoundarySegments_ = 0; // TODO: really required here(?)

        // compute the local and global coordinates of the element
        const Geometry& geometry = e.geometry();
        geometryType_ = geometry.type();
        const auto& referenceElement = Dune::ReferenceElements<CoordScalar,dim>::general(geometryType_);
        for (unsigned vertexIdx = 0; vertexIdx < numVertices; vertexIdx++) {
            subContVol[vertexIdx].local = referenceElement.position(static_cast<int>(vertexIdx), dim);
            subContVol[vertexIdx].global = geometry.corner(static_cast<int>(vertexIdx));
        }
    }

    void updatePrimaryTopology(const Element& element)
    {
        // since all degrees of freedom in a stencil are "primary" DOFs for the
        // vertex-centered finite volume method, there's no difference to
        // updateTopology()
        updateTopology(element);
    }

    void update(const Element& e)
    {
        updateTopology(e);

        const Geometry& geometry = e.geometry();
        geometryType_ = geometry.type();

        const auto& referenceElement = Dune::ReferenceElements<CoordScalar,dim>::general(geometryType_);

        elementVolume = geometry.volume();
        elementLocal = referenceElement.position(0,0);
        elementGlobal = geometry.global(elementLocal);

        // corners:
        for (unsigned vert = 0; vert < numVertices; vert++) {
            subContVol[vert].local = referenceElement.position(static_cast<int>(vert), dim);
            subContVol[vert].global = geometry.global(subContVol[vert].local);
        }

        // edges:
        for (unsigned edge = 0; edge < numEdges; edge++) {
            edgeCoord[edge] = geometry.global(referenceElement.position(static_cast<int>(edge), dim-1));
        }

        // faces:
        for (unsigned face = 0; face < numFaces; face++) {
            faceCoord[face] = geometry.global(referenceElement.position(static_cast<int>(face), 1));
        }

        // fill sub control volume data use specialization for this
        // \todo maybe it would be a good idea to migrate everything
        // which is dependend of the grid's dimension to
        // _VcfvFVElemGeomHelper in order to benefit from more aggressive
        // compiler optimizations...
        fillSubContVolData_();

        // fill sub control volume face data:
        for (unsigned k = 0; k < numEdges; k++) { // begin loop over edges / sub control volume faces
            unsigned short i = static_cast<unsigned short>(referenceElement.subEntity(static_cast<int>(k), dim-1, 0, dim));
            unsigned short j = static_cast<unsigned short>(referenceElement.subEntity(static_cast<int>(k), dim-1, 1, dim));
            if (numEdges == 4 && (i == 2 || j == 2))
                std::swap(i, j);
            subContVolFace[k].i = i;
            subContVolFace[k].j = j;

            // calculate the local integration point and
            // the face normal. note that since dim is a
            // constant which is known at compile time
            // the compiler can optimize away all if
            // cases which don't apply.
            LocalPosition ipLocal_;
            DimVector diffVec;
            if (dim==1) {
                subContVolFace[k].ipLocal_ = 0.5;
                subContVolFace[k].normal_ = 1.0;
                subContVolFace[k].area_ = 1.0;
                ipLocal_ = subContVolFace[k].ipLocal_;
            }
            else if (dim==2) {
                ipLocal_ = referenceElement.position(static_cast<int>(k), dim-1) + elementLocal;
                ipLocal_ *= 0.5;
                subContVolFace[k].ipLocal_ = ipLocal_;
                for (unsigned m = 0; m < dimWorld; ++m)
                    diffVec[m] = elementGlobal[m] - edgeCoord[k][m];
                subContVolFace[k].normal_[0] = diffVec[1];
                subContVolFace[k].normal_[1] = -diffVec[0];

                for (unsigned m = 0; m < dimWorld; ++m)
                    diffVec[m] = subContVol[j].global[m] - subContVol[i].global[m];
                // make sure the normal points to the right direction
                if (subContVolFace[k].normal_ * diffVec < 0)
                    subContVolFace[k].normal_ *= -1;

                subContVolFace[k].area_ = subContVolFace[k].normal_.two_norm();
                subContVolFace[k].normal_ /= subContVolFace[k].area_;
            }
            else if (dim==3) {
                unsigned leftFace;
                unsigned rightFace;
                getFaceIndices(numVertices, k, leftFace, rightFace);
                ipLocal_ = referenceElement.position(static_cast<int>(k), dim-1) + elementLocal
                    + referenceElement.position(static_cast<int>(leftFace), 1)
                    + referenceElement.position(static_cast<int>(rightFace), 1);
                ipLocal_ *= 0.25;
                subContVolFace[k].ipLocal_ = ipLocal_;
                normalOfQuadrilateral3D(subContVolFace[k].normal_,
                                        edgeCoord[k], faceCoord[rightFace],
                                        elementGlobal, faceCoord[leftFace]);
                subContVolFace[k].area_ = subContVolFace[k].normal_.two_norm();
                subContVolFace[k].normal_ /= subContVolFace[k].area_;
            }

            // get the global integration point and the Jacobian inverse
            subContVolFace[k].ipGlobal_ = geometry.global(ipLocal_);
        } // end loop over edges / sub control volume faces

        // fill boundary face data:
        IntersectionIterator endit = gridView_.iend(e);
        for (IntersectionIterator it = gridView_.ibegin(e); it != endit; ++it) {
            const typename IntersectionIterator::Intersection& intersection = *it ;

            if ( ! intersection.boundary())
                continue;

            unsigned face = static_cast<unsigned>(intersection.indexInInside());
            unsigned numVerticesOfFace = static_cast<unsigned>(referenceElement.size(static_cast<int>(face), 1, dim));
            for (unsigned vertInFace = 0; vertInFace < numVerticesOfFace; vertInFace++)
            {
                unsigned short vertInElement = static_cast<unsigned short>(referenceElement.subEntity(static_cast<int>(face), 1, static_cast<int>(vertInFace), dim));
                unsigned bfIdx = numBoundarySegments_;
                ++numBoundarySegments_;

                if (dim == 1) {
                    boundaryFace_[bfIdx].ipLocal_ = referenceElement.position(static_cast<int>(vertInElement), dim);
                    boundaryFace_[bfIdx].area_ = 1.0;
                }
                else if (dim == 2) {
                    boundaryFace_[bfIdx].ipLocal_ = referenceElement.position(static_cast<int>(vertInElement), dim)
                        + referenceElement.position(static_cast<int>(face), 1);
                    boundaryFace_[bfIdx].ipLocal_ *= 0.5;
                    boundaryFace_[bfIdx].area_ = 0.5 * intersection.geometry().volume();
                }
                else if (dim == 3) {
                    unsigned leftEdge;
                    unsigned rightEdge;
                    getEdgeIndices(numVertices, face, vertInElement, leftEdge, rightEdge);
                    boundaryFace_[bfIdx].ipLocal_ = referenceElement.position(static_cast<int>(vertInElement), dim)
                        + referenceElement.position(static_cast<int>(face), 1)
                        + referenceElement.position(static_cast<int>(leftEdge), dim-1)
                        + referenceElement.position(static_cast<int>(rightEdge), dim-1);
                    boundaryFace_[bfIdx].ipLocal_ *= 0.25;
                    boundaryFace_[bfIdx].area_ =
                        quadrilateralArea3D(subContVol[vertInElement].global,
                                            edgeCoord[rightEdge],
                                            faceCoord[face],
                                            edgeCoord[leftEdge]);
                }
                else
                    throw std::logic_error("Not implemented:VcfvStencil for dim = "+std::to_string(dim));

                boundaryFace_[bfIdx].ipGlobal_ = geometry.global(boundaryFace_[bfIdx].ipLocal_);
                boundaryFace_[bfIdx].i = vertInElement;
                boundaryFace_[bfIdx].j = vertInElement;

                // ASSUME constant normal on the segment of the boundary face
                boundaryFace_[bfIdx].normal_ = intersection.centerUnitOuterNormal();
            }
        }

        updateScvGeometry(e);
    }

    void updateScvGeometry(const Element& element)
    {
        auto geomType = element.geometry().type();

        // get the local geometries of the sub control volumes
        if (geomType.isTriangle() || geomType.isTetrahedron()) {
            for (unsigned vertIdx = 0; vertIdx < numVertices; ++vertIdx) {
                subContVol[vertIdx].geometry_.element_ = &element;
                subContVol[vertIdx].geometry_.localGeometry_ =
                    &VcfvScvGeometries<Scalar, dim, ElementType::simplex>::get(vertIdx);
            }
        }
        else if (geomType.isLine() || geomType.isQuadrilateral() || geomType.isHexahedron()) {
            for (unsigned vertIdx = 0; vertIdx < numVertices; ++vertIdx) {
                subContVol[vertIdx].geometry_.element_ = &element;
                subContVol[vertIdx].geometry_.localGeometry_ =
                    &VcfvScvGeometries<Scalar, dim, ElementType::cube>::get(vertIdx);
            }
        }
        else
            throw std::logic_error("Not implemented: SCV geometries for non hexahedron elements");
    }

#if HAVE_DUNE_LOCALFUNCTIONS
    void updateCenterGradients()
    {
        const auto& localFiniteElement = feCache_.get(element_.type());
        const auto& geom = element_.geometry();

        std::vector<ShapeJacobian> localJac;

        for (unsigned scvIdx = 0; scvIdx < numVertices; ++ scvIdx) {
            const auto& localCenter = subContVol[scvIdx].localGeometry().center();
            localFiniteElement.localBasis().evaluateJacobian(localCenter, localJac);
            const auto& globalPos = subContVol[scvIdx].geometry().center();

            const auto jacInvT = geom.jacobianInverseTransposed(globalPos);
            for (unsigned vert = 0; vert < numVertices; vert++) {
                jacInvT.mv(localJac[vert][0], subContVol[scvIdx].gradCenter[vert]);
            }
        }
    }
#endif

    unsigned numDof() const
    { return numVertices; }

    unsigned numPrimaryDof() const
    { return numDof(); }

    Dune::PartitionType partitionType(unsigned scvIdx) const
    { return element_.template subEntity</*codim=*/dim>(scvIdx)->partitionType(); }

    const SubControlVolume& subControlVolume(unsigned scvIdx) const
    {
        assert(0 <= scvIdx && scvIdx < numDof());
        return subContVol[scvIdx];
    }

    unsigned numInteriorFaces() const
    { return numEdges; }

    unsigned numBoundaryFaces() const
    { return numBoundarySegments_; }

    const SubControlVolumeFace& interiorFace(unsigned faceIdx) const
    { return subContVolFace[faceIdx]; }

    const BoundaryFace& boundaryFace(unsigned bfIdx) const
    { return boundaryFace_[bfIdx]; }

    /*!
     * \brief Return the global space index given the index of a degree of
     *        freedom.
     */
    unsigned globalSpaceIndex(unsigned dofIdx) const
    {
        assert(0 <= dofIdx && dofIdx < numDof());

        return static_cast<unsigned>(vertexMapper_.subIndex(element_, static_cast<int>(dofIdx), /*codim=*/dim));
    }

    /*!
     * \brief Return the global space index given the index of a degree of
     *        freedom.
     */
    Entity entity(unsigned dofIdx) const
    {
        assert(0 <= dofIdx && dofIdx < numDof());
        return element_.template subEntity<dim>(static_cast<int>(dofIdx));
    }

private:
#if __GNUC__ || __clang__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wpragmas"
#pragma GCC diagnostic ignored "-Warray-bounds"
#endif
    void fillSubContVolData_()
    {
        if (dim == 1) {
            // 1D
            subContVol[0].volume_ = 0.5*elementVolume;
            subContVol[1].volume_ = 0.5*elementVolume;
        }
        else if (dim == 2) {
            switch (numVertices) {
            case 3: // 2D, triangle
                subContVol[0].volume_ = elementVolume/3;
                subContVol[1].volume_ = elementVolume/3;
                subContVol[2].volume_ = elementVolume/3;
                break;
            case 4: // 2D, quadrilinear
                subContVol[0].volume_ =
                    quadrilateralArea(subContVol[0].global,
                                      edgeCoord[2],
                                      elementGlobal,
                                      edgeCoord[0]);
                subContVol[1].volume_ =
                    quadrilateralArea(subContVol[1].global,
                                      edgeCoord[1],
                                      elementGlobal,
                                      edgeCoord[2]);
                subContVol[2].volume_ =
                    quadrilateralArea(subContVol[2].global,
                                      edgeCoord[0],
                                      elementGlobal,
                                      edgeCoord[3]);
                subContVol[3].volume_ =
                    quadrilateralArea(subContVol[3].global,
                                      edgeCoord[3],
                                      elementGlobal,
                                      edgeCoord[1]);
                break;
            default:
                throw std::logic_error("Not implemented:VcfvStencil dim = "+std::to_string(dim)
                                       +", numVertices = "+std::to_string(numVertices));
            }
        }
        else if (dim == 3) {
            switch (numVertices) {
            case 4: // 3D, tetrahedron
                for (unsigned k = 0; k < numVertices; k++)
                    subContVol[k].volume_ = elementVolume / 4.0;
                break;
            case 5: // 3D, pyramid
                subContVol[0].volume_ =
                    hexahedronVolume(subContVol[0].global,
                                     edgeCoord[2],
                                     faceCoord[0],
                                     edgeCoord[0],
                                     edgeCoord[4],
                                     faceCoord[3],
                                     elementGlobal,
                                     faceCoord[1]);
                subContVol[1].volume_ =
                    hexahedronVolume(subContVol[1].global,
                                     edgeCoord[1],
                                     faceCoord[0],
                                     edgeCoord[2],
                                     edgeCoord[5],
                                     faceCoord[2],
                                     elementGlobal,
                                     faceCoord[3]);
                subContVol[2].volume_ =
                    hexahedronVolume(subContVol[2].global,
                                     edgeCoord[0],
                                     faceCoord[0],
                                     edgeCoord[3],
                                     edgeCoord[6],
                                     faceCoord[1],
                                     elementGlobal,
                                     faceCoord[4]);
                subContVol[3].volume_ =
                    hexahedronVolume(subContVol[3].global,
                                     edgeCoord[3],
                                     faceCoord[0],
                                     edgeCoord[1],
                                     edgeCoord[7],
                                     faceCoord[4],
                                     elementGlobal,
                                     faceCoord[2]);
                subContVol[4].volume_ = elementVolume -
                    subContVol[0].volume_ -
                    subContVol[1].volume_ -
                    subContVol[2].volume_ -
                    subContVol[3].volume_;
                break;
            case 6: // 3D, prism
                subContVol[0].volume_ =
                    hexahedronVolume(subContVol[0].global,
                                     edgeCoord[3],
                                     faceCoord[3],
                                     edgeCoord[4],
                                     edgeCoord[0],
                                     faceCoord[0],
                                     elementGlobal,
                                     faceCoord[1]);
                subContVol[1].volume_ =
                    hexahedronVolume(subContVol[1].global,
                                     edgeCoord[5],
                                     faceCoord[3],
                                     edgeCoord[3],
                                     edgeCoord[1],
                                     faceCoord[2],
                                     elementGlobal,
                                     faceCoord[0]);
                subContVol[2].volume_ =
                    hexahedronVolume(subContVol[2].global,
                                     edgeCoord[4],
                                     faceCoord[3],
                                     edgeCoord[5],
                                     edgeCoord[2],
                                     faceCoord[1],
                                     elementGlobal,
                                     faceCoord[2]);
                subContVol[3].volume_ =
                    hexahedronVolume(edgeCoord[0],
                                     faceCoord[0],
                                     elementGlobal,
                                     faceCoord[1],
                                     subContVol[3].global,
                                     edgeCoord[6],
                                     faceCoord[4],
                                     edgeCoord[7]);
                subContVol[4].volume_ =
                    hexahedronVolume(edgeCoord[1],
                                     faceCoord[2],
                                     elementGlobal,
                                     faceCoord[0],
                                     subContVol[4].global,
                                     edgeCoord[8],
                                     faceCoord[4],
                                     edgeCoord[6]);
                subContVol[5].volume_ =
                    hexahedronVolume(edgeCoord[2],
                                     faceCoord[1],
                                     elementGlobal,
                                     faceCoord[2],
                                     subContVol[5].global,
                                     edgeCoord[7],
                                     faceCoord[4],
                                     edgeCoord[8]);
                break;
            case 8: // 3D, hexahedron
                subContVol[0].volume_ =
                    hexahedronVolume(subContVol[0].global,
                                     edgeCoord[6],
                                     faceCoord[4],
                                     edgeCoord[4],
                                     edgeCoord[0],
                                     faceCoord[2],
                                     elementGlobal,
                                     faceCoord[0]);
                subContVol[1].volume_ =
                    hexahedronVolume(subContVol[1].global,
                                     edgeCoord[5],
                                     faceCoord[4],
                                     edgeCoord[6],
                                     edgeCoord[1],
                                     faceCoord[1],
                                     elementGlobal,
                                     faceCoord[2]);
                subContVol[2].volume_ =
                    hexahedronVolume(subContVol[2].global,
                                     edgeCoord[4],
                                     faceCoord[4],
                                     edgeCoord[7],
                                     edgeCoord[2],
                                     faceCoord[0],
                                     elementGlobal,
                                     faceCoord[3]);
                subContVol[3].volume_ =
                    hexahedronVolume(subContVol[3].global,
                                     edgeCoord[7],
                                     faceCoord[4],
                                     edgeCoord[5],
                                     edgeCoord[3],
                                     faceCoord[3],
                                     elementGlobal,
                                     faceCoord[1]);
                subContVol[4].volume_ =
                    hexahedronVolume(edgeCoord[0],
                                     faceCoord[2],
                                     elementGlobal,
                                     faceCoord[0],
                                     subContVol[4].global,
                                     edgeCoord[10],
                                     faceCoord[5],
                                     edgeCoord[8]);
                subContVol[5].volume_ =
                    hexahedronVolume(edgeCoord[1],
                                     faceCoord[1],
                                     elementGlobal,
                                     faceCoord[2],
                                     subContVol[5].global,
                                     edgeCoord[9],
                                     faceCoord[5],
                                     edgeCoord[10]);
                subContVol[6].volume_ =
                    hexahedronVolume(edgeCoord[2],
                                     faceCoord[0],
                                     elementGlobal,
                                     faceCoord[3],
                                     subContVol[6].global,
                                     edgeCoord[8],
                                     faceCoord[5],
                                     edgeCoord[11]);
                subContVol[7].volume_ =
                    hexahedronVolume(edgeCoord[3],
                                     faceCoord[3],
                                     elementGlobal,
                                     faceCoord[1],
                                     subContVol[7].global,
                                     edgeCoord[11],
                                     faceCoord[5],
                                     edgeCoord[9]);
                break;
            default:
                throw std::logic_error("Not implemented:VcfvStencil for dim = "+std::to_string(dim)
                                       +", numVertices = "+std::to_string(numVertices));
            }
        }
        else
            throw std::logic_error("Not implemented:VcfvStencil for dim = "+std::to_string(dim));
    }
#if __GNUC__ || __clang__
#pragma GCC diagnostic pop
#endif

    const GridView&     gridView_;
    const Mapper& vertexMapper_;

    Element element_;

#if HAVE_DUNE_LOCALFUNCTIONS
    static LocalFiniteElementCache feCache_;
#endif // HAVE_DUNE_LOCALFUNCTIONS

    //! local coordinate of element center
    LocalPosition elementLocal;
    //! global coordinate of element center
    GlobalPosition elementGlobal;
    //! element volume
    Scalar elementVolume;
    //! data of the sub control volumes
    SubControlVolume subContVol[maxNC];
    //! data of the sub control volume faces
    SubControlVolumeFace subContVolFace[maxNE];
    //! data of the boundary faces
    BoundaryFace boundaryFace_[maxBF];
    unsigned numBoundarySegments_;
    //! global coordinates of the edge centers
    GlobalPosition edgeCoord[maxNE];
    //! global coordinates of the face centers
    GlobalPosition faceCoord[maxNF];
    //! number of verts
    unsigned numVertices;
    //! number of edges
    unsigned numEdges;
    //! number of faces (0 in < 3D)
    unsigned numFaces;
    Dune::GeometryType geometryType_;
};

#if HAVE_DUNE_LOCALFUNCTIONS
template<class Scalar, class GridView>
typename VcfvStencil<Scalar, GridView>::LocalFiniteElementCache
VcfvStencil<Scalar, GridView>::feCache_;
#endif // HAVE_DUNE_LOCALFUNCTIONS

} // namespace Opm

#endif

