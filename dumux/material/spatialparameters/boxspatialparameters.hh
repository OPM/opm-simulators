// $Id$
/*****************************************************************************
 *   Copyright (C) 2010 by Andreas Lauser                                    *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 *
 * \ingroup SpatialParameters
 * \brief The base class for spatial parameters of problems using the
 *        box method.
 */
#ifndef DUMUX_BOX_SPATIAL_PARAMETERS_HH
#define DUMUX_BOX_SPATIAL_PARAMETERS_HH

#include <dumux/common/propertysystem.hh>
#include <dumux/common/math.hh>

#include <dumux/boxmodels/common/boxproperties.hh>

#include <dune/common/fmatrix.hh>

namespace Dumux {
// forward declation of property tags
namespace Properties {
NEW_PROP_TAG(SpatialParameters);
};

/*!
 * \ingroup SpatialParameters
 */


/**
 * \brief The base class for spatial parameters of problems using the
 *        box method.
 */
template<class TypeTag>
class BoxSpatialParameters
{
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(SpatialParameters)) Implementation;

    enum {
        dimWorld = GridView::dimensionworld
    };

    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FVElementGeometry)) FVElementGeometry;

    typedef typename GridView::ctype CoordScalar;
    typedef Dune::FieldMatrix<CoordScalar, dimWorld, dimWorld> Tensor;

public:
    BoxSpatialParameters(const GridView &gv)
    { }

    ~BoxSpatialParameters()
    {}

    /*!
     * \brief Returns the factor by which the volume of a sub control
     *        volume needs to be multiplied in order to get cubic
     *        meters.
     *
     * \param element The current finite element
     * \param fvElemGeom The current finite volume geometry of the element
     * \param scvIdx The index sub-control volume face where the
     *                      factor ought to be calculated.
     *
     * By default that's just 1.0
     */
    Scalar extrusionFactorScv(const Element &element,
                              const FVElementGeometry &fvElemGeom,
                              int scvIdx) const
    { return 1.0; }

    /*!
     * \brief Returns the factor by which the area of a sub control
     *        volume face needs to be multiplied in order to get
     *        square meters.
     *
     * \param element The current finite element
     * \param fvElemGeom The current finite volume geometry of the element
     * \param scvfIdx The index sub-control volume face where the
     *                      factor ought to be calculated.
     *
     * By default it is the arithmetic mean of the extrusion factor of
     * the face's two sub-control volumes.
     */
    Scalar extrusionFactorScvf(const Element &element,
                              const FVElementGeometry &fvElemGeom,
                              int scvfIdx) const
    {
        return
            0.5 *
            (asImp_().extrusionFactorScv(element,
                                         fvElemGeom,
                                         fvElemGeom.subContVolFace[scvfIdx].i)
             +
             asImp_().extrusionFactorScv(element,
                                         fvElemGeom,
                                         fvElemGeom.subContVolFace[scvfIdx].j));
    }

    /*!
     * \brief Averages the intrinsic permeability (Scalar).
     * \param result averaged intrinsic permeability
     * \param K1 intrinsic permeability of the first node
     * \param K2 intrinsic permeability of the second node
     */
    void meanK(Tensor &result,
               Scalar K1,
               Scalar K2) const
    {
        const Scalar K = Dumux::harmonicMean(K1, K2);
        for (int i = 0; i < dimWorld; ++i) {
            for (int j = 0; j < dimWorld; ++j)
                result[i][j] = 0;
            result[i][i] = K;
        }
    }

    /*!
     * \brief Averages the intrinsic permeability (Tensor).
     * \param result averaged intrinsic permeability
     * \param K1 intrinsic permeability of the first node
     * \param K2 intrinsic permeability of the second node
     */
    void meanK(Tensor &result,
               const Tensor &K1,
               const Tensor &K2) const
    {
        // entry-wise harmonic mean. this is almost certainly wrong if
        // you have off-main diagonal entries in your permeabilities!
        for (int i = 0; i < dimWorld; ++i)
            for (int j = 0; j < dimWorld; ++j)
                result[i][j] = harmonicMean(K1[i][j], K2[i][j]);
    }

protected:
    Implementation &asImp_()
    { return *static_cast<Implementation*>(this); }

    const Implementation &asImp_() const
    { return *static_cast<const Implementation*>(this); }
};

} // namespace Dumux

#endif
