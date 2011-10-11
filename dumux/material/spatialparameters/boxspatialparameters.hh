// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2010 by Andreas Lauser                                    *
 *   Institute for Modelling Hydraulic and Environmental Systems             *
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

#include "boxspatialparameters1p.hh"

namespace Dumux {
// forward declation of property tags
namespace Properties {
NEW_PROP_TAG(SpatialParameters);
NEW_PROP_TAG(MaterialLaw);
NEW_PROP_TAG(MaterialLawParams);
}

/*!
 * \ingroup SpatialParameters
 */


/**
 * \brief The base class for spatial parameters of problems using the
 *        box method.
 */
template<class TypeTag>
class BoxSpatialParameters: public BoxSpatialParametersOneP<TypeTag>
{
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, SpatialParameters) Implementation;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;

    typedef typename GET_PROP_TYPE(TypeTag, MaterialLawParams) MaterialLawParams;

    enum {
        dimWorld = GridView::dimensionworld
    };

    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, FluxVariables) FluxVariables;

    typedef typename GridView::ctype CoordScalar;
    typedef Dune::FieldVector<CoordScalar,dimWorld> GlobalPosition;
    typedef Dune::FieldVector<Scalar,dimWorld> Vector;

public:
    BoxSpatialParameters(const GridView &gv)
        : BoxSpatialParametersOneP<TypeTag>(gv)
    { }

    /*!
     * \brief Function for defining the parameters needed by constitutive relationships (kr-Sw, pc-Sw, etc.).
     *
     * \param element The current element
     * \param fvElemGeom The current finite volume geometry of the element
     * \param scvIdx The index of the sub-control volume.
     * \return the material parameters object
     */
    template <class Context>
    DUMUX_DEPRECATED_MSG("Old problem API used. Please use context objects for your problem!")
    const MaterialLawParams& materialLawParams(const Context &context, int localIdx) const
    {
        return asImp_().materialLawParams(context.element(),
                                          context.fvElemGeom(),
                                          localIdx);
    }

    const MaterialLawParams& materialLawParams(const Element &element,
                                               const FVElementGeometry &fvElemGeom,
                                               int localIdx) const
    {
        return asImp_().materialLawParamsAtPos(fvElemGeom.subContVol[localIdx].global);
    }

    /*!
     * \brief Function for defining the parameters needed by constitutive relationships (kr-Sw, pc-Sw, etc.).
     *
     * \param globalPos The position of the center of the element
     * \return the material parameters object
     */
    const MaterialLawParams& materialLawParamsAtPos(const GlobalPosition& globalPos) const
    {
        DUNE_THROW(Dune::InvalidStateException,
                   "The spatial parameters do not provide "
                   "a materialLawParams() method.");
    }

    /*!
     * \brief Returns the heat capacity [J/(kg K)] of the pure solid.
     */
    template <class Context>
    DUMUX_DEPRECATED_MSG("Old problem API used. Please use context objects for your problem!")
    Scalar heatCapacity(const Context &context, int localIdx) const
    {
        return asImp_().heatCapacity(context.element(),
                                          context.fvElemGeom(),
                                          localIdx);
    }

    Scalar heatCapacity(const Element &element,
                        const FVElementGeometry &fvElemGeom,
                        int localIdx) const
    {
        DUNE_THROW(Dune::InvalidStateException,
                   "The spatial parameters do not provide "
                   "a heatCapacity() method.");
    }

private:
    Implementation &asImp_()
    { return *static_cast<Implementation*>(this); }

    const Implementation &asImp_() const
    { return *static_cast<const Implementation*>(this); }
};

} // namespace Dumux

#endif
