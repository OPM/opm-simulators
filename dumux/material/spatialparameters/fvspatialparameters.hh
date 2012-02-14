// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2011 by Markus Wolff                                      *
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
 *        fv method.
 */
#ifndef DUMUX_FV_SPATIAL_PARAMETERS_HH
#define DUMUX_FV_SPATIAL_PARAMETERS_HH

#include "fvspatialparameters1p.hh"

namespace Dumux
{
// forward declation of property tags
namespace Properties {
NEW_PROP_TAG(MaterialLaw);
}

/*!
 * \ingroup SpatialParameters
 */

/**
 * \brief The base class for spatial parameters of a multi-phase problem using the
 *        fv method.
 */
template<class TypeTag>
class FVSpatialParameters: public FVSpatialParametersOneP<TypeTag>
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, SpatialParameters) Implementation;

    enum
    {
        dimWorld = GridView::dimensionworld
    };

    typedef typename GridView::template Codim<0>::Entity Element;
    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;
    /// @cond 0
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw)::Params MaterialLawParams;
    /// @endcond

public:
    FVSpatialParameters(const GridView &gv)
    :FVSpatialParametersOneP<TypeTag>(gv)
    {
    }

    /*!
     * \brief Function for defining the parameters needed by constitutive relationships (kr-Sw, pc-Sw, etc.).
     *
     * \return the material parameters object
     * \param element The element
     */
    const MaterialLawParams& materialLawParams(const Element &element) const
    {
            return asImp_().materialLawParamsAtPos(element.geometry().center());
    }

    /*!
     * \brief Function for defining the parameters needed by constitutive relationships (kr-Sw, pc-Sw, etc.).
     *
     * \return the material parameters object
     * \param globalPos The position of the center of the element
     */
    const MaterialLawParams& materialLawParamsAtPos(const GlobalPosition& globalPos) const
    {
        DUNE_THROW(Dune::InvalidStateException,
                   "The spatial parameters do not provide "
                   "a materialLawParamsAtPos() method.");
    }

private:
    Implementation &asImp_()
    {
        return *static_cast<Implementation*> (this);
    }

    const Implementation &asImp_() const
    {
        return *static_cast<const Implementation*> (this);
    }
};

} // namespace Dumux

#endif
