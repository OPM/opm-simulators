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

#include "boxspatialparameters1p.hh"

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
class BoxSpatialParameters: public BoxSpatialParametersOneP<TypeTag>
{
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(SpatialParameters)) Implementation;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(MaterialLawParams)) MaterialLawParams;

    enum {
        dimWorld = GridView::dimensionworld
    };

    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FVElementGeometry)) FVElementGeometry;

    typedef typename GridView::ctype CoordScalar;
    typedef Dune::FieldVector<CoordScalar,dimWorld> GlobalPosition;

public:
    BoxSpatialParameters(const GridView &gv)
    :BoxSpatialParametersOneP<TypeTag>(gv)
    { }

    /*!
     * \brief Function for defining the parameters needed by constitutive relationships (kr-Sw, pc-Sw, etc.).
     *
     * \return the material parameters object
     * \param element The element
     */
    const MaterialLawParams& materialLawParams(const Element &element,
            const FVElementGeometry &fvElemGeom,
            int scvIdx) const
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
    { return *static_cast<Implementation*>(this); }

    const Implementation &asImp_() const
    { return *static_cast<const Implementation*>(this); }
};

} // namespace Dumux

#endif
