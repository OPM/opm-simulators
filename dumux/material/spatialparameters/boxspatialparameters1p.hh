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
#ifndef DUMUX_BOX_SPATIAL_PARAMETERS_ONE_P_HH
#define DUMUX_BOX_SPATIAL_PARAMETERS_ONE_P_HH

#include <dumux/common/propertysystem.hh>
#include <dumux/common/math.hh>

#include <dumux/boxmodels/common/boxproperties.hh>

#include <dune/common/fmatrix.hh>

namespace Dumux {
// forward declation of property tags
namespace Properties {
NEW_PROP_TAG(SpatialParameters);
}

/*!
 * \ingroup SpatialParameters
 */


/**
 * \brief The base class for spatial parameters of problems using the
 *        box method.
 */
template<class TypeTag>
class BoxSpatialParametersOneP
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, SpatialParameters) Implementation;
    typedef typename GET_PROP_TYPE(TypeTag, HeatConductionLawParams) HeatConductionLawParams;

    enum {
        dimWorld = GridView::dimensionworld
    };

    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;

    typedef typename GridView::ctype CoordScalar;
    typedef Dune::FieldMatrix<Scalar, dimWorld, dimWorld> Tensor;
    typedef Dune::FieldVector<CoordScalar,dimWorld> GlobalPosition;
    typedef Dune::FieldVector<Scalar,dimWorld> Vector;

public:
    BoxSpatialParametersOneP(const GridView &gv)
    { }

    ~BoxSpatialParametersOneP()
    {}

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

    template <class Context>
    DUMUX_DEPRECATED_MSG("Old problem API used. Please use context objects for your problem!")
    const Tensor intrinsicPermeability(const Context &context,
                                       int spaceIdx, int timeIdx) const
    {
        return toTensor_(asImp_().intrinsicPermeability(context.element(),
                                                        context.fvElemGeom(timeIdx),
                                                        spaceIdx));
    }

    /*!
     * \brief Function for defining the intrinsic (absolute) permeability.
     *
     * \param element The current element
     * \param fvElemGeom The current finite volume geometry of the element
     * \param scvIdx The index of the sub-control volume.
     * \return the intrinsic permeability
     */
    const Tensor intrinsicPermeability(const Element &element,
                                        const FVElementGeometry &fvElemGeom,
                                        int scvIdx) const
    {
        return toTensor_(asImp_().intrinsicPermeabilityAtPos(element.geometry().center()));
    }

    /*!
     * \brief Function for defining the intrinsic (absolute) permeability.
     *
     * \return intrinsic (absolute) permeability
     * \param globalPos The position of the center of the element
     */
    const Tensor& intrinsicPermeabilityAtPos (const GlobalPosition& globalPos) const
    {
        DUNE_THROW(Dune::InvalidStateException,
                   "The spatial parameters do not provide "
                   "a intrinsicPermeabilityAtPos() method.");
    }

    template <class Context>
    DUMUX_DEPRECATED_MSG("Old problem API used. Please use context objects for your problem!")
    Scalar porosity(const Context &context,
                    int spaceIdx, int timeIdx) const
    {
        return asImp_().porosity(context.element(),
                                 context.fvElemGeom(timeIdx),
                                 spaceIdx);
    }

    /*!
     * \brief Function for defining the porosity.
     *
     * \param element The current element
     * \param fvElemGeom The current finite volume geometry of the element
     * \param scvIdx The index of the sub-control volume.
     * \return porosity
     */
    Scalar porosity(const Element &element,
                    const FVElementGeometry &fvElemGeom,
                    int scvIdx) const
    {
        return asImp_().porosityAtPos(element.geometry().center());
    }

    /*!
     * \brief Function for defining the porosity.
     *
     * \return porosity
     * \param globalPos The position of the center of the element
     */
    Scalar porosityAtPos(const GlobalPosition& globalPos) const
    {
        DUNE_THROW(Dune::InvalidStateException,
                   "The spatial parameters do not provide "
                   "a porosity() method.");
    }

    /*!
     * \brief Returns the heat capacity [J/(K m^3)] of the solid phase
     *        with no pores in the sub-control volume.
     */
    template <class Context>
    DUMUX_DEPRECATED_MSG("Old problem API used. Please use context objects for your problem!")
    Scalar heatCapacitySolid(const Context &context,
                             int spaceIdx, int timeIdx) const
    {
        return asImp_().heatCapacitySolid(context.element(),
                                          context.fvElemGeom(timeIdx),
                                          spaceIdx);
    }

    /*!
     * \brief Returns the heat capacity [J/(K m^3)] of the solid phase
     *        with no pores in the sub-control volume.
     */
    Scalar heatCapacitySolid(const Element &element,
                             const FVElementGeometry &fvElemGeom,
                             int scvIdx) const
    {
        return asImp_().heatCapacitySolidAtPos(element.geometry().center());
    }

    /*!
     * \brief Returns the heat capacity [J/(K m^3)] of the solid phase
     *        with no pores in the sub-control volume.
     */
    Scalar heatCapacitySolidAtPos(const GlobalPosition& globalPos) const
    {
        DUNE_THROW(Dune::InvalidStateException,
                   "The spatial parameters do not provide "
                   "a heatCapacitySolid() method.");
    }

    /*!
     * \brief Returns the thermal conductivity [W / (K m)] of the solid
     *        phase disregarding the pores in a sub-control volume.
     */
    template <class Context>
    const HeatConductionLawParams&
    heatConducionParams(const Context &context, int spaceIdx, int timeIdx) const
    {
        DUNE_THROW(Dune::InvalidStateException,
                   "Spatial parameters does not provide "
                   "a heatConducionParams() method!");
    }

    /*!
     * \brief Returns the thermal conductivity [W / (K m)] of the solid
     *        phase disregarding the pores in a sub-control volume.
     */
    Scalar thermalConductivitySolidAtPos(const GlobalPosition& globalPos) const
    {
        DUNE_THROW(Dune::InvalidStateException,
                   "The spatial parameters do not provide "
                   "a thermalConductivitySolid() method.");
    }

    /*!
     * \brief Used by the 1p2c model.
     */
    template <class Context>
    bool useTwoPointGradient(const Context &context,
                             int spaceIdx,
                             int timeIdx) const
    {
        return false; // use finite element gradient!
    }

private:
    Implementation &asImp_()
    { return *static_cast<Implementation*>(this); }

    const Implementation &asImp_() const
    { return *static_cast<const Implementation*>(this); }

    const Tensor &toTensor_(const Tensor &val) const
    { return val; };

    Tensor toTensor_(Scalar val) const
    {
        Tensor ret(0.0);
        for (int i = 0; i < Tensor::rows; ++i)
            ret[i][i] = val;
        return ret;
    };
};

} // namespace Dumux

#endif
