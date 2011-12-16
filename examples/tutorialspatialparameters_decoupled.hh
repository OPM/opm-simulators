// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2008-2009 by Markus Wolff                                 *
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
 * \brief spatial parameters for the sequential tutorial
 */
#ifndef DUMUX_TUTORIAL_SPATIAL_PARAMETERS_DECOUPLED_HH
#define DUMUX_TUTORIAL_SPATIAL_PARAMETERS_DECOUPLED_HH


#include <dumux/material/spatialparameters/fvspatialparameters.hh>
#include <dumux/material/fluidmatrixinteractions/2p/linearmaterial.hh>
#include <dumux/material/fluidmatrixinteractions/2p/regularizedbrookscorey.hh>
#include <dumux/material/fluidmatrixinteractions/2p/efftoabslaw.hh>

namespace Dumux
{

//forward declaration
template<class TypeTag>
class TutorialSpatialParametersDecoupled;

namespace Properties
{
// The spatial parameters TypeTag
NEW_TYPE_TAG(TutorialSpatialParametersDecoupled);

// Set the spatial parameters
SET_TYPE_PROP(TutorialSpatialParametersDecoupled, SpatialParameters, Dumux::TutorialSpatialParametersDecoupled<TypeTag>); /*@\label{tutorial-decoupled:set-spatialparameters}@*/

// Set the material law
SET_PROP(TutorialSpatialParametersDecoupled, MaterialLaw)
{
private:
    // material law typedefs
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef RegularizedBrooksCorey<Scalar> RawMaterialLaw;
public:
    typedef EffToAbsLaw<RawMaterialLaw> type;
};
}

//! Definition of the spatial parameters for the decoupled tutorial

template<class TypeTag>
class TutorialSpatialParametersDecoupled: public FVSpatialParameters<TypeTag>
{
    typedef FVSpatialParameters<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Grid)) Grid;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef typename Grid::ctype CoordScalar;

    enum
        {dim=Grid::dimension, dimWorld=Grid::dimensionworld, numEq=1};
    typedef typename Grid::Traits::template Codim<0>::Entity Element;

    typedef Dune::FieldVector<CoordScalar, dimWorld> GlobalPosition;
    typedef Dune::FieldVector<CoordScalar, dim> LocalPosition;
    typedef Dune::FieldMatrix<Scalar,dim,dim> FieldMatrix;

public:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(MaterialLaw)) MaterialLaw;
    typedef typename MaterialLaw::Params MaterialLawParams;

    //! Intrinsic permeability tensor K \f$[m^2]\f$ depending
    /*! on the position in the domain
     *
     *  \param element The finite volume element
     *
     *  Alternatively, the function intrinsicPermeabilityAtPos(const GlobalPosition& globalPos) could be defined, where globalPos
     *  is the vector including the global coordinates of the finite volume.
     */
    const FieldMatrix& intrinsicPermeability (const Element& element) const
    {
            return K_;
    }

    //! Define the porosity \f$[-]\f$ of the porous medium depending
    /*! on the position in the domain
     *
     *  \param element The finite volume element
     *
     *  Alternatively, the function porosityAtPos(const GlobalPosition& globalPos) could be defined, where globalPos
     *  is the vector including the global coordinates of the finite volume.
     */
    double porosity(const Element& element) const
    {
        return 0.2;
    }

    /*! Return the parameter object for the material law (i.e. Brooks-Corey)
     *  depending on the position in the domain
     *
     *  \param element The finite volume element
     *
     *  Alternatively, the function materialLawParamsAtPos(const GlobalPosition& globalPos) could be defined, where globalPos
     *  is the vector including the global coordinates of the finite volume.
     */
    const MaterialLawParams& materialLawParams(const Element &element) const
    {
            return materialLawParams_;
    }

    //! Constructor
    TutorialSpatialParametersDecoupled(const GridView& gridView)
    : ParentType(gridView), K_(0)
    {
        for (int i = 0; i < dim; i++)
                K_[i][i] = 1e-7;

        // residual saturations
        materialLawParams_.setSwr(0);
        materialLawParams_.setSnr(0);

        // parameters for the Brooks-Corey Law
        // entry pressures
        materialLawParams_.setPe(500);

        // Brooks-Corey shape parameters
        materialLawParams_.setLambda(2);
    }

private:
    MaterialLawParams materialLawParams_;
    FieldMatrix K_;
};

} // end namespace
#endif
