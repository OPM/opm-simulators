// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2008-2009 by Melanie Darcis                               *
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
 * \brief The spatial parameters for the fully coupled tutorial problem
 *        which uses the twophase box model.
 */
#ifndef DUMUX_TUTORIAL_SPATIAL_PARAMETERS_COUPLED_HH
#define DUMUX_TUTORIAL_SPATIAL_PARAMETERS_COUPLED_HH

// include parent spatialparameters
#include <dumux/material/spatialparameters/boxspatialparameters.hh>

// include material laws
#include <dumux/material/fluidmatrixinteractions/2p/regularizedbrookscorey.hh> /*@\label{tutorial-coupled:rawLawInclude}@*/
#include <dumux/material/fluidmatrixinteractions/2p/efftoabslaw.hh>

namespace Dumux
{

//forward declaration
template<class TypeTag>
class TutorialSpatialParametersCoupled;

namespace Properties
{
// The spatial parameters TypeTag
NEW_TYPE_TAG(TutorialSpatialParametersCoupled);/*@\label{tutorial-coupled:define-spatialparameters-typetag}@*/

// Set the spatial parameters
SET_TYPE_PROP(TutorialSpatialParametersCoupled, SpatialParameters, Dumux::TutorialSpatialParametersCoupled<TypeTag>); /*@\label{tutorial-coupled:set-spatialparameters}@*/

// Set the material law
SET_PROP(TutorialSpatialParametersCoupled, MaterialLaw)
{
private:
    // material law typedefs
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    // select material law to be used
    typedef RegularizedBrooksCorey<Scalar> RawMaterialLaw;     /*@\label{tutorial-coupled:rawlaw}@*/
public:
    // adapter for absolute law
    typedef EffToAbsLaw<RawMaterialLaw> type;   /*@\label{tutorial-coupled:eff2abs}@*/
};
}

/*!
 * \ingroup TwoPBoxModel
 *
 * \brief The spatial parameters for the fully coupled tutorial problem
 *        which uses the twophase box model.
 */
template<class TypeTag>
class TutorialSpatialParametersCoupled: public BoxSpatialParameters<TypeTag> /*@\label{tutorial-coupled:tutorialSpatialParameters}@*/
{
    // Get informations for current implementation via property system
    typedef typename GET_PROP_TYPE(TypeTag, Grid) Grid;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    enum
    {
        dim = Grid::dimension,
        dimWorld = Grid::dimensionworld
    };

    // Get object types for function arguments
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename Grid::Traits::template Codim<0>::Entity Element;

public:
    // get material law from property system
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;
    // determine appropriate parameters depening on selected materialLaw
    typedef typename MaterialLaw::Params MaterialLawParams;    /*@\label{tutorial-coupled:matLawObjectType}@*/


    //! Intrinsic permeability tensor K \f$[m^2]\f$ depending
    /*! on the position in the domain
     *
     *  \param element The finite volume element
     *  \param fvElemGeom The finite-volume geometry in the box scheme
     *  \param scvIdx The local vertex index
     *
     *  Alternatively, the function intrinsicPermeabilityAtPos(const GlobalPosition& globalPos) could be defined, where globalPos
     *  is the vector including the global coordinates of the finite volume.
     */
    const Dune::FieldMatrix<Scalar, dim, dim> &intrinsicPermeability(const Element &element, /*@\label{tutorial-coupled:permeability}@*/
                                                    const FVElementGeometry &fvElemGeom,
                                                    int scvIdx) const
    {
        return K_;
    }

    //! Define the porosity \f$[-]\f$ of the porous medium depending
    /*! on the position in the domain
     *
     *  \param element The finite volume element
     *  \param fvElemGeom The finite-volume geometry in the box scheme
     *  \param scvIdx The local vertex index
     *
     *  Alternatively, the function porosityAtPos(const GlobalPosition& globalPos) could be defined, where globalPos
     *  is the vector including the global coordinates of the finite volume.
     */
    double porosity(const Element &element,                    /*@\label{tutorial-coupled:porosity}@*/
                    const FVElementGeometry &fvElemGeom,
                    int scvIdx) const
    {
        return 0.2;
    }

    /*! Return the parameter object for the material law (i.e. Brooks-Corey)
     *  depending on the position in the domain
     *
     *  \param element The finite volume element
     *  \param fvElemGeom The finite-volume geometry in the box scheme
     *  \param scvIdx The local vertex index
     *
     *  Alternatively, the function materialLawParamsAtPos(const GlobalPosition& globalPos) could be defined, where globalPos
     *  is the vector including the global coordinates of the finite volume.
     */
    const MaterialLawParams& materialLawParams(const Element &element,            /*@\label{tutorial-coupled:matLawParams}@*/
                                               const FVElementGeometry &fvElemGeom,
                                               int scvIdx) const
    {
        return materialParams_;
    }

    // constructor
    TutorialSpatialParametersCoupled(const GridView& gridView) :
        BoxSpatialParameters<TypeTag>(gridView),
        K_(0)
    {
        //set main diagonal entries of the permeability tensor to a value
        //setting to one value means: isotropic, homogeneous
        for (int i = 0; i < dim; i++)
            K_[i][i] = 1e-7;

        //set residual saturations
        materialParams_.setSwr(0.0);                /*@\label{tutorial-coupled:setLawParams}@*/
        materialParams_.setSnr(0.0);

        //parameters of Brooks & Corey Law
        materialParams_.setPe(500.0);
        materialParams_.setLambda(2);
    }

private:
    Dune::FieldMatrix<Scalar, dim, dim> K_;
    // Object that holds the values/parameters of the selected material law.
    MaterialLawParams materialParams_;                 /*@\label{tutorial-coupled:matParamsObject}@*/
};
} // end namespace
#endif
