// $Id$
/*****************************************************************************
 *   Copyright (C) 2008-2009 by Melanie Darcis                               *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 2 of the License, or       *
 *   (at your option) any later version, as long as this copyright notice    *
 *   is included in its original form.                                       *
 *                                                                           *
 *   This program is distributed WITHOUT ANY WARRANTY.                       *
 *****************************************************************************/
#ifndef TUTORIALSPATIALPARAMETERS_COUPLED_HH
#define TUTORIALSPATIALPARAMETERS_COUPLED_HH

// include parent spatialparameters
#include <dumux/material/spatialparameters/boxspatialparameters.hh>

// include material laws
#include <dumux/material/fluidmatrixinteractions/2p/linearmaterial.hh> /*@\label{tutorial-coupled:rawLawInclude}@*/
#include <dumux/material/fluidmatrixinteractions/2p/efftoabslaw.hh>

namespace Dumux
{

template<class TypeTag>
class TutorialSpatialParametersCoupled: public BoxSpatialParameters<TypeTag> /*@\label{tutorial-coupled:tutorialSpatialParameters}@*/
{
    // Get informations for current implementation via property system
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Grid)) Grid;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;

    enum
    {
        dim = Grid::dimension,
        dimWorld = Grid::dimensionworld,
    };

    // Get object types for function arguments
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FVElementGeometry)) FVElementGeometry;
    typedef typename Grid::Traits::template Codim<0>::Entity Element;

    // select material law to be used
    typedef LinearMaterial<Scalar>        RawMaterialLaw;     /*@\label{tutorial-coupled:rawlaw}@*/

public:
    // adapter for absolute law
    typedef EffToAbsLaw<RawMaterialLaw> MaterialLaw;        /*@\label{tutorial-coupled:eff2abs}@*/
    // determine appropriate parameters depening on selected materialLaw
    typedef typename MaterialLaw::Params MaterialLawParams;    /*@\label{tutorial-coupled:matLawObjectType}@*/


    // method returning the intrinsic permeability tensor K depending
    // on the position within the domain
    const Dune::FieldMatrix<Scalar, dim, dim> &intrinsicPermeability(const Element &element, /*@\label{tutorial-coupled:permeability}@*/
                                                    const FVElementGeometry &fvElemGeom,
                                                    int scvIdx) const
    {
        return K_;
    }

    // method returning the porosity of the porous matrix depending on
    // the position within the domain
    double porosity(const Element &element,                    /*@\label{tutorial-coupled:porosity}@*/
                    const FVElementGeometry &fvElemGeom,
                    int scvIdx) const
    {
        return 0.2;
    }

    // return the materialLaw context (i.e. BC, regularizedVG, etc) depending on the position
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
        for (int i = 0; i < dim; i++)
            K_[i][i] = 1e-7;

        //set residual saturations
        materialParams_.setSwr(0.0);                /*@\label{tutorial-coupled:setLawParams}@*/
        materialParams_.setSnr(0.0);

        //linear material law
        materialParams_.setEntryPC(0.0);
        materialParams_.setMaxPC(0.0);
    }

private:
    Dune::FieldMatrix<Scalar, dim, dim> K_;
    // Object that helds the values/parameters of the selected material law.
    MaterialLawParams materialParams_;                 /*@\label{tutorial-coupled:matParamsObject}@*/
};
} // end namespace
#endif
