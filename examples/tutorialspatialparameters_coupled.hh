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
#include <dumux/new_material/spatialparameters/boxspatialparameters.hh>

// include material laws
#include <dumux/new_material/fluidmatrixinteractions/2p/linearmaterial.hh>
#include <dumux/new_material/fluidmatrixinteractions/2p/efftoabslaw.hh>

namespace Dumux
{

template<class TypeTag>
class TutorialSpatialParameters: public BoxSpatialParameters<TypeTag> /*@\label{tutorial-coupled:tutorialSpatialParameters}@*/
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

	// select materialLaw to be used
	//TODO: enter materialLaw stuff here









	// method returning the intrinsic permeability tensor K depending
	// on the position within the domain
	const Dune::FieldMatrix<Scalar, dim, dim> &intrinsicPermeability(const Element &element,
													const FVElementGeometry &fvElemGeom,
													int scvIdx) const
	{
		return K_;
	}

	// method returning the porosity of the porous matrix depending on
	// the position within the domain
	double porosity(const Element &element,
					const FVElementGeometry &fvElemGeom,
					int scvIdx) const
	{
		return 0.2;
	}

	// return the materialLaw context (i.e. BC, regularizedVG, etc) depending on the position
	const MaterialLawParams& materialLawParams(const Element &element,
											const FVElementGeometry &fvElemGeom,
											int scvIdx) const
	{
		//TODO: return current material law Object and delete the relPermFlag stuff!


	}




    // method returning the residual saturation of the wetting fluid
    // depending on the position within the domain and on the
    // temperature
    double Sr_w(const Dune::FieldVector<Scalar,dim>& globalPos, const Element& element, /*@\label{tutorial-coupled:srw}@*/
                const Dune::FieldVector<Scalar,dim>& localPos, const double T = 283.15) const
    {
        return 0;
    }

    // method returning the residual saturation of the non-wetting
    // fluid depending on the position within the domain and on the
    // temperature
    double Sr_n(const Dune::FieldVector<Scalar,dim>& globalPos, const Element& element, /*@\label{tutorial-coupled:srn}@*/
                const Dune::FieldVector<Scalar,dim>& localPos, const double T = 283.15) const
    {
        return 0;
    }

    // method returning the parameters of the capillary pressure and
    // the relative permeability functionms depending on the position
    // within the domain and on the temperature
    std::vector<double> paramRelPerm(const Dune::FieldVector<Scalar,dim>& globalPos, const Element& element, /*@\label{tutorial-coupled:parameters}@*/
                                     const Dune::FieldVector<Scalar,dim>& localPos, const double T = 283.15) const
    {
        std::vector<double> param(2);

        //linear law parameters
        param[0] = 0; // minimal capillary pressure
        param[1] = 0; // maximal capillary pressure

        //Brooks-Corey parameters
        //        param[0] = 2; // lambda
        //        param[1] = 0.; // entry-pressure

        return param;
    }

    // method returning the kind of relation used for the calculation
    // of the capillary pressure and the relative permeabilities
    // depending on the position within the domain
    typename Matrix2p<Grid,Scalar>::modelFlag relPermFlag(const Dune::FieldVector<Scalar,dim>& globalPos, const Element& element, /*@\label{tutorial-coupled:flags}@*/
                                                   const Dune::FieldVector<Scalar,dim>& localPos) const
    {
        return Matrix2p<Grid,Scalar>::linear; //flag types defined in
    }                                   //dumux/material/property_baseclasses.hh


	// constructor
	TutorialSpatialParameters(const GridView& gridView) :
		BoxSpatialParameters<TypeTag>(gridView), K_(0)
	{
		for (int i = 0; i < dim; i++)
			K_[i][i] = 1e-7;

		//TODO: set the actual values for the respective parameters that are
		// of interest in the applied materialLaw. e.g. in a linar law, these are
		// residual saturations, a minimum value and a maximum value.
		// Afterwards, please delete the paramRelPerm and Sr_n, Sr_w functions above.








	}

private:
	Dune::FieldMatrix<Scalar, dim, dim> K_;
	// Object that helds the values/parameters of the selected material law.
	//TODO: add something here!
};
} // end namespace
#endif
