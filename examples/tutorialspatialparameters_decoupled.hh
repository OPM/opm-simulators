// $Id: test_2p_spatialparamsinjection.hh 3456 2010-04-09 12:11:51Z mwolff $
/*****************************************************************************
 *   Copyright (C) 2008-2009 by Markus Wolff                                 *
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
#ifndef TUTORIALSPATIALPARAMETERS_DECOUPLED_HH
#define TUTORIALSPATIALPARAMETERS_DECOUPLED_HH


//#include <dumux/new_material/fluidmatrixinteractions/2p/linearmaterial.hh>
#include <dumux/new_material/fluidmatrixinteractions/2p/regularizedbrookscorey.hh>
#include <dumux/new_material/fluidmatrixinteractions/2p/efftoabslaw.hh>

namespace Dumux
{

/** \todo Please doc me! */

template<class TypeTag>
class TutorialSpatialParametersDecoupled
{
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Grid))     Grid;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar))   Scalar;
    typedef typename Grid::ctype                            CoordScalar;

    enum
        {dim=Grid::dimension, dimWorld=Grid::dimensionworld, numEq=1};
    typedef    typename Grid::Traits::template Codim<0>::Entity Element;

    typedef Dune::FieldVector<CoordScalar, dimWorld> GlobalPosition;
    typedef Dune::FieldVector<CoordScalar, dim> LocalPosition;
    typedef Dune::FieldMatrix<Scalar,dim,dim> FieldMatrix;

    typedef RegularizedBrooksCorey<Scalar>                RawMaterialLaw;
//    typedef LinearMaterial<Scalar>                        RawMaterialLaw;
public:
    typedef EffToAbsLaw<RawMaterialLaw>               MaterialLaw;
    typedef typename MaterialLaw::Params MaterialLawParams;

    void update (Scalar saturationW, const Element& element)
    {

    }

    const FieldMatrix& intrinsicPermeability  (const GlobalPosition& globalPos, const Element& element) const
    {
            return K_;
    }

    double porosity(const GlobalPosition& globalPos, const Element& element) const
    {
        return 0.2;
    }


    // return the brooks-corey context depending on the position
    const MaterialLawParams& materialLawParams(const GlobalPosition& globalPos, const Element &element) const
    {
            return materialLawParams_;
    }


    TutorialSpatialParametersDecoupled(const GridView& gridView)
    : K_(0)
    {
        for (int i = 0; i < dim; i++)
                K_[i][i] = 1e-7;

        // residual saturations
        materialLawParams_.setSwr(0);
        materialLawParams_.setSnr(0);

        // parameters for the Brooks-Corey Law
        // entry pressures
        materialLawParams_.setPe(10000);

        // Brooks-Corey shape parameters
        materialLawParams_.setAlpha(2);
    }

private:
    MaterialLawParams materialLawParams_;
    FieldMatrix K_;
};

} // end namespace
#endif
