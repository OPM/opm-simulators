// $Id: test_2p_spatialparamsinjection.hh 3456 2010-04-09 12:11:51Z mwolff $
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
#ifndef TUTORIALSPATIALPARAMETERS_DECOUPLED_HH
#define TUTORIALSPATIALPARAMETERS_DECOUPLED_HH


//#include <dumux/material/fluidmatrixinteractions/2p/linearmaterial.hh>
#include <dumux/material/fluidmatrixinteractions/2p/regularizedbrookscorey.hh>
#include <dumux/material/fluidmatrixinteractions/2p/efftoabslaw.hh>

namespace Dumux
{
//! Definition of the spatial parameters for the decoupled tutorial

template<class TypeTag>
class TutorialSpatialParametersDecoupled
{
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

    // material law typedefs
    typedef RegularizedBrooksCorey<Scalar> EffectiveMaterialLaw;
//    typedef LinearMaterial<Scalar> EffectiveMaterialLaw;
public:
    typedef EffToAbsLaw<EffectiveMaterialLaw> MaterialLaw;
    typedef typename MaterialLaw::Params MaterialLawParams;

    //! Update the spatial parameters with the flow solution after a timestep.
    /*! Function left blank as there is nothing to do for the tutorial.
     */
    void update (Scalar saturationW, const Element& element)
    {    }
    //! Intrinsic permeability tensor
    /*! Apply the intrinsic permeability tensor \f$[m^2]\f$ to a
     *  pressure potential gradient.
     */
    const FieldMatrix& intrinsicPermeability (const GlobalPosition& globalPos,
                                              const Element& element) const
    {
            return K_;
    }

    //! Define the porosity \f$[-]\f$ of the spatial parameters
    double porosity(const GlobalPosition& globalPos, const Element& element) const
    {
        return 0.2;
    }

    //! return the parameter object for the material law (i.e. Brooks-Corey)
    //! which may vary with the spatial position
    const MaterialLawParams& materialLawParams(const GlobalPosition& globalPos,
                                               const Element &element) const
    {
            return materialLawParams_;
    }

    //! Constructor
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
        materialLawParams_.setPe(500);

        // Brooks-Corey shape parameters
        materialLawParams_.setAlpha(2);
    }

private:
    MaterialLawParams materialLawParams_;
    FieldMatrix K_;
};

} // end namespace
#endif
