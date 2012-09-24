// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2009-2012 by Andreas Lauser                               *
 *   Copyright (C) 2010-2012 by Klaus Mosthaf                                *
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
 * \brief Tutorial problem for a fully coupled twophase box model.
 */
#ifndef DUMUX_TUTORIAL_PROBLEM_COUPLED_HH    // guardian macro /*@\label{tutorial-coupled:guardian1}@*/
#define DUMUX_TUTORIAL_PROBLEM_COUPLED_HH    // guardian macro /*@\label{tutorial-coupled:guardian2}@*/

// The numerical model
#include <dumux/boxmodels/immiscible/immisciblemodel.hh>

// The components that are used
#include <dumux/material/components/h2o.hh>
#include <dumux/material/components/lnapl.hh>

// The material laws
#include <dumux/material/fluidmatrixinteractions/2p/regularizedbrookscorey.hh> /*@\label{tutorial-coupled:rawLawInclude}@*/
#include <dumux/material/fluidmatrixinteractions/2p/efftoabslaw.hh>
#include <dumux/material/fluidmatrixinteractions/mp/2padapter.hh>

// The DUNE grid used
#include <dune/grid/yaspgrid.hh>
#include <dumux/common/cubegridcreator.hh>

// Dune::FieldVector and Dune::FieldMatrix
#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

namespace Dumux {

// forward declaration of the problem class
template <class TypeTag>
class TutorialProblemCoupled;

namespace Properties {
// Create a new type tag for the problem
NEW_TYPE_TAG(TutorialProblemCoupled, INHERITS_FROM(BoxImmiscibleTwoPhase)); /*@\label{tutorial-coupled:create-type-tag}@*/

// Set the "Problem" property
SET_PROP(TutorialProblemCoupled, Problem) /*@\label{tutorial-coupled:set-problem}@*/
{ typedef Dumux::TutorialProblemCoupled<TypeTag> type;};

// Set grid and the grid creator to be used
SET_TYPE_PROP(TutorialProblemCoupled, Grid, Dune::YaspGrid</*dim=*/2>); /*@\label{tutorial-coupled:set-grid}@*/
SET_TYPE_PROP(TutorialProblemCoupled, GridCreator, Dumux::CubeGridCreator<TypeTag>); /*@\label{tutorial-coupled:set-gridcreator}@*/

// Set the wetting phase
SET_PROP(TutorialProblemCoupled, WettingPhase) /*@\label{tutorial-coupled:2p-system-start}@*/
{
private: typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
public: typedef Dumux::LiquidPhase<Scalar, Dumux::H2O<Scalar> > type; /*@\label{tutorial-coupled:wettingPhase}@*/
};

// Set the non-wetting phase
SET_PROP(TutorialProblemCoupled, NonwettingPhase)
{
private: typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
public: typedef Dumux::LiquidPhase<Scalar, Dumux::LNAPL<Scalar> > type; /*@\label{tutorial-coupled:nonwettingPhase}@*/
}; /*@\label{tutorial-coupled:2p-system-end}@*/

// Set the material law
SET_PROP(TutorialProblemCoupled, MaterialLaw)
{
private:
    // material law typedefs
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    // select material law to be used
    typedef RegularizedBrooksCorey<Scalar> RawMaterialLaw;     /*@\label{tutorial-coupled:rawlaw}@*/
    // adapter for absolute law
    typedef EffToAbsLaw<RawMaterialLaw> TwoPMaterialLaw;   /*@\label{tutorial-coupled:eff2abs}@*/

    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    enum { wPhaseIdx = FluidSystem::wPhaseIdx };

public:
    typedef TwoPAdapter<wPhaseIdx, TwoPMaterialLaw> type;
};

// Disable gravity
SET_BOOL_PROP(TutorialProblemCoupled, EnableGravity, false); /*@\label{tutorial-coupled:gravity}@*/

// define how long the simulation should run [s]
SET_SCALAR_PROP(TutorialProblemCoupled, EndTime, 100e3);

// define the size of the initial time step [s]
SET_SCALAR_PROP(TutorialProblemCoupled, InitialTimeStepSize, 500.0);

// define the properties required by the cube grid creator
SET_SCALAR_PROP(TutorialProblemCoupled, DomainSizeX, 300.0);
SET_SCALAR_PROP(TutorialProblemCoupled, DomainSizeY, 60.0);
SET_SCALAR_PROP(TutorialProblemCoupled, DomainSizeZ, 0.0);

SET_INT_PROP(TutorialProblemCoupled, CellsX, 100);
SET_INT_PROP(TutorialProblemCoupled, CellsY, 1);
SET_INT_PROP(TutorialProblemCoupled, CellsZ, 0);
}

/*!
 * \ingroup TwoPBoxModel
 *
 * \brief  Tutorial problem for a fully coupled twophase box model.
 */
template <class TypeTag>
class TutorialProblemCoupled
    : public GET_PROP_TYPE(TypeTag, BaseProblem) /*@\label{tutorial-coupled:def-problem}@*/
{
    typedef typename GET_PROP_TYPE(TypeTag, BaseProblem) ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;

    // Grid dimension
    enum { dimWorld = GridView::dimensionworld };

    typedef typename GridView::ctype CoordScalar;
    typedef Dune::FieldVector<CoordScalar, dimWorld> GlobalPosition;
    typedef Dune::FieldMatrix<Scalar, dimWorld, dimWorld> DimMatrix;

    // Dumux specific types
    typedef typename GET_PROP_TYPE(TypeTag, TimeManager) TimeManager;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, RateVector) RateVector;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryRateVector) BoundaryRateVector;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    
    // get material law from property system
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;
    // determine type of the parameter objects depening on selected material law
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLawParams) MaterialLawParams;    /*@\label{tutorial-coupled:matLawObjectType}@*/

    // phase indices
    enum { numPhases = FluidSystem::numPhases };
    enum { wPhaseIdx = FluidSystem::wPhaseIdx };
    enum { nPhaseIdx = FluidSystem::nPhaseIdx };

    // indices of the conservation equations
    enum { contiWEqIdx = Indices::conti0EqIdx + wPhaseIdx };
    enum { contiNEqIdx = Indices::conti0EqIdx + nPhaseIdx };

public:
    //! The constructor of the problem
    TutorialProblemCoupled(TimeManager &timeManager)
        : ParentType(timeManager, GET_PROP_TYPE(TypeTag, GridCreator)::grid().leafView())
        , eps_(3e-6)
    {
        // set main diagonal entries of the permeability tensor to a value
        // setting to a single value means: isotropic, homogeneous
        K_ = this->toDimMatrix_(1e-7);

        //set residual saturations
        materialParams_.setSwr(0.0);                /*@\label{tutorial-coupled:setLawParams}@*/
        materialParams_.setSnr(0.0);

        //parameters of Brooks & Corey Law
        materialParams_.setPe(500.0);
        materialParams_.setLambda(2);
    }

    //! Specifies the problem name. This is used as a prefix for files
    //! generated by the simulation.
    const char *name() const
    { return "tutorial_coupled"; }

    //! Returns true if a restart file should be written.
    bool shouldWriteRestartFile() const /*@\label{tutorial-coupled:restart}@*/
    { return false; }

    //! Returns true if the current solution should be written to disk
    //! as a VTK file
    bool shouldWriteOutput() const /*@\label{tutorial-coupled:output}@*/
    { 
        return (this->timeManager().timeStepIndex() % 5 == 0) ||
            this->timeManager().willBeFinished() ;
    }

    //! Returns the temperature within a finite volume. We use constant
    //! 10 degrees Celsius.
    template <class Context>
    Scalar temperature(const Context &context, int spaceIdx, int timeIdx) const
    { return 283.15; }

    //! Returns the intrinsic permeability tensor K \f$[m^2]\f$
    //!  depending on the position in the domain.
    template <class Context>
    const DimMatrix &intrinsicPermeability(const Context &context, /*@\label{tutorial-coupled:permeability}@*/
                                           int spaceIdx, int timeIdx) const
    { return K_; }
    
    //! Defines the porosity \f$[-]\f$ of the porous medium depending
    //!  on the position in the domain.
    template <class Context>
    Scalar porosity(const Context &context,                    /*@\label{tutorial-coupled:porosity}@*/
                    int spaceIdx, int timeIdx) const
    { return 0.2; }

    //! Returns the parameter object for the material law (i.e. Brooks-Corey)
    //! depending on the position in the domain
    template <class Context>
    const MaterialLawParams& materialLawParams(const Context &context,            /*@\label{tutorial-coupled:matLawParams}@*/
                                               int spaceIdx, int timeIdx) const
    { return materialParams_; }

    //! Evaluates the boundary conditions.
    template <class Context>
    void boundary(BoundaryRateVector &values, 
                  const Context &context, int spaceIdx, int timeIdx) const
    {
        const GlobalPosition &pos = context.pos(spaceIdx, timeIdx);
        if (pos[0] < eps_) {
            // Free-flow conditions on left boundary
            const auto &materialParams = this->materialLawParams(context, spaceIdx, timeIdx);
            
            Scalar Sw = 1.0;
            ImmiscibleFluidState<Scalar, FluidSystem> fs;
            fs.setSaturation(wPhaseIdx, Sw);
            fs.setSaturation(nPhaseIdx, 1.0 - Sw);
            fs.setTemperature(temperature(context, spaceIdx, timeIdx));
            
            Scalar pC[numPhases];
            MaterialLaw::capillaryPressures(pC, materialParams, fs);
            fs.setPressure(wPhaseIdx, 200e3);
            fs.setPressure(nPhaseIdx, 200e3 + pC[nPhaseIdx] - pC[nPhaseIdx]);

            values.setFreeFlow(context, spaceIdx, timeIdx, fs);
        }
        else if (pos[0] > this->bboxMax()[0] - eps_) {
            // forced outflow at the right boundary
            RateVector massRate(0.0);

            massRate[contiWEqIdx] = 0.0; // [kg / (s m^2)]
            massRate[contiNEqIdx] = 3e-2; // [kg / (s m^2)]

            values.setMassRate(massRate);
        }
        else // no flow at the remaining boundaries
            values.setNoFlow();
    }

    //! Evaluates the source term for all conserved quantities at a
    //! given position in the pysical domain [(m^3 * s)]. Positive
    //! values mean that mass is created.
    template <class Context>
    void source(RateVector &values, const Context &context, int spaceIdx, int timeIdx) const
    {
        values[contiWEqIdx] = 0.0;
        values[contiNEqIdx]= 0.0;
    }

    //! Evaluates the initial value at a given position in the domain.
    template <class Context>
    void initial(PrimaryVariables &values, const Context &context, int spaceIdx, int timeIdx) const
    {
        ImmiscibleFluidState<Scalar, FluidSystem> fs;

        // the domain is initially fully saturated by LNAPL
        Scalar Sw = 0.0;
        fs.setSaturation(wPhaseIdx, Sw);
        fs.setSaturation(nPhaseIdx, 1.0 - Sw);

        // the temperature is given by the temperature() method
        fs.setTemperature(temperature(context, spaceIdx, timeIdx));
        
        // set pressure of the wetting phase to 200 kPa = 2 bar
        Scalar pC[numPhases];
        MaterialLaw::capillaryPressures(pC, materialLawParams(context, spaceIdx, timeIdx), fs);
        fs.setPressure(wPhaseIdx, 200e3);
        fs.setPressure(nPhaseIdx, 200e3 + pC[nPhaseIdx] - pC[nPhaseIdx]);
        
        values.assignNaive(fs);
    }

private:
    DimMatrix K_;
    // Object that holds the values/parameters of the selected material law.
    MaterialLawParams materialParams_;                 /*@\label{tutorial-coupled:matParamsObject}@*/

    // small epsilon value
    Scalar eps_;
};
}

#endif
