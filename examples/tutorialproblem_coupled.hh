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
 * \copydoc Dumux::TutorialProblemCoupled
 */
#ifndef DUMUX_TUTORIAL_PROBLEM_COUPLED_HH    // guardian macro /*@\label{tutorial-coupled:guardian1}@*/
#define DUMUX_TUTORIAL_PROBLEM_COUPLED_HH    // guardian macro /*@\label{tutorial-coupled:guardian2}@*/

// The numerical model
#include <dumux/boxmodels/immiscible/immisciblemodel.hh>

// The chemical species that are used
#include <dumux/material/components/h2o.hh>
#include <dumux/material/components/lnapl.hh>

// The material laws
#include <dumux/material/fluidmatrixinteractions/2p/regularizedbrookscorey.hh> /*@\label{tutorial-coupled:rawLawInclude}@*/
#include <dumux/material/fluidmatrixinteractions/2p/efftoabslaw.hh>
#include <dumux/material/fluidmatrixinteractions/mp/2padapter.hh>

// For the DUNE grid
#include <dune/grid/yaspgrid.hh>
#include <dumux/common/cubegridcreator.hh> /*@\label{tutorial-coupled:include-grid-creator}@*/

// For Dune::FieldMatrix
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

// Set the wetting phase /*@\label{tutorial-coupled:2p-system-start}@*/
SET_TYPE_PROP(TutorialProblemCoupled, WettingPhase,   /*@\label{tutorial-coupled:wettingPhase}@*/
              Dumux::LiquidPhase<typename GET_PROP_TYPE(TypeTag, Scalar),
              Dumux::H2O<typename GET_PROP_TYPE(TypeTag, Scalar)> >); 

// Set the non-wetting phase
SET_TYPE_PROP(TutorialProblemCoupled, NonwettingPhase,  /*@\label{tutorial-coupled:nonwettingPhase}@*/
              Dumux::LiquidPhase<typename GET_PROP_TYPE(TypeTag, Scalar), 
              Dumux::LNAPL<typename GET_PROP_TYPE(TypeTag, Scalar)> >);  /*@\label{tutorial-coupled:2p-system-end}@*/

// Set the material law
SET_PROP(TutorialProblemCoupled, MaterialLaw)
{
private:
    // Retrieve the C++ type used to represent scalar values
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    // Select the base material law to be used
    typedef RegularizedBrooksCorey<Scalar> RawMaterialLaw;     /*@\label{tutorial-coupled:rawlaw}@*/
    // Converts absolute saturations into effective ones before
    // passing it to the base material law
    typedef EffToAbsLaw<RawMaterialLaw> TwoPMaterialLaw;   /*@\label{tutorial-coupled:eff2abs}@*/

    // Retrieve the index of the wetting phase
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    enum { wPhaseIdx = FluidSystem::wPhaseIdx };

public:
    // Convert two-phase material law into a general M-phase one.
    typedef TwoPAdapter<wPhaseIdx, TwoPMaterialLaw> type;
};

// Disable gravity
SET_BOOL_PROP(TutorialProblemCoupled, EnableGravity, false); /*@\label{tutorial-coupled:gravity}@*/

// define how long the simulation should run [s]  /*@\label{tutorial-coupled:default-params-begin}@*/
SET_SCALAR_PROP(TutorialProblemCoupled, EndTime, 100e3);

// define the size of the initial time step [s]
SET_SCALAR_PROP(TutorialProblemCoupled, InitialTimeStepSize, 500.0);

// define the physical size of the problem's domain [m]
SET_SCALAR_PROP(TutorialProblemCoupled, DomainSizeX, 300.0);
SET_SCALAR_PROP(TutorialProblemCoupled, DomainSizeY, 60.0);
SET_SCALAR_PROP(TutorialProblemCoupled, DomainSizeZ, 0.0);

// // define the number of cells used for discretizing the physical domain
SET_INT_PROP(TutorialProblemCoupled, CellsX, 100);
SET_INT_PROP(TutorialProblemCoupled, CellsY, 1);
SET_INT_PROP(TutorialProblemCoupled, CellsZ, 1);   /*@\label{tutorial-coupled:default-params-end}@*/
} // namespace Properties

//! Tutorial problem using the fully-implicit immiscible model.
template <class TypeTag>
class TutorialProblemCoupled
    : public GET_PROP_TYPE(TypeTag, BaseProblem) /*@\label{tutorial-coupled:def-problem}@*/
{
    typedef typename GET_PROP_TYPE(TypeTag, BaseProblem) ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;

    // Grid dimension
    enum { dimWorld = GridView::dimensionworld };

    // The type of the intrinsic permeability tensor
    typedef Dune::FieldMatrix<Scalar, dimWorld, dimWorld> DimMatrix;

    // eWoms specific types are specified via the property system
    typedef typename GET_PROP_TYPE(TypeTag, TimeManager) TimeManager;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, RateVector) RateVector;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryRateVector) BoundaryRateVector;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLawParams) MaterialLawParams;    /*@\label{tutorial-coupled:matLawObjectType}@*/

    // phase indices
    enum { numPhases = FluidSystem::numPhases };
    enum { wPhaseIdx = FluidSystem::wPhaseIdx };
    enum { nPhaseIdx = FluidSystem::nPhaseIdx };

    // Indices of the conservation equations
    enum { contiWEqIdx = Indices::conti0EqIdx + wPhaseIdx };
    enum { contiNEqIdx = Indices::conti0EqIdx + nPhaseIdx };

public:
    //! The constructor of the problem
    TutorialProblemCoupled(TimeManager &timeManager)
        : ParentType(timeManager, GET_PROP_TYPE(TypeTag, GridCreator)::grid().leafView())
        , eps_(3e-6)
    {
        // Use an isotropic and homogeneous intrinsic permeability 
        K_ = this->toDimMatrix_(1e-7);
        
        // Parameters of the Brooks-Corey law
        materialParams_.setPe(500.0); // entry pressure [Pa]  /*@\label{tutorial-coupled:setLawParams}@*/
        materialParams_.setLambda(2); // shape parameter

        // Set the residual saturations
        materialParams_.setSwr(0.0);
        materialParams_.setSnr(0.0);
    }

    //! Specifies the problem name. This is used for files generated by the simulation.
    const char *name() const
    { return "tutorial_coupled"; }

    //! Returns the temperature at a given position.
    template <class Context>
    Scalar temperature(const Context &context, int spaceIdx, int timeIdx) const
    { return 283.15; }

    //! Returns the intrinsic permeability tensor [m^2] at a position.
    template <class Context>
    const DimMatrix &intrinsicPermeability(const Context &context, /*@\label{tutorial-coupled:permeability}@*/
                                           int spaceIdx, int timeIdx) const
    { return K_; }
    
    //! Defines the porosity [-] of the medium at a given position
    template <class Context>
    Scalar porosity(const Context &context, int spaceIdx, int timeIdx) const  /*@\label{tutorial-coupled:porosity}@*/
    { return 0.2; }

    //! Returns the parameter object for the material law at a given position
    template <class Context>
    const MaterialLawParams& materialLawParams(const Context &context,            /*@\label{tutorial-coupled:matLawParams}@*/
                                               int spaceIdx, int timeIdx) const
    { return materialParams_; }

    //! Evaluates the boundary conditions.
    template <class Context>
    void boundary(BoundaryRateVector &values, 
                  const Context &context, int spaceIdx, int timeIdx) const
    {
        const auto &pos = context.pos(spaceIdx, timeIdx);
        if (pos[0] < eps_) {
            // Free-flow conditions on left boundary
            const auto &materialParams = this->materialLawParams(context, spaceIdx, timeIdx);
            
            ImmiscibleFluidState<Scalar, FluidSystem> fs;
            Scalar Sw = 1.0;
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

    //! Evaluates the source term for all conserved quantities at a given position
    //! of the domain [kg/(m^3 * s)]. Positive values mean that mass is created.
    template <class Context>
    void source(RateVector &source, const Context &context, int spaceIdx, int timeIdx) const
    {
        source[contiWEqIdx] = 0.0;
        source[contiNEqIdx] = 0.0;
    }

    //! Evaluates the initial value at a given position in the domain.
    template <class Context>
    void initial(PrimaryVariables &values,
                 const Context &context, int spaceIdx, int timeIdx) const
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
} // namespace Dumux

#endif
