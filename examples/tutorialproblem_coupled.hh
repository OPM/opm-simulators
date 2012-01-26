// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2008-2009 by Melanie Darcis, Klaus Mosthaf                *
 *   Copyright (C) 2009-2012 by Andreas Lauser                               *
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
 * \brief Tutorial problem for a fully coupled twophase box model.
 */
#ifndef DUMUX_TUTORIAL_PROBLEM_COUPLED_HH    // guardian macro /*@\label{tutorial-coupled:guardian1}@*/
#define DUMUX_TUTORIAL_PROBLEM_COUPLED_HH    // guardian macro /*@\label{tutorial-coupled:guardian2}@*/

// The numerical model
#include <dumux/boxmodels/2p/2pmodel.hh>

// The DUNE grid used
#include <dune/grid/yaspgrid.hh>

// Spatially dependent parameters
// The components that are used
#include <dumux/material/components/h2o.hh>
#include <dumux/material/components/oil.hh>
#include <dumux/common/cubegridcreator.hh>

// include material laws
#include <dumux/material/fluidmatrixinteractions/2p/regularizedbrookscorey.hh> /*@\label{tutorial-coupled:rawLawInclude}@*/
#include <dumux/material/fluidmatrixinteractions/2p/efftoabslaw.hh>
#include <dumux/material/fluidmatrixinteractions/mp/2padapter.hh>

namespace Dumux {

// forward declaration of the problem class
template <class TypeTag>
class TutorialProblemCoupled;

namespace Properties {
// Create a new type tag for the problem
NEW_TYPE_TAG(TutorialProblemCoupled, INHERITS_FROM(BoxTwoP)); /*@\label{tutorial-coupled:create-type-tag}@*/

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
public: typedef Dumux::LiquidPhase<Scalar, Dumux::Oil<Scalar> > type; /*@\label{tutorial-coupled:nonwettingPhase}@*/
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

// define the properties required by the cube grid creator
SET_SCALAR_PROP(TutorialProblemCoupled, GridSizeX, 300.0);
SET_SCALAR_PROP(TutorialProblemCoupled, GridSizeY, 60.0);
SET_SCALAR_PROP(TutorialProblemCoupled, GridSizeZ, 0.0);

SET_INT_PROP(TutorialProblemCoupled, GridCellsX, 100);
SET_INT_PROP(TutorialProblemCoupled, GridCellsY, 1);
SET_INT_PROP(TutorialProblemCoupled, GridCellsZ, 0);
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
    enum { dim = GridView::dimension };

    typedef Dune::FieldVector<Scalar, dim> GlobalPosition;
    typedef Dune::FieldMatrix<Scalar, dim, dim> Tensor;

    // Dumux specific types
    typedef typename GET_PROP_TYPE(TypeTag, TimeManager) TimeManager;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, RateVector) RateVector;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryTypes) BoundaryTypes;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, TwoPIndices) Indices;
    
    // get material law from property system
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;
    // determine type of the parameter objects depening on selected material law
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLawParams) MaterialLawParams;    /*@\label{tutorial-coupled:matLawObjectType}@*/

    // indices of the conservation equations
    enum { contiWEqIdx = Indices::conti0EqIdx + FluidSystem::wPhaseIdx };
    enum { contiNEqIdx = Indices::conti0EqIdx + FluidSystem::nPhaseIdx };

    // indices of the primary variables
    enum { pwIdx = Indices::pwIdx };
    enum { SnIdx = Indices::SnIdx };

public:
    TutorialProblemCoupled(TimeManager &timeManager)
        : ParentType(timeManager, GET_PROP_TYPE(TypeTag, GridCreator)::grid().leafView())
        , eps_(3e-6)
    {
        //set main diagonal entries of the permeability tensor to a value
        //setting to one value means: isotropic, homogeneous
        K_ = 0;
        for (int i = 0; i < dim; i++)
            K_[i][i] = 1e-7;

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
        return
            this->timeManager().timeStepIndex() > 0 &&
            (this->timeManager().timeStepIndex() % 1 == 0);
    }

    //! Returns the temperature within a finite volume. We use constant
    //! 10 degrees Celsius.
    template <class Context>
    Scalar temperature(const Context &context, int spaceIdx, int timeIdx) const
    { return 283.15; }

    /*! Intrinsic permeability tensor K \f$[m^2]\f$ depending
     *  on the position in the domain
     *
     *  \param context The execution context
     *  \param scvIdx The local index of the degree of freedom
     *
     *  Alternatively, the function intrinsicPermeabilityAtPos(const GlobalPosition& globalPos) could be defined, where globalPos
     *  is the vector including the global coordinates of the finite volume.
     */
    template <class Context>
    const Tensor &intrinsicPermeability(const Context &context, /*@\label{tutorial-coupled:permeability}@*/
                                        int spaceIdx, int timeIdx) const
    { return K_; }

    /*! Define the porosity \f$[-]\f$ of the porous medium depending
     *  on the position in the domain
     *
     *  \param context The execution context
     *  \param scvIdx The local index of the degree of freedom
     *
     *  Alternatively, the function porosityAtPos(const GlobalPosition& globalPos) could be defined, where globalPos
     *  is the vector including the global coordinates of the finite volume.
     */
    template <class Context>
    Scalar porosity(const Context &context,                    /*@\label{tutorial-coupled:porosity}@*/
                    int spaceIdx, int timeIdx) const
    { return 0.2; }

    /*! Return the parameter object for the material law (i.e. Brooks-Corey)
     *  depending on the position in the domain
     *
     *  \param context The execution context
     *  \param scvIdx The local index of the degree of freedom
     *
     *  Alternatively, the function materialLawParamsAtPos(const GlobalPosition& globalPos) could be defined, where globalPos
     *  is the vector including the global coordinates of the finite volume.
     */
    template <class Context>
    const MaterialLawParams& materialLawParams(const Context &context,            /*@\label{tutorial-coupled:matLawParams}@*/
                                               int spaceIdx, int timeIdx) const
    { return materialParams_; }

    template <class Context>
    void boundaryTypes(BoundaryTypes &bcTypes, const Context &context, int spaceIdx, int timeIdx) const
    {
        const GlobalPosition &pos = context.pos(spaceIdx, timeIdx);
        if (pos[0] < eps_) // Dirichlet conditions on left boundary
           bcTypes.setAllDirichlet();
        else // neuman for the remaining boundaries
           bcTypes.setAllNeumann();

    }

    //! Evaluates the Dirichlet boundary conditions for a finite volume
    //! on the grid boundary. Here, the 'values' parameter stores
    //! primary variables.
    template <class Context>
    void dirichlet(PrimaryVariables &values, const Context &context, int spaceIdx, int timeIdx) const
    {
        values[Indices::pwIdx] = 200.0e3; // 200 kPa = 2 bar
        values[Indices::SnIdx] = 0.0; // 0 % oil saturation on left boundary
    }

    //! Evaluates the boundary conditions for a Neumann boundary
    //! segment. Here, the 'values' parameter stores the mass flux in
    //! [kg/(m^2 * s)] in normal direction of each phase. Negative
    template <class Context>
    void neumann(RateVector &values, const Context &context, int spaceIdx, int timeIdx) const
    {
        const GlobalPosition &pos = context.pos(spaceIdx, timeIdx);
        Scalar right = this->bboxMax()[0];
        // extraction of oil on the right boundary for approx. 1.e6 seconds
        if (pos[0] > right - eps_) {
            // oil outflux of 30 g/(m * s) on the right boundary.
            values[contiWEqIdx] = 0;
            values[contiNEqIdx] = 3e-2;
        } else {
            // no-flow on the remaining Neumann-boundaries.
            values[contiWEqIdx] = 0;
            values[contiNEqIdx] = 0;
        }
    }

    //! Evaluates the source term for all phases within a given
    //! sub-control-volume. In this case, the 'values' parameter
    //! stores the rate mass generated or annihilated per volume unit
    //! in [kg / (m^3 * s)]. Positive values mean that mass is created.
    template <class Context>
    void source(RateVector &values, const Context &context, int spaceIdx, int timeIdx) const
    {
        values[contiWEqIdx] = 0.0;
        values[contiNEqIdx]= 0.0;
    }

    // Evaluates the initial value for a control volume. For this
    // method, the 'values' parameter stores primary variables.
    template <class Context>
    void initial(PrimaryVariables &values, const Context &context, int spaceIdx, int timeIdx) const
    {
        values[pwIdx] = 200.0e3; // 200 kPa = 2 bar
        values[SnIdx] = 1.0;
    }

private:
    Tensor K_;
    // Object that holds the values/parameters of the selected material law.
    MaterialLawParams materialParams_;                 /*@\label{tutorial-coupled:matParamsObject}@*/

    // small epsilon value
    Scalar eps_;
};
}

#endif
