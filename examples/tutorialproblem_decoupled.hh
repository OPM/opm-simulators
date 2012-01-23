// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
*   Copyright (C) 2007-2008 by Klaus Mosthaf                                *
*   Copyright (C) 2007-2008 by Bernd Flemisch                               *
*   Copyright (C) 2008-2009 by Andreas Lauser                               *
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
 * \brief problem for the sequential tutorial
 */
#ifndef DUMUX_TUTORIALPROBLEM_DECOUPLED_HH // guardian macro /*@\label{tutorial-decoupled:guardian1}@*/
#define DUMUX_TUTORIALPROBLEM_DECOUPLED_HH // guardian macro /*@\label{tutorial-decoupled:guardian2}@*/

// the grid includes
#include <dune/grid/sgrid.hh>

// dumux 2p-decoupled environment
#include <dumux/decoupled/2p/diffusion/fv/fvpressureproperties2p.hh>
#include <dumux/decoupled/2p/transport/fv/fvtransportproperties2p.hh>
#include <dumux/decoupled/2p/impes/impesproblem2p.hh>
#include <dumux/decoupled/2p/impes/impesproblem2p.hh> /*@\label{tutorial-decoupled:parent-problem}@*/

// assign parameters dependent on space (e.g. spatial parameters)
#include "tutorialspatialparameters_decoupled.hh" /*@\label{tutorial-decoupled:spatialparameters}@*/

// the components that are used
#include <dumux/material/components/h2o.hh>
#include <dumux/material/components/oil.hh>

namespace Dumux
{

template<class TypeTag>
class TutorialProblemDecoupled;

//////////
// Specify the properties for the lens problem
//////////
namespace Properties
{
// create a new type tag for the problem
NEW_TYPE_TAG(TutorialProblemDecoupled, INHERITS_FROM(FVPressureTwoP, FVTransportTwoP, IMPESTwoP, TutorialSpatialParametersDecoupled)); /*@\label{tutorial-decoupled:create-type-tag}@*/

// Set the problem property
SET_PROP(TutorialProblemDecoupled, Problem) /*@\label{tutorial-decoupled:set-problem}@*/
{
    typedef Dumux::TutorialProblemDecoupled<TypeTag> type;
};

// Set the grid type
SET_PROP(TutorialProblemDecoupled, Grid) /*@\label{tutorial-decoupled:grid-begin}@*/
{
    typedef Dune::SGrid<2, 2> type; /*@\label{tutorial-decoupled:set-grid-type}@*/
    static type *create() /*@\label{tutorial-decoupled:create-grid-method}@*/
    {
        typedef typename type::ctype ctype;
        Dune::FieldVector<int, 2> cellRes;  // vector holding resolution of the grid
        Dune::FieldVector<ctype, 2> lowerLeft(0.0); // Coordinate of lower left corner of the grid
        Dune::FieldVector<ctype, 2> upperRight; // Coordinate of upper right corner of the grid
        cellRes[0] = 100;
        cellRes[1] = 1;
        upperRight[0] = 300;
        upperRight[1] = 60;
        return new Dune::SGrid<2,2>(cellRes,
                                    lowerLeft,
                                    upperRight);
    } /*@\label{tutorial-decoupled:grid-end}@*/
};

// Set the wetting phase
SET_PROP(TutorialProblemDecoupled, WettingPhase) /*@\label{tutorial-decoupled:2p-system-start}@*/
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
public:
    typedef Dumux::LiquidPhase<Scalar, Dumux::H2O<Scalar> > type; /*@\label{tutorial-decoupled:wettingPhase}@*/
};

// Set the non-wetting phase
SET_PROP(TutorialProblemDecoupled, NonwettingPhase)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
public:
    typedef Dumux::LiquidPhase<Scalar, Dumux::Oil<Scalar> > type; /*@\label{tutorial-decoupled:nonwettingPhase}@*/
}; /*@\label{tutorial-decoupled:2p-system-end}@*/

SET_SCALAR_PROP(TutorialProblemDecoupled, CFLFactor, 0.5); /*@\label{tutorial-decoupled:cfl}@*/

// Disable gravity
SET_BOOL_PROP(TutorialProblemDecoupled, EnableGravity, false); /*@\label{tutorial-decoupled:gravity}@*/
} /*@\label{tutorial-decoupled:propertysystem-end}@*/

/*! \ingroup DecoupledProblems
 * @brief Problem class for the decoupled tutorial
*/
template<class TypeTag>
class TutorialProblemDecoupled: public IMPESProblem2P<TypeTag> /*@\label{tutorial-decoupled:def-problem}@*/
{
    typedef IMPESProblem2P<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, TimeManager) TimeManager;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, FluidState) FluidState;

    typedef typename GET_PROP_TYPE(TypeTag, BoundaryTypes) BoundaryTypes;
    typedef typename GET_PROP(TypeTag, SolutionTypes) SolutionTypes;
    typedef typename SolutionTypes::PrimaryVariables PrimaryVariables;

    enum
    {
        dim = GridView::dimension, dimWorld = GridView::dimensionworld
    };

    enum
    {
        wPhaseIdx = Indices::wPhaseIdx,
        nPhaseIdx = Indices::nPhaseIdx,
        pWIdx = Indices::pwIdx,
        SwIdx = Indices::SwIdx,
        pressEqIdx = Indices::pressEqIdx,
        satEqIdx = Indices::satEqIdx
    };

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;

    typedef typename GridView::Traits::template Codim<0>::Entity Element;
    typedef typename GridView::Intersection Intersection;
    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;

public:
    TutorialProblemDecoupled(TimeManager &timeManager, const GridView &gridView)
        : ParentType(timeManager, gridView) /*@\label{tutorial-decoupled:constructor-problem}@*/
    {    }

    //! The problem name.
    /*! This is used as a prefix for files generated by the simulation.
    */
    const char *name() const    /*@\label{tutorial-decoupled:name}@*/
    {
        return "tutorial_decoupled";
    }

    //!  Returns true if a restart file should be written.
    /* The default behaviour is to write no restart file.
     */
    bool shouldWriteRestartFile() const /*@\label{tutorial-decoupled:restart}@*/
    {
        return false;
    }

    //! Returns true if the current solution should be written to disk (i.e. as a VTK file)
    /*! The default behaviour is to write out every the solution for
     *  very time step. Else, change divisor.
     */
    bool shouldWriteOutput() const /*@\label{tutorial-decoupled:output}@*/
    {
        return this->timeManager().timeStepIndex() > 0 &&
        (this->timeManager().timeStepIndex() % 1 == 0);
    }

    //! Returns the temperature within the domain at position globalPos.
    /*! This problem assumes a temperature of 10 degrees Celsius.
     *
     *  \param element The finite volume element
     *
     * Alternatively, the function temperatureAtPos(const GlobalPosition& globalPos) could be defined, where globalPos
     * is the vector including the global coordinates of the finite volume.
     */
    Scalar temperature(const Element& element) const /*@\label{tutorial-decoupled:temperature}@*/
    {
        return 273.15 + 10; // -> 10°C
    }

    //! Returns a constant pressure to enter material laws at position globalPos.
    /* For incrompressible simulations, a constant pressure is necessary
     * to enter the material laws to gain a constant density etc. In the compressible
     * case, the pressure is used for the initialization of material laws.
     *
     * \param element The finite volume element
     *
     * Alternatively, the function referencePressureAtPos(const GlobalPosition& globalPos) could be defined, where globalPos
     * is the vector including the global coordinates of the finite volume.
     */
    Scalar referencePressure(const Element& element) const /*@\label{tutorial-decoupled:refPressure}@*/
    {
        return 2e5;
    }

    //! Source of mass \f$ [\frac{kg}{m^3 \cdot s}] \f$ of a finite volume.
    /*! Evaluate the source term for all phases within a given
     *  volume.
     *
     *  \param values Includes sources for the two phases
     *  \param element The finite volume element
     *
     *  The method returns the mass generated (positive) or
     *  annihilated (negative) per volume unit.
     *
     * Alternatively, the function sourceAtPos(PrimaryVariables &values, const GlobalPosition& globalPos) could be defined, where globalPos
     * is the vector including the global coordinates of the finite volume.
     */
    void source(PrimaryVariables &values, const Element& element) const /*@\label{tutorial-decoupled:source}@*/
    {
        values = 0;
    }

    //! Type of boundary conditions at position globalPos.
    /*! Defines the type the boundary condition for the pressure equation,
     *  either pressure (dirichlet) or flux (neumann),
     *  and for the transport equation,
     *  either saturation (dirichlet) or flux (neumann).
     *
     *  \param bcTypes Includes the types of boundary conditions
     *  \param globalPos The position of the center of the finite volume
     *
     *  Alternatively, the function boundaryTypes(PrimaryVariables &values, const Intersection& intersection) could be defined,
     *  where intersection is the boundary intersection.
     */
    void boundaryTypesAtPos(BoundaryTypes &bcTypes, const GlobalPosition& globalPos) const /*@\label{tutorial-decoupled:bctype}@*/
    {
            if (globalPos[0] < this->bboxMin()[0] + eps_)
            {
                bcTypes.setDirichlet(pressEqIdx);
                bcTypes.setDirichlet(satEqIdx);
//                bcTypes.setAllDirichlet(); // alternative if the same BC is used for both types of equations
            }
            // all other boundaries
            else
            {
                bcTypes.setNeumann(pressEqIdx);
                bcTypes.setNeumann(satEqIdx);
//                bcTypes.setAllNeumann(); // alternative if the same BC is used for both types of equations
            }
    }
    //! Value for dirichlet boundary condition at position globalPos.
    /*! In case of a dirichlet BC for the pressure equation the pressure \f$ [Pa] \f$, and for the transport equation the saturation [-]
     *  have to be defined on boundaries.
     *
     *  \param values Values of primary variables at the boundary
     *  \param intersection The boundary intersection
     *
     *  Alternatively, the function dirichletAtPos(PrimaryVariables &values, const GlobalPosition& globalPos) could be defined, where globalPos
     *  is the vector including the global coordinates of the finite volume.
     */
    void dirichlet(PrimaryVariables &values, const Intersection& intersection) const /*@\label{tutorial-decoupled:dirichlet}@*/
    {
        values[pWIdx] = 2e5;
        values[SwIdx] = 1.0;
    }
    //! Value for neumann boundary condition \f$ [\frac{kg}{m^3 \cdot s}] \f$ at position globalPos.
    /*! In case of a neumann boundary condition, the flux of matter
     *  is returned as a vector.
     *
     *  \param values Boundary flux values for the different phases
     *  \param globalPos The position of the center of the finite volume
     *
     *  Alternatively, the function neumann(PrimaryVariables &values, const Intersection& intersection) could be defined,
     *  where intersection is the boundary intersection.
     */
    void neumannAtPos(PrimaryVariables &values, const GlobalPosition& globalPos) const /*@\label{tutorial-decoupled:neumann}@*/
    {
        values = 0;
        if (globalPos[0] > this->bboxMax()[0] - eps_)
        {
            values[nPhaseIdx] = 3e-2;
        }
    }
    //! Initial condition at position globalPos.
    /*! Only initial values for saturation have to be given!
     *
     *  \param values Values of primary variables
     *  \param element The finite volume element
     *
     *  Alternatively, the function initialAtPos(PrimaryVariables &values, const GlobalPosition& globalPos) could be defined, where globalPos
     *  is the vector including the global coordinates of the finite volume.
     */
    void initial(PrimaryVariables &values,
            const Element &element) const /*@\label{tutorial-decoupled:initial}@*/
    {
        values = 0;
    }

private:
    static constexpr Scalar eps_ = 1e-6;
};
} //end namespace

#endif
