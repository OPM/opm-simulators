// $Id$
/*****************************************************************************
 *   Copyright (C) 2008-2009 by Melanie Darcis, Klaus Mosthaf                *
 *   Copyright (C) 2009 by Andreas Lauser                                    *
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
 * \brief Tutorial problem for a fully coupled twophase box model.
 */
#ifndef DUMUX_TUTORIALPROBLEM_COUPLED_HH    // guardian macro /*@\label{tutorial-coupled:guardian1}@*/
#define DUMUX_TUTORIALPROBLEM_COUPLED_HH    // guardian macro /*@\label{tutorial-coupled:guardian2}@*/

// the numerical model
#include <dumux/boxmodels/2p/2pmodel.hh>

// the DUNE grid used
#include <dune/grid/sgrid.hh>

// spatialy dependent parameters
#include "tutorialspatialparameters_coupled.hh"

// the components that are used
#include <dumux/material/components/h2o.hh>
#include <dumux/material/components/oil.hh>

namespace Dumux
{

// forward declaration of the problem class
template <class TypeTag>
class TutorialProblemCoupled;

namespace Properties
{
// create a new type tag for the problem
NEW_TYPE_TAG(TutorialProblemCoupled, INHERITS_FROM(BoxTwoP)); /*@\label{tutorial-coupled:create-type-tag}@*/

// Set the "Problem" property
SET_PROP(TutorialProblemCoupled, Problem) /*@\label{tutorial-coupled:set-problem}@*/
{ typedef Dumux::TutorialProblemCoupled<TypeTag> type;};

// Set the grid
SET_PROP(TutorialProblemCoupled, Grid) /*@\label{tutorial-coupled:set-grid}@*/
{
    typedef Dune::SGrid<2,2> type;
    static type *create() /*@\label{tutorial-coupled:create-grid-method}@*/
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
    }
};

// Set the wetting phase
SET_PROP(TutorialProblemCoupled, WettingPhase) /*@\label{tutorial-coupled:2p-system-start}@*/
{
private: typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
public: typedef Dumux::LiquidPhase<Scalar, Dumux::H2O<Scalar> > type; /*@\label{tutorial-coupled:wettingPhase}@*/
};

// Set the non-wetting phase
SET_PROP(TutorialProblemCoupled, NonwettingPhase)
{
private: typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
public: typedef Dumux::LiquidPhase<Scalar, Dumux::Oil<Scalar> > type; /*@\label{tutorial-coupled:nonwettingPhase}@*/
}; /*@\label{tutorial-coupled:2p-system-end}@*/


// Set the spatial parameters
SET_PROP(TutorialProblemCoupled, SpatialParameters) /*@\label{tutorial-coupled:set-spatialparameters}@*/
{ typedef Dumux::TutorialSpatialParametersCoupled<TypeTag> type; };

// Disable gravity
SET_BOOL_PROP(TutorialProblemCoupled, EnableGravity, false); /*@\label{tutorial-coupled:gravity}@*/
}

/*!
* \ingroup TwoPBoxModel
*
* \brief Tutorial problem for a fully coupled twophase box model.
*/

// Definition of the actual problem
template <class TypeTag>
class TutorialProblemCoupled : public TwoPProblem<TypeTag> /*@\label{tutorial-coupled:def-problem}@*/
{
    typedef TwoPProblem<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;

    // Grid dimension
    enum { dim = GridView::dimension };

    // Types from DUNE-Grid
    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::template Codim<dim>::Entity Vertex;
    typedef typename GridView::Intersection Intersection;
    typedef Dune::FieldVector<Scalar, dim> GlobalPosition;

    // Dumux specific types
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(TimeManager)) TimeManager;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(TwoPIndices)) Indices;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(PrimaryVariables)) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(BoundaryTypes)) BoundaryTypes;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FVElementGeometry)) FVElementGeometry;

public:
    TutorialProblemCoupled(TimeManager &timeManager,
                           const GridView &gridView)
        : ParentType(timeManager, gridView)
    {
    }

    // Specifies the problem name. This is used as a prefix for files
    // generated by the simulation.
    const char *name() const
    { return "tutorial_coupled"; }

    //!  Returns true if a restart file should be written.
    /* The default behaviour is to write no restart file.
     */
    bool shouldWriteRestartFile() const /*@\label{tutorial-coupled:restart}@*/
    {
        return false;
    }

    //! Returns true if the current solution should be written to disk (i.e. as a VTK file)
    /*! The default behaviour is to write out the solution for
     *  every time step. Else, the user has to change the divisor in this function.
     */
    bool shouldWriteOutput() const /*@\label{tutorial-coupled:output}@*/
    {
        return this->timeManager().timeStepIndex() > 0 &&
        (this->timeManager().timeStepIndex() % 1 == 0);
    }

    // Return the temperature within a finite volume. We use constant
    // 10 degrees Celsius.
    Scalar temperature(const Element &element,
                       const FVElementGeometry &fvElemGeom,
                       int scvIdx) const
    { return 283.15; };

    // Specifies which kind of boundary condition should be used for
    // which equation for a finite volume on the boundary.
    void boundaryTypes(BoundaryTypes &BCtypes, const Vertex &vertex) const
    {
        const GlobalPosition &pos = vertex.geometry().center();
        if (pos[0] < eps_) // Dirichlet conditions on left boundary
           BCtypes.setAllDirichlet();
        else // neuman for the remaining boundaries
           BCtypes.setAllNeumann();

    }

    // Evaluate the Dirichlet boundary conditions for a finite volume
    // on the grid boundary. Here, the 'values' parameter stores
    // primary variables.
    void dirichlet(PrimaryVariables &values, const Vertex &vertex) const
    {
        values[Indices::pwIdx] = 200.0e3; // 200 kPa = 2 bar
        values[Indices::SnIdx] = 0.0; // 0 % oil saturation on left boundary
    }

    // Evaluate the boundary conditions for a Neumann boundary
    // segment. Here, the 'values' parameter stores the mass flux in
    // normal direction of each phase. Negative values mean influx.
    void neumann(PrimaryVariables &values,
                 const Element &element,
                 const FVElementGeometry &fvElemGeom,
                 const Intersection &isIt,
                 int scvIdx,
                 int boundaryFaceIdx) const
    {
        const GlobalPosition &pos =
            fvElemGeom.boundaryFace[boundaryFaceIdx].ipGlobal;
        Scalar right = this->bboxMax()[0];
        // extraction of oil on the right boundary for approx. 1.e6 seconds
        if (pos[0] > right - eps_) {
            // oil outflux of 30 g/(m * s) on the right boundary.
            values[Indices::contiWEqIdx] = 0;
            values[Indices::contiNEqIdx] = 3e-2;
        } else {
            // no-flow on the remaining Neumann-boundaries.
            values[Indices::contiWEqIdx] = 0;
            values[Indices::contiNEqIdx] = 0;
        }
    }

    // Evaluate the initial value for a control volume. For this
    // method, the 'values' parameter stores primary variables.
    void initial(PrimaryVariables &values,
                 const Element &element,
                 const FVElementGeometry &fvElemGeom,
                 int scvIdx) const
    {
        values[Indices::pwIdx] = 200.0e3; // 200 kPa = 2 bar
        values[Indices::SnIdx] = 1.0;
    }

    // Evaluate the source term for all phases within a given
    // sub-control-volume. In this case, the 'values' parameter stores
    // the rate mass generated or annihilate per volume unit. Positive
    // values mean that mass is created.
    void source(PrimaryVariables &values,
                const Element &element,
                const FVElementGeometry &fvElemGeom,
                int scvIdx) const
    {
        values[Indices::contiWEqIdx] = 0.0;
        values[Indices::contiNEqIdx]= 0.0;
    }

private:
    // small epsilon value
    static const Scalar eps_ = 3e-6;
};
}

#endif
