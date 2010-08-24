// $Id$
/*****************************************************************************
 *   Copyright (C) 2008-2009 by Melanie Darcis                               *
 *   Copyright (C) 2009 by Andreas Lauser                                    *
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
#ifndef DUMUX_TUTORIALPROBLEM_COUPLED_HH
#define DUMUX_TUTORIALPROBLEM_COUPLED_HH

// fluid properties
#include <dumux/material/fluidsystems/h2o_n2_system.hh>

// the numerical model
#include <dumux/boxmodels/2p/2pmodel.hh>

// the grid used
#include <dune/grid/yaspgrid.hh>
#include <dune/grid/io/file/dgfparser/dgfs.hh>

// assign parameters dependent on space (e.g. spatial parameters)
#include "tutorialspatialparameters_coupled.hh"

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
{
    typedef Dumux::TutorialProblemCoupled<TTAG(TutorialProblemCoupled)> type;
};

// Set the grid
SET_PROP(TutorialProblemCoupled, Grid) /*@\label{tutorial-coupled:set-grid}@*/
{
    typedef Dune::SGrid<2,2> type;
    static type *create() /*@\label{tutorial-coupled:create-grid-method}@*/
    {
        typedef typename type::ctype ctype;
        Dune::FieldVector<int, 2> cellRes;
        Dune::FieldVector<ctype, 2> lowerLeft(0.0);
        Dune::FieldVector<ctype, 2> upperRight;
        cellRes[0] = 30;
        cellRes[1] = 10;
        upperRight[0] = 300;
        upperRight[1] = 60;
        return new Dune::SGrid<2,2>(cellRes,
                                    lowerLeft,
                                    upperRight);
    }
};

// Select fluid system
SET_PROP(TutorialProblemCoupled,   FluidSystem) /*@\label{tutorial-coupled:set-fluidsystem}@*/
{
    typedef Dumux::H2O_N2_System<TypeTag> type;
};

// Set the spatial parameters
SET_PROP(TutorialProblemCoupled, SpatialParameters) /*@\label{tutorial-coupled:set-spatialparameters}@*/
{
    typedef Dumux::TutorialSpatialParametersCoupled<TypeTag> type;
};

// Disable gravity
SET_BOOL_PROP(TutorialProblemCoupled, EnableGravity, false); /*@\label{tutorial-coupled:gravity}@*/
}

// Definition of the actual problem
template <class TypeTag = TTAG(TutorialProblemCoupled) >
class TutorialProblemCoupled : public TwoPProblem<TypeTag> /*@\label{tutorial-coupled:def-problem}@*/
{
    typedef TutorialProblemCoupled<TypeTag> ThisType;
    typedef TwoPProblem<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(TimeManager)) TimeManager;

    // Grid and world dimension
    enum {
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld,
    };

    typedef typename GridView::Grid::ctype CoordScalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(TwoPIndices)) Indices;
    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::template Codim<dim>::Entity Vertex;
    typedef typename GridView::Intersection Intersection;
    typedef Dune::FieldVector<CoordScalar, dim> LocalPosition;
    typedef Dune::FieldVector<CoordScalar, dimWorld> GlobalPosition;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(PrimaryVariables)) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(BoundaryTypes)) BoundaryTypes;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FVElementGeometry)) FVElementGeometry;

public:
    TutorialProblemCoupled(TimeManager &timeManager,
                           const GridView &gridView)
        : ParentType(timeManager, gridView)
    {}

    // Return the temperature within the domain. We use 10 degrees Celsius.
    Scalar temperature(const Element &element,
                       const FVElementGeometry &fvElemGeom,
                       int scvIdx) const
    { return 283.15; };

    // Specifies which kind of boundary condition should be used for
    // which equation on a given boundary segment.
    void boundaryTypes(BoundaryTypes &BCtype,
                       const Element &element,
                       const FVElementGeometry &fvElemGeom,
                       const Intersection &isIt,
                       int scvIdx,
                       int boundaryFaceIdx) const
    {
        const GlobalPosition &pos = element.geometry().corner(scvIdx);
        if (pos[0] < eps_) // dirichlet conditions on left boundary
           BCtype.setAllDirichlet();
        else // neuman for the remaining boundaries
           BCtype.setAllNeumann();

    }

    // Evaluate the boundary conditions for a dirichlet boundary
    // segment.  For this method, the 'values' parameter stores
    // primary variables.
    void dirichlet(PrimaryVariables &values,
                   const Element &element,
                   const FVElementGeometry &fvElemGeom,
                   const Intersection &isIt,
                   int scvIdx,
                   int boundaryFaceIdx) const
    {
        values[Indices::pwIdx] = 200.0e3; // 200 kPa = 2 bar
        values[Indices::SnIdx] = 0.0; // 0 % oil saturation on left boundary
    }

    // Evaluate the boundary conditions for a neumann boundary
    // segment. For this method, the 'values' parameter stores the
    // mass flux in normal direction of each phase. Negative values
    // mean influx.
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
        if (pos[0] > right - eps_) {
            // oil outflux of 0.3 g/(m * s) on the right boundary of
            // the domain.
            values[Indices::contiWEqIdx] = 0;
            values[Indices::contiNEqIdx] = 0.3e-3;
        } else {
            // no-flow on the remaining neumann-boundaries
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
    // sub-control-volume. For this method, the \a values parameter
    // stores the rate mass generated or annihilate per volume
    // unit. Positive values mean that mass is created, negative ones
    // mean that it vanishes.
    void source(PrimaryVariables &values,
                const Element &element,
                const FVElementGeometry &fvElemGeom,
                int scvIdx) const
    {
        values[Indices::contiWEqIdx] = 0.0;
        values[Indices::contiNEqIdx]= 0.0;
    }

private:
    static const Scalar eps_ = 3e-6;
};
}

#endif
