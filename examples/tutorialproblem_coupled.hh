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
#include <dumux/material/fluids/water.hh>
#include <dumux/material/fluids/lowviscosityoil.hh>

// the numerical model
#include <dumux/boxmodels/2p/2pboxmodel.hh>

// the grid used
#include <dune/grid/yaspgrid.hh>
#include <dune/grid/io/file/dgfparser/dgfs.hh>

// the soil to be used
#include "tutorialsoil_coupled.hh"

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
        typedef typename SGrid<2,2>::ctype ctype;
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

// Set the wetting and non-wetting phases
SET_TYPE_PROP(TutorialProblemCoupled, WettingPhase, Dumux::Water); /*@\label{tutorial-coupled:set-wetting}@*/
SET_TYPE_PROP(TutorialProblemCoupled, NonwettingPhase, Dumux::LowViscosityOil);/*@\label{tutorial-coupled:set-nonwetting}@*/

// Set the soil properties
SET_PROP(TutorialProblemCoupled, Soil) /*@\label{tutorial-coupled:set-soil}@*/
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Grid)) Grid;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    
public:
    typedef Dumux::TutorialSoil<Grid, Scalar> type;
};

// Disable gravity
SET_BOOL_PROP(TutorialProblemCoupled, EnableGravity, false); /*@\label{tutorial-coupled:gravity}@*/
}

// Definition of the actual problem
template <class TypeTag = TTAG(TutorialProblemCoupled) >
class TutorialProblemCoupled : public TwoPBoxProblem<TypeTag, /*@\label{tutorial-coupled:def-problem}@*/
                                                     TutorialProblemCoupled<TypeTag> >
{
    typedef TutorialProblemCoupled<TypeTag>   ThisType;
    typedef TwoPBoxProblem<TypeTag, ThisType> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView))    GridView;

    // Grid and world dimension
    enum {       
        dim         = GridView::dimension,
        dimWorld    = GridView::dimensionworld,
    };

    typedef typename GridView::Grid::ctype                    CoordScalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar))     Scalar;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(TwoPIndices)) Indices;
    typedef typename GridView::template Codim<0>::Entity    Element;
    typedef typename GridView::template Codim<dim>::Entity  Vertex;
    typedef typename GridView::Intersection                 Intersection;
    typedef Dune::FieldVector<CoordScalar, dim>             LocalPosition;
    typedef Dune::FieldVector<CoordScalar, dimWorld>        GlobalPosition;
  
    typedef typename GET_PROP(TypeTag, PTAG(SolutionTypes))          SolutionTypes;
    typedef typename SolutionTypes::PrimaryVarVector                 PrimaryVarVector;
    typedef typename SolutionTypes::BoundaryTypeVector               BoundaryTypeVector;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FVElementGeometry)) FVElementGeometry;


public:
    TutorialProblemCoupled(const GridView &gridView)
        : ParentType(gridView)
    {}

    // Return the temperature within the domain. We use 10 degrees Celsius.
    Scalar temperature(const Element           &element,
                       const FVElementGeometry &fvElemGeom,
                       int                      scvIdx) const
    { return 283.15; };

    // Specifies which kind of boundary condition should be used for
    // which equation on a given boundary segment.
    void boundaryTypes(BoundaryTypeVector         &values,
                       const Element              &element,
                       const FVElementGeometry    &fvElemGeom,
                       const Intersection         &isIt,
                       int                         scvIdx,
                       int                         boundaryFaceIdx) const
    {    
        const GlobalPosition &pos = element.geometry().corner(scvIdx);
        if (pos[0] < eps_) // dirichlet conditions on left boundary
           values = BoundaryConditions::dirichlet;
        else // neuman for the remaining boundaries
            values = BoundaryConditions::neumann;

    }

    // Evaluate the boundary conditions for a dirichlet boundary
    // segment.  For this method, the 'values' parameter stores
    // primary variables.
    void dirichlet(PrimaryVarVector           &values,
                   const Element              &element,
                   const FVElementGeometry    &fvElemGeom,
                   const Intersection         &isIt,
                   int                         scvIdx,
                   int                         boundaryFaceIdx) const
    {
        values[Indices::pW] = 200.0e3; // 200 kPa = 2 bar
        values[Indices::sN] = 0.0; // 0 % oil saturation on left boundary
    }

    // Evaluate the boundary conditions for a neumann boundary
    // segment. For this method, the 'values' parameter stores the
    // mass flux in normal direction of each phase. Negative values
    // mean influx.
    void neumann(PrimaryVarVector           &values,
                 const Element              &element,
                 const FVElementGeometry    &fvElemGeom,
                 const Intersection         &isIt,
                 int                         scvIdx,
                 int                         boundaryFaceIdx) const
    {
        const GlobalPosition &pos = 
            fvElemGeom.boundaryFace[boundaryFaceIdx].ipGlobal;
        Scalar right = this->bboxMax()[0];
        if (pos[0] > right - eps_) {
            // oil outflux of 0.3 g/(m * s) on the right boundary of
            // the domain.
            values[Indices::phase2Mass(Indices::wPhase)] = 0;
            values[Indices::phase2Mass(Indices::nPhase)] = 0.3e-3;
        } else {
            // no-flow on the remaining neumann-boundaries
            values[Indices::phase2Mass(Indices::wPhase)] = 0;
            values[Indices::phase2Mass(Indices::nPhase)] = 0;
        }
    }

    // Evaluate the initial value for a control volume. For this
    // method, the 'values' parameter stores primary variables.
    void initial(PrimaryVarVector        &values,
                 const Element           &element,
                 const FVElementGeometry &fvElemGeom,
                 int                      scvIdx) const
    {
        values[Indices::pW] = 200.0e3; // 200 kPa = 2 bar
        values[Indices::sN] = 1.0;
    }

    // Evaluate the source term for all phases within a given
    // sub-control-volume. For this method, the \a values parameter
    // stores the rate mass generated or annihilate per volume
    // unit. Positive values mean that mass is created, negative ones
    // mean that it vanishes.
    void source(PrimaryVarVector        &values,
                const Element           &element,
                const FVElementGeometry &fvElemGeom,
                int                      scvIdx) const
    {
        values[Indices::phase2Mass(Indices::wPhase)] = 0.0;
        values[Indices::phase2Mass(Indices::nPhase)] = 0.0;
    }

private:
    static const Scalar eps_ = 3e-6;
};
}

#endif
