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
#ifndef DUNE_TUTORIALPROBLEM_COUPLED_HH
#define DUNE_TUTORIALPROBLEM_COUPLED_HH

// fluid properties
#include <dumux/material/fluids/water.hh>
#include <dumux/material/fluids/oil.hh>

// the numerical model
#include <dumux/boxmodels/2p/2pboxmodel.hh>

// the grid used
#include <dune/grid/yaspgrid.hh>
#include <dune/grid/io/file/dgfparser/dgfyasp.hh>

// the soil to be used
#include "tutorialsoil_coupled.hh"

namespace Dune
{

// forward declaration of the problem class
template <class TypeTag>
class TutorialProblemCoupled;

//////////
// Specify the properties of the problem
//////////
namespace Properties
{
// create a new type tag for the problem
NEW_TYPE_TAG(TutorialProblemCoupled, INHERITS_FROM(BoxTwoP));

// Set the "Problem" property
SET_PROP(TutorialProblemCoupled, Problem)
{
    typedef Dune::TutorialProblemCoupled<TTAG(TutorialProblemCoupled)> type;
};

// Set the grid type
SET_TYPE_PROP(TutorialProblemCoupled, Grid, Dune::YaspGrid<2>);

// Set the wetting phase
SET_TYPE_PROP(TutorialProblemCoupled, WettingPhase, Dune::Water);

// Set the non-wetting phase
SET_TYPE_PROP(TutorialProblemCoupled, NonwettingPhase, Dune::Oil);

// Set the soil properties
SET_PROP(TutorialProblemCoupled, Soil)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Grid)) Grid;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    
public:
    typedef Dune::TutorialSoil<Grid, Scalar> type;
};

// Enable gravity
SET_BOOL_PROP(TutorialProblemCoupled, EnableGravity, true);
}

/*!
 * \ingroup TwoPBoxProblems
 * \brief The problem used for the tutorial of the coupled models
 */
template <class TypeTag = TTAG(TutorialProblemCoupled) >
class TutorialProblemCoupled : public TwoPBoxProblem<TypeTag, 
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
    typedef typename GridView::IntersectionIterator         IntersectionIterator;
    typedef Dune::FieldVector<CoordScalar, dim>             LocalPosition;
    typedef Dune::FieldVector<CoordScalar, dimWorld>        GlobalPosition;
  
    typedef typename GET_PROP(TypeTag, PTAG(SolutionTypes))          SolutionTypes;
    typedef typename SolutionTypes::PrimaryVarVector                 PrimaryVarVector;
    typedef typename SolutionTypes::BoundaryTypeVector               BoundaryTypeVector;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FVElementGeometry)) FVElementGeometry;


public:
    TutorialProblemCoupled(const GridView &gridView)
        : ParentType(gridView)
    {
    }


    /*!
     * \name Problem parameters
     */
    // \{

    /*!
     * \brief Returns the temperature within the domain.
     * 
     * We use 10°C...
     */
    Scalar temperature() const
    { return 283.15; };

    // \}

    /*!
     * \name Boundary conditions
     */
    // \{

    /*! 
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary segment.
     */
    void boundaryTypes(BoundaryTypeVector         &values,
                       const Element              &element,
                       const FVElementGeometry    &fvElemGeom,
                       const IntersectionIterator &isIt,
                       int                         scvIdx,
                       int                         boundaryFaceIdx) const
    {    
        const GlobalPosition &pos = element.geometry().corner(scvIdx);
        if (pos[0] < eps_)
            // dirichlet conditions on left boundary
           values = BoundaryConditions::dirichlet;
        else 
            // neuman for the remaining boundaries
            values = BoundaryConditions::neumann;

    }

    /*! 
     * \brief Evaluate the boundary conditions for a dirichlet
     *        boundary segment.
     *
     * For this method, the \a values parameter stores primary variables.
     */
    void dirichlet(PrimaryVarVector           &values,
                   const Element              &element,
                   const FVElementGeometry    &fvElemGeom,
                   const IntersectionIterator &isIt,
                   int                         scvIdx,
                   int                         boundaryFaceIdx) const
    {
        values[Indices::pW] = 200.0e3; // 200 000 Pa = 2 bar
        values[Indices::sN] = 1.0; // 100 % oil saturation
    }

    /*! 
     * \brief Evaluate the boundary conditions for a neumann
     *        boundary segment.
     *
     * For this method, the \a values parameter stores the mass flux
     * in normal direction of each phase. Negative values mean influx.
     */
    void neumann(PrimaryVarVector           &values,
                 const Element              &element,
                 const FVElementGeometry    &fvElemGeom,
                 const IntersectionIterator &isIt,
                 int                         scvIdx,
                 int                         boundaryFaceIdx) const
    {
        const GlobalPosition &pos = element.geometry().corner(scvIdx);
        if (pos[0]> right_ - eps_) {
            // outflow of 0.3 g/(m * s) oil on the right boundary of the
            // domain
            values[Indices::phase2Mass(Indices::wPhase)] = 0;
            values[Indices::phase2Mass(Indices::nPhase)] = 0.3e-3;
        } else {
            // no-flow on the remaining neumann-boundaries
            values[Indices::phase2Mass(Indices::wPhase)] = 0;
            values[Indices::phase2Mass(Indices::nPhase)] = 0;
        }
    }
    // \}

    /*!
     * \name Volume terms
     */
    // \{

    /*! 
     * \brief Evaluate the initial value for a control volume.
     *
     * For this method, the \a values parameter stores primary
     * variables.
     */
    void initial(PrimaryVarVector        &values,
                 const Element           &element,
                 const FVElementGeometry &fvElemGeom,
                 int                      scvIdx) const
    {
        values[Indices::pW] = 200.0e3; // 200 000 Pa = 2 bar
        values[Indices::sN] = 1.0;
    }

    /*! 
     * \brief Evaluate the source term for all phases within a given
     *        sub-control-volume.
     *
     * For this method, the \a values parameter stores the rate mass
     * generated or annihilate per volume unit. Positive values mean
     * that mass is created, negative ones mean that it vanishes.
     */
    void source(PrimaryVarVector        &values,
                const Element           &element,
                const FVElementGeometry &,
                int subControlVolumeIdx) const
    {
        values[Indices::phase2Mass(Indices::wPhase)] = 0.0;
        values[Indices::phase2Mass(Indices::nPhase)] = 0.0;
    }
    // \}

private:
    static const Scalar eps_ = 3e-6;

    static const Scalar right_ = 5.0;
};
} //end namespace

#endif
