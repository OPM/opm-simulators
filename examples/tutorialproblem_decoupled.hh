// $Id$
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
#ifndef TUTORIALPROBLEM_DECOUPLED_HH
#define TUTORIALPROBLEM_DECOUPLED_HH

#include "dumux/fractionalflow/fractionalflowproblem.hh"

namespace Dumux
{

/** \todo Please doc me! */

template<class GridView, class Scalar, class VariableClass> class TutorialProblemDecoupled /*@\label{tutorial-decoupled:tutorialproblem}@*/
    : public FractionalFlowProblem<GridView, Scalar, VariableClass>
{
    enum
        {dim=GridView::dimension, dimWorld = GridView::dimensionworld};
    enum{wetting = 0, nonwetting = 1};
    typedef typename GridView::Grid Grid;
    typedef typename GridView::Intersection Intersection;
    typedef typename GridView::Traits::template Codim<0>::Entity Element;
    typedef Dune::FieldVector<Scalar,dim> LocalPosition;
    typedef Dune::FieldVector<Scalar,dimWorld> GlobalPosition;

public:
    TutorialProblemDecoupled(VariableClass& variables, Fluid& wettingphase, Fluid& nonwettingphase, Matrix2p<Grid, Scalar>& soil,
                             TwoPhaseRelations<Grid, Scalar>& materialLaw = *(new TwoPhaseRelations<Grid,Scalar>),
                             const Dune::FieldVector<Scalar,dim> Left = 0, const Dune::FieldVector<Scalar,dim> Right = 0)
        : FractionalFlowProblem<GridView, Scalar, VariableClass>(variables, wettingphase, nonwettingphase, soil, materialLaw),
          Left_(Left[0]), Right_(Right[0]), eps_(1e-8)
    {}

    // function returning source/sink terms for the pressure equation
    // depending on the position within the domain
    virtual std::vector<Scalar> source(const GlobalPosition& globalPos,
                          const Element& e, /*@\label{tutorial-decoupled:qpress}@*/
                          const LocalPosition& localPos)
    {
        return std::vector<Scalar>(2,0.0);
    }

    using FractionalFlowProblem<GridView, Scalar, VariableClass>::bctypePress;

    // function returning the boundary condition type for solution
    // of the pressure equation depending on the position within the domain
    typename BoundaryConditions::Flags bctypePress(const GlobalPosition& globalPos, const Intersection& intersection) const /*@\label{tutorial-decoupled:bctypepress}@*/
    {
        if (globalPos[0] < eps_)
        {
            return BoundaryConditions::dirichlet;
        }
        // all other boundaries
        return BoundaryConditions::neumann;
    }

    // function returning the boundary condition type for solution
    // of the saturation equation depending on the position within the domain
    BoundaryConditions::Flags bctypeSat (const GlobalPosition& globalPos, const Intersection& intersection) const /*@\label{tutorial-decoupled:bctypesat}@*/
    {
        if (globalPos[0]> (Right_ - eps_) || globalPos[0] < eps_)
        {
            return Dumux::BoundaryConditions::dirichlet;
        }
        // all other boundaries
        return Dumux::BoundaryConditions::neumann;
    }

    // function returning the Dirichlet boundary condition for the solution
    // of the pressure equation depending on the position within the domain
    Scalar dirichletPress(const GlobalPosition& globalPos, const Intersection& intersection) const /*@\label{tutorial-decoupled:gpress}@*/
    {
        return 1e6;
    }

    // function returning the Dirichlet boundary condition for the solution
    // of the saturation equation depending on the position within the domain
    Scalar dirichletSat(const GlobalPosition& globalPos, const Intersection& intersection) const /*@\label{tutorial-decoupled:gsat}@*/
    {
        if (globalPos[0] < eps_)
        {
            return 1.0;
        }
        // all other boundaries
        return 0.0;
    }

    using FractionalFlowProblem<GridView, Scalar, VariableClass>::neumannPress;

    // function returning the Neumann boundary condition for the solution
    // of the pressure equation depending on the position within the domain
    std::vector<Scalar> neumannPress(const GlobalPosition& globalPos, const Intersection& intersection) const /*@\label{tutorial-decoupled:jpress}@*/
    {
        std::vector<Scalar> neumannFlux(2, 0.0);
        if (globalPos[0]> Right_ - eps_)
        {
            neumannFlux[nonwetting] = 3e-4;
        }
        // all other boundaries
        return neumannFlux;
    }

    // function returning the initial saturation
    // depending on the position within the domain
    Scalar initSat (const GlobalPosition& globalPos, const Element& e, /*@\label{tutorial-decoupled:initsat}@*/
                    const Dune::FieldVector<Scalar,dim>& xi) const
    {
        return 0.0;
    }

private:
    Scalar Left_;
    Scalar Right_;

    Scalar eps_;
};
} // end namespace
#endif
