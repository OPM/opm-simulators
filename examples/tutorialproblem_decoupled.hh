// $Id$
/*****************************************************************************
*   Copyright (C) 2007-2008 by Klaus Mosthaf                                *
*   Copyright (C) 2007-2008 by Bernd Flemisch                               *
*   Copyright (C) 2008-2009 by Andreas Lauser                               *
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
#ifndef DUMUX_TUTORIALPROBLEM_DECOUPLED_HH
#define DUMUX_TUTORIALPROBLEM_DECOUPLED_HH

#if HAVE_UG
#include <dune/grid/uggrid.hh>
#endif

#include <dune/grid/yaspgrid.hh>
#include <dune/grid/sgrid.hh>

// fluid properties
#include <dumux/material/fluidsystems/2p_system.hh>

#include <dumux/decoupled/2p/impes/impesproblem2p.hh>
#include <dumux/decoupled/2p/diffusion/fv/fvvelocity2p.hh>
#include <dumux/decoupled/2p/transport/fv/fvsaturation2p.hh>
#include <dumux/decoupled/2p/transport/fv/capillarydiffusion.hh>
#include <dumux/decoupled/2p/transport/fv/gravitypart.hh>

#include "tutorialspatialparameters_decoupled.hh"

namespace Dumux
{

template<class TypeTag>
class TutorialProblemDecoupled;

//////////
// Specify the properties for the lens problem
//////////
namespace Properties
{
NEW_TYPE_TAG(TutorialProblemDecoupled, INHERITS_FROM(DecoupledTwoP, Transport));

// Set the grid type
SET_PROP(TutorialProblemDecoupled, Grid)
{
    typedef Dune::SGrid<2, 2> type;
    static type *create() /*@\label{tutorial-coupled:create-grid-method}@*/
    {
        typedef typename type::ctype ctype;
        Dune::FieldVector<int, 2> cellRes;
        Dune::FieldVector<ctype, 2> lowerLeft(0.0);
        Dune::FieldVector<ctype, 2> upperRight;
        upperRight[0] = 300;
        upperRight[1] = 60;
        cellRes[0] = 100;
        cellRes[1] = 1;
        return new Dune::SGrid<2,2>(cellRes,
                                    lowerLeft,
                                    upperRight);
    } /*@\label{tutorial-decoupled:grid-end}@*/
};

// Set the problem property
SET_PROP(TutorialProblemDecoupled, Problem)
{
public:
    typedef Dumux::TutorialProblemDecoupled<TTAG(TutorialProblemDecoupled)> type;
};

// Set the model properties
SET_PROP(TutorialProblemDecoupled, TransportModel)
{
    typedef Dumux::FVSaturation2P<TTAG(TutorialProblemDecoupled)> type;
};

SET_PROP(TutorialProblemDecoupled, PressureModel)
{
    typedef Dumux::FVVelocity2P<TTAG(TutorialProblemDecoupled)> type;
};

SET_INT_PROP(TutorialProblemDecoupled, VelocityFormulation,
        GET_PROP_TYPE(TypeTag, PTAG(TwoPIndices))::velocityW);

//SET_INT_PROP(TutorialProblemDecoupled, PressureFormulation,
//        GET_PROP_TYPE(TypeTag, PTAG(TwoPIndices))::pressureGlobal);

// Set the wetting phase
SET_PROP(TutorialProblemDecoupled, WettingPhase)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
public:
    typedef Dumux::LiquidPhase<Scalar, Dumux::H2O<Scalar> > type;
};

// Set the non-wetting phase
SET_PROP(TutorialProblemDecoupled, NonwettingPhase)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
public:
    typedef Dumux::LiquidPhase<Scalar, Dumux::Oil<Scalar> > type;
};

// Set the spatial parameters
SET_PROP(TutorialProblemDecoupled, SpatialParameters)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Grid)) Grid;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;

public:
    typedef Dumux::TutorialSpatialParametersDecoupled<TypeTag> type;
};

SET_TYPE_PROP(TutorialProblemDecoupled, DiffusivePart, Dumux::CapillaryDiffusion<TypeTag>);

// Disable gravity
SET_BOOL_PROP(TutorialProblemDecoupled, EnableGravity, false);

SET_SCALAR_PROP(TutorialProblemDecoupled, CFLFactor, 0.3);
}

/*!
* \ingroup DecoupledProblems
*/
template<class TypeTag = TTAG(TutorialProblemDecoupled)>
class TutorialProblemDecoupled: public IMPESProblem2P<TypeTag, TutorialProblemDecoupled<TypeTag> >
{
    typedef TutorialProblemDecoupled<TypeTag> ThisType;
    typedef IMPESProblem2P<TypeTag, ThisType> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(TwoPIndices)) Indices;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluidSystem)) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluidState)) FluidState;

    enum
    {
        dim = GridView::dimension, dimWorld = GridView::dimensionworld
    };

    enum
    {
        wPhaseIdx = Indices::wPhaseIdx, nPhaseIdx = Indices::nPhaseIdx
    };

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;

    typedef typename GridView::Traits::template Codim<0>::Entity Element;
    typedef typename GridView::Intersection Intersection;
    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;
    typedef Dune::FieldVector<Scalar, dim> LocalPosition;

public:
    TutorialProblemDecoupled(const GridView &gridView, const GlobalPosition lowerLeft = GlobalPosition(0.), const GlobalPosition upperRight = GlobalPosition(0.)) :
        ParentType(gridView), lowerLeft_(lowerLeft), upperRight_(upperRight)
    {
    }

    /*!
    * \name Problem parameters
    */
    // \{

    /*!
    * \brief The problem name.
    *
    * This is used as a prefix for files generated by the simulation.
    */
    const char *name() const
    {
        return "tutorial_decoupled";
    }

    bool shouldWriteRestartFile() const
    {
        return false;
    }

    /*!
    * \brief Returns the temperature within the domain.
    *
    * This problem assumes a temperature of 10 degrees Celsius.
    */
    Scalar temperature(const GlobalPosition& globalPos, const Element& element) const
    {
        return 273.15 + 10; // -> 10°C
    }

    // \}

    Scalar referencePressure(const GlobalPosition& globalPos, const Element& element) const
    {
        return 1e5; // -> 10°C
    }

    std::vector<Scalar> source(const GlobalPosition& globalPos, const Element& element)
        {
        return std::vector<Scalar>(2, 0.0);
        }

    typename BoundaryConditions::Flags bctypePress(const GlobalPosition& globalPos, const Intersection& intersection) const
    {
        if ((globalPos[0] < lowerLeft_[0] + eps_))
            return BoundaryConditions::dirichlet;
        // all other boundaries
        return BoundaryConditions::neumann;
    }

    BoundaryConditions::Flags bctypeSat(const GlobalPosition& globalPos, const Intersection& intersection) const
    {
        if (globalPos[0] < lowerLeft_[0] + eps_)
            return Dumux::BoundaryConditions::dirichlet;
        else
            return Dumux::BoundaryConditions::neumann;
    }

    Scalar dirichletPress(const GlobalPosition& globalPos, const Intersection& intersection) const
    {
        if (globalPos[0] < lowerLeft_[0] + eps_)
            return 2e5;
        // all other boundaries
        return 0;
    }

    Scalar dirichletSat(const GlobalPosition& globalPos, const Intersection& intersection) const
    {
        if (globalPos[0] < lowerLeft_[0] + eps_)
            return 1;
        // all other boundaries
        return 0;
    }

    std::vector<Scalar> neumannPress(const GlobalPosition& globalPos, const Intersection& intersection) const
    {
        std::vector<Scalar> neumannFlux(2,0.0);
        if (globalPos[0] > upperRight_[0] - eps_)
        {
            neumannFlux[nPhaseIdx] = 3e-4;
        }
        return neumannFlux;
    }

    Scalar neumannSat(const GlobalPosition& globalPos, const Intersection& intersection, Scalar factor) const
    {
        return 0;
    }

    Scalar initSat(const GlobalPosition& globalPos, const Element& element) const
    {
        return 0;
    }

private:
    GlobalPosition lowerLeft_;
    GlobalPosition upperRight_;

    static const Scalar eps_ = 1e-6;
};
} //end namespace

#endif
