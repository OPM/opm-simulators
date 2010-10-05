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

// the grid includes
#if HAVE_UG
#include <dune/grid/uggrid.hh>
#endif
#include <dune/grid/yaspgrid.hh>
#include <dune/grid/sgrid.hh>

// dumux 2p-decoupled environment
#include <dumux/decoupled/2p/impes/impesproblem2p.hh> /*@\label{tutorial-decoupled:parent-problem}@*/
#include <dumux/decoupled/2p/diffusion/fv/fvvelocity2p.hh>
#include <dumux/decoupled/2p/transport/fv/fvsaturation2p.hh>
#include <dumux/decoupled/2p/transport/fv/capillarydiffusion.hh>

// assign parameters dependent on space (e.g. spatial parameters)
#include "tutorialspatialparameters_decoupled.hh" /*@\label{tutorial-decoupled:spatialparameters}@*/

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
NEW_TYPE_TAG(TutorialProblemDecoupled, INHERITS_FROM(DecoupledTwoP)); /*@\label{tutorial-decoupled:create-type-tag}@*/

// Set the problem property
SET_PROP(TutorialProblemDecoupled, Problem) /*@\label{tutorial-decoupled:set-problem}@*/
{
public:
    typedef Dumux::TutorialProblemDecoupled<TTAG(TutorialProblemDecoupled)> type;
};

// Set the grid type
SET_PROP(TutorialProblemDecoupled, Grid) /*@\label{tutorial-decoupled:grid-begin}@*/
{
    typedef Dune::SGrid<2, 2> type;
    static type *create() /*@\label{tutorial-decoupled:create-grid-method}@*/
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

// Set the wetting phase
SET_PROP(TutorialProblemDecoupled, WettingPhase) /*@\label{tutorial-decoupled:2p-system-start}@*/
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
public:
    typedef Dumux::LiquidPhase<Scalar, Dumux::H2O<Scalar> > type; /*@\label{tutorial-decoupled:wettingPhase}@*/
};

// Set the non-wetting phase
SET_PROP(TutorialProblemDecoupled, NonwettingPhase)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
public:
    typedef Dumux::LiquidPhase<Scalar, Dumux::Oil<Scalar> > type; /*@\label{tutorial-decoupled:nonwettingPhase}@*/
}; /*@\label{tutorial-decoupled:2p-system-end}@*/

// Set the spatial parameters
SET_PROP(TutorialProblemDecoupled, SpatialParameters) /*@\label{tutorial-decoupled:set-spatialparameters}@*/
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Grid)) Grid;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;

public:
    typedef Dumux::TutorialSpatialParametersDecoupled<TypeTag> type;
};

// Set the model properties
SET_PROP(TutorialProblemDecoupled, TransportModel) /*@\label{tutorial-decoupled:TransportModel}@*/
{
    typedef Dumux::FVSaturation2P<TTAG(TutorialProblemDecoupled)> type;
};

SET_PROP(TutorialProblemDecoupled, PressureModel) /*@\label{tutorial-decoupled:PressureModel}@*/
{
    typedef Dumux::FVVelocity2P<TTAG(TutorialProblemDecoupled)> type;
};

// model-specific settings
SET_INT_PROP(TutorialProblemDecoupled, VelocityFormulation,
        GET_PROP_TYPE(TypeTag, PTAG(TwoPIndices))::velocityW); /*@\label{tutorial-decoupled:velocityFormulation}@*/


SET_TYPE_PROP(TutorialProblemDecoupled, DiffusivePart,
        Dumux::CapillaryDiffusion<TypeTag>); /*@\label{tutorial-decoupled:DiffusivePart}@*/

SET_SCALAR_PROP(TutorialProblemDecoupled, CFLFactor, 0.3); /*@\label{tutorial-decoupled:cfl}@*/

// Disable gravity
SET_BOOL_PROP(TutorialProblemDecoupled, EnableGravity, false); /*@\label{tutorial-decoupled:gravity}@*/
} /*@\label{tutorial-decoupled:propertysystem-end}@*/

/*!
* \ingroup DecoupledProblems
*/
template<class TypeTag = TTAG(TutorialProblemDecoupled)>
class TutorialProblemDecoupled: public IMPESProblem2P<TypeTag, TutorialProblemDecoupled<TypeTag> > /*@\label{tutorial-decoupled:def-problem}@*/
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
    TutorialProblemDecoupled(const GridView &gridView, const GlobalPosition lowerLeft = GlobalPosition(0.),
            const GlobalPosition upperRight = GlobalPosition(0.)) : ParentType(gridView) /*@\label{tutorial-decoupled:constructor-problem}@*/
    {    }

    /*!
    * \brief The problem name.
    *
    * This is used as a prefix for files generated by the simulation.
    */
    const char *name() const    /*@\label{tutorial-decoupled:name}@*/
    {
        return "tutorial_decoupled";
    }

    /*!
     * \brief Returns true if a restart file should be written.
     *
     * The default behaviour is to write no restart file.
     */
    bool shouldWriteRestartFile() const /*@\label{tutorial-decoupled:restart}@*/
    {
        return false;
    }

    /*!
     * \brief Returns true if the current solution should be written to
     *        disk (i.e. as a VTK file)
     *
     * The default behaviour is to write out every the solution for
     * very time step. Else, change divisor.
     */
    bool shouldWriteOutput() const /*@\label{tutorial-decoupled:output}@*/
    {
        return this->timeManager().timeStepIndex() > 0 &&
        (this->timeManager().timeStepIndex() % 1 == 0);
    }

    /*!
    * \brief Returns the temperature within the domain.
    *
    * This problem assumes a temperature of 10 degrees Celsius.
    */
    Scalar temperature(const GlobalPosition& globalPos, const Element& element) const /*@\label{tutorial-decoupled:temperature}@*/
    {
        return 273.15 + 10; // -> 10°C
    }

    /*!
    * \brief Returns a constant pressure to enter material laws
    *
    * For incrompressible simulations, a constant pressure is necessary
    * to enter the material laws to gain a constant density etc.
    */
    Scalar referencePressure(const GlobalPosition& globalPos, const Element& element) const /*@\label{tutorial-decoupled:refPressure}@*/
    {
        return 2e5;
    }
    /*!
    * \brief Source of mass \f$ [\frac{kg}{m^3 \cdot s}] \f$
    *
    *     Evaluate the source term for all phases within a given
    *     volume. The method returns the mass generated (positive) or
    *     annihilated (negative) per volume unit.
    */
    std::vector<Scalar> source(const GlobalPosition& globalPos, const Element& element) /*@\label{tutorial-decoupled:source}@*/
        {
        return std::vector<Scalar>(2, 0.);
        }

    /*!
     * \brief Type of pressure boundary condition.
     *
     * Defines the type the boundary condition for the pressure equation,
     * either pressure (dirichlet) or flux (neumann).
     */
    typename BoundaryConditions::Flags bctypePress(const GlobalPosition& globalPos, const Intersection& intersection) const /*@\label{tutorial-decoupled:bctypePress}@*/
    {
        if ((globalPos[0] < lowerLeft_[0] + eps_))
            return BoundaryConditions::dirichlet;
        // all other boundaries
        return BoundaryConditions::neumann;
    }

    /*!
     * \brief Type of Transport boundary condition.
     *
     * Defines the type the boundary condition for the transport equation,
     * either saturation (dirichlet) or flux (neumann).
     */
    BoundaryConditions::Flags bctypeSat(const GlobalPosition& globalPos, const Intersection& intersection) const /*@\label{tutorial-decoupled:bctypeSat}@*/
    {
        if (globalPos[0] < lowerLeft_[0] + eps_)
            return Dumux::BoundaryConditions::dirichlet;
        else
            return Dumux::BoundaryConditions::neumann;
    }
    /*!
     * \brief Value for dirichlet pressure boundary condition \f$ [Pa] \f$.
     *
     * In case of a dirichlet BC for the pressure equation, the pressure
     * have to be defined on boundaries.
     */
    Scalar dirichletPress(const GlobalPosition& globalPos, const Intersection& intersection) const /*@\label{tutorial-decoupled:dirichletPress}@*/
    {
        if (globalPos[0] < lowerLeft_[0] + eps_)
            return 2e5;
        // all other boundaries
        return 0;
    }
    /*!
     * \brief Value for transport dirichlet boundary condition (dimensionless).
     *
     * In case of a dirichlet BC for the transport equation, a saturation
     * have to be defined on boundaries.
     */
    Scalar dirichletSat(const GlobalPosition& globalPos, const Intersection& intersection) const /*@\label{tutorial-decoupled:dirichletSat}@*/
    {
        if (globalPos[0] < lowerLeft_[0] + eps_)
            return 1;
        // all other boundaries
        return 0;
    }
    //! Value for pressure neumann boundary condition \f$ [\frac{kg}{m^3 \cdot s}] \f$.
     /** In case of a neumann boundary condition, the flux of matter
      * is returned as a vector.
      */
    std::vector<Scalar> neumannPress(const GlobalPosition& globalPos, const Intersection& intersection) const /*@\label{tutorial-decoupled:neumannPress}@*/
    {
        std::vector<Scalar> neumannFlux(2,0.0);
        if (globalPos[0] > upperRight_[0] - eps_)
        {
            neumannFlux[nPhaseIdx] = 3e-4;
        }
        return neumannFlux;
    }
    //! Value for transport neumann boundary condition \f$ [\frac{kg}{m^3 \cdot s}] \f$.
     /** In case of a neumann boundary condition for the transport equation
      * the flux of matter for the primary variable is returned as a scalar.
      */
    Scalar neumannSat(const GlobalPosition& globalPos, const Intersection& intersection, Scalar factor) const /*@\label{tutorial-decoupled:neumannSat}@*/
    {
        return 0;
    }
    //! Saturation initial condition (dimensionless)
    /*
     * @param element reference to the cell for which the function is to be evaluated
     * @param localPos local coordinates inside element
     */
    Scalar initSat(const GlobalPosition& globalPos, const Element& element) const /*@\label{tutorial-decoupled:initSat}@*/
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
