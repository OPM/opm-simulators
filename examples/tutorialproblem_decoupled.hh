#ifndef TUTORIALPROBLEM_DECOUPLED_HH
#define TUTORIALPROBLEM_DECOUPLED_HH

#include "dumux/fractionalflow/fractionalflowproblem.hh"

namespace Dune
{

/** \todo Please doc me! */

template<class GridView, class Scalar, class VC> class TutorialProblemDecoupled /*@\label{tutorial-decoupled:tutorialproblem}@*/
    : public FractionalFlowProblem<GridView, Scalar, VC>
{
    enum
        {dim=GridView::dimension, dimWorld = GridView::dimensionworld};
    typedef typename GridView::Grid Grid;
    typedef typename GridView::Traits::template Codim<0>::Entity Element;
    typedef Dune::FieldVector<Scalar,dim> LocalPosition;
    typedef Dune::FieldVector<Scalar,dimWorld> GlobalPosition;

public:
    TutorialProblemDecoupled(VC& variables, Fluid& wettingphase, Fluid& nonwettingphase, Matrix2p<Grid, Scalar>& soil,
                             TwoPhaseRelations<Grid, Scalar>& materialLaw = *(new TwoPhaseRelations<Grid,Scalar>),
                             const FieldVector<Scalar,dim> Left = 0, const FieldVector<Scalar,dim> Right = 0)
        : FractionalFlowProblem<GridView, Scalar, VC>(variables, wettingphase, nonwettingphase, soil, materialLaw),
          Left_(Left[0]), Right_(Right[0]), eps_(1e-8)
    {}

    // function returning source/sink terms for the pressure equation
    // depending on the position within the domain
    virtual Scalar sourcePress (const GlobalPosition& globalPos, const Element& e, /*@\label{tutorial-decoupled:qpress}@*/
                                const LocalPosition& localPos)
    {
        return 0;
    }

    // function returning the boundary condition type for solution
    // of the pressure equation depending on the position within the domain
    typename BoundaryConditions::Flags bctypePress(const GlobalPosition& globalPos, const Element& e, /*@\label{tutorial-decoupled:bctypepress}@*/
                                                   const LocalPosition& localPos) const
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
    BoundaryConditions::Flags bctypeSat (const GlobalPosition& globalPos, const Element& e, /*@\label{tutorial-decoupled:bctypesat}@*/
                                         const LocalPosition& localPos) const
    {
        if (globalPos[0]> (Right_ - eps_) || globalPos[0] < eps_)
        {
            return Dune::BoundaryConditions::dirichlet;
        }
        // all other boundaries
        return Dune::BoundaryConditions::neumann;
    }

    // function returning the Dirichlet boundary condition for the solution
    // of the pressure equation depending on the position within the domain
    Scalar dirichletPress(const GlobalPosition& globalPos, const Element& e, /*@\label{tutorial-decoupled:gpress}@*/
                          const LocalPosition& localPos) const
    {
        return 2e5;
    }

    // function returning the Dirichlet boundary condition for the solution
    // of the saturation equation depending on the position within the domain
    Scalar dirichletSat(const GlobalPosition& globalPos, const Element& e, /*@\label{tutorial-decoupled:gsat}@*/
                        const LocalPosition& localPos) const
    {
        if (globalPos[0] < eps_)
        {
            return 1;
        }
        // all other boundaries
        return 0;
    }

    // function returning the Neumann boundary condition for the solution
    // of the pressure equation depending on the position within the domain
    Scalar neumannPress(const GlobalPosition& globalPos, const Element& e, /*@\label{tutorial-decoupled:jpress}@*/
                        const LocalPosition& localPos) const
    {
        if (globalPos[0]> Right_ - eps_)
        {
            return 3e-7;
        }
        // all other boundaries
        return 0;
    }

    // function returning the initial saturation
    // depending on the position within the domain
    Scalar initSat (const GlobalPosition& globalPos, const Element& e, /*@\label{tutorial-decoupled:initsat}@*/
                    const FieldVector<Scalar,dim>& xi) const
    {
        return 0;
    }

private:
    Scalar Left_;
    Scalar Right_;

    Scalar eps_;
};
} // end namespace
#endif
