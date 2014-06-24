/*
  Copyright (C) 2012-2013 by Andreas Lauser
  Copyright (C) 2012 by Klaus Mosthaf

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 2 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/
/*!
 * \file
 * \copydoc Ewoms::NavierStokesTestProblem
 */
#ifndef EWOMS_NAVIER_STOKES_TEST_PROBLEM_HH
#define EWOMS_NAVIER_STOKES_TEST_PROBLEM_HH

#include <ewoms/models/stokes/stokesmodel.hh>

#include <opm/material/fluidsystems/H2ON2FluidSystem.hpp>
#include <opm/material/fluidsystems/GasPhase.hpp>
#include <opm/material/components/N2.hpp>

#include <dune/grid/alugrid/2d/alugrid.hh>

#include <dune/common/version.hh>
#include <dune/common/fvector.hh>

namespace Ewoms {
template <class TypeTag>
class NavierStokesTestProblem;
}

namespace Opm {
namespace Properties {
NEW_TYPE_TAG(NavierStokesTestProblem, INHERITS_FROM(NavierStokesModel));

// Set the grid type
SET_TYPE_PROP(NavierStokesTestProblem, Grid, Dune::ALUGrid<2, 2, Dune::cube, Dune::nonconforming>);

// Set the property which defines the type of the physical problem
SET_TYPE_PROP(NavierStokesTestProblem, Problem,
              Ewoms::NavierStokesTestProblem<TypeTag>);

SET_PROP(NavierStokesTestProblem, Fluid)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;

public:
    typedef Opm::GasPhase<Scalar, Opm::N2<Scalar> > type;
};

// Disable gravity
SET_BOOL_PROP(NavierStokesTestProblem, EnableGravity, false);

// Enable constraints
SET_BOOL_PROP(NavierStokesTestProblem, EnableConstraints, true);

// Default simulation end time [s]
SET_SCALAR_PROP(NavierStokesTestProblem, EndTime, 1e-3);

// Default initial time step size [s]
SET_SCALAR_PROP(NavierStokesTestProblem, InitialTimeStepSize, 1e-3);

// Default grid file to load
SET_STRING_PROP(NavierStokesTestProblem, GridFile,
                "data/test_navierstokes.dgf");
}
}

namespace Ewoms {
/*!
 * \ingroup StokesModel
 * \ingroup VcfvTestProblems
 * \brief Stokes flow problem with modified nitrogen (N2) circulating in
 *        a cavity. (lid-driven cavity-flow)
 *
 * The example is taken from Ghia, Ghia, and Shin (1982), "High-Re solutions
 * for incompressible flow using the Navier-Stokes equations and a multigrid
 * method", Journal of Computational Physics, Vol. 48, pp. 387-411.
 *
 * The domain is two-dimensional and sized 1m times 1m. The boundary
 * conditions for the momentum balances are no-flow boundary
 * conditions except for the top, which is floating from left to right
 * with 1 m/s. The mass balance features outflow boundary
 * conditions. All vertices at the bottom, left and right boundaries
 * are constraint to a constant pressure level and zero velocity.
 */
template <class TypeTag>
class NavierStokesTestProblem : public StokesProblem<TypeTag>
{
    typedef StokesProblem<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Simulator) Simulator;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, RateVector) RateVector;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryRateVector) BoundaryRateVector;
    typedef typename GET_PROP_TYPE(TypeTag, Constraints) Constraints;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

    enum {
        // Number of equations and grid dimension
        dimWorld = GridView::dimensionworld,

        // copy some indices for convenience
        pressureIdx = Indices::pressureIdx,
        velocity0Idx = Indices::velocity0Idx,
        conti0EqIdx = Indices::conti0EqIdx,
        momentum0EqIdx = Indices::momentum0EqIdx
    };

    typedef typename GridView::ctype CoordScalar;
    typedef Dune::FieldVector<CoordScalar, dimWorld> GlobalPosition;

public:
    /*!
     * \copydoc Doxygen::defaultProblemConstructor
     */
    NavierStokesTestProblem(Simulator &simulator)
        : ParentType(simulator)
    { eps_ = 1e-6; }

    /*!
     * \name Problem parameters
     */
    //! \{

    /*!
     * \copydoc VcfvProblem::name
     */
    std::string name() const
    { return "navierstokes"; }

    /*!
     * \brief StokesProblem::temperature
     *
     * This problem assumes a constant temperature of 10 degrees Celsius.
     */
    template <class Context>
    Scalar temperature(const Context &context, int spaceIdx, int timeIdx) const
    { return 273.15 + 10; }

    //! \}

    /*!
     * \name Boundary conditions
     */
    //! \{

    /*!
     * \copydoc VcfvProblem::boundary
     */
    template <class Context>
    void boundary(BoundaryRateVector &values, const Context &context,
                  int spaceIdx, int timeIdx) const
    {
        /*        const GlobalPosition &pos = context.pos(spaceIdx, timeIdx);

                values.setOutflow(massBalanceIdx);
                values.setDirichlet(momentumXIdx);
                values.setDirichlet(momentumYIdx);
                // set pressure for all vertices at the bottom
                if (onLowerBoundary_(pos)) {
                    values.setDirichlet(massBalanceIdx);
                }
        */
        values.setNoFlow(context, spaceIdx, timeIdx);
    }

    //! \}

    /*!
     * \name Volumetric terms
     */
    //! \{

    /*!
     * \copydoc VcfvProblem::initial
     */
    template <class Context>
    void initial(PrimaryVariables &values, const Context &context, int spaceIdx,
                 int timeIdx) const
    { initial_(values); }

    /*!
     * \copydoc VcfvProblem::constraints
     *
     * For this problem, we fix the velocity of upper boundary.
     */
    template <class Context>
    void constraints(Constraints &constraints, const Context &context,
                     int spaceIdx, int timeIdx) const
    {
        const auto &pos = context.pos(spaceIdx, timeIdx);

        if (onUpperBoundary_(pos)) {
            // lid moves from left to right
            const Scalar lidVelocity = 1.0;
            constraints.setConstraint(momentum0EqIdx, velocity0Idx + 0,
                                      lidVelocity);
            constraints.setConstraint(momentum0EqIdx + 1, velocity0Idx + 1, 0);
            constraints.setConstraint(conti0EqIdx, pressureIdx, 1e5);
        }
    }

    /*!
     * \copydoc VcfvProblem::source
     */
    template <class Context>
    void source(RateVector &rate, const Context &context, int spaceIdx,
                int timeIdx) const
    { rate = Scalar(0.0); }

    //! \}

private:
    // internal method for the initial condition
    void initial_(PrimaryVariables &priVars) const
    {
        priVars[pressureIdx] = 1e5;
        priVars[velocity0Idx + 0] = 0.0;
        priVars[velocity0Idx + 1] = 0.0;
    }

    bool onLeftBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[0] < this->boundingBoxMin()[0] + eps_; }

    bool onRightBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[0] > this->boundingBoxMax()[0] - eps_; }

    bool onLowerBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[1] < this->boundingBoxMin()[1] + eps_; }

    bool onUpperBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[1] > this->boundingBoxMax()[1] - eps_; }

    Scalar eps_;
};

} // namespace Ewoms

#endif
