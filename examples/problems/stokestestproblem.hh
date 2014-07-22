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
 *
 * \copydoc Ewoms::StokesTestProblem
 */
#ifndef EWOMS_STOKES_TEST_PROBLEM_HH
#define EWOMS_STOKES_TEST_PROBLEM_HH

#include <ewoms/models/stokes/stokesmodel.hh>

#include <opm/material/fluidsystems/H2ON2FluidSystem.hpp>
#include <opm/material/fluidsystems/GasPhase.hpp>

#include <dune/grid/yaspgrid.hh>

#include <dune/common/version.hh>
#include <dune/common/fvector.hh>

namespace Ewoms {
template <class TypeTag>
class StokesTestProblem;
}

namespace Opm {
namespace Properties {
NEW_TYPE_TAG(StokesTestProblem, INHERITS_FROM(StokesModel));

// Set the grid type
SET_TYPE_PROP(StokesTestProblem, Grid, Dune::YaspGrid<2>);

// Set the problem property
SET_TYPE_PROP(StokesTestProblem, Problem, Ewoms::StokesTestProblem<TypeTag>);

// Use the default fluid system of the Stokes model. It requires to
// specify a fluid, though.
SET_PROP(StokesTestProblem, Fluid)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;

public:
    typedef Opm::GasPhase<Scalar, Opm::N2<Scalar> > type;
};

// Disable gravity
SET_BOOL_PROP(StokesTestProblem, EnableGravity, false);

// Enable constraints
SET_BOOL_PROP(StokesTestProblem, EnableConstraints, true);

// Default simulation end time [s]
SET_SCALAR_PROP(StokesTestProblem, EndTime, 10.0);

// Default initial time step size [s]
SET_SCALAR_PROP(StokesTestProblem, InitialTimeStepSize, 10.0);

// Default grid file to load
SET_STRING_PROP(StokesTestProblem, GridFile, "data/test_stokes.dgf");
} // namespace Properties
} // namespace Opm

namespace Ewoms {
/*!
 * \ingroup StokesModel
 * \ingroup TestProblems
 *
 * \brief Stokes flow problem with nitrogen (\f$N_2\f$) flowing
 *        from the left to the right.
 *
 * The domain is sized 1m times 1m. The boundary conditions for the
 * momentum balances are set to outflow on the right boundary and to
 * no-flow at the top and bottom of the domain. For the mass balance
 * equation, outflow boundary conditions are assumed on the right,
 * free-flow on the left and no-flow at the top and bottom boundaries.
 */
template <class TypeTag>
class StokesTestProblem : public GET_PROP_TYPE(TypeTag, BaseProblem)
{
    typedef typename GET_PROP_TYPE(TypeTag, BaseProblem) ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Simulator) Simulator;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    typedef typename GET_PROP_TYPE(TypeTag, EqVector) EqVector;
    typedef typename GET_PROP_TYPE(TypeTag, RateVector) RateVector;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryRateVector) BoundaryRateVector;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, Fluid) Fluid;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Constraints) Constraints;
    typedef typename GET_PROP_TYPE(TypeTag, Model) Model;

    enum {
        // Number of equations and grid dimension
        dimWorld = GridView::dimensionworld,

        // equation indices
        conti0EqIdx = Indices::conti0EqIdx,
        momentum0EqIdx = Indices::momentum0EqIdx,

        // primary variable indices
        velocity0Idx = Indices::velocity0Idx,
        pressureIdx = Indices::pressureIdx
    };

    typedef typename GridView::ctype CoordScalar;
    typedef Dune::FieldVector<CoordScalar, dimWorld> GlobalPosition;
    typedef Dune::FieldVector<Scalar, dimWorld> DimVector;

public:
    /*!
     * \copydoc Doxygen::defaultProblemConstructor
     */
    StokesTestProblem(Simulator &simulator)
        : ParentType(simulator)
    { eps_ = 1e-6; }

    /*!
     * \name Problem parameters
     */
    //! \{

    /*!
     * \copydoc FvBaseProblem::name
     */
    std::string name() const
    { return "stokestest"; }


    /*!
     * \copydoc FvBaseProblem::endTimeStep
     */
    void endTimeStep()
    {
#ifndef NDEBUG
        // checkConservativeness() does not include the effect of constraints, so we
        // disable it for this problem...
        //this->model().checkConservativeness();

        // Calculate storage terms
        EqVector storage;
        this->model().globalStorage(storage);

        // Write mass balance information for rank 0
        if (this->gridView().comm().rank() == 0) {
            std::cout << "Storage: " << storage << std::endl << std::flush;
        }
#endif // NDEBUG
    }

    /*!
     * \brief StokesProblem::temperature
     *
     * This problem assumes a constant temperature of 10 degrees Celsius.
     */
    template <class Context>
    Scalar temperature(const Context &context, int spaceIdx, int timeIdx) const
    { return 273.15 + 10; } // -> 10 deg C

    //! \}

    /*!
     * \name Boundary conditions
     */
    //! \{

    /*!
     * \copydoc FvBaseProblem::boundary
     *
     * For this problem, we use an out-flow boundary on the right,
     * no-flow at the top and at the bottom and the left boundary gets
     * a parabolic velocity profile via constraints.
     */
    template <class Context>
    void boundary(BoundaryRateVector &values, const Context &context,
                  int spaceIdx, int timeIdx) const
    {
        const GlobalPosition &pos = context.pos(spaceIdx, timeIdx);

        Scalar y = pos[1] - this->boundingBoxMin()[1];
        Scalar height = this->boundingBoxMax()[1] - this->boundingBoxMin()[1];

        // parabolic velocity profile
        const Scalar maxVelocity = 1.0;

        Scalar a = -4 * maxVelocity / (height * height);
        Scalar b = -a * height;
        Scalar c = 0;

        DimVector velocity(0.0);
        velocity[0] = a * y * y + b * y + c;

        if (onRightBoundary_(pos))
            values.setOutFlow(context, spaceIdx, timeIdx);
        else if (onLeftBoundary_(pos)) {
            // left boundary is constraint!
            values = 0.0;
        }
        else {
            // top and bottom
            values.setNoFlow(context, spaceIdx, timeIdx);
        }
    }

    //! \}

    /*!
     * \name Volumetric terms
     */
    //! \{

    /*!
     * \copydoc FvBaseProblem::initial
     */
    template <class Context>
    void initial(PrimaryVariables &values, const Context &context, int spaceIdx,
                 int timeIdx) const
    {
        const auto &pos = context.pos(spaceIdx, timeIdx);

        Scalar y = pos[1] - this->boundingBoxMin()[1];
        Scalar height = this->boundingBoxMax()[1] - this->boundingBoxMin()[1];

        // parabolic velocity profile on boundaries
        const Scalar maxVelocity = 1.0;

        Scalar a = -4 * maxVelocity / (height * height);
        Scalar b = -a * height;
        Scalar c = 0;

        DimVector velocity(0.0);
        velocity[0] = a * y * y + b * y + c;

        for (int axisIdx = 0; axisIdx < dimWorld; ++axisIdx)
            values[velocity0Idx + axisIdx] = velocity[axisIdx];
        values[pressureIdx] = 1e5;
    }

    /*!
     * \copydoc FvBaseProblem::source
     *
     * For this problem, the source term of all conserved quantities
     * is 0 everywhere.
     */
    template <class Context>
    void source(RateVector &rate, const Context &context, int spaceIdx,
                int timeIdx) const
    { rate = Scalar(0.0); }

    /*!
     * \copydoc FvBaseProblem::constraints
     *
     * For this problem, the left side of the domain gets a parabolic
     * velocity profile using constraints.
     */
    template <class Context>
    void constraints(Constraints &constraints, const Context &context,
                     int spaceIdx, int timeIdx) const
    {
        const auto &pos = context.pos(spaceIdx, timeIdx);

        if (onLeftBoundary_(pos) || onRightBoundary_(pos)) {
            PrimaryVariables initCond;
            initial(initCond, context, spaceIdx, timeIdx);

            constraints.setConstraint(pressureIdx, conti0EqIdx,
                                      initCond[pressureIdx]);
            ;
            for (int axisIdx = 0; axisIdx < dimWorld; ++axisIdx)
                constraints.setConstraint(velocity0Idx + axisIdx,
                                          momentum0EqIdx + axisIdx,
                                          initCond[velocity0Idx + axisIdx]);
        }
    }

    //! \}

private:
    bool onLeftBoundary_(const GlobalPosition &pos) const
    { return pos[0] < this->boundingBoxMin()[0] + eps_; }

    bool onRightBoundary_(const GlobalPosition &pos) const
    { return pos[0] > this->boundingBoxMax()[0] - eps_; }

    bool onLowerBoundary_(const GlobalPosition &pos) const
    { return pos[1] < this->boundingBoxMin()[1] + eps_; }

    bool onUpperBoundary_(const GlobalPosition &pos) const
    { return pos[1] > this->boundingBoxMax()[1] - eps_; }

    bool onBoundary_(const GlobalPosition &pos) const
    {
        return onLeftBoundary_(pos) || onRightBoundary_(pos)
               || onLowerBoundary_(pos) || onUpperBoundary_(pos);
    }

    Scalar eps_;
};

} // namespace Ewoms

#endif
