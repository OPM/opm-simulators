// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
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

  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the list of
  copyright holders.
*/
/*!
 * \file
 * \copydoc Ewoms::StokesNiTestProblem
 */
#ifndef EWOMS_STOKES_NI_TEST_PROBLEM_HH
#define EWOMS_STOKES_NI_TEST_PROBLEM_HH

#include <ewoms/models/stokes/stokesmodel.hh>
#include <ewoms/io/simplexgridmanager.hh>

#include <opm/material/fluidsystems/H2OAirFluidSystem.hpp>
#include <opm/material/common/Unused.hpp>

#include <dune/grid/yaspgrid.hh>
#include <dune/grid/io/file/dgfparser/dgfyasp.hh>

#include <dune/common/version.hh>
#include <dune/common/fvector.hh>

namespace Ewoms {
template <class TypeTag>
class StokesNiTestProblem;
}

namespace Ewoms {
namespace Properties {
NEW_TYPE_TAG(StokesNiTestProblem, INHERITS_FROM(StokesModel));

// Set the grid type
SET_TYPE_PROP(StokesNiTestProblem, Grid, Dune::YaspGrid<2>);

// Set the problem property
SET_TYPE_PROP(StokesNiTestProblem, Problem, Ewoms::StokesNiTestProblem<TypeTag>);

//! Select the fluid system
SET_TYPE_PROP(StokesNiTestProblem, FluidSystem,
              Opm::FluidSystems::H2OAir<typename GET_PROP_TYPE(TypeTag, Scalar)>);

//! Select the phase to be considered
SET_INT_PROP(StokesNiTestProblem, StokesPhaseIndex,
             GET_PROP_TYPE(TypeTag, FluidSystem)::gasPhaseIdx);

// Enable gravity
SET_BOOL_PROP(StokesNiTestProblem, EnableGravity, true);

// Enable the energy equation
SET_BOOL_PROP(StokesNiTestProblem, EnableEnergy, true);

// Enable constraints
SET_BOOL_PROP(StokesNiTestProblem, EnableConstraints, true);

// Default simulation end time [s]
SET_SCALAR_PROP(StokesNiTestProblem, EndTime, 3.0);

// Default initial time step size [s]
SET_SCALAR_PROP(StokesNiTestProblem, InitialTimeStepSize, 0.1);

// Increase the default raw tolerance of the Newton-Raphson method to 10^-4
SET_SCALAR_PROP(StokesNiTestProblem, NewtonRawTolerance, 1e-4);

// Default grid file to load
SET_STRING_PROP(StokesNiTestProblem, GridFile, "data/test_stokes2cni.dgf");
} // namespace Properties
} // namespace Ewoms

namespace Ewoms {
/*!
 * \ingroup StokesNiModel
 * \ingroup TestProblems
 * \brief Non-isothermal test problem for the Stokes model with a gas
 *        (N2) flowing from the left to the right.
 *
 * The domain of this problem is 1m times 1m. The upper and the lower
 * boundaries are fixed to the initial condition by means of
 * constraints, the left and the right boundaries are no-slip
 * conditions.
 */
template <class TypeTag>
class StokesNiTestProblem : public GET_PROP_TYPE(TypeTag, BaseProblem)
{
    typedef typename GET_PROP_TYPE(TypeTag, BaseProblem) ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Simulator) Simulator;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, Constraints) Constraints;
    typedef typename GET_PROP_TYPE(TypeTag, EqVector) EqVector;
    typedef typename GET_PROP_TYPE(TypeTag, RateVector) RateVector;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryRateVector) BoundaryRateVector;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;

    enum {
        // Number of equations and grid dimension

        dimWorld = GridView::dimensionworld,

        // primary variable indices
        pressureIdx = Indices::pressureIdx,
        moleFrac1Idx = Indices::moleFrac1Idx,
        velocity0Idx = Indices::velocity0Idx,
        temperatureIdx = Indices::temperatureIdx,

        // equation indices
        conti0EqIdx = Indices::conti0EqIdx,
        momentum0EqIdx = Indices::momentum0EqIdx,
        energyEqIdx = Indices::energyEqIdx
    };
    enum { numComponents = FluidSystem::numComponents };
    enum { H2OIdx = FluidSystem::H2OIdx };
    enum { AirIdx = FluidSystem::AirIdx };

    typedef typename GridView::ctype CoordScalar;
    typedef Dune::FieldVector<CoordScalar, dimWorld> GlobalPosition;
    typedef Dune::FieldVector<Scalar, dimWorld> DimVector;

public:
    /*!
     * \copydoc Doxygen::defaultProblemConstructor
     */
    StokesNiTestProblem(Simulator& simulator)
        : ParentType(simulator)
    { }

    /*!
     * \copydoc FvBaseProblem::finishInit
     */
    void finishInit()
    {
        ParentType::finishInit();

        eps_ = 1e-6;

        // initialize the tables of the fluid system
        FluidSystem::init(/*Tmin=*/280.0, /*Tmax=*/285, /*nT=*/10,
                          /*pmin=*/1e5, /*pmax=*/1e5 + 100, /*np=*/200);
    }

    /*!
     * \name Problem parameters
     */
    //! \{

    /*!
     * \copydoc FvBaseProblem::name
     */
    std::string name() const
    { return "stokestest_ni"; }

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

    //! \}

    /*!
     * \name Boundary conditions
     */
    //! \{

    /*!
     * \copydoc FvBaseProblem::boundary
     */
    template <class Context>
    void boundary(BoundaryRateVector& values, const Context& context,
                  unsigned spaceIdx, unsigned timeIdx) const
    {
        const GlobalPosition& pos = context.pos(spaceIdx, timeIdx);

        if (onUpperBoundary_(pos))
            values.setOutFlow(context, spaceIdx, timeIdx);
        else if (onLowerBoundary_(pos)) {
            // lower boundary is constraint!
            values = 0.0;
        }
        else {
            // left and right
            values.setNoFlow(context, spaceIdx, timeIdx);
        }
    }

    //! \}

    /*!
     * \name Volumetric terms
     */
    // \{

    /*!
     * \copydoc FvBaseProblem::initial
     */
    template <class Context>
    void initial(PrimaryVariables& values, const Context& context, unsigned spaceIdx,
                 unsigned timeIdx) const
    {
        const GlobalPosition& pos = context.pos(spaceIdx, timeIdx);

        Scalar moleFrac[numComponents];

        moleFrac[H2OIdx] = 1e-4;
        Scalar temperature = 283.15;
        if (inLens_(pos)) {
            moleFrac[H2OIdx] = 0.9e-4;
            temperature = 284.15;
        }
        moleFrac[AirIdx] = 1 - moleFrac[H2OIdx];

        // parabolic velocity profile
        Scalar y = this->boundingBoxMax()[1] - pos[1];
        Scalar x = pos[0] - this->boundingBoxMin()[0];
        Scalar width = this->boundingBoxMax()[0] - this->boundingBoxMin()[0];

        // parabolic velocity profile
        const Scalar maxVelocity = 1.0;

        Scalar a = -4 * maxVelocity / (width * width);
        Scalar b = -a * width;
        Scalar c = 0;

        DimVector velocity(0.0);
        velocity[1] = a * x * x + b * x + c;

        // hydrostatic pressure
        Scalar rho = 1.189;
        Scalar pressure = 1e5 - rho * this->gravity()[1] * y;

        for (unsigned axisIdx = 0; axisIdx < dimWorld; ++axisIdx)
            values[velocity0Idx + axisIdx] = velocity[axisIdx];

        values[pressureIdx] = pressure;
        values[moleFrac1Idx] = moleFrac[1];
        values[temperatureIdx] = temperature;
    }

    /*!
     * \copydoc FvBaseProblem::source
     *
     * For this problem, the source term of all conserved quantities
     * is 0 everywhere.
     */
    template <class Context>
    void source(RateVector& rate,
                const Context& OPM_UNUSED context,
                unsigned OPM_UNUSED spaceIdx,
                unsigned OPM_UNUSED timeIdx) const
    { rate = Scalar(0.0); }

    /*!
     * \copydoc FvBaseProblem::constraints
     *
     * This problem sets temperature constraints for the finite volumes
     * adjacent to the inlet.
     */
    template <class Context>
    void constraints(Constraints& constraints,
                     const Context& context,
                     unsigned spaceIdx,
                     unsigned timeIdx) const
    {
        const auto& pos = context.pos(spaceIdx, timeIdx);

        if (onLowerBoundary_(pos) || onUpperBoundary_(pos)) {
            constraints.setActive(true);
            initial(constraints, context, spaceIdx, timeIdx);
        }
    }

    //! \}

private:
    bool onLeftBoundary_(const GlobalPosition& pos) const
    { return pos[0] < this->boundingBoxMin()[0] + eps_; }

    bool onRightBoundary_(const GlobalPosition& pos) const
    { return pos[0] > this->boundingBoxMax()[0] - eps_; }

    bool onLowerBoundary_(const GlobalPosition& pos) const
    { return pos[1] < this->boundingBoxMin()[1] + eps_; }

    bool onUpperBoundary_(const GlobalPosition& pos) const
    { return pos[1] > this->boundingBoxMax()[1] - eps_; }

    bool onBoundary_(const GlobalPosition& pos) const
    {
        return onLeftBoundary_(pos) || onRightBoundary_(pos)
               || onLowerBoundary_(pos) || onUpperBoundary_(pos);
    }

    bool inLens_(const GlobalPosition& pos) const
    { return pos[0] < 0.75 && pos[0] > 0.25 && pos[1] < 0.75 && pos[1] > 0.25; }

    Scalar eps_;
};
} // namespace Ewoms

#endif
