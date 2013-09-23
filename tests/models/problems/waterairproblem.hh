// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2008-2012 by Andreas Lauser                               *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 *
 * \copydoc Ewoms::WaterAirProblem
 */
#ifndef EWOMS_WATER_AIR_PROBLEM_HH
#define EWOMS_WATER_AIR_PROBLEM_HH

#include <ewoms/models/pvs/pvsproperties.hh>

#include <opm/material/fluidsystems/H2OAirFluidSystem.hpp>
#include <opm/material/fluidstates/ImmiscibleFluidState.hpp>
#include <opm/material/fluidstates/CompositionalFluidState.hpp>
#include <opm/material/fluidmatrixinteractions/2p/LinearMaterial.hpp>
#include <opm/material/fluidmatrixinteractions/2p/RegularizedBrooksCorey.hpp>
#include <opm/material/fluidmatrixinteractions/2p/EffToAbsLaw.hpp>
#include <opm/material/fluidmatrixinteractions/mp/2pAdapter.hpp>
#include <opm/material/heatconduction/Somerton.hpp>
#include <opm/material/constraintsolvers/ComputeFromReferencePhase.hpp>

#include <dune/grid/io/file/dgfparser/dgfug.hh>
#include <dune/grid/io/file/dgfparser/dgfs.hh>
#include <dune/grid/io/file/dgfparser/dgfyasp.hh>

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

#include <sstream>
#include <string>

namespace Ewoms {

template <class TypeTag>
class WaterAirProblem;

namespace Properties {

NEW_TYPE_TAG(WaterAirBaseProblem);

// Set the grid type
SET_TYPE_PROP(WaterAirBaseProblem, Grid, Dune::YaspGrid<2>);

// Set the problem property
SET_TYPE_PROP(WaterAirBaseProblem, Problem, Ewoms::WaterAirProblem<TypeTag>);

// Set the material Law
SET_PROP(WaterAirBaseProblem, MaterialLaw)
{
private:
    // define the material law which is parameterized by effective
    // saturations
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef Opm::RegularizedBrooksCorey<Scalar> EffMaterialLaw;

    // define the material law parameterized by absolute saturations
    // which uses the two-phase API
    typedef Opm::EffToAbsLaw<EffMaterialLaw> TwoPMaterialLaw;

    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    enum { lPhaseIdx = FluidSystem::lPhaseIdx };

public:
    // define the type of the generic material law
    typedef Opm::TwoPAdapter<lPhaseIdx, TwoPMaterialLaw> type;
};

// Set the heat conduction law
SET_PROP(WaterAirBaseProblem, HeatConductionLaw)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

public:
    // define the material law parameterized by absolute saturations
    typedef Opm::Somerton<FluidSystem, Scalar> type;
};

// Set the fluid system. in this case, we use the one which describes
// air and water
SET_TYPE_PROP(WaterAirBaseProblem, FluidSystem,
              Opm::FluidSystems::H2OAir<typename GET_PROP_TYPE(TypeTag, Scalar)>);

// Enable gravity
SET_BOOL_PROP(WaterAirBaseProblem, EnableGravity, true);

// Enable constraints
SET_BOOL_PROP(WaterAirBaseProblem, EnableConstraints, true);

// Use forward differences instead of central differences
SET_INT_PROP(WaterAirBaseProblem, NumericDifferenceMethod, +1);

// Write newton convergence
SET_BOOL_PROP(WaterAirBaseProblem, NewtonWriteConvergence, false);

// The default for the end time of the simulation
SET_SCALAR_PROP(WaterAirBaseProblem, EndTime, 5e3);

// The default for the initial time step size of the simulation
SET_SCALAR_PROP(WaterAirBaseProblem, InitialTimeStepSize, 250);

// The default DGF file to load
SET_STRING_PROP(WaterAirBaseProblem, GridFile, "./grids/waterair.dgf");
}


/*!
 * \ingroup VcfvTestProblems
 * \brief Non-isothermal gas injection problem where a air
 *        is injected into a fully water saturated medium.
 *
 * During buoyancy driven upward migration, the gas passes a
 * rectangular high temperature area. This decreases the temperature
 * of the high-temperature area and accelerates gas infiltration due
 * to the lower viscosity of the gas. (Be aware that the pressure of
 * the gas is approximately constant within the lens, so the density
 * of the gas is reduced. This more than off-sets the viscosity
 * increase of the gas at constant density.)
 *
 * The domain is sized 40 m times 40 m. The rectangular area with
 * increased temperature (380 K) starts at (20 m, 5 m) and ends at (30
 * m, 35 m).
 *
 * For the mass conservation equation, no-flow boundary conditions are
 * used on the top and on the bottom of the domain, while free-flow
 * conditions apply on the left and the right boundary. Gas is
 * injected at bottom from 15 m to 25 m at a rate of 0.001 kg/(s m^2)
 * by means if a forced inflow boundary condition.
 *
 * At the free-flow boundaries, the initial condition for the bulk
 * part of the domain is assumed, i. e.  hydrostatic pressure, a gas
 * saturation of zero and a geothermal temperature gradient of 0.03
 * K/m.
 */
template <class TypeTag >
class WaterAirProblem
    : public GET_PROP_TYPE(TypeTag, BaseProblem)
{
    typedef typename GET_PROP_TYPE(TypeTag, BaseProblem) ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GridView::Grid Grid;

    // copy some indices for convenience
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    enum {
        numPhases = FluidSystem::numPhases,

        // energy related indices
        temperatureIdx = Indices::temperatureIdx,
        energyEqIdx = Indices::energyEqIdx,

        // component indices
        H2OIdx = FluidSystem::H2OIdx,
        AirIdx = FluidSystem::AirIdx,

        // phase indices
        lPhaseIdx = FluidSystem::lPhaseIdx,
        gPhaseIdx = FluidSystem::gPhaseIdx,

        // equation indices
        conti0EqIdx = Indices::conti0EqIdx,

        // Grid and world dimension
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld
    };

    typedef typename GET_PROP_TYPE(TypeTag, RateVector) RateVector;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryRateVector) BoundaryRateVector;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, Constraints) Constraints;
    typedef typename GET_PROP_TYPE(TypeTag, TimeManager) TimeManager;

    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLawParams) MaterialLawParams;

    typedef typename GET_PROP_TYPE(TypeTag, HeatConductionLaw) HeatConductionLaw;
    typedef typename GET_PROP_TYPE(TypeTag, HeatConductionLawParams) HeatConductionLawParams;

    typedef typename GridView::ctype CoordScalar;
    typedef Dune::FieldVector<CoordScalar, dimWorld> GlobalPosition;

    typedef Dune::FieldMatrix<Scalar, dimWorld, dimWorld> DimMatrix;

public:
    /*!
     * \copydoc Doxygen::defaultProblemConstructor
     */
    WaterAirProblem(TimeManager &timeManager)
        : ParentType(timeManager, GET_PROP_TYPE(TypeTag, GridCreator)::grid().leafView())
    {
        maxDepth_ = 1000.0; // [m]
        eps_ = 1e-6;

        FluidSystem::init(/*Tmin=*/275, /*Tmax=*/600, /*nT=*/100,
                          /*pmin=*/9.5e6, /*pmax=*/10.5e6, /*np=*/200);

        layerBottom_ = 22.0;

        // intrinsic permeabilities
        fineK_ = this->toDimMatrix_(1e-13);
        coarseK_ = this->toDimMatrix_(1e-12);

        // porosities
        finePorosity_ = 0.3;
        coarsePorosity_ = 0.3;

        // residual saturations
        fineMaterialParams_.setSwr(0.2);
        fineMaterialParams_.setSnr(0.0);
        coarseMaterialParams_.setSwr(0.2);
        coarseMaterialParams_.setSnr(0.0);

        // parameters for the Brooks-Corey law
        fineMaterialParams_.setPe(1e4);
        coarseMaterialParams_.setPe(1e4);
        fineMaterialParams_.setLambda(2.0);
        coarseMaterialParams_.setLambda(2.0);

        // parameters for the somerton law of heat conduction
        computeHeatCondParams_(fineHeatCondParams_, finePorosity_);
        computeHeatCondParams_(coarseHeatCondParams_, coarsePorosity_);
    }

    /*!
     * \name Problem parameters
     */
    //! \{

    /*!
     * \copydoc VcfvProblem::name
     */
    std::string name() const
    {
        std::ostringstream oss;
        oss << "waterair_" << this->model().name();
        if (GET_PROP_VALUE(TypeTag, EnableEnergy))
            oss << "_ni";

        return oss.str();
    }

    /*!
     * \copydoc VcfvMultiPhaseProblem::intrinsicPermeability
     *
     * In this problem, the upper part of the domain is sightly less
     * permeable than the lower one.
     */
    template <class Context>
    const DimMatrix &intrinsicPermeability(const Context &context, int spaceIdx, int timeIdx) const
    {
        const GlobalPosition &pos = context.pos(spaceIdx, timeIdx);
        if (isFineMaterial_(pos))
            return fineK_;
        return coarseK_;
    }

    /*!
     * \copydoc VcfvMultiPhaseProblem::porosity
     */
    template <class Context>
    Scalar porosity(const Context &context, int spaceIdx, int timeIdx) const
    {
        const GlobalPosition &pos = context.pos(spaceIdx, timeIdx);
        if (isFineMaterial_(pos))
            return finePorosity_;
        else
            return coarsePorosity_;
    }

    /*!
     * \copydoc VcfvMultiPhaseProblem::materialLawParams
     */
    template <class Context>
    const MaterialLawParams& materialLawParams(const Context &context, int spaceIdx, int timeIdx) const
    {
        const GlobalPosition &pos = context.pos(spaceIdx, timeIdx);
        if (isFineMaterial_(pos))
            return fineMaterialParams_;
        else
            return coarseMaterialParams_;
    }

    /*!
     * \copydoc VcfvMultiPhaseProblem::heatCapacitySolid
     *
     * In this case, we assume the rock-matrix to be granite.
     */
    template <class Context>
    Scalar heatCapacitySolid(const Context &context, int spaceIdx, int timeIdx) const
    {
        return
            790 // specific heat capacity of granite [J / (kg K)]
            * 2700; // density of granite [kg/m^3]
    }

    /*!
     * \copydoc VcfvMultiPhaseProblem::heatConductionParams
     */
    template <class Context>
    const HeatConductionLawParams&
    heatConductionParams(const Context &context, int spaceIdx, int timeIdx) const
    {
        const GlobalPosition &pos = context.pos(spaceIdx, timeIdx);
        if (isFineMaterial_(pos))
            return fineHeatCondParams_;
        return coarseHeatCondParams_;
    }

    //! \}

    /*!
     * \name Boundary conditions
     */
    //! \{

    /*!
     * \copydoc VcfvProblem::boundary
     *
     * For this problem, we inject air at the inlet on the center of
     * the lower domain boundary and use a no-flow condition on the
     * top boundary and a and a free-flow condition on the left and
     * right boundaries of the domain.
     */
    template <class Context>
    void boundary(BoundaryRateVector &values,
                  const Context &context,
                  int spaceIdx, int timeIdx) const
    {
        const auto &pos = context.cvCenter(spaceIdx, timeIdx);
        assert(onLeftBoundary_(pos) ||
               onLowerBoundary_(pos) ||
               onRightBoundary_(pos) ||
               onUpperBoundary_(pos));

        if (onInlet_(pos)) {
            RateVector massRate(0.0);
            massRate[conti0EqIdx + AirIdx] = -1e-3; // [kg/(m^2 s)]

            // impose an forced inflow boundary condition on the inlet
            values.setMassRate(massRate);
        }
        else if (onLeftBoundary_(pos) || onRightBoundary_(pos)) {
            //int globalIdx = context.elemContext().globalSpaceIndex(context.insideScvIndex(spaceIdx,timeIdx), timeIdx);

            Opm::CompositionalFluidState<Scalar, FluidSystem> fs;
            initialFluidState_(fs, context, spaceIdx, timeIdx);

            // impose an freeflow boundary condition
            values.setFreeFlow(context, spaceIdx, timeIdx, fs);
        }
        else
            // no flow on top and bottom
            values.setNoFlow();
    }

    //! \}

    /*!
     * \name Volume terms
     */
    //! \{

    /*!
     * \copydoc VcfvProblem::initial
     *
     * For this problem, we set the medium to be fully saturated by
     * liquid water and assume hydrostatic pressure.
     */
    template <class Context>
    void initial(PrimaryVariables &values, const Context &context, int spaceIdx, int timeIdx) const
    {
        //int globalIdx = context.globalSpaceIndex(spaceIdx, timeIdx);

        Opm::CompositionalFluidState<Scalar, FluidSystem> fs;
        initialFluidState_(fs, context, spaceIdx, timeIdx);

        const auto &matParams = materialLawParams(context, spaceIdx, timeIdx);
        values.assignMassConservative(fs, matParams, /*inEquilibrium=*/true);
    }

    /*!
     * \copydoc VcfvProblem::constraints
     *
     * In this problem, constraints are used to keep the temperature
     * of the finite-volumes which are closest to the inlet constant.
     */
    template <class Context>
    void constraints(Constraints &constraints,
                     const Context &context,
                     int spaceIdx, int timeIdx) const
    {
        const auto &pos = context.pos(spaceIdx, timeIdx);

        if (onInlet_(pos)) {
            constraints.setConstraint(temperatureIdx, energyEqIdx, 380);;
        }
    }

    /*!
     * \copydoc VcfvProblem::source
     *
     * For this problem, the source term of all components is 0
     * everywhere.
     */
    template <class Context>
    void source(RateVector &rate,
                const Context &context, int spaceIdx, int timeIdx) const
    { rate = 0; }

    //! \}


private:
    bool onLeftBoundary_(const GlobalPosition &pos) const
    { return pos[0] < eps_; }

    bool onRightBoundary_(const GlobalPosition &pos) const
    { return pos[0] > this->bboxMax()[0] - eps_; }

    bool onLowerBoundary_(const GlobalPosition &pos) const
    { return pos[1] < eps_; }

    bool onUpperBoundary_(const GlobalPosition &pos) const
    { return pos[1] > this->bboxMax()[1] - eps_; }

    bool onInlet_(const GlobalPosition &pos) const
    { return onLowerBoundary_(pos) && (15.0 < pos[0]) && (pos[0] < 25.0); }

    bool inHighTemperatureRegion_(const GlobalPosition &pos) const
    { return (20 < pos[0]) && (pos[0] < 30) && (pos[1] < 30); }

    template <class Context, class FluidState>
    void initialFluidState_(FluidState &fs, const Context &context, int spaceIdx, int timeIdx) const
    {
        const GlobalPosition &pos = context.pos(spaceIdx, timeIdx);

        Scalar densityW = 1000.0;
        fs.setPressure(lPhaseIdx, 1e5 + (maxDepth_ - pos[1])*densityW*9.81);
        fs.setSaturation(lPhaseIdx, 1.0);
        fs.setMoleFraction(lPhaseIdx, H2OIdx, 1.0);
        fs.setMoleFraction(lPhaseIdx, AirIdx, 0.0);

        if (inHighTemperatureRegion_(pos))
            fs.setTemperature(380);
        else
            fs.setTemperature(283.0 + (maxDepth_ - pos[1])*0.03);

        // set the gas saturation and pressure
        fs.setSaturation(gPhaseIdx, 0);
        Scalar pc[numPhases];
        const auto &matParams = materialLawParams(context, spaceIdx, timeIdx);
        MaterialLaw::capillaryPressures(pc, matParams, fs);
        fs.setPressure(gPhaseIdx, fs.pressure(lPhaseIdx) + (pc[gPhaseIdx] - pc[lPhaseIdx]));

        typename FluidSystem::ParameterCache paramCache;
        typedef Opm::ComputeFromReferencePhase<Scalar, FluidSystem> CFRP;
        CFRP::solve(fs, paramCache, lPhaseIdx, /*setViscosity=*/false,  /*setEnthalpy=*/true);
    }

    void computeHeatCondParams_(HeatConductionLawParams &params, Scalar poro)
    {
        Scalar lambdaGranite = 2.8; // [W / (K m)]

        // create a Fluid state which has all phases present
        Opm::ImmiscibleFluidState<Scalar, FluidSystem> fs;
        fs.setTemperature(293.15);
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            fs.setPressure(phaseIdx, 1.0135e5);
        }

        typename FluidSystem::ParameterCache paramCache;
        paramCache.updateAll(fs);
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            Scalar rho = FluidSystem::density(fs, paramCache, phaseIdx);
            fs.setDensity(phaseIdx, rho);
        }

        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            Scalar lambdaSaturated;
            if (FluidSystem::isLiquid(phaseIdx)) {
                Scalar lambdaFluid =
                    FluidSystem::thermalConductivity(fs, paramCache, phaseIdx);
                lambdaSaturated = std::pow(lambdaGranite, (1-poro)) + std::pow(lambdaFluid, poro);
            }
            else
                lambdaSaturated = std::pow(lambdaGranite, (1-poro));

            params.setFullySaturatedLambda(phaseIdx, lambdaSaturated);
            if (!FluidSystem::isLiquid(phaseIdx))
                params.setVacuumLambda(lambdaSaturated);
        }
    }

    bool isFineMaterial_(const GlobalPosition &pos) const
    { return pos[dim-1] > layerBottom_; }

    DimMatrix fineK_;
    DimMatrix coarseK_;
    Scalar layerBottom_;

    Scalar finePorosity_;
    Scalar coarsePorosity_;

    MaterialLawParams fineMaterialParams_;
    MaterialLawParams coarseMaterialParams_;

    HeatConductionLawParams fineHeatCondParams_;
    HeatConductionLawParams coarseHeatCondParams_;

    Scalar maxDepth_;
    Scalar eps_;
};
} //end namespace

#endif
