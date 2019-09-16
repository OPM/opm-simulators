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
 *
 * \copydoc Opm::WaterAirProblem
 */
#ifndef EWOMS_WATER_AIR_PROBLEM_HH
#define EWOMS_WATER_AIR_PROBLEM_HH

#include <opm/models/pvs/pvsproperties.hh>
#include <opm/simulators/linalg/parallelistlbackend.hh>

#include <opm/material/fluidsystems/H2OAirFluidSystem.hpp>
#include <opm/material/fluidstates/ImmiscibleFluidState.hpp>
#include <opm/material/fluidstates/CompositionalFluidState.hpp>
#include <opm/material/fluidmatrixinteractions/LinearMaterial.hpp>
#include <opm/material/fluidmatrixinteractions/RegularizedBrooksCorey.hpp>
#include <opm/material/fluidmatrixinteractions/EffToAbsLaw.hpp>
#include <opm/material/fluidmatrixinteractions/MaterialTraits.hpp>
#include <opm/material/thermal/ConstantSolidHeatCapLaw.hpp>
#include <opm/material/thermal/SomertonThermalConductionLaw.hpp>
#include <opm/material/constraintsolvers/ComputeFromReferencePhase.hpp>
#include <opm/material/common/Unused.hpp>

#include <dune/grid/yaspgrid.hh>
#include <dune/grid/io/file/dgfparser/dgfyasp.hh>

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/version.hh>

#include <sstream>
#include <string>

namespace Opm {
template <class TypeTag>
class WaterAirProblem;
}

BEGIN_PROPERTIES

NEW_TYPE_TAG(WaterAirBaseProblem);

// Set the grid type
SET_TYPE_PROP(WaterAirBaseProblem, Grid, Dune::YaspGrid<2>);

// Set the problem property
SET_TYPE_PROP(WaterAirBaseProblem, Problem, Opm::WaterAirProblem<TypeTag>);

// Set the material Law
SET_PROP(WaterAirBaseProblem, MaterialLaw)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef Opm::TwoPhaseMaterialTraits<Scalar,
                                        /*wettingPhaseIdx=*/FluidSystem::liquidPhaseIdx,
                                        /*nonWettingPhaseIdx=*/FluidSystem::gasPhaseIdx> Traits;

    // define the material law which is parameterized by effective
    // saturations
    typedef Opm::RegularizedBrooksCorey<Traits> EffMaterialLaw;

public:
    // define the material law parameterized by absolute saturations
    // which uses the two-phase API
    typedef Opm::EffToAbsLaw<EffMaterialLaw> type;
};

// Set the thermal conduction law
SET_PROP(WaterAirBaseProblem, ThermalConductionLaw)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

public:
    // define the material law parameterized by absolute saturations
    typedef Opm::SomertonThermalConductionLaw<FluidSystem, Scalar> type;
};

// set the energy storage law for the solid phase
SET_TYPE_PROP(WaterAirBaseProblem, SolidEnergyLaw,
              Opm::ConstantSolidHeatCapLaw<typename GET_PROP_TYPE(TypeTag, Scalar)>);

// Set the fluid system. in this case, we use the one which describes
// air and water
SET_TYPE_PROP(WaterAirBaseProblem, FluidSystem,
              Opm::H2OAirFluidSystem<typename GET_PROP_TYPE(TypeTag, Scalar)>);

// Enable gravity
SET_BOOL_PROP(WaterAirBaseProblem, EnableGravity, true);

// Use forward differences instead of central differences
SET_INT_PROP(WaterAirBaseProblem, NumericDifferenceMethod, +1);

// Write newton convergence
SET_BOOL_PROP(WaterAirBaseProblem, NewtonWriteConvergence, false);

// The default for the end time of the simulation (1 year)
SET_SCALAR_PROP(WaterAirBaseProblem, EndTime, 1.0 * 365 * 24 * 60 * 60);

// The default for the initial time step size of the simulation
SET_SCALAR_PROP(WaterAirBaseProblem, InitialTimeStepSize, 250);

// The default DGF file to load
SET_STRING_PROP(WaterAirBaseProblem, GridFile, "./data/waterair.dgf");

// Use the restarted GMRES linear solver with the ILU-2 preconditioner from dune-istl
SET_TAG_PROP(WaterAirBaseProblem, LinearSolverSplice, ParallelIstlLinearSolver);
SET_TYPE_PROP(WaterAirBaseProblem, LinearSolverWrapper,
              Opm::Linear::SolverWrapperRestartedGMRes<TypeTag>);
#if DUNE_VERSION_NEWER(DUNE_ISTL, 2,7)
SET_TYPE_PROP(WaterAirBaseProblem, PreconditionerWrapper,
              Opm::Linear::PreconditionerWrapperILU<TypeTag>);
#else
SET_TYPE_PROP(WaterAirBaseProblem, PreconditionerWrapper,
              Opm::Linear::PreconditionerWrapperILUn<TypeTag>);
#endif
SET_INT_PROP(WaterAirBaseProblem, PreconditionerOrder, 2);

END_PROPERTIES

namespace Opm {
/*!
 * \ingroup TestProblems
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
class WaterAirProblem : public GET_PROP_TYPE(TypeTag, BaseProblem)
{
    typedef typename GET_PROP_TYPE(TypeTag, BaseProblem) ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;

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
        liquidPhaseIdx = FluidSystem::liquidPhaseIdx,
        gasPhaseIdx = FluidSystem::gasPhaseIdx,

        // equation indices
        conti0EqIdx = Indices::conti0EqIdx,

        // Grid and world dimension
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld
    };

    static const bool enableEnergy = GET_PROP_VALUE(TypeTag, EnableEnergy);

    typedef typename GET_PROP_TYPE(TypeTag, EqVector) EqVector;
    typedef typename GET_PROP_TYPE(TypeTag, RateVector) RateVector;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryRateVector) BoundaryRateVector;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, Constraints) Constraints;
    typedef typename GET_PROP_TYPE(TypeTag, Simulator) Simulator;
    typedef typename GET_PROP_TYPE(TypeTag, Model) Model;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLawParams) MaterialLawParams;
    typedef typename GET_PROP_TYPE(TypeTag, ThermalConductionLawParams) ThermalConductionLawParams;
    typedef typename GET_PROP_TYPE(TypeTag, SolidEnergyLawParams) SolidEnergyLawParams;

    typedef typename GridView::ctype CoordScalar;
    typedef Dune::FieldVector<CoordScalar, dimWorld> GlobalPosition;

    typedef Dune::FieldMatrix<Scalar, dimWorld, dimWorld> DimMatrix;

public:
    /*!
     * \copydoc Doxygen::defaultProblemConstructor
     */
    WaterAirProblem(Simulator& simulator)
        : ParentType(simulator)
    { }

    /*!
     * \copydoc FvBaseProblem::finishInit
     */
    void finishInit()
    {
        ParentType::finishInit();

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
        fineMaterialParams_.setResidualSaturation(liquidPhaseIdx, 0.2);
        fineMaterialParams_.setResidualSaturation(gasPhaseIdx, 0.0);
        coarseMaterialParams_.setResidualSaturation(liquidPhaseIdx, 0.2);
        coarseMaterialParams_.setResidualSaturation(gasPhaseIdx, 0.0);

        // parameters for the Brooks-Corey law
        fineMaterialParams_.setEntryPressure(1e4);
        coarseMaterialParams_.setEntryPressure(1e4);
        fineMaterialParams_.setLambda(2.0);
        coarseMaterialParams_.setLambda(2.0);

        fineMaterialParams_.finalize();
        coarseMaterialParams_.finalize();

        // parameters for the somerton law of thermal conduction
        computeThermalCondParams_(fineThermalCondParams_, finePorosity_);
        computeThermalCondParams_(coarseThermalCondParams_, coarsePorosity_);

        // assume constant volumetric heat capacity and granite
        solidEnergyLawParams_.setSolidHeatCapacity(790.0 // specific heat capacity of granite [J / (kg K)]
                                                   * 2700.0); // density of granite [kg/m^3]
        solidEnergyLawParams_.finalize();
    }

    /*!
     * \name Problem parameters
     */
    //! \{

    /*!
     * \copydoc FvBaseProblem::name
     */
    std::string name() const
    {
        std::ostringstream oss;
        oss << "waterair_" << Model::name();
        if (GET_PROP_VALUE(TypeTag, EnableEnergy))
            oss << "_ni";

        return oss.str();
    }

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
     * \copydoc FvBaseMultiPhaseProblem::intrinsicPermeability
     *
     * In this problem, the upper part of the domain is sightly less
     * permeable than the lower one.
     */
    template <class Context>
    const DimMatrix& intrinsicPermeability(const Context& context, unsigned spaceIdx, unsigned timeIdx) const
    {
        const GlobalPosition& pos = context.pos(spaceIdx, timeIdx);
        if (isFineMaterial_(pos))
            return fineK_;
        return coarseK_;
    }

    /*!
     * \copydoc FvBaseMultiPhaseProblem::porosity
     */
    template <class Context>
    Scalar porosity(const Context& context, unsigned spaceIdx, unsigned timeIdx) const
    {
        const GlobalPosition& pos = context.pos(spaceIdx, timeIdx);
        if (isFineMaterial_(pos))
            return finePorosity_;
        else
            return coarsePorosity_;
    }

    /*!
     * \copydoc FvBaseMultiPhaseProblem::materialLawParams
     */
    template <class Context>
    const MaterialLawParams& materialLawParams(const Context& context,
                                               unsigned spaceIdx,
                                               unsigned timeIdx) const
    {
        const GlobalPosition& pos = context.pos(spaceIdx, timeIdx);
        if (isFineMaterial_(pos))
            return fineMaterialParams_;
        else
            return coarseMaterialParams_;
    }

    /*!
     * \brief Return the parameters for the energy storage law of the rock
     *
     * In this case, we assume the rock-matrix to be granite.
     */
    template <class Context>
    const SolidEnergyLawParams&
    solidEnergyLawParams(const Context& context OPM_UNUSED,
                         unsigned spaceIdx OPM_UNUSED,
                         unsigned timeIdx OPM_UNUSED) const
    { return solidEnergyLawParams_; }

    /*!
     * \copydoc FvBaseMultiPhaseProblem::thermalConductionParams
     */
    template <class Context>
    const ThermalConductionLawParams&
    thermalConductionLawParams(const Context& context,
                            unsigned spaceIdx,
                            unsigned timeIdx) const
    {
        const GlobalPosition& pos = context.pos(spaceIdx, timeIdx);
        if (isFineMaterial_(pos))
            return fineThermalCondParams_;
        return coarseThermalCondParams_;
    }

    //! \}

    /*!
     * \name Boundary conditions
     */
    //! \{

    /*!
     * \copydoc FvBaseProblem::boundary
     *
     * For this problem, we inject air at the inlet on the center of
     * the lower domain boundary and use a no-flow condition on the
     * top boundary and a and a free-flow condition on the left and
     * right boundaries of the domain.
     */
    template <class Context>
    void boundary(BoundaryRateVector& values,
                  const Context& context,
                  unsigned spaceIdx, unsigned timeIdx) const
    {
        const auto& pos = context.cvCenter(spaceIdx, timeIdx);
        assert(onLeftBoundary_(pos) ||
               onLowerBoundary_(pos) ||
               onRightBoundary_(pos) ||
               onUpperBoundary_(pos));

        if (onInlet_(pos)) {
            RateVector massRate(0.0);
            massRate[conti0EqIdx + AirIdx] = -1e-3; // [kg/(m^2 s)]

            // impose an forced inflow boundary condition on the inlet
            values.setMassRate(massRate);

            if (enableEnergy) {
                Opm::CompositionalFluidState<Scalar, FluidSystem> fs;
                initialFluidState_(fs, context, spaceIdx, timeIdx);

                Scalar hl = fs.enthalpy(liquidPhaseIdx);
                Scalar hg = fs.enthalpy(gasPhaseIdx);
                values.setEnthalpyRate(values[conti0EqIdx + AirIdx] * hg +
                                       values[conti0EqIdx + H2OIdx] * hl);
            }
        }
        else if (onLeftBoundary_(pos) || onRightBoundary_(pos)) {
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
     * \name Volumetric terms
     */
    //! \{

    /*!
     * \copydoc FvBaseProblem::initial
     *
     * For this problem, we set the medium to be fully saturated by
     * liquid water and assume hydrostatic pressure.
     */
    template <class Context>
    void initial(PrimaryVariables& values,
                 const Context& context,
                 unsigned spaceIdx,
                 unsigned timeIdx) const
    {
        Opm::CompositionalFluidState<Scalar, FluidSystem> fs;
        initialFluidState_(fs, context, spaceIdx, timeIdx);

        const auto& matParams = materialLawParams(context, spaceIdx, timeIdx);
        values.assignMassConservative(fs, matParams, /*inEquilibrium=*/true);
    }

    /*!
     * \copydoc FvBaseProblem::source
     *
     * For this problem, the source term of all components is 0
     * everywhere.
     */
    template <class Context>
    void source(RateVector& rate,
                const Context& context OPM_UNUSED,
                unsigned spaceIdx OPM_UNUSED,
                unsigned timeIdx OPM_UNUSED) const
    { rate = 0; }

    //! \}

private:
    bool onLeftBoundary_(const GlobalPosition& pos) const
    { return pos[0] < eps_; }

    bool onRightBoundary_(const GlobalPosition& pos) const
    { return pos[0] > this->boundingBoxMax()[0] - eps_; }

    bool onLowerBoundary_(const GlobalPosition& pos) const
    { return pos[1] < eps_; }

    bool onUpperBoundary_(const GlobalPosition& pos) const
    { return pos[1] > this->boundingBoxMax()[1] - eps_; }

    bool onInlet_(const GlobalPosition& pos) const
    { return onLowerBoundary_(pos) && (15.0 < pos[0]) && (pos[0] < 25.0); }

    bool inHighTemperatureRegion_(const GlobalPosition& pos) const
    { return (20 < pos[0]) && (pos[0] < 30) && (pos[1] < 30); }

    template <class Context, class FluidState>
    void initialFluidState_(FluidState& fs,
                            const Context& context,
                            unsigned spaceIdx,
                            unsigned timeIdx) const
    {
        const GlobalPosition& pos = context.pos(spaceIdx, timeIdx);

        Scalar densityW = 1000.0;
        fs.setPressure(liquidPhaseIdx, 1e5 + (maxDepth_ - pos[1])*densityW*9.81);
        fs.setSaturation(liquidPhaseIdx, 1.0);
        fs.setMoleFraction(liquidPhaseIdx, H2OIdx, 1.0);
        fs.setMoleFraction(liquidPhaseIdx, AirIdx, 0.0);

        if (inHighTemperatureRegion_(pos))
            fs.setTemperature(380);
        else
            fs.setTemperature(283.0 + (maxDepth_ - pos[1])*0.03);

        // set the gas saturation and pressure
        fs.setSaturation(gasPhaseIdx, 0);
        Scalar pc[numPhases];
        const auto& matParams = materialLawParams(context, spaceIdx, timeIdx);
        MaterialLaw::capillaryPressures(pc, matParams, fs);
        fs.setPressure(gasPhaseIdx, fs.pressure(liquidPhaseIdx) + (pc[gasPhaseIdx] - pc[liquidPhaseIdx]));

        typename FluidSystem::template ParameterCache<Scalar> paramCache;
        typedef Opm::ComputeFromReferencePhase<Scalar, FluidSystem> CFRP;
        CFRP::solve(fs, paramCache, liquidPhaseIdx, /*setViscosity=*/true,  /*setEnthalpy=*/true);
    }

    void computeThermalCondParams_(ThermalConductionLawParams& params, Scalar poro)
    {
        Scalar lambdaGranite = 2.8; // [W / (K m)]

        // create a Fluid state which has all phases present
        Opm::ImmiscibleFluidState<Scalar, FluidSystem> fs;
        fs.setTemperature(293.15);
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            fs.setPressure(phaseIdx, 1.0135e5);
        }

        typename FluidSystem::template ParameterCache<Scalar> paramCache;
        paramCache.updateAll(fs);
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            Scalar rho = FluidSystem::density(fs, paramCache, phaseIdx);
            fs.setDensity(phaseIdx, rho);
        }

        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
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

    bool isFineMaterial_(const GlobalPosition& pos) const
    { return pos[dim-1] > layerBottom_; }

    DimMatrix fineK_;
    DimMatrix coarseK_;
    Scalar layerBottom_;

    Scalar finePorosity_;
    Scalar coarsePorosity_;

    MaterialLawParams fineMaterialParams_;
    MaterialLawParams coarseMaterialParams_;

    ThermalConductionLawParams fineThermalCondParams_;
    ThermalConductionLawParams coarseThermalCondParams_;
    SolidEnergyLawParams solidEnergyLawParams_;

    Scalar maxDepth_;
    Scalar eps_;
};
} // namespace Opm

#endif
