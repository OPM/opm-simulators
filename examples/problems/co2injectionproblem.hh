/*
  Copyright (C) 2008-2013 by Andreas Lauser

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
 * \copydoc Ewoms::Co2InjectionProblem
 */
#ifndef EWOMS_CO2_INJECTION_PROBLEM_HH
#define EWOMS_CO2_INJECTION_PROBLEM_HH

#include <ewoms/models/immiscible/immisciblemodel.hh>

#include <opm/material/fluidsystems/H2ON2FluidSystem.hpp>
#include <opm/material/fluidsystems/BrineCO2FluidSystem.hpp>
#include <opm/material/fluidstates/CompositionalFluidState.hpp>
#include <opm/material/fluidstates/ImmiscibleFluidState.hpp>
#include <opm/material/constraintsolvers/ComputeFromReferencePhase.hpp>
#include <opm/material/fluidmatrixinteractions/LinearMaterial.hpp>
#include <opm/material/fluidmatrixinteractions/RegularizedBrooksCorey.hpp>
#include <opm/material/fluidmatrixinteractions/EffToAbsLaw.hpp>
#include <opm/material/fluidmatrixinteractions/MaterialTraits.hpp>
#include <opm/material/heatconduction/Somerton.hpp>
#include <opm/material/binarycoefficients/Brine_CO2.hpp>
#include <opm/material/StaticTabulated2dFunction.hpp>

#include <dune/grid/yaspgrid.hh>

#include <dune/common/version.hh>
#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

#include <sstream>
#include <iostream>
#include <string>

namespace Ewoms {
template <class TypeTag>
class Co2InjectionProblem;

namespace Co2Injection {
#include <opm/material/components/co2tables.inc>
}
} // namespace Ewoms

namespace Opm {
namespace Properties {
NEW_TYPE_TAG(Co2InjectionBaseProblem);

// declare the CO2 injection problem specific property tags
NEW_PROP_TAG(FluidSystemPressureLow);
NEW_PROP_TAG(FluidSystemPressureHigh);
NEW_PROP_TAG(FluidSystemNumPressure);
NEW_PROP_TAG(FluidSystemTemperatureLow);
NEW_PROP_TAG(FluidSystemTemperatureHigh);
NEW_PROP_TAG(FluidSystemNumTemperature);

NEW_PROP_TAG(MaxDepth);
NEW_PROP_TAG(Temperature);
NEW_PROP_TAG(SimulationName);

// Set the grid type
SET_TYPE_PROP(Co2InjectionBaseProblem, Grid, Dune::YaspGrid<2>);

// Set the problem property
SET_TYPE_PROP(Co2InjectionBaseProblem, Problem,
              Ewoms::Co2InjectionProblem<TypeTag>);

// Set fluid configuration
SET_PROP(Co2InjectionBaseProblem, FluidSystem)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef Ewoms::Co2Injection::CO2Tables CO2Tables;

    static const bool useComplexRelations = false;

public:
    typedef Opm::FluidSystems::BrineCO2<Scalar, CO2Tables> type;
    // typedef Opm::FluidSystems::H2ON2<Scalar, useComplexRelations> type;
};

// Set the material Law
SET_PROP(Co2InjectionBaseProblem, MaterialLaw)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    enum { lPhaseIdx = FluidSystem::lPhaseIdx };
    enum { gPhaseIdx = FluidSystem::gPhaseIdx };

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef Opm::TwoPhaseMaterialTraits<Scalar,
                                        /*wettingPhaseIdx=*/FluidSystem::lPhaseIdx,
                                        /*nonWettingPhaseIdx=*/FluidSystem::gPhaseIdx>
    Traits;

    // define the material law which is parameterized by effective
    // saturations
    typedef Opm::RegularizedBrooksCorey<Traits> EffMaterialLaw;

public:
    // define the material law parameterized by absolute saturations
    typedef Opm::EffToAbsLaw<EffMaterialLaw> type;
};

// Set the heat conduction law
SET_PROP(Co2InjectionBaseProblem, HeatConductionLaw)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

public:
    // define the material law parameterized by absolute saturations
    typedef Opm::Somerton<FluidSystem, Scalar> type;
};

// Write the Newton convergence behavior to disk?
SET_BOOL_PROP(Co2InjectionBaseProblem, NewtonWriteConvergence, false);

// Enable gravity
SET_BOOL_PROP(Co2InjectionBaseProblem, EnableGravity, true);

// Reuse Jacobian matrices if possible?
SET_BOOL_PROP(Co2InjectionBaseProblem, EnableJacobianRecycling, true);

// Smoothen the upwinding method?
SET_BOOL_PROP(Co2InjectionBaseProblem, EnableSmoothUpwinding, false);

// set the defaults for the problem specific properties
SET_SCALAR_PROP(Co2InjectionBaseProblem, FluidSystemPressureLow, 3e7);
SET_SCALAR_PROP(Co2InjectionBaseProblem, FluidSystemPressureHigh, 4e7);
SET_INT_PROP(Co2InjectionBaseProblem, FluidSystemNumPressure, 100);
SET_SCALAR_PROP(Co2InjectionBaseProblem, FluidSystemTemperatureLow, 290);
SET_SCALAR_PROP(Co2InjectionBaseProblem, FluidSystemTemperatureHigh, 500);
SET_INT_PROP(Co2InjectionBaseProblem, FluidSystemNumTemperature, 100);

SET_SCALAR_PROP(Co2InjectionBaseProblem, MaxDepth, 2500);
SET_SCALAR_PROP(Co2InjectionBaseProblem, Temperature, 293.15);
SET_STRING_PROP(Co2InjectionBaseProblem, SimulationName, "co2injection");

// The default for the end time of the simulation
SET_SCALAR_PROP(Co2InjectionBaseProblem, EndTime, 1e4);

// The default for the initial time step size of the simulation
SET_SCALAR_PROP(Co2InjectionBaseProblem, InitialTimeStepSize, 250);

// The default DGF file to load
SET_STRING_PROP(Co2InjectionBaseProblem, GridFile, "grids/co2injection.dgf");
} // namespace Properties
} // namespace Opm

namespace Ewoms {
/*!
 * \ingroup VcfvTestProblems
 *
 * \brief Problem where \f$CO_2\f$ is injected under a low permeable
 *        layer at a depth of 2700m.
 *
 * The domain is sized 60m times 40m and consists of two layers, one
 * which is moderately permeable (\f$K = 10^{-12}\;m^2\f$) for \f$ y >
 * 22\; m\f$ and one with a lower intrinsic permeablility (\f$
 * K=10^{-13}\;m^2\f$) in the rest of the domain.
 *
 * \f$CO_2\f$ gets injected by means of a forced-flow boundary
 * condition into water-filled aquifer, which is situated 2700m below
 * sea level, at the lower-right boundary (\f$5m<y<15m\f$) and
 * migrates upwards due to buoyancy. It accumulates and eventually
 * enters the lower permeable aquitard.
 *
 * The boundary conditions applied by this problem are no-flow
 * conditions on the top bottom and right boundaries and a free-flow
 * boundary condition on the left. For the free-flow condition,
 * hydrostatic pressure is assumed.
 */
template <class TypeTag>
class Co2InjectionProblem : public GET_PROP_TYPE(TypeTag, BaseProblem)
{
    typedef typename GET_PROP_TYPE(TypeTag, BaseProblem) ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

    enum { dim = GridView::dimension };
    enum { dimWorld = GridView::dimensionworld };

    // copy some indices for convenience
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    enum { numPhases = FluidSystem::numPhases };
    enum { gPhaseIdx = FluidSystem::gPhaseIdx };
    enum { lPhaseIdx = FluidSystem::lPhaseIdx };
    enum { CO2Idx = FluidSystem::CO2Idx };
    enum { BrineIdx = FluidSystem::BrineIdx };
    enum { conti0EqIdx = Indices::conti0EqIdx };
    enum { contiCO2EqIdx = conti0EqIdx + CO2Idx };

    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, RateVector) RateVector;
    typedef typename GET_PROP_TYPE(TypeTag,
                                   BoundaryRateVector) BoundaryRateVector;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;
    typedef typename GET_PROP_TYPE(TypeTag, TimeManager) TimeManager;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLawParams) MaterialLawParams;
    typedef typename GET_PROP_TYPE(TypeTag, HeatConductionLaw) HeatConductionLaw;
    typedef typename HeatConductionLaw::Params HeatConductionLawParams;

    typedef typename GridView::ctype CoordScalar;
    typedef Dune::FieldVector<CoordScalar, dimWorld> GlobalPosition;
    typedef Dune::FieldMatrix<Scalar, dimWorld, dimWorld> DimMatrix;

public:
    /*!
     * \copydoc Doxygen::defaultProblemConstructor
     */
    Co2InjectionProblem(TimeManager &timeManager)
#if DUNE_VERSION_NEWER(DUNE_COMMON, 2, 3)
        : ParentType(timeManager,
                     GET_PROP_TYPE(TypeTag, GridCreator)::grid().leafGridView())
#else
        : ParentType(timeManager,
                     GET_PROP_TYPE(TypeTag, GridCreator)::grid().leafView())
#endif
    {
        eps_ = 1e-6;

        temperatureLow_ = EWOMS_GET_PARAM(TypeTag, Scalar, FluidSystemTemperatureLow);
        temperatureHigh_ = EWOMS_GET_PARAM(TypeTag, Scalar, FluidSystemTemperatureHigh);
        nTemperature_ = EWOMS_GET_PARAM(TypeTag, int, FluidSystemNumTemperature);

        nPressure_ = EWOMS_GET_PARAM(TypeTag, int, FluidSystemNumPressure);
        pressureLow_ = EWOMS_GET_PARAM(TypeTag, Scalar, FluidSystemPressureLow);
        pressureHigh_ = EWOMS_GET_PARAM(TypeTag, Scalar, FluidSystemPressureHigh);

        maxDepth_ = EWOMS_GET_PARAM(TypeTag, Scalar, MaxDepth);
        temperature_ = EWOMS_GET_PARAM(TypeTag, Scalar, Temperature);
        name_ = EWOMS_GET_PARAM(TypeTag, std::string, SimulationName);

        // initialize the tables of the fluid system
        // FluidSystem::init();
        FluidSystem::init(/*Tmin=*/temperatureLow_,
                          /*Tmax=*/temperatureHigh_,
                          /*nT=*/nTemperature_,
                          /*pmin=*/pressureLow_,
                          /*pmax=*/pressureHigh_,
                          /*np=*/nPressure_);

        fineLayerBottom_ = 22.0;

        // intrinsic permeabilities
        fineK_ = this->toDimMatrix_(1e-13);
        coarseK_ = this->toDimMatrix_(1e-12);

        // porosities
        finePorosity_ = 0.3;
        coarsePorosity_ = 0.3;

        // residual saturations
        fineMaterialParams_.setResidualSaturation(lPhaseIdx, 0.2);
        fineMaterialParams_.setResidualSaturation(gPhaseIdx, 0.0);
        coarseMaterialParams_.setResidualSaturation(lPhaseIdx, 0.2);
        coarseMaterialParams_.setResidualSaturation(gPhaseIdx, 0.0);

        // parameters for the Brooks-Corey law
        fineMaterialParams_.setEntryPressure(1e4);
        coarseMaterialParams_.setEntryPressure(5e3);
        fineMaterialParams_.setLambda(2.0);
        coarseMaterialParams_.setLambda(2.0);

        fineMaterialParams_.finalize();
        coarseMaterialParams_.finalize();

        // parameters for the somerton law of heat conduction
        computeHeatCondParams_(fineHeatCondParams_, finePorosity_);
        computeHeatCondParams_(coarseHeatCondParams_, coarsePorosity_);
    }

    /*!
     * \copydoc FvBaseMultiPhaseProblem::registerParameters
     */
    static void registerParameters()
    {
        ParentType::registerParameters();

        EWOMS_REGISTER_PARAM(TypeTag, Scalar, FluidSystemTemperatureLow,
                             "The lower temperature [K] for tabulation of the "
                             "fluid system");
        EWOMS_REGISTER_PARAM(TypeTag, Scalar, FluidSystemTemperatureHigh,
                             "The upper temperature [K] for tabulation of the "
                             "fluid system");
        EWOMS_REGISTER_PARAM(TypeTag, int, FluidSystemNumTemperature,
                             "The number of intervals between the lower and "
                             "upper temperature");

        EWOMS_REGISTER_PARAM(TypeTag, Scalar, FluidSystemPressureLow,
                             "The lower pressure [Pa] for tabulation of the "
                             "fluid system");
        EWOMS_REGISTER_PARAM(TypeTag, Scalar, FluidSystemPressureHigh,
                             "The upper pressure [Pa] for tabulation of the "
                             "fluid system");
        EWOMS_REGISTER_PARAM(TypeTag, int, FluidSystemNumPressure,
                             "The number of intervals between the lower and "
                             "upper pressure");

        EWOMS_REGISTER_PARAM(TypeTag, Scalar, Temperature,
                             "The temperature [K] in the reservoir");
        EWOMS_REGISTER_PARAM(TypeTag, Scalar, MaxDepth,
                             "The maximum depth [m] of the reservoir");
        EWOMS_REGISTER_PARAM(TypeTag, std::string, SimulationName,
                             "The name of the simulation used for the output "
                             "files");
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
        oss << name_ << "_" << this->model().name();
        if (GET_PROP_VALUE(TypeTag, EnableEnergy))
            oss << "_ni";
        oss << "_" << this->model().discretizationName();
        return oss.str();
    }

    /*!
     * \copydoc VcfvProblem::postTimeStep
     */
    void postTimeStep()
    {
        // Calculate storage terms
        PrimaryVariables storageL, storageG;
        this->model().globalPhaseStorage(storageL, /*phaseIdx=*/0);
        this->model().globalPhaseStorage(storageG, /*phaseIdx=*/1);

        // Write mass balance information for rank 0
        if (this->gridView().comm().rank() == 0) {
            std::cout << "Storage: liquid=[" << storageL << "]"
                      << " gas=[" << storageG << "]\n" << std::flush;
        }
    }

    /*!
     * \copydoc FvBaseMultiPhaseProblem::temperature
     */
    template <class Context>
    Scalar temperature(const Context &context, int spaceIdx, int timeIdx) const
    {
        const auto &pos = context.pos(spaceIdx, timeIdx);
        if (inHighTemperatureRegion_(pos))
            return temperature_ + 100;
        return temperature_;
    }

    /*!
     * \copydoc FvBaseMultiPhaseProblem::intrinsicPermeability
     */
    template <class Context>
    const DimMatrix &intrinsicPermeability(const Context &context, int spaceIdx,
                                           int timeIdx) const
    {
        const GlobalPosition &pos = context.pos(spaceIdx, timeIdx);
        if (isFineMaterial_(pos))
            return fineK_;
        return coarseK_;
    }

    /*!
     * \copydoc FvBaseMultiPhaseProblem::porosity
     */
    template <class Context>
    Scalar porosity(const Context &context, int spaceIdx, int timeIdx) const
    {
        const GlobalPosition &pos = context.pos(spaceIdx, timeIdx);
        if (isFineMaterial_(pos))
            return finePorosity_;
        return coarsePorosity_;
    }

    /*!
     * \copydoc FvBaseMultiPhaseProblem::materialLawParams
     */
    template <class Context>
    const MaterialLawParams &materialLawParams(const Context &context,
                                               int spaceIdx, int timeIdx) const
    {
        const GlobalPosition &pos = context.pos(spaceIdx, timeIdx);
        if (isFineMaterial_(pos))
            return fineMaterialParams_;
        return coarseMaterialParams_;
    }

    /*!
     * \copydoc FvBaseMultiPhaseProblem::heatCapacitySolid
     *
     * In this case, we assume the rock-matrix to be granite.
     */
    template <class Context>
    Scalar heatCapacitySolid(const Context &context, int spaceIdx,
                             int timeIdx) const
    {
        return 790     // specific heat capacity of granite [J / (kg K)]
               * 2700; // density of granite [kg/m^3]
    }

    /*!
     * \copydoc FvBaseMultiPhaseProblem::heatConductionParams
     */
    template <class Context>
    const HeatConductionLawParams &
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
     */
    template <class Context>
    void boundary(BoundaryRateVector &values, const Context &context,
                  int spaceIdx, int timeIdx) const
    {
        const auto &pos = context.pos(spaceIdx, timeIdx);

        if (onLeftBoundary_(pos)) {
            Opm::CompositionalFluidState<Scalar, FluidSystem> fs;
            initialFluidState_(fs, context, spaceIdx, timeIdx);
            fs.checkDefined();

            // impose an freeflow boundary condition
            values.setFreeFlow(context, spaceIdx, timeIdx, fs);
        }
        else if (onInlet_(pos)) {
            RateVector massRate(0.0);
            massRate[contiCO2EqIdx] = -1e-3; // [kg/(m^3 s)]

            Opm::ImmiscibleFluidState<Scalar, FluidSystem> fs;
            fs.setSaturation(gPhaseIdx, 1.0);
            fs.setPressure(gPhaseIdx,
                           context.volVars(spaceIdx, timeIdx).fluidState().pressure(
                               gPhaseIdx));
            fs.setTemperature(temperature(context, spaceIdx, timeIdx));
            typename FluidSystem::ParameterCache paramCache;
            paramCache.updatePhase(fs, gPhaseIdx);
            Scalar h = FluidSystem::enthalpy(fs, paramCache, gPhaseIdx);

            // impose an forced inflow boundary condition
            values.setMassRate(massRate);
            values.setEnthalpyRate(massRate[contiCO2EqIdx] * h);
        }
        else
            // no flow on top and bottom
            values.setNoFlow();
    }

    // \}

    /*!
     * \name Volume terms
     */
    //! \{

    /*!
     * \copydoc VcfvProblem::initial
     */
    template <class Context>
    void initial(PrimaryVariables &values, const Context &context, int spaceIdx,
                 int timeIdx) const
    {
        Opm::CompositionalFluidState<Scalar, FluidSystem> fs;
        initialFluidState_(fs, context, spaceIdx, timeIdx);

        // const auto &matParams = this->materialLawParams(context, spaceIdx,
        // timeIdx);
        // values.assignMassConservative(fs, matParams, /*inEquilibrium=*/true);
        values.assignNaive(fs);
    }

    /*!
     * \copydoc VcfvProblem::source
     *
     * For this problem, the source term of all components is 0
     * everywhere.
     */
    template <class Context>
    void source(RateVector &rate, const Context &context, int spaceIdx,
                int timeIdx) const
    { rate = Scalar(0.0); }

    //! \}

private:
    template <class Context, class FluidState>
    void initialFluidState_(FluidState &fs, const Context &context,
                            int spaceIdx, int timeIdx) const
    {
        const GlobalPosition &pos = context.pos(spaceIdx, timeIdx);

        //////
        // set temperature
        //////
        fs.setTemperature(temperature(context, spaceIdx, timeIdx));

        //////
        // set saturations
        //////
        fs.setSaturation(FluidSystem::lPhaseIdx, 1.0);
        fs.setSaturation(FluidSystem::gPhaseIdx, 0.0);

        //////
        // set pressures
        //////
        Scalar densityL = FluidSystem::Brine::liquidDensity(temperature_, 1e5);
        Scalar depth = maxDepth_ - pos[dim - 1];
        Scalar pl = 1e5 - densityL * this->gravity()[dim - 1] * depth;

        Scalar pC[numPhases];
        const auto &matParams
            = this->materialLawParams(context, spaceIdx, timeIdx);
        MaterialLaw::capillaryPressures(pC, matParams, fs);

        fs.setPressure(lPhaseIdx, pl + (pC[lPhaseIdx] - pC[lPhaseIdx]));
        fs.setPressure(gPhaseIdx, pl + (pC[gPhaseIdx] - pC[lPhaseIdx]));

        //////
        // set composition of the liquid phase
        //////
        fs.setMoleFraction(lPhaseIdx, CO2Idx, 0.005);
        fs.setMoleFraction(lPhaseIdx, BrineIdx,
                           1.0 - fs.moleFraction(lPhaseIdx, CO2Idx));

        typename FluidSystem::ParameterCache paramCache;
        typedef Opm::ComputeFromReferencePhase<Scalar, FluidSystem> CFRP;
        CFRP::solve(fs, paramCache,
                    /*refPhaseIdx=*/lPhaseIdx,
                    /*setViscosity=*/true,
                    /*setEnthalpy=*/true);
    }

    bool onLeftBoundary_(const GlobalPosition &pos) const
    { return pos[0] < eps_; }

    bool onRightBoundary_(const GlobalPosition &pos) const
    { return pos[0] > this->boundingBoxMax()[0] - eps_; }

    bool onInlet_(const GlobalPosition &pos) const
    { return onRightBoundary_(pos) && (5 < pos[1]) && (pos[1] < 15); }

    bool inHighTemperatureRegion_(const GlobalPosition &pos) const
    { return (pos[0] > 20) && (pos[0] < 30) && (pos[1] > 5) && (pos[1] < 35); }

    void computeHeatCondParams_(HeatConductionLawParams &params, Scalar poro)
    {
        Scalar lambdaWater = 0.6;
        Scalar lambdaGranite = 2.8;

        Scalar lambdaWet = std::pow(lambdaGranite, (1 - poro))
                           * std::pow(lambdaWater, poro);
        Scalar lambdaDry = std::pow(lambdaGranite, (1 - poro));

        params.setFullySaturatedLambda(gPhaseIdx, lambdaDry);
        params.setFullySaturatedLambda(lPhaseIdx, lambdaWet);
        params.setVacuumLambda(lambdaDry);
    }

    bool isFineMaterial_(const GlobalPosition &pos) const
    { return pos[dim - 1] > fineLayerBottom_; }

    DimMatrix fineK_;
    DimMatrix coarseK_;
    Scalar fineLayerBottom_;

    Scalar finePorosity_;
    Scalar coarsePorosity_;

    MaterialLawParams fineMaterialParams_;
    MaterialLawParams coarseMaterialParams_;

    HeatConductionLawParams fineHeatCondParams_;
    HeatConductionLawParams coarseHeatCondParams_;

    Scalar temperature_;
    Scalar maxDepth_;
    Scalar eps_;

    int nTemperature_;
    int nPressure_;

    std::string name_;

    Scalar pressureLow_, pressureHigh_;
    Scalar temperatureLow_, temperatureHigh_;
};
} // namespace Ewoms

#endif
