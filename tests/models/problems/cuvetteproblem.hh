// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2008-2012 by Andreas Lauser                               *
 *   Copyright (C) 2009-2012 by Bernd Flemisch                               *
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
 * \copydoc Ewoms::CuvetteProblem
 */
#ifndef EWOMS_CUVETTE_PROBLEM_HH
#define EWOMS_CUVETTE_PROBLEM_HH

#include <opm/material/fluidstates/CompositionalFluidState.hpp>
#include <opm/material/fluidstates/ImmiscibleFluidState.hpp>
#include <opm/material/fluidsystems/H2OAirMesityleneFluidSystem.hpp>
#include <opm/material/fluidmatrixinteractions/3p/3pParkerVanGenuchten.hpp>
#include <opm/material/fluidmatrixinteractions/3pAdapter.hpp>
#include <opm/material/fluidmatrixinteractions/MpLinearMaterial.hpp>
#include <opm/material/heatconduction/Somerton.hpp>
#include <opm/material/constraintsolvers/MiscibleMultiPhaseComposition.hpp>

#include <ewoms/models/pvs/pvsproperties.hh>

#include <dune/grid/io/file/dgfparser/dgfyasp.hh>

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

#include <string>

namespace Ewoms {
template <class TypeTag>
class CuvetteProblem;
}

namespace Opm {
namespace Properties {

// create a new type tag for the cuvette steam injection problem
NEW_TYPE_TAG(CuvetteBaseProblem);

SET_BOOL_PROP(CuvetteBaseProblem, EnablePartialReassemble, true);
SET_BOOL_PROP(CuvetteBaseProblem, EnableJacobianRecycling, true);

// Set the grid type
SET_TYPE_PROP(CuvetteBaseProblem, Grid, Dune::YaspGrid<2>);

// Set the problem property
SET_TYPE_PROP(CuvetteBaseProblem, Problem, Ewoms::CuvetteProblem<TypeTag>);

// Set the fluid system
SET_TYPE_PROP(CuvetteBaseProblem,
              FluidSystem,
              Opm::FluidSystems::H2OAirMesitylene<typename GET_PROP_TYPE(TypeTag, Scalar)>);

// Enable gravity
SET_BOOL_PROP(CuvetteBaseProblem, EnableGravity, true);

// Set the maximum time step
SET_SCALAR_PROP(CuvetteBaseProblem, MaxTimeStepSize, 600.);

// Set the material Law
SET_PROP(CuvetteBaseProblem, MaterialLaw)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

    enum { wPhaseIdx = FluidSystem::wPhaseIdx };
    enum { nPhaseIdx = FluidSystem::nPhaseIdx };
    enum { gPhaseIdx = FluidSystem::gPhaseIdx };

    // define the three-phase material law
    typedef Opm::ThreePParkerVanGenuchten<Scalar> ThreePLaw;

public:
    // wrap the three-phase law in an adaptor to make use the generic
    // material law API
    typedef Opm::ThreePAdapter<wPhaseIdx, nPhaseIdx, gPhaseIdx, ThreePLaw> type;

    //typedef Opm::MpLinearMaterial<FluidSystem::numPhases, Scalar> type;
};

// Set the heat conduction law
SET_PROP(CuvetteBaseProblem, HeatConductionLaw)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

public:
    // define the material law parameterized by absolute saturations
    typedef Opm::Somerton<FluidSystem, Scalar> type;
};

// The default for the end time of the simulation
SET_SCALAR_PROP(CuvetteBaseProblem, EndTime, 180);

// The default for the initial time step size of the simulation
SET_SCALAR_PROP(CuvetteBaseProblem, InitialTimeStepSize, 1);

// The default DGF file to load
SET_STRING_PROP(CuvetteBaseProblem, GridFile, "./grids/cuvette_11x4.dgf");
} // namespace Properties
} // namespace Opm

namespace Ewoms {
/*!
 * \ingroup VcfvTestProblems
 *
 * \brief Non-isothermal three-phase gas injection problem where a hot gas
 *        is injected into a unsaturated porous medium with a residually
 *        trapped NAPL contamination.
 *
 * The domain is a quasi-two-dimensional container (cuvette). Its
 * dimensions are 1.5 m x 0.74 m. The top and bottom boundaries are
 * closed, the right boundary is a free-flow boundary allowing fluids
 * to escape. From the left, an injection of a hot water-air mixture
 * is injected. The set-up is aimed at remediating an initial NAPL
 * (Non-Aquoeus Phase Liquid) contamination in the domain.  The
 * contamination is initially placed partly into the ambient coarse
 * sand and partly into a fine sand lens.
 *
 * This simulation can be varied through assigning different boundary conditions
 * at the left boundary as described in Class (2001):
 * Theorie und numerische Modellierung nichtisothermer Mehrphasenprozesse in
 * NAPL-kontaminierten poroesen Medien, Dissertation, Eigenverlag des Instituts
 * fuer Wasserbau
 *
 * To see the basic effect and the differences to scenarios with pure
 * steam or pure air injection, it is sufficient to simulate this
 * problem to about 2-3 hours simulation time.  Complete remediation
 * of the domain requires much longer (about 10 days simulated time).
 */
template <class TypeTag >
class CuvetteProblem : public GET_PROP_TYPE(TypeTag, BaseProblem)
{
    typedef typename GET_PROP_TYPE(TypeTag, BaseProblem) ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLawParams) MaterialLawParams;
    typedef typename GET_PROP_TYPE(TypeTag, HeatConductionLaw) HeatConductionLaw;
    typedef typename GET_PROP_TYPE(TypeTag, HeatConductionLawParams) HeatConductionLawParams;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, RateVector) RateVector;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryRateVector) BoundaryRateVector;
    typedef typename GET_PROP_TYPE(TypeTag, TimeManager) TimeManager;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

    // copy some indices for convenience
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    enum {
        numPhases = FluidSystem::numPhases,
        numComponents = FluidSystem::numComponents,

        wPhaseIdx = FluidSystem::wPhaseIdx,
        nPhaseIdx = FluidSystem::nPhaseIdx,
        gPhaseIdx = FluidSystem::gPhaseIdx,

        H2OIdx = FluidSystem::H2OIdx,
        airIdx = FluidSystem::airIdx,
        NAPLIdx = FluidSystem::NAPLIdx,

        conti0EqIdx = Indices::conti0EqIdx,

        // Grid and world dimension
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld
    };

    typedef typename GridView::ctype CoordScalar;
    typedef Dune::FieldVector<CoordScalar, dimWorld> GlobalPosition;
    typedef Dune::FieldMatrix<Scalar, dimWorld, dimWorld> DimMatrix;

public:
    /*!
     * \copydoc Doxygen::defaultProblemConstructor
     */
    CuvetteProblem(TimeManager &timeManager)
        : ParentType(timeManager,
                     GET_PROP_TYPE(TypeTag, GridCreator)::grid().leafView())
        , eps_(1e-6)
    {
        if (Valgrind::IsRunning())
            FluidSystem::init(/*minT=*/283.15, /*maxT=*/500.0, /*nT=*/20,
                              /*minp=*/0.8e5, /*maxp=*/2e5, /*np=*/10);
        else
            FluidSystem::init(/*minT=*/283.15, /*maxT=*/500.0, /*nT=*/200,
                              /*minp=*/0.8e5, /*maxp=*/2e5, /*np=*/100);

        // intrinsic permeabilities
        fineK_ = this->toDimMatrix_(6.28e-12);
        coarseK_ = this->toDimMatrix_(9.14e-10);

        // porosities
        finePorosity_ = 0.42;
        coarsePorosity_ = 0.42;

        // parameters for the capillary pressure law
#if 1
        // three-phase van Genuchten law
        fineMaterialParams_.setVgAlpha(0.0005);
        coarseMaterialParams_.setVgAlpha(0.005);
        fineMaterialParams_.setVgN(4.0);
        coarseMaterialParams_.setVgN(4.0);

        coarseMaterialParams_.setkrRegardsSnr(true);
        fineMaterialParams_.setkrRegardsSnr(true);

        coarseMaterialParams_.setKdNAPL(0.);
        coarseMaterialParams_.setRhoBulk(1500.);
        fineMaterialParams_.setKdNAPL(0.);
        fineMaterialParams_.setRhoBulk(1500.);

        // residual saturations
        fineMaterialParams_.setSwr(0.1201);
        fineMaterialParams_.setSwrx(0.1201);
        fineMaterialParams_.setSnr(0.0701);
        fineMaterialParams_.setSgr(0.0101);
        coarseMaterialParams_.setSwr(0.1201);
        coarseMaterialParams_.setSwrx(0.1201);
        coarseMaterialParams_.setSnr(0.0701);
        coarseMaterialParams_.setSgr(0.0101);
#else
        // linear material law
        fineMaterialParams_.setPcMinSat(gPhaseIdx, 0);
        fineMaterialParams_.setPcMaxSat(gPhaseIdx, 0);
        fineMaterialParams_.setPcMinSat(nPhaseIdx, 0);
        fineMaterialParams_.setPcMaxSat(nPhaseIdx, -1000);
        fineMaterialParams_.setPcMinSat(wPhaseIdx, 0);
        fineMaterialParams_.setPcMaxSat(wPhaseIdx, -10000);

        coarseMaterialParams_.setPcMinSat(gPhaseIdx, 0);
        coarseMaterialParams_.setPcMaxSat(gPhaseIdx, 0);
        coarseMaterialParams_.setPcMinSat(nPhaseIdx, 0);
        coarseMaterialParams_.setPcMaxSat(nPhaseIdx, -100);
        coarseMaterialParams_.setPcMinSat(wPhaseIdx, 0);
        coarseMaterialParams_.setPcMaxSat(wPhaseIdx, -1000);

        // residual saturations
        fineMaterialParams_.setResidSat(wPhaseIdx, 0.1201);
        fineMaterialParams_.setResidSat(nPhaseIdx, 0.0701);
        fineMaterialParams_.setResidSat(gPhaseIdx, 0.0101);

        coarseMaterialParams_.setResidSat(wPhaseIdx, 0.1201);
        coarseMaterialParams_.setResidSat(nPhaseIdx, 0.0701);
        coarseMaterialParams_.setResidSat(gPhaseIdx, 0.0101);
#endif

        // initialize parameters for the heat conduction law
        computeHeatCondParams_(heatCondParams_, finePorosity_);

        initInjectFluidState_();
    }

    /*!
     * \name Auxiliary methods
     */
    //! \{

    /*!
     * \copydoc VcfvProblem::shouldWriteRestartFile
     *
     * This problem writes a restart file after every time step.
     */
    bool shouldWriteRestartFile() const
    { return true; }

    /*!
     * \copydoc VcfvProblem::name
     */
    const char *name() const
    {
        static std::string tmp = std::string("cuvette_")+this->model().name();
        return tmp.c_str();
    }

    //! \}

    /*!
     * \name Soil parameters
     */
    //! \{

    /*!
     * \copydoc VcfvMultiPhaseProblem::temperature
     */
    template <class Context>
    Scalar temperature(const Context &context, int spaceIdx, int timeIdx) const
    { return 293.15; /* [K] */ }

    /*!
     * \copydoc VcfvMultiPhaseProblem::intrinsicPermeability
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
     * \copydoc VcfvMultiPhaseProblem::heatConductionParams
     */
    template <class Context>
    const HeatConductionLawParams&
    heatConductionParams(const Context &context, int spaceIdx, int timeIdx) const
    { return heatCondParams_; }

    /*!
     * \copydoc VcfvMultiPhaseProblem::heatCapacitySolid
     */
    template <class Context>
    Scalar heatCapacitySolid(const Context &context, int spaceIdx, int timeIdx) const
    {
        return
            850 // specific heat capacity [J / (kg K)]
            * 2650; // density of sand [kg/m^3]
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
    void boundary(BoundaryRateVector &values, const Context &context, int spaceIdx, int timeIdx) const
    {
        const auto &pos = context.pos(spaceIdx, timeIdx);

        if (onRightBoundary_(pos)) {
            Opm::CompositionalFluidState<Scalar, FluidSystem> fs;

            initialFluidState_(fs, context, spaceIdx, timeIdx);

            values.setFreeFlow(context, spaceIdx, timeIdx, fs);
            values.setNoFlow();
        }
        else  if (onLeftBoundary_(pos))
        {
            // injection
            RateVector molarRate;

            // inject with the same composition as the gas phase of
            // the injection fluid state
            Scalar molarInjectionRate = 0.3435; // [mol/(m^2 s)]
            for (int compIdx = 0; compIdx < numComponents; ++ compIdx)
                molarRate[conti0EqIdx + compIdx] =
                    - molarInjectionRate
                    * injectFluidState_.moleFraction(gPhaseIdx, compIdx);

            // calculate the total mass injection rate [kg / (m^2 s)
            Scalar massInjectionRate =
                molarInjectionRate
                * injectFluidState_.averageMolarMass(gPhaseIdx);

            // set the boundary rate vector
            values.setMolarRate(molarRate);
            values.setEnthalpyRate(- injectFluidState_.enthalpy(gPhaseIdx)
                                   * massInjectionRate); // [J / (m^2 s)]
        }
        else
            values.setNoFlow();
    }

    //! \}

    /*!
     * \name Volume terms
     */
    //! \{

    /*!
     * \copydoc VcfvProblem::initial
     */
    template <class Context>
    void initial(PrimaryVariables &values, const Context &context, int spaceIdx, int timeIdx) const
    {
        Opm::CompositionalFluidState<Scalar, FluidSystem> fs;

        initialFluidState_(fs, context, spaceIdx, timeIdx);

        const auto &matParams = materialLawParams(context, spaceIdx, timeIdx);
        values.assignMassConservative(fs, matParams, /*inEquilibrium=*/false);
    }

    /*!
     * \copydoc VcfvProblem::source
     *
     * For this problem, the source term of all components is 0
     * everywhere.
     */
    template <class Context>
    void source(RateVector &rate, const Context &context, int spaceIdx, int timeIdx) const
    { rate = Scalar(0.0); }

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

    bool isContaminated_(const GlobalPosition &pos) const
    {
        return
            (0.20 <= pos[0]) && (pos[0] <= 0.80)
            && (0.4 <= pos[1]) && (pos[1] <= 0.65);
    }

    bool isFineMaterial_(const GlobalPosition &pos) const
    {
        if (0.13 <= pos[0] && 1.20 >= pos[0] && 0.32 <= pos[1] && pos[1] <= 0.57)
            return true;
        else if (pos[1] <= 0.15 && 1.20 <= pos[0])
            return true;
        else return false;
    }

    template <class FluidState, class Context>
    void initialFluidState_(FluidState &fs,
                            const Context &context,
                            int spaceIdx,
                            int timeIdx) const
    {
        const GlobalPosition &pos = context.pos(spaceIdx, timeIdx);

        fs.setTemperature(293.0 /*[K]*/);

        Scalar pw = 1e5;

        if(isContaminated_(pos)) {
            fs.setSaturation(wPhaseIdx, 0.12);
            fs.setSaturation(nPhaseIdx, 0.07);
            fs.setSaturation(gPhaseIdx, 1 - 0.12 - 0.07);

            // set the capillary pressures
            const auto &matParams = materialLawParams(context, spaceIdx, timeIdx);
            Scalar pc[numPhases];
            MaterialLaw::capillaryPressures(pc, matParams, fs);
            for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
                fs.setPressure(phaseIdx,
                               pw + (pc[phaseIdx] - pc[wPhaseIdx]));

            // compute the phase compositions
            typedef Opm::MiscibleMultiPhaseComposition<Scalar, FluidSystem> MMPC;
            typename FluidSystem::ParameterCache paramCache;
            MMPC::solve(fs, paramCache, /*setViscosity=*/true, /*setEnthalpy=*/true);
        }
        else {
            fs.setSaturation(wPhaseIdx, 0.12);
            fs.setSaturation(gPhaseIdx, 1 - fs.saturation(wPhaseIdx));
            fs.setSaturation(nPhaseIdx, 0);

            // set the capillary pressures
            const auto &matParams = materialLawParams(context, spaceIdx, timeIdx);
            Scalar pc[numPhases];
            MaterialLaw::capillaryPressures(pc, matParams, fs);
            for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
                fs.setPressure(phaseIdx,
                               pw + (pc[phaseIdx] - pc[wPhaseIdx]));

            // compute the phase compositions
            typedef Opm::MiscibleMultiPhaseComposition<Scalar, FluidSystem> MMPC;
            typename FluidSystem::ParameterCache paramCache;
            MMPC::solve(fs, paramCache, /*setViscosity=*/true, /*setEnthalpy=*/true);

            // set the contaminant mole fractions to zero. this is a
            // little bit hacky...
            for (int phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx) {
                fs.setMoleFraction(phaseIdx, NAPLIdx, 0.0);

                if (phaseIdx == nPhaseIdx)
                    continue;

                Scalar sumx = 0;
                for (int compIdx = 0; compIdx < numComponents; ++ compIdx)
                    sumx += fs.moleFraction(phaseIdx, compIdx);

                for (int compIdx = 0; compIdx < numComponents; ++ compIdx)
                    fs.setMoleFraction(phaseIdx,
                                       compIdx,
                                       fs.moleFraction(phaseIdx, compIdx) / sumx);
            }
        }
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

    void initInjectFluidState_()
    {
        injectFluidState_.setTemperature(383.0); // [K]
        injectFluidState_.setPressure(gPhaseIdx, 1e5); // [Pa]
        injectFluidState_.setSaturation(gPhaseIdx, 1.0); // [-]

        Scalar xgH2O = 0.417;
        injectFluidState_.setMoleFraction(gPhaseIdx,
                                          H2OIdx,
                                          xgH2O); // [-]
        injectFluidState_.setMoleFraction(gPhaseIdx,
                                          airIdx,
                                          1 - xgH2O); // [-]
        injectFluidState_.setMoleFraction(gPhaseIdx,
                                          NAPLIdx,
                                          0.0); // [-]

        // set the specific enthalpy of the gas phase
        typename FluidSystem::ParameterCache paramCache;
        paramCache.updatePhase(injectFluidState_, gPhaseIdx);

        Scalar h = FluidSystem::enthalpy(injectFluidState_, paramCache, gPhaseIdx);
        injectFluidState_.setEnthalpy(gPhaseIdx, h);
    }

    DimMatrix fineK_;
    DimMatrix coarseK_;

    Scalar finePorosity_;
    Scalar coarsePorosity_;

    MaterialLawParams fineMaterialParams_;
    MaterialLawParams coarseMaterialParams_;

    HeatConductionLawParams heatCondParams_;

    Opm::CompositionalFluidState<Scalar, FluidSystem> injectFluidState_;

    const Scalar eps_;
};
} // namespace Ewoms

#endif
