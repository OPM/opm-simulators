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
 * \copydoc Opm::CuvetteProblem
 */
#ifndef EWOMS_CUVETTE_PROBLEM_HH
#define EWOMS_CUVETTE_PROBLEM_HH

#include <ewoms/models/pvs/pvsproperties.hh>

#include <opm/material/fluidstates/CompositionalFluidState.hpp>
#include <opm/material/fluidstates/ImmiscibleFluidState.hpp>
#include <opm/material/fluidsystems/H2OAirMesityleneFluidSystem.hpp>
#include <opm/material/fluidmatrixinteractions/ThreePhaseParkerVanGenuchten.hpp>
#include <opm/material/fluidmatrixinteractions/LinearMaterial.hpp>
#include <opm/material/thermal/ConstantSolidHeatCapLaw.hpp>
#include <opm/material/thermal/SomertonThermalConductionLaw.hpp>
#include <opm/material/constraintsolvers/MiscibleMultiPhaseComposition.hpp>
#include <opm/material/fluidmatrixinteractions/MaterialTraits.hpp>
#include <opm/material/common/Valgrind.hpp>
#include <opm/material/common/Unused.hpp>

#include <dune/grid/yaspgrid.hh>
#include <dune/grid/io/file/dgfparser/dgfyasp.hh>

#include <dune/common/version.hh>
#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

#include <string>

namespace Opm {
template <class TypeTag>
class CuvetteProblem;
}

BEGIN_PROPERTIES


// create a new type tag for the cuvette steam injection problem
NEW_TYPE_TAG(CuvetteBaseProblem);

// Set the grid type
SET_TYPE_PROP(CuvetteBaseProblem, Grid, Dune::YaspGrid<2>);

// Set the problem property
SET_TYPE_PROP(CuvetteBaseProblem, Problem, Opm::CuvetteProblem<TypeTag>);

// Set the fluid system
SET_TYPE_PROP(
    CuvetteBaseProblem, FluidSystem,
    Opm::H2OAirMesityleneFluidSystem<typename GET_PROP_TYPE(TypeTag, Scalar)>);

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

    typedef Opm::ThreePhaseMaterialTraits<
        Scalar,
        /*wettingPhaseIdx=*/FluidSystem::waterPhaseIdx,
        /*nonWettingPhaseIdx=*/FluidSystem::naplPhaseIdx,
        /*gasPhaseIdx=*/FluidSystem::gasPhaseIdx> Traits;

public:
    typedef Opm::ThreePhaseParkerVanGenuchten<Traits> type;
};

// set the energy storage law for the solid phase
SET_TYPE_PROP(CuvetteBaseProblem, SolidEnergyLaw,
              Opm::ConstantSolidHeatCapLaw<typename GET_PROP_TYPE(TypeTag, Scalar)>);

// Set the thermal conduction law
SET_PROP(CuvetteBaseProblem, ThermalConductionLaw)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

public:
    // define the material law parameterized by absolute saturations
    typedef Opm::SomertonThermalConductionLaw<FluidSystem, Scalar> type;
};

// The default for the end time of the simulation
SET_SCALAR_PROP(CuvetteBaseProblem, EndTime, 180);

// The default for the initial time step size of the simulation
SET_SCALAR_PROP(CuvetteBaseProblem, InitialTimeStepSize, 1);

// The default DGF file to load
SET_STRING_PROP(CuvetteBaseProblem, GridFile, "./data/cuvette_11x4.dgf");

END_PROPERTIES

namespace Opm {
/*!
 * \ingroup TestProblems
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
template <class TypeTag>
class CuvetteProblem : public GET_PROP_TYPE(TypeTag, BaseProblem)
{
    typedef typename GET_PROP_TYPE(TypeTag, BaseProblem) ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLawParams) MaterialLawParams;
    typedef typename GET_PROP_TYPE(TypeTag, ThermalConductionLawParams) ThermalConductionLawParams;
    typedef typename GET_PROP_TYPE(TypeTag, SolidEnergyLawParams) SolidEnergyLawParams;
    typedef typename GET_PROP_TYPE(TypeTag, EqVector) EqVector;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, RateVector) RateVector;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryRateVector) BoundaryRateVector;
    typedef typename GET_PROP_TYPE(TypeTag, Simulator) Simulator;
    typedef typename GET_PROP_TYPE(TypeTag, Model) Model;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

    // copy some indices for convenience
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    enum { numPhases = FluidSystem::numPhases };
    enum { numComponents = FluidSystem::numComponents };
    enum { waterPhaseIdx = FluidSystem::waterPhaseIdx };
    enum { naplPhaseIdx = FluidSystem::naplPhaseIdx };
    enum { gasPhaseIdx = FluidSystem::gasPhaseIdx };
    enum { H2OIdx = FluidSystem::H2OIdx };
    enum { airIdx = FluidSystem::airIdx };
    enum { NAPLIdx = FluidSystem::NAPLIdx };
    enum { conti0EqIdx = Indices::conti0EqIdx };

    // Grid and world dimension
    enum { dimWorld = GridView::dimensionworld };

    typedef typename GridView::ctype CoordScalar;
    typedef Dune::FieldVector<CoordScalar, dimWorld> GlobalPosition;
    typedef Dune::FieldMatrix<Scalar, dimWorld, dimWorld> DimMatrix;

public:
    /*!
     * \copydoc Doxygen::defaultProblemConstructor
     */
    CuvetteProblem(Simulator& simulator)
        : ParentType(simulator)
        , eps_(1e-6)
    { }

    /*!
     * \copydoc FvBaseProblem::finishInit
     */
    void finishInit()
    {
        ParentType::finishInit();

        if (Opm::Valgrind::IsRunning())
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
        // three-phase Parker -- van Genuchten law
        fineMaterialParams_.setVgAlpha(0.0005);
        coarseMaterialParams_.setVgAlpha(0.005);
        fineMaterialParams_.setVgN(4.0);
        coarseMaterialParams_.setVgN(4.0);

        coarseMaterialParams_.setkrRegardsSnr(true);
        fineMaterialParams_.setkrRegardsSnr(true);

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
        fineMaterialParams_.setPcMinSat(gasPhaseIdx, 0);
        fineMaterialParams_.setPcMaxSat(gasPhaseIdx, 0);
        fineMaterialParams_.setPcMinSat(naplPhaseIdx, 0);
        fineMaterialParams_.setPcMaxSat(naplPhaseIdx, -1000);
        fineMaterialParams_.setPcMinSat(waterPhaseIdx, 0);
        fineMaterialParams_.setPcMaxSat(waterPhaseIdx, -10000);

        coarseMaterialParams_.setPcMinSat(gasPhaseIdx, 0);
        coarseMaterialParams_.setPcMaxSat(gasPhaseIdx, 0);
        coarseMaterialParams_.setPcMinSat(naplPhaseIdx, 0);
        coarseMaterialParams_.setPcMaxSat(naplPhaseIdx, -100);
        coarseMaterialParams_.setPcMinSat(waterPhaseIdx, 0);
        coarseMaterialParams_.setPcMaxSat(waterPhaseIdx, -1000);

        // residual saturations
        fineMaterialParams_.setResidSat(waterPhaseIdx, 0.1201);
        fineMaterialParams_.setResidSat(naplPhaseIdx, 0.0701);
        fineMaterialParams_.setResidSat(gasPhaseIdx, 0.0101);

        coarseMaterialParams_.setResidSat(waterPhaseIdx, 0.1201);
        coarseMaterialParams_.setResidSat(naplPhaseIdx, 0.0701);
        coarseMaterialParams_.setResidSat(gasPhaseIdx, 0.0101);
#endif

        fineMaterialParams_.finalize();
        coarseMaterialParams_.finalize();

        // initialize parameters for the thermal conduction law
        computeThermalCondParams_(thermalCondParams_, finePorosity_);

        // assume constant volumetric heat capacity and granite
        solidEnergyLawParams_.setSolidHeatCapacity(790.0 // specific heat capacity of granite [J / (kg K)]
                                                   * 2700.0); // density of granite [kg/m^3]
        solidEnergyLawParams_.finalize();

        initInjectFluidState_();
    }

    /*!
     * \name Auxiliary methods
     */
    //! \{

    /*!
     * \copydoc FvBaseProblem::shouldWriteRestartFile
     *
     * This problem writes a restart file after every time step.
     */
    bool shouldWriteRestartFile() const
    { return true; }

    /*!
     * \copydoc FvBaseProblem::name
     */
    std::string name() const
    { return std::string("cuvette_") + Model::name(); }

    /*!
     * \copydoc FvBaseProblem::endTimeStep
     */
    void endTimeStep()
    {
#ifndef NDEBUG
        this->model().checkConservativeness();

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
     * \name Soil parameters
     */
    //! \{

    /*!
     * \copydoc FvBaseMultiPhaseProblem::temperature
     */
    template <class Context>
    Scalar temperature(const Context& context OPM_UNUSED,
                       unsigned spaceIdx OPM_UNUSED,
                       unsigned timeIdx OPM_UNUSED) const
    { return 293.15; /* [K] */ }

    /*!
     * \copydoc FvBaseMultiPhaseProblem::intrinsicPermeability
     */
    template <class Context>
    const DimMatrix& intrinsicPermeability(const Context& context, unsigned spaceIdx,
                                           unsigned timeIdx) const
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
                                               unsigned spaceIdx, unsigned timeIdx) const
    {
        const GlobalPosition& pos = context.pos(spaceIdx, timeIdx);
        if (isFineMaterial_(pos))
            return fineMaterialParams_;
        else
            return coarseMaterialParams_;
    }

    /*!
     * \copydoc FvBaseMultiPhaseProblem::thermalConductionParams
     */
    template <class Context>
    const ThermalConductionLawParams &
    thermalConductionParams(const Context& context OPM_UNUSED,
                         unsigned spaceIdx OPM_UNUSED,
                         unsigned timeIdx OPM_UNUSED) const
    { return thermalCondParams_; }

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
        const auto& pos = context.pos(spaceIdx, timeIdx);

        if (onRightBoundary_(pos)) {
            Opm::CompositionalFluidState<Scalar, FluidSystem> fs;

            initialFluidState_(fs, context, spaceIdx, timeIdx);

            values.setFreeFlow(context, spaceIdx, timeIdx, fs);
            values.setNoFlow();
        }
        else if (onLeftBoundary_(pos)) {
            // injection
            RateVector molarRate;

            // inject with the same composition as the gas phase of
            // the injection fluid state
            Scalar molarInjectionRate = 0.3435; // [mol/(m^2 s)]
            for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx)
                molarRate[conti0EqIdx + compIdx] =
                    -molarInjectionRate
                    * injectFluidState_.moleFraction(gasPhaseIdx, compIdx);

            // calculate the total mass injection rate [kg / (m^2 s)
            Scalar massInjectionRate =
                molarInjectionRate
                * injectFluidState_.averageMolarMass(gasPhaseIdx);

            // set the boundary rate vector [J / (m^2 s)]
            values.setMolarRate(molarRate);
            values.setEnthalpyRate(-injectFluidState_.enthalpy(gasPhaseIdx) * massInjectionRate);
        }
        else
            values.setNoFlow();
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
    void initial(PrimaryVariables& values, const Context& context, unsigned spaceIdx,
                 unsigned timeIdx) const
    {
        Opm::CompositionalFluidState<Scalar, FluidSystem> fs;

        initialFluidState_(fs, context, spaceIdx, timeIdx);

        const auto& matParams = materialLawParams(context, spaceIdx, timeIdx);
        values.assignMassConservative(fs, matParams, /*inEquilibrium=*/false);
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
    { rate = Scalar(0.0); }

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

    bool isContaminated_(const GlobalPosition& pos) const
    {
        return (0.20 <= pos[0]) && (pos[0] <= 0.80) && (0.4 <= pos[1])
               && (pos[1] <= 0.65);
    }

    bool isFineMaterial_(const GlobalPosition& pos) const
    {
        if (0.13 <= pos[0] && 1.20 >= pos[0] && 0.32 <= pos[1] && pos[1] <= 0.57)
            return true;
        else if (pos[1] <= 0.15 && 1.20 <= pos[0])
            return true;
        else
            return false;
    }

    template <class FluidState, class Context>
    void initialFluidState_(FluidState& fs, const Context& context,
                            unsigned spaceIdx, unsigned timeIdx) const
    {
        const GlobalPosition& pos = context.pos(spaceIdx, timeIdx);

        fs.setTemperature(293.0 /*[K]*/);

        Scalar pw = 1e5;

        if (isContaminated_(pos)) {
            fs.setSaturation(waterPhaseIdx, 0.12);
            fs.setSaturation(naplPhaseIdx, 0.07);
            fs.setSaturation(gasPhaseIdx, 1 - 0.12 - 0.07);

            // set the capillary pressures
            const auto& matParams = materialLawParams(context, spaceIdx, timeIdx);
            Scalar pc[numPhases];
            MaterialLaw::capillaryPressures(pc, matParams, fs);
            for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
                fs.setPressure(phaseIdx, pw + (pc[phaseIdx] - pc[waterPhaseIdx]));

            // compute the phase compositions
            typedef Opm::MiscibleMultiPhaseComposition<Scalar, FluidSystem> MMPC;
            typename FluidSystem::template ParameterCache<Scalar> paramCache;
            MMPC::solve(fs, paramCache, /*setViscosity=*/true, /*setEnthalpy=*/true);
        }
        else {
            fs.setSaturation(waterPhaseIdx, 0.12);
            fs.setSaturation(gasPhaseIdx, 1 - fs.saturation(waterPhaseIdx));
            fs.setSaturation(naplPhaseIdx, 0);

            // set the capillary pressures
            const auto& matParams = materialLawParams(context, spaceIdx, timeIdx);
            Scalar pc[numPhases];
            MaterialLaw::capillaryPressures(pc, matParams, fs);
            for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
                fs.setPressure(phaseIdx, pw + (pc[phaseIdx] - pc[waterPhaseIdx]));

            // compute the phase compositions
            typedef Opm::MiscibleMultiPhaseComposition<Scalar, FluidSystem> MMPC;
            typename FluidSystem::template ParameterCache<Scalar> paramCache;
            MMPC::solve(fs, paramCache, /*setViscosity=*/true, /*setEnthalpy=*/true);

            // set the contaminant mole fractions to zero. this is a little bit hacky...
            for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                fs.setMoleFraction(phaseIdx, NAPLIdx, 0.0);

                if (phaseIdx == naplPhaseIdx)
                    continue;

                Scalar sumx = 0;
                for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx)
                    sumx += fs.moleFraction(phaseIdx, compIdx);

                for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx)
                    fs.setMoleFraction(phaseIdx, compIdx,
                                       fs.moleFraction(phaseIdx, compIdx) / sumx);
            }
        }
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
                Scalar lambdaFluid = FluidSystem::thermalConductivity(fs, paramCache, phaseIdx);
                lambdaSaturated =
                    std::pow(lambdaGranite, (1 - poro))
                    +
                    std::pow(lambdaFluid, poro);
            }
            else
                lambdaSaturated = std::pow(lambdaGranite, (1 - poro));

            params.setFullySaturatedLambda(phaseIdx, lambdaSaturated);
            if (!FluidSystem::isLiquid(phaseIdx))
                params.setVacuumLambda(lambdaSaturated);
        }
    }

    void initInjectFluidState_()
    {
        injectFluidState_.setTemperature(383.0);         // [K]
        injectFluidState_.setPressure(gasPhaseIdx, 1e5);   // [Pa]
        injectFluidState_.setSaturation(gasPhaseIdx, 1.0); // [-]

        Scalar xgH2O = 0.417;
        injectFluidState_.setMoleFraction(gasPhaseIdx, H2OIdx, xgH2O);     // [-]
        injectFluidState_.setMoleFraction(gasPhaseIdx, airIdx, 1 - xgH2O); // [-]
        injectFluidState_.setMoleFraction(gasPhaseIdx, NAPLIdx, 0.0);      // [-]

        // set the specific enthalpy of the gas phase
        typename FluidSystem::template ParameterCache<Scalar> paramCache;
        paramCache.updatePhase(injectFluidState_, gasPhaseIdx);

        Scalar h = FluidSystem::enthalpy(injectFluidState_, paramCache, gasPhaseIdx);
        injectFluidState_.setEnthalpy(gasPhaseIdx, h);
    }

    DimMatrix fineK_;
    DimMatrix coarseK_;

    Scalar finePorosity_;
    Scalar coarsePorosity_;

    MaterialLawParams fineMaterialParams_;
    MaterialLawParams coarseMaterialParams_;

    ThermalConductionLawParams thermalCondParams_;
    SolidEnergyLawParams solidEnergyLawParams_;

    Opm::CompositionalFluidState<Scalar, FluidSystem> injectFluidState_;

    const Scalar eps_;
};
} // namespace Opm

#endif
