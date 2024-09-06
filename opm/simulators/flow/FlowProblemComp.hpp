// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  Copyright 2024 SINTEF Digital
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
 * \copydoc Opm::FlowProblemComp
 */
#ifndef OPM_FLOW_PROBLEM_COMP_HPP
#define OPM_FLOW_PROBLEM_COMP_HPP


#include <opm/simulators/flow/FlowProblem.hpp>
#include <opm/material/fluidstates/CompositionalFluidState.hpp>

#include <opm/material/thermal/EclThermalLawManager.hpp>


#include <algorithm>
#include <functional>
#include <set>
#include <string>
#include <vector>

namespace Opm {

/*!
 * \ingroup CompositionalSimulator
 *
 * \brief This problem simulates an input file given in the data format used by the
 *        commercial ECLiPSE simulator.
 */
template <class TypeTag>
class FlowProblemComp : public FlowProblem<TypeTag>
{
    // TODO: the naming of the Types will be adjusted
    using FlowProblemType = FlowProblem<TypeTag>;

    using typename FlowProblemType::Scalar;
    using typename FlowProblemType::Simulator;
    using typename FlowProblemType::GridView;
    using typename FlowProblemType::FluidSystem;
    using typename FlowProblemType::Vanguard;

    // might not be needed
    using FlowProblemType::dim;
    using FlowProblemType::dimWorld;

    using InitialFluidState = CompositionalFluidState<Scalar, FluidSystem>;

public:
    using FlowProblemType::porosity;
    using FlowProblemType::pvtRegionIndex;

    /*!
     * \copydoc FvBaseProblem::registerParameters
     */
    static void registerParameters()
    {
        FlowProblemType::registerParameters();
    }


    /*!
     * \copydoc Doxygen::defaultProblemConstructor
     */
    explicit FlowProblemComp(Simulator& simulator)
        : FlowProblemType(simulator)
    {
    }

    /*!
     * \copydoc FvBaseProblem::finishInit
     */
    void finishInit()
    {
        // TODO: we need to do differently from the black oil version
        // because the way we handle the transmissiblity
        // TODO: we come back to refine this function
        FlowProblemType::finishInit();

        auto& simulator = this->simulator();

        auto finishTransmissibilities = [updated = false, this]() mutable
        {
            if (updated) { return; }
            this->transmissibilities_.finishInit([&vg = this->simulator().vanguard()](const unsigned int it) {
                return vg.gridIdxToEquilGridIdx(it);
            });
            updated = true;
        };

        // TODO: we call this for compostional
        finishTransmissibilities();

        // TODO: test whether we need to remove the following
        const auto& eclState = simulator.vanguard().eclState();
        const auto& schedule = simulator.vanguard().schedule();

        // Set the start time of the simulation
        simulator.setStartTime(schedule.getStartTime());
        simulator.setEndTime(schedule.simTime(schedule.size() - 1));

        // We want the episode index to be the same as the report step index to make
        // things simpler, so we have to set the episode index to -1 because it is
        // incremented by endEpisode(). The size of the initial time step and
        // length of the initial episode is set to zero for the same reason.
        simulator.setEpisodeIndex(-1);
        simulator.setEpisodeLength(0.0);

        // the "NOGRAV" keyword from Frontsim or setting the EnableGravity to false
        // disables gravity, else the standard value of the gravity constant at sea level
        // on earth is used
        this->gravity_ = 0.0;
        if (Parameters::Get<Parameters::EnableGravity>())
            this->gravity_[dim - 1] = 9.80665;
        if (!eclState.getInitConfig().hasGravity())
            this->gravity_[dim - 1] = 0.0;

        if (this->enableTuning_) {
            // if support for the TUNING keyword is enabled, we get the initial time
            // steping parameters from it instead of from command line parameters
            const auto& tuning = schedule[0].tuning();
            this->initialTimeStepSize_ = tuning.TSINIT.has_value() ? tuning.TSINIT.value() : -1.0;
            this->maxTimeStepAfterWellEvent_ = tuning.TMAXWC;
        }

        this->initFluidSystem_();

        if (FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx) &&
            FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)) {
            this->maxOilSaturation_.resize(this->model().numGridDof(), 0.0);
        }

        this->readRockParameters_(simulator.vanguard().cellCenterDepths(),
                                  [&simulator](const unsigned idx)
                                  {
                                      std::array<int,dim> coords;
                                      simulator.vanguard().cartesianCoordinate(idx, coords);
                                      for (auto& c : coords) {
                                          ++c;
                                      }
                                      return coords;
                                  });
        readMaterialParameters_();
        readThermalParameters_();

        const auto& initconfig = eclState.getInitConfig();
        if (initconfig.restartRequested())
            readEclRestartSolution_();
        else
            readInitialCondition_();

        updatePffDofData_();

        if constexpr (getPropValue<TypeTag, Properties::EnablePolymer>()) {
            const auto& vanguard = this->simulator().vanguard();
            const auto& gridView = vanguard.gridView();
            int numElements = gridView.size(/*codim=*/0);
            this->polymer_.maxAdsorption.resize(numElements, 0.0);
        }

        readBoundaryConditions_();

        // compute and set eq weights based on initial b values
        computeAndSetEqWeights_();

        if (enableDriftCompensation_) {
            drift_.resize(this->model().numGridDof());
            drift_ = 0.0;
        }

	// TODO: check wether the following can work with compostional
        if (enableVtkOutput_ && eclState.getIOConfig().initOnly()) {
            simulator.setTimeStepSize(0.0);
            FlowProblemType::writeOutput(true);
        }

        // after finishing the initialization and writing the initial solution, we move
        // to the first "real" episode/report step
        // for restart the episode index and start is already set
        if (!initconfig.restartRequested()) {
            simulator.startNextEpisode(schedule.seconds(1));
            simulator.setEpisodeIndex(0);
            simulator.setTimeStepIndex(0);
        }
    }

    /* void finalizeOutput() {
        OPM_TIMEBLOCK(finalizeOutput);
    } */

    /*!
     * \copydoc EclTransmissiblity::transmissibilityBoundary
     */
    template <class Context>
    Scalar transmissibilityBoundary(const Context& elemCtx,
                                    unsigned boundaryFaceIdx) const
    {
        unsigned elemIdx = elemCtx.globalSpaceIndex(/*dofIdx=*/0, /*timeIdx=*/0);
        return transmissibilities_.transmissibilityBoundary(elemIdx, boundaryFaceIdx);
    }

    /*!
     * \brief Direct access to a boundary transmissibility.
     */
    Scalar transmissibilityBoundary(const unsigned globalSpaceIdx,
                                    const unsigned boundaryFaceIdx) const
    {
        return transmissibilities_.transmissibilityBoundary(globalSpaceIdx, boundaryFaceIdx);
    }


    /*!
     * \copydoc EclTransmissiblity::thermalHalfTransmissibility
     */
    Scalar thermalHalfTransmissibility(const unsigned globalSpaceIdxIn,
                                       const unsigned globalSpaceIdxOut) const
    {
        return transmissibilities_.thermalHalfTrans(globalSpaceIdxIn,globalSpaceIdxOut);
    }

    /*!
     * \copydoc EclTransmissiblity::thermalHalfTransmissibility
     */
    template <class Context>
    Scalar thermalHalfTransmissibilityIn(const Context& context,
                                         unsigned faceIdx,
                                         unsigned timeIdx) const
    {
        const auto& face = context.stencil(timeIdx).interiorFace(faceIdx);
        unsigned toDofLocalIdx = face.exteriorIndex();
        return *pffDofData_.get(context.element(), toDofLocalIdx).thermalHalfTransIn;
    }

    /*!
     * \copydoc EclTransmissiblity::thermalHalfTransmissibility
     */
    template <class Context>
    Scalar thermalHalfTransmissibilityOut(const Context& context,
                                          unsigned faceIdx,
                                          unsigned timeIdx) const
    {
        const auto& face = context.stencil(timeIdx).interiorFace(faceIdx);
        unsigned toDofLocalIdx = face.exteriorIndex();
        return *pffDofData_.get(context.element(), toDofLocalIdx).thermalHalfTransOut;
    }

    /*!
     * \copydoc EclTransmissiblity::thermalHalfTransmissibility
     */
    template <class Context>
    Scalar thermalHalfTransmissibilityBoundary(const Context& elemCtx,
                                               unsigned boundaryFaceIdx) const
    {
        unsigned elemIdx = elemCtx.globalSpaceIndex(/*dofIdx=*/0, /*timeIdx=*/0);
        return transmissibilities_.thermalHalfTransBoundary(elemIdx, boundaryFaceIdx);
    }

    /*!
     * \brief Return a reference to the object that handles the "raw" transmissibilities.
     */
    const typename Vanguard::TransmissibilityType& eclTransmissibilities() const
    { return transmissibilities_; }

    /*!
     * \copydoc BlackOilBaseProblem::thresholdPressure
     */
    Scalar thresholdPressure(unsigned elem1Idx, unsigned elem2Idx) const
    { return thresholdPressures_.thresholdPressure(elem1Idx, elem2Idx); }

    const FlowThresholdPressure<TypeTag>& thresholdPressure() const
    { return thresholdPressures_; }

    FlowThresholdPressure<TypeTag>& thresholdPressure()
    { return thresholdPressures_; }

    /*!
     * \copydoc FvBaseMultiPhaseProblem::porosity
     *
     * For the FlowProblemBlackoil, this method is identical to referencePorosity(). The intensive
     * quantities object may apply various multipliers (e.g. ones which model rock
     * compressibility and water induced rock compaction) to it which depend on the
     * current physical conditions.
     */
    template <class Context>
    Scalar porosity(const Context& context, unsigned spaceIdx, unsigned timeIdx) const
    {
        unsigned globalSpaceIdx = context.globalSpaceIndex(spaceIdx, timeIdx);
        return this->porosity(globalSpaceIdx, timeIdx);
    }

    /*!
     * \brief Returns the depth of an degree of freedom [m]
     *
     * For ECL problems this is defined as the average of the depth of an element and is
     * thus slightly different from the depth of an element's centroid.
     */
    template <class Context>
    Scalar dofCenterDepth(const Context& context, unsigned spaceIdx, unsigned timeIdx) const
    {
        unsigned globalSpaceIdx = context.globalSpaceIndex(spaceIdx, timeIdx);
        return this->dofCenterDepth(globalSpaceIdx);
    }

    /*!
     * \brief Direct indexed acces to the depth of an degree of freedom [m]
     *
     * For ECL problems this is defined as the average of the depth of an element and is
     * thus slightly different from the depth of an element's centroid.
     */
    Scalar dofCenterDepth(unsigned globalSpaceIdx) const
    {
        return this->simulator().vanguard().cellCenterDepth(globalSpaceIdx);
    }

    /*!
     * \copydoc BlackoilProblem::rockCompressibility
     */
    template <class Context>
    Scalar rockCompressibility(const Context& context, unsigned spaceIdx, unsigned timeIdx) const
    {
        unsigned globalSpaceIdx = context.globalSpaceIndex(spaceIdx, timeIdx);
        return this->rockCompressibility(globalSpaceIdx);
    }

    /*!
     * \copydoc BlackoilProblem::rockReferencePressure
     */
    template <class Context>
    Scalar rockReferencePressure(const Context& context, unsigned spaceIdx, unsigned timeIdx) const
    {
        unsigned globalSpaceIdx = context.globalSpaceIndex(spaceIdx, timeIdx);
        return this->rockReferencePressure(globalSpaceIdx);
    }

    /*!
     * \copydoc FvBaseMultiPhaseProblem::materialLawParams
     */
    template <class Context>
    const MaterialLawParams& materialLawParams(const Context& context,
                                               unsigned spaceIdx, unsigned timeIdx) const
    {
        unsigned globalSpaceIdx = context.globalSpaceIndex(spaceIdx, timeIdx);
        return this->materialLawParams(globalSpaceIdx);
    }

    const MaterialLawParams& materialLawParams(unsigned globalDofIdx) const
    {
        return materialLawManager_->materialLawParams(globalDofIdx);
    }

    const MaterialLawParams& materialLawParams(unsigned globalDofIdx, FaceDir::DirEnum facedir) const
    {
        return materialLawManager_->materialLawParams(globalDofIdx, facedir);
    }

    /*!
     * \brief Return the parameters for the energy storage law of the rock
     */
    template <class Context>
    const SolidEnergyLawParams&
    solidEnergyLawParams(const Context& context,
                         unsigned spaceIdx,
                         unsigned timeIdx) const
    {
        unsigned globalSpaceIdx = context.globalSpaceIndex(spaceIdx, timeIdx);
        return thermalLawManager_->solidEnergyLawParams(globalSpaceIdx);
    }

    /*!
     * \copydoc FvBaseMultiPhaseProblem::thermalConductionParams
     */
    template <class Context>
    const ThermalConductionLawParams &
    thermalConductionLawParams(const Context& context, unsigned spaceIdx, unsigned timeIdx) const
    {
        unsigned globalSpaceIdx = context.globalSpaceIndex(spaceIdx, timeIdx);
        return thermalLawManager_->thermalConductionLawParams(globalSpaceIdx);
    }

    /*!
     * \brief Returns the ECL material law manager
     *
     * Note that this method is *not* part of the generic eWoms problem API because it
     * would force all problens use the ECL material laws.
     */
    std::shared_ptr<const EclMaterialLawManager> materialLawManager() const
    { return materialLawManager_; }

    template <class FluidState>
    void updateRelperms(
        std::array<Evaluation,numPhases> &mobility,
        DirectionalMobilityPtr &dirMob,
        FluidState &fluidState,
        unsigned globalSpaceIdx) const
    {
        OPM_TIMEBLOCK_LOCAL(updateRelperms);
        {
            // calculate relative permeabilities. note that we store the result into the
            // mobility_ class attribute. the division by the phase viscosity happens later.
            const auto& materialParams = materialLawParams(globalSpaceIdx);
            MaterialLaw::relativePermeabilities(mobility, materialParams, fluidState);
            Valgrind::CheckDefined(mobility);
        }
        if (materialLawManager_->hasDirectionalRelperms()
               || materialLawManager_->hasDirectionalImbnum())
        {
            using Dir = FaceDir::DirEnum;
            constexpr int ndim = 3;
            dirMob = std::make_unique<DirectionalMobility<TypeTag, Evaluation>>();
            Dir facedirs[ndim] = {Dir::XPlus, Dir::YPlus, Dir::ZPlus};
            for (int i = 0; i<ndim; i++) {
                const auto& materialParams = materialLawParams(globalSpaceIdx, facedirs[i]);
                auto& mob_array = dirMob->getArray(i);
                MaterialLaw::relativePermeabilities(mob_array, materialParams, fluidState);
            }
        }
    }

    /*!
     * \copydoc materialLawManager()
     */
    std::shared_ptr<EclMaterialLawManager> materialLawManager()
    { return materialLawManager_; }

    using BaseType::pvtRegionIndex;
    /*!
     * \brief Returns the index of the relevant region for thermodynmic properties
     */
    template <class Context>
    unsigned pvtRegionIndex(const Context& context, unsigned spaceIdx, unsigned timeIdx) const
    { return pvtRegionIndex(context.globalSpaceIndex(spaceIdx, timeIdx)); }

    using BaseType::satnumRegionIndex;
    /*!
     * \brief Returns the index of the relevant region for thermodynmic properties
     */
    template <class Context>
    unsigned satnumRegionIndex(const Context& context, unsigned spaceIdx, unsigned timeIdx) const
    { return this->satnumRegionIndex(context.globalSpaceIndex(spaceIdx, timeIdx)); }

    using BaseType::miscnumRegionIndex;
    /*!
     * \brief Returns the index of the relevant region for thermodynmic properties
     */
    template <class Context>
    unsigned miscnumRegionIndex(const Context& context, unsigned spaceIdx, unsigned timeIdx) const
    { return this->miscnumRegionIndex(context.globalSpaceIndex(spaceIdx, timeIdx)); }

    using BaseType::plmixnumRegionIndex;
    /*!
     * \brief Returns the index of the relevant region for thermodynmic properties
     */
    template <class Context>
    unsigned plmixnumRegionIndex(const Context& context, unsigned spaceIdx, unsigned timeIdx) const
    { return this->plmixnumRegionIndex(context.globalSpaceIndex(spaceIdx, timeIdx)); }

    using BaseType::maxPolymerAdsorption;
    /*!
     * \brief Returns the max polymer adsorption value
     */
    template <class Context>
    Scalar maxPolymerAdsorption(const Context& context, unsigned spaceIdx, unsigned timeIdx) const
    { return this->maxPolymerAdsorption(context.globalSpaceIndex(spaceIdx, timeIdx)); }

    /*!
     * \copydoc FvBaseProblem::name
     */
    std::string name() const
    { return this->simulator().vanguard().caseName(); }

    /*!
     * \copydoc FvBaseMultiPhaseProblem::temperature
     */
    template <class Context>
    Scalar temperature(const Context& context, unsigned spaceIdx, unsigned timeIdx) const
    {
        // use the initial temperature of the DOF if temperature is not a primary
        // variable
        unsigned globalDofIdx = context.globalSpaceIndex(spaceIdx, timeIdx);
        return initialFluidStates_[globalDofIdx].temperature(/*phaseIdx=*/0);
    }


    Scalar temperature(unsigned globalDofIdx, unsigned /*timeIdx*/) const
    {
        // use the initial temperature of the DOF if temperature is not a primary
        // variable
         return initialFluidStates_[globalDofIdx].temperature(/*phaseIdx=*/0);
    }

    const SolidEnergyLawParams&
    solidEnergyLawParams(unsigned globalSpaceIdx,
                         unsigned /*timeIdx*/) const
    {
        return this->thermalLawManager_->solidEnergyLawParams(globalSpaceIdx);
    }
    const ThermalConductionLawParams &
    thermalConductionLawParams(unsigned globalSpaceIdx,
                               unsigned /*timeIdx*/)const
    {
        return this->thermalLawManager_->thermalConductionLawParams(globalSpaceIdx);
    }

    /*!
     * \copydoc FvBaseProblem::boundary
     *
     * Reservoir simulation uses no-flow conditions as default for all boundaries.
     */
     // TODO: the boundary function might need to be emptied for compositional
    template <class Context>
    void boundary(BoundaryRateVector& values,
                  const Context& context,
                  unsigned spaceIdx,
                  unsigned timeIdx) const
    {
    }


    /*!
     * \brief Returns an element's historic maximum oil phase saturation that was
     *        observed during the simulation.
     *
     * In this context, "historic" means the the time before the current timestep began.
     *
     * This is a bit of a hack from the conceptional point of view, but it is required to
     * match the results of the 'flow' and ECLIPSE 100 simulators.
     */
    Scalar maxOilSaturation(unsigned globalDofIdx) const
    {
        if (!this->vapparsActive(this->episodeIndex()))
            return 0.0;

        return this->maxOilSaturation_[globalDofIdx];
    }

    /*!
     * \brief Sets an element's maximum oil phase saturation observed during the
     *        simulation.
     *
     * In this context, "historic" means the the time before the current timestep began.
     *
     * This a hack on top of the maxOilSaturation() hack but it is currently required to
     * do restart externally. i.e. from the flow code.
     */
    void setMaxOilSaturation(unsigned globalDofIdx, Scalar value)
    {
        if (!this->vapparsActive(this->episodeIndex()))
            return;

        this->maxOilSaturation_[globalDofIdx] = value;
    }

    /*!
     * \copydoc FvBaseProblem::initial
     *
     * The reservoir problem uses a constant boundary condition for
     * the whole domain.
     */
    template <class Context>
    void initial(PrimaryVariables& values, const Context& context, unsigned spaceIdx, unsigned timeIdx) const
    {
        // TODO: checking why we need the following one
	    const unsigned globalDofIdx = context.globalSpaceIndex(spaceIdx, timeIdx);
        const auto& initial_fs = initialFluidStates_[globalDofIdx];
        Opm::CompositionalFluidState<Evaluation, FluidSystem> fs;
        using ComponentVector = Dune::FieldVector<Evaluation, numComponents>;
        for (unsigned p = 0; p < numPhases; ++p) { // TODO: assuming the phaseidx continuous
            ComponentVector evals;
            auto& last_eval = evals[numComponents - 1];
            last_eval = 1.;
            for (unsigned c = 0; c < numComponents - 1; ++c) {
                const auto val = initial_fs.moleFraction(p, c);
                const Evaluation eval = Evaluation::createVariable(val, c+1);
                evals[c] = eval;
                last_eval -= eval;
            }
            for (unsigned c = 0; c < numComponents; ++c) {
                fs.setMoleFraction(p, c, evals[c]);
            }

            // pressure
            const auto p_val = initial_fs.pressure(p);
            fs.setPressure(p, Evaluation::createVariable(p_val, 0));

            const auto sat_val = initial_fs.saturation(p);
            fs.setSaturation(p, sat_val);

            const auto temp_val = initial_fs.temperature(p);
            fs.setTemperature(temp_val);
        }

        {
            typename FluidSystem::template ParameterCache<Evaluation> paramCache;
            paramCache.updatePhase(fs, FluidSystem::oilPhaseIdx);
            paramCache.updatePhase(fs, FluidSystem::gasPhaseIdx);
            fs.setDensity(FluidSystem::oilPhaseIdx, FluidSystem::density(fs, paramCache, FluidSystem::oilPhaseIdx));
            fs.setDensity(FluidSystem::gasPhaseIdx, FluidSystem::density(fs, paramCache, FluidSystem::gasPhaseIdx));
            fs.setViscosity(FluidSystem::oilPhaseIdx, FluidSystem::viscosity(fs, paramCache, FluidSystem::oilPhaseIdx));
            fs.setViscosity(FluidSystem::gasPhaseIdx, FluidSystem::viscosity(fs, paramCache, FluidSystem::gasPhaseIdx));
        }

        // Set initial K and L
        for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
            const Evaluation Ktmp = fs.wilsonK_(compIdx);
            fs.setKvalue(compIdx, Ktmp);
        }

        const Evaluation& Ltmp = -1.0;
        fs.setLvalue(Ltmp);

        values.assignNaive(fs);
    }

    /*!
     * \copydoc FvBaseProblem::initialSolutionApplied()
     */
    void initialSolutionApplied()
    {
        // Calculate all intensive quantities.
        this->model().invalidateAndUpdateIntensiveQuantities(/*timeIdx*/0);

        // We also need the intensive quantities for timeIdx == 1
        // corresponding to the start of the current timestep, if we
        // do not use the storage cache, or if we cannot recycle the
        // first iteration storage.
        if (!this->model().enableStorageCache() || !this->recycleFirstIterationStorage()) {
            this->model().invalidateAndUpdateIntensiveQuantities(/*timeIdx*/1);
        }

        // initialize the wells. Note that this needs to be done after initializing the
        // intrinsic permeabilities and the after applying the initial solution because
        // the well model uses these...
        wellModel_.init();

        // let the object for threshold pressures initialize itself. this is done only at
        // this point, because determining the threshold pressures may require to access
        // the initial solution.
	// TODO: comment out the following for compositional
        // thresholdPressures_.finishInit();

        aquiferModel_.initialSolutionApplied();


        const bool invalidateFromHyst = updateHysteresis_();
        if (invalidateFromHyst) {
            OPM_TIMEBLOCK(beginTimeStepInvalidateIntensiveQuantities);
            this->model().invalidateAndUpdateIntensiveQuantities(/*timeIdx=*/0);
        }
    }

    /*!
     * \copydoc FvBaseProblem::source
     *
     * For this problem, the source term of all components is 0 everywhere.
     */
    template <class Context>
    void source(RateVector& rate,
                const Context& context,
                unsigned spaceIdx,
                unsigned timeIdx) const
    {
        const unsigned globalDofIdx = context.globalSpaceIndex(spaceIdx, timeIdx);
        source(rate, globalDofIdx, timeIdx);
    }

    void source(RateVector& rate,
                unsigned globalDofIdx,
                unsigned timeIdx) const
    {
        OPM_TIMEBLOCK_LOCAL(eclProblemSource);
        rate = 0.0;

        // Add well contribution to source here.
        wellModel_.computeTotalRatesForDof(rate, globalDofIdx);

        // convert the source term from the total mass rate of the
        // cell to the one per unit of volume as used by the model.
        for (unsigned eqIdx = 0; eqIdx < numEq; ++ eqIdx) {
            rate[eqIdx] /= this->model().dofTotalVolume(globalDofIdx);

            Valgrind::CheckDefined(rate[eqIdx]);
            assert(isfinite(rate[eqIdx]));
        }

        // Add non-well sources.
        addToSourceDense(rate, globalDofIdx, timeIdx);
    }

    void addToSourceDense(RateVector& rate,
                          unsigned globalDofIdx,
                          unsigned timeIdx) const
    {
    }

    /*!
     * \brief Returns a reference to the ECL well manager used by the problem.
     *
     * This can be used for inspecting wells outside of the problem.
     */
    const WellModel& wellModel() const
    { return wellModel_; }

    WellModel& wellModel()
    { return wellModel_; }

    const AquiferModel& aquiferModel() const
    { return aquiferModel_; }

    AquiferModel& mutableAquiferModel()
    { return aquiferModel_; }

    // temporary solution to facilitate output of initial state from flow
    const InitialFluidState& initialFluidState(unsigned globalDofIdx) const
    { return initialFluidStates_[globalDofIdx]; }

    bool nonTrivialBoundaryConditions() const
    { return nonTrivialBoundaryConditions_; }

    const InitialFluidState boundaryFluidState(unsigned globalDofIdx, const int directionId) const
    {
	// TODO: we might need to make this function empty for compositional (except the last line)
        OPM_TIMEBLOCK_LOCAL(boundaryFluidState);
        const auto& bcprop = this->simulator().vanguard().schedule()[this->episodeIndex()].bcprop;
        if (bcprop.size() > 0) {
            FaceDir::DirEnum dir = FaceDir::FromIntersectionIndex(directionId);

            // index == 0: no boundary conditions for this
            // global cell and direction
            if (bcindex_(dir)[globalDofIdx] == 0)
                return initialFluidStates_[globalDofIdx];

            const auto& bc = bcprop[bcindex_(dir)[globalDofIdx]];
            if (bc.bctype == BCType::DIRICHLET )
            {
                InitialFluidState fluidState;
                const int pvtRegionIdx = this->pvtRegionIndex(globalDofIdx);
                fluidState.setPvtRegionIndex(pvtRegionIdx);

                switch (bc.component) {
                    case BCComponent::OIL:
                        if (!FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx))
                            throw std::logic_error("oil is not active and you're trying to add oil BC");

                        fluidState.setSaturation(FluidSystem::oilPhaseIdx, 1.0);
                        break;
                    case BCComponent::GAS:
                        if (!FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx))
                            throw std::logic_error("gas is not active and you're trying to add gas BC");

                        fluidState.setSaturation(FluidSystem::gasPhaseIdx, 1.0);
                        break;
                        case BCComponent::WATER:
                        if (!FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx))
                            throw std::logic_error("water is not active and you're trying to add water BC");

                        fluidState.setSaturation(FluidSystem::waterPhaseIdx, 1.0);
                        break;
                    case BCComponent::SOLVENT:
                    case BCComponent::POLYMER:
                    case BCComponent::NONE:
                        throw std::logic_error("you need to specify a valid component (OIL, WATER or GAS) when DIRICHLET type is set in BC");
                        break;
                }
                fluidState.setTotalSaturation(1.0);
                double pressure = initialFluidStates_[globalDofIdx].pressure(refPressurePhaseIdx_());
                const auto pressure_input = bc.pressure;
                if (pressure_input) {
                    pressure = *pressure_input;
                }

                std::array<Scalar, numPhases> pc = {0};
                const auto& matParams = materialLawParams(globalDofIdx);
                MaterialLaw::capillaryPressures(pc, matParams, fluidState);
                Valgrind::CheckDefined(pressure);
                Valgrind::CheckDefined(pc);
                for (unsigned activePhaseIdx = 0; activePhaseIdx < FluidSystem::numActivePhases(); ++activePhaseIdx) {
                    const auto phaseIdx = FluidSystem::activeToCanonicalPhaseIdx(activePhaseIdx);

                    fluidState.setPc(phaseIdx, pc[phaseIdx]);
                    if (Indices::oilEnabled)
                        fluidState.setPressure(phaseIdx, pressure + (pc[phaseIdx] - pc[oilPhaseIdx]));
                    else if (Indices::gasEnabled)
                        fluidState.setPressure(phaseIdx, pressure + (pc[phaseIdx] - pc[gasPhaseIdx]));
                    else if (Indices::waterEnabled)
                        //single (water) phase
                        fluidState.setPressure(phaseIdx, pressure);
                }
                
                double temperature = initialFluidStates_[globalDofIdx].temperature(0); // we only have one temperature
                const auto temperature_input = bc.temperature;
                if(temperature_input)
                    temperature = *temperature_input;
                fluidState.setTemperature(temperature);

                if (FluidSystem::enableDissolvedGas()) {
                    fluidState.setRs(0.0);
                    fluidState.setRv(0.0);
                }
                if (FluidSystem::enableDissolvedGasInWater()) {
                    fluidState.setRsw(0.0);
                }
                if (FluidSystem::enableVaporizedWater())
                    fluidState.setRvw(0.0);

                for (unsigned activePhaseIdx = 0; activePhaseIdx < FluidSystem::numActivePhases(); ++activePhaseIdx) {
                    const auto phaseIdx = FluidSystem::activeToCanonicalPhaseIdx(activePhaseIdx);

                    const auto& b = FluidSystem::inverseFormationVolumeFactor(fluidState, phaseIdx, pvtRegionIdx);
                    fluidState.setInvB(phaseIdx, b);

                    const auto& rho = FluidSystem::density(fluidState, phaseIdx, pvtRegionIdx);
                    fluidState.setDensity(phaseIdx, rho);
                    if (enableEnergy) {
                        const auto& h = FluidSystem::enthalpy(fluidState, phaseIdx, pvtRegionIdx);
                        fluidState.setEnthalpy(phaseIdx, h);
                    }
                }
                fluidState.checkDefined();
                return fluidState;
            }
        }
        return initialFluidStates_[globalDofIdx];
    }

    /*!
     * \brief Propose the size of the next time step to the simulator.
     *
     * This method is only called if the Newton solver does converge, the simulator
     * automatically cuts the time step in half without consultating this method again.
     */
    Scalar nextTimeStepSize() const
    {
        OPM_TIMEBLOCK(nexTimeStepSize);
        // allow external code to do the timestepping
        if (this->nextTimeStepSize_ > 0.0)
            return this->nextTimeStepSize_;

        const auto& simulator = this->simulator();
        int episodeIdx = simulator.episodeIndex();

        // for the initial episode, we use a fixed time step size
        if (episodeIdx < 0)
            return this->initialTimeStepSize_;

        // ask the newton method for a suggestion. This suggestion will be based on how
        // well the previous time step converged. After that, apply the runtime time
        // stepping constraints.
        const auto& newtonMethod = this->model().newtonMethod();
        return limitNextTimeStepSize_(newtonMethod.suggestTimeStepSize(simulator.timeStepSize()));
    }

    /*!
     * \brief Calculate the porosity multiplier due to water induced rock compaction.
     *
     * TODO: The API of this is a bit ad-hoc, it would be better to use context objects.
     */
    template <class LhsEval>
    LhsEval rockCompPoroMultiplier(const IntensiveQuantities& intQuants, unsigned elementIdx) const
    {
        OPM_TIMEBLOCK_LOCAL(rockCompPoroMultiplier);
        if (this->rockCompPoroMult_.empty() && this->rockCompPoroMultWc_.empty())
            return 1.0;

        unsigned tableIdx = 0;
        if (!this->rockTableIdx_.empty())
            tableIdx = this->rockTableIdx_[elementIdx];

        const auto& fs = intQuants.fluidState();
        LhsEval effectivePressure = decay<LhsEval>(fs.pressure(refPressurePhaseIdx_()));
        if (!this->minRefPressure_.empty())
            // The pore space change is irreversible
            effectivePressure =
                min(decay<LhsEval>(fs.pressure(refPressurePhaseIdx_())),
                                   this->minRefPressure_[elementIdx]);

        if (!this->overburdenPressure_.empty())
            effectivePressure -= this->overburdenPressure_[elementIdx];


        if (!this->rockCompPoroMult_.empty()) {
            return this->rockCompPoroMult_[tableIdx].eval(effectivePressure, /*extrapolation=*/true);
        }

        // water compaction
        assert(!this->rockCompPoroMultWc_.empty());
        LhsEval SwMax = max(decay<LhsEval>(fs.saturation(waterPhaseIdx)), this->maxWaterSaturation_[elementIdx]);
        LhsEval SwDeltaMax = SwMax - initialFluidStates_[elementIdx].saturation(waterPhaseIdx);

        return this->rockCompPoroMultWc_[tableIdx].eval(effectivePressure, SwDeltaMax, /*extrapolation=*/true);
    }

    /*!
     * \brief Calculate the transmissibility multiplier due to water induced rock compaction.
     *
     * TODO: The API of this is a bit ad-hoc, it would be better to use context objects.
     */
    template <class LhsEval>
    LhsEval rockCompTransMultiplier(const IntensiveQuantities& intQuants, unsigned elementIdx) const
    {
        const bool implicit = !this->explicitRockCompaction_;
        return implicit ? this->simulator().problem().template computeRockCompTransMultiplier_<LhsEval>(intQuants, elementIdx)
                        : this->simulator().problem().getRockCompTransMultVal(elementIdx);
    }


    /*!
     * \brief Return the well transmissibility multiplier due to rock changues.
     */
    template <class LhsEval>
    LhsEval wellTransMultiplier(const IntensiveQuantities& intQuants, unsigned elementIdx) const
    {
        OPM_TIMEBLOCK_LOCAL(wellTransMultiplier);
        
        const bool implicit = !this->explicitRockCompaction_;
        double trans_mult = implicit ? this->simulator().problem().template computeRockCompTransMultiplier_<double>(intQuants, elementIdx)
                                     : this->simulator().problem().getRockCompTransMultVal(elementIdx);
        trans_mult *= this->simulator().problem().template permFactTransMultiplier<double>(intQuants);
    
        return trans_mult;
    }

//    std::pair<BCType, RateVector> boundaryCondition(const unsigned int globalSpaceIdx, const int directionId) const
//    {
//        // TODO: we might need to empty this function for compositional
//        OPM_TIMEBLOCK_LOCAL(boundaryCondition);
//        if (!nonTrivialBoundaryConditions_) {
//            return { BCType::NONE, RateVector(0.0) };
//        }
//        FaceDir::DirEnum dir = FaceDir::FromIntersectionIndex(directionId);
//        const auto& schedule = this->simulator().vanguard().schedule();
//        if (bcindex_(dir)[globalSpaceIdx] == 0) {
//            return { BCType::NONE, RateVector(0.0) };
//        }
//        if (schedule[this->episodeIndex()].bcprop.size() == 0) {
//            return { BCType::NONE, RateVector(0.0) };
//        }
//        const auto& bc = schedule[this->episodeIndex()].bcprop[bcindex_(dir)[globalSpaceIdx]];
//        if (bc.bctype!=BCType::RATE) {
//            return { bc.bctype, RateVector(0.0) };
//        }
//
//        RateVector rate = 0.0;
//        switch (bc.component) {
//        case BCComponent::OIL:
//            rate[Indices::canonicalToActiveComponentIndex(oilCompIdx)] = bc.rate;
//            break;
//        case BCComponent::GAS:
//            rate[Indices::canonicalToActiveComponentIndex(gasCompIdx)] = bc.rate;
//            break;
//        case BCComponent::WATER:
//            rate[Indices::canonicalToActiveComponentIndex(waterCompIdx)] = bc.rate;
//            break;
//        case BCComponent::SOLVENT:
//            if constexpr (!enableSolvent)
//                throw std::logic_error("solvent is disabled and you're trying to add solvent to BC");
//
//            rate[Indices::solventSaturationIdx] = bc.rate;
//            break;
//        case BCComponent::POLYMER:
//            if constexpr (!enablePolymer)
//                throw std::logic_error("polymer is disabled and you're trying to add polymer to BC");
//
//            rate[Indices::polymerConcentrationIdx] = bc.rate;
//            break;
//        case BCComponent::NONE:
//            throw std::logic_error("you need to specify the component when RATE type is set in BC");
//            break;
//        }
//        //TODO add support for enthalpy rate
//        return {bc.bctype, rate};
//	//         return { BCType::NONE, RateVector(0.0) };
//    }



    template<class Serializer>
    void serializeOp(Serializer& serializer)
    {
        serializer(static_cast<BaseType&>(*this));
        serializer(drift_);
        serializer(wellModel_);
        serializer(aquiferModel_);
        serializer(*materialLawManager_);
    }
private:
    Implementation& asImp_()
    { return *static_cast<Implementation *>(this); }
protected:
    void updateExplicitQuantities_()
    {
        // TODO: we might need to empty the this function for compositional
        OPM_TIMEBLOCK(updateExplicitQuantities);
        const bool invalidateFromMaxWaterSat = updateMaxWaterSaturation_();
        const bool invalidateFromMinPressure = updateMinPressure_();

        // update hysteresis and max oil saturation used in vappars
        const bool invalidateFromHyst = updateHysteresis_();
        const bool invalidateFromMaxOilSat = updateMaxOilSaturation_();


        // deal with DRSDT and DRVDT
        // const bool invalidateDRDT = this->asImp_().updateCompositionChangeLimits_();

        // the derivatives may have changed
        bool invalidateIntensiveQuantities
            = invalidateFromMaxWaterSat || invalidateFromMinPressure || invalidateFromHyst || invalidateFromMaxOilSat;// || invalidateDRDT;
        if (invalidateIntensiveQuantities) {
            OPM_TIMEBLOCK(beginTimeStepInvalidateIntensiveQuantities);
            this->model().invalidateAndUpdateIntensiveQuantities(/*timeIdx=*/0);
        }

        updateRockCompTransMultVal_();
    }

    template<class UpdateFunc>
    void updateProperty_(const std::string& failureMsg,
                         UpdateFunc func)
    {
        OPM_TIMEBLOCK(updateProperty);
        const auto& model = this->simulator().model();
        const auto& primaryVars = model.solution(/*timeIdx*/0);
        const auto& vanguard = this->simulator().vanguard();
        std::size_t numGridDof = primaryVars.size();
        OPM_BEGIN_PARALLEL_TRY_CATCH();
#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (unsigned dofIdx = 0; dofIdx < numGridDof; ++dofIdx) {
                const auto& iq = *model.cachedIntensiveQuantities(dofIdx, /*timeIdx=*/ 0);
                func(dofIdx, iq);
        }
        OPM_END_PARALLEL_TRY_CATCH(failureMsg, vanguard.grid().comm());
    }

    // update the parameters needed for DRSDT and DRVDT
    bool updateCompositionChangeLimits_()
    {
            return true;
    }


    bool updateMaxOilSaturation_()
    {
        OPM_TIMEBLOCK(updateMaxOilSaturation);
        int episodeIdx = this->episodeIndex();

        // we use VAPPARS
        if (this->vapparsActive(episodeIdx)) {
            this->updateProperty_("FlowProblemBlackoil::updateMaxOilSaturation_() failed:",
                                  [this](unsigned compressedDofIdx, const IntensiveQuantities& iq)
                                  {
                                      this->updateMaxOilSaturation_(compressedDofIdx,iq);
                                  });
            return true;
        }

        return false;
    }

    bool updateMaxOilSaturation_(unsigned compressedDofIdx, const IntensiveQuantities& iq)
    {
        OPM_TIMEBLOCK_LOCAL(updateMaxOilSaturation);
        const auto& fs = iq.fluidState();
        const Scalar So = decay<Scalar>(fs.saturation(refPressurePhaseIdx_()));
        auto& mos = this->maxOilSaturation_;
        if(mos[compressedDofIdx] < So){
            mos[compressedDofIdx] = So;
            return true;
        }else{
            return false;
        }
    }

    bool updateMaxWaterSaturation_()
    {
        OPM_TIMEBLOCK(updateMaxWaterSaturation);
        // water compaction is activated in ROCKCOMP
        if (this->maxWaterSaturation_.empty())
            return false;

        this->maxWaterSaturation_[/*timeIdx=*/1] = this->maxWaterSaturation_[/*timeIdx=*/0];
        this->updateProperty_("FlowProblemBlackoil::updateMaxWaterSaturation_() failed:",
                              [this](unsigned compressedDofIdx, const IntensiveQuantities& iq)
                              {
                                  this->updateMaxWaterSaturation_(compressedDofIdx,iq);
                               });
        return true;
    }


    bool updateMaxWaterSaturation_(unsigned compressedDofIdx, const IntensiveQuantities& iq)
    {
        OPM_TIMEBLOCK_LOCAL(updateMaxWaterSaturation);
        const auto& fs = iq.fluidState();
        const Scalar Sw = decay<Scalar>(fs.saturation(waterPhaseIdx));
        auto& mow = this->maxWaterSaturation_;
        if(mow[compressedDofIdx]< Sw){
            mow[compressedDofIdx] = Sw;
            return true;
        }else{
            return false;
        }
    }

    bool updateMinPressure_()
    {
        OPM_TIMEBLOCK(updateMinPressure);
        // IRREVERS option is used in ROCKCOMP
        if (this->minRefPressure_.empty())
            return false;

        this->updateProperty_("FlowProblemBlackoil::updateMinPressure_() failed:",
                              [this](unsigned compressedDofIdx, const IntensiveQuantities& iq)
                              {
                                  this->updateMinPressure_(compressedDofIdx,iq);
                              });
        return true;
    }

    bool updateMinPressure_(unsigned compressedDofIdx, const IntensiveQuantities& iq){
        OPM_TIMEBLOCK_LOCAL(updateMinPressure);
        const auto& fs = iq.fluidState();
        const Scalar min_pressure = getValue(fs.pressure(refPressurePhaseIdx_()));
        auto& min_pressures = this->minRefPressure_;
        if(min_pressures[compressedDofIdx]> min_pressure){
            min_pressures[compressedDofIdx] = min_pressure;
            return true;
        }else{
            return false;
        }
    }

    // \brief Function to assign field properties of type double, on the leaf grid view.
    //
    // For CpGrid with local grid refinement, the field property of a cell on the leaf
    // is inherited from its parent or equivalent (when has no parent) cell on level zero.
    std::function<std::vector<double>(const FieldPropsManager&, const std::string&)>
    fieldPropDoubleOnLeafAssigner_()
    {
        const auto& lookup = this->lookUpData_;
        return [&lookup](const FieldPropsManager& fieldPropManager, const std::string& propString)
        {
            return lookup.assignFieldPropsDoubleOnLeaf(fieldPropManager, propString);
        };
    }

    // \brief Function to assign field properties of type int, unsigned int, ..., on the leaf grid view.
    //
    // For CpGrid with local grid refinement, the field property of a cell on the leaf
    // is inherited from its parent or equivalent (when has no parent) cell on level zero.
    template<typename IntType>
    std::function<std::vector<IntType>(const FieldPropsManager&, const std::string&, bool)>
    fieldPropIntTypeOnLeafAssigner_()
    {
        const auto& lookup = this->lookUpData_;
        return [&lookup](const FieldPropsManager& fieldPropManager, const std::string& propString, bool needsTranslation)
        {
            return lookup.template assignFieldPropsIntOnLeaf<IntType>(fieldPropManager, propString, needsTranslation);
        };
    }

    void readMaterialParameters_()
    {
        OPM_TIMEBLOCK(readMaterialParameters);
        const auto& simulator = this->simulator();
        const auto& vanguard = simulator.vanguard();
        const auto& eclState = vanguard.eclState();

        // the PVT and saturation region numbers
        OPM_BEGIN_PARALLEL_TRY_CATCH();
        this->updatePvtnum_();
        this->updateSatnum_();

        // the MISC region numbers (solvent model)
        this->updateMiscnum_();
        // the PLMIX region numbers (polymer model)
        this->updatePlmixnum_();

        OPM_END_PARALLEL_TRY_CATCH("Invalid region numbers: ", vanguard.gridView().comm());
        ////////////////////////////////
        // porosity
        updateReferencePorosity_();
        this->referencePorosity_[1] = this->referencePorosity_[0];
        ////////////////////////////////

        ////////////////////////////////
        // fluid-matrix interactions (saturation functions; relperm/capillary pressure)
        materialLawManager_ = std::make_shared<EclMaterialLawManager>();
        materialLawManager_->initFromState(eclState);
        materialLawManager_->initParamsForElements(eclState, this->model().numGridDof(),
                                                   this-> template fieldPropIntTypeOnLeafAssigner_<int>(),
                                                   this-> lookupIdxOnLevelZeroAssigner_());
        ////////////////////////////////
    }

    void readThermalParameters_()
    {
        if constexpr (enableEnergy)
        {
            const auto& simulator = this->simulator();
            const auto& vanguard = simulator.vanguard();
            const auto& eclState = vanguard.eclState();

            // fluid-matrix interactions (saturation functions; relperm/capillary pressure)
            thermalLawManager_ = std::make_shared<EclThermalLawManager>();
            thermalLawManager_->initParamsForElements(eclState, this->model().numGridDof(),
                                                      this-> fieldPropDoubleOnLeafAssigner_(),
                                                      this-> template fieldPropIntTypeOnLeafAssigner_<unsigned int>());
        }
    }

    void updateReferencePorosity_()
    {
        const auto& simulator = this->simulator();
        const auto& vanguard = simulator.vanguard();
        const auto& eclState = vanguard.eclState();

        std::size_t numDof = this->model().numGridDof();

        this->referencePorosity_[/*timeIdx=*/0].resize(numDof);

        const auto& fp = eclState.fieldProps();
        const std::vector<double> porvData = this -> fieldPropDoubleOnLeafAssigner_()(fp, "PORV");
        for (std::size_t dofIdx = 0; dofIdx < numDof; ++dofIdx) {
            Scalar poreVolume = porvData[dofIdx];

            // we define the porosity as the accumulated pore volume divided by the
            // geometric volume of the element. Note that -- in pathetic cases -- it can
            // be larger than 1.0!
            Scalar dofVolume = simulator.model().dofTotalVolume(dofIdx);
            assert(dofVolume > 0.0);
            this->referencePorosity_[/*timeIdx=*/0][dofIdx] = poreVolume/dofVolume;
        }
    }

    void readInitialCondition_()
    {
        // TODO: empty this funcition and use the following two lines for this function for compositional
	// this->readExplicitInitialCondition_();
        // this->model().applyInitialSolution();
        const auto& simulator = this->simulator();
        const auto& vanguard = simulator.vanguard();
        const auto& eclState = vanguard.eclState();

        if (eclState.getInitConfig().hasEquil())
            readEquilInitialCondition_();
        else
            readExplicitInitialCondition_();

        if constexpr (enableSolvent || enablePolymer || enablePolymerMolarWeight || enableMICP)
            this->readBlackoilExtentionsInitialConditions_(this->model().numGridDof(),
                                                           enableSolvent,
                                                           enablePolymer,
                                                           enablePolymerMolarWeight,
                                                           enableMICP);

        //initialize min/max values
        std::size_t numElems = this->model().numGridDof();
        for (std::size_t elemIdx = 0; elemIdx < numElems; ++elemIdx) {
            const auto& fs = initialFluidStates_[elemIdx];
            if (!this->maxWaterSaturation_.empty())
                this->maxWaterSaturation_[elemIdx] = std::max(this->maxWaterSaturation_[elemIdx], fs.saturation(waterPhaseIdx));
            if (!this->maxOilSaturation_.empty())
                this->maxOilSaturation_[elemIdx] = std::max(this->maxOilSaturation_[elemIdx], fs.saturation(oilPhaseIdx));
            if (!this->minRefPressure_.empty())
                this->minRefPressure_[elemIdx] = std::min(this->minRefPressure_[elemIdx], fs.pressure(refPressurePhaseIdx_()));
        }


    }

    void readEquilInitialCondition_()
    {
    }

    void readEclRestartSolution_()
    {
    }

    void readExplicitInitialCondition_()
    {
        readExplicitInitialConditionCompositional_();
    }
    
    void readExplicitInitialConditionCompositional_()
    {
        const auto& simulator = this->simulator();
        const auto& vanguard = simulator.vanguard();
        const auto& eclState = vanguard.eclState();
        const auto& fp = eclState.fieldProps();
        // const bool has_xmf = fp.has_double("XMF");
        assert(fp.has_double("XMF"));
        // const bool has_ymf = fp.has_double("YMF");
        assert(fp.has_double("YMF"));
        const bool has_temp = fp.has_double("TEMPI");
        const bool has_pressure = fp.has_double("PRESSURE");

        // const bool has_gas = fp.has_double("SGAS");
        assert(fp.has_double("SGAS"));
        if (!has_pressure)
            throw std::runtime_error("The ECL input file requires the presence of the PRESSURE "
                                      "keyword if the model is initialized explicitly");

        std::size_t numDof = this->model().numGridDof();

        initialFluidStates_.resize(numDof);

        std::vector<double> xmfData;
        std::vector<double> ymfData;
        std::vector<double> waterSaturationData;
        std::vector<double> gasSaturationData;
        std::vector<double> soilData;
        std::vector<double> pressureData;
        std::vector<double> tempiData;

        if (waterPhaseIdx > 0 && Indices::numPhases > 1)
            waterSaturationData = fp.get_double("SWAT");
        else
            waterSaturationData.resize(numDof);

        pressureData = fp.get_double("PRESSURE");

        assert(waterPhaseIdx < 0);

        xmfData = fp.get_double("XMF");
        ymfData = fp.get_double("YMF");
        if (has_temp) {
            tempiData = fp.get_double("TEMPI");
        } else {
            ; // TODO: throw?
        }

        if (gasPhaseIdx > 0) // && FluidSystem::phaseIsActive(oilPhaseIdx))
            gasSaturationData = fp.get_double("SGAS");
        else
            gasSaturationData.resize(numDof);

        // constexpr std::size_t num_comps = 3;

        for (std::size_t dofIdx = 0; dofIdx < numDof; ++dofIdx) {
            auto& dofFluidState = initialFluidStates_[dofIdx];
            // dofFluidState.setPvtRegionIndex(pvtRegionIndex(dofIdx));

            Scalar temperatureLoc = tempiData[dofIdx];
            assert(std::isfinite(temperatureLoc) && temperatureLoc > 0);
            dofFluidState.setTemperature(temperatureLoc);

            if (gasPhaseIdx > 0) {
                    dofFluidState.setSaturation(FluidSystem::gasPhaseIdx,
                                                gasSaturationData[dofIdx]);
            }
            if (oilPhaseIdx > 0) {
                dofFluidState.setSaturation(FluidSystem::oilPhaseIdx,
                                            1.0
                                            - waterSaturationData[dofIdx]
                                            - gasSaturationData[dofIdx]);
            }

            //////
            // set phase pressures
            //////
            const Scalar pressure = pressureData[dofIdx]; // oil pressure (or gas pressure for water-gas system or water pressure for single phase)

            // TODO: zero capillary pressure for now
            const std::array<Scalar, numPhases> pc = {0};
            for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                if (!FluidSystem::phaseIsActive(phaseIdx))
                    continue;

                if (Indices::oilEnabled)
                    dofFluidState.setPressure(phaseIdx, pressure + (pc[phaseIdx] - pc[oilPhaseIdx]));
                else if (Indices::gasEnabled)
                    dofFluidState.setPressure(phaseIdx, pressure + (pc[phaseIdx] - pc[gasPhaseIdx]));
                else if (Indices::waterEnabled)
                    //single (water) phase
                    dofFluidState.setPressure(phaseIdx, pressure);
            }

            for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
                const std::size_t data_idx = compIdx * numDof + dofIdx;
                const Scalar xmf = xmfData[data_idx];
                const Scalar ymf = xmfData[data_idx];

                dofFluidState.setMoleFraction(FluidSystem::oilPhaseIdx, compIdx, xmf);
                dofFluidState.setMoleFraction(FluidSystem::gasPhaseIdx, compIdx, ymf);
            }
        }
    }

    // update the hysteresis parameters of the material laws for the whole grid
    bool updateHysteresis_()
    {
        if (!materialLawManager_->enableHysteresis())
            return false;

        // we need to update the hysteresis data for _all_ elements (i.e., not just the
        // interior ones) to avoid desynchronization of the processes in the parallel case!
        this->updateProperty_("FlowProblemBlackoil::updateHysteresis_() failed:",
                              [this](unsigned compressedDofIdx, const IntensiveQuantities& iq)
                              {
                                  materialLawManager_->updateHysteresis(iq.fluidState(), compressedDofIdx);
                              });
        return true;
    }


    bool updateHysteresis_(unsigned compressedDofIdx, const IntensiveQuantities& iq)
    {
        OPM_TIMEBLOCK_LOCAL(updateHysteresis_);
        materialLawManager_->updateHysteresis(iq.fluidState(), compressedDofIdx);
        //TODO change materials to give a bool
        return true;
    }

    Scalar getRockCompTransMultVal(std::size_t dofIdx) const
    {
        if (this->rockCompTransMultVal_.empty())
            return 1.0;

        return this->rockCompTransMultVal_[dofIdx];
    }


private:
    struct PffDofData_
    {
        ConditionalStorage<enableEnergy, Scalar> thermalHalfTransIn;
        ConditionalStorage<enableEnergy, Scalar> thermalHalfTransOut;
        ConditionalStorage<enableDiffusion, Scalar> diffusivity;
        ConditionalStorage<enableDispersion, Scalar> dispersivity;
        Scalar transmissibility;
    };

    // update the prefetch friendly data object
    void updatePffDofData_()
    {
        const auto& distFn =
            [this](PffDofData_& dofData,
                   const Stencil& stencil,
                   unsigned localDofIdx)
            -> void
        {
            const auto& elementMapper = this->model().elementMapper();

            unsigned globalElemIdx = elementMapper.index(stencil.entity(localDofIdx));
            if (localDofIdx != 0) {
                unsigned globalCenterElemIdx = elementMapper.index(stencil.entity(/*dofIdx=*/0));
                dofData.transmissibility = transmissibilities_.transmissibility(globalCenterElemIdx, globalElemIdx);

                if constexpr (enableEnergy) {
                    *dofData.thermalHalfTransIn = transmissibilities_.thermalHalfTrans(globalCenterElemIdx, globalElemIdx);
                    *dofData.thermalHalfTransOut = transmissibilities_.thermalHalfTrans(globalElemIdx, globalCenterElemIdx);
                }
                if constexpr (enableDiffusion)
                    *dofData.diffusivity = transmissibilities_.diffusivity(globalCenterElemIdx, globalElemIdx);
                if (enableDispersion)
                    dofData.dispersivity = transmissibilities_.dispersivity(globalCenterElemIdx, globalElemIdx);
            }
        };

        pffDofData_.update(distFn);
    }

    void readBoundaryConditions_()
    {
        const auto& simulator = this->simulator();
        const auto& vanguard = simulator.vanguard();
        const auto& bcconfig = vanguard.eclState().getSimulationConfig().bcconfig();
        if (bcconfig.size() > 0) {
            nonTrivialBoundaryConditions_ = true;

            std::size_t numCartDof = vanguard.cartesianSize();
            unsigned numElems = vanguard.gridView().size(/*codim=*/0);
            std::vector<int> cartesianToCompressedElemIdx(numCartDof, -1);

            for (unsigned elemIdx = 0; elemIdx < numElems; ++elemIdx)
                cartesianToCompressedElemIdx[vanguard.cartesianIndex(elemIdx)] = elemIdx;

            bcindex_.resize(numElems, 0);
            auto loopAndApply = [&cartesianToCompressedElemIdx,
                                 &vanguard](const auto& bcface,
                                            auto apply)
            {
                for (int i = bcface.i1; i <= bcface.i2; ++i) {
                    for (int j = bcface.j1; j <= bcface.j2; ++j) {
                        for (int k = bcface.k1; k <= bcface.k2; ++k) {
                            std::array<int, 3> tmp = {i,j,k};
                            auto elemIdx = cartesianToCompressedElemIdx[vanguard.cartesianIndex(tmp)];
                            if (elemIdx >= 0)
                                apply(elemIdx);
                        }
                    }
                }
            };
            for (const auto& bcface : bcconfig) {
                std::vector<int>& data = bcindex_(bcface.dir);
                const int index = bcface.index;
                    loopAndApply(bcface,
                                 [&data,index](int elemIdx)
                                 { data[elemIdx] = index; });
            }
        }
    }

    // this method applies the runtime constraints specified via the deck and/or command
    // line parameters for the size of the next time step.
    Scalar limitNextTimeStepSize_(Scalar dtNext) const
    {
        if constexpr (enableExperiments) {
            const auto& simulator = this->simulator();
            const auto& schedule = simulator.vanguard().schedule();
            int episodeIdx = simulator.episodeIndex();

            // first thing in the morning, limit the time step size to the maximum size
            Scalar maxTimeStepSize = Parameters::Get<Parameters::SolverMaxTimeStepInDays<Scalar>>() * 24 * 60 * 60;
            int reportStepIdx = std::max(episodeIdx, 0);
            if (this->enableTuning_) {
                const auto& tuning = schedule[reportStepIdx].tuning();
                maxTimeStepSize = tuning.TSMAXZ;
            }

            dtNext = std::min(dtNext, maxTimeStepSize);

            Scalar remainingEpisodeTime =
                simulator.episodeStartTime() + simulator.episodeLength()
                - (simulator.startTime() + simulator.time());
            assert(remainingEpisodeTime >= 0.0);

            // if we would have a small amount of time left over in the current episode, make
            // two equal time steps instead of a big and a small one
            if (remainingEpisodeTime/2.0 < dtNext && dtNext < remainingEpisodeTime*(1.0 - 1e-5))
                // note: limiting to the maximum time step size here is probably not strictly
                // necessary, but it should not hurt and is more fool-proof
                dtNext = std::min(maxTimeStepSize, remainingEpisodeTime/2.0);

            if (simulator.episodeStarts()) {
                // if a well event occurred, respect the limit for the maximum time step after
                // that, too
                const auto& events = simulator.vanguard().schedule()[reportStepIdx].events();
                bool wellEventOccured =
                        events.hasEvent(ScheduleEvents::NEW_WELL)
                        || events.hasEvent(ScheduleEvents::PRODUCTION_UPDATE)
                        || events.hasEvent(ScheduleEvents::INJECTION_UPDATE)
                        || events.hasEvent(ScheduleEvents::WELL_STATUS_CHANGE);
                if (episodeIdx >= 0 && wellEventOccured && this->maxTimeStepAfterWellEvent_ > 0)
                    dtNext = std::min(dtNext, this->maxTimeStepAfterWellEvent_);
            }
        }

        return dtNext;
    }

    void computeAndSetEqWeights_()
    {
    }

    int refPressurePhaseIdx_() const {
        if (FluidSystem::phaseIsActive(oilPhaseIdx)) {
            return oilPhaseIdx;
        }
        else if (FluidSystem::phaseIsActive(gasPhaseIdx)) {
            return gasPhaseIdx;
        }
        else {
            return waterPhaseIdx;
        }
    }

    void updateRockCompTransMultVal_()
    {
        const auto& model = this->simulator().model();
        std::size_t numGridDof = this->model().numGridDof();
        this->rockCompTransMultVal_.resize(numGridDof, 1.0);
        for (std::size_t elementIdx = 0; elementIdx < numGridDof; ++elementIdx) {
            const auto& iq = *model.cachedIntensiveQuantities(elementIdx, /*timeIdx=*/ 0);
            Scalar trans_mult = computeRockCompTransMultiplier_<Scalar>(iq, elementIdx);
            this->rockCompTransMultVal_[elementIdx] = trans_mult;
        }
    }

    /*!
     * \brief Calculate the transmissibility multiplier due to water induced rock compaction.
     *
     * TODO: The API of this is a bit ad-hoc, it would be better to use context objects.
     */
    template <class LhsEval>
    LhsEval computeRockCompTransMultiplier_(const IntensiveQuantities& intQuants, unsigned elementIdx) const
    {
        OPM_TIMEBLOCK_LOCAL(computeRockCompTransMultiplier);
        if (this->rockCompTransMult_.empty() && this->rockCompTransMultWc_.empty())
            return 1.0;
	// TODO: we might need to comment out the following lines
        unsigned tableIdx = 0;
        if (!this->rockTableIdx_.empty())
            tableIdx = this->rockTableIdx_[elementIdx];

        const auto& fs = intQuants.fluidState();
        LhsEval effectivePressure = decay<LhsEval>(fs.pressure(refPressurePhaseIdx_()));

        if (!this->minRefPressure_.empty())
            // The pore space change is irreversible
            effectivePressure =
                min(decay<LhsEval>(fs.pressure(refPressurePhaseIdx_())),
                    this->minRefPressure_[elementIdx]);

        if (!this->overburdenPressure_.empty())
            effectivePressure -= this->overburdenPressure_[elementIdx];

        if (!this->rockCompTransMult_.empty())
            return this->rockCompTransMult_[tableIdx].eval(effectivePressure, /*extrapolation=*/true);

        // water compaction
        assert(!this->rockCompTransMultWc_.empty());
        LhsEval SwMax = max(decay<LhsEval>(fs.saturation(waterPhaseIdx)), this->maxWaterSaturation_[elementIdx]);
        LhsEval SwDeltaMax = SwMax - initialFluidStates_[elementIdx].saturation(waterPhaseIdx);

        return this->rockCompTransMultWc_[tableIdx].eval(effectivePressure, SwDeltaMax, /*extrapolation=*/true);
    }

    typename Vanguard::TransmissibilityType transmissibilities_;

    std::shared_ptr<EclMaterialLawManager> materialLawManager_;
    std::shared_ptr<EclThermalLawManager> thermalLawManager_;

    FlowThresholdPressure<TypeTag> thresholdPressures_;

    std::vector<InitialFluidState> initialFluidStates_;

    bool enableDriftCompensation_;
    GlobalEqVector drift_;

    WellModel wellModel_;
    AquiferModel aquiferModel_;

    bool enableEclOutput_;
    bool enableVtkOutput_;

    PffGridVector<GridView, Stencil, PffDofData_, DofMapper> pffDofData_;

    template<class T>
    struct BCData
    {
        std::array<std::vector<T>,6> data;

        void resize(std::size_t size, T defVal)
        {
            for (auto& d : data)
                d.resize(size, defVal);
        }

        const std::vector<T>& operator()(FaceDir::DirEnum dir) const
        {
            if (dir == FaceDir::DirEnum::Unknown)
                throw std::runtime_error("Tried to access BC data for the 'Unknown' direction");
            int idx = 0;
            int div = static_cast<int>(dir);
            while ((div /= 2) >= 1)
              ++idx;
            assert(idx >= 0 && idx <= 5);
            return data[idx];
        }

        std::vector<T>& operator()(FaceDir::DirEnum dir)
        {
            return const_cast<std::vector<T>&>(std::as_const(*this)(dir));
        }
    };

    BCData<int> bcindex_;
    bool nonTrivialBoundaryConditions_ = false;
    bool explicitRockCompaction_ = false;
};

} // namespace Opm

#endif // OPM_FLOW_PROBLEM_COMP_HPP
