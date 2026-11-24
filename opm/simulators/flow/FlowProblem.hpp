// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  Copyright 2023 INRIA

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
 * \copydoc Opm::FlowProblem
 */
#ifndef OPM_FLOW_PROBLEM_HPP
#define OPM_FLOW_PROBLEM_HPP

#include <dune/common/version.hh>
#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

#include <opm/common/utility/TimeService.hpp>

#include <opm/input/eclipse/EclipseState/EclipseState.hpp>
#include <opm/input/eclipse/Schedule/Schedule.hpp>
#include <opm/input/eclipse/Units/Units.hpp>

#include <opm/material/common/ConditionalStorage.hpp>
#include <opm/material/common/Valgrind.hpp>
#include <opm/material/densead/Evaluation.hpp>
#include <opm/material/fluidmatrixinteractions/EclMaterialLawManager.hpp>
#include <opm/material/thermal/EclThermalLawManager.hpp>

#include <opm/models/blackoil/blackoilconvectivemixingmoduleparam.hpp>
#include <opm/models/blackoil/blackoilmoduleparams.hh>
#include <opm/models/common/directionalmobility.hh>
#include <opm/models/utils/pffgridvector.hh>
#include <opm/models/discretization/ecfv/ecfvdiscretization.hh>

#include <opm/output/eclipse/EclipseIO.hpp>

#include <opm/simulators/flow/CpGridVanguard.hpp>
#include <opm/simulators/flow/DummyGradientCalculator.hpp>
#include <opm/simulators/flow/EclWriter.hpp>
#include <opm/simulators/flow/EquilInitializer.hpp>
#include <opm/simulators/flow/FlowGenericProblem.hpp>
// TODO: maybe we can name it FlowProblemProperties.hpp
#include <opm/simulators/flow/FlowAuxCellModule.hpp>
#include <opm/simulators/flow/FlowBaseProblemProperties.hpp>
#include <opm/simulators/flow/FlowProblemParameters.hpp>
#include <opm/simulators/flow/FlowThresholdPressure.hpp>
#include <opm/simulators/flow/flux/FluxBoundary.hpp>
#include <opm/simulators/flow/flux/ParentSummary.hpp>
#include <opm/simulators/flow/FlowUtils.hpp>
#include <opm/simulators/flow/LgrOutputTransGather.hpp>
#include <opm/simulators/flow/TracerModel.hpp>
#include <opm/simulators/flow/TemperatureModel.hpp>
#include <opm/simulators/flow/Transmissibility.hpp>
#include <opm/simulators/timestepping/AdaptiveTimeStepping.hpp>
#include <opm/simulators/timestepping/SimulatorReport.hpp>

#include <opm/simulators/utils/DeferredLoggingErrorHelpers.hpp>
#include <opm/simulators/utils/ParallelSerialization.hpp>
#include <opm/simulators/utils/satfunc/RelpermDiagnostics.hpp>

#include <opm/utility/CopyablePtr.hpp>

#include <algorithm>
#include <cstddef>
#include <cmath>
#include <filesystem>
#include <functional>
#include <sstream>
#include <set>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

namespace Opm {

/*!
 * \ingroup BlackOilSimulator
 *
 * \brief This problem simulates an input file given in the data format used by the
 *        commercial ECLiPSE simulator.
 */
template <class TypeTag>
class FlowProblem : public GetPropType<TypeTag, Properties::BaseProblem>
                  , public FlowGenericProblem<GetPropType<TypeTag, Properties::GridView>,
                                              GetPropType<TypeTag, Properties::FluidSystem>>
{
protected:
    using BaseType = FlowGenericProblem<GetPropType<TypeTag, Properties::GridView>,
                                        GetPropType<TypeTag, Properties::FluidSystem>>;
    using ParentType = GetPropType<TypeTag, Properties::BaseProblem>;
    using Implementation = GetPropType<TypeTag, Properties::Problem>;

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using Stencil = GetPropType<TypeTag, Properties::Stencil>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using GlobalEqVector = GetPropType<TypeTag, Properties::GlobalEqVector>;
    using EqVector = GetPropType<TypeTag, Properties::EqVector>;
    using Vanguard = GetPropType<TypeTag, Properties::Vanguard>;
    using Indices = GetPropType<TypeTag, Properties::Indices>;

    // Grid and world dimension
    enum { dim = GridView::dimension };
    enum { dimWorld = GridView::dimensionworld };

    // copy some indices for convenience
    enum { numEq = getPropValue<TypeTag, Properties::NumEq>() };
    enum { numPhases = FluidSystem::numPhases };
    enum { numComponents = FluidSystem::numComponents };

    static constexpr bool enableBioeffects = getPropValue<TypeTag, Properties::EnableBioeffects>();
    static constexpr bool enableBrine = getPropValue<TypeTag, Properties::EnableBrine>();
    static constexpr bool enableConvectiveMixing = getPropValue<TypeTag, Properties::EnableConvectiveMixing>();
    static constexpr bool enableDiffusion = getPropValue<TypeTag, Properties::EnableDiffusion>();
    static constexpr bool enableDispersion = getPropValue<TypeTag, Properties::EnableDispersion>();
    static constexpr bool enableExtbo = getPropValue<TypeTag, Properties::EnableExtbo>();
    static constexpr bool enableFoam = getPropValue<TypeTag, Properties::EnableFoam>();
    static constexpr bool enablePolymer = getPropValue<TypeTag, Properties::EnablePolymer>();
    static constexpr bool enablePolymerMolarWeight = getPropValue<TypeTag, Properties::EnablePolymerMW>();
    static constexpr bool enableSolvent = getPropValue<TypeTag, Properties::EnableSolvent>();

    static constexpr EnergyModules energyModuleType = getPropValue<TypeTag, Properties::EnergyModuleType>();
    enum { enableFullyImplicitThermal = getPropValue<TypeTag, Properties::EnergyModuleType>() == EnergyModules::FullyImplicitThermal };
    enum { enableExperiments = getPropValue<TypeTag, Properties::EnableExperiments>() };
    enum { enableMICP = Indices::enableMICP };
    enum { enableSaltPrecipitation = getPropValue<TypeTag, Properties::EnableSaltPrecipitation>() };
    enum { enableThermalFluxBoundaries = getPropValue<TypeTag, Properties::EnableThermalFluxBoundaries>() };

    enum { gasPhaseIdx = FluidSystem::gasPhaseIdx };
    enum { oilPhaseIdx = FluidSystem::oilPhaseIdx };
    enum { waterPhaseIdx = FluidSystem::waterPhaseIdx };

    // TODO: later, gasCompIdx, oilCompIdx and waterCompIdx should go to the FlowProblemBlackoil in the future
    // we do not want them in the compositional setting
    enum { gasCompIdx = FluidSystem::gasCompIdx };
    enum { oilCompIdx = FluidSystem::oilCompIdx };
    enum { waterCompIdx = FluidSystem::waterCompIdx };

    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using RateVector = GetPropType<TypeTag, Properties::RateVector>;
    using Simulator = GetPropType<TypeTag, Properties::Simulator>;
    using Element = typename GridView::template Codim<0>::Entity;
    using ElementContext = GetPropType<TypeTag, Properties::ElementContext>;
    using EclMaterialLawManager = typename GetProp<TypeTag, Properties::MaterialLaw>::EclMaterialLawManager;
    using EclThermalLawManager = typename GetProp<TypeTag, Properties::SolidEnergyLaw>::EclThermalLawManager;
    using MaterialLawParams = typename EclMaterialLawManager::MaterialLawParams;
    using SolidEnergyLawParams = typename EclThermalLawManager::SolidEnergyLawParams;
    using ThermalConductionLawParams = typename EclThermalLawManager::ThermalConductionLawParams;
    using MaterialLaw = GetPropType<TypeTag, Properties::MaterialLaw>;
    using DofMapper = GetPropType<TypeTag, Properties::DofMapper>;
    using Evaluation = GetPropType<TypeTag, Properties::Evaluation>;
    using IntensiveQuantities = GetPropType<TypeTag, Properties::IntensiveQuantities>;
    using WellModel = GetPropType<TypeTag, Properties::WellModel>;
    using AquiferModel = GetPropType<TypeTag, Properties::AquiferModel>;

    using Toolbox = MathToolbox<Evaluation>;
    using DimMatrix = Dune::FieldMatrix<Scalar, dimWorld, dimWorld>;

    using TemperatureModel = GetPropType<TypeTag, Properties::TemperatureModel>;
    using TracerModel = GetPropType<TypeTag, Properties::TracerModel>;
    using DirectionalMobilityPtr = Utility::CopyablePtr<DirectionalMobility<TypeTag>>;
    using ModuleParams = BlackoilModuleParams<ConvectiveMixingModuleParam<Scalar>>;

public:

    /// Extra module-specific lines appended to the final simulation trailer
    /// (e.g. geomechanics timing). Derived problems may shadow this.
    std::string extraTrailerSummary() const
    { return {}; }
    using BaseType::briefDescription;
    using BaseType::helpPreamble;
    using BaseType::shouldWriteOutput;
    using BaseType::shouldWriteRestartFile;
    using BaseType::rockBiotComp;
    using BaseType::rockBiotTemp;
    using BaseType::lame;
    using BaseType::biotCoeff;
    using BaseType::biotTemp;
    using BaseType::rockCompressibility;
    using BaseType::porosity;

    /*!
     * \copydoc FvBaseProblem::registerParameters
     */
    static void registerParameters()
    {
        ParentType::registerParameters();

        registerFlowProblemParameters<Scalar>();
    }

    /*!
     * \copydoc FvBaseProblem::handlePositionalParameter
     *
     * \param addKey Callback for adding a key-value pair
     * \param seenParams Already seen parameters
     * \param errorMsg Error message
     * \param argv Command line parameters
     * \param paramIdx Index of parameter to handle
     */
    static int handlePositionalParameter(std::function<void(const std::string&,
                                                            const std::string&)> addKey,
                                         std::set<std::string>& seenParams,
                                         std::string& errorMsg,
                                         int,
                                         const char** argv,
                                         int paramIdx,
                                         int)
    {
        return detail::eclPositionalParameter(addKey,
                                              seenParams,
                                              errorMsg,
                                              argv,
                                              paramIdx);
    }

    /*!
     * \copydoc Doxygen::defaultProblemConstructor
     */
    explicit FlowProblem(Simulator& simulator)
        : ParentType(simulator)
        , BaseType(simulator.vanguard().eclState(),
                   simulator.vanguard().schedule(),
                   simulator.vanguard().gridView())
        , transmissibilities_(simulator.vanguard().eclState(),
                              simulator.vanguard().gridView(),
                              simulator.vanguard().cartesianIndexMapper(),
                              simulator.vanguard().grid(),
                              simulator.vanguard().cellCentroids(),
                              (energyModuleType == EnergyModules::FullyImplicitThermal ||
                               energyModuleType == EnergyModules::SequentialImplicitThermal),
                              enableDiffusion,
                              enableDispersion)
        , wellModel_(simulator, this->iterationContext())
        , aquiferModel_(simulator)
        , pffDofData_(simulator.gridView(), this->elementMapper())
        , tracerModel_(simulator)
        , temperatureModel_(simulator)
        , thresholdPressures_(simulator)
        , enable_state_rollback_(Parameters::Get<Parameters::EnableStateRollback>())
    {
        if (! Parameters::Get<Parameters::CheckSatfuncConsistency>()) {
            // User did not enable the "new" saturation function consistency
            // check module.  Run the original checker instead.  This is a
            // temporary measure.
            RelpermDiagnostics relpermDiagnostics{};
            relpermDiagnostics.diagnosis(simulator.vanguard().eclState(),
                                         simulator.vanguard().levelCartesianIndexMapper());
        }

        if (energyModuleType == EnergyModules::SequentialImplicitThermal) {
            this->enableDriftCompensationTemp_ = Parameters::Get<Parameters::EnableDriftCompensationTemp>();
        }

    }

    virtual ~FlowProblem() = default;

    void prefetch(const Element& elem) const
    { this->pffDofData_.prefetch(elem); }

    /*!
     * \brief This method restores the complete state of the problem and its sub-objects
     *        from disk.
     *
     * The serialization format used by this method is ad-hoc. It is the inverse of the
     * serialize() method.
     *
     * \tparam Restarter The deserializer type
     *
     * \param res The deserializer object
     */
    template <class Restarter>
    void deserialize(Restarter& res)
    {
        // reload the current episode/report step from the deck
        this->beginEpisode();

        // deserialize the wells
        wellModel_.deserialize(res);

        // deserialize the aquifer
        aquiferModel_.deserialize(res);
    }

    /*!
     * \brief This method writes the complete state of the problem and its subobjects to
     *        disk.
     *
     * The file format used here is ad-hoc.
     */
    template <class Restarter>
    void serialize(Restarter& res)
    {
        wellModel_.serialize(res);

        aquiferModel_.serialize(res);
    }

    int episodeIndex() const
    {
        return std::max(this->simulator().episodeIndex(), 0);
    }

    /*!
     * \brief Called by the simulator before an episode begins.
     */
    virtual void beginEpisode()
    {
        OPM_TIMEBLOCK(beginEpisode);
        // Proceed to the next report step
        auto& simulator = this->simulator();
        int episodeIdx = simulator.episodeIndex();
        auto& eclState = simulator.vanguard().eclState();
        const auto& schedule = simulator.vanguard().schedule();
        const auto& events = schedule[episodeIdx].events();

        if (episodeIdx >= 0 && events.hasEvent(ScheduleEvents::GEO_MODIFIER)) {
            // bring the contents of the keywords to the current state of the SCHEDULE
            // section.
            //
            // TODO (?): make grid topology changes possible (depending on what exactly
            // has changed, the grid may need be re-created which has some serious
            // implications on e.g., the solution of the simulation.)
            const auto& miniDeck = schedule[episodeIdx].geo_keywords();
            const auto& cc = simulator.vanguard().grid().comm();
            eclState.apply_schedule_keywords( miniDeck );
            eclBroadcast(cc, eclState.getTransMult() );

            // Re-ordering in case of ALUGrid
            std::function<unsigned int(unsigned int)> equilGridToGrid = [&simulator](unsigned int i) {
                  return simulator.vanguard().gridEquilIdxToGridIdx(i);
            };

            // re-compute all quantities which may possibly be affected.
            using TransUpdateQuantities = typename Vanguard::TransmissibilityType::TransUpdateQuantities;
            transmissibilities_.update(true, TransUpdateQuantities::All, equilGridToGrid);
            this->referencePorosity_[1] = this->referencePorosity_[0];
            updateReferencePorosity_();
            this->rockFraction_[1] = this->rockFraction_[0];
            updateRockFraction_();
            updatePffDofData_();
            this->applyFluxBoundaryTransmissibilityOverrides_();
            this->model().linearizer().updateDiscretizationParameters();
        }

        bool tuningEvent = this->beginEpisode_(enableExperiments, this->episodeIndex());

        // set up the wells for the next episode.
        wellModel_.beginEpisode();

        // set up the aquifers for the next episode.
        aquiferModel_.beginEpisode();

        // set the size of the initial time step of the episode
        Scalar dt = limitNextTimeStepSize_(simulator.episodeLength());
        // negative value of initialTimeStepSize_ indicates no active limit from TSINIT or NEXTSTEP
        if ( (episodeIdx == 0 || tuningEvent) && this->initialTimeStepSize_ > 0)
            // allow the size of the initial time step to be set via an external parameter
            // if TUNING is enabled, also limit the time step size after a tuning event to TSINIT
            dt = std::min(dt, this->initialTimeStepSize_);
        simulator.setTimeStepSize(dt);
    }

    /*!
     * \brief Called by the simulator before each time integration.
     */
    virtual void beginTimeStep()
    {
        OPM_TIMEBLOCK(beginTimeStep);
        const int episodeIdx = this->episodeIndex();
        const int timeStepSize = this->simulator().timeStepSize();

        if (enable_state_rollback_) {
            this->captureBeginTimeStepState_();
        }

        // Any auxiliary cell still without a start-of-step accumulation gets one, before
        // anything assembles.  A cell not yet handed out never earns one on its own: the
        // linearizer fills the entry from the previous step's, and a cell that has been
        // sitting idle has no previous step to have moved in.
        this->seedAuxCellStorageCache_();

        this->beginTimeStep_(enableExperiments,
                             episodeIdx,
                             this->simulator().timeStepIndex(),
                             this->simulator().startTime(),
                             this->simulator().time(),
                             timeStepSize,
                             this->simulator().endTime());

        // update maximum water saturation and minimum pressure
        // used when ROCKCOMP is activated
        // Do not update max RS first step after a restart
        this->updateExplicitQuantities_(episodeIdx, timeStepSize, first_step_ && (episodeIdx > 0));
        first_step_ = false;

        if (nonTrivialBoundaryConditions()) {
            this->model().linearizer().updateBoundaryConditionData();
        }

        wellModel_.beginTimeStep();
        aquiferModel_.beginTimeStep();
        tracerModel_.beginTimeStep();
        temperatureModel_.beginTimeStep();
    }

    /*!
     * \brief Called by the simulator to restore the state captured at the
     *        beginning of the timestep after a failed timestep.
     */
    void updateFailed()
    {
        if (enable_state_rollback_) {
            this->restoreBeginTimeStepState_();
        }
        wellModel_.updateFailed();
        this->model().updateFailed();
    }

    /*!
     * \brief Called by the simulator to accept the current state as the new
     *        time level after a successful timestep.
     */
    void advanceTimeLevel()
    {
        this->model().advanceTimeLevel();
        wellModel_.advanceTimeLevel();
    }

    /*!
     * \brief Called by the simulator before each Newton-Raphson iteration.
     */
    void beginIteration()
    {
        OPM_TIMEBLOCK(beginIteration);
        wellModel_.beginIteration();
        aquiferModel_.beginIteration();
    }

    /*!
     * \brief Called by the simulator after each Newton-Raphson iteration.
     */
    void endIteration()
    {
        OPM_TIMEBLOCK(endIteration);
        wellModel_.endIteration();
        aquiferModel_.endIteration();
    }

    /*!
     * \brief Called by the simulator after each time integration.
     */
    virtual void endTimeStep()
    {
        OPM_TIMEBLOCK(endTimeStep);

#ifndef NDEBUG
        if constexpr (getPropValue<TypeTag, Properties::EnableDebuggingChecks>()) {
            // in debug mode, we don't care about performance, so we check
            // if the model does the right thing (i.e., the mass change
            // inside the whole reservoir must be equivalent to the fluxes
            // over the grid's boundaries plus the source rates specified by
            // the problem).
            const int rank = this->simulator().gridView().comm().rank();
            if (rank == 0) {
                std::cout << "checking conservativeness of solution\n";
            }

            this->model().checkConservativeness(/*tolerance=*/-1, /*verbose=*/true);
            if (rank == 0) {
                std::cout << "solution is sufficiently conservative\n";
            }
        }
#endif // NDEBUG

        auto& simulator = this->simulator();
        simulator.setTimeStepIndex(simulator.timeStepIndex()+1);

        this->wellModel_.endTimeStep();
        this->aquiferModel_.endTimeStep();
        this->tracerModel_.endTimeStep();

        // Compute flux for output
        this->model().linearizer().updateFlowsInfo();

        if (this->enableDriftCompensation_ || this->enableDriftCompensationTemp_) {
            OPM_TIMEBLOCK(driftCompansation);

            const auto& residual = this->model().linearizer().residual();

            for (unsigned globalDofIdx = 0; globalDofIdx < residual.size(); globalDofIdx ++) {
                int sfcdofIdx = simulator.vanguard().gridEquilIdxToGridIdx(globalDofIdx);
                this->drift_[sfcdofIdx] = residual[sfcdofIdx] * simulator.timeStepSize();

                if constexpr (getPropValue<TypeTag, Properties::UseVolumetricResidual>()) {
                    this->drift_[sfcdofIdx] *= this->model().dofTotalVolume(sfcdofIdx);
                }
            }
        }

        // For sequential implicit thermal: finalize TEMP solve,
        // then refresh cached intensive quantities (temperature-dependent).
        if constexpr(energyModuleType == EnergyModules::SequentialImplicitThermal) {
            this->temperatureModel_.endTimeStep(wellModel_.wellState());
            this->model().invalidateAndUpdateIntensiveQuantities(/*timeIdx=*/0);
        }
    }

    /*!
     * \brief Called by the simulator after the end of an episode.
     */
    virtual void endEpisode()
    {
        const int episodeIdx = this->episodeIndex();

        this->wellModel_.endEpisode();
        this->aquiferModel_.endEpisode();

        const auto& schedule = this->simulator().vanguard().schedule();

        // End simulation when completed.
        if (episodeIdx + 1 >= static_cast<int>(schedule.size()) - 1) {
            this->simulator().setFinished(true);
            return;
        }

        // Otherwise, start next episode (report step).
        this->simulator().startNextEpisode(schedule.stepLength(episodeIdx + 1));
    }

    /*!
     * \brief Write the requested quantities of the current solution into the output
     *        files.
     */
    virtual void writeOutput(bool verbose)
    {
        OPM_TIMEBLOCK(problemWriteOutput);

        if (Parameters::Get<Parameters::EnableWriteAllSolutions>() ||
            this->episodeWillBeOver())
        {
            // Create VTK output as needed.
            ParentType::writeOutput(verbose);
        }
    }

    /*!
     * \copydoc FvBaseMultiPhaseProblem::intrinsicPermeability
     */
    template <class Context>
    const DimMatrix& intrinsicPermeability(const Context& context,
                                           unsigned spaceIdx,
                                           unsigned timeIdx) const
    {
        unsigned globalSpaceIdx = context.globalSpaceIndex(spaceIdx, timeIdx);
        return transmissibilities_.permeability(globalSpaceIdx);
    }

    /*!
     * \brief This method returns the intrinsic permeability tensor
     *        given a global element index.
     *
     * Its main (only?) usage is the ECL transmissibility calculation code...
     */
    const DimMatrix& intrinsicPermeability(unsigned globalElemIdx) const
    { return transmissibilities_.permeability(globalElemIdx); }

    /*!
     * \copydoc EclTransmissiblity::transmissibility
     */
    template <class Context>
    Scalar transmissibility(const Context& context,
                            [[maybe_unused]] unsigned fromDofLocalIdx,
                            unsigned toDofLocalIdx) const
    {
        assert(fromDofLocalIdx == 0);
        return pffDofData_.get(context.element(), toDofLocalIdx).transmissibility;
    }

    /*!
     * \brief Direct access to the transmissibility between two elements.
     */
    Scalar transmissibility(unsigned globalCenterElemIdx, unsigned globalElemIdx) const
    {
        return transmissibilities_.transmissibility(globalCenterElemIdx, globalElemIdx);
    }

    /*!
     * \brief Threshold pressure [Pa] for the intersection between two elements.
     */
    Scalar thresholdPressure(unsigned elem1Idx, unsigned elem2Idx) const
    {
        // An auxiliary cell is in no equilibration region or fault; the tables are grid-sized.
        const auto numGridDof = this->model().numGridDof();
        if (elem1Idx >= numGridDof || elem2Idx >= numGridDof) {
            return 0.0;
        }
        return thresholdPressures_.thresholdPressure(elem1Idx, elem2Idx);
    }

    const FlowThresholdPressure<TypeTag>& thresholdPressure() const
    { return thresholdPressures_; }

    FlowThresholdPressure<TypeTag>& thresholdPressure()
    { return thresholdPressures_; }

    const ModuleParams& moduleParams() const
    { return moduleParams_; }

    /*!
     * \copydoc Transmissibility::transmissibilityOrZero
     */
    Scalar transmissibilityOrZero(unsigned globalCenterElemIdx, unsigned globalElemIdx) const
    {
        return transmissibilities_.transmissibilityOrZero(globalCenterElemIdx, globalElemIdx);
    }

    /*!
     * \copydoc EclTransmissiblity::diffusivity
     */
    template <class Context>
    Scalar diffusivity(const Context& context,
                       [[maybe_unused]] unsigned fromDofLocalIdx,
                       unsigned toDofLocalIdx) const
    {
        assert(fromDofLocalIdx == 0);
        return *pffDofData_.get(context.element(), toDofLocalIdx).diffusivity;
    }

    /*!
     * give the transmissibility for a face i.e. pair. should be symmetric?
     */
    Scalar diffusivity(const unsigned globalCellIn, const unsigned globalCellOut) const{
        return transmissibilities_.diffusivity(globalCellIn, globalCellOut);
    }

    /*!
     * give the dispersivity for a face i.e. pair.
     */
    Scalar dispersivity(const unsigned globalCellIn, const unsigned globalCellOut) const{
        return transmissibilities_.dispersivity(globalCellIn, globalCellOut);
    }

    /*!
     * \brief Direct access to a boundary transmissibility.
     */
    Scalar thermalTransmissibilityBoundary(const unsigned globalSpaceIdx,
                                    const unsigned boundaryFaceIdx) const
    {
        return transmissibilities_.thermalTransmissibilityBoundary(globalSpaceIdx, boundaryFaceIdx);
    }




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
     * \brief Pore-volume weighted sums for the cells a USEFLUX run does not
     *        have, as recorded by the run that produced its boundary data.
     *
     * \details Returns false when there are none, which is every run that is
     *   not reading a .FLUX file. See FluxFile::ReportStep::externalRegionSums
     *   for the layout.
     */
    bool fluxConverterExternalSums(std::array<Scalar, 8>& hydrocarbonPvWeighted,
                                   std::array<Scalar, 8>& poreVolumeWeighted) const
    {
        const auto* fluxStep = this->fluxBoundaryReportStep_();
        if (fluxStep == nullptr || fluxStep->externalRegionSums.size() != 16) {
            return false;
        }

        for (std::size_t i = 0; i < 8; ++i) {
            hydrocarbonPvWeighted[i] = static_cast<Scalar>(fluxStep->externalRegionSums[i]);
            poreVolumeWeighted[i] = static_cast<Scalar>(fluxStep->externalRegionSums[i + 8]);
        }

        return true;
    }

    /*!
     * \brief Threshold pressure across a boundary face.
     *
     * \details Zero unless the run is reading a sector boundary out of a .FLUX
     *   file and the producing run had threshold pressures across it. The
     *   ordering matters when THPRES is irreversible, hence the direction flag:
     *   \p interiorToExterior selects the threshold resisting flow out of this
     *   cell.
     */
    Scalar thresholdPressureBoundary(unsigned /*globalSpaceIdx*/,
                                     unsigned /*boundaryFaceIdx*/,
                                     bool /*interiorToExterior*/) const
    { return 0.0; }

    /*!
     * \brief Relative permeabilities of the fluid entering through a boundary
     *        face.
     *
     * \details A USEFLUX run takes them from the .FLUX file, where the
     *   producing run recorded the exterior cell's values, so that the
     *   saturation-function region, the scaled end points and the hysteresis
     *   state of the cell the flow comes from are all honoured. Everything else
     *   falls back to evaluating this cell's saturation functions at the
     *   exterior saturations, which is what the code did unconditionally
     *   before.
     */
    template <class FluidState, class Container>
    void boundaryRelativePermeabilities(Container& kr,
                                        const unsigned globalSpaceIdx,
                                        const unsigned boundaryFaceIdx,
                                        const FluidState& exFluidState) const
    {
        // boundaryFaceIdx is the face's POSITION in this cell's boundary-face
        // list, which is how the discretisation addresses it, not a direction.
        // The two only agree while a cell's boundary faces occupy the leading
        // direction slots, and a cell can have more than six of them: a
        // corner-point face split by a fault contributes several, and in
        // parallel every intersection that reaches off this rank counts too.
        // Reading it as a direction then asks for a seventh one.
        const auto dir = this->fluxBoundaryDirAt_(globalSpaceIdx, boundaryFaceIdx);
        const auto* fluxData = this->fluxBoundaryData_();
        const auto* fluxStep = this->fluxBoundaryReportStep_();
        const auto fluxSlot = (dir == FaceDir::Unknown)
            ? 0
            : this->fluxBoundaryFaceSlot_(globalSpaceIdx, dir);

        if (fluxData != nullptr && fluxStep != nullptr && fluxSlot > 0) {
            const auto numPhaseSlots = static_cast<std::size_t>(fluxData->header.numPhases);
            const auto base = static_cast<std::size_t>(fluxSlot - 1) * numPhaseSlots;

            if (fluxStep->relPerm.size() >= base + numPhaseSlots) {
                for (unsigned phaseIdx = 0; phaseIdx < FluidSystem::numPhases; ++phaseIdx) {
                    kr[phaseIdx] = 0.0;
                }

                std::size_t slot = 0;
                const auto load = [&](const unsigned phaseIdx)
                {
                    if (FluidSystem::phaseIsActive(phaseIdx)) {
                        kr[phaseIdx] = fluxStep->relPerm[base + slot];
                        ++slot;
                    }
                };

                load(oilPhaseIdx);
                load(waterPhaseIdx);
                load(gasPhaseIdx);
                return;
            }
        }

        MaterialLaw::relativePermeabilities(kr, this->materialLawParams(globalSpaceIdx), exFluidState);
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


    const TracerModel& tracerModel() const
    { return tracerModel_; }

    TracerModel& tracerModel()
    { return tracerModel_; }

    TemperatureModel& temperatureModel() // need for restart
    { return temperatureModel_; }

    /*!
     * \copydoc FvBaseMultiPhaseProblem::porosity
     *
     * For the FlowProblem, this method is identical to referencePorosity(). The intensive
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
        if (globalSpaceIdx >= this->model().numGridDof()) {
            return this->auxCellDepth_(globalSpaceIdx);
        }

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
     * \copydoc BlackoilProblem::rockBiotComp
     */
    template <class Context>
    Scalar rockBiotComp(const Context& context, unsigned spaceIdx, unsigned timeIdx) const
    {
        unsigned globalSpaceIdx = context.globalSpaceIndex(spaceIdx, timeIdx);
        return this->rockBiotComp(globalSpaceIdx);
    }

    /*!
     * \copydoc BlackoilProblem::rockBiotTemp
     */
    template <class Context>
    Scalar rockBiotTemp(const Context& context, unsigned spaceIdx, unsigned timeIdx) const
    {
        unsigned globalSpaceIdx = context.globalSpaceIndex(spaceIdx, timeIdx);
        return this->rockBiotTemp(globalSpaceIdx);
    }

    /*!
     * \copydoc BlackoilProblem::lame
     */
    template <class Context>
    Scalar lame(const Context& context, unsigned spaceIdx, unsigned timeIdx) const
    {
        unsigned globalSpaceIdx = context.globalSpaceIndex(spaceIdx, timeIdx);
        return this->lame(globalSpaceIdx);
    }

    /*!
     * \copydoc BlackoilProblem::biotCoeff
     */
    template <class Context>
    Scalar biotCoeff(const Context& context, unsigned spaceIdx, unsigned timeIdx) const
    {
        unsigned globalSpaceIdx = context.globalSpaceIndex(spaceIdx, timeIdx);
        return this->biotCoeff(globalSpaceIdx);
    }

    /*!
     * \copydoc BlackoilProblem::biotTemp
     */
    template <class Context>
    Scalar biotTemp(const Context& context, unsigned spaceIdx, unsigned timeIdx) const
    {
        unsigned globalSpaceIdx = context.globalSpaceIndex(spaceIdx, timeIdx);
        return this->biotTemp(globalSpaceIdx);
    }

    /*!
     * \copydoc BlackoilProblem::rockReferencePressure
     */
    template <class Context>
    Scalar rockReferencePressure(const Context& context, unsigned spaceIdx, unsigned timeIdx) const
    {
        unsigned globalSpaceIdx = context.globalSpaceIndex(spaceIdx, timeIdx);
        return rockReferencePressure(globalSpaceIdx);
    }

    /*!
    * \copydoc BlackoilProblem::rockReferencePressure
    */
    Scalar rockReferencePressure(unsigned globalSpaceIdx) const
    {
        const auto& rock_config = this->simulator().vanguard().eclState().getSimulationConfig().rock_config();
        if (rock_config.store()) {
            return asImp_().initialFluidState(globalSpaceIdx).pressure(refPressurePhaseIdx_());
        }
        else {
            if (this->rockParams_.empty())
                return 1e5;

            unsigned tableIdx = 0;
            if (!this->rockTableIdx_.empty()) {
                // Auxiliary DOFs have no entry here; they take the first ROCK region.
                if (globalSpaceIdx < this->rockTableIdx_.size()) {
                    tableIdx = this->rockTableIdx_[globalSpaceIdx];
                }
            }
            return this->rockParams_[tableIdx].referencePressure;
        }
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
        return materialLawManager_->materialLawParams(auxCellSaturationProxy_(globalDofIdx));
    }

    //! Grid cell whose saturation functions a DOF uses; identity for grid cells.
    unsigned saturationFunctionCell(unsigned globalDofIdx) const
    { return auxCellSaturationProxy_(globalDofIdx); }

    const MaterialLawParams& materialLawParams(unsigned globalDofIdx, FaceDir::DirEnum facedir) const
    {
        return materialLawManager_->materialLawParams(auxCellSaturationProxy_(globalDofIdx), facedir);
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
    const EclMaterialLawManager& materialLawManager() const
    { return *materialLawManager_; }

    std::shared_ptr<const EclThermalLawManager> thermalLawManager() const
    { return thermalLawManager_; }

    template <class FluidState, class ...Args>
    void updateRelperms(
        std::array<Evaluation,numPhases> &mobility,
        DirectionalMobilityPtr &dirMob,
        FluidState &fluidState,
        unsigned globalSpaceIdx) const
    {
        using ContainerT = std::array<Evaluation, numPhases>;
        OPM_TIMEBLOCK_LOCAL(updateRelperms, Subsystem::SatProps);
        {
            // calculate relative permeabilities. note that we store the result into the
            // mobility_ class attribute. the division by the phase viscosity happens later.
            const auto& materialParams = materialLawParams(globalSpaceIdx);
            MaterialLaw::template relativePermeabilities<ContainerT, FluidState, Args...>(mobility, materialParams, fluidState);
            Valgrind::CheckDefined(mobility);
        }
        if (materialLawManager_->hasDirectionalRelperms()
               || materialLawManager_->hasDirectionalImbnum())
        {
            using Dir = FaceDir::DirEnum;
            constexpr int ndim = 3;
            dirMob = std::make_unique<DirectionalMobility<TypeTag>>();
            Dir facedirs[ndim] = {Dir::XPlus, Dir::YPlus, Dir::ZPlus};
            for (int i = 0; i<ndim; i++) {
                const auto& materialParams = materialLawParams(globalSpaceIdx, facedirs[i]);
                auto& mob_array = dirMob->getArray(i);
                MaterialLaw::template relativePermeabilities<ContainerT, FluidState, Args...>(mob_array, materialParams, fluidState);
            }
        }
    }

    /*!
     * \copydoc materialLawManager()
     */
    EclMaterialLawManager& materialLawManager()
    { return *materialLawManager_; }

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

    // TODO: polymer related might need to go to the blackoil side
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
        if constexpr (energyModuleType == EnergyModules::SequentialImplicitThermal)
            return temperatureModel_.temperature(globalDofIdx);

        return asImp_().initialFluidState(globalDofIdx).temperature(/*phaseIdx=*/0);
    }


    Scalar temperature(unsigned globalDofIdx, unsigned /*timeIdx*/) const
    {
        // use the initial temperature of the DOF if temperature is not a primary
        // variable
        if constexpr (energyModuleType == EnergyModules::SequentialImplicitThermal)
            return temperatureModel_.temperature(globalDofIdx);

        return asImp_().initialFluidState(globalDofIdx).temperature(/*phaseIdx=*/0);
    }

    const SolidEnergyLawParams&
    solidEnergyLawParams(unsigned globalSpaceIdx,
                         unsigned /*timeIdx*/) const
    {
        return this->thermalLawManager_->solidEnergyLawParams(auxCellSaturationProxy_(globalSpaceIdx));
    }
    const ThermalConductionLawParams &
    thermalConductionLawParams(unsigned globalSpaceIdx,
                               unsigned /*timeIdx*/)const
    {
        return this->thermalLawManager_->thermalConductionLawParams(auxCellSaturationProxy_(globalSpaceIdx));
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
     * \copydoc FvBaseProblem::initialSolutionApplied()
     */
    virtual void initialSolutionApplied()
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

        // An auxiliary cell goes through the ordinary assembly whether or not it is in
        // use, and that reaches for the accumulation of the previous time level.  The
        // linearizer fills it on the first Newton iteration of a step, which is soon
        // enough for a cell that is there from the start, but not for one still waiting to
        // be handed out -- nor for one handed out after that iteration has gone by, which
        // is when a fracture opens its cells.  Seed them all once here so that no
        // auxiliary row can reach for an entry that was never written.
        this->seedAuxCellStorageCache_();

        // initialize the wells. Note that this needs to be done after initializing the
        // intrinsic permeabilities and the after applying the initial solution because
        // the well model uses these...
        wellModel_.init();

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
        OPM_TIMEBLOCK_LOCAL(eclProblemSource, Subsystem::Assembly);
        rate = 0.0;

        // Add well contribution to source here.
        wellModel_.computeTotalRatesForDof(rate, globalDofIdx);

        // A degree of freedom with no volume holds nothing and so can be the source of
        // nothing -- an auxiliary cell that has been reserved but not yet put to use.
        // Dividing by its volume below would turn a zero rate into a NaN.
        const auto volume = this->model().dofTotalVolume(globalDofIdx);
        if (volume <= 0.0) {
            return;
        }

        // convert the source term from the total mass rate of the
        // cell to the one per unit of volume as used by the model.
        for (unsigned eqIdx = 0; eqIdx < numEq; ++ eqIdx) {
            rate[eqIdx] /= volume;

            Valgrind::CheckDefined(rate[eqIdx]);
            assert(isfinite(rate[eqIdx]));
        }

        // Add non-well sources.
        addToSourceDense(rate, globalDofIdx, timeIdx);
    }

    virtual void addToSourceDense(RateVector& rate,
                                  unsigned globalDofIdx,
                                  unsigned timeIdx) const = 0;

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

    bool nonTrivialBoundaryConditions() const
    { return nonTrivialBoundaryConditions_; }

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
        OPM_TIMEBLOCK_LOCAL(rockCompPoroMultiplier, Subsystem::PvtProps);
        if (this->rockCompPoroMult_.empty() && this->rockCompPoroMultWc_.empty())
            return 1.0;

        unsigned tableIdx = 0;
        if (!this->rockTableIdx_.empty())
            // Auxiliary DOFs have no entry here; they take the first ROCK region.
            if (elementIdx < this->rockTableIdx_.size()) {
                tableIdx = this->rockTableIdx_[elementIdx];
            }

        const auto& fs = intQuants.fluidState();
        const auto& rock_config = this->simulator().vanguard().eclState().getSimulationConfig().rock_config();

        if (!this->rockCompPoroMultElastic_.empty()) {
            // ROCKCOMP HYSTERESIS=HYSTER: follow the deflation (virgin/plastic) curve
            // at or below the lowest pressure ever reached by this cell, and the
            // reversible elastic reload curve, anchored at that turning pressure,
            // above it.
            LhsEval effectivePressure = decay<LhsEval>(fs.pressure(refPressurePhaseIdx_()));
            LhsEval turningPressure = this->minRefPressure_[elementIdx];

            if (!this->overburdenPressure_.empty()) {
                effectivePressure -= this->overburdenPressure_[elementIdx];
                turningPressure -= this->overburdenPressure_[elementIdx];
            }

            if (rock_config.store()) {
                const auto& initialPressure = asImp_().initialFluidState(elementIdx).pressure(refPressurePhaseIdx_());
                effectivePressure -= initialPressure;
                turningPressure -= initialPressure;
            }

            if (effectivePressure <= turningPressure)
                return this->rockCompPoroMult_[tableIdx].eval(effectivePressure, /*extrapolation=*/true);

            return this->rockCompPoroMultElastic_[tableIdx].eval(turningPressure, effectivePressure, /*extrapolation=*/true);
        }

        LhsEval effectivePressure = decay<LhsEval>(fs.pressure(refPressurePhaseIdx_()));
        if (!this->minRefPressure_.empty())
            // The pore space change is irreversible
            effectivePressure =
                min(decay<LhsEval>(fs.pressure(refPressurePhaseIdx_())),
                                   this->minRefPressure_[elementIdx]);

        if (!this->overburdenPressure_.empty())
            effectivePressure -= this->overburdenPressure_[elementIdx];

        if (rock_config.store()) {
            effectivePressure -= asImp_().initialFluidState(elementIdx).pressure(refPressurePhaseIdx_());
        }

        if (!this->rockCompPoroMult_.empty()) {
            return this->rockCompPoroMult_[tableIdx].eval(effectivePressure, /*extrapolation=*/true);
        }

        // water compaction
        assert(!this->rockCompPoroMultWc_.empty());
        LhsEval SwMax = max(decay<LhsEval>(fs.saturation(waterPhaseIdx)), this->maxWaterSaturation_[elementIdx]);
        LhsEval SwDeltaMax = SwMax - asImp_().initialFluidStates()[elementIdx].saturation(waterPhaseIdx);

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
        auto obtain = [](const auto& value)
                      {
                          if constexpr (std::is_same_v<LhsEval, Scalar>) {
                              return getValue(value);
                          } else {
                              return value;
                          }
                      };
        return rockCompTransMultiplier<LhsEval>(intQuants, elementIdx, obtain);
    }

    template <class LhsEval, class Callback>
    LhsEval rockCompTransMultiplier(const IntensiveQuantities& intQuants, unsigned elementIdx, Callback& obtain) const
    {
        const bool implicit = !this->explicitRockCompaction_;
        return implicit ? this->simulator().problem().template computeRockCompTransMultiplier_<LhsEval>(intQuants, elementIdx, obtain)
                        : this->simulator().problem().getRockCompTransMultVal(elementIdx);
    }

    template <class LhsEval, class Callback>
    LhsEval wellTransMultiplier(const IntensiveQuantities& intQuants, unsigned elementIdx, Callback& obtain) const
    {
        OPM_TIMEBLOCK_LOCAL(wellTransMultiplier, Subsystem::Wells);

        const bool implicit = !this->explicitRockCompaction_;
        LhsEval trans_mult = implicit ? this->simulator().problem().template computeRockCompTransMultiplier_<LhsEval>(intQuants, elementIdx, obtain)
                                      : this->simulator().problem().getRockCompTransMultVal(elementIdx);
        trans_mult *= this->simulator().problem().template permFactTransMultiplier<LhsEval>(intQuants, elementIdx, obtain);

        return trans_mult;
    }

    std::pair<BCType, RateVector> boundaryCondition(const unsigned int globalSpaceIdx, const int directionId) const
    {
        OPM_TIMEBLOCK_LOCAL(boundaryCondition, Subsystem::Assembly);
        if (!nonTrivialBoundaryConditions_) {
            return { BCType::NONE, RateVector(0.0) };
        }
        const FaceDir::DirEnum dir = FaceDir::FromIntersectionIndex(directionId);

        const auto& bcDirIndex = bcindex_(dir);
        const bool hasDeckBc = globalSpaceIdx < bcDirIndex.size() && bcDirIndex[globalSpaceIdx] != 0;
        if (!hasDeckBc) {
            RateVector phaseVolRate = 0.0;
            bool isMassRate = false;
            if (this->fluxBoundaryRate_(globalSpaceIdx, dir, phaseVolRate, isMassRate)) {
                if (isMassRate) {
                    // The producing run already formed the mass rate, using the
                    // density of the cell the flow comes from. Nothing to do.
                    return { BCType::RATE, phaseVolRate };
                }

                // The FLUX payload holds reservoir-condition volumetric rates, so
                // the conversion to a mass rate has to use the in-situ density.
                // Using the surface reference density would overstate the
                // exchanged mass by the formation volume factor.
                const auto pvtRegionIdx = this->pvtRegionIndex(globalSpaceIdx);
                const auto* intQuants =
                    this->simulator().model().cachedIntensiveQuantities(globalSpaceIdx, /*timeIdx=*/0);

                auto phaseDensity = [&intQuants, pvtRegionIdx](const unsigned phaseIdx)
                {
                    if (intQuants != nullptr) {
                        const auto rho = getValue(intQuants->fluidState().density(phaseIdx));
                        if (std::isfinite(rho) && rho > 0.0) {
                            return static_cast<Scalar>(rho);
                        }
                    }

                    return FluidSystem::referenceDensity(phaseIdx, pvtRegionIdx);
                };

                RateVector rate = 0.0;
                if (FluidSystem::phaseIsActive(oilPhaseIdx)) {
                    rate[FluidSystem::canonicalToActiveCompIdx(oilCompIdx)] =
                        phaseVolRate[FluidSystem::canonicalToActiveCompIdx(oilCompIdx)]
                        * phaseDensity(oilPhaseIdx);
                }
                if (FluidSystem::phaseIsActive(waterPhaseIdx)) {
                    rate[FluidSystem::canonicalToActiveCompIdx(waterCompIdx)] =
                        phaseVolRate[FluidSystem::canonicalToActiveCompIdx(waterCompIdx)]
                        * phaseDensity(waterPhaseIdx);
                }
                if (FluidSystem::phaseIsActive(gasPhaseIdx)) {
                    rate[FluidSystem::canonicalToActiveCompIdx(gasCompIdx)] =
                        phaseVolRate[FluidSystem::canonicalToActiveCompIdx(gasCompIdx)]
                        * phaseDensity(gasPhaseIdx);
                }

                return { BCType::RATE, rate };
            }

            const auto* fluxData = this->fluxBoundaryData_();
            const auto mode = fluxData ? fluxData->header.mode : EclIO::FluxFile::Mode::Flux;
            const auto pressureEnabled = (static_cast<int>(mode) & static_cast<int>(EclIO::FluxFile::Mode::Pressure)) != 0;

            // Only the faces that the FLUX file actually describes may act as a
            // pressure boundary. The remaining outer faces of the sector are
            // ordinary no-flow boundaries; treating them as Dirichlet would
            // surround the region with a constant-pressure source taken from the
            // initial state.
            if (pressureEnabled && this->fluxBoundaryFaceSlot_(globalSpaceIdx, dir) > 0) {
                return { BCType::FREE, RateVector(0.0) };
            }

            return { BCType::NONE, RateVector(0.0) };
        }

        const auto& schedule = this->simulator().vanguard().schedule();
        if (schedule[this->episodeIndex()].bcstate.size() == 0) {
            return { BCType::NONE, RateVector(0.0) };
        }
        const auto& bc = schedule[this->episodeIndex()].bcstate[bcDirIndex[globalSpaceIdx]];
        if (bc.bctype!=BCType::RATE) {
            return { bc.bctype, RateVector(0.0) };
        }

        RateVector rate = 0.0;
        switch (bc.component) {
        case BCComponent::OIL:
            rate[FluidSystem::canonicalToActiveCompIdx(oilCompIdx)] = bc.rate;
            break;
        case BCComponent::GAS:
            rate[FluidSystem::canonicalToActiveCompIdx(gasCompIdx)] = bc.rate;
            break;
        case BCComponent::WATER:
            rate[FluidSystem::canonicalToActiveCompIdx(waterCompIdx)] = bc.rate;
            break;
        case BCComponent::SOLVENT:
            this->handleSolventBC(bc, rate);
            break;
        case BCComponent::POLYMER:
            this->handlePolymerBC(bc, rate);
            break;
        case BCComponent::MICR:
            this->handleMicrBC(bc, rate);
            break;
        case BCComponent::OXYG:
            this->handleOxygBC(bc, rate);
            break;
        case BCComponent::UREA:
            this->handleUreaBC(bc, rate);
            break;
        case BCComponent::NONE:
            throw std::logic_error("you need to specify the component when RATE type is set in BC");
            break;
        }
        //TODO add support for enthalpy rate
        return {bc.bctype, rate};
    }

protected:
    //! \brief Boundary rate for a face, as stored by the producing run.
    //!
    //! \param[out] isMassRate Set when the returned vector holds component mass
    //!   rates, which can be imposed directly. Otherwise it holds reservoir
    //!   volumetric rates, which still have to be converted through a density.
    bool fluxBoundaryRate_(const unsigned int globalSpaceIdx,
                           const FaceDir::DirEnum dir,
                           RateVector& phaseVolRate,
                           bool& isMassRate) const
    {
        const auto& bcDirIndex = bcindex_(dir);
        const bool hasDeckBc = globalSpaceIdx < bcDirIndex.size() && bcDirIndex[globalSpaceIdx] != 0;
        if (hasDeckBc) {
            return false;
        }

        const auto* fluxData = this->fluxBoundaryData_();
        const auto* fluxStep = this->fluxBoundaryReportStep_();
        const auto fluxSlot = this->fluxBoundaryFaceSlot_(globalSpaceIdx, dir);
        if (!fluxData || !fluxStep || fluxSlot <= 0) {
            return false;
        }

        const auto mode = fluxData->header.mode;
        const auto fluxEnabled = (static_cast<int>(mode) & static_cast<int>(EclIO::FluxFile::Mode::Flux)) != 0;
        if (!fluxEnabled) {
            return false;
        }

        const auto faceIndex = static_cast<std::size_t>(fluxSlot - 1);
        const auto faceCount = static_cast<std::size_t>(fluxData->header.numBoundaryFaces);
        const auto phaseCount = static_cast<std::size_t>(fluxData->header.numPhases);
        const auto expectedSize = faceCount * phaseCount;
        if (phaseCount == 0 || faceIndex >= faceCount
            || fluxStep->massRates.size() != expectedSize)
        {
            return false;
        }

        phaseVolRate = 0.0;
        const auto base = faceIndex * phaseCount;
        std::size_t phaseSlot = 0;

        // The file stores component masses, which is the only thing a reduced
        // run can impose: for flow entering the sector the upstream cell lies
        // outside its grid, so it has no state with which to turn a phase
        // volumetric flux into components.
        isMassRate = true;
        const auto& storedRates = fluxStep->massRates;

        // The stored rates are totals across the face, positive into the
        // sector. The linearizer multiplies the boundary rate vector by the
        // face area, and opm-models orients boundary rates along the outer
        // normal, i.e. a flux into the domain is negative. Convert accordingly.
        const auto faceArea = this->fluxBoundaryFaceAreaAt_(globalSpaceIdx, dir);
        if (!(faceArea > 0.0)) {
            return false;
        }

        const auto areaScale = Scalar{-1} / faceArea;

        if (fluxData->header.hasPhase(EclIO::FluxFile::Phase::Oil)) {
            if (FluidSystem::phaseIsActive(oilPhaseIdx)) {
                phaseVolRate[FluidSystem::canonicalToActiveCompIdx(oilCompIdx)] =
                    storedRates[base + phaseSlot] * areaScale;
            }
            ++phaseSlot;
        }
        if (fluxData->header.hasPhase(EclIO::FluxFile::Phase::Water)) {
            if (FluidSystem::phaseIsActive(waterPhaseIdx)) {
                phaseVolRate[FluidSystem::canonicalToActiveCompIdx(waterCompIdx)] =
                    storedRates[base + phaseSlot] * areaScale;
            }
            ++phaseSlot;
        }
        if (fluxData->header.hasPhase(EclIO::FluxFile::Phase::Gas)) {
            if (FluidSystem::phaseIsActive(gasPhaseIdx)) {
                phaseVolRate[FluidSystem::canonicalToActiveCompIdx(gasCompIdx)] =
                    storedRates[base + phaseSlot] * areaScale;
            }
        }


        return true;
    }

    int fluxBoundaryFaceSlot_(const unsigned int globalSpaceIdx, const FaceDir::DirEnum dir) const
    {
        const auto& dirData = fluxBoundaryFaceIndex_(dir);
        if (globalSpaceIdx >= dirData.size()) {
            return 0;
        }

        return dirData[globalSpaceIdx];
    }

    //! \brief Sector-boundary faces reaching this cell through a non-neighbour
    //!   connection.
    //!
    //! \details Empty for all but a handful of cells, and empty everywhere when
    //!   the sector boundary follows the grid axes. What crosses one of these
    //!   has to enter as a source, because the reduced run has no such
    //!   connection of its own to carry it.
    const std::vector<int>& fluxNncFacesAt_(const unsigned int globalSpaceIdx) const
    {
        static const std::vector<int> none{};

        if (globalSpaceIdx >= this->fluxNncFaceIndex_.size()) {
            return none;
        }

        return this->fluxNncFaceIndex_[globalSpaceIdx];
    }

    const FluxBoundary::Face* fluxBoundaryFace_(const unsigned int globalSpaceIdx,
                                                const FaceDir::DirEnum dir) const
    {
        if (!this->fluxBoundary_) {
            return nullptr;
        }

        return this->fluxBoundary_->faceFromSlot(this->fluxBoundaryFaceSlot_(globalSpaceIdx, dir));
    }

    //! \brief The sector-boundary face a cell's \p bfIdx'th boundary face is,
    //!        or nullptr if it is not one.
    //!
    //! \details The discretisation numbers a cell's boundary faces by position
    //!   and knows nothing of their direction, while everything read out of the
    //!   .FLUX file is keyed by direction.
    const FluxBoundary::Face* fluxBoundaryFaceAtOrdinal_(const unsigned int globalSpaceIdx,
                                                         const unsigned int bfIdx) const
    {
        const auto dir = this->fluxBoundaryDirAt_(globalSpaceIdx, bfIdx);

        return (dir == FaceDir::Unknown)
            ? nullptr
            : this->fluxBoundaryFace_(globalSpaceIdx, dir);
    }

    const EclIO::FluxFile::Data* fluxBoundaryData_() const
    {
        if (!this->fluxBoundary_) {
            return nullptr;
        }

        return &this->fluxBoundary_->data();
    }

    //! \brief The boundary record that applies to the current time step.
    //!
    //! \details Resolved once per time step by refreshFluxBoundaryRecord_().
    //!          This is read per boundary face per Newton iteration during
    //!          assembly, so it must not search the record list.
    const EclIO::FluxFile::ReportStep* fluxBoundaryReportStep_() const
    {
        return this->fluxBoundaryActiveRecord_;
    }

    //! \brief Point the cached boundary record at the record covering the step
    //!        (\p time, \p stepEnd].
    //!
    //! \details A record holds the parent's state at the END of the interval
    //!   it covers. Held over the whole interval, that state reaches the
    //!   sector early: a front that crosses the boundary late in a long
    //!   interval is imposed from its start, and when the parent wrote records
    //!   only every year or two, as a restart-based file often has, that is a
    //!   year or two early. Interpolating between the previous record's state
    //!   and this one's, at the end of the step as an implicit step wants it,
    //!   keeps the boundary on the parent's trajectory, and gives back the
    //!   recorded state exactly whenever a step ends on a record's time.
    //!
    //!   Only the state is interpolated. Flux-mode mass rates are averages
    //!   over the record's interval and stay as recorded.
    void refreshFluxBoundaryRecord_(const double time, const double stepEnd)
    {
        const auto* data = this->fluxBoundaryData_();

        this->fluxBoundaryActiveRecord_ = (data == nullptr)
            ? nullptr
            : FluxBoundary::selectRecordAt(*data, time);

        if (!this->fluxBoundaryInterpolate_
            || (this->fluxBoundaryActiveRecord_ == nullptr)
            || (this->fluxBoundaryActiveRecord_ == &data->reportSteps.front())
            || !(this->fluxBoundaryActiveRecord_->stepLength > 0.0))
        {
            return;
        }

        const auto& record = *this->fluxBoundaryActiveRecord_;
        const auto& previous = *(this->fluxBoundaryActiveRecord_ - 1);
        const auto weight = std::clamp((stepEnd - record.startTime) / record.stepLength, 0.0, 1.0);
        if (!(weight < 1.0)) {
            return;
        }

        auto& blended = this->fluxBoundaryBlendedRecord_;
        blended = record;

        const auto blend = [weight](std::vector<double>& values, const std::vector<double>& before)
        {
            if (before.size() != values.size()) {
                return;
            }
            for (std::size_t i = 0; i < values.size(); ++i) {
                if (std::isfinite(values[i]) && std::isfinite(before[i])) {
                    values[i] = before[i] + weight * (values[i] - before[i]);
                }
            }
        };

        blend(blended.pressures, previous.pressures);
        blend(blended.swat, previous.swat);
        blend(blended.sgas, previous.sgas);
        blend(blended.rs, previous.rs);
        blend(blended.rv, previous.rv);
        blend(blended.temperature, previous.temperature);
        blend(blended.relPerm, previous.relPerm);
        blend(blended.capPressure, previous.capPressure);
        blend(blended.externalRegionSums, previous.externalRegionSums);

        this->fluxBoundaryActiveRecord_ = &blended;
    }

    //! \brief The parent run's summary vectors, or nullptr when unavailable.
    const ParentSummary* fluxParentSummary_() const
    {
        return this->fluxParentSummaryData_.get();
    }

    //! \brief Simulation time, in seconds, at the end of the step that is
    //!        about to be taken, where the boundary state and the parent's
    //!        summary vectors are sampled for it.
    double fluxStepEndTime_() const
    {
        return static_cast<double>(this->simulator().time())
            + static_cast<double>(this->simulator().timeStepSize());
    }

public:

    template<class Serializer>
    void serializeOp(Serializer& serializer)
    {
        serializer(static_cast<BaseType&>(*this));
        serializer(first_step_);
        serializer(drift_);
        serializer(wellModel_);
        serializer(aquiferModel_);
        serializer(tracerModel_);
        serializer(temperatureModel_);
        serializer(*materialLawManager_);
    }

    const GlobalEqVector& drift() const
    {
        return drift_;
    }

    //! For reporting, which needs to see what lives outside the grid.
    const std::vector<std::unique_ptr<FlowAuxCellModule<TypeTag>>>& auxCellModules() const
    { return auxCellModules_; }

    /*!
     * \brief Give the auxiliary cell modules their turn at the linear system.
     *
     * Their degrees of freedom are assembled by the model's own local residual, like grid
     * cells, so most of them have nothing to add here.  What they do own is the rows of
     * the cells they have preallocated but not yet put to use: those have no volume and
     * no connections, so nothing else writes to them, and an all-zero row is a singular
     * matrix.  Conditioning them is the module's business, since only it knows which they
     * are.
     *
     * Deliberately not the linearizer's linearizeAuxiliaryEquations(), which drives
     * *every* auxiliary module -- the well model among them, whose contribution is
     * applied separately and would be counted twice.
     */
    void linearizeAuxCellModules(GetPropType<TypeTag, Properties::SparseMatrixAdapter>& matrix,
                                 GlobalEqVector& residual)
    {
        for (const auto& module : auxCellModules_) {
            module->linearize(matrix, residual);
        }
    }

    //! \brief The parent run's summary vectors, or nullptr when unavailable.
    //!
    //! \details Only populated in a USEFLUX run for which a parent summary
    //!          source could be located.
    const ParentSummary* fluxParentSummary() const
    {
        return this->fluxParentSummaryData_.get();
    }

private:
    Implementation& asImp_()
    { return *static_cast<Implementation *>(this); }

    const Implementation& asImp_() const
    { return *static_cast<const Implementation *>(this); }

protected:
    //! \brief Snapshot the explicit quantities before the timestep runs.
    //!
    //! updateExplicitQuantities_() advances all of these from the current
    //! solution at the start of every timestep.  They are running extrema, so
    //! a step that later fails would otherwise leave the maxima raised (and the
    //! minima lowered) by a solution that is being thrown away.
    virtual void captureBeginTimeStepState_()
    {
        prev_timestep_state_.first_step = first_step_;
        prev_timestep_state_.max_polymer_adsorption = this->polymer_.maxAdsorption;
        prev_timestep_state_.max_oil_saturation = this->maxOilSaturation_;
        prev_timestep_state_.max_water_saturation = this->maxWaterSaturation_;
        prev_timestep_state_.min_ref_pressure = this->minRefPressure_;
        prev_timestep_state_.rock_comp_trans_mult_val = this->rockCompTransMultVal_;
        if (materialLawManager_ && materialLawManager_->enableHysteresis()) {
            materialLawManager_->captureBeginTimeStepState();
        }
    }

    //! \brief Put the explicit quantities back as they were before the step.
    virtual void restoreBeginTimeStepState_()
    {
        first_step_ = prev_timestep_state_.first_step;
        this->polymer_.maxAdsorption = prev_timestep_state_.max_polymer_adsorption;
        this->maxOilSaturation_ = prev_timestep_state_.max_oil_saturation;
        this->maxWaterSaturation_ = prev_timestep_state_.max_water_saturation;
        this->minRefPressure_ = prev_timestep_state_.min_ref_pressure;
        this->rockCompTransMultVal_ = prev_timestep_state_.rock_comp_trans_mult_val;
        if (materialLawManager_ && materialLawManager_->enableHysteresis()) {
            materialLawManager_->restoreBeginTimeStepState();
        }
    }

    void finishTransmissibilities_()
    {
        this->transmissibilities_.finishInit(
            [&vg = this->simulator().vanguard()](const unsigned int index) {
                return vg.gridIdxToEquilGridIdx(index);
            });
    }

    template<class EclWriterType>
    bool prepareTransmissibilityOutput_(EclWriterType& eclWriter,
                                        const bool enableEclOutput)
    {
        auto& simulator = this->simulator();
        bool localTransmissibilitiesFinished = false;

        if (enableEclOutput) {
            // Parallel TRANX, TRANY, TRANZ and NNC output requires the global
            // grid's transmissibilities on the I/O rank -- with gathered LGR
            // records, only for the connections that no rank records, below.
            if (simulator.vanguard().grid().comm().size() > 1) {
                bool wholeGridTransNeeded = simulator.vanguard().grid().comm().rank() == 0;
                // A parallel Conforming LGR run writes from the I/O rank's refined
                // output grid; the transmissibility must live there too.  The
                // backend decides, not that grid: it exists on the I/O rank only,
                // and the gathered path below is collective.
                bool refinedOutputGrid = false;
                if constexpr (requires { simulator.vanguard().eclOutputTransmissibility(); }) {
                    refinedOutputGrid = simulator.vanguard().lgrBackend() == "conforming";
                    if (refinedOutputGrid) {
                        if (wholeGridTransNeeded) {
                            eclWriter.setTransmissibilities(&simulator.vanguard().eclOutputTransmissibility());
                        }
                        wholeGridTransNeeded = false;
                    }
                }
                // Parallel LGR: reuse the simulator's own (distributed) transmissibilities for the
                // INIT output -- each rank contributes its interior connections, gathered on the
                // I/O rank and keyed by level-Cartesian indices so the output walk over the global
                // (equil) grid can look them up directly. This reuses the values already computed
                // in parallel for the simulation itself instead of recomputing a whole-grid
                // transmissibility. The local transmissibilities are built here, before the INIT
                // write, and reported as finished so that they are not built again.
                if constexpr (std::is_same_v<GetPropType<TypeTag, Properties::Grid>, Dune::CpGrid>) {
                    // Gate on the deck's LGRs -- the same condition the writer keys its LGR
                    // output walk on (writeInit / extractOutputTransAndNNC). Grid refinement
                    // alone (grid().maxLevel() > 0, e.g. after adapt()) must not enable the
                    // gathered mode: the writer would look up different keys than the ones
                    // the leaf-grid walk records.
                    if (!refinedOutputGrid && simulator.vanguard().eclState().getLgrs().size() > 0) {
                        this->finishTransmissibilities_();
                        localTransmissibilitiesFinished = true;
                        const auto& localTrans = simulator.problem().eclTransmissibilities();
                        eclWriter.setGatheredLgrTrans(
                            gatherLgrOutputTrans(simulator.vanguard().grid(),
                                                 simulator.vanguard().gridView(),
                                                 [&localTrans](unsigned c1, unsigned c2)
                                                 { return static_cast<double>(localTrans.transmissibility(c1, c2)); }));
                        // The output values (TRANX/Y/Z and NNC) come from the gathered records.
                        // The whole-grid transmissibilities set below tell, for a connection
                        // that no rank records, whether its face is sealed.
                    }
                }
                if (wholeGridTransNeeded) {
                    eclWriter.setTransmissibilities(&simulator.vanguard().globalTransmissibility());
                }
            }
            else {
                this->finishTransmissibilities_();
                localTransmissibilitiesFinished = true;
                eclWriter.setTransmissibilities(&simulator.problem().eclTransmissibilities());
            }

            std::function<unsigned int(unsigned int)> equilGridToGrid =
                [&simulator](const unsigned int index) {
                    return simulator.vanguard().gridEquilIdxToGridIdx(index);
                };
            eclWriter.extractOutputTransAndNNC(equilGridToGrid);
        }

        simulator.vanguard().releaseGlobalTransmissibilities();
        return localTransmissibilitiesFinished;
    }

    void initializeSimulatorTime_()
    {
        auto& simulator = this->simulator();
        const auto& eclState = simulator.vanguard().eclState();
        const auto& schedule = simulator.vanguard().schedule();

        simulator.setStartTime(TimeService::to_time_t(schedule.getStartTime()));
        simulator.setEndTime(TimeService::to_time_t(schedule.simTime(schedule.size() - 1)));

        // Represent initialization as a zero-length episode before report step zero.
        simulator.setEpisodeIndex(-1);
        simulator.setEpisodeLength(0.0);

        this->initGravity_(eclState);

        if (this->enableTuning_) {
            // If support for the TUNING keyword is enabled, get the initial time
            // stepping parameters from it instead of from command line parameters.
            const auto& tuning = schedule[0].tuning();
            this->initialTimeStepSize_ = tuning.TSINIT.value_or(-1.0);
            this->maxTimeStepAfterWellEvent_ = tuning.TMAXWC;
        }
    }

    void initializeModelProperties_()
    {
        auto& simulator = this->simulator();

        if (FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx) &&
            FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)) {
            // Indexed by DOF, auxiliary cells included.
            this->maxOilSaturation_.resize(this->model().numTotalDof(), 0.0);
        }

        this->readRockParameters_(
            simulator.vanguard().cellCenterDepths(),
            [&simulator](const unsigned int index) {
                std::array<int, dim> coords;
                simulator.vanguard().cartesianCoordinate(index, coords);
                std::ranges::transform(coords, coords.begin(),
                                       [](const auto coordinate) { return coordinate + 1; });
                return coords;
            });

        this->readMaterialParameters_();
        this->readThermalParameters_();
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
                const auto& iq = model.intensiveQuantities(dofIdx, /*timeIdx=*/ 0);
                func(dofIdx, iq);
        }
        OPM_END_PARALLEL_TRY_CATCH(failureMsg, vanguard.grid().comm());
    }

    bool updateMaxOilSaturation_()
    {
        OPM_TIMEBLOCK(updateMaxOilSaturation);
        int episodeIdx = this->episodeIndex();

        // we use VAPPARS
        if (this->vapparsActive(episodeIdx)) {
            this->updateProperty_("FlowProblem::updateMaxOilSaturation_() failed:",
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
        OPM_TIMEBLOCK_LOCAL(updateMaxOilSaturation, Subsystem::SatProps);
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

        this->updateProperty_("FlowProblem::updateMaxWaterSaturation_() failed:",
                              [this](unsigned compressedDofIdx, const IntensiveQuantities& iq)
                              {
                                  this->updateMaxWaterSaturation_(compressedDofIdx,iq);
                               });
        return true;
    }


    bool updateMaxWaterSaturation_(unsigned compressedDofIdx, const IntensiveQuantities& iq)
    {
        OPM_TIMEBLOCK_LOCAL(updateMaxWaterSaturation, Subsystem::SatProps);
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

        this->updateProperty_("FlowProblem::updateMinPressure_() failed:",
                              [this](unsigned compressedDofIdx, const IntensiveQuantities& iq)
                              {
                                  this->updateMinPressure_(compressedDofIdx,iq);
                              });
        return true;
    }

    bool updateMinPressure_(unsigned compressedDofIdx, const IntensiveQuantities& iq){
        OPM_TIMEBLOCK_LOCAL(updateMinPressure, Subsystem::PvtProps);
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
        this->authorAuxCellRegions_();
        ////////////////////////////////
        // porosity
        updateReferencePorosity_();
        this->referencePorosity_[1] = this->referencePorosity_[0];
        ////////////////////////////////

        ////////////////////////////////
        // rock fraction
        updateRockFraction_();
        this->rockFraction_[1] = this->rockFraction_[0];
        ////////////////////////////////

        ////////////////////////////////
        // fluid-matrix interactions (saturation functions; relperm/capillary pressure)
        materialLawManager_ = std::make_shared<EclMaterialLawManager>();
        materialLawManager_->initFromState(eclState);
        // Grid DOFs only: endpoint scaling/hysteresis are keyed on grid entities.
        // Auxiliary cells go through auxCellSaturationProxy_() instead.
        materialLawManager_->initParamsForElements(eclState, this->model().numGridDof(),
                                                   this-> template fieldPropIntTypeOnLeafAssigner_<int>(),
                                                   this-> lookupIdxOnLevelZeroAssigner_());
        ////////////////////////////////
    }

    void readThermalParameters_()
    {
        if constexpr (energyModuleType == EnergyModules::FullyImplicitThermal ||
                      energyModuleType == EnergyModules::SequentialImplicitThermal )
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

        this->referencePorosity_[/*timeIdx=*/0].resize(this->model().numTotalDof());
        this->authorAuxCellPorosity_();

        const auto& fp = eclState.fieldProps();
        const std::vector<double> porvData = this -> fieldPropDoubleOnLeafAssigner_()(fp, "PORV");
        for (std::size_t dofIdx = 0; dofIdx < numDof; ++dofIdx) {
            int sfcdofIdx = simulator.vanguard().gridEquilIdxToGridIdx(dofIdx);
            Scalar poreVolume = porvData[dofIdx];

            // we define the porosity as the accumulated pore volume divided by the
            // geometric volume of the element. Note that -- in pathetic cases -- it can
            // be larger than 1.0!
            Scalar dofVolume = simulator.model().dofTotalVolume(sfcdofIdx);
            assert(dofVolume > 0.0);
            this->referencePorosity_[/*timeIdx=*/0][sfcdofIdx] = poreVolume/dofVolume;
        }
    }

    void updateRockFraction_()    {

        const bool solveEnergyEquation = (energyModuleType == EnergyModules::FullyImplicitThermal ||
        energyModuleType == EnergyModules::SequentialImplicitThermal);
        if (!solveEnergyEquation)
            return;

        const auto& simulator = this->simulator();
        const auto& vanguard = simulator.vanguard();
        const auto& eclState = vanguard.eclState();

        std::size_t numDof = this->model().numGridDof();
        this->rockFraction_[/*timeIdx=*/0].resize(this->model().numTotalDof(), 0.0);

        // Zero for auxiliary cells: they store no rock heat.
        // For the energy equation, we need the volume of the rock.
        // The volume of the rock is computed by rockFraction * geometric volume of the element.
        // The reference porosity is defined as porosity * ntg * pore-volume-multiplier.
        // A common practice in reservoir simulation is to use large pore-volume-multipliers in boundary cells
        // to model boundary conditions other than no-flow. This may result in reference porosities that are larger than 1.
        // A simple (1-reference porosity) * geometric volume of the element may give unphysical results.
        // We therefore instead consider the pore-volume-multiplier as a volume multiplier. The rock fraction is thus given by
        // (1 - porosity * ntg) * pore-volume-multiplier = (1 - porosity * ntg) * reference porosity / (porosity * ntg)
        const auto& fp = eclState.fieldProps();
        const std::vector<double> poroData = this->fieldPropDoubleOnLeafAssigner_()(fp, "PORO");
        const std::vector<double> ntgData = this->fieldPropDoubleOnLeafAssigner_()(fp, "NTG");

        for (std::size_t dofIdx = 0; dofIdx < numDof; ++dofIdx) {
            const auto ntg = ntgData[dofIdx];
            const auto poro_eff = ntg * poroData[dofIdx];
            const int sfcdofIdx = simulator.vanguard().gridEquilIdxToGridIdx(dofIdx);
            const auto rock_fraction = (1 - poro_eff) * this->referencePorosity_[/*timeIdx=*/0][sfcdofIdx] / poro_eff;
            this->rockFraction_[/*timeIdx=*/0][sfcdofIdx] = rock_fraction;
        }
    }

    virtual void readInitialCondition_()
    {
        // TODO: whether we should move this to FlowProblemBlackoil
        const auto& simulator = this->simulator();
        const auto& vanguard = simulator.vanguard();
        const auto& eclState = vanguard.eclState();

        if (eclState.getInitConfig().hasEquil())
            readEquilInitialCondition_();
        else
            readExplicitInitialCondition_();

        //initialize min/max values
        std::size_t numElems = this->model().numGridDof();
        for (std::size_t elemIdx = 0; elemIdx < numElems; ++elemIdx) {
            const auto& fs = asImp_().initialFluidStates()[elemIdx];
            if (!this->maxWaterSaturation_.empty() && waterPhaseIdx > -1)
                this->maxWaterSaturation_[elemIdx] = std::max(this->maxWaterSaturation_[elemIdx], fs.saturation(waterPhaseIdx));
            if (!this->maxOilSaturation_.empty() && oilPhaseIdx > -1)
                this->maxOilSaturation_[elemIdx] = std::max(this->maxOilSaturation_[elemIdx], fs.saturation(oilPhaseIdx));
            if (!this->minRefPressure_.empty() && refPressurePhaseIdx_() > -1)
                this->minRefPressure_[elemIdx] = std::min(this->minRefPressure_[elemIdx], fs.pressure(refPressurePhaseIdx_()));
        }
    }

    virtual void readEquilInitialCondition_() = 0;
    virtual void readExplicitInitialCondition_() = 0;

    // update the hysteresis parameters of the material laws for the whole grid
    bool updateHysteresis_()
    {
        if (!materialLawManager_->enableHysteresis())
            return false;

        // we need to update the hysteresis data for _all_ elements (i.e., not just the
        // interior ones) to avoid desynchronization of the processes in the parallel case!
        this->updateProperty_("FlowProblem::updateHysteresis_() failed:",
                              [this](unsigned compressedDofIdx, const IntensiveQuantities& iq)
                              {
                                  materialLawManager_->updateHysteresis(iq.fluidState(), compressedDofIdx);
                              });
        return true;
    }


    bool updateHysteresis_(unsigned compressedDofIdx, const IntensiveQuantities& iq)
    {
        OPM_TIMEBLOCK_LOCAL(updateHysteresis_, Subsystem::SatProps);
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

protected:
    /*!
     * \brief Set the gravity vector from the run's configuration.
     *
     * Gravity is switched off by the NOGRAV keyword from Frontsim or by
     * setting EnableGravity to false, else the standard value of the gravity
     * constant at sea level on earth is used.  The vertical component is
     * positive, which is the magnitude the equilibration expects.
     */
    void initGravity_(const EclipseState& eclState)
    {
        this->gravity_ = 0.0;

        if (Parameters::Get<Parameters::EnableGravity>() &&
            eclState.getInitConfig().hasGravity())
        {
            // unit::gravity is 9.80665 m^2/s--i.e., standard measure at Tellus equator.
            this->gravity_[dim - 1] = unit::gravity;
        }
    }

    //! Must run before the model sizes its per-DOF containers; connections are built
    //! after, since the global DOF offset is assigned on registration.
    template <class Module>
    Module& registerAuxCellModule_(std::unique_ptr<Module> module)
    {
        // Aux DOFs are not in the communication index set yet, so each rank would
        // solve its own copy of them.
        if (this->simulator().gridView().comm().size() > 1) {
            OPM_THROW(std::runtime_error,
                      "Degrees of freedom outside the grid (auxiliary cells: numerical "
                      "aquifers or embedded fractures) are not supported in parallel "
                      "yet. Run on one process; for numerical aquifers "
                      "--numerical-aquifer-mode=grid also avoids them.");
        }

        auto& ref = *module;
        this->simulator().model().addAuxiliaryModule(&ref);
        this->auxCellModules_.push_back(std::move(module));
        ref.buildConnections();
        return ref;
    }

    //! Auxiliary cells borrow their initialisation partner's per-element tables; exact
    //! only if both share a saturation region.
    unsigned auxCellSaturationProxy_(unsigned globalDofIdx) const
    {
        if (globalDofIdx < this->model().numGridDof()) {
            return globalDofIdx;
        }

        for (const auto& module : this->auxCellModules_) {
            const auto begin = static_cast<unsigned>(module->dofOffset());
            if ((globalDofIdx >= begin) && (globalDofIdx < begin + module->numDofs())) {
                return module->initialisationPartner(globalDofIdx - begin);
            }
        }

        return 0;
    }

public:
    //! Newton temperature step limit of an auxiliary cell, or <= 0 for the global one.
    Scalar auxMaxTemperatureChange(unsigned globalDofIdx) const
    {
        if (globalDofIdx < this->model().numGridDof()) {
            return 0.0;
        }
        for (const auto& module : this->auxCellModules_) {
            const auto begin = static_cast<unsigned>(module->dofOffset());
            if ((globalDofIdx >= begin) && (globalDofIdx < begin + module->numDofs())) {
                return static_cast<Scalar>(module->maxTemperatureChange());
            }
        }
        return 0.0;
    }

protected:
    Scalar auxCellDepth_(unsigned globalSpaceIdx) const
    {
        for (const auto& module : this->auxCellModules_) {
            const auto begin = static_cast<unsigned>(module->dofOffset());
            if ((globalSpaceIdx >= begin) && (globalSpaceIdx < begin + module->numDofs())) {
                return module->depth(globalSpaceIdx - begin);
            }
        }

        return 0.0;
    }

    //! Region arrays are read per DOF, so they must cover aux DOFs too; empty means one region.
    void authorAuxCellRegions_()
    {
        const auto numTotalDof = this->model().numTotalDof();
        if (numTotalDof == this->model().numGridDof()) {
            return;
        }

        const auto extend = [numTotalDof](auto& numbers, auto&& regionOf) {
            if (numbers.empty()) {
                return;
            }

            numbers.resize(numTotalDof, 0);
            regionOf(numbers);
        };

        extend(this->pvtnum_, [this](auto& numbers) {
            for (const auto& module : this->auxCellModules_) {
                for (unsigned localIdx = 0; localIdx < module->numDofs(); ++localIdx) {
                    numbers[module->localToGlobalDof(localIdx)] = module->pvtRegionIndex(localIdx);
                }
            }
        });

        extend(this->satnum_, [this](auto& numbers) {
            for (const auto& module : this->auxCellModules_) {
                for (unsigned localIdx = 0; localIdx < module->numDofs(); ++localIdx) {
                    numbers[module->localToGlobalDof(localIdx)] = module->satRegionIndex(localIdx);
                }
            }
        });

        // No solvent/polymer support for aux cells yet; region zero as placeholder.
        extend(this->miscnum_, [](auto&) {});
        extend(this->plmixnum_, [](auto&) {});
    }

    //! Re-read what the auxiliary cell modules author; a new connection
    //! (topologyChanged) also needs the sparsity pattern rebuilt.
    void refreshAuxCellModules_(const bool topologyChanged)
    {
        if (this->auxCellModules_.empty()) {
            return;
        }

        this->model().updateAuxiliaryDofVolumes();
        this->authorAuxCellPorosity_();
        this->referencePorosity_[1] = this->referencePorosity_[0];
        this->applyAuxCellTransmissibilities_();

        if (topologyChanged) {
            this->model().linearizer().eraseMatrix();

            // A cell that has only just appeared has no history: it did not exist at the
            // start of the time step, so there is no earlier state to have changed from.
            // Give it the state it starts with, which is what says its mass has not moved
            // by coming into existence.  The grid cells keep theirs -- that is the state
            // the step is measuring the change from.
            auto& previous = this->model().solution(/*timeIdx=*/1);
            const auto& current = this->model().solution(/*timeIdx=*/0);
            for (unsigned dofIdx = this->model().numGridDof();
                 dofIdx < this->model().numTotalDof(); ++dofIdx)
            {
                previous[dofIdx] = current[dofIdx];
            }

        }

        this->model().invalidateAndUpdateIntensiveQuantities(/*timeIdx=*/0);

        if (topologyChanged) {
            this->seedAuxCellStorageCache_();
        }
    }

    /*!
     * \brief Give the auxiliary cells a start-of-step storage term.
     *
     * The linearizer keeps the accumulation of the previous time level in a cache, and
     * fills it on the first Newton iteration of a step -- which a cell that appears part
     * way through one has already missed.  Reaching for it then throws rather than
     * quietly using a stale number, which is the right behaviour and exactly why it has
     * to be filled here.
     *
     * The value is the accumulation the cell's own starting state implies, so that its
     * mass is unchanged by the cell having come into existence.
     */
    void seedAuxCellStorageCache_()
    {
        using LocalResidual = GetPropType<TypeTag, Properties::LocalResidual>;
        using EqVector = GetPropType<TypeTag, Properties::EqVector>;

        using Linearizer = GetPropType<TypeTag, Properties::Linearizer>;

        // Auxiliary cells carrying the model's equations exist only under a linearizer
        // that assembles per degree of freedom, and only that one pairs with the local
        // residual whose accumulation can be formed from a degree-of-freedom index alone.
        if constexpr (!Linearizer::assemblesAuxiliaryDofEquations) {
            return;
        }
        else {

        auto& model = this->model();
        if (!model.enableStorageCache()) {
            return;
        }

        // Formed from the current state rather than from the previous time level's
        // intensive quantities, which need not be kept at all when the first iteration's
        // storage is recycled.  For these cells the two are the same: a cell that has
        // just appeared starts the step where it is, and one that has not been handed out
        // has not moved.
        for (unsigned dofIdx = model.numGridDof(); dofIdx < model.numTotalDof(); ++dofIdx) {
            if (model.storageCacheIsUpToDate(dofIdx, /*timeIdx=*/1)) {
                continue;
            }

            EqVector storage;
            LocalResidual::template computeStorage<Scalar>
                (storage, model.intensiveQuantities(dofIdx, /*timeIdx=*/0));

            model.updateCachedStorage(dofIdx, /*timeIdx=*/1, storage);
        }

        } // if constexpr
    }

    //! Porosity times the DOF volume must give back the module's pore volume.
    void authorAuxCellPorosity_()
    {
        for (const auto& module : this->auxCellModules_) {
            for (unsigned localIdx = 0; localIdx < module->numDofs(); ++localIdx) {
                const auto globalIdx = static_cast<unsigned>(module->localToGlobalDof(localIdx));
                const auto bulkVolume = module->bulkVolume(localIdx);

                this->referencePorosity_[/*timeIdx=*/0][globalIdx] = (bulkVolume > 0.0)
                    ? module->poreVolume(localIdx) / bulkVolume
                    : 0.0;
            }
        }
    }

    //! Call once the grid's transmissibilities are final.
    void applyAuxCellTransmissibilities_()
    {
        using ConnectionVector = std::vector<typename FlowAuxCellModule<TypeTag>::Connection>;

        for (const auto& module : this->auxCellModules_) {
            ConnectionVector conns;
            module->connections(conns);

            for (const auto& conn : conns) {
                this->transmissibilities_.setTransmissibility(conn.dof1, conn.dof2, conn.trans);

                if constexpr (enableFullyImplicitThermal) {
                    this->transmissibilities_
                        .setThermalHalfTrans(conn.dof1, conn.dof2, conn.thermalHalfTrans12);
                    this->transmissibilities_
                        .setThermalHalfTrans(conn.dof2, conn.dof1, conn.thermalHalfTrans21);
                }
            }
        }
    }

    struct PffDofData_
    {
        ConditionalStorage<enableFullyImplicitThermal, Scalar> thermalHalfTransIn;
        ConditionalStorage<enableFullyImplicitThermal, Scalar> thermalHalfTransOut;
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

                if constexpr (enableFullyImplicitThermal) {
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

    virtual void updateExplicitQuantities_(int episodeIdx, int timeStepSize, bool first_step_after_restart) = 0;

    void readBoundaryConditions_()
    {
        const auto& simulator = this->simulator();
        const auto& vanguard = simulator.vanguard();
        this->readFluxBoundary_();

        // Sized even without BCCON: a USEFLUX boundary also makes the boundary
        // conditions non-trivial, and boundaryConditionIndex() is then read for
        // every boundary face, e.g. by the tracer model.
        const unsigned numElems = vanguard.gridView().size(/*codim=*/0);
        bcindex_.resize(numElems, 0);

        const auto& bcconfig = vanguard.eclState().getSimulationConfig().bcconfig();
        if (bcconfig.size() > 0) {
            nonTrivialBoundaryConditions_ = true;

            std::size_t numCartDof = vanguard.cartesianSize();
            std::vector<int> cartesianToCompressedElemIdx(numCartDof, -1);

            for (unsigned elemIdx = 0; elemIdx < numElems; ++elemIdx)
                cartesianToCompressedElemIdx[vanguard.cartesianIndex(elemIdx)] = elemIdx;

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

    void readFluxBoundary_()
    {
        const auto& ioConfig = this->simulator().vanguard().eclState().getIOConfig();
        if (!ioConfig.getUseFlux()) {
            this->fluxBoundary_.reset();
            this->fluxBoundaryFaceIndex_.data = {};
            this->fluxBoundaryFaceArea_.data = {};
            this->fluxBoundaryFaceOrdinal_.data = {};
            this->fluxBoundaryFaceDepth_.data = {};
            return;
        }

        const auto inputDir = std::filesystem::path{ioConfig.getInputDir()};
        const auto baseName = ioConfig.getUseFluxInputBaseName();
        const auto selectedFluxPath = FluxBoundary::selectInputPath(inputDir, baseName);
        if (!std::filesystem::exists(selectedFluxPath)) {
            std::ostringstream msg;
            msg << "USEFLUX requested but no FLUX file found: expected '" << baseName
                << ".FLUX' or '" << baseName << ".FLUXdddd' in '" << inputDir.string() << "'";
            throw std::runtime_error(msg.str());
        }

        const auto fluxData = EclIO::FluxFile::read(selectedFluxPath.string());
        const auto fluxModeEnabled = (static_cast<int>(fluxData.header.mode)
                                      & static_cast<int>(EclIO::FluxFile::Mode::Flux)) != 0;
        if (fluxModeEnabled) {
            for (std::size_t stepIdx = 0; stepIdx < fluxData.reportSteps.size(); ++stepIdx) {
                const auto& masses = fluxData.reportSteps[stepIdx].massRates;
                for (std::size_t rateIdx = 0; rateIdx < masses.size(); ++rateIdx) {
                    const auto value = masses[rateIdx];
                    if (!std::isfinite(value)) {
                        std::ostringstream msg;
                        msg << "Invalid USEFLUX FLUX-mode component mass rate in '"
                            << selectedFluxPath.string() << "' at report-step index "
                            << stepIdx << ", rate index " << rateIdx
                            << ": " << value;
                        throw std::runtime_error(msg.str());
                    }
                }
            }
        }

        const auto& vanguard = this->simulator().vanguard();

        // The file names the region's cells by their Cartesian index in the
        // PARENT grid. A reduced run on the parent deck shares that grid, but
        // one on a model cut out to the region's box has a grid of its own,
        // numbered from the box corner, and a parent index read as one of its
        // own lands on some unrelated cell or none at all -- which on one field
        // left 553 of 583 boundary faces without a grid face and the sector
        // all but sealed.
        const auto& header = fluxData.header;
        const auto& dims = vanguard.cartesianDimensions();
        const bool sameGrid = (dims[0] == header.parentNx)
            && (dims[1] == header.parentNy)
            && (dims[2] == header.parentNz);
        const bool boxGrid = (dims[0] == header.boxNx)
            && (dims[1] == header.boxNy)
            && (dims[2] == header.boxNz);

        if (!sameGrid && !boxGrid) {
            throw std::runtime_error(fmt::format(
                "USEFLUX: '{}' describes a region of a {}x{}x{} grid, boxed in {}x{}x{} "
                "cells from ({},{},{}), and this run's grid is {}x{}x{}, which is "
                "neither",
                selectedFluxPath.string(),
                header.parentNx, header.parentNy, header.parentNz,
                header.boxNx, header.boxNy, header.boxNz,
                header.boxI1, header.boxJ1, header.boxK1,
                dims[0], dims[1], dims[2]));
        }

        const auto toThisGrid = [&header, sameGrid](const int parentCell)
        {
            if (sameGrid) {
                return parentCell;
            }

            const int i = parentCell % header.parentNx - (header.boxI1 - 1);
            const int j = (parentCell / header.parentNx) % header.parentNy - (header.boxJ1 - 1);
            const int k = parentCell / (header.parentNx * header.parentNy) - (header.boxK1 - 1);

            if ((i < 0) || (i >= header.boxNx)
                || (j < 0) || (j >= header.boxNy)
                || (k < 0) || (k >= header.boxNz))
            {
                return -1;
            }

            return i + header.boxNx * (j + header.boxNy * k);
        };

        // Where each of the region's cells sits in THIS rank's grid. The file
        // numbers them in region order, which a serial run happens to activate
        // them in, but a parallel rank holds an arbitrary subset under its own
        // numbering and has to look each one up. Cells that belong to another
        // rank come back as -1 and their faces are left to whoever holds them.
        std::vector<int> localToActive(fluxData.localToGlobal.size(), -1);
        for (std::size_t local = 0; local < fluxData.localToGlobal.size(); ++local) {
            const auto globalCell = fluxData.localToGlobal[local];
            if (globalCell >= 0) {
                const auto cell = toThisGrid(globalCell);
                if (cell >= 0) {
                    localToActive[local] = vanguard.compressedIndex(cell);
                }
            }
        }

        auto fluxBoundary = FluxBoundary::fromData(fluxData, localToActive);

        const auto numElems = this->simulator().vanguard().gridView().size(/*codim=*/0);
        this->fluxBoundaryFaceIndex_.resize(numElems, 0);
        this->fluxBoundaryFaceIndex_.data = fluxBoundary.buildDirectionalFaceIndices(numElems);
        this->fluxNncFaceIndex_ = fluxBoundary.buildNncFaceIndices(numElems);
        this->fluxBoundary_ = std::make_shared<FluxBoundary>(std::move(fluxBoundary));
        this->collectFluxBoundaryFaceAreas_(numElems);
        this->applyFluxBoundaryTransmissibilityOverrides_();
        this->reportUnmappedFluxBoundaryFaces_();

        this->fluxBoundaryInterpolate_ = Parameters::Get<Parameters::FluxBoundaryInterpolate>();

        OpmLog::info(fmt::format("USEFLUX: {} boundary record(s) read from '{}' ({}), {}",
                                 fluxData.reportSteps.size(),
                                 selectedFluxPath.string(),
                                 fluxData.header.boundaryPerTimestep
                                 ? "per time step"
                                 : "per report step",
                                 this->fluxBoundaryInterpolate_
                                 ? "interpolated in time between records"
                                 : "each held over the interval it covers"));

        // Start out on the record covering the initial time, so that the very
        // first assembly does not run without one.
        this->refreshFluxBoundaryRecord_(0.0, 0.0);

        if (!this->fluxBoundary_->faces().empty()) {
            this->nonTrivialBoundaryConditions_ = true;
        }

        this->loadFluxParentSummary_(fluxData, selectedFluxPath, inputDir, baseName);
    }

    //! \brief Establish where the reduced run reads the parent run's summary
    //!        vectors from.
    //! \details The vectors embedded in the .FLUX file are preferred because
    //!          they are sampled at sub-report-step resolution and their rate
    //!          entries are interval averages. Otherwise fall back to the
    //!          parent run's own summary output. Either way, say so in the log:
    //!          the two sources have different time resolution and the choice
    //!          affects the accuracy of the reduced run.
    void loadFluxParentSummary_(const EclIO::FluxFile::Data&  fluxData,
                                const std::filesystem::path&  fluxPath,
                                const std::filesystem::path&  inputDir,
                                const std::string&            baseName)
    {
        this->fluxParentSummaryData_.reset();

        if (auto embedded = ParentSummary::fromFluxFile(fluxData, fluxPath.string());
            embedded.has_value())
        {
            OpmLog::info(fmt::format("USEFLUX: summary vectors loaded from flux file '{}' "
                                     "({} keys, {} samples)",
                                     fluxPath.string(),
                                     embedded->keys().size(),
                                     embedded->times().size()));

            this->fluxParentSummaryData_ =
                std::make_shared<ParentSummary>(std::move(*embedded));
            return;
        }

        for (const auto* extension : {".ESMRY", ".SMSPEC", ".FSMSPEC"}) {
            const auto candidate = inputDir / (baseName + extension);

            if (auto external = ParentSummary::fromSummaryFile(candidate);
                external.has_value())
            {
                OpmLog::info(fmt::format("USEFLUX: summary vectors loaded from parent "
                                         "summary file '{}' ({} keys, {} samples)",
                                         candidate.string(),
                                         external->keys().size(),
                                         external->times().size()));

                this->fluxParentSummaryData_ =
                    std::make_shared<ParentSummary>(std::move(*external));
                return;
            }
        }

        OpmLog::warning("USEFLUX: no parent summary vectors available, neither embedded "
                        "in the flux file nor as a parent summary file next to it. "
                        "Field and group aggregates will only reflect the wells inside "
                        "the sector.");
    }

    // The FLUX file stores total volumetric rates across each boundary face,
    // while the linearizer multiplies the boundary rate vector by the face
    // area. Record the areas so the prescribed rates can be normalised.
    //! \brief Face areas and stencil ordinals of the sector boundary, per
    //!        (cell, direction).
    //!
    //! \details The ordinal is where the face sits in the cell's list of
    //!   boundary faces, which is the index the discretisation addresses it by.
    //!   It is NOT the direction index: the stencil appends a boundary face for
    //!   every intersection that has no neighbour, in iteration order, so a cell
    //!   whose I- neighbour lies inside the sector has its I+ face at ordinal
    //!   zero. Anything handing the discretisation a per-face quantity has to
    //!   go through here or it lands on the wrong face.
    void collectFluxBoundaryFaceAreas_(const std::size_t numElems)
    {
        this->fluxBoundaryFaceArea_.resize(numElems, 0.0);
        this->fluxBoundaryFaceOrdinal_.resize(numElems, -1);
        this->fluxBoundaryFaceDepth_.resize(numElems, 0.0);

        if (!this->fluxBoundary_) {
            return;
        }

        const auto& gridView = this->simulator().vanguard().gridView();
        const auto& elementMapper = this->model().elementMapper();

        for (const auto& elem : elements(gridView)) {
            const auto elemIdx = elementMapper.index(elem);

            int ordinal = 0;
            for (const auto& is : intersections(gridView, elem)) {
                // Counted exactly as EcfvStencil::updateTopology() does, which
                // means counting intersections without a neighbour rather than
                // boundary intersections: on a process boundary the two differ,
                // and it is the stencil's numbering the solver will use.
                if (is.neighbor()) {
                    continue;
                }

                const auto thisOrdinal = ordinal++;

                if (!is.boundary()) {
                    continue;
                }

                const auto dir = FaceDir::FromIntersectionIndex(is.indexInInside());
                auto& areas = this->fluxBoundaryFaceArea_(dir);
                if (static_cast<std::size_t>(elemIdx) < areas.size()) {
                    areas[elemIdx] = is.geometry().volume();
                    this->fluxBoundaryFaceOrdinal_(dir)[elemIdx] = thisOrdinal;

                    // Where the boundary condition is imposed. The pressures in
                    // the file belong to the cell on the far side, taken at its
                    // centre, and have to be carried here before they are used.
                    this->fluxBoundaryFaceDepth_(dir)[elemIdx] =
                        is.geometry().center()[dimWorld - 1];
                }
            }
        }
    }

    //! \brief Warn about boundary faces the .FLUX file describes but this grid
    //!        cannot carry.
    //!
    //! \details A face is applied by dividing the stored total rate by the area
    //!   of the matching grid intersection, so a face with no such intersection
    //!   silently drops whatever flux the producing run recorded for it. That
    //!   is a mass leak rather than a small inaccuracy, so it is worth saying
    //!   out loud.
    void reportUnmappedFluxBoundaryFaces_()
    {
        if (!this->fluxBoundary_) {
            return;
        }

        std::size_t unmapped = 0;
        std::size_t total = 0;
        for (const auto& face : this->fluxBoundary_->faces()) {
            if (face.isNnc || face.direction == FaceDir::Unknown) {
                continue;
            }

            // Faces belonging to another rank are not this one's to account
            // for. Counting them here would report most of the region as
            // unmapped on every rank of a parallel run.
            if (face.interiorActiveCell < 0) {
                continue;
            }

            // Nothing crosses a face with zero transmissibility, so it cannot leak.
            if ((face.transmissibility == 0.0) && (face.thermalTransmissibility == 0.0)) {
                continue;
            }

            ++total;
            const auto area =
                this->fluxBoundaryFaceAreaAt_(static_cast<unsigned>(face.interiorActiveCell),
                                              face.direction);
            if (!(area > 0.0)) {
                ++unmapped;
            }
        }

        if (unmapped > 0) {
            OpmLog::warning(fmt::format("USEFLUX: {} of {} boundary faces have no matching grid "
                                        "face in this run, so any flow the producing run recorded "
                                        "across them is lost. The sector grid and the region the "
                                        "FLUX file describes do not agree.",
                                        unmapped, total));
        }
    }

    Scalar fluxBoundaryFaceAreaAt_(const unsigned int globalSpaceIdx, const FaceDir::DirEnum dir) const
    {
        const auto& areas = this->fluxBoundaryFaceArea_(dir);
        if (globalSpaceIdx >= areas.size()) {
            return 0.0;
        }

        return areas[globalSpaceIdx];
    }

    //! \brief Which direction the discretisation's \p boundaryFaceIdx refers to
    //!        on this cell, or Unknown when it is not a sector boundary face.
    //!
    //! \details The inverse of the ordinal collected alongside the face areas.
    //!   Six comparisons rather than a table, because a cell may hold more
    //!   boundary faces than there are directions and only those the sector
    //!   boundary uses can be answered for.
    FaceDir::DirEnum fluxBoundaryDirAt_(const unsigned globalSpaceIdx,
                                        const unsigned boundaryFaceIdx) const
    {
        static constexpr auto directions = std::array {
            FaceDir::XMinus, FaceDir::XPlus,
            FaceDir::YMinus, FaceDir::YPlus,
            FaceDir::ZMinus, FaceDir::ZPlus,
        };

        for (const auto dir : directions) {
            const auto& ordinals = this->fluxBoundaryFaceOrdinal_(dir);
            if ((globalSpaceIdx < ordinals.size())
                && (ordinals[globalSpaceIdx] >= 0)
                && (static_cast<unsigned>(ordinals[globalSpaceIdx]) == boundaryFaceIdx))
            {
                return dir;
            }
        }

        return FaceDir::Unknown;
    }

    //! \brief Depth at which a sector boundary condition is imposed.
    //!
    //! \details The face the discretisation uses, which is not where the file's
    //!   pressures were measured: those belong to the cell on the far side,
    //!   taken at its centre.
    Scalar fluxBoundaryFaceDepthAt_(unsigned globalDofIdx, const FaceDir::DirEnum dir) const
    {
        const auto& depths = this->fluxBoundaryFaceDepth_(dir);

        return (globalDofIdx < depths.size())
            ? depths[globalDofIdx]
            : Scalar{0};
    }

    void applyFluxBoundaryTransmissibilityOverrides_()
    {
        if (!this->fluxBoundary_) {
            return;
        }

        // The discretisation addresses a boundary face by its position in the
        // cell's boundary-face list, not by its direction, so the override has
        // to be told where each face sits or it writes onto a different face.
        const auto ordinal = [this](const int activeCell, const FaceDir::DirEnum dir)
        {
            const auto& ordinals = this->fluxBoundaryFaceOrdinal_(dir);
            const auto cell = static_cast<std::size_t>(activeCell);

            return (activeCell >= 0) && (cell < ordinals.size()) ? ordinals[cell] : -1;
        };

        this->fluxBoundary_->applyTransmissibilityOverrides(transmissibilities_, ordinal);
        this->fluxBoundary_->applyThermalTransmissibilityOverrides(transmissibilities_, ordinal);
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
            const auto& iq = model.intensiveQuantities(elementIdx, /*timeIdx=*/ 0);
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
        auto obtain = [](const auto& value)
                      {
                          if constexpr (std::is_same_v<LhsEval, Scalar>) {
                              return getValue(value);
                          } else {
                              return value;
                          }
                      };

        return computeRockCompTransMultiplier_<LhsEval>(intQuants, elementIdx, obtain);
    }

    template <class LhsEval, class Callback>
    LhsEval computeRockCompTransMultiplier_(const IntensiveQuantities& intQuants, unsigned elementIdx, Callback& obtain) const
    {
        OPM_TIMEBLOCK_LOCAL(computeRockCompTransMultiplier, Subsystem::PvtProps);
        if (this->rockCompTransMult_.empty() && this->rockCompTransMultWc_.empty())
            return 1.0;

        unsigned tableIdx = 0;
        if (!this->rockTableIdx_.empty())
            // Auxiliary DOFs have no entry here; they take the first ROCK region.
            if (elementIdx < this->rockTableIdx_.size()) {
                tableIdx = this->rockTableIdx_[elementIdx];
            }

        const auto& fs = intQuants.fluidState();
        const auto& rock_config = this->simulator().vanguard().eclState().getSimulationConfig().rock_config();

        if (!this->rockCompTransMultElastic_.empty()) {
            // ROCKCOMP HYSTERESIS=HYSTER: see rockCompPoroMultiplier() above.
            LhsEval effectivePressure = obtain(fs.pressure(refPressurePhaseIdx_()));
            LhsEval turningPressure = this->minRefPressure_[elementIdx];

            if (!this->overburdenPressure_.empty()) {
                effectivePressure -= this->overburdenPressure_[elementIdx];
                turningPressure -= this->overburdenPressure_[elementIdx];
            }

            if (rock_config.store()) {
                const auto& initialPressure = asImp_().initialFluidState(elementIdx).pressure(refPressurePhaseIdx_());
                effectivePressure -= initialPressure;
                turningPressure -= initialPressure;
            }

            if (effectivePressure <= turningPressure)
                return this->rockCompTransMult_[tableIdx].eval(effectivePressure, /*extrapolation=*/true);

            return this->rockCompTransMultElastic_[tableIdx].eval(turningPressure, effectivePressure, /*extrapolation=*/true);
        }

        LhsEval effectivePressure = obtain(fs.pressure(refPressurePhaseIdx_()));
        if (!this->minRefPressure_.empty())
            // The pore space change is irreversible
            effectivePressure =
                min(obtain(fs.pressure(refPressurePhaseIdx_())),
                    this->minRefPressure_[elementIdx]);

        if (!this->overburdenPressure_.empty())
            effectivePressure -= this->overburdenPressure_[elementIdx];

        if (rock_config.store()) {
            effectivePressure -= asImp_().initialFluidState(elementIdx).pressure(refPressurePhaseIdx_());
        }

        if (!this->rockCompTransMult_.empty())
            return this->rockCompTransMult_[tableIdx].eval(effectivePressure, /*extrapolation=*/true);

        // water compaction
        assert(!this->rockCompTransMultWc_.empty());
        LhsEval SwMax = max(obtain(fs.saturation(waterPhaseIdx)), this->maxWaterSaturation_[elementIdx]);
        LhsEval SwDeltaMax = SwMax - asImp_().initialFluidStates()[elementIdx].saturation(waterPhaseIdx);

        return this->rockCompTransMultWc_[tableIdx].eval(effectivePressure, SwDeltaMax, /*extrapolation=*/true);
    }

    typename Vanguard::TransmissibilityType transmissibilities_;

    //! Owned here because their data feeds the problem's per-DOF tables.
    std::vector<std::unique_ptr<FlowAuxCellModule<TypeTag>>> auxCellModules_;

    std::shared_ptr<EclMaterialLawManager> materialLawManager_;
    std::shared_ptr<EclThermalLawManager> thermalLawManager_;

    GlobalEqVector drift_;

    WellModel wellModel_;
    AquiferModel aquiferModel_;

    PffGridVector<GridView, Stencil, PffDofData_, DofMapper> pffDofData_;
    TracerModel tracerModel_;
    TemperatureModel temperatureModel_;

    FlowThresholdPressure<TypeTag> thresholdPressures_;
    ModuleParams moduleParams_;

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

    virtual void handleSolventBC(const BCState::BCFace&, RateVector&) const = 0;

    virtual void handlePolymerBC(const BCState::BCFace&, RateVector&) const = 0;

    virtual void handleMicrBC(const BCState::BCFace&, RateVector&) const = 0;

    virtual void handleOxygBC(const BCState::BCFace&, RateVector&) const = 0;

    virtual void handleUreaBC(const BCState::BCFace&, RateVector&) const = 0;

    BCData<int> bcindex_;
    BCData<int> fluxBoundaryFaceIndex_;

    //! \brief Sector-boundary NNC faces per interior cell, indexed by active cell.
    std::vector<std::vector<int>> fluxNncFaceIndex_;
    BCData<Scalar> fluxBoundaryFaceArea_;
    BCData<int> fluxBoundaryFaceOrdinal_;

    //! \brief Depth of each sector boundary face, per (cell, direction).
    BCData<Scalar> fluxBoundaryFaceDepth_;
    std::shared_ptr<FluxBoundary> fluxBoundary_;
    std::shared_ptr<ParentSummary> fluxParentSummaryData_;
    const EclIO::FluxFile::ReportStep* fluxBoundaryActiveRecord_ = nullptr;
    EclIO::FluxFile::ReportStep fluxBoundaryBlendedRecord_;
    bool fluxBoundaryInterpolate_ = true;
    bool nonTrivialBoundaryConditions_ = false;
    bool first_step_ = true;
    bool enable_state_rollback_ = false;

    //! \brief Explicit quantities as they stood at the start of the current
    //!        timestep, restored if that timestep fails.
    //!
    //! Only written by captureBeginTimeStepState_() and only read by
    //! restoreBeginTimeStepState_(), both of which run inside a single
    //! timestep, so this never has to survive a restart.
    struct PrevTimestepState
    {
        bool first_step = true;                            //!< first step of the episode
        std::vector<Scalar> max_polymer_adsorption;        //!< POLYMER
        std::vector<Scalar> max_oil_saturation;            //!< VAPPARS / DRSDT saturation history
        std::vector<Scalar> max_water_saturation;          //!< ROCKCOMP water-induced compaction
        std::vector<Scalar> min_ref_pressure;              //!< ROCKCOMP irreversible compaction
        std::vector<Scalar> rock_comp_trans_mult_val;      //!< ROCKCOMP transmissibility multiplier
    };

    PrevTimestepState prev_timestep_state_;

    /// Whether or not the current episode will end at the end of the
    /// current time step.
    virtual bool episodeWillBeOver() const
    {
        return this->simulator().episodeWillBeOver();
    }
};

} // namespace Opm

#endif // OPM_FLOW_PROBLEM_HPP
