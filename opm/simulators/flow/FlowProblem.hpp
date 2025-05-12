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
#include <opm/simulators/flow/FlowBaseProblemProperties.hpp>
#include <opm/simulators/flow/FlowUtils.hpp>
#include <opm/simulators/flow/TracerModel.hpp>
#include <opm/simulators/flow/Transmissibility.hpp>
#include <opm/simulators/timestepping/AdaptiveTimeStepping.hpp>
#include <opm/simulators/timestepping/SimulatorReport.hpp>

#include <opm/simulators/utils/DeferredLoggingErrorHelpers.hpp>
#include <opm/simulators/utils/ParallelSerialization.hpp>
#include <opm/simulators/utils/satfunc/RelpermDiagnostics.hpp>

#include <opm/utility/CopyablePtr.hpp>

#include <algorithm>
#include <cstddef>
#include <functional>
#include <set>
#include <stdexcept>
#include <string>
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

    // Grid and world dimension
    enum { dim = GridView::dimension };
    enum { dimWorld = GridView::dimensionworld };

    // copy some indices for convenience
    enum { numEq = getPropValue<TypeTag, Properties::NumEq>() };
    enum { numPhases = FluidSystem::numPhases };
    enum { numComponents = FluidSystem::numComponents };

    enum { enableConvectiveMixing = getPropValue<TypeTag, Properties::EnableConvectiveMixing>() };
    enum { enableBrine = getPropValue<TypeTag, Properties::EnableBrine>() };
    enum { enableDiffusion = getPropValue<TypeTag, Properties::EnableDiffusion>() };
    enum { enableDispersion = getPropValue<TypeTag, Properties::EnableDispersion>() };
    enum { enableEnergy = getPropValue<TypeTag, Properties::EnableEnergy>() };
    enum { enableExperiments = getPropValue<TypeTag, Properties::EnableExperiments>() };
    enum { enableExtbo = getPropValue<TypeTag, Properties::EnableExtbo>() };
    enum { enableFoam = getPropValue<TypeTag, Properties::EnableFoam>() };
    enum { enableMICP = getPropValue<TypeTag, Properties::EnableMICP>() };
    enum { enablePolymer = getPropValue<TypeTag, Properties::EnablePolymer>() };
    enum { enablePolymerMolarWeight = getPropValue<TypeTag, Properties::EnablePolymerMW>() };
    enum { enableSaltPrecipitation = getPropValue<TypeTag, Properties::EnableSaltPrecipitation>() };
    enum { enableSolvent = getPropValue<TypeTag, Properties::EnableSolvent>() };
    enum { enableTemperature = getPropValue<TypeTag, Properties::EnableTemperature>() };
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
    using Indices = GetPropType<TypeTag, Properties::Indices>;
    using IntensiveQuantities = GetPropType<TypeTag, Properties::IntensiveQuantities>;
    using WellModel = GetPropType<TypeTag, Properties::WellModel>;
    using AquiferModel = GetPropType<TypeTag, Properties::AquiferModel>;

    using Toolbox = MathToolbox<Evaluation>;
    using DimMatrix = Dune::FieldMatrix<Scalar, dimWorld, dimWorld>;


    using TracerModel = GetPropType<TypeTag, Properties::TracerModel>;
    using DirectionalMobilityPtr = Utility::CopyablePtr<DirectionalMobility<TypeTag>>;

public:
    using BaseType::briefDescription;
    using BaseType::helpPreamble;
    using BaseType::shouldWriteOutput;
    using BaseType::shouldWriteRestartFile;
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
                              enableEnergy,
                              enableDiffusion,
                              enableDispersion)
        , wellModel_(simulator)
        , aquiferModel_(simulator)
        , pffDofData_(simulator.gridView(), this->elementMapper())
        , tracerModel_(simulator)
    {
        if (! Parameters::Get<Parameters::CheckSatfuncConsistency>()) {
            // User did not enable the "new" saturation function consistency
            // check module.  Run the original checker instead.  This is a
            // temporary measure.
            RelpermDiagnostics relpermDiagnostics{};
            relpermDiagnostics.diagnosis(simulator.vanguard().eclState(),
                                         simulator.vanguard().levelCartesianIndexMapper());
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
            updatePffDofData_();
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
    void beginTimeStep()
    {
        OPM_TIMEBLOCK(beginTimeStep);
        const int episodeIdx = this->episodeIndex();
        const int timeStepSize = this->simulator().timeStepSize();

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

        if (this->enableDriftCompensation_) {
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
        // use the generic code to prepare the output fields and to
        // write the desired VTK files.
        if (Parameters::Get<Parameters::EnableWriteAllSolutions>() ||
            this->simulator().episodeWillBeOver()) {
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
                tableIdx = this->rockTableIdx_[globalSpaceIdx];
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
            dirMob = std::make_unique<DirectionalMobility<TypeTag>>();
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
        return asImp_().initialFluidState(globalDofIdx).temperature(/*phaseIdx=*/0);
    }


    Scalar temperature(unsigned globalDofIdx, unsigned /*timeIdx*/) const
    {
        // use the initial temperature of the DOF if temperature is not a primary
        // variable
         return asImp_().initialFluidState(globalDofIdx).temperature(/*phaseIdx=*/0);
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
        OPM_TIMEBLOCK_LOCAL(rockCompPoroMultiplier);
        if (this->rockCompPoroMult_.empty() && this->rockCompPoroMultWc_.empty())
            return 1.0;

        unsigned tableIdx = 0;
        if (!this->rockTableIdx_.empty())
            tableIdx = this->rockTableIdx_[elementIdx];

        const auto& fs = intQuants.fluidState();
        LhsEval effectivePressure = decay<LhsEval>(fs.pressure(refPressurePhaseIdx_()));
        const auto& rock_config = this->simulator().vanguard().eclState().getSimulationConfig().rock_config();
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
        const bool implicit = !this->explicitRockCompaction_;
        return implicit ? this->simulator().problem().template computeRockCompTransMultiplier_<LhsEval>(intQuants, elementIdx)
                        : this->simulator().problem().getRockCompTransMultVal(elementIdx);
    }


    /*!
     * \brief Return the well transmissibility multiplier due to rock changes.
     */
    template <class LhsEval>
    LhsEval wellTransMultiplier(const IntensiveQuantities& intQuants, unsigned elementIdx) const
    {
        OPM_TIMEBLOCK_LOCAL(wellTransMultiplier);
        
        const bool implicit = !this->explicitRockCompaction_;
        double trans_mult = implicit ? this->simulator().problem().template computeRockCompTransMultiplier_<double>(intQuants, elementIdx)
                                     : this->simulator().problem().getRockCompTransMultVal(elementIdx);
        trans_mult *= this->simulator().problem().template permFactTransMultiplier<double>(intQuants, elementIdx);
    
        return trans_mult;
    }

    std::pair<BCType, RateVector> boundaryCondition(const unsigned int globalSpaceIdx, const int directionId) const
    {
        OPM_TIMEBLOCK_LOCAL(boundaryCondition);
        if (!nonTrivialBoundaryConditions_) {
            return { BCType::NONE, RateVector(0.0) };
        }
        FaceDir::DirEnum dir = FaceDir::FromIntersectionIndex(directionId);
        const auto& schedule = this->simulator().vanguard().schedule();
        if (bcindex_(dir)[globalSpaceIdx] == 0) {
            return { BCType::NONE, RateVector(0.0) };
        }
        if (schedule[this->episodeIndex()].bcprop.size() == 0) {
            return { BCType::NONE, RateVector(0.0) };
        }
        const auto& bc = schedule[this->episodeIndex()].bcprop[bcindex_(dir)[globalSpaceIdx]];
        if (bc.bctype!=BCType::RATE) {
            return { bc.bctype, RateVector(0.0) };
        }

        RateVector rate = 0.0;
        switch (bc.component) {
        case BCComponent::OIL:
            rate[Indices::canonicalToActiveComponentIndex(oilCompIdx)] = bc.rate;
            break;
        case BCComponent::GAS:
            rate[Indices::canonicalToActiveComponentIndex(gasCompIdx)] = bc.rate;
            break;
        case BCComponent::WATER:
            rate[Indices::canonicalToActiveComponentIndex(waterCompIdx)] = bc.rate;
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


    template<class Serializer>
    void serializeOp(Serializer& serializer)
    {
        serializer(static_cast<BaseType&>(*this));
        serializer(drift_);
        serializer(wellModel_);
        serializer(aquiferModel_);
        serializer(tracerModel_);
        serializer(*materialLawManager_);
    }

private:
    Implementation& asImp_()
    { return *static_cast<Implementation *>(this); }

    const Implementation& asImp_() const
    { return *static_cast<const Implementation *>(this); }

protected:
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
        this->updateProperty_("FlowProblem::updateMaxWaterSaturation_() failed:",
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

        this->updateProperty_("FlowProblem::updateMinPressure_() failed:",
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

protected:
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

    virtual void updateExplicitQuantities_(int episodeIdx, int timeStepSize, bool first_step_after_restart) = 0;

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

        unsigned tableIdx = 0;
        if (!this->rockTableIdx_.empty())
            tableIdx = this->rockTableIdx_[elementIdx];

        const auto& fs = intQuants.fluidState();
        LhsEval effectivePressure = decay<LhsEval>(fs.pressure(refPressurePhaseIdx_()));
        const auto& rock_config = this->simulator().vanguard().eclState().getSimulationConfig().rock_config();
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

        if (!this->rockCompTransMult_.empty())
            return this->rockCompTransMult_[tableIdx].eval(effectivePressure, /*extrapolation=*/true);

        // water compaction
        assert(!this->rockCompTransMultWc_.empty());
        LhsEval SwMax = max(decay<LhsEval>(fs.saturation(waterPhaseIdx)), this->maxWaterSaturation_[elementIdx]);
        LhsEval SwDeltaMax = SwMax - asImp_().initialFluidStates()[elementIdx].saturation(waterPhaseIdx);

        return this->rockCompTransMultWc_[tableIdx].eval(effectivePressure, SwDeltaMax, /*extrapolation=*/true);
    }

    typename Vanguard::TransmissibilityType transmissibilities_;

    std::shared_ptr<EclMaterialLawManager> materialLawManager_;
    std::shared_ptr<EclThermalLawManager> thermalLawManager_;

    GlobalEqVector drift_;

    WellModel wellModel_;
    AquiferModel aquiferModel_;

    PffGridVector<GridView, Stencil, PffDofData_, DofMapper> pffDofData_;
    TracerModel tracerModel_;

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

    virtual void handleSolventBC(const BCProp::BCFace&, RateVector&) const = 0;

    virtual void handlePolymerBC(const BCProp::BCFace&, RateVector&) const = 0;

    virtual void handleMicrBC(const BCProp::BCFace&, RateVector&) const = 0;

    virtual void handleOxygBC(const BCProp::BCFace&, RateVector&) const = 0;

    virtual void handleUreaBC(const BCProp::BCFace&, RateVector&) const = 0;

    BCData<int> bcindex_;
    bool nonTrivialBoundaryConditions_ = false;
    bool first_step_ = true;
};

} // namespace Opm

#endif // OPM_FLOW_PROBLEM_HPP
