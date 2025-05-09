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
/**
 * \file
 *
 * \copydoc Opm::TracerModel
 */
#ifndef OPM_TRACER_MODEL_HPP
#define OPM_TRACER_MODEL_HPP

#include <opm/common/OpmLog/OpmLog.hpp>
#include <opm/common/TimingMacros.hpp>

#include <opm/input/eclipse/Schedule/Well/WellConnections.hpp>

#include <opm/grid/utility/ElementChunks.hpp>

#include <opm/models/parallel/threadmanager.hpp>
#include <opm/models/utils/propertysystem.hh>

#include <opm/simulators/flow/GenericTracerModel.hpp>
#include <opm/simulators/utils/DeferredLoggingErrorHelpers.hpp>
#include <opm/simulators/utils/VectorVectorDataHandle.hpp>

#include <array>
#include <cstddef>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace Opm::Properties {

template<class TypeTag, class MyTypeTag>
struct EnableTracerModel {
    using type = UndefinedProperty;
};

} // namespace Opm::Properties

namespace Opm {

/*!
 * \ingroup BlackOilSimulator
 *
 * \brief A class which handles tracers as specified in by ECL
 */
template <class TypeTag>
class TracerModel : public GenericTracerModel<GetPropType<TypeTag, Properties::Grid>,
                                              GetPropType<TypeTag, Properties::GridView>,
                                              GetPropType<TypeTag, Properties::DofMapper>,
                                              GetPropType<TypeTag, Properties::Stencil>,
                                              GetPropType<TypeTag, Properties::FluidSystem>,
                                              GetPropType<TypeTag, Properties::Scalar>>
{
    using BaseType = GenericTracerModel<GetPropType<TypeTag, Properties::Grid>,
                                        GetPropType<TypeTag, Properties::GridView>,
                                        GetPropType<TypeTag, Properties::DofMapper>,
                                        GetPropType<TypeTag, Properties::Stencil>,
                                        GetPropType<TypeTag, Properties::FluidSystem>,
                                        GetPropType<TypeTag, Properties::Scalar>>;
    using Simulator = GetPropType<TypeTag, Properties::Simulator>;
    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using Grid = GetPropType<TypeTag, Properties::Grid>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Stencil = GetPropType<TypeTag, Properties::Stencil>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using ElementContext = GetPropType<TypeTag, Properties::ElementContext>;
    using RateVector = GetPropType<TypeTag, Properties::RateVector>;
    using Indices = GetPropType<TypeTag, Properties::Indices>;

    using TracerEvaluation = DenseAd::Evaluation<Scalar,1>;

    using TracerMatrix = typename BaseType::TracerMatrix;
    using TracerVector = typename BaseType::TracerVector;
    using TracerVectorSingle = typename BaseType::TracerVectorSingle;

    enum { numEq = getPropValue<TypeTag, Properties::NumEq>() };
    enum { numPhases = FluidSystem::numPhases };
    enum { waterPhaseIdx = FluidSystem::waterPhaseIdx };
    enum { oilPhaseIdx = FluidSystem::oilPhaseIdx };
    enum { gasPhaseIdx = FluidSystem::gasPhaseIdx };

public:
    explicit TracerModel(Simulator& simulator)
        : BaseType(simulator.vanguard().gridView(),
                   simulator.vanguard().eclState(),
                   simulator.vanguard().cartesianIndexMapper(),
                   simulator.model().dofMapper(),
                   simulator.vanguard().cellCentroids())
        , simulator_(simulator)
        , tbatch({waterPhaseIdx, oilPhaseIdx, gasPhaseIdx})
        , wat_(tbatch[0])
        , oil_(tbatch[1])
        , gas_(tbatch[2])
        , element_chunks_(simulator.gridView(), Dune::Partitions::all, ThreadManager::maxThreads())
    { }


    /*
      The initialization of the tracer model is a three step process:

      1. The init() method is called. This will allocate buffers and initialize
         some phase index stuff. If this is a normal run the initial tracer
         concentrations will be assigned from the TBLK or TVDPF keywords.

      2. [Restart only:] The tracer concentration are read from the restart
         file and the concentrations are applied with repeated calls to the
         setTracerConcentration() method. This is currently done in the
         eclwriter::beginRestart() method.

      3. Internally the tracer model manages the concentrations in "batches" for
         the oil, water and gas tracers respectively. The batches should be
         initialized with the initial concentration, that must be performed
         after the concentration values have been assigned. This is done in
         method prepareTracerBatches() called from eclproblem::finishInit().
    */
    void init(bool rst)
    {
        this->doInit(rst, simulator_.model().numGridDof(),
                     gasPhaseIdx, oilPhaseIdx, waterPhaseIdx);
    }

    void prepareTracerBatches()
    {
        for (std::size_t tracerIdx = 0; tracerIdx < this->tracerPhaseIdx_.size(); ++tracerIdx) {
            if (this->tracerPhaseIdx_[tracerIdx] == FluidSystem::waterPhaseIdx) {
                if (! FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)){
                    throw std::runtime_error("Water tracer specified for non-water fluid system: " +
                                             this->name(tracerIdx));
                }

                wat_.addTracer(tracerIdx, this->tracerConcentration_[tracerIdx]);
            }
            else if (this->tracerPhaseIdx_[tracerIdx] == FluidSystem::oilPhaseIdx) {
                if (! FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx)){
                    throw std::runtime_error("Oil tracer specified for non-oil fluid system: " +
                                             this->name(tracerIdx));
                }

                oil_.addTracer(tracerIdx, this->tracerConcentration_[tracerIdx]);
            }
            else if (this->tracerPhaseIdx_[tracerIdx] == FluidSystem::gasPhaseIdx) {
                if (! FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)){
                    throw std::runtime_error("Gas tracer specified for non-gas fluid system: " +
                                             this->name(tracerIdx));
                }

                gas_.addTracer(tracerIdx, this->tracerConcentration_[tracerIdx]);
            }

            // resize free and solution volume storages
            vol1_[0][this->tracerPhaseIdx_[tracerIdx]].
                resize(this->freeTracerConcentration_[tracerIdx].size());
            vol1_[1][this->tracerPhaseIdx_[tracerIdx]].
                resize(this->freeTracerConcentration_[tracerIdx].size());
            dVol_[0][this->tracerPhaseIdx_[tracerIdx]].
                resize(this->solTracerConcentration_[tracerIdx].size());
            dVol_[1][this->tracerPhaseIdx_[tracerIdx]].
                resize(this->solTracerConcentration_[tracerIdx].size());
        }

        // will be valid after we move out of tracerMatrix_
        TracerMatrix* base = this->tracerMatrix_.get();
        for (auto& tr : this->tbatch) {
            if (tr.numTracer() != 0) {
                if (this->tracerMatrix_) {
                    tr.mat = std::move(this->tracerMatrix_);
                }
                else {
                    tr.mat = std::make_unique<TracerMatrix>(*base);
                }
            }
        }
    }

    void beginTimeStep()
    {
        if (this->numTracers() == 0) {
            return;
        }

        OPM_TIMEBLOCK(tracerUpdateCache);
        updateStorageCache();
    }

    /*!
     * \brief Informs the tracer model that a time step has just been finished.
     */
    void endTimeStep()
    {
        if (this->numTracers() == 0) {
            return;
        }

        OPM_TIMEBLOCK(tracerAdvance);
        advanceTracerFields();
    }

    /*!
     * \brief This method writes the complete state of all tracer
     *        to the hard disk.
     */
    template <class Restarter>
    void serialize(Restarter&)
    { /* not implemented */ }

    /*!
     * \brief This method restores the complete state of the tracer
     *        from disk.
     *
     * It is the inverse of the serialize() method.
     */
    template <class Restarter>
    void deserialize(Restarter&)
    { /* not implemented */ }

    template<class Serializer>
    void serializeOp(Serializer& serializer)
    {
        serializer(static_cast<BaseType&>(*this));
        serializer(tbatch);
    }

protected:
    using TracerTypeIdx = typename BaseType::TracerTypeIdx;
    using BaseType::Free;
    using BaseType::Solution;

    // compute volume associated with free/solution concentration
    template<TracerTypeIdx Index>
    Scalar computeVolume_(const int tracerPhaseIdx,
                          const unsigned globalDofIdx,
                          const unsigned timeIdx) const
    {
        const auto& intQuants = simulator_.model().intensiveQuantities(globalDofIdx, timeIdx);
        const auto& fs = intQuants.fluidState();
        constexpr Scalar min_volume = 1e-10;

        if constexpr (Index == Free) {
            return std::max(decay<Scalar>(fs.saturation(tracerPhaseIdx)) *
                            decay<Scalar>(fs.invB(tracerPhaseIdx)) *
                            decay<Scalar>(intQuants.porosity()),
                            min_volume);
        } else {
            // vaporized oil
            if (tracerPhaseIdx == FluidSystem::oilPhaseIdx && FluidSystem::enableVaporizedOil()) {
                return std::max(decay<Scalar>(fs.saturation(FluidSystem::gasPhaseIdx)) *
                                decay<Scalar>(fs.invB(FluidSystem::gasPhaseIdx)) *
                                decay<Scalar>(fs.Rv()) *
                                decay<Scalar>(intQuants.porosity()),
                                min_volume);
            }

            // dissolved gas
            else if (tracerPhaseIdx == FluidSystem::gasPhaseIdx && FluidSystem::enableDissolvedGas()) {
                return std::max(decay<Scalar>(fs.saturation(FluidSystem::oilPhaseIdx)) *
                                decay<Scalar>(fs.invB(FluidSystem::oilPhaseIdx)) *
                                decay<Scalar>(fs.Rs()) *
                                decay<Scalar>(intQuants.porosity()),
                                min_volume);
            }

            return min_volume;
        }
    }

    template<TracerTypeIdx Index>
    std::pair<TracerEvaluation, bool>
    computeFlux_(const int tracerPhaseIdx,
                 const ElementContext& elemCtx,
                 const unsigned scvfIdx,
                 const unsigned timeIdx) const
    {
        const auto& stencil = elemCtx.stencil(timeIdx);
        const auto& scvf = stencil.interiorFace(scvfIdx);

        const auto& extQuants = elemCtx.extensiveQuantities(scvfIdx, timeIdx);
        const unsigned inIdx = extQuants.interiorIndex();

        Scalar v;
        unsigned upIdx;

        if constexpr (Index == Free) {
            upIdx = extQuants.upstreamIndex(tracerPhaseIdx);
            const auto& intQuants = elemCtx.intensiveQuantities(upIdx, timeIdx);
            const auto& fs = intQuants.fluidState();
            v = decay<Scalar>(extQuants.volumeFlux(tracerPhaseIdx)) *
                decay<Scalar>(fs.invB(tracerPhaseIdx));
        } else {
            if (tracerPhaseIdx == FluidSystem::oilPhaseIdx && FluidSystem::enableVaporizedOil()) {
                upIdx = extQuants.upstreamIndex(FluidSystem::gasPhaseIdx);

                const auto& intQuants = elemCtx.intensiveQuantities(upIdx, timeIdx);
                const auto& fs = intQuants.fluidState();
                v = decay<Scalar>(fs.invB(FluidSystem::gasPhaseIdx)) *
                    decay<Scalar>(extQuants.volumeFlux(FluidSystem::gasPhaseIdx)) *
                    decay<Scalar>(fs.Rv());
            }
            // dissolved gas
            else if (tracerPhaseIdx == FluidSystem::gasPhaseIdx && FluidSystem::enableDissolvedGas()) {
                upIdx = extQuants.upstreamIndex(FluidSystem::oilPhaseIdx);

                const auto& intQuants = elemCtx.intensiveQuantities(upIdx, timeIdx);
                const auto& fs = intQuants.fluidState();
                v = decay<Scalar>(fs.invB(FluidSystem::oilPhaseIdx)) *
                    decay<Scalar>(extQuants.volumeFlux(FluidSystem::oilPhaseIdx)) *
                    decay<Scalar>(fs.Rs());
            }
            else {
                upIdx = 0;
                v = 0.0;
            }
        }

        const Scalar A = scvf.area();
        return inIdx == upIdx
            ? std::pair{A * v * variable<TracerEvaluation>(1.0, 0), true}
            : std::pair{A * v, false};
    }

    template<TracerTypeIdx Index, class TrRe>
    Scalar storage1_(const TrRe& tr,
                     const unsigned tIdx,
                     const unsigned I,
                     const unsigned I1,
                     const bool cache)
    {
        if (cache) {
            return tr.storageOfTimeIndex1_[tIdx][I][Index];
        } else {
            return computeVolume_<Index>(tr.phaseIdx_, I1, 1) *
                   tr.concentration_[tIdx][I1][Index];
        }
    }

    template<class TrRe>
    void assembleTracerEquationVolume(TrRe& tr,
                                      const ElementContext& elemCtx,
                                      const Scalar scvVolume,
                                      const Scalar dt,
                                      unsigned I,
                                      unsigned I1)

    {
        if (tr.numTracer() == 0) {
            return;
        }

        const TracerEvaluation fVol = computeVolume_<Free>(tr.phaseIdx_, I, 0) * variable<TracerEvaluation>(1.0, 0);
        const TracerEvaluation sVol = computeVolume_<Solution>(tr.phaseIdx_, I, 0) * variable<TracerEvaluation>(1.0, 0);
        dVol_[Solution][tr.phaseIdx_][I] += sVol.value() * scvVolume - vol1_[1][tr.phaseIdx_][I];
        dVol_[Free][tr.phaseIdx_][I] += fVol.value() * scvVolume - vol1_[0][tr.phaseIdx_][I];
        for (int tIdx = 0; tIdx < tr.numTracer(); ++tIdx) {
            // Free part
            const Scalar fStorageOfTimeIndex0 = fVol.value() * tr.concentration_[tIdx][I][Free];
            const Scalar fLocalStorage = (fStorageOfTimeIndex0 - storage1_<Free>(tr, tIdx, I, I1,
                                                                                 elemCtx.enableStorageCache())) * scvVolume / dt;
            tr.residual_[tIdx][I][Free] += fLocalStorage; // residual + flux

            // Solution part
            const Scalar sStorageOfTimeIndex0 = sVol.value() * tr.concentration_[tIdx][I][Solution];
            const Scalar sLocalStorage = (sStorageOfTimeIndex0 - storage1_<Solution>(tr, tIdx, I, I1,
                                                                                     elemCtx.enableStorageCache())) * scvVolume / dt;
            tr.residual_[tIdx][I][Solution] += sLocalStorage; // residual + flux
        }

        // Derivative matrix
        (*tr.mat)[I][I][Free][Free] += fVol.derivative(0) * scvVolume/dt;
        (*tr.mat)[I][I][Solution][Solution] += sVol.derivative(0) * scvVolume/dt;
    }

    template<class TrRe>
    void assembleTracerEquationFlux(TrRe& tr,
                                    const ElementContext& elemCtx,
                                    unsigned scvfIdx,
                                    unsigned I,
                                    unsigned J,
                                    const Scalar dt)
    {
        if (tr.numTracer() == 0) {
            return;
        }

        const auto& [fFlux, isUpF] = computeFlux_<Free>(tr.phaseIdx_, elemCtx, scvfIdx, 0);
        const auto& [sFlux, isUpS] = computeFlux_<Solution>(tr.phaseIdx_, elemCtx, scvfIdx, 0);
        dVol_[Solution][tr.phaseIdx_][I] += sFlux.value() * dt;
        dVol_[Free][tr.phaseIdx_][I] += fFlux.value() * dt;
        const int fGlobalUpIdx = isUpF ? I : J;
        const int sGlobalUpIdx = isUpS ? I : J;
        for (int tIdx = 0; tIdx < tr.numTracer(); ++tIdx) {
            // Free and solution fluxes
            tr.residual_[tIdx][I][Free] += fFlux.value()*tr.concentration_[tIdx][fGlobalUpIdx][Free]; // residual + flux
            tr.residual_[tIdx][I][Solution] += sFlux.value()*tr.concentration_[tIdx][sGlobalUpIdx][Solution]; // residual + flux
        }

        // Derivative matrix
        if (isUpF){
            (*tr.mat)[J][I][Free][Free] = -fFlux.derivative(0);
            (*tr.mat)[I][I][Free][Free] += fFlux.derivative(0);
        }
        if (isUpS) {
            (*tr.mat)[J][I][Solution][Solution] = -sFlux.derivative(0);
            (*tr.mat)[I][I][Solution][Solution] += sFlux.derivative(0);
        }
    }

    template<class TrRe, class Well>
    void assembleTracerEquationWell(TrRe& tr,
                                    const Well& well)
    {
        if (tr.numTracer() == 0) {
            return;
        }

        const auto& eclWell = well.wellEcl();

        // Init. well output to zero
        auto& tracerRate = this->wellTracerRate_[eclWell.seqIndex()];
        tracerRate.reserve(tr.numTracer());
        auto& solTracerRate = this->wellSolTracerRate_[eclWell.seqIndex()];
        solTracerRate.reserve(tr.numTracer());
        auto& freeTracerRate = this->wellFreeTracerRate_[eclWell.seqIndex()];
        freeTracerRate.reserve(tr.numTracer());
        auto* mswTracerRate = eclWell.isMultiSegment()
            ? &this->mSwTracerRate_[eclWell.seqIndex()]
            : nullptr;
        if (mswTracerRate) {
            mswTracerRate->reserve(tr.numTracer());
        }
        for (int tIdx = 0; tIdx < tr.numTracer(); ++tIdx) {
            tracerRate.emplace_back(this->name(tr.idx_[tIdx]), 0.0);
            freeTracerRate.emplace_back(this->wellfname(tr.idx_[tIdx]), 0.0);
            solTracerRate.emplace_back(this->wellsname(tr.idx_[tIdx]), 0.0);
            if (eclWell.isMultiSegment()) {
                auto& wtr = mswTracerRate->emplace_back(this->name(tr.idx_[tIdx]));
                wtr.rate.reserve(eclWell.getConnections().size());
                for (std::size_t i = 0; i < eclWell.getConnections().size(); ++i) {
                    wtr.rate.emplace(eclWell.getConnections().get(i).segment(), 0.0);
                }
            }
        }

        std::vector<Scalar> wtracer(tr.numTracer());
        for (int tIdx = 0; tIdx < tr.numTracer(); ++tIdx) {
            wtracer[tIdx] = this->currentConcentration_(eclWell, this->name(tr.idx_[tIdx]));
        }

        const Scalar dt = simulator_.timeStepSize();
        const auto& ws = simulator_.problem().wellModel().wellState().well(well.name());
        for (std::size_t i = 0; i < ws.perf_data.size(); ++i) {
            const auto I = ws.perf_data.cell_index[i];
            const Scalar rate = well.volumetricSurfaceRateForConnection(I, tr.phaseIdx_);
            Scalar rate_s;
            if (tr.phaseIdx_ == FluidSystem::oilPhaseIdx && FluidSystem::enableVaporizedOil()) {
                rate_s = ws.perf_data.phase_mixing_rates[i][ws.vaporized_oil];
            }
            else if (tr.phaseIdx_ == FluidSystem::gasPhaseIdx && FluidSystem::enableDissolvedGas()) {
                rate_s = ws.perf_data.phase_mixing_rates[i][ws.dissolved_gas];
            }
            else {
                rate_s = 0.0;
            }

            const Scalar rate_f = rate - rate_s;
            if (rate_f > 0) {
                for (int tIdx = 0; tIdx < tr.numTracer(); ++tIdx) {
                    const Scalar delta = rate_f * wtracer[tIdx];
                    // Injection of free tracer only
                    tr.residual_[tIdx][I][Free] -= delta;

                    // Store _injector_ tracer rate for reporting
                    // (can be done here since WTRACER is constant)
                    tracerRate[tIdx].rate += delta;
                    freeTracerRate[tIdx].rate += delta;
                    if (eclWell.isMultiSegment()) {
                        (*mswTracerRate)[tIdx].rate[eclWell.getConnections().get(i).segment()] += delta;
                    }
                }
                dVol_[Free][tr.phaseIdx_][I] -= rate_f * dt;
            }
            else if (rate_f < 0) {
                for (int tIdx = 0; tIdx < tr.numTracer(); ++tIdx) {
                    const Scalar delta = rate_f * wtracer[tIdx];
                    // Store _injector_ tracer rate for cross-flowing well connections
                    // (can be done here since WTRACER is constant)
                    tracerRate[tIdx].rate += delta;
                    freeTracerRate[tIdx].rate += delta;

                    // Production of free tracer
                    tr.residual_[tIdx][I][Free] -= rate_f * tr.concentration_[tIdx][I][Free];
                }
                dVol_[Free][tr.phaseIdx_][I] -= rate_f * dt;

                // Derivative matrix for free tracer producer
                (*tr.mat)[I][I][Free][Free] -= rate_f * variable<TracerEvaluation>(1.0, 0).derivative(0);
            }
            if (rate_s < 0) {
                for (int tIdx = 0; tIdx < tr.numTracer(); ++tIdx) {
                    // Production of solution tracer
                    tr.residual_[tIdx][I][Solution] -= rate_s * tr.concentration_[tIdx][I][Solution];
                }
                dVol_[Solution][tr.phaseIdx_][I] -= rate_s * dt;

                // Derivative matrix for solution tracer producer
                (*tr.mat)[I][I][Solution][Solution] -= rate_s * variable<TracerEvaluation>(1.0, 0).derivative(0);
            }
        }
    }

    template<class TrRe>
    void assembleTracerEquationSource(TrRe& tr,
                                      const Scalar dt,
                                      unsigned I)
    {
        if (tr.numTracer() == 0) {
            return;
        }

        // Skip if solution tracers do not exist
        if (tr.phaseIdx_ ==  FluidSystem::waterPhaseIdx ||
            (tr.phaseIdx_ ==  FluidSystem::gasPhaseIdx && !FluidSystem::enableDissolvedGas()) ||
            (tr.phaseIdx_ ==  FluidSystem::oilPhaseIdx && !FluidSystem::enableVaporizedOil()))
        {
            return;
        }

        const Scalar& dsVol = dVol_[Solution][tr.phaseIdx_][I];
        const Scalar& dfVol = dVol_[Free][tr.phaseIdx_][I];

        // Source term determined by sign of dsVol: if dsVol > 0 then ms -> mf, else mf -> ms
        for (int tIdx = 0; tIdx < tr.numTracer(); ++tIdx) {
            if (dsVol >= 0) {
                const auto delta = (dfVol / dt) * tr.concentration_[tIdx][I][Free];
                tr.residual_[tIdx][I][Free] -= delta;
                tr.residual_[tIdx][I][Solution] += delta;
            }
            else {
                const auto delta = (dsVol / dt) * tr.concentration_[tIdx][I][Solution];
                tr.residual_[tIdx][I][Free] += delta;
                tr.residual_[tIdx][I][Solution] -= delta;
            }
        }

        // Derivative matrix
        if (dsVol >= 0) {
            const auto delta = (dfVol / dt) * variable<TracerEvaluation>(1.0, 0).derivative(0);
            (*tr.mat)[I][I][Free][Free] -= delta;
            (*tr.mat)[I][I][Solution][Free] += delta;
        }
        else {
            const auto delta = (dsVol / dt) * variable<TracerEvaluation>(1.0, 0).derivative(0);
            (*tr.mat)[I][I][Free][Solution] += delta;
            (*tr.mat)[I][I][Solution][Solution] -= delta;
        }
    }

    void assembleTracerEquations_()
    {
        // Note that we formulate the equations in terms of a concentration update
        // (compared to previous time step) and not absolute concentration.
        // This implies that current concentration (tr.concentration_[][]) contributes
        // to the rhs both through storage and flux terms.
        // Compare also advanceTracerFields(...) below.

        DeferredLogger local_deferredLogger{};
        OPM_BEGIN_PARALLEL_TRY_CATCH()
        {
            OPM_TIMEBLOCK(tracerAssemble);
            for (auto& tr : tbatch) {
                if (tr.numTracer() != 0) {
                    (*tr.mat) = 0.0;
                    for (int tIdx = 0; tIdx < tr.numTracer(); ++tIdx) {
                        tr.residual_[tIdx] = 0.0;
                    }
                }
            }

            this->wellTracerRate_.clear();
            this->wellFreeTracerRate_.clear();
            this->wellSolTracerRate_.clear();

             // educated guess for new container size
            const auto num_msw = this->mSwTracerRate_.size();
            this->mSwTracerRate_.clear();

            // Well terms
            const auto& wellPtrs = simulator_.problem().wellModel().localNonshutWells();
            this->wellTracerRate_.reserve(wellPtrs.size());
            this->wellFreeTracerRate_.reserve(wellPtrs.size());
            this->wellSolTracerRate_.reserve(wellPtrs.size());
            this->mSwTracerRate_.reserve(num_msw);
            for (const auto& wellPtr : wellPtrs) {
                for (auto& tr : tbatch) {
                    this->assembleTracerEquationWell(tr, *wellPtr);
                }
            }

            // Parallel loop over element chunks
            #ifdef _OPENMP
            #pragma omp parallel for
            #endif
            for (const auto& chunk : element_chunks_) {
                ElementContext elemCtx(simulator_);
                const Scalar dt = elemCtx.simulator().timeStepSize();

                for (const auto& elem : chunk) {
                    elemCtx.updateStencil(elem);
                    const std::size_t I = elemCtx.globalSpaceIndex(/*dofIdx=*/ 0, /*timeIdx=*/0);

                    if (elem.partitionType() != Dune::InteriorEntity) {
                        // Dirichlet boundary conditions for parallel matrix
                        // This is safe as each element has a unique I. So each thread
                        // always writes to different memory locations in the shared arrays.
                        for (const auto& tr : tbatch) {
                            if (tr.numTracer() != 0) {
                                (*tr.mat)[I][I][0][0] = 1.;
                                (*tr.mat)[I][I][1][1] = 1.;
                            }
                        }
                        continue;
                    }
                    elemCtx.updateAllIntensiveQuantities();
                    elemCtx.updateAllExtensiveQuantities();

                    const Scalar extrusionFactor =
                        elemCtx.intensiveQuantities(/*dofIdx=*/ 0, /*timeIdx=*/0).extrusionFactor();
                    Valgrind::CheckDefined(extrusionFactor);
                    assert(isfinite(extrusionFactor));
                    assert(extrusionFactor > 0.0);
                    const Scalar scvVolume =
                        elemCtx.stencil(/*timeIdx=*/0).subControlVolume(/*dofIdx=*/ 0).volume() * extrusionFactor;
                    const std::size_t I1 = elemCtx.globalSpaceIndex(/*dofIdx=*/ 0, /*timeIdx=*/1);

                    // This is safe as each element has a unique I. So each thread
                    // always writes to different memory locations in the shared arrays.
                    for (auto& tr : tbatch) {
                        if (tr.numTracer() == 0) {
                            continue;
                        }
                        this->assembleTracerEquationVolume(tr, elemCtx, scvVolume, dt, I, I1);
                    }

                    const std::size_t numInteriorFaces = elemCtx.numInteriorFaces(/*timIdx=*/0);
                    for (unsigned scvfIdx = 0; scvfIdx < numInteriorFaces; scvfIdx++) {
                        const auto& face = elemCtx.stencil(0).interiorFace(scvfIdx);
                        const unsigned j = face.exteriorIndex();
                        const unsigned J = elemCtx.globalSpaceIndex(/*dofIdx=*/ j, /*timIdx=*/0);
                        for (auto& tr : tbatch) {
                            if (tr.numTracer() == 0) {
                                continue;
                            }
                            this->assembleTracerEquationFlux(tr, elemCtx, scvfIdx, I, J, dt);
                        }
                    }

                     // Source terms (mass transfer between free and solution tracer)
                    for (auto& tr : tbatch) {
                        if (tr.numTracer() == 0) {
                            continue;
                        }
                        this->assembleTracerEquationSource(tr, dt, I);
                    }
                }
            }
        }
        OPM_END_PARALLEL_TRY_CATCH_LOG(local_deferredLogger,
                                       "assembleTracerEquations() failed: ",
                                       true, simulator_.gridView().comm())

        // Communicate overlap using grid Communication
        for (auto& tr : tbatch) {
            if (tr.numTracer() == 0) {
                continue;
            }
            auto handle = VectorVectorDataHandle<GridView, std::vector<TracerVector>>(tr.residual_,
                                                                                      simulator_.gridView());
            simulator_.gridView().communicate(handle, Dune::InteriorBorder_All_Interface,
                                              Dune::ForwardCommunication);
        }
    }

    template<TracerTypeIdx Index, class TrRe>
    void updateElem(TrRe& tr,
                    const Scalar scvVolume,
                    const unsigned globalDofIdx)
    {
        const Scalar vol1 = computeVolume_<Index>(tr.phaseIdx_, globalDofIdx, 0);
        vol1_[Index][tr.phaseIdx_][globalDofIdx] = vol1 * scvVolume;
        dVol_[Index][tr.phaseIdx_][globalDofIdx] = 0.0;
        for (int tIdx = 0; tIdx < tr.numTracer(); ++tIdx) {
            tr.storageOfTimeIndex1_[tIdx][globalDofIdx][Index] =
                vol1 * tr.concentrationInitial_[tIdx][globalDofIdx][Index];
        }
    }

    void updateStorageCache()
    {
        for (auto& tr : tbatch) {
            if (tr.numTracer() != 0) {
                tr.concentrationInitial_ = tr.concentration_;
            }
        }

        // Parallel loop over element chunks
        #ifdef _OPENMP
        #pragma omp parallel for
        #endif
        for (const auto& chunk : element_chunks_) {
            ElementContext elemCtx(simulator_);

            for (const auto& elem : chunk) {
                elemCtx.updatePrimaryStencil(elem);
                elemCtx.updatePrimaryIntensiveQuantities(/*timeIdx=*/0);
                const Scalar extrusionFactor = elemCtx.intensiveQuantities(/*dofIdx=*/ 0, /*timeIdx=*/0).extrusionFactor();
                const Scalar scvVolume = elemCtx.stencil(/*timeIdx=*/0).subControlVolume(/*dofIdx=*/ 0).volume() * extrusionFactor;
                const unsigned globalDofIdx = elemCtx.globalSpaceIndex(0, /*timeIdx=*/0);

                for (auto& tr : tbatch) {
                    if (tr.numTracer() == 0) {
                        continue;
                    }
                    // This is safe as each element has a unique globalDofIdx. So each thread
                    // always writes to different memory locations in the shared arrays.
                    updateElem<Free>(tr, scvVolume, globalDofIdx);
                    updateElem<Solution>(tr, scvVolume, globalDofIdx);
                }
            }
        }
    }

    template<TracerTypeIdx Index, class TrRe>
    void copyForOutput(TrRe& tr,
                       const std::vector<TracerVector>& dx,
                       const Scalar S,
                       const unsigned tIdx,
                       const unsigned globalDofIdx,
                       std::vector<TracerVectorSingle>& sc)
    {
        constexpr Scalar tol_gas_sat = 1e-6;
        tr.concentration_[tIdx][globalDofIdx][Index] -= dx[tIdx][globalDofIdx][Index];
        if (tr.concentration_[tIdx][globalDofIdx][Index] < 0.0 || S < tol_gas_sat) {
            tr.concentration_[tIdx][globalDofIdx][Index] = 0.0;
        }
        sc[tr.idx_[tIdx]][globalDofIdx] = tr.concentration_[tIdx][globalDofIdx][Index];
    }

    template<TracerTypeIdx Index, class TrRe>
    void assignRates(const TrRe& tr,
                     const Well& eclWell,
                     const std::size_t i,
                     const std::size_t I,
                     const Scalar rate,
                     std::vector<WellTracerRate<Scalar>>& tracerRate,
                     std::vector<MSWellTracerRate<Scalar>>* mswTracerRate,
                     std::vector<WellTracerRate<Scalar>>& splitRate)
    {
        if (rate < 0) {
            for (int tIdx = 0; tIdx < tr.numTracer(); ++tIdx) {
                // Store _producer_ free tracer rate for reporting
                const Scalar delta = rate * tr.concentration_[tIdx][I][Index];
                tracerRate[tIdx].rate += delta;
                splitRate[tIdx].rate += delta;
                if (eclWell.isMultiSegment()) {
                    (*mswTracerRate)[tIdx].rate[eclWell.getConnections().get(i).segment()] += delta;
                }
            }
        }
    }

    void advanceTracerFields()
    {
        assembleTracerEquations_();

        for (auto& tr : tbatch) {
            if (tr.numTracer() == 0) {
                continue;
            }

            // Note that we solve for a concentration update (compared to previous time step)
            // Confer also assembleTracerEquations_(...) above.
            std::vector<TracerVector> dx(tr.concentration_);
            for (int tIdx = 0; tIdx < tr.numTracer(); ++tIdx) {
                dx[tIdx] = 0.0;
            }

            const bool converged = this->linearSolveBatchwise_(*tr.mat, dx, tr.residual_);
            if (!converged) {
                OpmLog::warning("### Tracer model: Linear solver did not converge. ###");
            }

            OPM_TIMEBLOCK(tracerPost);

            for (int tIdx = 0; tIdx < tr.numTracer(); ++tIdx) {
                for (std::size_t globalDofIdx = 0; globalDofIdx < tr.concentration_[tIdx].size(); ++globalDofIdx) {
                    // New concetration. Concentrations that are negative or where free/solution phase is not
                    // present are set to zero
                    const auto& intQuants = simulator_.model().intensiveQuantities(globalDofIdx, 0);
                    const auto& fs = intQuants.fluidState();
                    const Scalar Sf = decay<Scalar>(fs.saturation(tr.phaseIdx_));
                    Scalar Ss = 0.0;

                    if (tr.phaseIdx_ == FluidSystem::gasPhaseIdx && FluidSystem::enableDissolvedGas()) {
                        Ss = decay<Scalar>(fs.saturation(FluidSystem::oilPhaseIdx));
                    }
                    else if (tr.phaseIdx_ == FluidSystem::oilPhaseIdx && FluidSystem::enableVaporizedOil()) {
                        Ss = decay<Scalar>(fs.saturation(FluidSystem::gasPhaseIdx));
                    }

                    copyForOutput<Free>(tr, dx, Sf, tIdx, globalDofIdx, this->freeTracerConcentration_);
                    copyForOutput<Solution>(tr, dx, Ss, tIdx, globalDofIdx, this->solTracerConcentration_);
                }
            }

            // Store _producer_ tracer rate for reporting
            const auto& wellPtrs = simulator_.problem().wellModel().localNonshutWells();
            for (const auto& wellPtr : wellPtrs) {
                const auto& eclWell = wellPtr->wellEcl();

                // Injection rates already reported during assembly
                if (!eclWell.isProducer()) {
                    continue;
                }

                Scalar rateWellPos = 0.0;
                Scalar rateWellNeg = 0.0;
                const std::size_t well_index = simulator_.problem().wellModel().wellState().index(eclWell.name()).value();
                const auto& ws = simulator_.problem().wellModel().wellState().well(well_index);
                auto& tracerRate = this->wellTracerRate_[eclWell.seqIndex()];
                auto& freeTracerRate = this->wellFreeTracerRate_[eclWell.seqIndex()];
                auto& solTracerRate = this->wellSolTracerRate_[eclWell.seqIndex()];
                auto* mswTracerRate = eclWell.isMultiSegment() ? &this->mSwTracerRate_[eclWell.seqIndex()] : nullptr;

                for (std::size_t i = 0; i < ws.perf_data.size(); ++i) {
                    const auto I = ws.perf_data.cell_index[i];
                    const Scalar rate = wellPtr->volumetricSurfaceRateForConnection(I, tr.phaseIdx_);

                    Scalar rate_s;
                    if (tr.phaseIdx_ == FluidSystem::oilPhaseIdx && FluidSystem::enableVaporizedOil()) {
                        rate_s = ws.perf_data.phase_mixing_rates[i][ws.vaporized_oil];
                    }
                    else if (tr.phaseIdx_ == FluidSystem::gasPhaseIdx && FluidSystem::enableDissolvedGas()) {
                        rate_s = ws.perf_data.phase_mixing_rates[i][ws.dissolved_gas];
                    }
                    else {
                        rate_s = 0.0;
                    }

                    const Scalar rate_f = rate - rate_s;
                    assignRates<Free>(tr, eclWell, i, I, rate_f,
                                      tracerRate, mswTracerRate, freeTracerRate);
                    assignRates<Solution>(tr, eclWell, i, I, rate_s,
                                          tracerRate, mswTracerRate, solTracerRate);

                    if (rate < 0) {
                        rateWellNeg += rate;
                    }
                    else {
                        rateWellPos += rate;
                    }
                }

                // TODO: Some inconsistencies here that perhaps should be clarified.
                // The "offical" rate as reported below is occasionally significant
                // different from the sum over connections (as calculated above). Only observed
                // for small values, neglible for the rate itself, but matters when used to
                // calculate tracer concentrations.
                const Scalar official_well_rate_total =
                    simulator_.problem().wellModel().wellState().well(well_index).surface_rates[tr.phaseIdx_];

                const Scalar rateWellTotal = official_well_rate_total;

                if (rateWellTotal > rateWellNeg) { // Cross flow
                    constexpr Scalar bucketPrDay = 10.0 / (1000. * 3600. * 24.); // ... keeps (some) trouble away
                    const Scalar factor = (rateWellTotal < -bucketPrDay) ? rateWellTotal / rateWellNeg : 0.0;
                    for (int tIdx = 0; tIdx < tr.numTracer(); ++tIdx) {
                        tracerRate[tIdx].rate *= factor;
                    }
                }
            }
        }
    }

    Simulator& simulator_;

    // This struct collects tracers of the same type (i.e, transported in same phase).
    // The idea being that, under the assumption of linearity, tracers of same type can
    // be solved in concert, having a common system matrix but separate right-hand-sides.

    // Since oil or gas tracers appears in dual compositions when VAPOIL respectively DISGAS
    // is active, the template argument is intended to support future extension to these
    // scenarios by supplying an extended vector type.

    template <typename TV>
    struct TracerBatch
    {
        std::vector<int> idx_;
        const int phaseIdx_;
        std::vector<TV> concentrationInitial_;
        std::vector<TV> concentration_;
        std::vector<TV> storageOfTimeIndex1_;
        std::vector<TV> residual_;
        std::unique_ptr<TracerMatrix> mat;

        bool operator==(const TracerBatch& rhs) const
        {
            return this->concentrationInitial_ == rhs.concentrationInitial_ &&
                   this->concentration_ == rhs.concentration_;
        }

        static TracerBatch serializationTestObject()
        {
            TracerBatch<TV> result(4);
            result.idx_ = {1,2,3};
            result.concentrationInitial_ = {5.0, 6.0};
            result.concentration_ = {7.0, 8.0};
            result.storageOfTimeIndex1_ = {9.0, 10.0, 11.0};
            result.residual_ = {12.0, 13.0};

            return result;
        }

        template<class Serializer>
        void serializeOp(Serializer& serializer)
        {
            serializer(concentrationInitial_);
            serializer(concentration_);
        }

        TracerBatch(int phaseIdx = 0) : phaseIdx_(phaseIdx) {}

        int numTracer() const
        { return idx_.size(); }

        void addTracer(const int idx, const TV& concentration)
        {
            const int numGridDof = concentration.size();
            idx_.emplace_back(idx);
            concentrationInitial_.emplace_back(concentration);
            concentration_.emplace_back(concentration);
            residual_.emplace_back(numGridDof);
            storageOfTimeIndex1_.emplace_back(numGridDof);
        }
    };

    std::array<TracerBatch<TracerVector>,numPhases> tbatch;
    TracerBatch<TracerVector>& wat_;
    TracerBatch<TracerVector>& oil_;
    TracerBatch<TracerVector>& gas_;
    std::array<std::array<std::vector<Scalar>,numPhases>,2> vol1_;
    std::array<std::array<std::vector<Scalar>,numPhases>,2> dVol_;
    ElementChunks<GridView, Dune::Partitions::All> element_chunks_;
};

} // namespace Opm

#endif // OPM_TRACER_MODEL_HPP
