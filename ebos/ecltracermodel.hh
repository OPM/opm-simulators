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
 * \copydoc Opm::EclTracerModel
 */
#ifndef EWOMS_ECL_TRACER_MODEL_HH
#define EWOMS_ECL_TRACER_MODEL_HH

#include <ebos/eclgenerictracermodel.hh>

#include <opm/common/OpmLog/OpmLog.hpp>

#include <opm/models/utils/propertysystem.hh>

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
 * \ingroup EclBlackOilSimulator
 *
 * \brief A class which handles tracers as specified in by ECL
 */
template <class TypeTag>
class EclTracerModel : public EclGenericTracerModel<GetPropType<TypeTag, Properties::Grid>,
                                                    GetPropType<TypeTag, Properties::GridView>,
                                                    GetPropType<TypeTag, Properties::DofMapper>,
                                                    GetPropType<TypeTag, Properties::Stencil>,
                                                    GetPropType<TypeTag, Properties::Scalar>>
{
    using BaseType = EclGenericTracerModel<GetPropType<TypeTag, Properties::Grid>,
                                           GetPropType<TypeTag, Properties::GridView>,
                                           GetPropType<TypeTag, Properties::DofMapper>,
                                           GetPropType<TypeTag, Properties::Stencil>,
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

    enum { numEq = getPropValue<TypeTag, Properties::NumEq>() };
    enum { numPhases = FluidSystem::numPhases };
    enum { waterPhaseIdx = FluidSystem::waterPhaseIdx };
    enum { oilPhaseIdx = FluidSystem::oilPhaseIdx };
    enum { gasPhaseIdx = FluidSystem::gasPhaseIdx };

public:
    EclTracerModel(Simulator& simulator)
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
    { }


    /*
      The initialization of the tracer model is a three step process:

      1. The init() method is called. This will allocate buffers and initialize
         some phase index stuff. If this is a normal run the initial tracer
         concentrations will be assigned from the TBLK or TVDPF keywords.

      2. [Restart only:] The tracer concenntration are read from the restart
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
                    throw std::runtime_error("Water tracer specified for non-water fluid system:" + this->name(tracerIdx));
                }

                wat_.addTracer(tracerIdx, this->tracerConcentration_[tracerIdx]);
            }
            else if (this->tracerPhaseIdx_[tracerIdx] == FluidSystem::oilPhaseIdx) {
                if (! FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx)){
                    throw std::runtime_error("Oil tracer specified for non-oil fluid system:" + this->name(tracerIdx));
                }
                if (FluidSystem::enableVaporizedOil()) {
                    throw std::runtime_error("Oil tracer in combination with kw VAPOIL is not supported: " + this->name(tracerIdx));
                }

                oil_.addTracer(tracerIdx, this->tracerConcentration_[tracerIdx]);
            }
            else if (this->tracerPhaseIdx_[tracerIdx] == FluidSystem::gasPhaseIdx) {
                if (! FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)){
                    throw std::runtime_error("Gas tracer specified for non-gas fluid system:" + this->name(tracerIdx));
                }
                if (FluidSystem::enableDissolvedGas()) {
                    throw std::runtime_error("Gas tracer in combination with kw DISGAS is not supported: " + this->name(tracerIdx));
                }

                gas_.addTracer(tracerIdx, this->tracerConcentration_[tracerIdx]);
            }
        }

        // will be valid after we move out of tracerMatrix_
        TracerMatrix* base = this->tracerMatrix_.get();
        for (auto& tr : this->tbatch) {
            if (tr.numTracer() != 0) {
                if (this->tracerMatrix_)
                    tr.mat = std::move(this->tracerMatrix_);
                else
                    tr.mat = std::make_unique<TracerMatrix>(*base);
            }
        }
    }

    void beginTimeStep()
    {
        if (this->numTracers() == 0)
            return;

        updateStorageCache();
    }

    /*!
     * \brief Informs the tracer model that a time step has just been finished.
     */
    void endTimeStep()
    {
        if (this->numTracers() == 0)
            return;

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

    // evaluate water storage volume(s) in a single cell
    template <class LhsEval>
    void computeVolume_(LhsEval& freeVolume,
                        const int tracerPhaseIdx,
                        const ElementContext& elemCtx,
                        unsigned scvIdx,
                        unsigned timeIdx)
    {
        const auto& intQuants = elemCtx.intensiveQuantities(scvIdx, timeIdx);
        const auto& fs = intQuants.fluidState();
        Scalar phaseVolume =
            decay<Scalar>(fs.saturation(tracerPhaseIdx))
            *decay<Scalar>(fs.invB(tracerPhaseIdx))
            *decay<Scalar>(intQuants.porosity());

        // avoid singular matrix if no water is present.
        phaseVolume = max(phaseVolume, 1e-10);

        if (std::is_same<LhsEval, Scalar>::value)
            freeVolume = phaseVolume;
        else
            freeVolume = phaseVolume * variable<LhsEval>(1.0, 0);

    }

    // evaluate the flux(es) over one face
    void computeFlux_(TracerEvaluation & freeFlux,
                      bool & isUpFree,
                      const int tracerPhaseIdx,
                      const ElementContext& elemCtx,
                      unsigned scvfIdx,
                      unsigned timeIdx)

    {
        const auto& stencil = elemCtx.stencil(timeIdx);
        const auto& scvf = stencil.interiorFace(scvfIdx);

        const auto& extQuants = elemCtx.extensiveQuantities(scvfIdx, timeIdx);
        unsigned inIdx = extQuants.interiorIndex();

        unsigned upIdx = extQuants.upstreamIndex(tracerPhaseIdx);

        const auto& intQuants = elemCtx.intensiveQuantities(upIdx, timeIdx);
        const auto& fs = intQuants.fluidState();

        Scalar A = scvf.area();
        Scalar v = decay<Scalar>(extQuants.volumeFlux(tracerPhaseIdx));
        Scalar b = decay<Scalar>(fs.invB(tracerPhaseIdx));

        if (inIdx == upIdx) {
            freeFlux = A*v*b*variable<TracerEvaluation>(1.0, 0);
            isUpFree = true;
        }
        else {
            freeFlux = A*v*b*1.0;
            isUpFree = false;
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
        if (tr.numTracer() == 0)
            return;

        std::vector<Scalar> storageOfTimeIndex1(tr.numTracer());
        if (elemCtx.enableStorageCache()) {
            for (int tIdx = 0; tIdx < tr.numTracer(); ++tIdx) {
                storageOfTimeIndex1[tIdx] = tr.storageOfTimeIndex1_[tIdx][I];
            }
        }
        else {
            Scalar fVolume1;
            computeVolume_(fVolume1, tr.phaseIdx_, elemCtx, 0, /*timeIdx=*/1);
            for (int tIdx = 0; tIdx < tr.numTracer(); ++tIdx) {
                storageOfTimeIndex1[tIdx] = fVolume1*tr.concentrationInitial_[tIdx][I1];
            }
        }

        TracerEvaluation fVolume;
        computeVolume_(fVolume, tr.phaseIdx_, elemCtx, 0, /*timeIdx=*/0);
        for (int tIdx = 0; tIdx < tr.numTracer(); ++tIdx) {
            Scalar storageOfTimeIndex0 = fVolume.value()*tr.concentration_[tIdx][I];
            Scalar localStorage = (storageOfTimeIndex0 - storageOfTimeIndex1[tIdx]) * scvVolume/dt;
            tr.residual_[tIdx][I][0] += localStorage; //residual + flux
        }
        (*tr.mat)[I][I][0][0] += fVolume.derivative(0) * scvVolume/dt;
    }

    template<class TrRe>
    void assembleTracerEquationFlux(TrRe& tr,
                                    const ElementContext& elemCtx,
                                    unsigned scvfIdx,
                                    unsigned I,
                                    unsigned J)
    {
        if (tr.numTracer() == 0)
            return;

        TracerEvaluation flux;
        bool isUpF;
        computeFlux_(flux, isUpF, tr.phaseIdx_, elemCtx, scvfIdx, 0);
        int globalUpIdx = isUpF ? I : J;
        for (int tIdx =0; tIdx < tr.numTracer(); ++tIdx) {
            tr.residual_[tIdx][I][0] += flux.value()*tr.concentration_[tIdx][globalUpIdx]; //residual + flux
        }
        if (isUpF) {
            (*tr.mat)[J][I][0][0] = -flux.derivative(0);
            (*tr.mat)[I][I][0][0] += flux.derivative(0);
        }
    }

    template<class TrRe, class Well>
    void assembleTracerEquationWell(TrRe& tr,
                                    const Well& well)
    {
        if (tr.numTracer() == 0)
            return;

        const auto& eclWell = well.wellEcl();
        for (int tIdx = 0; tIdx < tr.numTracer(); ++tIdx) {
            this->wellTracerRate_[std::make_pair(eclWell.name(), this->name(tr.idx_[tIdx]))] = 0.0;
        }

        std::vector<double> wtracer(tr.numTracer());
        for (int tIdx = 0; tIdx < tr.numTracer(); ++tIdx) {
            wtracer[tIdx] = this->currentConcentration_(eclWell, this->name(tr.idx_[tIdx]));
        }

        for (auto& perfData : well.perforationData()) {
            const auto I = perfData.cell_index;
            Scalar rate = well.volumetricSurfaceRateForConnection(I, tr.phaseIdx_);
            if (rate > 0) {
                for (int tIdx = 0; tIdx < tr.numTracer(); ++tIdx) {
                    tr.residual_[tIdx][I][0] -= rate*wtracer[tIdx];
                    // Store _injector_ tracer rate for reporting
                    this->wellTracerRate_.at(std::make_pair(eclWell.name(),this->name(tr.idx_[tIdx]))) += rate*wtracer[tIdx];
                }
            }
            else if (rate < 0) {
                for (int tIdx = 0; tIdx < tr.numTracer(); ++tIdx) {
                    tr.residual_[tIdx][I][0] -= rate*tr.concentration_[tIdx][I];
                }
                (*tr.mat)[I][I][0][0] -= rate*variable<TracerEvaluation>(1.0, 0).derivative(0);
            }
        }
    }

    void assembleTracerEquations_()
    {
        // Note that we formulate the equations in terms of a concentration update
        // (compared to previous time step) and not absolute concentration.
        // This implies that current concentration (tr.concentration_[][]) contributes
        // to the rhs both through storrage and flux terms.
        // Compare also advanceTracerFields(...) below.

        for (auto& tr : tbatch) {
            if (tr.numTracer() != 0) {
                (*tr.mat) = 0.0;
                for (int tIdx = 0; tIdx < tr.numTracer(); ++tIdx)
                    tr.residual_[tIdx] = 0.0;
            }
        }

        ElementContext elemCtx(simulator_);
        for (const auto& elem : elements(simulator_.gridView())) {
            elemCtx.updateStencil(elem);

            std::size_t I = elemCtx.globalSpaceIndex(/*dofIdx=*/ 0, /*timeIdx=*/0);

            if (elem.partitionType() != Dune::InteriorEntity)
            {
                // Dirichlet boundary conditions needed for the parallel matrix
                for (auto& tr : tbatch) {
                    if (tr.numTracer() != 0)
                        (*tr.mat)[I][I][0][0] = 1.;
                }
                continue;
            }
            elemCtx.updateAllIntensiveQuantities();
            elemCtx.updateAllExtensiveQuantities();

            Scalar extrusionFactor =
                    elemCtx.intensiveQuantities(/*dofIdx=*/ 0, /*timeIdx=*/0).extrusionFactor();
            Valgrind::CheckDefined(extrusionFactor);
            assert(isfinite(extrusionFactor));
            assert(extrusionFactor > 0.0);
            Scalar scvVolume =
                    elemCtx.stencil(/*timeIdx=*/0).subControlVolume(/*dofIdx=*/ 0).volume()
                    * extrusionFactor;
            Scalar dt = elemCtx.simulator().timeStepSize();

            std::size_t I1 = elemCtx.globalSpaceIndex(/*dofIdx=*/ 0, /*timeIdx=*/1);

            for (auto& tr : tbatch) {
                this->assembleTracerEquationVolume(tr, elemCtx, scvVolume, dt, I, I1);
            }

            std::size_t numInteriorFaces = elemCtx.numInteriorFaces(/*timIdx=*/0);
            for (unsigned scvfIdx = 0; scvfIdx < numInteriorFaces; scvfIdx++) {
                const auto& face = elemCtx.stencil(0).interiorFace(scvfIdx);
                unsigned j = face.exteriorIndex();
                unsigned J = elemCtx.globalSpaceIndex(/*dofIdx=*/ j, /*timIdx=*/0);
                for (auto& tr : tbatch) {
                    this->assembleTracerEquationFlux(tr, elemCtx, scvfIdx, I, J);
                }
            }
        }

        // Well terms
        const auto& wellPtrs = simulator_.problem().wellModel().localNonshutWells();
        for (const auto& wellPtr : wellPtrs) {
            for (auto& tr : tbatch) {
                this->assembleTracerEquationWell(tr, *wellPtr);
            }
        }

        // Communicate overlap using grid Communication
        for (auto& tr : tbatch) {
            if (tr.numTracer() == 0)
                continue;
            auto handle = VectorVectorDataHandle<GridView, std::vector<TracerVector>>(tr.residual_,
                                                                                      simulator_.gridView());
            simulator_.gridView().communicate(handle, Dune::InteriorBorder_All_Interface,
                                              Dune::ForwardCommunication);
        }
    }

    void updateStorageCache()
    {
        for (auto& tr : tbatch) {
            if (tr.numTracer() != 0)
                tr.concentrationInitial_ = tr.concentration_;
        }

        ElementContext elemCtx(simulator_);
        for (const auto& elem : elements(simulator_.gridView())) {
            elemCtx.updatePrimaryStencil(elem);
            elemCtx.updatePrimaryIntensiveQuantities(/*timeIdx=*/0);
            int globalDofIdx = elemCtx.globalSpaceIndex(0, /*timeIdx=*/0);
            for (auto& tr : tbatch) {
                if (tr.numTracer() == 0)
                    continue;
                Scalar fVolume;
                computeVolume_(fVolume, tr.phaseIdx_, elemCtx, 0, /*timeIdx=*/0);
                for (int tIdx = 0; tIdx < tr.numTracer(); ++tIdx) {
                    tr.storageOfTimeIndex1_[tIdx][globalDofIdx] = fVolume*tr.concentrationInitial_[tIdx][globalDofIdx];
                }
            }
        }
    }

    void advanceTracerFields()
    {
        assembleTracerEquations_();

        for (auto& tr : tbatch) {
            if (tr.numTracer() == 0)
                continue;

            // Note that we solve for a concentration update (compared to previous time step)
            // Confer also assembleTracerEquations_(...) above.
            std::vector<TracerVector> dx(tr.concentration_);
            for (int tIdx = 0; tIdx < tr.numTracer(); ++tIdx)
                dx[tIdx] = 0.0;

            bool converged = this->linearSolveBatchwise_(*tr.mat, dx, tr.residual_);
            if (!converged) {
                OpmLog::warning("### Tracer model: Linear solver did not converge. ###");
            }

            for (int tIdx = 0; tIdx < tr.numTracer(); ++tIdx) {
                tr.concentration_[tIdx] -= dx[tIdx];
                // Tracer concentrations for restart report
                this->tracerConcentration_[tr.idx_[tIdx]] = tr.concentration_[tIdx];
            }

            // Store _producer_ tracer rate for reporting
            const auto& wellPtrs = simulator_.problem().wellModel().localNonshutWells();
            for (const auto& wellPtr : wellPtrs) {
                const auto& well = wellPtr->wellEcl();

                if (!well.isProducer()) //Injection rates already reported during assembly
                    continue;

                Scalar rateWellPos = 0.0;
                Scalar rateWellNeg = 0.0;
                for (auto& perfData : wellPtr->perforationData()) {
                    const int I = perfData.cell_index;
                    Scalar rate = wellPtr->volumetricSurfaceRateForConnection(I, tr.phaseIdx_);
                    if (rate < 0) {
                        rateWellNeg += rate;
                        for (int tIdx = 0; tIdx < tr.numTracer(); ++tIdx) {
                            this->wellTracerRate_.at(std::make_pair(well.name(),this->name(tr.idx_[tIdx]))) += rate*tr.concentration_[tIdx][I];
                        }
                    }
                    else {
                        rateWellPos += rate;
                    }
                }

                Scalar rateWellTotal = rateWellNeg + rateWellPos;

                // TODO: Some inconsistencies here that perhaps should be clarified. The "offical" rate as reported below is
                //  occasionally significant different from the sum over connections (as calculated above). Only observed
                //  for small values, neglible for the rate itself, but matters when used to calculate tracer concentrations.
                std::size_t well_index = simulator_.problem().wellModel().wellState().index(well.name()).value();
                Scalar official_well_rate_total = simulator_.problem().wellModel().wellState().well(well_index).surface_rates[tr.phaseIdx_];

                rateWellTotal = official_well_rate_total;

                if (rateWellTotal > rateWellNeg) { // Cross flow
                    const Scalar bucketPrDay = 10.0/(1000.*3600.*24.); // ... keeps (some) trouble away
                    const Scalar factor = (rateWellTotal < -bucketPrDay) ? rateWellTotal/rateWellNeg : 0.0;
                    for (int tIdx = 0; tIdx < tr.numTracer(); ++tIdx) {
                        this->wellTracerRate_.at(std::make_pair(well.name(),this->name(tr.idx_[tIdx]))) *= factor;
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
    struct TracerBatch {
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

      int numTracer() const {return idx_.size(); }

      void addTracer(const int idx, const TV & concentration)
      {
          int numGridDof = concentration.size();
          idx_.emplace_back(idx);
          concentrationInitial_.emplace_back(concentration);
          concentration_.emplace_back(concentration);
          storageOfTimeIndex1_.emplace_back(numGridDof);
          residual_.emplace_back(numGridDof);
      }
    };

    std::array<TracerBatch<TracerVector>,3> tbatch;
    TracerBatch<TracerVector>& wat_;
    TracerBatch<TracerVector>& oil_;
    TracerBatch<TracerVector>& gas_;
};

} // namespace Opm

#endif
