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

#include <opm/models/utils/propertysystem.hh>

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
 *
 * TODO: MPI parallelism.
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
                   simulator.model().dofMapper())
        , simulator_(simulator)
    { }


    /*!
     * \brief Initialize all internal data structures needed by the tracer module
     */
    void init()
    {
        bool enabled = EWOMS_GET_PARAM(TypeTag, bool, EnableTracerModel);
        this->doInit(enabled, simulator_.model().numGridDof(),
                     gasPhaseIdx, oilPhaseIdx, waterPhaseIdx);
    }

    void beginTimeStep()
    {
        if (this->numTracers()==0)
            return;

        this->tracerConcentrationInitial_ = this->tracerConcentration_;

        // compute storageCache
        ElementContext elemCtx(simulator_);
        auto elemIt = simulator_.gridView().template begin</*codim=*/0>();
        auto elemEndIt = simulator_.gridView().template end</*codim=*/0>();
        for (; elemIt != elemEndIt; ++ elemIt) {
            elemCtx.updateAll(*elemIt);
            int globalDofIdx = elemCtx.globalSpaceIndex(0, 0);
            for (int tracerIdx = 0; tracerIdx < this->numTracers(); ++ tracerIdx){
                Scalar storageOfTimeIndex1;
                computeStorage_(storageOfTimeIndex1, elemCtx, 0, /*timIdx=*/0, tracerIdx);
                this->storageOfTimeIndex1_[tracerIdx][globalDofIdx] = storageOfTimeIndex1;
            }
        }
    }

    /*!
     * \brief Informs the tracer model that a time step has just been finished.
     */
    void endTimeStep()
    {
        if (this->numTracers()==0)
            return;

        for (int tracerIdx = 0; tracerIdx < this->numTracers(); ++ tracerIdx){

            typename BaseType::TracerVector dx(this->tracerResidual_.size());
            // Newton step (currently the system is linear, converge in one iteration)
            for (int iter = 0; iter < 5; ++ iter){
                linearize_(tracerIdx);
                this->linearSolve_(*this->tracerMatrix_, dx, this->tracerResidual_);
                this->tracerConcentration_[tracerIdx] -= dx;

                if (dx.two_norm()<1e-2)
                    break;
            }
        }
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

protected:
    // evaluate storage term for all tracers in a single cell
    template <class LhsEval>
    void computeStorage_(LhsEval& tracerStorage,
                         const ElementContext& elemCtx,
                         unsigned scvIdx,
                         unsigned timeIdx,
                         const int tracerIdx)
    {
        int globalDofIdx = elemCtx.globalSpaceIndex(scvIdx, timeIdx);

        const auto& intQuants = elemCtx.intensiveQuantities(scvIdx, timeIdx);
        const auto& fs = intQuants.fluidState();
        Scalar phaseVolume =
            decay<Scalar>(fs.saturation(this->tracerPhaseIdx_[tracerIdx]))
            *decay<Scalar>(fs.invB(this->tracerPhaseIdx_[tracerIdx]))
            *decay<Scalar>(intQuants.porosity());

        // avoid singular matrix if no water is present.
        phaseVolume = max(phaseVolume, 1e-10);

        if (std::is_same<LhsEval, Scalar>::value)
            tracerStorage = phaseVolume * this->tracerConcentrationInitial_[tracerIdx][globalDofIdx];
        else
            tracerStorage =
                    phaseVolume
                    * variable<LhsEval>(this->tracerConcentration_[tracerIdx][globalDofIdx][0], 0);
    }

    // evaluate the tracerflux over one face
    void computeFlux_(TracerEvaluation & tracerFlux,
                      const ElementContext& elemCtx,
                      unsigned scvfIdx,
                      unsigned timeIdx,
                      const int tracerIdx)

    {
        const auto& stencil = elemCtx.stencil(timeIdx);
        const auto& scvf = stencil.interiorFace(scvfIdx);

        const auto& extQuants = elemCtx.extensiveQuantities(scvfIdx, timeIdx);
        unsigned inIdx = extQuants.interiorIndex();

        const int tracerPhaseIdx = this->tracerPhaseIdx_[tracerIdx];

        unsigned upIdx = extQuants.upstreamIndex(tracerPhaseIdx);
        int globalUpIdx = elemCtx.globalSpaceIndex(upIdx, timeIdx);

        const auto& intQuants = elemCtx.intensiveQuantities(upIdx, timeIdx);
        const auto& fs = intQuants.fluidState();

        Scalar A = scvf.area();
        Scalar v = decay<Scalar>(extQuants.volumeFlux(tracerPhaseIdx));
        Scalar b = decay<Scalar>(fs.invB(this->tracerPhaseIdx_[tracerIdx]));
        Scalar c = this->tracerConcentration_[tracerIdx][globalUpIdx];

        if (inIdx == upIdx)
            tracerFlux = A*v*b*variable<TracerEvaluation>(c, 0);
        else
            tracerFlux = A*v*b*c;

    }

    void linearize_(int tracerIdx)
    {
        (*this->tracerMatrix_) = 0.0;
        this->tracerResidual_ = 0.0;

        size_t numGridDof =  simulator_.model().numGridDof();
        std::vector<double> volumes(numGridDof, 0.0);
        ElementContext elemCtx(simulator_);
        auto elemIt = simulator_.gridView().template begin</*codim=*/0>();
        auto elemEndIt = simulator_.gridView().template end</*codim=*/0>();
        for (; elemIt != elemEndIt; ++ elemIt) {
            elemCtx.updateAll(*elemIt);

            Scalar extrusionFactor =
                    elemCtx.intensiveQuantities(/*dofIdx=*/ 0, /*timeIdx=*/0).extrusionFactor();
            Valgrind::CheckDefined(extrusionFactor);
            assert(isfinite(extrusionFactor));
            assert(extrusionFactor > 0.0);
            Scalar scvVolume =
                    elemCtx.stencil(/*timeIdx=*/0).subControlVolume(/*dofIdx=*/ 0).volume()
                    * extrusionFactor;
            Scalar dt = elemCtx.simulator().timeStepSize();

            size_t I = elemCtx.globalSpaceIndex(/*dofIdx=*/ 0, /*timIdx=*/0);
            volumes[I] = scvVolume;
            TracerEvaluation localStorage;
            TracerEvaluation storageOfTimeIndex0;
            Scalar storageOfTimeIndex1;
            computeStorage_(storageOfTimeIndex0, elemCtx, 0, /*timIdx=*/0, tracerIdx);
            if (elemCtx.enableStorageCache())
                storageOfTimeIndex1 = this->storageOfTimeIndex1_[tracerIdx][I];
            else
                computeStorage_(storageOfTimeIndex1, elemCtx, 0, /*timIdx=*/1, tracerIdx);

            localStorage = (storageOfTimeIndex0 - storageOfTimeIndex1) * scvVolume/dt;
            this->tracerResidual_[I][0] += localStorage.value(); //residual + flux
            (*this->tracerMatrix_)[I][I][0][0] = localStorage.derivative(0);
            size_t numInteriorFaces = elemCtx.numInteriorFaces(/*timIdx=*/0);
            for (unsigned scvfIdx = 0; scvfIdx < numInteriorFaces; scvfIdx++) {
                TracerEvaluation flux;
                const auto& face = elemCtx.stencil(0).interiorFace(scvfIdx);
                unsigned j = face.exteriorIndex();
                unsigned J = elemCtx.globalSpaceIndex(/*dofIdx=*/ j, /*timIdx=*/0);
                computeFlux_(flux, elemCtx, scvfIdx, 0, tracerIdx);
                this->tracerResidual_[I][0] += flux.value(); //residual + flux
                (*this->tracerMatrix_)[J][I][0][0] = -flux.derivative(0);
                (*this->tracerMatrix_)[I][J][0][0] = flux.derivative(0);
            }

        }

        // Wells
        const int episodeIdx = simulator_.episodeIndex();
        const auto& wells = simulator_.vanguard().schedule().getWells(episodeIdx);
        for (const auto& well : wells) {

            if (well.getStatus() == Well::Status::SHUT)
                continue;

            const double wtracer = well.getTracerProperties().getConcentration(this->tracerNames_[tracerIdx]);
            std::array<int, 3> cartesianCoordinate;
            for (auto& connection : well.getConnections()) {

                if (connection.state() == Connection::State::SHUT)
                    continue;

                cartesianCoordinate[0] = connection.getI();
                cartesianCoordinate[1] = connection.getJ();
                cartesianCoordinate[2] = connection.getK();
                const size_t cartIdx = simulator_.vanguard().cartesianIndex(cartesianCoordinate);
                const int I = this->cartToGlobal_[cartIdx];
                Scalar rate = simulator_.problem().wellModel().well(well.name())->volumetricSurfaceRateForConnection(I, this->tracerPhaseIdx_[tracerIdx]);
                if (rate > 0)
                    this->tracerResidual_[I][0] -= rate*wtracer;
                else if (rate < 0)
                    this->tracerResidual_[I][0] -= rate*this->tracerConcentration_[tracerIdx][I];
            }
        }
    }

    Simulator& simulator_;
};

} // namespace Opm

#endif
