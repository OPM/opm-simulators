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
 * \copydoc Ewoms::EclTracerModel
 */
#ifndef EWOMS_ECL_TRACER_MODEL_HH
#define EWOMS_ECL_TRACER_MODEL_HH

#include "tracervdtable.hh"

#include <ewoms/models/blackoil/blackoilmodel.hh>

#include <dune/istl/operators.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/preconditioners.hh>

#include <dune/common/version.hh>

#include <string>
#include <vector>
#include <iostream>

BEGIN_PROPERTIES

NEW_PROP_TAG(EnableTracerModel);

END_PROPERTIES

namespace Ewoms {

/*!
 * \ingroup EclBlackOilSimulator
 *
 * \brief A class which handles tracers as specified in by ECL
 *
 * TODO: MPI parallelism.
 */
template <class TypeTag>
class EclTracerModel
{
    typedef typename GET_PROP_TYPE(TypeTag, Simulator) Simulator;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Grid) Grid;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Stencil) Stencil;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP_TYPE(TypeTag, RateVector) RateVector;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

    typedef Opm::DenseAd::Evaluation<Scalar,1> TracerEvaluation;

    enum { numEq = GET_PROP_VALUE(TypeTag, NumEq) };
    enum { numPhases = FluidSystem::numPhases };
    enum { waterPhaseIdx = FluidSystem::waterPhaseIdx };
    enum { oilPhaseIdx = FluidSystem::oilPhaseIdx };
    enum { gasPhaseIdx = FluidSystem::gasPhaseIdx };

    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;

    typedef Dune::BCRSMatrix<Dune::FieldMatrix<Scalar, 1, 1>> TracerMatrix;
    typedef Dune::BlockVector<Dune::FieldVector<Scalar,1>> TracerVector;

public:
    EclTracerModel(Simulator& simulator)
        : simulator_(simulator)
    { }


    /*!
     * \brief Initialize all internal data structures needed by the tracer module
     */
    void init()
    {
        const Opm::Deck& deck = simulator_.vanguard().deck();

        if (!deck.hasKeyword("TRACERS"))
            return; // tracer treatment is supposed to be disabled

        if (!EWOMS_GET_PARAM(TypeTag, bool, EnableTracerModel)) {
            if (simulator_.gridView().comm().rank() == 0) {
                std::cout << "Warning: Tracer model is disabled but the deck contains the TRACERS keyword\n"
                          << "The tracer model must be explictly activated using --enable-tracer-model=true\n"
                          << std::flush;
            }
            return; // Tracer transport must be enabled by the user
        }

        if (!deck.hasKeyword("TRACER"))
            throw std::runtime_error("The deck does not contain the TRACER keyword");

        if (simulator_.gridView().comm().size() > 1) {
            tracerNames_.resize(0);
            if (simulator_.gridView().comm().rank() == 0)
                std::cout << "Warning: The tracer model currently does not work for parallel runs\n"
                          << std::flush;
            return;
        }

        // retrieve the number of tracers from the deck
        const int numTracers = deck.getKeyword("TRACER").size();
        tracerNames_.resize(numTracers);
        tracerConcentration_.resize(numTracers);
        storageOfTimeIndex1_.resize(numTracers);

        // the phase where the tracer is
        tracerPhaseIdx_.resize(numTracers);
        size_t numGridDof =  simulator_.model().numGridDof();
        for (int tracerIdx = 0;  tracerIdx < numTracers; ++tracerIdx) {
            const auto& tracerRecord = deck.getKeyword("TRACER").getRecord(tracerIdx);
            tracerNames_[tracerIdx] = tracerRecord.getItem("NAME").template get<std::string>(0);
            const std::string& fluidName = tracerRecord.getItem("FLUID").template get<std::string>(0);
            if (fluidName == "WAT")
                tracerPhaseIdx_[tracerIdx] = waterPhaseIdx;
            else if (fluidName == "OIL")
                tracerPhaseIdx_[tracerIdx] = oilPhaseIdx;
            else if (fluidName == "GAS")
                tracerPhaseIdx_[tracerIdx] = gasPhaseIdx;
            else
                throw std::invalid_argument("Tracer: invalid fluid name "
                                            +fluidName+" for "+tracerNames_[tracerIdx]);

            tracerConcentration_[tracerIdx].resize(numGridDof);
            storageOfTimeIndex1_[tracerIdx].resize(numGridDof);
            std::string tmp = "TVDPF" +tracerNames_[tracerIdx];


            //TBLK keyword
            if (deck.hasKeyword("TBLKF" +tracerNames_[tracerIdx])){
                const auto& cartMapper = simulator_.vanguard().cartesianIndexMapper();
                const auto& tblkData =
                        deck.getKeyword("TBLKF"
                                        +tracerNames_
                                        [tracerIdx]).getRecord(0).getItem(0).getSIDoubleData();
                int tblkDatasize = tblkData.size();
                if (tblkDatasize < simulator_.vanguard().cartesianSize()){
                    throw std::runtime_error("Uninitialized tracer concentration (TBLKF) for tracer "
                                             + tracerName(tracerIdx));
                }
                for (size_t globalDofIdx = 0; globalDofIdx < numGridDof; ++globalDofIdx){
                    int cartDofIdx = cartMapper.cartesianIndex(globalDofIdx);
                    tracerConcentration_[tracerIdx][globalDofIdx] = tblkData[cartDofIdx];
                }
            }
            //TVDPF keyword
            else if (deck.hasKeyword(tmp)){
                TracerVdTable dtable(deck.getKeyword(tmp).getRecord(0).getItem(0));
                const auto& eclGrid = simulator_.vanguard().eclState().getInputGrid();
                const auto& cartMapper = simulator_.vanguard().cartesianIndexMapper();

                for (size_t globalDofIdx = 0; globalDofIdx < numGridDof; ++globalDofIdx){
                    int cartDofIdx = cartMapper.cartesianIndex(globalDofIdx);
                    const auto& center = eclGrid.getCellCenter(cartDofIdx);
                    tracerConcentration_[tracerIdx][globalDofIdx]
                            = dtable.evaluate("TRACER_CONCENTRATION", center[2]);
                }
            }
            else {
                throw std::runtime_error("Uninitialized tracer concentration for tracer "
                                         + tracerName(tracerIdx));
            }

        }

        // initial tracer concentration
        tracerConcentrationInitial_ = tracerConcentration_;

        // residual of tracers
        tracerResidual_.resize(numGridDof);

        // allocate matrix for storing the Jacobian of the tracer residual
        tracerMatrix_ = new TracerMatrix(numGridDof, numGridDof, TracerMatrix::random);

        // find the sparsity pattern of the tracer matrix
        typedef std::set<unsigned> NeighborSet;
        std::vector<NeighborSet> neighbors(numGridDof);

        Stencil stencil(simulator_.gridView(), simulator_.model().dofMapper() );
        ElementIterator elemIt = simulator_.gridView().template begin<0>();
        const ElementIterator elemEndIt = simulator_.gridView().template end<0>();
        for (; elemIt != elemEndIt; ++elemIt) {
            const Element& elem = *elemIt;
            stencil.update(elem);

            for (unsigned primaryDofIdx = 0; primaryDofIdx < stencil.numPrimaryDof(); ++primaryDofIdx) {
                unsigned myIdx = stencil.globalSpaceIndex(primaryDofIdx);

                for (unsigned dofIdx = 0; dofIdx < stencil.numDof(); ++dofIdx) {
                    unsigned neighborIdx = stencil.globalSpaceIndex(dofIdx);
                    neighbors[myIdx].insert(neighborIdx);
                }
            }
        }

        // allocate space for the rows of the matrix
        for (unsigned dofIdx = 0; dofIdx < numGridDof; ++ dofIdx)
            tracerMatrix_->setrowsize(dofIdx, neighbors[dofIdx].size());
        tracerMatrix_->endrowsizes();

        // fill the rows with indices. each degree of freedom talks to
        // all of its neighbors. (it also talks to itself since
        // degrees of freedom are sometimes quite egocentric.)
        for (unsigned dofIdx = 0; dofIdx < numGridDof; ++ dofIdx) {
            typename NeighborSet::iterator nIt = neighbors[dofIdx].begin();
            typename NeighborSet::iterator nEndIt = neighbors[dofIdx].end();
            for (; nIt != nEndIt; ++nIt)
                tracerMatrix_->addindex(dofIdx, *nIt);
        }
        tracerMatrix_->endindices();

        const int sizeCartGrid = simulator_.vanguard().cartesianSize();
        cartToGlobal_.resize(sizeCartGrid);
        for (unsigned i = 0; i < numGridDof; ++i) {
            int cartIdx = simulator_.vanguard().cartesianIndex(i);
            cartToGlobal_[cartIdx] = i;
        }

    }

    /*!
     * \brief Return the number of tracers considered by the tracerModel.
     */
    int numTracers() const
    { return tracerNames_.size(); }

    /*!
     * \brief Return the tracer name
     */
    const std::string& tracerName(int tracerIdx) const
    {
        if (numTracers()==0)
            throw std::logic_error("This method should never be called when there are no tracers in the model");

        return tracerNames_[tracerIdx];
    }

    /*!
     * \brief Return the tracer concentration for tracer index and global DofIdx
     */
    Scalar tracerConcentration(int tracerIdx, int globalDofIdx) const
    {
        if (numTracers()==0)
            return 0.0;

        return tracerConcentration_[tracerIdx][globalDofIdx];
    }

    void beginTimeStep()
    {
        if (numTracers()==0)
            return;

        tracerConcentrationInitial_ = tracerConcentration_;

        // compute storageCache
        ElementContext elemCtx(simulator_);
        auto elemIt = simulator_.gridView().template begin</*codim=*/0>();
        auto elemEndIt = simulator_.gridView().template end</*codim=*/0>();
        for (; elemIt != elemEndIt; ++ elemIt) {
            elemCtx.updateAll(*elemIt);
            int globalDofIdx = elemCtx.globalSpaceIndex(0, 0);
            for (int tracerIdx = 0; tracerIdx < numTracers(); ++ tracerIdx){
                Scalar storageOfTimeIndex1;
                computeStorage_(storageOfTimeIndex1, elemCtx, 0, /*timIdx=*/0, tracerIdx);
                storageOfTimeIndex1_[tracerIdx][globalDofIdx] = storageOfTimeIndex1;
            }
        }
    }

    /*!
     * \brief Informs the tracer model that a time step has just been finished.
     */
    void endTimeStep()
    {
        if (numTracers()==0)
            return;

        for (int tracerIdx = 0; tracerIdx < numTracers(); ++ tracerIdx){

            TracerVector dx(tracerResidual_.size());
            // Newton step (currently the system is linear, converge in one iteration)
            for (int iter = 0; iter < 5; ++ iter){
                linearize_(tracerIdx);
                linearSolve_(*tracerMatrix_, dx, tracerResidual_);
                tracerConcentration_[tracerIdx] -= dx;

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
    void serialize(Restarter& res OPM_UNUSED)
    { /* not implemented */ }

    /*!
     * \brief This method restores the complete state of the tracer
     *        from disk.
     *
     * It is the inverse of the serialize() method.
     */
    template <class Restarter>
    void deserialize(Restarter& res OPM_UNUSED)
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
            Opm::decay<Scalar>(fs.saturation(tracerPhaseIdx_[tracerIdx]))
            *Opm::decay<Scalar>(fs.invB(tracerPhaseIdx_[tracerIdx]))
            *Opm::decay<Scalar>(intQuants.porosity());

        // avoid singular matrix if no water is present.
        phaseVolume = Opm::max(phaseVolume, 1e-10);

        if (std::is_same<LhsEval, Scalar>::value)
            tracerStorage = phaseVolume * tracerConcentrationInitial_[tracerIdx][globalDofIdx];
        else
            tracerStorage =
                    phaseVolume
                    * Opm::variable<LhsEval>(tracerConcentration_[tracerIdx][globalDofIdx][0], 0);
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

        const int tracerPhaseIdx = tracerPhaseIdx_[tracerIdx];

        unsigned upIdx = extQuants.upstreamIndex(tracerPhaseIdx);
        int globalUpIdx = elemCtx.globalSpaceIndex(upIdx, timeIdx);

        const auto& intQuants = elemCtx.intensiveQuantities(upIdx, timeIdx);
        const auto& fs = intQuants.fluidState();

        Scalar A = scvf.area();
        Scalar v = Opm::decay<Scalar>(extQuants.volumeFlux(tracerPhaseIdx));
        Scalar b = Opm::decay<Scalar>(fs.invB(tracerPhaseIdx_[tracerIdx]));
        Scalar c = tracerConcentration_[tracerIdx][globalUpIdx];

        if (inIdx == upIdx)
            tracerFlux = A*v*b*Opm::variable<TracerEvaluation>(c, 0);
        else
            tracerFlux = A*v*b*c;

    }

    bool linearSolve_(const TracerMatrix& M, TracerVector& x, TracerVector& b)
    {
#if ! DUNE_VERSION_NEWER(DUNE_COMMON, 2,7)
        Dune::FMatrixPrecision<Scalar>::set_singular_limit(1.e-30);
        Dune::FMatrixPrecision<Scalar>::set_absolute_limit(1.e-30);
#endif
        x = 0.0;
        Scalar tolerance = 1e-2;
        int maxIter = 100;

        int verbosity = 0;
        typedef Dune::BiCGSTABSolver<TracerVector> TracerSolver;
        typedef Dune::MatrixAdapter<TracerMatrix, TracerVector , TracerVector > TracerOperator;
        typedef Dune::SeqScalarProduct< TracerVector > TracerScalarProduct ;
#if DUNE_VERSION_NEWER(DUNE_ISTL, 2,6)
        typedef Dune::SeqILU< TracerMatrix, TracerVector, TracerVector  > TracerPreconditioner;
#else
        typedef Dune::SeqILUn< TracerMatrix, TracerVector, TracerVector  > TracerPreconditioner;
#endif

        TracerOperator tracerOperator(M);
        TracerScalarProduct tracerScalarProduct;
        TracerPreconditioner tracerPreconditioner(M, 0, 1); // results in ILU0

        TracerSolver solver (tracerOperator, tracerScalarProduct,
                             tracerPreconditioner, tolerance, maxIter,
                             verbosity);

        Dune::InverseOperatorResult result;
        solver.apply(x, b, result);

        // return the result of the solver
        return result.converged;
    }

    void linearize_(int tracerIdx)
    {
        (*tracerMatrix_) = 0.0;
        tracerResidual_ = 0.0;

        size_t numGridDof =  simulator_.model().numGridDof();
        std::vector<double> volumes(numGridDof, 0.0);
        ElementContext elemCtx(simulator_);
        auto elemIt = simulator_.gridView().template begin</*codim=*/0>();
        auto elemEndIt = simulator_.gridView().template end</*codim=*/0>();
        for (; elemIt != elemEndIt; ++ elemIt) {
            elemCtx.updateAll(*elemIt);

            Scalar extrusionFactor =
                    elemCtx.intensiveQuantities(/*dofIdx=*/ 0, /*timeIdx=*/0).extrusionFactor();
            Opm::Valgrind::CheckDefined(extrusionFactor);
            assert(Opm::isfinite(extrusionFactor));
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
                storageOfTimeIndex1 = storageOfTimeIndex1_[tracerIdx][I];
            else
                computeStorage_(storageOfTimeIndex1, elemCtx, 0, /*timIdx=*/1, tracerIdx);

            localStorage = (storageOfTimeIndex0 - storageOfTimeIndex1) * scvVolume/dt;
            tracerResidual_[I][0] += localStorage.value(); //residual + flux
            (*tracerMatrix_)[I][I][0][0] = localStorage.derivative(0);
            size_t numInteriorFaces = elemCtx.numInteriorFaces(/*timIdx=*/0);
            for (unsigned scvfIdx = 0; scvfIdx < numInteriorFaces; scvfIdx++) {
                TracerEvaluation flux;
                const auto& face = elemCtx.stencil(0).interiorFace(scvfIdx);
                unsigned j = face.exteriorIndex();
                unsigned J = elemCtx.globalSpaceIndex(/*dofIdx=*/ j, /*timIdx=*/0);
                computeFlux_(flux, elemCtx, scvfIdx, 0, tracerIdx);
                tracerResidual_[I][0] += flux.value(); //residual + flux
                (*tracerMatrix_)[J][I][0][0] = -flux.derivative(0);
                (*tracerMatrix_)[I][J][0][0] = flux.derivative(0);
            }

        }

        // Wells
        const int episodeIdx = simulator_.episodeIndex();
        const auto& wells = simulator_.vanguard().schedule().getWells2(episodeIdx);
        for (const auto& well : wells) {

            if (well.getStatus() == Opm::WellCommon::SHUT)
                continue;

            const double wtracer = well.getTracerProperties().getConcentration(tracerNames_[tracerIdx]);
            std::array<int, 3> cartesianCoordinate;
            for (auto& connection : well.getConnections()) {

                if (connection.state() == Opm::WellCompletion::SHUT)
                    continue;

                cartesianCoordinate[0] = connection.getI();
                cartesianCoordinate[1] = connection.getJ();
                cartesianCoordinate[2] = connection.getK();
                const size_t cartIdx = simulator_.vanguard().cartesianIndex(cartesianCoordinate);
                const int I = cartToGlobal_[cartIdx];
                Scalar rate = simulator_.problem().wellModel().well(well.name())->volumetricSurfaceRateForConnection(I, tracerPhaseIdx_[tracerIdx]);
                if (rate > 0)
                    tracerResidual_[I][0] -= rate*wtracer;
                else if (rate < 0)
                    tracerResidual_[I][0] -= rate*tracerConcentration_[tracerIdx][I];
            }
        }
    }

    Simulator& simulator_;

    std::vector<std::string> tracerNames_;
    std::vector<int> tracerPhaseIdx_;
    std::vector<Dune::BlockVector<Dune::FieldVector<Scalar, 1>>> tracerConcentration_;
    std::vector<Dune::BlockVector<Dune::FieldVector<Scalar, 1>>> tracerConcentrationInitial_;
    TracerMatrix *tracerMatrix_;
    TracerVector tracerResidual_;
    std::vector<int> cartToGlobal_;
    std::vector<Dune::BlockVector<Dune::FieldVector<Scalar, 1>>> storageOfTimeIndex1_;

};
} // namespace Ewoms

#endif
