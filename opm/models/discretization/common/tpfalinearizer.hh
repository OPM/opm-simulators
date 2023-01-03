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
 * \copydoc Opm::FvBaseLinearizer
 */
#ifndef TPFA_LINEARIZER_HH
#define TPFA_LINEARIZER_HH

#include "fvbaseproperties.hh"
#include "linearizationtype.hh"

#include <opm/common/Exceptions.hpp>
#include <opm/common/TimingMacros.hpp>

#include <opm/models/discretization/common/baseauxiliarymodule.hh>

#include <opm/grid/utility/SparseTable.hpp>
#include <opm/input/eclipse/EclipseState/Grid/FaceDir.hpp>

#include <dune/common/version.hh>
#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

#include <type_traits>
#include <iostream>
#include <vector>
#include <thread>
#include <set>
#include <exception>   // current_exception, rethrow_exception
#include <mutex>

namespace Opm::Properties {
    template<class TypeTag, class MyTypeTag>
    struct SeparateSparseSourceTerms {
        using type = bool;
        static constexpr type value = true;
    };
}

namespace Opm {

// forward declarations
template<class TypeTag>
class EcfvDiscretization;

/*!
 * \ingroup FiniteVolumeDiscretizations
 *
 * \brief The common code for the linearizers of non-linear systems of equations
 *
 * This class assumes that these system of equations to be linearized are stemming from
 * models that use an finite volume scheme for spatial discretization and an Euler
 * scheme for time discretization.
 */
template<class TypeTag>
class TpfaLinearizer
{
//! \cond SKIP_THIS
    using Model = GetPropType<TypeTag, Properties::Model>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using Simulator = GetPropType<TypeTag, Properties::Simulator>;
    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Evaluation = GetPropType<TypeTag, Properties::Evaluation>;

    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    using GlobalEqVector = GetPropType<TypeTag, Properties::GlobalEqVector>;
    using SparseMatrixAdapter = GetPropType<TypeTag, Properties::SparseMatrixAdapter>;
    using EqVector = GetPropType<TypeTag, Properties::EqVector>;
    using Constraints = GetPropType<TypeTag, Properties::Constraints>;
    using Stencil = GetPropType<TypeTag, Properties::Stencil>;
    using LocalResidual = GetPropType<TypeTag, Properties::LocalResidual>;
    using IntensiveQuantities = GetPropType<TypeTag, Properties::IntensiveQuantities>;

    using Element = typename GridView::template Codim<0>::Entity;
    using ElementIterator = typename GridView::template Codim<0>::Iterator;

    using Vector = GlobalEqVector;

    enum { numEq = getPropValue<TypeTag, Properties::NumEq>() };
    enum { historySize = getPropValue<TypeTag, Properties::TimeDiscHistorySize>() };
    enum { dimWorld = GridView::dimensionworld };

    using MatrixBlock = typename SparseMatrixAdapter::MatrixBlock;
    using VectorBlock = Dune::FieldVector<Scalar, numEq>;
    using ADVectorBlock = GetPropType<TypeTag, Properties::RateVector>;

    static const bool linearizeNonLocalElements = getPropValue<TypeTag, Properties::LinearizeNonLocalElements>();

    // copying the linearizer is not a good idea
    TpfaLinearizer(const TpfaLinearizer&);
//! \endcond

public:
    TpfaLinearizer()
        : jacobian_()
    {
        simulatorPtr_ = 0;
        separateSparseSourceTerms_ = EWOMS_GET_PARAM(TypeTag, bool, SeparateSparseSourceTerms);
    }

    ~TpfaLinearizer()
    {
    }

    /*!
     * \brief Register all run-time parameters for the Jacobian linearizer.
     */
    static void registerParameters()
    {
        EWOMS_REGISTER_PARAM(TypeTag, bool, SeparateSparseSourceTerms,
                             "Treat well source terms all in one go, instead of on a cell by cell basis.");
    }

    /*!
     * \brief Initialize the linearizer.
     *
     * At this point we can assume that all objects in the simulator
     * have been allocated. We cannot assume that they are fully
     * initialized, though.
     *
     * \copydetails Doxygen::simulatorParam
     */
    void init(Simulator& simulator)
    {
        simulatorPtr_ = &simulator;
        eraseMatrix();
    }

    /*!
     * \brief Causes the Jacobian matrix to be recreated from scratch before the next
     *        iteration.
     *
     * This method is usally called if the sparsity pattern has changed for some
     * reason. (e.g. by modifications of the grid or changes of the auxiliary equations.)
     */
    void eraseMatrix()
    {
        jacobian_.reset();
    }

    /*!
     * \brief Linearize the full system of non-linear equations.
     *
     * The linearizationType() controls the scheme used and the focus
     * time index. The default is fully implicit scheme, and focus index
     * equal to 0, i.e. current time (end of step).
     *
     * This linearizes the spatial domain and all auxiliary equations.
     */
    void linearize()
    {
        linearizeDomain();
        linearizeAuxiliaryEquations();
    }

    /*!
     * \brief Linearize the part of the non-linear system of equations that is associated
     *        with the spatial domain.
     *
     * That means that the global Jacobian of the residual is assembled and the residual
     * is evaluated for the current solution.
     *
     * The current state of affairs (esp. the previous and the current solutions) is
     * represented by the model object.
     */
    void linearizeDomain()
    {
        OPM_TIMEBLOCK(linearizeDomain);
        // we defer the initialization of the Jacobian matrix until here because the
        // auxiliary modules usually assume the problem, model and grid to be fully
        // initialized...
        if (!jacobian_)
            initFirstIteration_();

        int succeeded;
        try {
            linearize_();
            succeeded = 1;
        }
        catch (const std::exception& e)
        {
            std::cout << "rank " << simulator_().gridView().comm().rank()
                      << " caught an exception while linearizing:" << e.what()
                      << "\n"  << std::flush;
            succeeded = 0;
        }
        catch (...)
        {
            std::cout << "rank " << simulator_().gridView().comm().rank()
                      << " caught an exception while linearizing"
                      << "\n"  << std::flush;
            succeeded = 0;
        }
        succeeded = gridView_().comm().min(succeeded);

        if (!succeeded)
            throw NumericalProblem("A process did not succeed in linearizing the system");
    }

    void finalize()
    { jacobian_->finalize(); }

    /*!
     * \brief Linearize the part of the non-linear system of equations that is associated
     *        with the spatial domain.
     */
    void linearizeAuxiliaryEquations()
    {
        OPM_TIMEBLOCK(linearizeAuxilaryEquations);
        // flush possible local caches into matrix structure
        jacobian_->commit();

        auto& model = model_();
        const auto& comm = simulator_().gridView().comm();
        for (unsigned auxModIdx = 0; auxModIdx < model.numAuxiliaryModules(); ++auxModIdx) {
            bool succeeded = true;
            try {
                model.auxiliaryModule(auxModIdx)->linearize(*jacobian_, residual_);
            }
            catch (const std::exception& e) {
                succeeded = false;

                std::cout << "rank " << simulator_().gridView().comm().rank()
                          << " caught an exception while linearizing:" << e.what()
                          << "\n"  << std::flush;
            }

            succeeded = comm.min(succeeded);

            if (!succeeded)
                throw NumericalProblem("linearization of an auxiliary equation failed");
        }
    }

    /*!
     * \brief Return constant reference to global Jacobian matrix backend.
     */
    const SparseMatrixAdapter& jacobian() const
    { return *jacobian_; }

    SparseMatrixAdapter& jacobian()
    { return *jacobian_; }

    /*!
     * \brief Return constant reference to global residual vector.
     */
    const GlobalEqVector& residual() const
    { return residual_; }

    GlobalEqVector& residual()
    { return residual_; }

    void setLinearizationType(LinearizationType linearizationType){
        linearizationType_ = linearizationType;
    };

    const LinearizationType& getLinearizationType() const{
        return linearizationType_;
    };

    /*!
     * \brief Return constant reference to the flowsInfo.
     *
     * (This object is only non-empty if the FLOWS keyword is true.)
     */
    const auto& getFlowsInfo() const{

        return flowsInfo_;
    }   

    /*!
     * \brief Return constant reference to the floresInfo.
     *
     * (This object is only non-empty if the FLORES keyword is true.)
     */
    const auto& getFloresInfo() const{

        return floresInfo_;
    }

    void updateDiscretizationParameters()
    {
        updateStoredTransmissibilities();
    }

    /*!
     * \brief Returns the map of constraint degrees of freedom.
     *
     * (This object is only non-empty if the EnableConstraints property is true.)
     */
    const std::map<unsigned, Constraints> constraintsMap() const
    { return {}; }

private:
    Simulator& simulator_()
    { return *simulatorPtr_; }
    const Simulator& simulator_() const
    { return *simulatorPtr_; }

    Problem& problem_()
    { return simulator_().problem(); }
    const Problem& problem_() const
    { return simulator_().problem(); }

    Model& model_()
    { return simulator_().model(); }
    const Model& model_() const
    { return simulator_().model(); }

    const GridView& gridView_() const
    { return problem_().gridView(); }

    void initFirstIteration_()
    {
        // initialize the BCRS matrix for the Jacobian of the residual function
        createMatrix_();

        // initialize the Jacobian matrix and the vector for the residual function
        residual_.resize(model_().numTotalDof());
        resetSystem_();

        // initialize the sparse tables for Flows and Flores
        createFlows_();
    }

    // Construct the BCRS matrix for the Jacobian of the residual function
    void createMatrix_()
    {
        OPM_TIMEBLOCK(createMatrix);
        if (!neighborInfo_.empty()) {
            // It is ok to call this function multiple times, but it
            // should not do anything if already called.
            return;
        }
        const auto& model = model_();
        Stencil stencil(gridView_(), model_().dofMapper());

        // for the main model, find out the global indices of the neighboring degrees of
        // freedom of each primary degree of freedom
        using NeighborSet = std::set< unsigned >;
        std::vector<NeighborSet> sparsityPattern(model.numTotalDof());

        unsigned numCells = model.numTotalDof();
        neighborInfo_.reserve(numCells, 6 * numCells);
        std::vector<NeighborInfo> loc_nbinfo;
        const auto& materialLawManager = problem_().materialLawManager();
        using FaceDirection = FaceDir::DirEnum;
        for (const auto& elem : elements(gridView_())) {
            stencil.update(elem);

            for (unsigned primaryDofIdx = 0; primaryDofIdx < stencil.numPrimaryDof(); ++primaryDofIdx) {
                unsigned myIdx = stencil.globalSpaceIndex(primaryDofIdx);
                loc_nbinfo.resize(stencil.numDof() - 1); // Do not include the primary dof in neighborInfo_

                for (unsigned dofIdx = 0; dofIdx < stencil.numDof(); ++dofIdx) {
                    unsigned neighborIdx = stencil.globalSpaceIndex(dofIdx);
                    sparsityPattern[myIdx].insert(neighborIdx);
                    if (dofIdx > 0) {
                        const double trans = problem_().transmissibility(myIdx, neighborIdx);
                        const auto scvfIdx = dofIdx - 1;
                        const auto& scvf = stencil.interiorFace(scvfIdx);
                        const double area = scvf.area();
                        FaceDirection dirId = FaceDirection::Unknown;
                        if (materialLawManager->hasDirectionalRelperms()) {
                            dirId = scvf.faceDirFromDirId();
                        }
                        loc_nbinfo[dofIdx - 1] = NeighborInfo{neighborIdx, trans, area, dirId, nullptr};
                    }
                }
                neighborInfo_.appendRow(loc_nbinfo.begin(), loc_nbinfo.end());
                for (unsigned bfIndex = 0; bfIndex < stencil.numBoundaryFaces(); ++bfIndex) {
                    const auto& bf = stencil.boundaryFace(bfIndex);
                    const int dir_id = bf.dirId();
                    const auto [free, massrateAD] = problem_().boundaryCondition(myIdx, dir_id);
                    // Strip the unnecessary (and zero anyway) derivatives off massrate.
                    VectorBlock massrate(0.0);
                    for (size_t ii = 0; ii < massrate.size(); ++ii) {
                        massrate[ii] = massrateAD[ii].value();
                    }
                    const bool nonzero_massrate = massrate != VectorBlock(0.0);
                    if (free || nonzero_massrate) {
                        const auto& exFluidState = problem_().boundaryFluidState(myIdx, dir_id);
                        BoundaryConditionData bcdata{free ? BCType::FREE : BCType::RATE,
                                                     massrate,
                                                     exFluidState.pvtRegionIndex(),
                                                     bfIndex,
                                                     bf.area(),
                                                     bf.integrationPos()[dimWorld - 1],
                                                     exFluidState};
                        boundaryInfo_.push_back({myIdx, bcdata});
                    }
                }
            }
        }

        // add the additional neighbors and degrees of freedom caused by the auxiliary
        // equations
        size_t numAuxMod = model.numAuxiliaryModules();
        for (unsigned auxModIdx = 0; auxModIdx < numAuxMod; ++auxModIdx)
            model.auxiliaryModule(auxModIdx)->addNeighbors(sparsityPattern);

        // allocate raw matrix
        jacobian_.reset(new SparseMatrixAdapter(simulator_()));
        diagMatAddress_.resize(numCells);
        // create matrix structure based on sparsity pattern
        jacobian_->reserve(sparsityPattern);
        for (unsigned globI = 0; globI < numCells; globI++) {
            const auto& nbInfos = neighborInfo_[globI];
            diagMatAddress_[globI] = jacobian_->blockAddress(globI, globI);
            for (auto& nbInfo : nbInfos) {
                nbInfo.matBlockAddress = jacobian_->blockAddress(nbInfo.neighbor, globI);
            }
        }
    }

    // reset the global linear system of equations.
    void resetSystem_()
    {
        residual_ = 0.0;
        // zero all matrix entries
        jacobian_->clear();
    }

    // Initialize the flows and flores sparse tables
    void createFlows_()
    {
        OPM_TIMEBLOCK(createFlows);
        // If FLOWS/FLORES is set in any RPTRST in the schedule, then we initializate the sparse tables
        const bool anyFlows = simulator_().problem().eclWriter()->eclOutputModule().anyFlows();
        const bool anyFlores = simulator_().problem().eclWriter()->eclOutputModule().anyFlores();
        if ((!anyFlows || !flowsInfo_.empty())  && (!anyFlores || !floresInfo_.empty())) {
            return;
        }
        const auto& model = model_();
        const auto& nncOutput = simulator_().problem().eclWriter()->getOutputNnc();
        Stencil stencil(gridView_(), model_().dofMapper());
        unsigned numCells = model.numTotalDof();
        std::unordered_multimap<int, std::pair<int, int>> nncIndices;
        std::vector<FlowInfo> loc_flinfo;
        unsigned int nncId = 0;
        VectorBlock flow(0.0);

        // Create a nnc structure to use fast lookup
        for (unsigned int nncIdx = 0; nncIdx < nncOutput.size(); ++nncIdx) {
            const int ci1 = nncOutput[nncIdx].cell1;
            const int ci2 = nncOutput[nncIdx].cell2;
            nncIndices.emplace(ci1, std::make_pair(ci2, nncIdx));
        }

        if (anyFlows) {
            flowsInfo_.reserve(numCells, 6 * numCells);
        }
        if (anyFlores) {
            floresInfo_.reserve(numCells, 6 * numCells);
        }

        for (const auto& elem : elements(gridView_())) {
            stencil.update(elem);
            for (unsigned primaryDofIdx = 0; primaryDofIdx < stencil.numPrimaryDof(); ++primaryDofIdx) {
                unsigned myIdx = stencil.globalSpaceIndex(primaryDofIdx);
                loc_flinfo.resize(stencil.numDof() - 1);
                for (unsigned dofIdx = 0; dofIdx < stencil.numDof(); ++dofIdx) {
                    unsigned neighborIdx = stencil.globalSpaceIndex(dofIdx);
                    if (dofIdx > 0) {
                        const auto scvfIdx = dofIdx - 1;
                        const auto& scvf = stencil.interiorFace(scvfIdx);
                        int faceId = scvf.dirId();
                        const int cartMyIdx = simulator_().vanguard().cartesianIndex(myIdx);
                        const int cartNeighborIdx = simulator_().vanguard().cartesianIndex(neighborIdx);
                        const auto& range = nncIndices.equal_range(cartMyIdx);
                        for (auto it = range.first; it != range.second; ++it) { 
                            if (it->second.first == cartNeighborIdx){
                                // -1 gives problem since is used for the nncInput from the deck
                                faceId = -2;
                                // the index is stored to be used for writting the outputs
                                nncId = it->second.second;
                            }
                        }
                        loc_flinfo[dofIdx - 1] = FlowInfo{faceId, flow, nncId};
                    }
                }
                if (anyFlows) {
                    flowsInfo_.appendRow(loc_flinfo.begin(), loc_flinfo.end());
                }
                if (anyFlores) {
                    floresInfo_.appendRow(loc_flinfo.begin(), loc_flinfo.end());
                }
            }
        }
    }

public:
    void setResAndJacobi(VectorBlock& res, MatrixBlock& bMat, const ADVectorBlock& resid) const
    {
        for (unsigned eqIdx = 0; eqIdx < numEq; eqIdx++)
            res[eqIdx] = resid[eqIdx].value();

        for (unsigned eqIdx = 0; eqIdx < numEq; eqIdx++) {
            for (unsigned pvIdx = 0; pvIdx < numEq; pvIdx++) {
                // A[dofIdx][focusDofIdx][eqIdx][pvIdx] is the partial derivative of
                // the residual function 'eqIdx' for the degree of freedom 'dofIdx'
                // with regard to the focus variable 'pvIdx' of the degree of freedom
                // 'focusDofIdx'
                bMat[eqIdx][pvIdx] = resid[eqIdx].derivative(pvIdx);
            }
        }
    }

private:
    void linearize_()
    {
        OPM_TIMEBLOCK(linearize);
        resetSystem_();
        unsigned numCells = model_().numTotalDof();
        const bool& enableFlows = simulator_().problem().eclWriter()->eclOutputModule().hasFlows();
        const bool& enableFlores = simulator_().problem().eclWriter()->eclOutputModule().hasFlores();
#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (unsigned globI = 0; globI < numCells; globI++) {
            OPM_TIMEBLOCK_LOCAL(linearizationForEachCell);
            const auto& nbInfos = neighborInfo_[globI]; // this is a set but should maybe be changed
            VectorBlock res(0.0);
            MatrixBlock bMat(0.0);
            ADVectorBlock adres(0.0);
            ADVectorBlock darcyFlux(0.0);
            const IntensiveQuantities* intQuantsInP = model_().cachedIntensiveQuantities(globI, /*timeIdx*/ 0);
            if (intQuantsInP == nullptr) {
                throw std::logic_error("Missing updated intensive quantities for cell " + std::to_string(globI));
            }
            const IntensiveQuantities& intQuantsIn = *intQuantsInP;

            // Flux term.
            {
            OPM_TIMEBLOCK_LOCAL(fluxCalculationForEachCell);    
            short loc = 0;
            for (const auto& nbInfo : nbInfos) {
                OPM_TIMEBLOCK_LOCAL(fluxCalculationForEachFace);
                unsigned globJ = nbInfo.neighbor;
                assert(globJ != globI);
                res = 0.0;
                bMat = 0.0;
                adres = 0.0;
                darcyFlux = 0.0;
                const IntensiveQuantities* intQuantsExP = model_().cachedIntensiveQuantities(globJ, /*timeIdx*/ 0);
                if (intQuantsExP == nullptr) {
                    throw std::logic_error("Missing updated intensive quantities for cell " + std::to_string(globJ) + " when assembling fluxes for cell " + std::to_string(globI));
                }
                const IntensiveQuantities& intQuantsEx = *intQuantsExP;
                LocalResidual::computeFlux(
                       adres, darcyFlux, problem_(), globI, globJ, intQuantsIn, intQuantsEx,
                           nbInfo.trans, nbInfo.faceArea, nbInfo.faceDirection);
                adres *= nbInfo.faceArea;
                if (enableFlows) {
                    for (unsigned phaseIdx = 0; phaseIdx < numEq; ++ phaseIdx) {
                        flowsInfo_[globI][loc].flow[phaseIdx] = adres[phaseIdx].value();
                    }
                }
                if (enableFlores) {
                    for (unsigned phaseIdx = 0; phaseIdx < numEq; ++ phaseIdx) {
                        floresInfo_[globI][loc].flow[phaseIdx] = darcyFlux[phaseIdx].value();
                    }
                }
                setResAndJacobi(res, bMat, adres);
                residual_[globI] += res;
                //SparseAdapter syntax:  jacobian_->addToBlock(globI, globI, bMat);
                *diagMatAddress_[globI] += bMat;
                bMat *= -1.0;
                //SparseAdapter syntax: jacobian_->addToBlock(globJ, globI, bMat);
                *nbInfo.matBlockAddress += bMat;
                ++loc;
            }
            }

            // Accumulation term.
            double dt = simulator_().timeStepSize();
            double volume = model_().dofTotalVolume(globI);
            Scalar storefac = volume / dt;
            adres = 0.0;
            {
                OPM_TIMEBLOCK_LOCAL(computeStorage);
                LocalResidual::computeStorage(adres, intQuantsIn);
            }
            setResAndJacobi(res, bMat, adres);
            // TODO: check recycleFirst etc.
            // first we use it as storage cache
            if (model_().newtonMethod().numIterations() == 0) {
                model_().updateCachedStorage(globI, /*timeIdx=*/1, res);
            }
            res -= model_().cachedStorage(globI, 1);
            res *= storefac;
            bMat *= storefac;
            // residual_[globI] -= model_().cachedStorage(globI, 1); //*storefac;
            residual_[globI] += res;
            //SparseAdapter syntax: jacobian_->addToBlock(globI, globI, bMat);
            *diagMatAddress_[globI] += bMat;

            // Cell-wise source terms.
            // This will include well sources if SeparateSparseSourceTerms is false.
            res = 0.0;
            bMat = 0.0;
            adres = 0.0;
            if (separateSparseSourceTerms_) {
                LocalResidual::computeSourceDense(adres, problem_(), globI, 0);
            } else {
                LocalResidual::computeSource(adres, problem_(), globI, 0);
            }
            adres *= -volume;
            setResAndJacobi(res, bMat, adres);
            residual_[globI] += res;
            //SparseAdapter syntax: jacobian_->addToBlock(globI, globI, bMat);
            *diagMatAddress_[globI] += bMat;
        } // end of loop for cell globI.

        // Add sparse source terms. For now only wells.
        if (separateSparseSourceTerms_) {
            problem_().wellModel().addReservoirSourceTerms(residual_, diagMatAddress_);
        }

        // Boundary terms. Only looping over cells with nontrivial bcs.
        for (const auto& bdyInfo : boundaryInfo_) {
            VectorBlock res(0.0);
            MatrixBlock bMat(0.0);
            ADVectorBlock adres(0.0);
            const unsigned globI = bdyInfo.cell;
            const IntensiveQuantities* insideIntQuants = model_().cachedIntensiveQuantities(globI, /*timeIdx*/ 0);
            if (insideIntQuants == nullptr) {
                throw std::logic_error("Missing updated intensive quantities for cell " + std::to_string(globI));
            }
            LocalResidual::computeBoundaryFlux(adres, problem_(), bdyInfo.bcdata, *insideIntQuants, globI);
            adres *= bdyInfo.bcdata.faceArea;
            setResAndJacobi(res, bMat, adres);
            residual_[globI] += res;
            ////SparseAdapter syntax: jacobian_->addToBlock(globI, globI, bMat);
            *diagMatAddress_[globI] += bMat;
        }
    }

    void updateStoredTransmissibilities()
    {
        if (neighborInfo_.empty()) {
            // This function was called before createMatrix_() was called.
            // We call initFirstIteration_(), not createMatrix_(), because
            // that will also initialize the residual consistently.
            initFirstIteration_();
        }
        unsigned numCells = model_().numTotalDof();
#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (unsigned globI = 0; globI < numCells; globI++) {
            auto nbInfos = neighborInfo_[globI]; // nbInfos will be a SparseTable<...>::mutable_iterator_range.
            for (auto& nbInfo : nbInfos) {
                unsigned globJ = nbInfo.neighbor;
                nbInfo.trans = problem_().transmissibility(globI, globJ);
            }
        }
    }


    Simulator *simulatorPtr_;

    // the jacobian matrix
    std::unique_ptr<SparseMatrixAdapter> jacobian_;

    // the right-hand side
    GlobalEqVector residual_;

    LinearizationType linearizationType_;

    struct NeighborInfo
    {
        unsigned int neighbor;
        double trans;
        double faceArea;
        FaceDir::DirEnum faceDirection;
        MatrixBlock* matBlockAddress;
    };
    SparseTable<NeighborInfo> neighborInfo_;
    std::vector<MatrixBlock*> diagMatAddress_;

    struct FlowInfo
    {
        int faceId;
        VectorBlock flow;
        unsigned int nncId;
    };
    SparseTable<FlowInfo> flowsInfo_;
    SparseTable<FlowInfo> floresInfo_;

    using ScalarFluidState = typename IntensiveQuantities::ScalarFluidState;
    struct BoundaryConditionData
    {
        BCType type;
        VectorBlock massRate;
        unsigned pvtRegionIdx;
        unsigned boundaryFaceIndex;
        double faceArea;
        double faceZCoord;
        ScalarFluidState exFluidState;
    };
    struct BoundaryInfo
    {
        unsigned int cell;
        BoundaryConditionData bcdata;
    };
    std::vector<BoundaryInfo> boundaryInfo_;
    bool separateSparseSourceTerms_ = false;
};

} // namespace Opm

#endif // TPFA_LINEARIZER
