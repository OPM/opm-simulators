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

#include <dune/common/version.hh>
#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

#include <opm/common/Exceptions.hpp>
#include <opm/common/TimingMacros.hpp>

#include <opm/grid/utility/SparseTable.hpp>

#include <opm/material/common/ConditionalStorage.hpp>

#include <opm/input/eclipse/EclipseState/Grid/FaceDir.hpp>
#include <opm/input/eclipse/Schedule/BCProp.hpp>

#include <opm/models/blackoil/blackoilproperties.hh>
#include <opm/models/common/multiphasebaseproperties.hh>
#include <opm/models/discretization/common/baseauxiliarymodule.hh>
#include <opm/models/discretization/common/fvbaseproperties.hh>
#include <opm/models/discretization/common/linearizationtype.hh>

#include <cassert>
#include <cstddef>
#include <exception>   // current_exception, rethrow_exception
#include <iostream>
#include <map>
#include <memory>
#include <numeric>
#include <set>
#include <stdexcept>
#include <unordered_map>
#include <vector>

namespace Opm::Parameters {

struct SeparateSparseSourceTerms { static constexpr bool value = false; };

} // namespace Opm::Parameters

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

    static constexpr bool linearizeNonLocalElements =
        getPropValue<TypeTag, Properties::LinearizeNonLocalElements>();
    static constexpr bool enableEnergy = getPropValue<TypeTag, Properties::EnableEnergy>();
    static constexpr bool enableDiffusion = getPropValue<TypeTag, Properties::EnableDiffusion>();
    static constexpr bool enableDispersion = getPropValue<TypeTag, Properties::EnableDispersion>();
    static constexpr bool enableMICP = getPropValue<TypeTag, Properties::EnableMICP>();

    // copying the linearizer is not a good idea
    TpfaLinearizer(const TpfaLinearizer&) = delete;
//! \endcond

public:
    TpfaLinearizer()
    {
        simulatorPtr_ = nullptr;
        separateSparseSourceTerms_ = Parameters::Get<Parameters::SeparateSparseSourceTerms>();
    }

    /*!
     * \brief Register all run-time parameters for the Jacobian linearizer.
     */
    static void registerParameters()
    {
        Parameters::Register<Parameters::SeparateSparseSourceTerms>
            ("Treat well source terms all in one go, instead of on a cell by cell basis.");
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
        int succeeded;
        try {
            linearizeDomain(fullDomain_);
            succeeded = 1;
        }
        catch (const std::exception& e) {
            std::cout << "rank " << simulator_().gridView().comm().rank()
                      << " caught an exception while linearizing:" << e.what()
                      << "\n"  << std::flush;
            succeeded = 0;
        }
        catch (...) {
            std::cout << "rank " << simulator_().gridView().comm().rank()
                      << " caught an exception while linearizing"
                      << "\n"  << std::flush;
            succeeded = 0;
        }
        succeeded = simulator_().gridView().comm().min(succeeded);

        if (!succeeded) {
            throw NumericalProblem("A process did not succeed in linearizing the system");
        }
    }

    /*!
     * \brief Linearize the part of the non-linear system of equations that is associated
     *        with a part of the spatial domain.
     *
     * That means that the Jacobian of the residual is assembled and the residual
     * is evaluated for the current solution, on the domain passed in as argument.
     *
     * The current state of affairs (esp. the previous and the current solutions) is
     * represented by the model object.
     */
    template <class SubDomainType>
    void linearizeDomain(const SubDomainType& domain)
    {
        OPM_TIMEBLOCK(linearizeDomain);
        // we defer the initialization of the Jacobian matrix until here because the
        // auxiliary modules usually assume the problem, model and grid to be fully
        // initialized...
        if (!jacobian_) {
            initFirstIteration_();
        }

        // Called here because it is no longer called from linearize_().
        if (domain.cells.size() == model_().numTotalDof()) {
            // We are on the full domain.
            resetSystem_();
        }
        else {
            resetSystem_(domain);
        }

        linearize_(domain);
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

            if (!succeeded) {
                throw NumericalProblem("linearization of an auxiliary equation failed");
            }
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

    void setLinearizationType(LinearizationType linearizationType)
    { linearizationType_ = linearizationType; }

    const LinearizationType& getLinearizationType() const
    { return linearizationType_; }

    /*!
     * \brief Return constant reference to the flowsInfo.
     *
     * (This object is only non-empty if the FLOWS keyword is true.)
     */
    const auto& getFlowsInfo() const
    { return flowsInfo_; }

    /*!
     * \brief Return constant reference to the floresInfo.
     *
     * (This object is only non-empty if the FLORES keyword is true.)
     */
    const auto& getFloresInfo() const
    { return floresInfo_; }

    /*!
     * \brief Return constant reference to the velocityInfo.
     *
     * (This object is only non-empty if the DISPERC keyword is true.)
     */
    const auto& getVelocityInfo() const
    { return velocityInfo_; }

    void updateDiscretizationParameters()
    {
        updateStoredTransmissibilities();
    }

    void updateBoundaryConditionData()
    {
        for (auto& bdyInfo : boundaryInfo_) {
            const auto [type, massrateAD] = problem_().boundaryCondition(bdyInfo.cell, bdyInfo.dir);

            // Strip the unnecessary (and zero anyway) derivatives off massrate.
            VectorBlock massrate(0.0);
            for (std::size_t ii = 0; ii < massrate.size(); ++ii) {
                massrate[ii] = massrateAD[ii].value();
            }
            if (type != BCType::NONE) {
                const auto& exFluidState = problem_().boundaryFluidState(bdyInfo.cell, bdyInfo.dir);
                bdyInfo.bcdata.type = type;
                bdyInfo.bcdata.massRate = massrate;
                bdyInfo.bcdata.exFluidState = exFluidState;
            }
        }
    }

    /*!
     * \brief Returns the map of constraint degrees of freedom.
     *
     * (This object is only non-empty if the EnableConstraints property is true.)
     */
    std::map<unsigned, Constraints> constraintsMap() const
    { return {}; }

    template <class SubDomainType>
    void resetSystem_(const SubDomainType& domain)
    {
        if (!jacobian_) {
            initFirstIteration_();
        }
        for (int globI : domain.cells) {
            residual_[globI] = 0.0;
            jacobian_->clearRow(globI, 0.0);
        }
    }

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
        using NeighborSet = std::set<unsigned>;
        std::vector<NeighborSet> sparsityPattern(model.numTotalDof());
        const Scalar gravity = problem_().gravity()[dimWorld - 1];
        unsigned numCells = model.numTotalDof();
        neighborInfo_.reserve(numCells, 6 * numCells);
        std::vector<NeighborInfo> loc_nbinfo;
        for (const auto& elem : elements(gridView_())) {
            stencil.update(elem);

            for (unsigned primaryDofIdx = 0; primaryDofIdx < stencil.numPrimaryDof(); ++primaryDofIdx) {
                const unsigned myIdx = stencil.globalSpaceIndex(primaryDofIdx);
                loc_nbinfo.resize(stencil.numDof() - 1); // Do not include the primary dof in neighborInfo_

                for (unsigned dofIdx = 0; dofIdx < stencil.numDof(); ++dofIdx) {
                    const unsigned neighborIdx = stencil.globalSpaceIndex(dofIdx);
                    sparsityPattern[myIdx].insert(neighborIdx);
                    if (dofIdx > 0) {
                        const Scalar trans = problem_().transmissibility(myIdx, neighborIdx);
                        const auto scvfIdx = dofIdx - 1;
                        const auto& scvf = stencil.interiorFace(scvfIdx);
                        const Scalar area = scvf.area();
                        const Scalar Vin = problem_().model().dofTotalVolume(myIdx);
                        const Scalar Vex = problem_().model().dofTotalVolume(neighborIdx);
                        const Scalar zIn = problem_().dofCenterDepth(myIdx);
                        const Scalar zEx = problem_().dofCenterDepth(neighborIdx);
                        const Scalar dZg = (zIn - zEx)*gravity;
                        const Scalar thpres = problem_().thresholdPressure(myIdx, neighborIdx);
                        const auto dirId = scvf.dirId();
                        auto faceDir = dirId < 0 ? FaceDir::DirEnum::Unknown
                                                 : FaceDir::FromIntersectionIndex(dirId);
                        ResidualNBInfo nbinfo{trans, area, thpres, dZg, faceDir, Vin, Vex, {}, {}, {}, {}};
                        if constexpr (enableEnergy) {
                            nbinfo.inAlpha = problem_().thermalHalfTransmissibility(myIdx, neighborIdx);
                            nbinfo.outAlpha = problem_().thermalHalfTransmissibility(neighborIdx, myIdx);
                        }
                        if constexpr (enableDiffusion) {
                            nbinfo.diffusivity = problem_().diffusivity(myIdx, neighborIdx);
                        }
                        if constexpr (enableDispersion) {
                            nbinfo.dispersivity = problem_().dispersivity(myIdx, neighborIdx);
                        }
                        loc_nbinfo[dofIdx - 1] = NeighborInfo{neighborIdx, nbinfo, nullptr};
                    }
                }
                neighborInfo_.appendRow(loc_nbinfo.begin(), loc_nbinfo.end());
                if (problem_().nonTrivialBoundaryConditions()) {
                    for (unsigned bfIndex = 0; bfIndex < stencil.numBoundaryFaces(); ++bfIndex) {
                        const auto& bf = stencil.boundaryFace(bfIndex);
                        const int dir_id = bf.dirId();
                        // not for NNCs
                        if (dir_id < 0) {
                            continue;
                        }
                        const auto [type, massrateAD] = problem_().boundaryCondition(myIdx, dir_id);
                        // Strip the unnecessary (and zero anyway) derivatives off massrate.
                        VectorBlock massrate(0.0);
                        for (std::size_t ii = 0; ii < massrate.size(); ++ii) {
                            massrate[ii] = massrateAD[ii].value();
                        }
                        const auto& exFluidState = problem_().boundaryFluidState(myIdx, dir_id);
                        BoundaryConditionData bcdata{type,
                                                     massrate,
                                                     exFluidState.pvtRegionIndex(),
                                                     bfIndex,
                                                     bf.area(),
                                                     bf.integrationPos()[dimWorld - 1],
                                                     exFluidState};
                        boundaryInfo_.push_back({myIdx, dir_id, bfIndex, bcdata});
                    }
                }
            }
        }

        // add the additional neighbors and degrees of freedom caused by the auxiliary
        // equations
        const std::size_t numAuxMod = model.numAuxiliaryModules();
        for (unsigned auxModIdx = 0; auxModIdx < numAuxMod; ++auxModIdx) {
            model.auxiliaryModule(auxModIdx)->addNeighbors(sparsityPattern);
        }

        // allocate raw matrix
        jacobian_ = std::make_unique<SparseMatrixAdapter>(simulator_());
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

        // Create dummy full domain.
        fullDomain_.cells.resize(numCells);
        std::iota(fullDomain_.cells.begin(), fullDomain_.cells.end(), 0);
    }

    // reset the global linear system of equations.
    void resetSystem_()
    {
        residual_ = 0.0;
        // zero all matrix entries
        jacobian_->clear();
    }

    // Initialize the flows, flores, and velocity sparse tables
    void createFlows_()
    {
        OPM_TIMEBLOCK(createFlows);
        // If FLOWS/FLORES is set in any RPTRST in the schedule, then we initializate the sparse tables
        // For now, do the same also if any block flows are requested (TODO: only save requested cells...)
        // If DISPERC is in the deck, we initialize the sparse table here as well.
        const bool anyFlows = simulator_().problem().eclWriter().outputModule().getFlows().anyFlows() ||
                              simulator_().problem().eclWriter().outputModule().getFlows().hasBlockFlows();
        const bool anyFlores = simulator_().problem().eclWriter().outputModule().getFlows().anyFlores();
        const bool dispersionActive = simulator_().vanguard().eclState().getSimulationConfig().rock_config().dispersion();
        if (((!anyFlows || !flowsInfo_.empty()) && (!anyFlores || !floresInfo_.empty())) && (!dispersionActive && !enableMICP)) {
            return;
        }
        const auto& model = model_();
        const auto& nncOutput = simulator_().problem().eclWriter().getOutputNnc();
        Stencil stencil(gridView_(), model_().dofMapper());
        const unsigned numCells = model.numTotalDof();
        std::unordered_multimap<int, std::pair<int, int>> nncIndices;
        std::vector<FlowInfo> loc_flinfo;
        std::vector<VelocityInfo> loc_vlinfo;
        unsigned int nncId = 0;
        VectorBlock flow(0.0);

        // Create a nnc structure to use fast lookup
        for (unsigned nncIdx = 0; nncIdx < nncOutput.size(); ++nncIdx) {
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
        if (dispersionActive || enableMICP) {
            velocityInfo_.reserve(numCells, 6 * numCells);
        }

        for (const auto& elem : elements(gridView_())) {
            stencil.update(elem);
            for (unsigned primaryDofIdx = 0; primaryDofIdx < stencil.numPrimaryDof(); ++primaryDofIdx) {
                const unsigned myIdx = stencil.globalSpaceIndex(primaryDofIdx);
                const int numFaces = stencil.numBoundaryFaces() + stencil.numInteriorFaces();
                loc_flinfo.resize(numFaces);
                loc_vlinfo.resize(stencil.numDof() - 1);

                for (unsigned dofIdx = 0; dofIdx < stencil.numDof(); ++dofIdx) {
                    const unsigned neighborIdx = stencil.globalSpaceIndex(dofIdx);
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
                        loc_vlinfo[dofIdx - 1] = VelocityInfo{flow};
                    }
                } 

                for (unsigned bdfIdx = 0; bdfIdx < stencil.numBoundaryFaces(); ++bdfIdx) {
                    const auto& scvf = stencil.boundaryFace(bdfIdx);
                    const int faceId = scvf.dirId();
                    loc_flinfo[stencil.numInteriorFaces() + bdfIdx] = FlowInfo{faceId, flow, nncId};
                }

                if (anyFlows) {
                    flowsInfo_.appendRow(loc_flinfo.begin(), loc_flinfo.end());
                }
                if (anyFlores) {
                    floresInfo_.appendRow(loc_flinfo.begin(), loc_flinfo.end());
                }
                if (dispersionActive || enableMICP) {
                    velocityInfo_.appendRow(loc_vlinfo.begin(), loc_vlinfo.end());
                }
            }
        }
    }

public:
    void setResAndJacobi(VectorBlock& res, MatrixBlock& bMat, const ADVectorBlock& resid) const
    {
        for (unsigned eqIdx = 0; eqIdx < numEq; ++eqIdx) {
            res[eqIdx] = resid[eqIdx].value();
        }

        for (unsigned eqIdx = 0; eqIdx < numEq; ++eqIdx) {
            for (unsigned pvIdx = 0; pvIdx < numEq; ++pvIdx) {
                // A[dofIdx][focusDofIdx][eqIdx][pvIdx] is the partial derivative of
                // the residual function 'eqIdx' for the degree of freedom 'dofIdx'
                // with regard to the focus variable 'pvIdx' of the degree of freedom
                // 'focusDofIdx'
                bMat[eqIdx][pvIdx] = resid[eqIdx].derivative(pvIdx);
            }
        }
    }

    void updateFlowsInfo()
    {
        OPM_TIMEBLOCK(updateFlows);
        const bool enableFlows = simulator_().problem().eclWriter().outputModule().getFlows().hasFlows() ||
                                 simulator_().problem().eclWriter().outputModule().getFlows().hasBlockFlows();
        const bool enableFlores = simulator_().problem().eclWriter().outputModule().getFlows().hasFlores();
        if (!enableFlows && !enableFlores) {
            return;
        }
        const unsigned int numCells = model_().numTotalDof();
#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (unsigned globI = 0; globI < numCells; ++globI) {
            OPM_TIMEBLOCK_LOCAL(linearizationForEachCell);
            const auto& nbInfos = neighborInfo_[globI];
            ADVectorBlock adres(0.0);
            ADVectorBlock darcyFlux(0.0);
            const IntensiveQuantities& intQuantsIn = model_().intensiveQuantities(globI, /*timeIdx*/ 0);
            // Flux term.
            {
                OPM_TIMEBLOCK_LOCAL(fluxCalculationForEachCell);
                short loc = 0;
                for (const auto& nbInfo : nbInfos) {
                    OPM_TIMEBLOCK_LOCAL(fluxCalculationForEachFace);
                    const unsigned globJ = nbInfo.neighbor;
                    assert(globJ != globI);
                    adres = 0.0;
                    darcyFlux = 0.0;
                    const IntensiveQuantities& intQuantsEx = model_().intensiveQuantities(globJ, /*timeIdx*/ 0);
                    LocalResidual::computeFlux(adres,darcyFlux, globI, globJ, intQuantsIn,
                                               intQuantsEx, nbInfo.res_nbinfo, problem_().moduleParams());
                    adres *= nbInfo.res_nbinfo.faceArea;
                    if (enableFlows) {
                        for (unsigned eqIdx = 0; eqIdx < numEq; ++eqIdx) {
                            flowsInfo_[globI][loc].flow[eqIdx] = adres[eqIdx].value();
                        }
                    }
                    if (enableFlores) {
                        for (unsigned eqIdx = 0; eqIdx < numEq; ++eqIdx) {
                            floresInfo_[globI][loc].flow[eqIdx] = darcyFlux[eqIdx].value();
                        }
                    }
                    ++loc;
                }
            }
        }

        // Boundary terms. Only looping over cells with nontrivial bcs.
        for (const auto& bdyInfo : boundaryInfo_) {
            if (bdyInfo.bcdata.type == BCType::NONE) {
                continue;
            }

            ADVectorBlock adres(0.0);
            const unsigned globI = bdyInfo.cell;
            const auto& nbInfos = neighborInfo_[globI];
            const IntensiveQuantities& insideIntQuants = model_().intensiveQuantities(globI, /*timeIdx*/ 0);
            LocalResidual::computeBoundaryFlux(adres, problem_(), bdyInfo.bcdata, insideIntQuants, globI);
            adres *= bdyInfo.bcdata.faceArea;
            const unsigned bfIndex = bdyInfo.bfIndex;
            if (enableFlows) {
                for (unsigned eqIdx = 0; eqIdx < numEq; ++eqIdx) {
                    flowsInfo_[globI][nbInfos.size() + bfIndex].flow[eqIdx] = adres[eqIdx].value();
                }
            }
            // TODO also store Flores?
        }
    }

private:
    template <class SubDomainType>
    void linearize_(const SubDomainType& domain)
    {
        // This check should be removed once this is addressed by
        // for example storing the previous timesteps' values for
        // rsmax (for DRSDT) and similar.
        if (!problem_().recycleFirstIterationStorage()) {
            if (!model_().storeIntensiveQuantities() && !model_().enableStorageCache()) {
                OPM_THROW(std::runtime_error, "Must have cached either IQs or storage when we cannot recycle.");
            }
        }

        OPM_TIMEBLOCK(linearize);

        // We do not call resetSystem_() here, since that will set
        // the full system to zero, not just our part.
        // Instead, that must be called before starting the linearization.
        const bool dispersionActive = simulator_().vanguard().eclState().getSimulationConfig().rock_config().dispersion();
        const unsigned int numCells = domain.cells.size();
        const bool on_full_domain = (numCells == model_().numTotalDof());

#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (unsigned ii = 0; ii < numCells; ++ii) {
            OPM_TIMEBLOCK_LOCAL(linearizationForEachCell);
            const unsigned globI = domain.cells[ii];
            const auto& nbInfos = neighborInfo_[globI];
            VectorBlock res(0.0);
            MatrixBlock bMat(0.0);
            ADVectorBlock adres(0.0);
            ADVectorBlock darcyFlux(0.0);
            const IntensiveQuantities& intQuantsIn = model_().intensiveQuantities(globI, /*timeIdx*/ 0);

            // Flux term.
            {
                OPM_TIMEBLOCK_LOCAL(fluxCalculationForEachCell);
                short loc = 0;
                for (const auto& nbInfo : nbInfos) {
                    OPM_TIMEBLOCK_LOCAL(fluxCalculationForEachFace);
                    const unsigned globJ = nbInfo.neighbor;
                    assert(globJ != globI);
                    res = 0.0;
                    bMat = 0.0;
                    adres = 0.0;
                    darcyFlux = 0.0;
                    const IntensiveQuantities& intQuantsEx = model_().intensiveQuantities(globJ, /*timeIdx*/ 0);
                    LocalResidual::computeFlux(adres,darcyFlux, globI, globJ, intQuantsIn, intQuantsEx,
                                               nbInfo.res_nbinfo,  problem_().moduleParams());
                    adres *= nbInfo.res_nbinfo.faceArea;
                    if (dispersionActive || enableMICP) {
                        for (unsigned phaseIdx = 0; phaseIdx < numEq; ++phaseIdx) {
                            velocityInfo_[globI][loc].velocity[phaseIdx] =
                                darcyFlux[phaseIdx].value() / nbInfo.res_nbinfo.faceArea;
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
            const double dt = simulator_().timeStepSize();
            const double volume = model_().dofTotalVolume(globI);
            const Scalar storefac = volume / dt;
            adres = 0.0;
            {
                OPM_TIMEBLOCK_LOCAL(computeStorage);
                LocalResidual::computeStorage(adres, intQuantsIn);
            }
            setResAndJacobi(res, bMat, adres);
            // Either use cached storage term, or compute it on the fly.
            if (model_().enableStorageCache()) {
                // The cached storage for timeIdx 0 (current time) is not
                // used, but after storage cache is shifted at the end of the
                // timestep, it will become cached storage for timeIdx 1.
                model_().updateCachedStorage(globI, /*timeIdx=*/0, res);
                if (model_().newtonMethod().numIterations() == 0) {
                    // Need to update the storage cache.
                    if (problem_().recycleFirstIterationStorage()) {
                        // Assumes nothing have changed in the system which
                        // affects masses calculated from primary variables.
                        if (on_full_domain) {
                            // This is to avoid resetting the start-of-step storage
                            // to incorrect numbers when we do local solves, where the iteration
                            // number will start from 0, but the starting state may not be identical
                            // to the start-of-step state.
                            // Note that a full assembly must be done before local solves
                            // otherwise this will be left un-updated.
                            model_().updateCachedStorage(globI, /*timeIdx=*/1, res);
                        }
                    }
                    else {
                        Dune::FieldVector<Scalar, numEq> tmp;
                        const IntensiveQuantities intQuantOld = model_().intensiveQuantities(globI, 1);
                        LocalResidual::computeStorage(tmp, intQuantOld);
                        model_().updateCachedStorage(globI, /*timeIdx=*/1, tmp);
                    }
                }
                res -= model_().cachedStorage(globI, 1);
            }
            else {
                OPM_TIMEBLOCK_LOCAL(computeStorage0);
                Dune::FieldVector<Scalar, numEq> tmp;
                const IntensiveQuantities intQuantOld = model_().intensiveQuantities(globI, 1);
                LocalResidual::computeStorage(tmp, intQuantOld);
                // assume volume do not change
                res -= tmp;
            }
            res *= storefac;
            bMat *= storefac;
            residual_[globI] += res;
            //SparseAdapter syntax: jacobian_->addToBlock(globI, globI, bMat);
            *diagMatAddress_[globI] += bMat;

            // Cell-wise source terms.
            // This will include well sources if SeparateSparseSourceTerms is false.
            res = 0.0;
            bMat = 0.0;
            adres = 0.0;
            if (separateSparseSourceTerms_) {
                LocalResidual::computeSourceDense(adres, problem_(), intQuantsIn, globI, 0);
            }
            else {
                LocalResidual::computeSource(adres, problem_(), intQuantsIn, globI, 0);
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
            if (bdyInfo.bcdata.type == BCType::NONE) {
                continue;
            }

            VectorBlock res(0.0);
            MatrixBlock bMat(0.0);
            ADVectorBlock adres(0.0);
            const unsigned globI = bdyInfo.cell;
            const IntensiveQuantities& insideIntQuants = model_().intensiveQuantities(globI, /*timeIdx*/ 0);
            LocalResidual::computeBoundaryFlux(adres, problem_(), bdyInfo.bcdata, insideIntQuants, globI);
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

        const unsigned numCells = model_().numTotalDof();
#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (unsigned globI = 0; globI < numCells; globI++) {
            auto nbInfos = neighborInfo_[globI]; // nbInfos will be a SparseTable<...>::mutable_iterator_range.
            for (auto& nbInfo : nbInfos) {
                const unsigned globJ = nbInfo.neighbor;
                nbInfo.res_nbinfo.trans = problem_().transmissibility(globI, globJ);
            }
        }
    }

    Simulator* simulatorPtr_{};

    // the jacobian matrix
    std::unique_ptr<SparseMatrixAdapter> jacobian_{};

    // the right-hand side
    GlobalEqVector residual_;

    LinearizationType linearizationType_{};

    using ResidualNBInfo = typename LocalResidual::ResidualNBInfo;
    struct NeighborInfo
    {
        unsigned int neighbor;
        ResidualNBInfo res_nbinfo;
        MatrixBlock* matBlockAddress;
    };
    SparseTable<NeighborInfo> neighborInfo_{};
    std::vector<MatrixBlock*> diagMatAddress_{};

    struct FlowInfo
    {
        int faceId;
        VectorBlock flow;
        unsigned int nncId;
    };
    SparseTable<FlowInfo> flowsInfo_;
    SparseTable<FlowInfo> floresInfo_;

    struct VelocityInfo
    {
        VectorBlock velocity;
    };
    SparseTable<VelocityInfo> velocityInfo_;

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
        int dir;
        unsigned int bfIndex;
        BoundaryConditionData bcdata;
    };
    std::vector<BoundaryInfo> boundaryInfo_;

    bool separateSparseSourceTerms_ = false;

    struct FullDomain
    {
        std::vector<int> cells;
        std::vector<bool> interior;
    };
    FullDomain fullDomain_;
};

} // namespace Opm

#endif // TPFA_LINEARIZER
