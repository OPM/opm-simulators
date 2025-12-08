// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  Copyright 2025 NORCE AS

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
#ifndef TPSA_LINEARIZER_HPP
#define TPSA_LINEARIZER_HPP

#include <dune/common/fvector.hh>

#include <opm/common/TimingMacros.hpp>

#include <opm/grid/utility/SparseTable.hpp>

#include <opm/input/eclipse/Schedule/BCProp.hpp>

#include <opm/material/materialstates/MaterialStateTPSA.hpp>

#include <opm/models/discretization/common/linearizationtype.hh>
#include <opm/models/tpsa/tpsabaseproperties.hpp>

#include <cassert>
#include <set>
#include <vector>


namespace Opm {

/*!
* \brief Linearizes TPSA equations and generates system matrix and residual for linear solver
*/
template<class TypeTag>
class TpsaLinearizer
{
    using Constraints = GetPropType<TypeTag, Properties::Constraints>;  // TODO: make TPSA constraints
    using Evaluation = GetPropType<TypeTag, Properties::EvaluationTPSA>;
    using FlowModel = GetPropType<TypeTag, Properties::Model>;
    using GeomechModel = GetPropType<TypeTag, Properties::ModelTPSA>;
    using GlobalEqVector = GetPropType<TypeTag, Properties::GlobalEqVectorTPSA>;
    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using LocalResidual = GetPropType<TypeTag, Properties::LocalResidualTPSA>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Simulator = GetPropType<TypeTag, Properties::Simulator>;
    using SparseMatrixAdapter = GetPropType<TypeTag, Properties::SparseMatrixAdapterTPSA>;
    using Stencil = GetPropType<TypeTag, Properties::Stencil>;

    using MaterialState = MaterialStateTPSA<Evaluation>;

    enum { numEq = getPropValue<TypeTag, Properties::NumEqTPSA>() };

    using ADVectorBlock = Dune::FieldVector<Evaluation, numEq>;
    using MatrixBlock = typename SparseMatrixAdapter::MatrixBlock;
    using VectorBlock = Dune::FieldVector<Scalar, numEq>;

public:
    // ///
    // Public functions
    // ///
    /*!
    * \brief Constructor
    */
    TpsaLinearizer()
    {
        simulatorPtr_ = nullptr;
    }

    /*!
    * \brief Initialize the linearizer
    *
    * \param simulator Simulator object
    *
    * At this point we can assume that all objects in the simulator have been allocated. We cannot assume that they are
    * fully initialized, though.
    */
    void init(Simulator& simulator)
    {
        simulatorPtr_ = &simulator;
        eraseMatrix();
    }

    /*!
    * \brief Causes the Jacobian matrix to be recreated from scratch before the next iteration.
    *
    * This method is usually called if the sparsity pattern has changed for some reason. (e.g. by modifications of the
    * grid or changes of the auxiliary equations.)
    */
    void eraseMatrix()
    {
        jacobian_.reset();
    }

    /*!
    * \brief Finalize creation of Jacobian matrix and make ready for linear solver
    */
    void finalize()
    {
        jacobian_->finalize();
    }

    /*!
    * \brief Initializing and/or reset residual and Jacobian
    *
    * \param domain (Sub-)domain to reset
    */
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

    /*!
    * \brief Linearize the non-linear system
    */
    void linearize()
    {
        linearizeDomain();
        linearizeAuxiliaryEquations();
    }

    /*!
    * \brief Linearize the non-linear system for the spatial domain
    *
    * The global Jacobian of the residual is assembled and the residual is evaluated for the current solution. The
    * current state of affairs (esp. the previous and the current solutions) is represented by the model object.
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
                      << " caught an exception while linearizing TPSA system:" << e.what()
                      << "\n"  << std::flush;
            succeeded = 0;
        }
        catch (...) {
            std::cout << "rank " << simulator_().gridView().comm().rank()
                      << " caught an exception while linearizing TPSA system"
                      << "\n"  << std::flush;
            succeeded = 0;
        }
        succeeded = simulator_().gridView().comm().min(succeeded);

        if (!succeeded) {
            throw NumericalProblem("A process did not succeed in linearizing the TPSA system");
        }
    }

    /*!
    * \brief Linearize the non-linear system for the spatial domain
    *
    * \param domain (Sub-)domain to linearize
    *
    * The global Jacobian of the residual is assembled and the residual is evaluated for the current solution.  The
    * current state of affairs (esp. the previous and the current solutions) is represented by the model object.
    */
    template <class SubDomainType>
    void linearizeDomain(const SubDomainType& domain)
    {
        OPM_TIMEBLOCK(linearizeDomain);

        // We defer the initialization of the Jacobian matrix until here because the auxiliary modules usually assume
        // the problem, model and grid to be fully initialized'
        if (!jacobian_) {
            initFirstIteration_();
        }

        if (domain.cells.size() == flowModel_().numTotalDof()) {
            // We are on the full domain.
            resetSystem_();
        }
        else {
            resetSystem_(domain);
        }

        linearize_(domain);
    }

    /*!
    * \brief Linearize auxillary equation
    */
    void linearizeAuxiliaryEquations()
    { }

    /*!
    * \brief Extract local residuals and sub-block Jacobians from locally computed AD residual
    *
    * \param res Residual to set
    * \param bMat Block matrix from Jacobian to set
    * \param resid Computed AD residual
    */
    void setResAndJacobi(VectorBlock& res, MatrixBlock& bMat, const ADVectorBlock& resid) const
    {
        // Scalar local residual
        for (unsigned eqIdx = 0; eqIdx < numEq; ++eqIdx) {
            res[eqIdx] = resid[eqIdx].value();
        }

        // A[dofIdx][focusDofIdx][eqIdx][pvIdx] is the partial derivative of the residual function 'eqIdx' for the
        // degree of freedom 'dofIdx' with regard to the focus variable 'pvIdx' of the degree of freedom 'focusDofIdx'
        for (unsigned eqIdx = 0; eqIdx < numEq; ++eqIdx) {
            for (unsigned pvIdx = 0; pvIdx < numEq; ++pvIdx) {
                bMat[eqIdx][pvIdx] = resid[eqIdx].derivative(pvIdx);
            }
        }
    }

    /*!
    * \brief Update boundary condition information
    */
    void updateBoundaryConditionData()
    {
        for (auto& bdyInfo : boundaryInfo_) {
            // Get boundary information from problem
            const auto [type, displacementAD] = problem_().mechBoundaryCondition(bdyInfo.cell, bdyInfo.dir);

            // Strip the unnecessary (and zero anyway) derivatives off displacement
            std::vector<double> displacement(3, 0.0);
            for (std::size_t ii = 0; ii < displacement.size(); ++ii) {
                displacement[ii] = displacementAD[ii].value();
            }

            // Update boundary information
            bdyInfo.bcdata.type = type;
            bdyInfo.bcdata.displacement = displacement;
        }
    }

    // ///
    // Public get and set functions
    // ///
    /*!
    * \brief Get Jacobian matrix
    *
    * \returns Reference to (sparse) Jacobian matrix object
    */
    const SparseMatrixAdapter& jacobian() const
    { return *jacobian_; }

    /*!
    * \brief Get Jacobian matrix
    *
    * \returns Reference to (sparse) Jacobian matrix object
    */
    SparseMatrixAdapter& jacobian()
    { return *jacobian_; }

    /*!
    * \brief Get residual vector
    *
    * \returns Reference to the residual
    */
    const GlobalEqVector& residual() const
    { return residual_; }

    /*!
    * \brief Get residual vector
    *
    * \returns Reference to the residual
    */
    GlobalEqVector& residual()
    { return residual_; }

    /*!
    * \brief Get linearization type
    *
    * \returns Reference to the linearization object
    *
    * The LinearizationType controls the scheme used and the focus time index. The default is fully implicit scheme, and
    * focus index equal to 0, i.e. current time (end of step).
    */
    const LinearizationType& getLinearizationType() const
    { return linearizationType_; }

    /*!
    * \brief Get constraints map
    *
    * \returns Constraints map
    *
    * \note No constraints implemented in TPSA
    */
    std::map<unsigned, Constraints> constraintsMap() const
    { return {}; }

    /*!
    * \brief Set linearization type
    *
    * \param linearizationType Linearization object
    */
    void setLinearizationType(LinearizationType linearizationType)
    { linearizationType_ = linearizationType; }

private:
    // ///
    // Private functions
    // ///
    /*!
    * \brief Construct the Jacobian matrix and create cell-neighbor information
    *
    * Sets up cell-neigbor and boundary information structs for easy calculation in linearize_(). In addition, sparsity
    * patterns for Jacobian matrix is set up,
    */
    void createMatrix_()
    {
        OPM_TIMEBLOCK(createMatrixTPSA);

        // If the Jacobian has been initialize before, we jump out
        if (!neighborInfo_.empty()) {
            return;
        }

        // Init. the stencil
        const auto& flowModel = flowModel_();
        Stencil stencil(gridView_(), flowModel.dofMapper());

        // Build up sparsity patterns and neighboring information for Jacobian and linearization
        std::vector<std::set<unsigned>> sparsityPattern(flowModel.numTotalDof());
        unsigned numCells = flowModel.numTotalDof();
        neighborInfo_.reserve(numCells, 6 * numCells);
        std::vector<NeighborInfo> loc_nbinfo;
        for (const auto& elem : elements(gridView_())) {
            // Loop over primary dofs in the element
            stencil.update(elem);

            for (unsigned primaryDofIdx = 0; primaryDofIdx < stencil.numPrimaryDof(); ++primaryDofIdx) {
                // Build up neighboring information for curret primary dof
                const unsigned myIdx = stencil.globalSpaceIndex(primaryDofIdx);
                loc_nbinfo.resize(stencil.numDof() - 1);

                for (unsigned dofIdx = 0; dofIdx < stencil.numDof(); ++dofIdx) {
                    // NOTE: NeighborInfo could/should be expanded with cell face parameters located in problem_()
                    // needed when computing face terms in LocalResidual
                    const unsigned neighborIdx = stencil.globalSpaceIndex(dofIdx);
                    sparsityPattern[myIdx].insert(neighborIdx);
                    if (dofIdx > 0) {
                        const auto scvfIdx = dofIdx - 1;
                        const auto& scvf = stencil.interiorFace(scvfIdx);
                        const Scalar area = scvf.area();
                        loc_nbinfo[dofIdx - 1] = NeighborInfo{ neighborIdx, area, nullptr };
                    }
                }

                // Insert local neighbor info
                neighborInfo_.appendRow(loc_nbinfo.begin(), loc_nbinfo.end());

                // Boundary condition information
                unsigned bfIndex = 0;
                for (const auto& intersection : intersections(gridView_(), elem)) {
                    if (intersection.boundary()) {
                        // Get boundary face direction
                        const auto& bf = stencil.boundaryFace(bfIndex);
                        const int dir_id = bf.dirId();

                        // Skip NNCs
                        if (dir_id < 0) {
                            continue;
                        }

                        // Get boundary information from problem()
                        const auto [type, displacementAD] = problem_().mechBoundaryCondition(myIdx, dir_id);

                        // Strip the unnecessary (and zero anyway) derivatives off displacement
                        std::vector<double> displacement(3, 0.0);
                        for (std::size_t ii = 0; ii < displacement.size(); ++ii) {
                            displacement[ii] = displacementAD[ii].value();
                        }

                        // Insert boundary condition data in container
                        BoundaryConditionData bcdata { type, displacement, bfIndex, bf.area() };
                        boundaryInfo_.push_back( { myIdx, dir_id, bfIndex, bcdata } );
                        ++bfIndex;
                        continue;
                    }
                    if (!intersection.neighbor()) {
                        ++bfIndex;
                        continue;
                    }
                }
            }
        }

        // Allocate Jacobian matrix and pointers to its sub-blocks
        jacobian_ = std::make_unique<SparseMatrixAdapter>(simulator_());
        diagMatAddress_.resize(numCells);
        jacobian_->reserve(sparsityPattern);
        for (unsigned globI = 0; globI < numCells; globI++) {
            const auto& nbInfos = neighborInfo_[globI];
            diagMatAddress_[globI] = jacobian_->blockAddress(globI, globI);
            for (auto& nbInfo : nbInfos) {
                nbInfo.matBlockAddress = jacobian_->blockAddress(nbInfo.neighbor, globI);
            }
        }

        // Create full domain
        fullDomain_.cells.resize(numCells);
        std::iota(fullDomain_.cells.begin(), fullDomain_.cells.end(), 0);
    }

    /*!
    * \brief Linearize the non-linear system of equation, computing residual and its derivatives
    *
    * Calculates TPSA equation terms for the AD residual, which in turn is used to get the derivates for the Jacobian
    */
    template <class SubDomainType>
    void linearize_(const SubDomainType& domain)
    {
        OPM_TIMEBLOCK(linearizeTPSA);

        // Extract misc. variables used in linearization
        const auto& flowModel = flowModel_();
        const auto& geoMechModel = geoMechModel_();
        auto& problem = problem_();
        const unsigned int numCells = domain.cells.size();

#ifdef _OPENMP
#pragma omp parallel for
#endif
        // Loop over cells in the domain and compute local residual and jacobian
        for (unsigned ii = 0; ii < numCells; ++ii) {
            OPM_TIMEBLOCK_LOCAL(linearizationForEachCellTPSA, Subsystem::Assembly);

            const unsigned globI = domain.cells[ii];
            const auto& nbInfos = neighborInfo_[globI];
            VectorBlock res(0.0);
            MatrixBlock bMat(0.0);
            ADVectorBlock adres(0.0);
            const MaterialState& materialStateIn = geoMechModel.materialState(globI, /*timeIdx=*/0);

            // ///
            // Face term
            // ///
            {
                OPM_TIMEBLOCK_LOCAL(faceCalculationForEachCellTPSA, Subsystem::Assembly);

                // Loop over neighboring cells
                for (const auto& nbInfo : nbInfos) {
                    OPM_TIMEBLOCK_LOCAL(calculationForEachFaceTPSA, Subsystem::Assembly);

                    // Reset local residual and Jacobian
                    res = 0.0;
                    bMat = 0.0;
                    adres = 0.0;

                    // Neighbor information
                    const unsigned globJ = nbInfo.neighbor;
                    assert(globJ != globI);
                    const MaterialState& materialStateEx = geoMechModel.materialState(globJ, /*timeIdx=*/0);

                    // Compute local face term
                    LocalResidual::computeFaceTerm(adres,
                                                   materialStateIn,
                                                   materialStateEx,
                                                   problem,
                                                   globI,
                                                   globJ);
                    adres *= nbInfo.faceArea;

                    // Extract residual and sub-block Jacobian entries from computed AD residual
                    setResAndJacobi(res, bMat, adres);

                    // Insert into global residual
                    residual_[globI] += res;

                    // Insert contribution to (globI, globI) sub-block
                    // SparseAdapter syntax: jacobian_->addToBlock(globI, globI, bMat);
                    *diagMatAddress_[globI] += bMat;

                    // Insert contribution to (globJ, globI) sub-block
                    // Note: since LocalResidual::computeFaceTerm have Evaluation on globI primary variables, it is
                    // natural to insert Jacobian entries for (globJ, globI) here, since we only need to flip the signs
                    // in the calculated face terms
                    // SparseAdapter syntax: jacobian_->addToBlock(globJ, globI, bMat);
                    bMat *= -1.0;
                    *nbInfo.matBlockAddress += bMat;
                }
            }

            // ///
            // Volume term
            // //
            adres = 0.0;
            {
                OPM_TIMEBLOCK_LOCAL(computeVolumeTerm, Subsystem::Assembly);

                // Compute local volume term
                LocalResidual::computeVolumeTerm(adres,
                                                 materialStateIn,
                                                 problem,
                                                 globI);
            }
            const double volume = flowModel.dofTotalVolume(globI);
            adres *= volume;

            // Extract residual and sub-block Jacobian entries from computed AD residual
            setResAndJacobi(res, bMat, adres);

            // Insert in global residual
            residual_[globI] += res;

            // Insert contribution to (globI, globI) sub-block
            // SparseAdapter syntax: jacobian_->addToBlock(globI, globI, bMat);
            *diagMatAddress_[globI] += bMat;

            // ///
            // Source term
            // ///
            res = 0.0;
            bMat = 0.0;
            adres = 0.0;

            // Compute local source term
            LocalResidual::computeSourceTerm(adres,
                                             problem,
                                             globI,
                                             0);
            adres *= -volume;

            // Extract residual and sub-block Jacobian entries from computed AD residual
            setResAndJacobi(res, bMat, adres);

            // Insert into global residual
            residual_[globI] += res;

            // Insert contribution to (globI, globI) sub-block
            // SparseAdapter syntax: jacobian_->addToBlock(globI, globI, bMat);
            *diagMatAddress_[globI] += bMat;
        }  // globI loop

        // ///
        // Boundary term
        // ///
        for (const auto& bdyInfo : boundaryInfo_) {
            VectorBlock res(0.0);
            MatrixBlock bMat(0.0);
            ADVectorBlock adres(0.0);
            const unsigned globI = bdyInfo.cell;
            const MaterialState& materialStateIn = geoMechModel.materialState(globI, /*timeIdx=*/0);

            // Compute local boundary condition
            LocalResidual::computeBoundaryTerm(adres,
                                               materialStateIn,
                                               bdyInfo.bcdata,
                                               problem,
                                               globI);
            adres *= bdyInfo.bcdata.faceArea;

            // Extract residual and sub-block Jacobian entries from computed AD residual
            setResAndJacobi(res, bMat, adres);

            // Insert in global residual
            residual_[globI] += res;

            // Insert contribution to (globI, globI) sub-block
            // SparseAdapter syntax: jacobian_->addToBlock(globI, globI, bMat);
            *diagMatAddress_[globI] += bMat;
        }
    }

    /*!
    * \brief Create residual and Jacobian, and reset both
    */
    void initFirstIteration_()
    {
        // initialize the BCRS matrix for the Jacobian of the residual function
        createMatrix_();

        // initialize the Jacobian matrix and the vector for the residual function
        residual_.resize(flowModel_().numTotalDof());
        resetSystem_();
    }

    /*!
    * \brief Reset Jacobian matrix and residuals to zero
    */
    void resetSystem_()
    {
        // Set residual vector entries to zero
        residual_ = 0.0;

        // Set all Jacobian matrix entries to zero
        jacobian_->clear();
    }

    // ///
    // Private get functions
    // ///
    /*!
    * \brief Return simulator object
    *
    * \returns Reference to simulator
    */
    Simulator& simulator_()
    { return *simulatorPtr_; }

    /*!
    * \brief Return simulator object
    *
    * \returns Reference to simulator
    */
    const Simulator& simulator_() const
    { return *simulatorPtr_; }

    /*!
    * \brief Return problem object
    *
    * \returns Reference to problem
    */
    Problem& problem_()
    { return simulator_().problem(); }

    /*!
    * \brief Return problem object
    *
    * \returns Reference to problem
    */
    const Problem& problem_() const
    { return simulator_().problem(); }

    /*!
    * \brief Return Flow model object
    *
    * \returns Reference to Flow model
    */
    FlowModel& flowModel_()
    { return simulator_().model(); }

    /*!
    * \brief Return Flow model object
    *
    * \returns Reference to Flow model
    */
    const FlowModel& flowModel_() const
    { return simulator_().model(); }

    /*!
    * \brief Return TPSA model object
    *
    * \returns Reference to geomechanics model
    */
    GeomechModel& geoMechModel_()
    { return problem_().geoMechModel(); }

    /*!
    * \brief Return TPSA model object
    *
    * \returns Reference to geomechanics model
    */
    const GeomechModel& geoMechModel_() const
    { return problem_().geoMechModel(); }

    /*!
    * \brief Return grid view
    *
    * \returns Reference to grid view
    */
    const GridView& gridView_() const
    { return problem_().gridView(); }

    Simulator* simulatorPtr_{};
    LinearizationType linearizationType_{};

    std::vector<MatrixBlock*> diagMatAddress_{};
    std::unique_ptr<SparseMatrixAdapter> jacobian_{};
    GlobalEqVector residual_;

    //
    // Helper structs
    //
    /*!
    * \brief Neighbor information for flux terms
    */
    struct NeighborInfo
    {
        unsigned int neighbor;
        double faceArea;
        MatrixBlock* matBlockAddress;
    };
    SparseTable<NeighborInfo> neighborInfo_{};

    /*!
    * \brief Information for boundary conditions
    */
    struct BoundaryConditionData
    {
        BCMECHType type;
        std::vector<double> displacement;
        unsigned boundaryFaceIndex;
        double faceArea;
    };

    /*!
    * \brief Boundary information
    */
    struct BoundaryInfo
    {
        unsigned int cell;
        int dir;
        unsigned int bfIndex;
        BoundaryConditionData bcdata;
    };
    std::vector<BoundaryInfo> boundaryInfo_;

    /*!
    * \brief Inforamtion of the domain to linearize
    */
    struct FullDomain
    {
        std::vector<int> cells;
        std::vector<bool> interior;
    };
    FullDomain fullDomain_;
};  // class TpsaLinearizer

}  // namespace Opm

#endif
