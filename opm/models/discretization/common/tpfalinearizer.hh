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

#include <opm/models/discretization/common/baseauxiliarymodule.hh>

#include <opm/material/common/Exceptions.hpp>
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
    }

    ~TpfaLinearizer()
    {
    }

    /*!
     * \brief Register all run-time parameters for the Jacobian linearizer.
     */
    static void registerParameters()
    { }

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
            throw NumericalIssue("A process did not succeed in linearizing the system");
    }

    void finalize()
    { jacobian_->finalize(); }

    /*!
     * \brief Linearize the part of the non-linear system of equations that is associated
     *        with the spatial domain.
     */
    void linearizeAuxiliaryEquations()
    {
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
                throw NumericalIssue("linearization of an auxiliary equation failed");
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
    }

    // Construct the BCRS matrix for the Jacobian of the residual function
    void createMatrix_()
    {
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
        ElementIterator elemIt = gridView_().template begin<0>();
        const ElementIterator elemEndIt = gridView_().template end<0>();
        for (; elemIt != elemEndIt; ++elemIt) {
            const Element& elem = *elemIt;
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
                            dirId = scvf.dirId();
                        }
                        loc_nbinfo[dofIdx - 1] = NeighborInfo{neighborIdx, trans, area, dirId};
                    }
                }
                neighborInfo_.appendRow(loc_nbinfo.begin(), loc_nbinfo.end());
            }
        }

        // add the additional neighbors and degrees of freedom caused by the auxiliary
        // equations
        size_t numAuxMod = model.numAuxiliaryModules();
        for (unsigned auxModIdx = 0; auxModIdx < numAuxMod; ++auxModIdx)
            model.auxiliaryModule(auxModIdx)->addNeighbors(sparsityPattern);

        // allocate raw matrix
        jacobian_.reset(new SparseMatrixAdapter(simulator_()));

        // create matrix structure based on sparsity pattern
        jacobian_->reserve(sparsityPattern);
    }

    // reset the global linear system of equations.
    void resetSystem_()
    {
        residual_ = 0.0;
        // zero all matrix entries
        jacobian_->clear();
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
        const bool well_local = true;
        resetSystem_();
        unsigned numCells = model_().numTotalDof();
#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (unsigned globI = 0; globI < numCells; globI++) {
            const auto& nbInfos = neighborInfo_[globI]; // this is a set but should maybe be changed
            VectorBlock res(0.0);
            MatrixBlock bMat(0.0);
            ADVectorBlock adres(0.0);
            const IntensiveQuantities* intQuantsInP = model_().cachedIntensiveQuantities(globI, /*timeIdx*/ 0);
            if (intQuantsInP == nullptr) {
                throw std::logic_error("Missing updated intensive quantities for cell " + std::to_string(globI));
            }
            const IntensiveQuantities& intQuantsIn = *intQuantsInP;

            // Flux term.
            short loc = 0;
            for (const auto& nbInfo : nbInfos) {
                unsigned globJ = nbInfo.neighbor;
                assert(globJ != globI);
                res = 0.0;
                bMat = 0.0;
                adres = 0.0;
                const IntensiveQuantities* intQuantsExP = model_().cachedIntensiveQuantities(globJ, /*timeIdx*/ 0);
                if (intQuantsExP == nullptr) {
                    throw std::logic_error("Missing updated intensive quantities for cell " + std::to_string(globJ) + " when assembling fluxes for cell " + std::to_string(globI));
                }
                const IntensiveQuantities& intQuantsEx = *intQuantsExP;
                LocalResidual::computeFlux(
                       adres, problem_(), globI, globJ, intQuantsIn, intQuantsEx,
                           nbInfo.trans, nbInfo.faceArea, nbInfo.faceDirection);
                adres *= nbInfo.faceArea;
                setResAndJacobi(res, bMat, adres);
                residual_[globI] += res;
                jacobian_->addToBlock(globI, globI, bMat);
                bMat *= -1.0;
                jacobian_->addToBlock(globJ, globI, bMat);
                ++loc;
            }

            // Accumulation term.
            double dt = simulator_().timeStepSize();
            double volume = model_().dofTotalVolume(globI);
            Scalar storefac = volume / dt;
            adres = 0.0;
            LocalResidual::computeStorage(adres, intQuantsIn);
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
            jacobian_->addToBlock(globI, globI, bMat);
            // wells sources for now (should be moved out)
            if (well_local) {
                res = 0.0;
                bMat = 0.0;
                adres = 0.0;
                LocalResidual::computeSource(adres, problem_(), globI, 0);
                adres *= -volume;
                setResAndJacobi(res, bMat, adres);
                residual_[globI] += res;
                jacobian_->addToBlock(globI, globI, bMat);
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
    };
    SparseTable<NeighborInfo> neighborInfo_;
};

} // namespace Opm

#endif // TPFA_LINEARIZER
