/*
  Copyright 2019 SINTEF Digital, Mathematics and Cybernetics.

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef OPM_ISTLSOLVEREBOSFLEXIBLE_HEADER_INCLUDED
#define OPM_ISTLSOLVEREBOSFLEXIBLE_HEADER_INCLUDED

#include <opm/simulators/linalg/matrixblock.hh>
#include <opm/simulators/linalg/findOverlapRowsAndColumns.hpp>
#include <opm/simulators/linalg/FlexibleSolver.hpp>
#include <opm/simulators/linalg/setupPropertyTree.hpp>
#include <opm/simulators/linalg/WriteSystemMatrixHelper.hpp>

#include <opm/common/ErrorMacros.hpp>

#include <boost/property_tree/json_parser.hpp>

#include <memory>
#include <utility>

namespace Opm::Properties {

namespace TTag {
struct FlowIstlSolverFlexible {
    using InheritsFrom = std::tuple<FlowIstlSolverParams>;
};
}

} // namespace Opm::Properties


namespace Opm
{

//=====================================================================
// Implementation for ISTL-matrix based operator
//=====================================================================
/// This class solves the fully implicit black-oil system by
/// solving the reduced system (after eliminating well variables)
/// as a block-structured matrix (one block for all cell variables) for a fixed
/// number of cell variables.
///
/// The solvers and preconditioners used are run-time configurable.
template <class TypeTag>
class ISTLSolverEbosFlexible
{
    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using SparseMatrixAdapter = GetPropType<TypeTag, Properties::SparseMatrixAdapter>;
    using VectorType = GetPropType<TypeTag, Properties::GlobalEqVector>;
    using Simulator = GetPropType<TypeTag, Properties::Simulator>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using MatrixType = typename SparseMatrixAdapter::IstlMatrix;
    using WellModel = GetPropType<TypeTag, Properties::EclWellModel>;
#if HAVE_MPI
    using Communication = Dune::OwnerOverlapCopyCommunication<int, int>;
#else
    using Communication = int; // Dummy type.
#endif
    using AbstractOperatorType = Dune::AssembledLinearOperator<MatrixType, VectorType, VectorType>;
    using WellModelOpType = WellModelAsLinearOperator<WellModel, VectorType, VectorType>;
    using SolverType = Dune::FlexibleSolver<MatrixType, VectorType>;

    // for quasiImpesWeights
    using Vector = GetPropType<TypeTag, Properties::GlobalEqVector>;
    using Indices = GetPropType<TypeTag, Properties::Indices>;
    typedef typename SparseMatrixAdapter::IstlMatrix Matrix;
    typedef typename SparseMatrixAdapter::MatrixBlock MatrixBlockType;
    typedef typename Vector::block_type BlockVector;
    using Evaluation = GetPropType<TypeTag, Properties::Evaluation>;
    using ThreadManager = GetPropType<TypeTag, Properties::ThreadManager>;
    typedef typename GridView::template Codim<0>::Entity Element;
    using ElementContext = GetPropType<TypeTag, Properties::ElementContext>;

public:
    static void registerParameters()
    {
        FlowLinearSolverParameters::registerParameters<TypeTag>();
    }

    explicit ISTLSolverEbosFlexible(const Simulator& simulator)
        : simulator_(simulator)
        , ownersFirst_(EWOMS_GET_PARAM(TypeTag, bool, OwnerCellsFirst))
        , matrixAddWellContributions_(EWOMS_GET_PARAM(TypeTag, bool, MatrixAddWellContributions))
        , interiorCellNum_(detail::numMatrixRowsToUseInSolver(simulator_.vanguard().grid(), ownersFirst_))
    {
        parameters_.template init<TypeTag>();
        prm_ = setupPropertyTree<TypeTag>(parameters_);
        extractParallelGridInformationToISTL(simulator_.vanguard().grid(), parallelInformation_);
        // For some reason simulator_.model().elementMapper() is not initialized at this stage
        // Hence const auto& elemMapper = simulator_.model().elementMapper(); does not work.
        // Set it up manually
        using ElementMapper =
            Dune::MultipleCodimMultipleGeomTypeMapper<GridView>;
        ElementMapper elemMapper(simulator_.vanguard().grid().leafGridView(), Dune::mcmgElementLayout());
        detail::findOverlapAndInterior(simulator_.vanguard().grid(), elemMapper, overlapRows_, interiorRows_);
#if HAVE_MPI
        if (parallelInformation_.type() == typeid(ParallelISTLInformation)) {
            // Parallel case.
            const ParallelISTLInformation* parinfo = std::any_cast<ParallelISTLInformation>(&parallelInformation_);
            assert(parinfo);
            comm_.reset(new Communication(parinfo->communicator()));
        }
#endif
        // Print parameters to PRT/DBG logs.
        if (simulator.gridView().comm().rank() == 0) {
            std::ostringstream os;
            os << "Property tree for linear solver:\n";
            boost::property_tree::write_json(os, prm_, true);
            OpmLog::note(os.str());
        }
    }

    void eraseMatrix()
    {
    }

    void prepare(SparseMatrixAdapter& mat, VectorType& b)
    {
#if HAVE_MPI
        static bool firstcall = true;
        if (firstcall && parallelInformation_.type() == typeid(ParallelISTLInformation)) {
            // Parallel case.
            const ParallelISTLInformation* parinfo = std::any_cast<ParallelISTLInformation>(&parallelInformation_);
            assert(parinfo);
            const size_t size = mat.istlMatrix().N();
            parinfo->copyValuesTo(comm_->indexSet(), comm_->remoteIndices(), size, 1);
            firstcall = false;
        }
        if (isParallel() && matrixAddWellContributions_) {
            makeOverlapRowsInvalid(mat.istlMatrix());
        }
#endif

        matrix_ = &mat.istlMatrix(); // Store pointer for output if needed.
        std::function<VectorType()> weightsCalculator = getWeightsCalculator(mat.istlMatrix(), b);

        if (shouldCreateSolver()) {
            if (isParallel()) {
#if HAVE_MPI
                if (matrixAddWellContributions_) {
                    using ParOperatorType = Dune::OverlappingSchwarzOperator<MatrixType, VectorType, VectorType, Communication>;
                    auto op = std::make_unique<ParOperatorType>(mat.istlMatrix(), *comm_);
                    auto sol = std::make_unique<SolverType>(*op, *comm_, prm_, weightsCalculator);
                    solver_ = std::move(sol);
                    linear_operator_ = std::move(op);
                } else {
                    if (!ownersFirst_) {
                        OPM_THROW(std::runtime_error, "In parallel, the flexible solver requires "
                                  "--owner-cells-first=true when --matrix-add-well-contributions=false is used.");
                    }
                    using ParOperatorType = WellModelGhostLastMatrixAdapter<MatrixType, VectorType, VectorType, true>;
                    auto well_op = std::make_unique<WellModelOpType>(simulator_.problem().wellModel());
                    auto op = std::make_unique<ParOperatorType>(mat.istlMatrix(), *well_op, interiorCellNum_);
                    auto sol = std::make_unique<SolverType>(*op, *comm_, prm_, weightsCalculator);
                    solver_ = std::move(sol);
                    linear_operator_ = std::move(op);
                    well_operator_ = std::move(well_op);
                }
#endif // HAVE_MPI
            } else {
                if (matrixAddWellContributions_) {
                    using SeqOperatorType = Dune::MatrixAdapter<MatrixType, VectorType, VectorType>;
                    auto op = std::make_unique<SeqOperatorType>(mat.istlMatrix());
                    auto sol = std::make_unique<SolverType>(*op, prm_, weightsCalculator);
                    solver_ = std::move(sol);
                    linear_operator_ = std::move(op);
                } else {
                    using SeqOperatorType = WellModelMatrixAdapter<MatrixType, VectorType, VectorType, false>;
                    auto well_op = std::make_unique<WellModelOpType>(simulator_.problem().wellModel());
                    auto op = std::make_unique<SeqOperatorType>(mat.istlMatrix(), *well_op);
                    auto sol = std::make_unique<SolverType>(*op, prm_, weightsCalculator);
                    solver_ = std::move(sol);
                    linear_operator_ = std::move(op);
                    well_operator_ = std::move(well_op);
                }
            }
            rhs_ = b;
        } else {
            solver_->preconditioner().update();
            rhs_ = b;
        }
    }

    bool solve(VectorType& x)
    {
        solver_->apply(x, rhs_, res_);
        this->writeMatrix();
        return res_.converged;
    }

    bool isParallel() const
    {
#if HAVE_MPI
        return parallelInformation_.type() == typeid(ParallelISTLInformation);
#else
        return false;
#endif
    }

    int iterations() const
    {
        return res_.iterations;
    }

    void setResidual(VectorType& /* b */)
    {
        // rhs_ = &b; // Must be handled in prepare() instead.
    }

    void setMatrix(const SparseMatrixAdapter& /* M */)
    {
        // matrix_ = &M.istlMatrix(); // Must be handled in prepare() instead.
    }

protected:

    bool shouldCreateSolver() const
    {
        // Decide if we should recreate the solver or just do
        // a minimal preconditioner update.
        if (!solver_) {
            return true;
        }
        const int newton_iteration = this->simulator_.model().newtonMethod().numIterations();
        bool recreate_solver = false;
        if (this->parameters_.cpr_reuse_setup_ == 0) {
            // Always recreate solver.
            recreate_solver = true;
        } else if (this->parameters_.cpr_reuse_setup_ == 1) {
            // Recreate solver on the first iteration of every timestep.
            if (newton_iteration == 0) {
                recreate_solver = true;
            }
        } else if (this->parameters_.cpr_reuse_setup_ == 2) {
            // Recreate solver if the last solve used more than 10 iterations.
            if (this->iterations() > 10) {
                recreate_solver = true;
            }
        } else {
            assert(this->parameters_.cpr_reuse_setup_ == 3);
            assert(recreate_solver == false);
            // Never recreate solver.
        }
        return recreate_solver;
    }

    std::function<VectorType()> getWeightsCalculator(const MatrixType& mat, const VectorType& b) const
    {
        std::function<VectorType()> weightsCalculator;

        auto preconditionerType = prm_.get("preconditioner.type", "cpr");
        if (preconditionerType == "cpr" || preconditionerType == "cprt") {
            const bool transpose = preconditionerType == "cprt";
            const auto weightsType = prm_.get("preconditioner.weight_type", "quasiimpes");
            const auto pressureIndex = this->prm_.get("preconditioner.pressure_var_index", 1);
            if (weightsType == "quasiimpes") {
                // weighs will be created as default in the solver
                weightsCalculator = [&mat, transpose, pressureIndex]() {
                    return Opm::Amg::getQuasiImpesWeights<MatrixType, VectorType>(mat, pressureIndex, transpose);
                };
            } else if (weightsType == "trueimpes") {
                weightsCalculator = [this, &b, pressureIndex]() {
                    return this->getTrueImpesWeights(b, pressureIndex);
                                    };
            } else if (weightsType == "quasiimpesmod") {
                // weighs will be created as default in the solver
                weightsCalculator = [this, transpose, pressureIndex]() {
                                        return Amg::getQuasiImpesWeightsMod<Matrix, Vector>(this->getMatrix(), pressureIndex, transpose);
                                    };
            } else if (weightsType == "trueimpesmod") {
                weightsCalculator = [this, pressureIndex]() {
                                        return this->getTrueImpesWeightsMod(pressureIndex);
                                    };
            } else if (weightsType == "trueimpesbasic") {
                weightsCalculator = [this, pressureIndex]() {
                                        return this->getTrueImpesWeightsBasic(pressureIndex);
                                    };
    
                    
    
            } else {
                OPM_THROW(std::invalid_argument,
                          "Weights type " << weightsType << "not implemented for cpr."
                                          << " Please use quasiimpes or trueimpes.");
            }
        }
        return weightsCalculator;
    }

    /// Zero out off-diagonal blocks on rows corresponding to overlap cells
    /// Diagonal blocks on ovelap rows are set to diag(1.0).
    void makeOverlapRowsInvalid(MatrixType& matrix) const
    {
        //value to set on diagonal
        const int numEq = MatrixType::block_type::rows;
        typename MatrixType::block_type diag_block(0.0);
        for (int eq = 0; eq < numEq; ++eq)
            diag_block[eq][eq] = 1.0;

        //loop over precalculated overlap rows and columns
        for (auto row = overlapRows_.begin(); row != overlapRows_.end(); row++ )
        {
            int lcell = *row;
            // Zero out row.
            matrix[lcell] = 0.0;

            //diagonal block set to diag(1.0).
            matrix[lcell][lcell] = diag_block;
        }
    }

    VectorType getTrueImpesWeights(const VectorType& b, const int pressureVarIndex) const
    {
        VectorType weights(b.size());
        ElementContext elemCtx(simulator_);
        Opm::Amg::getTrueImpesWeights(pressureVarIndex, weights, simulator_.vanguard().gridView(),
                                      elemCtx, simulator_.model(),
                                      ThreadManager::threadId());
        return weights;
    }

    Vector getTrueImpesWeightsMod(int pressureVarIndex) const
    {
        Vector weights(rhs_->size());
        ElementContext elemCtx(simulator_);
        Amg::getTrueImpesWeightsMod(pressureVarIndex, weights, simulator_.vanguard().gridView(),
                                    elemCtx, simulator_.model(),
                                    ThreadManager::threadId());
        return weights;
    }

    Vector getTrueImpesWeightsBasic(int pressureVarIndex) const
    {
        Vector weights(rhs_->size());
        ElementContext elemCtx(simulator_);
        Amg::getTrueImpesWeightsBasic(pressureVarIndex, weights, simulator_.vanguard().gridView(),
                                      elemCtx, simulator_.model(),
                                      ThreadManager::threadId());
        return weights;
    }

    
    void writeMatrix()
    {
        const int verbosity = prm_.get<int>("verbosity");
        const bool write_matrix = verbosity > 10;
        if (write_matrix) {
            Opm::Helper::writeSystem(this->simulator_, //simulator is only used to get names
                                     *(this->matrix_),
                                     this->rhs_,
                                     comm_.get());
        }
    }


    const Simulator& simulator_;
    MatrixType* matrix_;
    std::unique_ptr<WellModelOpType> well_operator_;
    std::unique_ptr<AbstractOperatorType> linear_operator_;
    std::unique_ptr<SolverType> solver_;
    FlowLinearSolverParameters parameters_;
    boost::property_tree::ptree prm_;
    VectorType rhs_;
    Dune::InverseOperatorResult res_;
    std::any parallelInformation_;
    bool ownersFirst_;
    bool matrixAddWellContributions_;
    int interiorCellNum_;
    std::unique_ptr<Communication> comm_;
    std::vector<int> overlapRows_;
    std::vector<int> interiorRows_;
}; // end ISTLSolverEbosFlexible

} // namespace Opm

#endif // OPM_ISTLSOLVEREBOSFLEXIBLE_HEADER_INCLUDED
