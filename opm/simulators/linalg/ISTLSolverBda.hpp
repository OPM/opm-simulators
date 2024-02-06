/*
  Copyright 2016 IRIS AS
  Copyright 2019, 2020 Equinor ASA
  Copyright 2020 SINTEF Digital, Mathematics and Cybernetics

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

#ifndef OPM_ISTLSOLVER_WITH_BDA_HEADER_INCLUDED
#define OPM_ISTLSOLVER_WITH_BDA_HEADER_INCLUDED

#include <opm/simulators/linalg/ISTLSolver.hpp>

#include <cstddef>
#include <memory>
#include <set>
#include <string>
#include <vector>

namespace Opm {

class Well;

template<class Matrix, class Vector, int block_size> class BdaBridge;
class WellContributions;
namespace detail {

template<class Matrix, class Vector>
struct BdaSolverInfo
{
  using WellContribFunc = std::function<void(WellContributions&)>;
  using Bridge = BdaBridge<Matrix,Vector,Matrix::block_type::rows>;

  BdaSolverInfo(const std::string& accelerator_mode,
                const int linear_solver_verbosity,
                const int maxit,
                const double tolerance,
                const int platformID,
                const int deviceID,
                const bool opencl_ilu_parallel,
                const std::string& linsolver);

  ~BdaSolverInfo();

  template<class Grid>
  void prepare(const Grid& grid,
               const Dune::CartesianIndexMapper<Grid>& cartMapper,
               const std::vector<Well>& wellsForConn,
               const std::vector<int>& cellPartition,
               const std::size_t nonzeroes,
               const bool useWellConn);

  bool apply(Vector& rhs,
             const bool useWellConn,
             WellContribFunc getContribs,
             const int rank,
             Matrix& matrix,
             Vector& x,
             Dune::InverseOperatorResult& result);

  bool gpuActive();

  int numJacobiBlocks_ = 0;

private:
  /// Create sparsity pattern for block-Jacobi matrix based on partitioning of grid.
  /// Do not initialize the values, that is done in copyMatToBlockJac()
  template<class Grid>
  void blockJacobiAdjacency(const Grid& grid,
                            const std::vector<int>& cell_part,
                            std::size_t nonzeroes);

  void copyMatToBlockJac(const Matrix& mat, Matrix& blockJac);

  std::unique_ptr<Bridge> bridge_;
  std::string accelerator_mode_;
  std::unique_ptr<Matrix> blockJacobiForGPUILU0_;
  std::vector<std::set<int>> wellConnectionsGraph_;
};

}

/// This class solves the fully implicit black-oil system by
/// solving the reduced system (after eliminating well variables)
/// as a block-structured matrix (one block for all cell variables) for a fixed
/// number of cell variables np .
template <class TypeTag>
class ISTLSolverBda : public ISTLSolver<TypeTag>
{
protected:
    using ParentType = ISTLSolver<TypeTag>;
    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using SparseMatrixAdapter = GetPropType<TypeTag, Properties::SparseMatrixAdapter>;
    using Vector = GetPropType<TypeTag, Properties::GlobalEqVector>;
    using Indices = GetPropType<TypeTag, Properties::Indices>;
    using WellModel = GetPropType<TypeTag, Properties::WellModel>;
    using Simulator = GetPropType<TypeTag, Properties::Simulator>;
    using Matrix = typename SparseMatrixAdapter::IstlMatrix;
    using ThreadManager = GetPropType<TypeTag, Properties::ThreadManager>;
    using ElementContext = GetPropType<TypeTag, Properties::ElementContext>;
    using AbstractSolverType = Dune::InverseOperator<Vector, Vector>;
    using AbstractOperatorType = Dune::AssembledLinearOperator<Matrix, Vector, Vector>;
    using AbstractPreconditionerType = Dune::PreconditionerWithUpdate<Vector, Vector>;
    using WellModelOperator = WellModelAsLinearOperator<WellModel, Vector, Vector>;
    using ElementMapper = GetPropType<TypeTag, Properties::ElementMapper>;
    constexpr static std::size_t pressureIndex = GetPropType<TypeTag, Properties::Indices>::pressureSwitchIdx;

#if HAVE_MPI
    using CommunicationType = Dune::OwnerOverlapCopyCommunication<int,int>;
#else
    using CommunicationType = Dune::CollectiveCommunication<int>;
#endif

public:
    using AssembledLinearOperatorType = Dune::AssembledLinearOperator< Matrix, Vector, Vector >;

    /// Construct a system solver.
    /// \param[in] simulator   The opm-models simulator object
    /// \param[in] parameters  Explicit parameters for solver setup, do not
    ///                        read them from command line parameters.
    ISTLSolverBda(const Simulator& simulator, const FlowLinearSolverParameters& parameters)
        : ParentType(simulator, parameters)
    {
        initializeBda();
    }

    /// Construct a system solver.
    /// \param[in] simulator   The opm-models simulator object
    explicit ISTLSolverBda(const Simulator& simulator)
        : ParentType(simulator)
    {
        initializeBda();
    }

    void initializeBda()
    {
        OPM_TIMEBLOCK(initializeBda);

        std::string accelerator_mode = EWOMS_GET_PARAM(TypeTag, std::string, AcceleratorMode);
        // Force accelerator mode to none if using MPI.
        if ((this->simulator_.vanguard().grid().comm().size() > 1) && (accelerator_mode != "none")) {
            const bool on_io_rank = (this->simulator_.gridView().comm().rank() == 0);
            if (on_io_rank) {
                OpmLog::warning("Cannot use AcceleratorMode feature with MPI, setting AcceleratorMode to 'none'.");
            }
            accelerator_mode = "none";
        }

        if (accelerator_mode == "none") {
            return;
        }

        // Initialize the BdaBridge
        const int platformID = EWOMS_GET_PARAM(TypeTag, int, OpenclPlatformId);
        const int deviceID = EWOMS_GET_PARAM(TypeTag, int, BdaDeviceId);
        const int maxit = EWOMS_GET_PARAM(TypeTag, int, LinearSolverMaxIter);
        const double tolerance = EWOMS_GET_PARAM(TypeTag, double, LinearSolverReduction);
        const bool opencl_ilu_parallel = EWOMS_GET_PARAM(TypeTag, bool, OpenclIluParallel);
        const int linear_solver_verbosity = this->parameters_[0].linear_solver_verbosity_;
        std::string linsolver = EWOMS_GET_PARAM(TypeTag, std::string, LinearSolver);
        bdaBridge_ = std::make_unique<detail::BdaSolverInfo<Matrix,Vector>>(accelerator_mode,
                                                                            linear_solver_verbosity,
                                                                            maxit,
                                                                            tolerance,
                                                                            platformID,
                                                                            deviceID,
                                                                            opencl_ilu_parallel,
                                                                            linsolver);
    }

    void prepare(const Matrix& M, Vector& b)
    {
        OPM_TIMEBLOCK(prepare);
        [[maybe_unused]] const bool firstcall = (this->matrix_ == nullptr);

        // Avoid performing the decomposition on CPU when we also do it on GPU,
        // but we do need to initialize the pointers.
        if (bdaBridge_) {
            ParentType::initPrepare(M,b);
        } else {
            ParentType::prepare(M,b);
        }

#if HAVE_OPENCL
        // update matrix entries for solvers.
        if (firstcall && bdaBridge_) {
            // model will not change the matrix object. Hence simply store a pointer
            // to the original one with a deleter that does nothing.
            // Outch! We need to be able to scale the linear system! Hence const_cast
            // setup sparsity pattern for jacobi matrix for preconditioner (only used for openclSolver)
            bdaBridge_->numJacobiBlocks_ = EWOMS_GET_PARAM(TypeTag, int, NumJacobiBlocks);
            bdaBridge_->prepare(this->simulator_.vanguard().grid(),
                               this->simulator_.vanguard().cartesianIndexMapper(),
                               this->simulator_.vanguard().schedule().getWellsatEnd(),
                               this->simulator_.vanguard().cellPartition(),
                               this->getMatrix().nonzeroes(), this->useWellConn_);
        }
#endif
    }


    void setResidual(Vector& /* b */)
    {
        // rhs_ = &b; // Must be handled in prepare() instead.
    }

    void getResidual(Vector& b) const
    {
        b = *(this->rhs_);
    }

    void setMatrix(const SparseMatrixAdapter& /* M */)
    {
        // matrix_ = &M.istlMatrix(); // Must be handled in prepare() instead.
    }

    bool solve(Vector& x)
    {
        if (!bdaBridge_) {
            return ParentType::solve(x);
        }

        OPM_TIMEBLOCK(istlSolverBdaSolve);
        this->solveCount_ += 1;
        // Write linear system if asked for.
        const int verbosity = this->prm_[this->activeSolverNum_].template get<int>("verbosity", 0);
        const bool write_matrix = verbosity > 10;
        if (write_matrix) {
            Helper::writeSystem(this->simulator_, //simulator is only used to get names
                                this->getMatrix(),
                                *(this->rhs_),
                                this->comm_.get());
        }

        // Solve system.
        Dune::InverseOperatorResult result;

        std::function<void(WellContributions&)> getContribs =
            [this](WellContributions& w)
            {
                this->simulator_.problem().wellModel().getWellContributions(w);
            };
        if (!bdaBridge_->apply(*(this->rhs_), this->useWellConn_, getContribs,
                              this->simulator_.gridView().comm().rank(),
                              const_cast<Matrix&>(this->getMatrix()),
                              x, result))
        {
            if(bdaBridge_->gpuActive()){
                // bda solve fails use istl solver setup need to be done since it is not setup in prepare
                ParentType::prepareFlexibleSolver();
            }
            assert(this->flexibleSolver_[this->activeSolverNum_].solver_);
            this->flexibleSolver_[this->activeSolverNum_].solver_->apply(x, *(this->rhs_), result);
        }

        // Check convergence, iterations etc.
        this->checkConvergence(result);

        return this->converged_;
    }

protected:
    std::unique_ptr<detail::BdaSolverInfo<Matrix, Vector>> bdaBridge_;
}; // end ISTLSolver

} // namespace Opm

#endif // OPM_ISTLSOLVER_WITH_BDA_HEADER_INCLUDED
