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

#ifndef OPM_ISTLSOLVER_EBOS_HEADER_INCLUDED
#define OPM_ISTLSOLVER_EBOS_HEADER_INCLUDED
#if COMPILE_BDA_BRIDGE

#include <opm/simulators/linalg/ISTLSolverEbos.hpp>

namespace Opm::Properties {


namespace Opm
{

template<class Matrix, class Vector, int block_size> class BdaBridge;
class WellContributions;
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
               const size_t nonzeroes,
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
                            size_t nonzeroes);

  void copyMatToBlockJac(const Matrix& mat, Matrix& blockJac);

  std::unique_ptr<Bridge> bridge_;
  std::string accelerator_mode_;
  std::unique_ptr<Matrix> blockJacobiForGPUILU0_;
  std::vector<std::set<int>> wellConnectionsGraph_;
};



/// Create sparsity pattern for block-Jacobi matrix based on partitioning of grid.
/// Do not initialize the values, that is done in copyMatToBlockJac()
}

    /// This class solves the fully implicit black-oil system by
    /// solving the reduced system (after eliminating well variables)
    /// as a block-structured matrix (one block for all cell variables) for a fixed
    /// number of cell variables np .
    template <class TypeTag>
    class ISTLSolverEbosWithGPU : public ISTLSolverEbos
    {
    protected:
        using GridView = GetPropType<TypeTag, Properties::GridView>;
        using Scalar = GetPropType<TypeTag, Properties::Scalar>;
        using SparseMatrixAdapter = GetPropType<TypeTag, Properties::SparseMatrixAdapter>;
        using Vector = GetPropType<TypeTag, Properties::GlobalEqVector>;
        using Indices = GetPropType<TypeTag, Properties::Indices>;
        using WellModel = GetPropType<TypeTag, Properties::EclWellModel>;
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
        ISTLSolverEbosGPU(const Simulator& simulator, const FlowLinearSolverParameters& parameters)
            : ISTLSolverEbos(simulator),
        {
            this->initialize();
        }

        /// Construct a system solver.
        /// \param[in] simulator   The opm-models simulator object
        explicit ISTLSolverEbos(const Simulator& simulator)
            : simulator_(simulator),
              iterations_( 0 ),
              calls_( 0 ),
              converged_(false),
              matrix_(nullptr)
        {
            parameters_.template init<TypeTag>(simulator_.vanguard().eclState().getSimulationConfig().useCPR());
            initialize();
        }

        void initialize()
        {
            OPM_TIMEBLOCK(IstlSolverEbos);
            IstlSolverEbos::initialize()
#if COMPILE_BDA_BRIDGE
            {
                std::string accelerator_mode = EWOMS_GET_PARAM(TypeTag, std::string, AcceleratorMode);
                if ((simulator_.vanguard().grid().comm().size() > 1) && (accelerator_mode != "none")) {
                    if (on_io_rank) {
                        OpmLog::warning("Cannot use GPU with MPI, GPU are disabled");
                    }
                    accelerator_mode = "none";
                }
                const int platformID = EWOMS_GET_PARAM(TypeTag, int, OpenclPlatformId);
                const int deviceID = EWOMS_GET_PARAM(TypeTag, int, BdaDeviceId);
                const int maxit = EWOMS_GET_PARAM(TypeTag, int, LinearSolverMaxIter);
                const double tolerance = EWOMS_GET_PARAM(TypeTag, double, LinearSolverReduction);
                const bool opencl_ilu_parallel = EWOMS_GET_PARAM(TypeTag, bool, OpenclIluParallel);
                const int linear_solver_verbosity = parameters_.linear_solver_verbosity_;
                std::string linsolver = EWOMS_GET_PARAM(TypeTag, std::string, LinearSolver);
                bdaBridge = std::make_unique<detail::BdaSolverInfo<Matrix,Vector>>(accelerator_mode,
                                                                                   linear_solver_verbosity,
                                                                                   maxit,
                                                                                   tolerance,
                                                                                   platformID,
                                                                                   deviceID,
                                                                                   opencl_ilu_parallel,
                                                                                   linsolver);
            }
#else
            if (EWOMS_GET_PARAM(TypeTag, std::string, AcceleratorMode) != "none") {
                OPM_THROW(std::logic_error,"Cannot use accelerated solver since CUDA, OpenCL and amgcl were not found by cmake");
            }
#endif
        }

        void prepare(const Matrix& M, Vector& b)
        {
            OPM_TIMEBLOCK(istlSolverEbosPrepare);
            IstlSolverEbos::prepare(M,b);
            const bool firstcall = (matrix_ == nullptr);
            // update matrix entries for solvers.
            if (firstcall) {
                // ebos will not change the matrix object. Hence simply store a pointer
                // to the original one with a deleter that does nothing.
                // Outch! We need to be able to scale the linear system! Hence const_cast
                // setup sparsity pattern for jacobi matrix for preconditioner (only used for openclSolver)
#if HAVE_OPENCL
                bdaBridge->numJacobiBlocks_ = EWOMS_GET_PARAM(TypeTag, int, NumJacobiBlocks);
                bdaBridge->prepare(simulator_.vanguard().grid(),
                                   simulator_.vanguard().cartesianIndexMapper(),
                                   simulator_.vanguard().schedule().getWellsatEnd(),
                                   simulator_.vanguard().cellPartition(),
                                   getMatrix().nonzeroes(), useWellConn_);
#endif
            }
        }


        void setResidual(Vector& /* b */)
        {
            // rhs_ = &b; // Must be handled in prepare() instead.
        }

        void getResidual(Vector& b) const
        {
            b = *rhs_;
        }

        void setMatrix(const SparseMatrixAdapter& /* M */)
        {
            // matrix_ = &M.istlMatrix(); // Must be handled in prepare() instead.
        }

        bool solve(Vector& x)
        {
            OPM_TIMEBLOCK(istlSolverEbosSolve);
            calls_ += 1;
            // Write linear system if asked for.
            const int verbosity = prm_.get<int>("verbosity", 0);
            const bool write_matrix = verbosity > 10;
            if (write_matrix) {
                Helper::writeSystem(simulator_, //simulator is only used to get names
                                    getMatrix(),
                                    *rhs_,
                                    comm_.get());
            }

            // Solve system.
            Dune::InverseOperatorResult result;

            std::function<void(WellContributions&)> getContribs =
                [this](WellContributions& w)
                {
                    this->simulator_.problem().wellModel().getWellContributions(w);
                };
            if (!bdaBridge->apply(*rhs_, useWellConn_, getContribs,
                                  simulator_.gridView().comm().rank(),
                                  const_cast<Matrix&>(this->getMatrix()),
                                  x, result))
            {
                OPM_TIMEBLOCK(flexibleSolverApply);
                if(bdaBridge->gpuActive()){
                    // bda solve fails use istl solver setup need to be done since it is not setup in prepeare
                    ISTLSolverEbos::prepareFlexibleSolver();
                }
                assert(flexibleSolver_.solver_);
                flexibleSolver_.solver_->apply(x, *rhs_, result);
            }

            // Check convergence, iterations etc.
            checkConvergence(result);

            return converged_;
        }
    protected:

        void prepareFlexibleSolver()
        {
            if(bdaBridge->gpuActive()){
                ISTLSolverEbos::prepareFlexibleSolver();
            }
        }
        std::unique_ptr<detail::BdaSolverInfo<Matrix, Vector>> bdaBridge;
    }; // end ISTLSolver

} // namespace Opm
#endif //COMPILE_BDA
#endif
