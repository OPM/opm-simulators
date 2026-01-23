#pragma once

#include <opm/simulators/linalg/ISTLSolver.hpp>
#include "WellMatrixMerger.hpp"
namespace Opm
{
        namespace SystemSolver
        {
                const int numResDofs = 3;
                const int numWellDofs = 4;

                // Define matrix and vector types
                // using RRMatrix = Dune::BCRSMatrix<Dune::FieldMatrix<double, numResDofs, numResDofs>>;
                using RRMatrix = Dune::BCRSMatrix<Opm::MatrixBlock<double, numResDofs, numResDofs>>;
                // using RWtype = Dune::FieldMatrix<double, numResDofs, numWellDofs>;
                using RWMatrix = Dune::BCRSMatrix<Dune::FieldMatrix<double, numResDofs, numWellDofs>>;
                using WRMatrix = Dune::BCRSMatrix<Dune::FieldMatrix<double, numWellDofs, numResDofs>>;
                using WWMatrix = Dune::BCRSMatrix<Dune::FieldMatrix<double, numWellDofs, numWellDofs>>;
                using RVector = Dune::BlockVector<Dune::FieldVector<double, numResDofs>>;
                using WVector = Dune::BlockVector<Dune::FieldVector<double, numWellDofs>>;
                // Define system matrix and vector types
                using SystemMatrix = Dune::MultiTypeBlockMatrix<
                    Dune::MultiTypeBlockVector<RRMatrix, RWMatrix>,
                    Dune::MultiTypeBlockVector<WRMatrix, WWMatrix>>;
                using SystemVector = Dune::MultiTypeBlockVector<RVector, WVector>;
                using Comm = Dune::OwnerOverlapCopyCommunication<int, int>;
                Dune::InverseOperatorResult solveSystem(const SystemMatrix &S, SystemVector &x, const SystemVector &b,  
                        const std::function<RVector()> &weightCalculator,
                        int pressureIndex, const Opm::PropertyTree &prm, const Comm &comm);
                Dune::InverseOperatorResult solveSystem(const SystemMatrix &S, SystemVector &x, const SystemVector &b,  
                        const std::function<RVector()> &weightCalculator,
                        int pressureIndex, const Opm::PropertyTree &prm);        
               
        } // namespace SystemSolver
        template <class TypeTag>
        class ISTLSolverExperiment : public ISTLSolver<TypeTag>
        {

        protected:
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
                using ElementChunksType = ElementChunks<GridView, Dune::Partitions::All>;

                constexpr static std::size_t pressureIndex = GetPropType<TypeTag, Properties::Indices>::pressureSwitchIdx;

                enum
                {
                        enablePolymerMolarWeight = getPropValue<TypeTag, Properties::EnablePolymerMW>()
                };
                constexpr static bool isIncompatibleWithCprw = enablePolymerMolarWeight;

#if HAVE_MPI
                using CommunicationType = Dune::OwnerOverlapCopyCommunication<int, int>;
#else
                using CommunicationType = Dune::Communication<int>;
#endif
                using Parent = ISTLSolver<TypeTag>;

        public:
                ISTLSolverExperiment(const Simulator &simulator,
                                     const FlowLinearSolverParameters &parameters,
                                     bool forceSerial = false)
                    : Parent(simulator, parameters, forceSerial)
                {
                }
                explicit ISTLSolverExperiment(const Simulator &simulator)
                    : Parent(simulator)
                {
                }
                bool solve(Vector &x) override
                {
                        OPM_TIMEBLOCK(ISTLSolverExperiment_solve);
                        // Here we could add experimental features before or after calling the base class solve.
                        const auto& prm = this->prm_[this->activeSolverNum_];
                        bool solve_system = prm.get("use_system_solver", true);
                        if (solve_system)
                        {
                                
                                // get reservoir matrix
                                const int numResDofs = 3;
                                const int numWellDofs = 4;
                                // using Indices = GetPropType<TypeTag, Properties::Indices>;
                                // static constexpr int numResDofs = Indices::numEq;
                                // static constexpr int numWellDofs = numWellDofs + 1;
                                //  Define matrix and vector types
                                // using RRMatrix = Dune::BCRSMatrix<Dune::FieldMatrix<double, numResDofs, numResDofs>>;
                                using RRMatrix = Dune::BCRSMatrix<Opm::MatrixBlock<double, numResDofs, numResDofs>>;
                                // using RWtype = Dune::FieldMatrix<double, numResDofs, numWellDofs>;
                                using RWMatrix = Dune::BCRSMatrix<Dune::FieldMatrix<double, numResDofs, numWellDofs>>;
                                using WRMatrix = Dune::BCRSMatrix<Dune::FieldMatrix<double, numWellDofs, numResDofs>>;
                                using WWMatrix = Dune::BCRSMatrix<Dune::FieldMatrix<double, numWellDofs, numWellDofs>>;
                                using RVector = Dune::BlockVector<Dune::FieldVector<double, numResDofs>>;
                                using WVector = Dune::BlockVector<Dune::FieldVector<double, numWellDofs>>;
                                // Define system matrix and vector types
                                using SystemMatrix = Dune::MultiTypeBlockMatrix<
                                    Dune::MultiTypeBlockVector<RRMatrix, RWMatrix>,
                                    Dune::MultiTypeBlockVector<WRMatrix, WWMatrix>>;
                                using SystemVector = Dune::MultiTypeBlockVector<RVector, WVector>;

                                // Define helper constants for accessing MultiTypeBlockMatrix elements
                                constexpr auto _0 = Dune::Indices::_0;
                                constexpr auto _1 = Dune::Indices::_1;
                                RRMatrix A_copy = *Parent::matrix_;
                                Opm::WellMatrixMerger merger(A_copy.N());
                                WRMatrix B1;
                                RWMatrix C1;
                                WWMatrix D1;
                                WVector w_res;
                                std::vector<int> cells1 = {}; // Cells for well 1
                                std::vector<WRMatrix> b_matrices;
                                std::vector<RWMatrix> c_matrices;
                                std::vector<WWMatrix> d_matrices;
                                std::vector<std::vector<int>> wcells;
                                std::vector<WVector> residual;

                                this->simulator_.problem().wellModel().addBCDMatrix(b_matrices, c_matrices, d_matrices, wcells, residual);
                                for (size_t i = 0; i < b_matrices.size(); ++i)
                                {
                                        merger.addWell(b_matrices[i], c_matrices[i], d_matrices[i], wcells[i], static_cast<int>(i), "Well" + std::to_string(i + 1));
                                }
                                //   merger.addWell(B1, C1, D1, cells1, 0, "Well1");

                                merger.finalize();
                                


                                // Get the merged matrices
                                const auto &mergedB = merger.getMergedB();
                                const auto &mergedC = merger.getMergedC();
                                const auto &mergedD = merger.getMergedD();
                                
                                RWMatrix C_copy = mergedC;
                                WRMatrix B_copy = mergedB;
                                WWMatrix D_copy = mergedD;
                                SystemMatrix S;
                                S[_0][_0] = A_copy;
                                S[_0][_1] = C_copy;
                                S[_1][_0] = B_copy;
                                S[_1][_1] = D_copy;

                                RVector r_res = *Parent::rhs_;
                                RVector x_r(A_copy.N());
                                x_r = 0.0;
                                WVector x_w(D_copy.N());
                                x_w = 0.0;
                                w_res = x_w;
                                // set well residual
                                size_t well_dof=0;
                                for (size_t i = 0; i < residual.size(); ++i)
                                {
                                        for (size_t j = 0; j < residual[i].size(); ++j)
                                        {
                                                w_res[well_dof] += residual[i][j];
                                                well_dof++;
                                        }
                                }
                                //std::cout << "Well residual norm: " << w_res.two_norm2() << std::endl;
                                w_res = 0.0;// this should be applied to reservoir at this point
                                SystemVector x_s{x_r, x_w};
                                SystemVector r_s{r_res, w_res};
                                const auto& prm_system = prm.get_child("system_solver");
                                std::function<RVector()> weightCalculator = this->getWeightsCalculator(prm_system.get_child("preconditioner.reservoir_solver"), this->getMatrix(), pressureIndex);
                                #if HAVE_MPI
                                        using Comm = Dune::OwnerOverlapCopyCommunication<int, int>;
                                #endif
                                        
                                const Comm& comm = *(this->comm_.get());        
                                if(this->comm_->communicator().size() > 1)
                                {
                                  //std::cout << "Parallel run with system solver..." << this->comm_->communicator().rank() <<  std::endl;
                                        const auto result = SystemSolver::solveSystem(S, x_s, r_s, weightCalculator, pressureIndex, prm_system, comm);
                                        this->iterations_ = result.iterations;
                                }else
                                {
                                   const auto result = SystemSolver::solveSystem(S, x_s, r_s, weightCalculator, pressureIndex, prm_system);        // Serial run
                                   this->iterations_ = result.iterations;
                                }
                                
                                x = x_s[_0];
                                
                        }else
                        {
                                // Call the base class solve method
                                return Parent::solve(x);
                        }
                        return true;
                }

        private:
        };

}
