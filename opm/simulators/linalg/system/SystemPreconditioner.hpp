#pragma once
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/bvector.hh>
#include <dune/istl/operators.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/multitypeblockmatrix.hh>
#include <dune/istl/multitypeblockvector.hh>
#include <opm/simulators/linalg/FlexibleSolver.hpp>
#include <opm/simulators/linalg/PropertyTree.hpp>
#include <opm/simulators/linalg/matrixblock.hh>
#include "MultiComm.hpp"
//#include "WellMatrixMerger.hpp"
namespace Opm {
// Simple diagonal preconditioner for the system
const int numResDofs = 3;
const int numWellDofs = 4;

using ResVector = Dune::BlockVector<Dune::FieldVector<double, numResDofs>>;
using WellVector = Dune::BlockVector<Dune::FieldVector<double, numWellDofs>>;
using SystemVector = Dune::MultiTypeBlockVector<ResVector, WellVector>;

// template<class Matrix>
// auto diagvec(const Matrix& matrix) {
//     using BlockType = typename Matrix::block_type;
//     using FieldType = typename BlockType::field_type;
//     using VectorBlock = Dune::FieldVector<FieldType, BlockType::rows>;
//     using VectorType = Dune::BlockVector<VectorBlock>;
    
//     VectorType diag(matrix.N());
//     for (size_t i = 0; i < matrix.N(); ++i) {
//         auto row = matrix[i];
//         auto col = row.find(i);
//         if (col != row.end()) {
//             // Extract diagonal elements
//             for (size_t j = 0; j < BlockType::rows; ++j) {
//                 diag[i][j] = (*col)[j][j];
//             }
//         } else {
//             // If diagonal element doesn't exist, set to 1.0
//             for (size_t j = 0; j < BlockType::rows; ++j) {
//                 diag[i][j] = 1.0;
//             }
//         }
//     }
//     return diag;
// }
class SystemPreconditioner : public Dune::Preconditioner<SystemVector, SystemVector>
{
public:
    using RRMatrix = Dune::BCRSMatrix<Opm::MatrixBlock<double, numResDofs, numResDofs>>;
    using RWMatrix = Dune::BCRSMatrix<Dune::FieldMatrix<double, numResDofs, numWellDofs>>;
    using WRMatrix = Dune::BCRSMatrix<Dune::FieldMatrix<double, numWellDofs, numResDofs>>;
    using WWMatrix = Dune::BCRSMatrix<Dune::FieldMatrix<double, numWellDofs, numWellDofs>>;
    using SystemMatrix = Dune::MultiTypeBlockMatrix<
    Dune::MultiTypeBlockVector<RRMatrix, RWMatrix>,
    Dune::MultiTypeBlockVector<WRMatrix, WWMatrix>
>;

    using ResVector = Dune::BlockVector<Dune::FieldVector<double, numResDofs>>;
    using WellVector = Dune::BlockVector<Dune::FieldVector<double, numWellDofs>>;
    using ResOperator = Dune::MatrixAdapter<RRMatrix, ResVector, ResVector>;
    using ResFlexibleSolverType = Dune::FlexibleSolver<ResOperator>;
    using WellOperator = Dune::MatrixAdapter<WWMatrix, WellVector, WellVector>;
    using WellFlexibleSolverType = Dune::FlexibleSolver<WellOperator>;        
    static constexpr auto _0 = Dune::Indices::_0;
    static constexpr auto _1 = Dune::Indices::_1;
    
    

    SystemPreconditioner(const SystemMatrix& S, 
        const std::function<ResVector()> &weightCalculator, 
        int pressureIndex, 
        const Opm::PropertyTree& prm);
    virtual void apply(SystemVector& v, const SystemVector& d) override;
    virtual void pre(SystemVector& /*x*/, SystemVector& /*b*/) override {}
    virtual void post(SystemVector& /*x*/) override {}
    virtual Dune::SolverCategory::Category category() const override {
        return Dune::SolverCategory::sequential;
        //return SolverCategory::overlapping;
    }
    
private:
    const SystemMatrix& S_;
    const Opm::PropertyTree& prm_;
    std::unique_ptr<ResOperator> rop_;
    std::unique_ptr<WellOperator> wop_;
    std::unique_ptr<Dune::InverseOperator<ResVector, ResVector>> resSolver_;
    std::unique_ptr<Dune::InverseOperator<ResVector, ResVector>> resSmoother_;
    std::unique_ptr<Dune::InverseOperator<WellVector, WellVector>> wellSolver_;
    
};

class SystemPreconditionerParallel : public Dune::Preconditioner<SystemVector, SystemVector>
{
public:
    
    using WellComm = Dune::JacComm;
    using ResComm = Dune::OwnerOverlapCopyCommunication<int, int>;
    using SystemComm = Dune::MultiCommunicator<const ResComm&, const WellComm&>;
    using RRMatrix = Dune::BCRSMatrix<Opm::MatrixBlock<double, numResDofs, numResDofs>>;
    using RWMatrix = Dune::BCRSMatrix<Dune::FieldMatrix<double, numResDofs, numWellDofs>>;
    using WRMatrix = Dune::BCRSMatrix<Dune::FieldMatrix<double, numWellDofs, numResDofs>>;
    using WWMatrix = Dune::BCRSMatrix<Dune::FieldMatrix<double, numWellDofs, numWellDofs>>;
    using SystemMatrix = Dune::MultiTypeBlockMatrix<
    Dune::MultiTypeBlockVector<RRMatrix, RWMatrix>,
    Dune::MultiTypeBlockVector<WRMatrix, WWMatrix>
>;

    using ResVector = Dune::BlockVector<Dune::FieldVector<double, numResDofs>>;
    using WellVector = Dune::BlockVector<Dune::FieldVector<double, numWellDofs>>;
    using ResOperator = Dune::OverlappingSchwarzOperator<RRMatrix, ResVector, ResVector, ResComm >; //
    //using ResOperator = Dune::MatrixAdapter<RRMatrix, ResVector, ResVector>;
    using ResFlexibleSolverType = Dune::FlexibleSolver<ResOperator>;
    using WellOperator = Dune::MatrixAdapter<WWMatrix, WellVector, WellVector>;//NB not parallel for now
    using WellFlexibleSolverType = Dune::FlexibleSolver<WellOperator>;        
    static constexpr auto _0 = Dune::Indices::_0;
    static constexpr auto _1 = Dune::Indices::_1;
    
    

    SystemPreconditionerParallel(const SystemMatrix& S, 
                                 const std::function<ResVector()> &weightCalculator, 
                                 int pressureIndex, 
                                 const Opm::PropertyTree& prm,
                                 const SystemComm& syscomm);
    virtual void apply(SystemVector& v, const SystemVector& d) override;
    virtual void pre(SystemVector& /*x*/, SystemVector& /*b*/) override {}
    virtual void post(SystemVector& /*x*/) override {}
    virtual Dune::SolverCategory::Category category() const override {
        //return Dune::SolverCategory::sequential;
        return Dune::SolverCategory::overlapping;
    }
    
private:
    const SystemMatrix& S_;
    const Opm::PropertyTree& prm_;
    const SystemComm& syscomm_;
    std::unique_ptr<ResOperator> rop_;
    std::unique_ptr<WellOperator> wop_;
    std::unique_ptr<Dune::InverseOperator<ResVector, ResVector>> resSolver_;
    std::unique_ptr<Dune::InverseOperator<ResVector, ResVector>> resSmoother_;
    std::unique_ptr<Dune::InverseOperator<WellVector, WellVector>> wellSolver_;
    
};
}
