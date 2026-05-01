#ifndef OPM_MIXED_ADAPTER_HEADER_INCLUDED
#define OPM_MIXED_ADAPTER_HEADER_INCLUDED

//#include <opm/simulators/flow/BlackoilModelParameters.hpp>
//#include <opm/simulators/linalg/FlowLinearSolverParameters.hpp>
#include <opm/simulators/linalg/PreconditionerWithUpdate.hpp>
//#include <dune/istl/solver.hh>

#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/bvector.hh>
#include <dune/istl/schwarz.hh>
#include <dune/istl/operators.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/owneroverlapcopy.hh>

#include <dune/istl/matrixindexset.hh>

/*
// Type definitions
typedef Dune::BCRSMatrix<Dune::FieldMatrix<double,1,1>> BCRSMat;
typedef Dune::BlockVector<Dune::FieldVector<double,1>> BVector;
typedef Dune::OwnerOverlapCopyCommunication<int,int> Communication;

//Assume parallel_A is the assembled BCRSMatrix and comm is the communicator
//1. Create the parallel operator
typedef Dune::OverlappingSchwarzOperator<BCRSMat,BVector,BVector,Communication> Operator;
Operator op(parallel_A,comm);

//2. Create the scalar product
typedef Dune::OverlappingSchwarzScalarProduct<BVector,Communication> ScalarProduct;
ScalarProduct sp(comm);

//3. Create the parallel preconditioner (e.g., usingSSOR)
typedef Dune::SeqSSOR<BCRSMat,BVector,BVector> Smoother;
typedef Dune::BlockPreconditioner<BVector,BVector,Communication,Smoother> ParPrec;
Smoother smoother(parallel_A,1,1.0);
ParPrec pprec(smoother,comm);
*/

namespace Dune
{
template <class Comm, class Operator, class Vector>
class MixedAdapter:public InverseOperator<Vector, Vector>
{

    public:

    using AbstractPrecondType = Dune::PreconditionerWithUpdate<Vector,Vector>;
    using AbstractScalarProductType = Dune::ScalarProduct<Vector>;

    using typename InverseOperator<Vector, Vector>::domain_type;
    static constexpr auto block_size = domain_type::block_type::dimension;

    using SingleMatrixType  = Dune::BCRSMatrix<Opm::MatrixBlock<float, block_size, block_size>>;
#if HAVE_MPI
    using MixedOperatorType = Dune::OverlappingSchwarzOperator<SingleMatrixType, Vector, Vector, Comm>;
#else // HAVE_MPI
    using MixedOperatorType = Dune::MatrixAdapter<SingleMatrixType, Vector, Vector>;
#endif // HAVE_MPI

    MixedAdapter(Operator *op,
                 std::shared_ptr<AbstractScalarProductType> sp,
                 std::shared_ptr<AbstractPrecondType> prec,
                 const double& tol,
                 const int& maxiter,
                 const int& verbosity,
                 const Comm &comm)
    {

        auto &A = op->getmat();
        int nrows = A.N();
        //int nnz = A.nonzeroes();
        //int b = A[0][0].N();

        single_matrix_ = std::make_shared<SingleMatrixType>();
        auto &B = *single_matrix_;

        //copy sparsity pattern
        int irow=0;
        int icol=0;
        MatrixIndexSet sparsity(nrows,nrows);
        for(auto row=A.begin();row!=A.end();row++)
        {
            for(unsigned int i=0;i<row->getsize();i++)
            {
                icol = row->getindexptr()[i];
                sparsity.add(irow,icol);
            }
            irow++;
        }
        sparsity.exportIdx(B);
        //B.compress();

        //initializemixedoperator
        double_operator_ = op;

#if HAVE_MPI
        mixed_operator_ = std::make_shared<MixedOperatorType>(B,comm);
#else // HAVE_MPI
        mixed_operator_ = std::make_shared<MixedOperatorType>(B);
#endif // HAVE_MPI

        //pointerstodataarrays
        //double_data_=&A[0][0][0][0];
        //single_data_=&B[0][0][0][0];

        //initialize bicgstab solver from Dune
        solver_ = std::make_shared<Dune::BiCGSTABSolver<Vector>>(
                                                              //*op,
                                                              *mixed_operator_,
                                                              *sp,
                                                              *prec,
                                                              tol, // desired residual reduction factor
                                                              maxiter, // maximum number of iterations
                                                              verbosity);
    }

    virtual void apply(Vector &x, Vector &b, InverseOperatorResult &res) override
    {
        //demote jacobian to single precision
        auto &A = double_operator_->getmat();
        auto &B = *single_matrix_;

        double const *dbl = &A[0][0][0][0];
        float *flt        = &B[0][0][0][0];
        int N = A.nonzeroes()*block_size*block_size;
        for(int k=0;k<N;k++) flt[k] = dbl[k];
/*
        //demote jacobian to single precision
        auto &A = double_operator_->getmat();
        auto &B = *single_matrix_;

        int irow=0;
        int icol=0;
        for(auto row=A.begin();row!=A.end();row++)
        {
            for(unsigned int i=0;i<row->getsize();i++)
            {
                icol = row->getindexptr()[i];
                for(int ii=0;ii<block_size;ii++)
                for(int jj=0;jj<block_size;jj++)
                B[irow][icol][ii][jj]=A[irow][icol][ii][jj];
            }
            irow++;
        }
*/
        //apply bicgstab solver from Dune
        solver_->apply(x,b,res);
    }

    virtual void apply(Vector &x, Vector &b, double reduction, InverseOperatorResult &res) override
    {
        x=0;
        b=0;
        res.reduction = reduction;
        OPM_THROW(std::invalid_argument, "MixedAdapter::apply(...) not implemented yet.");
    }

    virtual Dune::SolverCategory::Category category() const override{return Dune::SolverCategory::overlapping;};

    private:
    using AbstractSolverType = Dune::InverseOperator<Vector,Vector>;


    //Operator* operator_;
    //std::shared_ptr<AbstractPrecondType> prec_;
    //std::shared_ptr<AbstractScalarProductType> product_;
    std::shared_ptr<AbstractSolverType> solver_;
    std::shared_ptr<MixedOperatorType> mixed_operator_;
    std::shared_ptr<SingleMatrixType> single_matrix_;
    Operator *double_operator_;
    //SingleMatrixType *single_matrix_;
    double const *double_data_;
    //float *single_data_;

};

}

#endif // OPM_MIXED_ADAPTER_HEADER_INCLUDED
