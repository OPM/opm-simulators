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

#include <opm/simulators/linalg/mixed/MatrixWrapper.hpp>

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

    //using SingleMatrixType  = Dune::BCRSMatrix<Opm::MatrixBlock<float, block_size, block_size>>;
    using SingleMatrixType  = MatrixWrapper<Vector, block_size>;
    using MixedOperatorType = Dune::OverlappingSchwarzOperator<SingleMatrixType, Vector, Vector, Comm>;
    //using MixedOperatorType = Dune::MatrixAdapter<SingleMatrixType, Vector, Vector>;

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
        int nnz = A.nonzeroes();
        //int b = A[0][0].N();

        double_data_ = &A[0][0][0][0];
        single_matrix_ = std::make_shared<SingleMatrixType>(nrows,nnz);
        auto &B = *single_matrix_;

        // copy sparsity pattern
        int *rows = single_matrix_->rowptr();
        int *cols = single_matrix_->colidx();

        int irow = 0;
        int icol = 0;
        rows[0]  = 0;
        for(auto row=A.begin(); row!=A.end(); row++)
        {
            for(unsigned int i=0; i<row->getsize(); i++)
            {
                cols[icol++] = row->getindexptr()[i];
            }
            rows[irow+1]     = rows[irow]+row->getsize();
            irow++;
        }

        //initializemixedoperator
        double_operator_ = op;

        mixed_operator_ = std::make_shared<MixedOperatorType>(B,comm);
        //mixed_operator_ = std::make_shared<MixedOperatorType>(B);

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
        //auto &A = double_operator_->getmat();
        single_matrix_->update(double_data_);
        //single_matrix_->update(&A[0][0][0][0]);

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


    Operator *double_operator_;
    std::shared_ptr<AbstractSolverType> solver_;
    std::shared_ptr<MixedOperatorType> mixed_operator_;
    std::shared_ptr<SingleMatrixType> single_matrix_;
    double const *double_data_;

};

}

#endif // OPM_MIXED_ADAPTER_HEADER_INCLUDED
