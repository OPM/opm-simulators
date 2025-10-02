
#include "bsr.h"
#include "bslv.h"

namespace Dune
{
template <class X, class M>
class MixedSolver : public InverseOperator<X,X>
{
    public:

    MixedSolver(M A, double tol, int maxiter, bool use_dilu)
    {
        int nrows = A.N();
        int nnz   = A.nonzeroes();
        int b     = A[0][0].N();

        // create jacobian matrix object and allocate various arrays
        jacobian_ = bsr_new();
        bsr_init(jacobian_,nrows,nnz,b);

        // initialize sparsity pattern
        int *rows = jacobian_->rowptr;
        int *cols = jacobian_->colidx;

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

        // allocate and intialize solver memory
        mem_ =bslv_new();
        bslv_init(mem_, tol, maxiter, jacobian_, use_dilu);

        //pointer to nonzero blocks
        data_ = &A[0][0][0][0];
    }

    virtual void apply (X& x, X& b, InverseOperatorResult& res) override
    {
        // transpose each dense block to make them column-major
        double B[9];
        for(int k=0;k<jacobian_->nnz;k++)
        {
            for(int i=0;i<3;i++) for(int j=0;j<3;j++) B[3*j+i] = data_[9*k + 3*i + j];
            for(int i=0;i<9;i++) jacobian_->dbl[9*k + i] = B[i];
        }

        // downcast to allow mixed precision
        bsr_downcast(jacobian_);

        // solve linear system
        //int count = bslv_pbicgstab3m(mem_, jacobian_, &b[0][0], &x[0][0]);
        int count = bslv_pbicgstab3d(mem_, jacobian_, &b[0][0], &x[0][0]);

        // return convergence information
        res.converged  = (mem_->e[count] < mem_->tol);
        res.reduction  = mem_->e[count];
        res.iterations = count;
    }

    virtual void apply (X& x, X& b, double reduction, InverseOperatorResult& res) override
    {
        x=0;
        b=0;
        res.reduction = reduction;
        OPM_THROW(std::invalid_argument, "MixedSolver::apply(...) not implemented yet.");

    }

    virtual Dune::SolverCategory::Category category() const override { return Dune::SolverCategory::sequential; };

    private:
    bsr_matrix  *jacobian_;
    bslv_memory *mem_;
    double const *data_;
};

}

