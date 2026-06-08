#ifndef OPM_MIXED_SOLVER_HEADER_INCLUDED
#define OPM_MIXED_SOLVER_HEADER_INCLUDED

#include <opm/simulators/flow/BlackoilModelParameters.hpp>
#include <opm/simulators/linalg/FlowLinearSolverParameters.hpp>

#include "bsr.h"
#include "bslv.h"

namespace Dune
{
template <class X, class M>
class MixedSolver : public InverseOperator<X,X>
{
    public:

    // extract block size
    static constexpr auto block_size = X::block_type::dimension;

    MixedSolver(const M &A, double tol, int maxiter, bool use_dilu)
    {
        // verify that well contributions are added to the matrix
        if (!Opm::Parameters::Get<Opm::Parameters::MatrixAddWellContributions>()) {
        OPM_THROW(std::logic_error, "Well operators are currently not supported for mixed precision. "
        "Use --matrix-add-well-contributions=true to add well contributions to the matrix instead.");}

        int nrows = A.N();
        int nnz   = A.nonzeroes();
        int b     = A[0][0].N();

        // verify that block size is 3x3 or 4x4
        if (b<3 || b>4) {OPM_THROW(std::logic_error, "Legacy mixed precision only supports 3x3 and 4x4 blocks.");}

        // create jacobian matrix object and allocate various arrays
        jacobian_ = bsr_alloc();
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
        mem_ = bslv_alloc();
        bslv_init(mem_, tol, maxiter, jacobian_, use_dilu);

        //pointer to nonzero blocks
        data_ = &A[0][0][0][0];
    }

    ~MixedSolver()
    {
        bsr_free(jacobian_);
        bslv_free(mem_);
    }

    virtual void apply (X& x, X& b, InverseOperatorResult& res) override
    {
        // transpose each dense block to make them column-major
        int const N = block_size;
        int const NN = N*N;
        double B[NN];
        for(int k=0;k<jacobian_->nnz;k++)
        {
            for(int i=0;i<N;i++) for(int j=0;j<N;j++) B[N*j+i] = data_[NN*k + N*i + j];
            for(int i=0;i<NN;i++) jacobian_->dbl[NN*k + i] = B[i];
        }

        // downcast to allow mixed precision
        bsr_downcast(jacobian_);

        // solve linear system
        int count = 0;
             if constexpr(N==3) count = bslv_pbicgstab3m(mem_, jacobian_, &b[0][0], &x[0][0]);
        else if constexpr(N==4) count = bslv_pbicgstab4m(mem_, jacobian_, &b[0][0], &x[0][0]);

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

#endif // OPM_MIXED_SOLVER_HEADER_INCLUDED

