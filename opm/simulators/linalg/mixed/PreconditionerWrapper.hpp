#ifndef OPM_MIXED_PREC_HEADER_INCLUDED
#define OPM_MIXED_PREC_HEADER_INCLUDED

#include <opm/simulators/linalg/mixed/prec.h>

namespace Opm {

template <class M, class X, class Y>
class MixedPreconditioner : public Dune::PreconditionerWithUpdate<X, Y>
{
    public:

    using domain_type = X;
    static constexpr auto block_size = domain_type::block_type::dimension;

    MixedPreconditioner(const M& A)
    {
        int nrows = A.N();
        int nnz = A.nonzeroes();

        double_data_  = &A[0][0][0][0];
        prec_         = prec_alloc();
        mixed_matrix_ = bsr_alloc();
        bsr_init(mixed_matrix_, nrows, nnz, block_size);

        // copy sparsity pattern
        int *rows = mixed_matrix_->rowptr;
        int *cols = mixed_matrix_->colidx;

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

        prec_init(prec_,mixed_matrix_);
        nnz_ = nnz;

        update();
    };
    virtual void update() override;
    virtual bool hasPerfectUpdate() const override {return true;}
    virtual void pre ([[maybe_unused]] X& x, [[maybe_unused]] Y& y) override {};
    virtual void post ([[maybe_unused]] X& x) override {};
    virtual void apply ([[maybe_unused]] X& x, [[maybe_unused]] const Y& y) override;
    virtual Dune::SolverCategory::Category category() const override { return Dune::SolverCategory::sequential; };

    private:
    double const *double_data_;
    bsr_matrix *mixed_matrix_;
    prec_t *prec_;
    int nnz_;
};

template <class M, class X, class Y>
void MixedPreconditioner<M,X,Y>::
update ()
{
    // transpose each dense block to make them column-major
    int b = block_size;
    int bb=b*b;
    double B[bb];
    for(int k=0;k<nnz_;k++)
    {
        for(int i=0;i<b;i++) for(int j=0;j<b;j++) B[b*j+i] = double_data_[bb*k + b*i + j];
        for(int i=0;i<bb;i++) mixed_matrix_->dbl[bb*k + i] = B[i];
    }

    prec_ilu0_factorize(prec_, mixed_matrix_);
    prec_downcast(prec_);
}

template <class M, class X, class Y>
void MixedPreconditioner<M,X,Y>::
apply ([[maybe_unused]] X& x, [[maybe_unused]] const Y& y)
{
    x=y;
    prec_mapply3c(prec_,&x[0][0]);
}

}
#endif // OPM_MIXED_PREC_HEADER_INCLUDED
