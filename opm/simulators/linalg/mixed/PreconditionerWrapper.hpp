#ifndef OPM_MIXED_PREC_HEADER_INCLUDED
#define OPM_MIXED_PREC_HEADER_INCLUDED

#include <opm/simulators/linalg/mixed/prec.h>

namespace Opm {

//! @brief Wraps c-implementation of mixed-precision preconditioner
//!
//! @tparam Vector the block-vector used by linear operator
//! @tparam M block-sparse matrix type
//! @tparam X block-vector input type
//! @tparam Y block vector output type
template <class M, class X, class Y>
class MixedPreconditioner : public Dune::PreconditionerWithUpdate<X, Y>
{
    public:

    using domain_type = X;
    static constexpr auto block_size = domain_type::block_type::dimension;

    //! @brief constructor
    //!
    //! @param A block-sparse matrix
    //! @param use_dilu toggle between dilu or ilu0 factorization
    MixedPreconditioner(const M& A, bool use_dilu = false)
    {
        // Access double precision matrix data
        int nrows = A.N();
        int nnz = A.nonzeroes();
        double_data_  = &A[0][0][0][0];
        prec_         = prec_alloc();

        // allocate and initialize mixed-precision matrix
        mixed_matrix_ = bsr_alloc();
        bsr_init(mixed_matrix_, nrows, nnz, block_size);

        // copy sparsity pattern from double-precision matrix
        int *rows = mixed_matrix_->rowptr;
        int *cols = mixed_matrix_->colidx;

        int irow = 0;
        int icol = 0;
        rows[0]  = 0;
        for(auto row=A.begin(); row!=A.end(); row++)
        {
            for(auto col = row->begin(); col != row->end(); ++col)
            {
                cols[icol++] = col.index();
            }
            rows[irow+1] = icol;
            irow++;
        }

        // allocate and initialize preconditioner
        prec_init(prec_,mixed_matrix_);

        // attribute initialization
        nnz_ = nnz;
        use_dilu_ = use_dilu;

        // perform matrix factorization
        update();
    };

    //! @brief destructor
    ~MixedPreconditioner()
    {
        bsr_free(mixed_matrix_);
        prec_free(prec_);
    }

    virtual void update() override;
    virtual bool hasPerfectUpdate() const override {return true;}
    virtual void pre ([[maybe_unused]] X& x, [[maybe_unused]] Y& y) override {};
    virtual void post ([[maybe_unused]] X& x) override {};
    virtual void apply ([[maybe_unused]] X& x, [[maybe_unused]] const Y& y) override;
    virtual Dune::SolverCategory::Category category() const override { return Dune::SolverCategory::sequential; };

    private:
    bool use_dilu_;
    double const *double_data_;
    bsr_matrix *mixed_matrix_;
    prec_t *prec_;
    int nnz_;


    void matvec_mul(double *y,    float const *A, double const * x);
    void matvec_mulsub(double *y, float const *A, double const * x);
    void mat_copy(double *C, double const * A);
    void mat_inv(double *invA, const double *A);
    void mat_mulsub(double *C, double const *A, double const * B);
    void mat_rmul(double *C, double const *A);
    void mat_lmul(double const *A, double *C);
};

template <class M, class X, class Y>
void MixedPreconditioner<M,X,Y>::
update ()
{
    // transpose each dense block to make them column-major
    int const N = block_size;
    int const NN=N*N;
    double B[NN];
    for(int k=0;k<nnz_;k++)
    {
        for(int i=0;i<N;i++) for(int j=0;j<N;j++) B[N*j+i] = double_data_[NN*k + N*i + j];
        for(int i=0;i<NN;i++) mixed_matrix_->dbl[NN*k + i] = B[i];
    }

    if      constexpr(N==1){printf("MixedPreconditioner::update does not support block size == 1!\n");getchar();}
    else if constexpr(N==2) use_dilu_ ? prec_dilu_factorize2(prec_, mixed_matrix_) : prec_ilu0_factorize2(prec_, mixed_matrix_);
    else if constexpr(N==3) use_dilu_ ? prec_dilu_factorize(prec_, mixed_matrix_) : prec_ilu0_factorize(prec_, mixed_matrix_);
    else if constexpr(N==4) use_dilu_ ? prec_dilu_factorize4(prec_, mixed_matrix_) : prec_ilu0_factorize4(prec_, mixed_matrix_);
    else
    {
        bsr_matrix const *A = mixed_matrix_;
        bsr_matrix *L=prec_->L;
        bsr_matrix *D=prec_->D;
        bsr_matrix *U=prec_->U;

        int const nrows = A->nrows;

        // Splitting values of A into L, D, and U, respectively
        int kU=0;
        for(int i=0;i<nrows;i++)
        {
            for (int k=A->rowptr[i];k<A->rowptr[i+1];k++)
            {
                int j=A->colidx[k];
                if(j<i)       // struct-transpose of L
                {
                    int kL = L->rowptr[j];
                    mat_copy(L->dbl + NN*kL, A->dbl + NN*k);
                    L->rowptr[j]++;
                }
                else if(j==i) // struct-copy of D
                {
                    mat_copy(D->dbl + NN*i, A->dbl + NN*k);
                }
                else if(j>i) // struct-copy of U
                {
                    mat_copy(U->dbl + NN*kU, A->dbl + NN*k);
                    kU++;
                }
            }
        }
        // reset rowptr of L
        for(int i=nrows;i>0;i--) L->rowptr[i]=L->rowptr[i-1];
        L->rowptr[0]=0;

        // Factorizing
        int idx=0;
        int next = prec_->offsets[idx][0];
        double scale[NN];
        for(int i=0;i<A->nrows;i++)
        {
            mat_inv(scale,D->dbl+i*NN);
            mat_copy(D->dbl+NN*i, scale); //store inverse instead to simplify application
            for(int k=L->rowptr[i];k<L->rowptr[i+1];k++)
            {
                //scale column i of L
                mat_rmul(L->dbl+k*NN,scale);

                //update diagonal D
                int j=L->colidx[k];
                mat_mulsub(D->dbl+j*NN,L->dbl+k*NN,U->dbl+k*NN);
            }

            if (!use_dilu_)
            while(next<U->rowptr[i+1])
            {
                int ij = prec_->offsets[idx][0];
                int ik = prec_->offsets[idx][1];
                int jk = prec_->offsets[idx][2];

                //update off-diagonals L and U
                mat_mulsub(U->dbl+jk*NN,L->dbl+ij*NN,U->dbl+ik*NN);
                mat_mulsub(L->dbl+jk*NN,L->dbl+ik*NN,U->dbl+ij*NN);

                //update marker
                next=prec_->offsets[++idx][0];
            }

            for(int k=L->rowptr[i];k<L->rowptr[i+1];k++)
            {
                //scale row i of U
                mat_lmul(scale,U->dbl+k*NN);
            }
        }
        //prec_test(); getchar();
    }

    prec_downcast(prec_);
}

template <class M, class X, class Y>
void MixedPreconditioner<M,X,Y>::
apply ([[maybe_unused]] X& x, [[maybe_unused]] const Y& y)
{
    x=y;

    int const b = block_size;
    if      constexpr(b==1){printf("MixedPreconditioner::apply does not support block size == 1!\n");getchar();}
    else if constexpr(b==2) prec_mapply2c(prec_,&x[0][0]);
    else if constexpr(b==3) prec_mapply3c(prec_,&x[0][0]);
    else if constexpr(b==4) prec_mapply4c(prec_,&x[0][0]);
    else //if constexpr(b==4)
    {
        bsr_matrix const *L  = prec_->L;
        bsr_matrix const *D  = prec_->D;
        bsr_matrix const *U  = prec_->U;

        int const N = block_size;
        int const NN = N*N;

        // Lower triangular solve assuming ones on diagonal
        for(int i=0;i<L->ncols;i++)
        {
            double *xi = &x[0][0]+N*i;
            for(int k=L->rowptr[i];k<L->rowptr[i+1];k++)
            {
                const float *A = L->flt+k*NN;
                int j=U->colidx[k]; // should be L
                double *xj = &x[0][0]+N*j;
                matvec_mulsub(xj,A,xi);
            }

            // Muliply by (inverse) diagonal block
            const float *A = D->flt+i*NN;
            matvec_mul(xi,A,xi);
        }

        // Upper triangular solve assuming ones on diagonal`
        for(int i=U->ncols;i>0;i--)
        {
            double *xi = &x[0][0]+N*(i-1);
            for(int k=U->rowptr[i]-1;k>U->rowptr[i-1]-1;k--)
            {
                const float *A = U->flt+k*NN;
                int j=U->colidx[k];
                double const *xj =&x[0][0]+N*j;
                matvec_mulsub(xi,A,xj);
            }
        }

    }
}

template <class M, class X, class Y>
void MixedPreconditioner<M,X,Y>::
matvec_mul(double *y, float const *A, double const * x)
{
    int const N = block_size;
    double z[N];
    for(int i=0;i<N;i++) z[i] = 0.0;
    for(int j=0;j<N;j++)
    {
        double xj = x[j];
        for(int i=0;i<N;i++) z[i] += A[i+N*j]*xj;
    }
    for(int i=0;i<N;i++) y[i] = z[i];
}

template <class M, class X, class Y>
void MixedPreconditioner<M,X,Y>::
matvec_mulsub(double *y, float const *A, double const * x)
{
    int const N = block_size;
    double z[N];
    for(int i=0;i<N;i++) z[i] = 0.0;
    for(int j=0;j<N;j++)
    {
        double xj = x[j];
        for(int i=0;i<N;i++) z[i] += A[i+N*j]*xj;
    }
    for(int i=0;i<N;i++) y[i] -= z[i];
}

template <class M, class X, class Y>
void MixedPreconditioner<M,X,Y>::
mat_copy(double *C, double const * A)
{
    int const N = block_size;
    int const NN =N*N;
    for(int i=0;i<NN;i++) C[i] = A[i];
}

template <class M, class X, class Y>
void MixedPreconditioner<M,X,Y>::
mat_inv(double *invA, const double *A)
{
    int const N = block_size;
    int const NN =N*N;
    double T[NN];
    mat_copy(T,A);

    for(int k=0;k<N;k++)
    {
        double scale=-1.0/T[(N+1)*k];
        for(int i=0;i<N;i++) T[i+N*k] *= i==k?0:scale; // scale column k
        for(int j=0;j<N;j++)
        {
            if (j==k) continue;
            for(int i=0;i<N;i++) T[i+N*j] += i==k?0:T[i+N*k]*T[k+N*j]; //sweep
        }
        scale=-scale;
        for(int j=0;j<N;j++) T[k+N*j] *= scale; // scale row k
        T[(N+1)*k] = scale;
    }
    mat_copy(invA,T);
}

template <class M, class X, class Y>
void MixedPreconditioner<M,X,Y>::
mat_mulsub(double *C, double const *A, double const * B)
{
    int const N = block_size;
    double z[N];
    for(int j=0;j<N;j++)
    {
        for(int k=0;k<N;k++) z[k] = 0.0;
        for(int k=0;k<N;k++)
        {
            double xk = B[k+N*j];
            for(int i=0;i<N;i++) z[i] += A[i+N*k]*xk;
        }
        for(int i=0;i<N;i++) C[i+N*j] -= z[i];
    }
}

template <class M, class X, class Y>
void MixedPreconditioner<M,X,Y>::
mat_rmul(double *C, double const *A)
{
    int const N = block_size;
    int const NN =N*N;
    double T[NN];
    for(int j=0;j<N;j++)
    {
        for(int k=0;k<N;k++) T[k+N*j] = 0.0;
        for(int k=0;k<N;k++)
        {
            double xk = A[k+N*j];
            for(int i=0;i<N;i++) T[i+N*j] += C[i+N*k]*xk;
        }
    }
    mat_copy(C,T);
}

template <class M, class X, class Y>
void MixedPreconditioner<M,X,Y>::
mat_lmul(double const *A, double *C)
{
    int const N = block_size;
    double z[N];
    for(int j=0;j<N;j++)
    {
        for(int k=0;k<N;k++) z[k] = 0.0;
        for(int k=0;k<N;k++)
        {
            double xk = C[k+N*j];
            for(int i=0;i<N;i++) z[i] += A[i+N*k]*xk;
        }
        for(int i=0;i<N;i++) C[i+N*j] = z[i];
    }
}


}
#endif // OPM_MIXED_PREC_HEADER_INCLUDED
