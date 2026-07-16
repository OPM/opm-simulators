#ifndef OPM_MIXED_MATRIX_HEADER_INCLUDED
#define OPM_MIXED_MATRIX_HEADER_INCLUDED


#include <opm/simulators/linalg/mixed/bsr.h>

namespace Opm
{

//! @brief Wraps c-implementation of mixed-precision matrix
//!
//! @note matrix is stored in single precision, but all
//!       operations are performed in double-precision
//!
//! @tparam Vector the block-vector used by linear operator
template <class Vector>
class MixedMatrixWrapper
{
    public:

    // extract block size
    static constexpr auto block_size = Vector::block_type::dimension;

    //! @brief constructor
    //!
    //! @param nrows number of block rows
    //! @param nnz number of nonzero blocks
    MixedMatrixWrapper(int nrows, int nnz)
    {
        nnz_=nnz;
        M_ = bsr_alloc();
        bsr_init(M_, nrows, nnz, block_size);
    }

    //! @brief destructor
    ~MixedMatrixWrapper() {bsr_free(M_);}

    //! @brief update matrix entries
    //!
    //! @note downcasts from double precision and transposes
    //! each non-zero block entry
    //!
    //! @param data pointer to double precision data
    void update(double const *data);

    //! @brief block-sparse matrix-vector multiplication (y = M.x)
    //!
    //! @param x input vector
    //! @param y output vector
    virtual void mv(const Vector& x, Vector& y) const;

    //! @brief block-sparse matrix-vector multiplication with
    //! update (y += M.x)
    //!
    //! @param x input vector
    //! @param y output vector
    virtual void umv(const Vector& x, Vector& y) const;

    //! @brief block-sparse matrix-vector multiplication with
    //! scaled update (y += alpha * M.x)
    //!
    //! @param alpha scaling factor
    //! @param x input vector
    //! @param y output vector
    virtual void usmv(double alpha, const Vector& x, Vector& y) const;

    //! @brief access row offset pointer
    int *rowptr(){return M_->rowptr;}

    //! @brief access column index pointer
    int *colidx(){return M_->colidx;}

    private:
    int nnz_;
    bsr_matrix  *M_;
};

//! @brief mixed-precision block-sparse matrix-vector multiplication
//! (y = M.x)
//!
//! @note hand-optimized versions are provided for block-sizes
//! 2,3, and 4. A generic implementation is provided for block-
//! sizes > 4
//!
//! @param x input vector
//! @param y output vector
template <class Vector>
void MixedMatrixWrapper<Vector>::
mv(const Vector& x, Vector& y) const
{
    int const b = block_size;
    if      constexpr(b==1){OPM_THROW(std::invalid_argument, "MixedMatrixWrapper::mv does not support block size == 1!\n");}
    else if constexpr(b==2) bsr_vmspmv2(M_, &x[0][0], &y[0][0]);
    else if constexpr(b==3) bsr_vmspmv3(M_, &x[0][0], &y[0][0]);
    else if constexpr(b==4) bsr_vmspmv4(M_, &x[0][0], &y[0][0]);
    else
    {
        int nrows = M_->nrows;
        int *rowptr=M_->rowptr;
        int *colidx=M_->colidx;
        const float *data=M_->flt;

        int bb = b*b;
        double yy[bb];
        for(int i=0;i<nrows;i++)
        {
            for(int k=0;k<bb;k++) yy[k] = 0.0;
            for(int k=rowptr[i];k<rowptr[i+1];k++)
            {
                const float *AA=data+bb*k;

                int j = colidx[k];
                double const *xx = &x[0][0]+b*j;
                for(int n=0;n<b;n++) for(int m=0;m<b;m++) yy[m+b*n] += AA[m+b*n]*xx[n];
            }

            // sum over columns
            double z[b];
            for(int m=0;m<b;m++) z[m] = yy[m];
            for(int n=1;n<b;n++) for(int m=0;m<b;m++) z[m] += yy[m+b*n];

            double *y_i = &y[0][0] + b*i;
            for(int m=0;m<b;m++) y_i[m] = z[m];
        }
    }
}

//! @brief mixed-precision block-sparse matrix-vector multiplication
//! with update (y += M.x)
//!
//! @note hand-optimized versions are provided for block-sizes
//! 2,3, and 4. A generic implementation for block-sizes > 4 is
//! NOT provided
//!
//! @param x input vector
//! @param y output vector
template <class Vector>
void MixedMatrixWrapper<Vector>::
umv(const Vector& x, Vector& y) const
{
    int const b = block_size;
    if      constexpr(b==1){OPM_THROW(std::invalid_argument, "MixedMatrixWrapper::umv does not support block size == 1!\n");}
    else if constexpr(b==2) bsr_vmspumv2(M_, &x[0][0], &y[0][0], 1.0);
    else if constexpr(b==3) bsr_vmspumv3(M_, &x[0][0], &y[0][0], 1.0);
    else if constexpr(b==4) bsr_vmspumv4(M_, &x[0][0], &y[0][0], 1.0);
    else {OPM_THROW(std::invalid_argument, "MixedMatrixWrapper::umv does not support block size == 1!\n");}
}

//! @brief mixed-precision block-sparse matrix-vector multiplication
//! with scaled update (y += alpha * M.x)
//!
//! @note hand-optimized versions are provided for block-sizes
//! 2,3, and 4. A generic implementation is provided for block-
//! sizes > 4
//!
//! @param alpha scaling factor
//! @param x input vector
//! @param y output vector
template <class Vector>
void MixedMatrixWrapper<Vector>::
usmv(double alpha, const Vector& x, Vector& y) const
{
    int const b = block_size;
    if      constexpr(b==1){OPM_THROW(std::invalid_argument, "MixedMatrixWrapper::usmv does not support block size == 1!\n");}
    else if constexpr(b==2) bsr_vmspumv2(M_, &x[0][0], &y[0][0], alpha);
    else if constexpr(b==3) bsr_vmspumv3(M_, &x[0][0], &y[0][0], alpha);
    else if constexpr(b==4) bsr_vmspumv4(M_, &x[0][0], &y[0][0], alpha);
    else
    {
        int nrows = M_->nrows;
        int *rowptr=M_->rowptr;
        int *colidx=M_->colidx;
        const float *data=M_->flt;

        int bb = b*b;
        double yy[bb];
        for(int i=0;i<nrows;i++)
        {
            for(int k=0;k<bb;k++) yy[k] = 0.0;
            for(int k=rowptr[i];k<rowptr[i+1];k++)
            {
                const float *AA=data+bb*k;

                int j = colidx[k];
                double const *xx = &x[0][0]+b*j;
                for(int n=0;n<b;n++) for(int m=0;m<b;m++) yy[m+b*n] += AA[m+b*n]*xx[n];
            }

            // sum over columns
            double z[b];
            for(int m=0;m<b;m++) z[m] = yy[m];
            for(int n=1;n<b;n++) for(int m=0;m<b;m++) z[m] += yy[m+b*n];

            double *y_i = &y[0][0] + b*i;
            for(int m=0;m<b;m++) y_i[m] += alpha*z[m];
        }

    }
}

//! @brief update matrix entries
//!
//! @note downcasts from double precision and transposes
//! each non-zero block entry
//!
//! @param data pointer to double precision data
template <class Vector>
void MixedMatrixWrapper<Vector>::
update(double const *data)
{
    // transpose each dense block to make them column-major
    int const b = block_size;
    int const bb=b*b;
    double B[bb];
    for(int k=0;k<nnz_;k++)
    {
        for(int i=0;i<b;i++) for(int j=0;j<b;j++) B[b*j+i] = data[bb*k + b*i + j];
        for(int i=0;i<bb;i++) M_->flt[bb*k + i] = B[i];
    }
}


//! @brief Adapter to take advantage of the fact that all matrix rows
//! associated with ghost cells are located at the end of the matrix
//!
//! @note The underlying mixed-matrix already ignores ghost rows.
//!
//! @param M matrix class
//! @param V vector class
//! @param C communicator class
template<class M, class V, class C>
class MixedGhostLastMatrixAdapter : public Dune::AssembledLinearOperator<M,V,V>
{
public:

    //! constructor: just store a reference to matrix and communicator
    MixedGhostLastMatrixAdapter (const M& A, const C& comm) : A_( A ), comm_(comm) {}

    // y = A * x
    virtual void apply( const V& x, V& y ) const override
    {
        A_.mv(x,y);
        comm_.project(y);
    }

    // y += \alpha * A * x
    virtual void applyscaleadd (double alpha, const V& x, V& y) const override
    {
        A_.usmv(alpha,x,y);
        comm_.project(y);
    }

    // accessor to matix object
    virtual const M& getmat() const override { return A_; }

    // solver category
    Dune::SolverCategory::Category category() const override
    {
        return Dune::SolverCategory::overlapping;
    }

private:
    const M& A_;
    const C& comm_;
};

//! @brief Adapter to combine a matrix with another linear operator while
//! taking advantage of the fact that all matrix rows associated with ghost
//! cells are located at the end of the matrix
//!
//! @note The underlying mixed-matrix already ignores ghost rows.
//!
//! @param M matrix class
//! @param V vector class
//! @param C communicator class
template<class M, class V, class C>
class WellModelMixedGhostLastMatrixAdapter : public Dune::AssembledLinearOperator<M,V,V>
{
public:
    using field_type = typename V::field_type;
    using PressureMatrix = Dune::BCRSMatrix<MatrixBlock<field_type, 1, 1>>;

    //! constructor: just store a reference to a matrix
    WellModelMixedGhostLastMatrixAdapter (const M& A,
                                     const LinearOperatorExtra<V, V>& wellOper,
                                     const C& comm
                                    )
        : A_( A ), wellOper_( wellOper ), comm_ ( comm )
    {}

    // y = A * x
    virtual void apply( const V& x, V& y ) const override
    {
        A_.mv(x,y);
        wellOper_.apply(x, y);
        comm_.project(y);
    }

    // y += \alpha * A * x
    virtual void applyscaleadd (double alpha, const V& x, V& y) const override
    {
        A_.usmv(alpha,x,y);
        wellOper_.applyscaleadd(alpha, x, y);
        comm_.project(y);
    }

    // accessor to matix object
    const M& getmat() const override { return A_; }

    void addWellPressureEquations(PressureMatrix& jacobian,
                                  const V& weights,
                                  const bool use_well_weights) const
    {
        OPM_TIMEBLOCK(addWellPressureEquations);
        wellOper_.addWellPressureEquations(jacobian, weights, use_well_weights);
    }

    void addWellPressureEquationsStruct(PressureMatrix& jacobian) const
    {
        OPM_TIMEBLOCK(addWellPressureEquationsStruct);
        wellOper_.addWellPressureEquationsStruct(jacobian);
    }

    int getNumberOfExtraEquations() const
    {
        return wellOper_.getNumberOfExtraEquations();
    }

    // solver category
    Dune::SolverCategory::Category category() const override
    {
        return Dune::SolverCategory::overlapping;
    }


protected:

    const M& A_ ;
    const C& comm_ ;
    const LinearOperatorExtra<V, V>& wellOper_;
};


} // namespace Opm
#endif // OPM_MIXED_MATRIX_HEADER_INCLUDED
