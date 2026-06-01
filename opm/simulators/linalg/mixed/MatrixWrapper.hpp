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
//! @tparam b block size
template <class Vector>
class MixedMatrixWrapper
{
    public:

    // extract block size
    static constexpr auto block_size = Vector::block_type::dimension;

    virtual void mv(const Vector& x, Vector& y) const;
    virtual void umv(const Vector& x, Vector& y) const;
    virtual void usmv(double alpha, const Vector& x, Vector& y) const;

    //! @brief constructor
    //!
    //! @param nrows number of block rows
    //! @param nnz number of nonzero blocks
    MixedMatrixWrapper(int nrows, int nnz)
    {
        if constexpr(block_size!=3) OPM_THROW(std::invalid_argument, "MixedMatrixWrapper only supports block size == 3! \n");

        nnz_=nnz;
        M_ = bsr_alloc();
        bsr_init(M_, nrows, nnz, block_size);
    }

    //! @brief destructor
    ~MixedMatrixWrapper() {bsr_free(M_);}

    void update(double const *data);

    int *rowptr(){return M_->rowptr;}
    int *colidx(){return M_->colidx;}

    private:
    int nnz_;
    bsr_matrix  *M_;
};

template <class Vector>
void MixedMatrixWrapper<Vector>::
mv(const Vector& x, Vector& y) const
{
    // mixed-precision block spmv (y = M.x)
    if      constexpr(b==1){printf("MixedMatrixWrapper::mv does not support block size == 1!\n");getchar();}
    else if constexpr(b==2){printf("MixedMatrixWrapper::mv does not support block size == 2!\n");getchar();}
    else if constexpr(b==3) bsr_vmspmv3(M_, &x[0][0], &y[0][0]);
    else if constexpr(b==4) bsr_vmspmv4(M_, &x[0][0], &y[0][0]);
    else
    {
        printf("MixedMatrixWrapper::mv only supports block sizes < 5!\n");
        getchar();
    }
}

template <class Vector>
void MixedMatrixWrapper<Vector>::
umv(const Vector& x, Vector& y) const
{
    // mixed-precision block spmv with update (y += M.x)
    if      constexpr(b==1){printf("MixedMatrixWrapper::umv does not support block size == 1!\n");getchar();}
    else if constexpr(b==2){printf("MixedMatrixWrapper::umv does not support block size == 2!\n");getchar();}
    else if constexpr(b==3) bsr_vmspumv3(M_, &x[0][0], &y[0][0], 1.0);
    else if constexpr(b==4) bsr_vmspumv4(M_, &x[0][0], &y[0][0], 1.0);
    else
    {
        printf("MixedMatrixWrapper::umv only supports block sizes < 4!\n");
        getchar();
    }
}

template <class Vector>
void MixedMatrixWrapper<Vector>::
usmv(double alpha, const Vector& x, Vector& y) const
{
    // scaled mixed-precision block spmv with update (y += alpha *  M.x)
    //bsr_vmspumv3(M_, &x[0][0], &y[0][0], alpha);
    if      constexpr(b==1){printf("MixedMatrixWrapper::usmv does not support block size == 1!\n");getchar();}
    else if constexpr(b==2){printf("MixedMatrixWrapper::usmv does not support block size == 2!\n");getchar();}
    else if constexpr(b==3) bsr_vmspumv3(M_, &x[0][0], &y[0][0], alpha);
    else if constexpr(b==4) bsr_vmspumv4(M_, &x[0][0], &y[0][0], alpha);
    //else if constexpr(b==4){}
    else
    {
        printf("MixedMatrixWrapper::usmv only supports block sizes < 5!\n");
        getchar();
    }
}

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
        for(int i=0;i<bb;i++) M_->dbl[bb*k + i] = B[i];
    }

    // downcast to single precision
    bsr_downcast(M_);
}
} // namespace Opm
#endif // OPM_MIXED_MATRIX_HEADER_INCLUDED
