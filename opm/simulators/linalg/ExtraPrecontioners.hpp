#ifndef OPM_EXTRAPRECONTIONERS_HPP
#define OPM_EXTRAPRECONTIONERS_HPP

#include <dune/common/unused.hh>
#include <dune/common/transpose.hh>
#include <dune/istl/matrixutils.hh>
#include <dune/istl/preconditioner.hh>
#include <vector>

namespace Dune
{
namespace Details
{
    template <class DenseMatrix>
    DenseMatrix transposeDenseMatrix(const DenseMatrix& M)
    {
        DenseMatrix tmp;
        for (int i = 0; i < M.rows; ++i)
            for (int j = 0; j < M.cols; ++j)
                tmp[j][i] = M[i][j];

        return tmp;
    }
} // namespace Details

/*! \brief The sequential jacobian preconditioner.
 * It is a reimplementation to prepare for the SPAI0 smoother

   Wraps the naked ISTL generic block Jacobi preconditioner into the
    solver framework.

   \tparam M The matrix type to operate on
   \tparam X Type of the update
   \tparam Y Type of the defect
   \tparam l The block level to invert. Default is 1
 */
template <class M, class X, class Y, int l = 1>
class SeqJacNew : public Preconditioner<X, Y>
{
public:
    //! \brief The matrix type the preconditioner is for.
    typedef M matrix_type;
    //! \brief The domain type of the preconditioner.
    typedef X domain_type;
    //! \brief The range type of the preconditioner.
    typedef Y range_type;
    //! \brief The field type of the preconditioner.
    typedef typename X::field_type field_type;
    //! \brief scalar type underlying the field_type
#if DUNE_VERSION_NEWER(DUNE_ISTL, 2, 7)
    typedef Simd::Scalar<field_type> scalar_field_type;
#else
    typedef SimdScalar<field_type> scalar_field_type;
#endif

    /*! \brief Constructor.

       Constructor gets all parameters to operate the prec.
       \param A The matrix to operate on.
       \param n The number of iterations to perform.
       \param w The relaxation factor.
     */
    SeqJacNew(const M& A, int n, scalar_field_type w)
        : _A_(A)
        , _n(n)
        , _w(w)
    {
        CheckIfDiagonalPresent<M, l>::check(_A_);
        // we build the inverseD matrix
        _invD_.resize(_A_.N());
        for (size_t i = 0; i < _A_.N(); ++i) {
            _invD_[i] = _A_[i][i];
            _invD_[i].invert();
        }
    }

    /*!
       \brief Prepare the preconditioner.

       \copydoc Preconditioner::pre(X&,Y&)
     */
    virtual void pre(X& x, Y& b)
    {
        DUNE_UNUSED_PARAMETER(x);
        DUNE_UNUSED_PARAMETER(b);
    }

    /*!
       \brief Apply the preconditioner.

       \copydoc Preconditioner::apply(X&,const Y&)
     */
    virtual void apply(X& v, const Y& d)
    {
        // we need to update the defect there
        Y dd(d);
        X vv(v.size());
        for (int i = 0; i < _n; ++i) {
            for (size_t ii = 0; ii < _invD_.size(); ++ii) {
                // vv = _invD_ * dd;
                _invD_[ii].mv(dd[ii], vv[ii]);
            }

            v.axpy(_w, vv);
            // update dd for next iteration
            // dd -= invD_* (_w * vv); or
            // dd = d - A_ * v;
            if (i < _n - 1) {
                dd = d;
                _A_.mmv(v, dd);
            }
        }
    }

    /*!
       \brief Clean up.

       \copydoc Preconditioner::post(X&)
     */
    virtual void post(X& x)
    {
        DUNE_UNUSED_PARAMETER(x);
    }

    //! Category of the preconditioner (see SolverCategory::Category)
    virtual SolverCategory::Category category() const
    {
        return SolverCategory::sequential;
    }

private:
    //! \brief The matrix we operate on.
    const M& _A_;
    //! \brief The inverse of the diagnal matrix
    typedef typename M::block_type matrix_block_type;
    std::vector<matrix_block_type> _invD_;
    //! \brief The number of steps to perform during apply.
    int _n;
    //! \brief The relaxation parameter to use.
    scalar_field_type _w;
};





template <class M, class X, class Y, int l = 1>
class SeqSpai0 : public Preconditioner<X, Y>
{
public:
    //! \brief The matrix type the preconditioner is for.
    typedef M matrix_type;
    //! \brief The domain type of the preconditioner.
    typedef X domain_type;
    //! \brief The range type of the preconditioner.
    typedef Y range_type;
    //! \brief The field type of the preconditioner.
    typedef typename X::field_type field_type;
    //! \brief The inverse of the diagnal matrix
    typedef typename M::block_type matrix_block_type;
    
    //! \brief scalar type underlying the field_type
#if DUNE_VERSION_NEWER(DUNE_ISTL, 2, 7)
    typedef Simd::Scalar<field_type> scalar_field_type;
#else
    typedef SimdScalar<field_type> scalar_field_type;
#endif
    /*! \brief Constructor.

       Constructor gets all parameters to operate the prec.
       \param A The matrix to operate on.
       \param n The number of iterations to perform.
       \param w The relaxation factor.
     */
    SeqSpai0(const M& A, int n, scalar_field_type w)
        : _A_(A)
        , _n(n)
        , _w(w)
    {
        // EASY_FUNCTION();
        CheckIfDiagonalPresent<M, l>::check(_A_);
        // we build the scaling matrix, for SPAI0, it is a diagonal matrix
        _M_.resize(_A_.N());
        matrix_block_type temp;
        constexpr int sz = matrix_block_type::rows;
        using dune_matrix = Dune::FieldMatrix<double,sz,sz>;
        // FIXME: without considering the block size
        // Assuming the block size to be 1
        
        if constexpr (sz == 1){
                  //scalar case             
            for (auto row = _A_.begin(); row != _A_.end(); ++row) {
                double den = 0.;
                double v = 0.;
                for (auto col = (*row).begin(); col != (*row).end(); ++col) {
                    const double tempv = (*col)[0][0];
                    den += tempv * tempv;
                    if (col.index() == row.index()) {
                        v = tempv;
                    }
                }
                _M_[row.index()][0][0] = v / den;
            }
        }else{
#if DUNE_VERSION_NEWER(DUNE_ISTL, 2, 7)
            // block matrix case
            for (auto row = _A_.begin(); row != _A_.end(); ++row) {
                dune_matrix den(0.0);
                dune_matrix v(0.0);
                for (auto col = (*row).begin(); col != (*row).end(); ++col) {
                    const dune_matrix tempv = (*col);
                    const dune_matrix tempvt = Details::transposeDenseMatrix(tempv);
                    den += tempv.template leftmultiplyany<sz>(tempvt);//tempvt * tempv;
                    if (col.index() == row.index()) {
                        v = tempv;
                    }
                }
                //NB better with LU factorization
                //dune_matrix invden =
                den.invert();
                _M_[row.index()] = v.template rightmultiplyany<sz>(den);// v* den.invert()
            }
#else
            OPM_THROW(std::invalid_argument, "Spai0 with blocksize>0 not suppoted for dune<=2.7 ");
#endif
        }
 

    }

    /*!
       \brief Prepare the preconditioner.

       \copydoc Preconditioner::pre(X&,Y&)
     */
    virtual void pre(X& x, Y& b)
    {
        DUNE_UNUSED_PARAMETER(x);
        DUNE_UNUSED_PARAMETER(b);
    }

    /*!
       \brief Apply the preconditioner.

       \copydoc Preconditioner::apply(X&,const Y&)
     */
    virtual void apply(X& v, const Y& d)
    {
        // EASY_FUNCTION(profiler::colors::Magenta);
        if (_n == 1) {
            this->applyOnce_(v, d);
        } else {
            this->applyMultiple_(v, d);
        }
    }

    /*!
       \brief Clean up.

       \copydoc Preconditioner::post(X&)
     */
    virtual void post(X& x)
    {
        DUNE_UNUSED_PARAMETER(x);
    }

    //! Category of the preconditioner (see SolverCategory::Category)
    virtual SolverCategory::Category category() const
    {
        return SolverCategory::sequential;
    }

private:
    //! \brief The matrix we operate on.
    const M& _A_;
    //! \brief the diagnal matrix handling the scaling
    //typedef typename M::block_type matrix_block_type;
    std::vector<matrix_block_type> _M_;
    //! \brief The number of steps to perform during apply.
    int _n;
    //! \brief The relaxation parameter to use.
    scalar_field_type _w;

    void applyOnce_(X& v, const Y& d)
    {
        // v = _M_ * d;
        for (size_t ii = 0; ii < _M_.size(); ++ii) {
            _M_[ii].mv(d[ii], v[ii]);
        }
    }

    void applyMultiple_(X& v, const Y& d)
    {
        // we need to update the defect there
        Y dd(d);
        X vv(v.size());
        for (int i = 0; i < _n; ++i) {
            // vv = _M_ * dd;
            for (size_t ii = 0; ii < _M_.size(); ++ii) {
                _M_[ii].mv(dd[ii], vv[ii]);
            }

            v.axpy(_w, vv);
            // update dd for next iteration
            // dd -= invD_* (_w * vv); or
            // dd = d - A_ * v;
            if (i < _n - 1) {
                dd = d;
                _A_.mmv(v, dd);
            }
        }
    }

};


} // namespace Dune

#endif // OPM_EXTRAPRECONTIONERS_HPP
