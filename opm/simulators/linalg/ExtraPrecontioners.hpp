#ifndef OPM_EXTRAPRECONTIONERS_HPP
#define OPM_EXTRAPRECONTIONERS_HPP

#include <dune/common/unused.hh>
#include <dune/common/transpose.hh>
#include <dune/istl/matrixutils.hh>
#include <dune/istl/preconditioner.hh>
#include <vector>
#include <dune/common/densematrix.hh>
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

    struct ElimPivot {
        using simd_index_type = std::size_t;
        using size_type = std::size_t;
        // ElimPivot(std::vector<simd_index_type> & pivot) : pivot_(pivot)
        ElimPivot(std::vector<std::size_t>& pivot)
            : pivot_(pivot)
        {
            // typedef typename std::vector<size_type>::size_type size_type;
            for (size_type i = 0; i < pivot_.size(); ++i)
                pivot_[i] = i;
        }

        void swap(std::size_t i, simd_index_type j)
        {
            pivot_[i] = Simd::cond(Simd::Scalar<simd_index_type>(i) == j, pivot_[i], j);
        }

        template <typename T>
        void operator()(const T&, int, int)
        {
        }

        std::vector<simd_index_type>& pivot_;
    };
    // template<typename MAT>
    // template<typename Func, class Mask>
    template <class DenseMatrix, typename Func, class Mask>
    inline void luDecomposition(DenseMatrix& A, Func func, Mask& nonsingularLanes, bool throwEarly, bool doPivoting)
    {
        using std::max;
        using std::swap;

        //typedef typename FieldTraits<value_type>::real_type real_type;
        using real_type = double;
	using size_type = std::size_t;
	using simd_index_type = std::size_t;
	using field_type = double;
        // LU decomposition of A in A
        for (size_type i = 0; i < A.N(); i++) // loop over all rows
        {
            real_type pivmax = fvmeta::absreal(A[i][i]);

            if (doPivoting) {
                // compute maximum of column
                simd_index_type imax = i;
                for (size_type k = i + 1; k < A.N(); k++) {
                    auto abs = fvmeta::absreal(A[k][i]);
                    auto mask = abs > pivmax;
                    pivmax = Simd::cond(mask, abs, pivmax);
                    imax = Simd::cond(mask, simd_index_type(k), imax);
                }
                // swap rows
                for (size_type j = 0; j < A.N(); j++) {
                    // This is a swap operation where the second operand is scattered,
                    // and on top of that is also extracted from deep within a
                    // moderately complicated data structure (a DenseMatrix), where we
                    // can't assume much on the memory layout.  On intel processors,
                    // the only instruction that might help us here is vgather, but it
                    // is unclear whether that is even faster than a software
                    // implementation, and we would also need vscatter which does not
                    // exist.  So break vectorization here and do it manually.
                    for (std::size_t l = 0; l < Simd::lanes(A[i][j]); ++l)
                        swap(Simd::lane(l, A[i][j]), Simd::lane(l, A[Simd::lane(l, imax)][j]));
                }
                func.swap(i, imax); // swap the pivot or rhs
            }

            // singular ?
            nonsingularLanes = nonsingularLanes && (pivmax != real_type(0));
            if (throwEarly) {
                if (!Simd::allTrue(nonsingularLanes))
                    DUNE_THROW(FMatrixError, "matrix is singular");
            } else { // !throwEarly
                if (!Simd::anyTrue(nonsingularLanes))
                    return;
            }

            // eliminate
            for (size_type k = i + 1; k < A.N(); k++) {
                // in the simd case, A[i][i] may be close to zero in some lanes.  Pray
                // that the result is no worse than a quiet NaN.
                field_type factor = A[k][i] / A[i][i];
                A[k][i] = factor;
                for (size_type j = i + 1; j < A.N(); j++)
                    A[k][j] -= factor * A[i][j];
                func(factor, k, i);
            }
        }
    }

    template <class DenseMatrix>
    DenseMatrix invertMatrix(DenseMatrix AIN, bool doPivoting)
    {
        using size_type = std::size_t;
        using simd_index_type = std::size_t;
        // copied from dune::common::densmatrix
        using std::swap;
        using MAT = DenseMatrix;
        AutonomousValue<MAT> A(AIN);
        MAT AOUT(AIN);
        std::vector<simd_index_type> pivot(AIN.N());
        // Simd::Mask<typename FieldTraits<value_type>::real_type> nonsingularLanes(true);
        Simd::Mask<double> nonsingularLanes(true);
        Dune::Details::luDecomposition(A, Dune::Details::ElimPivot(pivot), nonsingularLanes, true, doPivoting);
        auto& L = A;
        auto& U = A;

        // initialize inverse
        AOUT = 0.0; // field_type();

        for (size_type i = 0; i < A.N(); ++i)
            AOUT[i][i] = 1;

        // L Y = I; multiple right hand sides
        for (size_type i = 0; i < A.N(); i++)
            for (size_type j = 0; j < i; j++)
                for (size_type k = 0; k < A.N(); k++)
                    AOUT[i][k] -= L[i][j] * AOUT[j][k];

        // U A^{-1} = Y
        for (size_type i = A.N(); i > 0;) {
            --i;
            for (size_type k = 0; k < A.N(); k++) {
                for (size_type j = i + 1; j < A.N(); j++)
                    AOUT[i][k] -= U[i][j] * AOUT[j][k];
                AOUT[i][k] /= U[i][i];
            }
        }

        for (size_type i = A.N(); i > 0;) {
            --i;
            for (std::size_t l = 0; l < Simd::lanes(A[0][0]); ++l) {
                std::size_t pi = Simd::lane(l, pivot[i]);
                if (i != pi)
                    for (size_type j = 0; j < A.N(); ++j)
                        swap(Simd::lane(l, AOUT[j][pi]),
			     Simd::lane(l, AOUT[j][i]));
            }
        }
        return AOUT;
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
  SeqSpai0(const M& A, int n, scalar_field_type w, bool left_precond = true)
        : _A_(A)
        , _n(n)
        , _w(w)
	, _left_precond(left_precond)  
    {
        // EASY_FUNCTION();
        CheckIfDiagonalPresent<M, l>::check(_A_);
        // we build the scaling matrix, for SPAI0, it is a diagonal matrix
        _M_.resize(_A_.N());
        matrix_block_type temp;
        constexpr int sz = matrix_block_type::rows;
        //using dune_matrix = Dune::FieldMatrix<double,sz,sz>;
        // FIXME: without considering the block size
        // Assuming the block size to be 1
        
        if constexpr (sz == 0){
             OPM_THROW(std::invalid_argument, "Now allways use invert branch ");    
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
	  // code based on f(M) = min || A M -Y ||
	  // f'(M) = A' (A M -Y)  W = (A'*A) \A' *Y(would give sparse right invers)  
	  // code based on g(M) = min || M A  -Y ||
          // g(M') = A (A' M' - Y') M' =(A*A')\ A*Y' -> M = Y*A'/(A*A')
	  // with blackoil matrices left is more stable.
            for (auto row = _A_.begin(); row != _A_.end(); ++row) {
                matrix_block_type den(0.0);
                matrix_block_type vt(0.0);
		//matrix_block_type v(0.0);
                for (auto col = (*row).begin(); col != (*row).end(); ++col) {
                    const matrix_block_type tempv = (*col);
                    const matrix_block_type tempvt = Details::transposeDenseMatrix(tempv);
                    if (col.index() == row.index()) {
		      // if SPA >0 has to be extended to an matrix matrix multiplication
		      if(_left_precond){
                        den += tempv.template rightmultiplyany<sz>(tempvt);//tempv * tempvt;
		      }else{
			den += tempvt.template rightmultiplyany<sz>(tempv);//tempvt * tempv;
		      }
		      vt = tempvt;
                    }
                }
                //NB better with LU factorization
                //matrix_block_type invden = den;
		//invden.invert();
		//LU for robustenes
		matrix_block_type invden = Dune::Details::invertMatrix(den,true);
		if(_left_precond){		  
		  _M_[row.index()] = vt.template rightmultiplyany<sz>(invden);// vt* inv(v*vt)
		  
		  //matrix_block_type Mtrans = den.solve(v);
		  //_M_[row.index()] = transposeDenseMatrix(Mtrans);// vt* den.invert()
		}else{
		  _M_[row.index()] = invden.template leftmultiplyany<sz>(vt);// inv(vt*v) *vt
		  //Best code probably
		  //_M_[row.index()] = den.solve(vt);// vt* den.invert()
		}
		// IF SPAI> 0 should be implemented more indexing
		
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
    bool _left_precond;

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
