/*
  Copyright 2015 Dr. Blatt - HPC-Simulation-Software & Services
  Copyright 2015 Statoil AS

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef OPM_PARALLELOVERLAPPINGILU0_HEADER_INCLUDED
#define OPM_PARALLELOVERLAPPINGILU0_HEADER_INCLUDED

#include <opm/autodiff/GraphColoring.hpp>
#include <opm/common/Exceptions.hpp>
#include <opm/common/ErrorMacros.hpp>
#include <dune/common/version.hh>
#include <dune/istl/preconditioner.hh>
#include <dune/istl/paamg/smoother.hh>
#include <dune/istl/paamg/graph.hh>
#include <dune/istl/paamg/pinfo.hh>

#include <type_traits>
#include <numeric>
#include <limits>
#include <cstddef>
#include <string>

namespace Opm
{

//template<class M, class X, class Y, class C>
//class ParallelOverlappingILU0;
template<class Matrix, class Domain, class Range, class ParallelInfo = Dune::Amg::SequentialInformation>
class ParallelOverlappingILU0;

enum class MILU_VARIANT{
    /// \brief Do not perform modified ILU
    ILU = 0,
    /// \brief \f$U_{ii} = U_{ii} +\f$  sum(dropped entries)
    MILU_1 = 1,
    /// \brief \f$U_{ii} = U_{ii} + sign(U_{ii}) * \f$ sum(dropped entries)
    MILU_2 = 2,
    /// \brief \f$U_{ii} = U_{ii} sign(U_{ii}) * \f$ sum(|dropped entries|)
    MILU_3 = 3,
    /// \brief \f$U_{ii} = U_{ii} + (U_{ii}>0?1:0) * \f$ sum(dropped entries)
    MILU_4 = 4
};

inline MILU_VARIANT convertString2Milu(std::string milu)
{
    if( 0 == milu.compare("MILU_1") )
    {
        return MILU_VARIANT::MILU_1;
    }
    if ( 0 == milu.compare("MILU_2") )
    {
        return MILU_VARIANT::MILU_2;
    }
    if ( 0 == milu.compare("MILU_3") )
    {
        return MILU_VARIANT::MILU_3;
    }
    return MILU_VARIANT::ILU;
}

template<class F>
class ParallelOverlappingILU0Args
    : public Dune::Amg::DefaultSmootherArgs<F>
{
 public:
    ParallelOverlappingILU0Args(MILU_VARIANT milu = MILU_VARIANT::ILU )
        : milu_(milu)
    {}
    void setMilu(MILU_VARIANT milu)
    {
        milu_ = milu;
    }
    MILU_VARIANT getMilu() const
    {
        return milu_;
    }
    void setN(int n)
    {
        n_ = n;
    }
    int getN() const
    {
        return n_;
    }
 private:
    MILU_VARIANT milu_;
    int n_;
};
} // end namespace Opm

namespace Dune
{

namespace Amg
{


template<class M, class X, class Y, class C>
struct SmootherTraits<Opm::ParallelOverlappingILU0<M,X,Y,C> >
{
    using Arguments = Opm::ParallelOverlappingILU0Args<typename M::field_type>;
};

/// \brief Tells AMG how to construct the Opm::ParallelOverlappingILU0 smoother
/// \tparam Matrix The type of the Matrix.
/// \tparam Domain The type of the Vector representing the domain.
/// \tparam Range The type of the Vector representing the range.
/// \tparam ParallelInfo The type of the parallel information object
///         used, e.g. Dune::OwnerOverlapCommunication
template<class Matrix, class Domain, class Range, class ParallelInfo>
struct ConstructionTraits<Opm::ParallelOverlappingILU0<Matrix,Domain,Range,ParallelInfo> >
{
    typedef Opm::ParallelOverlappingILU0<Matrix,Domain,Range,ParallelInfo> T;
    typedef DefaultParallelConstructionArgs<T,ParallelInfo> Arguments;
    static inline Opm::ParallelOverlappingILU0<Matrix,Domain,Range,ParallelInfo>* construct(Arguments& args)
    {
        return new T(args.getMatrix(),
                     args.getComm(),
                     args.getArgs().getN(),
                     args.getArgs().relaxationFactor,
                     args.getArgs().getMilu());
    }

    static inline void deconstruct(T* bp)
    {
        delete bp;
    }

};

} // end namespace Amg



} // end namespace Dune

namespace Opm
{
    namespace detail
    {
    struct Reorderer
    {
        virtual std::size_t operator[](std::size_t i) const = 0;
        virtual ~Reorderer() {}
    };

    struct NoReorderer : public Reorderer
    {
        virtual std::size_t operator[](std::size_t i) const
        {
            return i;
        }
    };

    struct RealReorderer : public Reorderer
    {
        RealReorderer(const std::vector<std::size_t>& ordering)
            : ordering_(&ordering)
        {}
        virtual std::size_t operator[](std::size_t i) const
        {
            return (*ordering_)[i];
        }
        const std::vector<std::size_t>* ordering_;
    };

    struct IdentityFunctor
    {
        template<class T>
        T operator()(const T& t)
        {
            return t;
        }
    };

    struct OneFunctor
    {
        template<class T>
        T operator()(const T&)
        {
            return 1.0;
        }
    };
    struct SignFunctor
    {
        template<class T>
        double operator()(const T& t)
        {
            if ( t < 0.0 )
            {
                return -1;
            }
            else
            {
                return 1.0;
            }
        }
    };

    struct IsPositiveFunctor
    {
        template<class T>
        double operator()(const T& t)
        {
            if ( t < 0.0 )
            {
                return 0;
            }
            else
            {
                return 1;
            }
        }
    };
    struct AbsFunctor
    {
        template<class T>
        T operator()(const T& t)
        {
            using std::abs;
            return abs(t);
        }
    };

    template<class M, class F1=detail::IdentityFunctor, class F2=detail::OneFunctor >
    void milu0_decomposition(M& A, F1 absFunctor = F1(), F2 signFunctor = F2(),
                             std::vector<typename M::block_type>* diagonal = nullptr)
    {
        if( diagonal )
        {
            diagonal->reserve(A.N());
        }

        for ( auto irow = A.begin(), iend = A.end(); irow != iend; ++irow)
        {
            auto a_i_end = irow->end();
            auto a_ik    = irow->begin();

            std::array<typename M::field_type, M::block_type::rows> sum_dropped{};

            // Eliminate entries in lower triangular matrix
            // and store factors for L
            for ( ; a_ik.index() < irow.index(); ++a_ik )
            {
                auto k = a_ik.index();
                auto a_kk = A[k].find(k);
                // L_ik = A_kk^-1 * A_ik
                a_ik->rightmultiply(*a_kk);

                // modify the rest of the row, everything right of a_ik
                // a_i* -=a_ik * a_k*
                auto a_k_end = A[k].end();
                auto a_kj = a_kk, a_ij = a_ik;
                ++a_kj; ++a_ij;

                while ( a_kj != a_k_end)
                {
                    auto modifier = *a_kj;
                    modifier.leftmultiply(*a_ik);

                    while( a_ij != a_i_end && a_ij.index() < a_kj.index())
                    {
                        ++a_ij;
                    }

                    if ( a_ij != a_i_end && a_ij.index() == a_kj.index() )
                    {
                        // Value is not dropped
                        *a_ij -= modifier;
                        ++a_ij; ++a_kj;
                    }
                    else
                    {
                        auto entry = sum_dropped.begin();
                        for( const auto& row: modifier )
                        {
                            for( const auto& colEntry: row )
                            {
                                *entry += absFunctor(-colEntry);
                            }
                            ++entry;
                        }
                        ++a_kj;
                    }
                }
            }

            if ( a_ik.index() != irow.index() )
                OPM_THROW(std::logic_error, "Matrix is missing diagonal for row " << irow.index());

            int index = 0;
            for(const auto& entry: sum_dropped)
            {
                auto& bdiag = (*a_ik)[index][index];
                bdiag += signFunctor(bdiag) * entry;
                ++index;
            }

            if ( diagonal )
            {
                diagonal->push_back(*a_ik);
            }
            a_ik->invert();   // compute inverse of diagonal block
        }
    }

    template<class M>
    void milu0_decomposition(M& A,
                             std::vector<typename M::block_type>* diagonal)
    {
        milu0_decomposition(A, detail::IdentityFunctor(), detail::OneFunctor(),
                            diagonal);
    }

    template<class M>
    void milun_decomposition(const M& A, int n, MILU_VARIANT milu, M& ILU,
                             Reorderer& ordering, Reorderer& inverseOrdering)
    {
        using Map = std::map<std::size_t, int>;

        auto iluRow = ILU.createbegin();

        for(std::size_t i = 0, iend = A.N(); i < iend; ++i)
        {
            auto& orow = A[inverseOrdering[i]];

            Map rowPattern;
            for ( auto col = orow.begin(), cend = orow.end(); col != cend; ++col)
            {
                rowPattern[ordering[col.index()]] = 0;
            }

            for(auto ik = rowPattern.begin(); ik->first < i; ++ik)
            {
                if ( ik->second < n )
                {
                    auto& rowk = ILU[ik->first];

                    for ( auto kj = rowk.find(ik->first), endk = rowk.end();
                          kj != endk; ++kj)
                    {
                        // Assume double and block_type FieldMatrix
                        // first element is misused to store generation number
                        int generation = (*kj)[0][0];
                        if(generation < n)
                        {
                            auto ij = rowPattern.find(kj.index());
                            if ( ij == rowPattern.end() )
                            {
                                rowPattern[ordering[kj.index()]] = generation + 1;
                            }
                        }
                    }
                }
            }
            // create the row
            for(const auto entry: rowPattern)
            {
                iluRow.insert(entry.first);
            }
            ++iluRow;

            // write generation to newly created row.
            auto generationPair = rowPattern.begin();
            for ( auto col = ILU[i].begin(), cend = ILU[i].end(); col != cend;
                  ++col, ++generationPair)
            {
                assert(col.index() == generationPair->first);
                (*col)[0][0] = generationPair->second;
            }
        }

        // copy Entries from A
        for(auto iter=A.begin(), iend = A.end(); iter != iend; ++iter)
        {
            auto& newRow = ILU[ordering[iter.index()]];
            // reset stored generation
            for ( auto& col: newRow)
            {
                col = 0;
            }
            // copy row.
            for(auto col = iter->begin(), cend = iter->end(); col != cend; ++col)
            {
                newRow[ordering[col.index()]] = *col;
            }
        }
        // call decomposition on pattern        
        switch ( milu )
        {
        case MILU_VARIANT::MILU_1:
            detail::milu0_decomposition ( ILU);
            break;
        case MILU_VARIANT::MILU_2:
            detail::milu0_decomposition ( ILU, detail::IdentityFunctor(),
                                          detail::SignFunctor() );
            break;
        case MILU_VARIANT::MILU_3:
            detail::milu0_decomposition ( ILU, detail::AbsFunctor(),
                                          detail::SignFunctor() );
            break;
        case MILU_VARIANT::MILU_4:
            detail::milu0_decomposition ( ILU, detail::IdentityFunctor(),
                                          detail::IsPositiveFunctor() );
            break;
        default:
            bilu0_decomposition( ILU );
            break;
        }
    }

      //! compute ILU decomposition of A. A is overwritten by its decomposition
      template<class M, class CRS, class InvVector>
      void convertToCRS(const M& A, CRS& lower, CRS& upper, InvVector& inv )
      {
        // No need to do anything for 0 rows. Return to prevent indexing a
        // a zero sized array.
        if ( A.N() == 0 )
        {
          return;
        }

        typedef typename M :: size_type size_type;

        lower.resize( A.N() );
        upper.resize( A.N() );
        inv.resize( A.N() );

        // Count the lower and upper matrix entries.
        size_type numLower = 0;
        size_type numUpper = 0;
        const auto endi = A.end();
        for (auto i = A.begin(); i != endi; ++i) {
          const size_type iIndex = i.index();
          size_type numLowerRow = 0;
          for (auto j = (*i).begin(); j.index() < iIndex; ++j) {
              ++numLowerRow;
          }
          numLower += numLowerRow;
          numUpper += (*i).size() - numLowerRow - 1;
        }
        assert(numLower + numUpper + A.N() == A.nonzeroes());

        lower.reserveAdditional( numLower );

        // implement left looking variant with stored inverse
        size_type row = 0;
        size_type colcount = 0;
        lower.rows_[ 0 ] = colcount;
        for (auto i=A.begin(); i!=endi; ++i, ++row)
        {
          const size_type iIndex  = i.index();

          // eliminate entries left of diagonal; store L factor
          for (auto j=(*i).begin(); j.index() < iIndex; ++j )
          {
            lower.push_back( (*j), j.index() );
            ++colcount;
          }
          lower.rows_[ iIndex+1 ] = colcount;
        }

        assert(colcount == numLower);

        const auto rendi = A.beforeBegin();
        row = 0;
        colcount = 0;
        upper.rows_[ 0 ] = colcount ;

        upper.reserveAdditional( numUpper );

        // NOTE: upper and inv store entries in reverse order, reverse here
        // relative to ILU
        for (auto i=A.beforeEnd(); i!=rendi; --i, ++ row )
        {
          const size_type iIndex = i.index();

          // store in reverse row order
          // eliminate entries left of diagonal; store L factor
          for (auto j=(*i).beforeEnd(); j.index()>=iIndex; --j )
          {
            const size_type jIndex = j.index();
            if( j.index() == iIndex )
            {
              inv[ row ] = (*j);
	      break;
            }
            else if ( j.index() >= i.index() )
            {
              upper.push_back( (*j), jIndex );
              ++colcount ;
            }
          }
          upper.rows_[ row+1 ] = colcount;
        }
        assert(colcount == numUpper);
      }
    } // end namespace detail


/// \brief A two-step version of an overlapping Schwarz preconditioner using one step ILU0 as
///
/// This preconditioner differs from a ParallelRestrictedOverlappingSchwarz with
/// Dune:SeqILU0 in the follwing way:
/// During apply we make sure that the current residual is consistent (i.e.
/// each process knows the same value for each index. The we solve
/// Ly= d for y and make y consistent again. Last we solve Ux = y and
/// make sure that x is consistent.
/// In contrast for ParallelRestrictedOverlappingSchwarz we solve (LU)x = d for x
/// without forcing consistency between the two steps.
/// \tparam Matrix The type of the Matrix.
/// \tparam Domain The type of the Vector representing the domain.
/// \tparam Range The type of the Vector representing the range.
/// \tparam ParallelInfo The type of the parallel information object
///         used, e.g. Dune::OwnerOverlapCommunication
template<class Matrix, class Domain, class Range, class ParallelInfoT>
class ParallelOverlappingILU0
    : public Dune::Preconditioner<Domain,Range>
{
    typedef ParallelInfoT ParallelInfo;


public:
    //! \brief The matrix type the preconditioner is for.
    typedef typename std::remove_const<Matrix>::type matrix_type;
    //! \brief The domain type of the preconditioner.
    typedef Domain domain_type;
    //! \brief The range type of the preconditioner.
    typedef Range range_type;
    //! \brief The field type of the preconditioner.
    typedef typename Domain::field_type field_type;

    typedef typename matrix_type::block_type  block_type;
    typedef typename matrix_type::size_type   size_type;

protected:
    struct CRS
    {
      CRS() : nRows_( 0 ) {}

      size_type rows() const { return nRows_; }

      size_type nonZeros() const
      {
        assert( rows_[ rows() ] != size_type(-1) );
        return rows_[ rows() ];
      }

      void resize( const size_type nRows )
      {
          if( nRows_ != nRows )
          {
            nRows_ = nRows ;
            rows_.resize( nRows_+1, size_type(-1) );
          }
      }

      void reserveAdditional( const size_type nonZeros )
      {
          const size_type needed = values_.size() + nonZeros ;
          if( values_.capacity() < needed )
          {
              const size_type estimate = needed * 1.1;
              values_.reserve( estimate );
              cols_.reserve( estimate );
          }
      }

      void push_back( const block_type& value, const size_type index )
      {
          values_.push_back( value );
          cols_.push_back( index );
      }

      std::vector< size_type  > rows_;
      std::vector< block_type > values_;
      std::vector< size_type  > cols_;
      size_type nRows_;
    };

public:
#if DUNE_VERSION_NEWER(DUNE_ISTL, 2, 6)
    Dune::SolverCategory::Category category() const override
    {
      return std::is_same<ParallelInfoT, Dune::Amg::SequentialInformation>::value ?
              Dune::SolverCategory::sequential : Dune::SolverCategory::overlapping;
    }

#else
    // define the category
    enum {
        //! \brief The category the preconditioner is part of.
        category = std::is_same<ParallelInfoT, Dune::Amg::SequentialInformation>::value ?
            Dune::SolverCategory::sequential : Dune::SolverCategory::overlapping
    };
#endif

    /*! \brief Constructor.

      Constructor gets all parameters to operate the prec.
      \param A The matrix to operate on.
      \param n ILU fill in level (for testing). This does not work in parallel.
      \param w The relaxation factor.
      \param milu The modified ILU variant to use. 0 means traditional ILU. \see MILU_VARIANT.
      \param redblack Whether to use a red-black ordering.
      \param reorder_sphere If true, we start the reordering at a root node.
                            The vertices on each layer aound it (same distance) are
                            ordered consecutivly. If false, we preserver the order of
                            the vertices with the same color.
    */
    template<class BlockType, class Alloc>
    ParallelOverlappingILU0 (const Dune::BCRSMatrix<BlockType,Alloc>& A,
                             const int n, const field_type w,
                             MILU_VARIANT milu, bool redblack=false,
                             bool reorder_sphere=true)
        : lower_(),
          upper_(),
          inv_(),
          comm_(nullptr), w_(w),
          relaxation_( std::abs( w - 1.0 ) > 1e-15 )
    {
        // BlockMatrix is a Subclass of FieldMatrix that just adds
        // methods. Therefore this cast should be safe.
        init( reinterpret_cast<const Matrix&>(A), n, milu, redblack,
              reorder_sphere  );
    }

    /*! \brief Constructor gets all parameters to operate the prec.
      \param A The matrix to operate on.
      \param comm   communication object, e.g. Dune::OwnerOverlapCopyCommunication
      \param n ILU fill in level (for testing). This does not work in parallel.
      \param w The relaxation factor.
      \param milu The modified ILU variant to use. 0 means traditional ILU. \see MILU_VARIANT.
      \param redblack Whether to use a red-black ordering.
      \param reorder_sphere If true, we start the reordering at a root node.
                            The vertices on each layer aound it (same distance) are
                            ordered consecutivly. If false, we preserver the order of
                            the vertices with the same color.
    */
    template<class BlockType, class Alloc>
    ParallelOverlappingILU0 (const Dune::BCRSMatrix<BlockType,Alloc>& A,
                             const ParallelInfo& comm, const int n, const field_type w,
                             MILU_VARIANT milu, bool redblack=false,
                             bool reorder_sphere=true)
        : lower_(),
          upper_(),
          inv_(),
          comm_(&comm), w_(w),
          relaxation_( std::abs( w - 1.0 ) > 1e-15 )
    {
        // BlockMatrix is a Subclass of FieldMatrix that just adds
        // methods. Therefore this cast should be safe.
        init( reinterpret_cast<const Matrix&>(A), n, milu, redblack,
              reorder_sphere );
    }

    /*! \brief Constructor.

      Constructor gets all parameters to operate the prec.
      \param A The matrix to operate on.
      \param w The relaxation factor.
      \param milu The modified ILU variant to use. 0 means traditional ILU. \see MILU_VARIANT.
      \param redblack Whether to use a red-black ordering.
      \param reorder_sphere If true, we start the reordering at a root node.
                  The vertices on each layer aound it (same distance) are
                  ordered consecutivly. If false, we preserver the order of
                  the vertices with the same color.
    */
    template<class BlockType, class Alloc>
    ParallelOverlappingILU0 (const Dune::BCRSMatrix<BlockType,Alloc>& A,
                             const field_type w, MILU_VARIANT milu, bool redblack=false,
                             bool reorder_sphere=true)
        : ParallelOverlappingILU0( A, 0, w, milu, redblack, reorder_sphere )
    {
    }

    /*! \brief Constructor.

      Constructor gets all parameters to operate the prec.
      \param A      The matrix to operate on.
      \param comm   communication object, e.g. Dune::OwnerOverlapCopyCommunication
      \param w      The relaxation factor.
      \param milu   The modified ILU variant to use. 0 means traditional ILU. \see MILU_VARIANT.
      \param redblack Whether to use a red-black ordering.
      \param reorder_sphere If true, we start the reordering at a root node.
                            The vertices on each layer aound it (same distance) are
                            ordered consecutivly. If false, we preserver the order of
                            the vertices with the same color.
    */
    template<class BlockType, class Alloc>
    ParallelOverlappingILU0 (const Dune::BCRSMatrix<BlockType,Alloc>& A,
                             const ParallelInfo& comm, const field_type w,
                             MILU_VARIANT milu, bool redblack=false,
                             bool reorder_sphere=true)
        : lower_(),
          upper_(),
          inv_(),
          comm_(&comm), w_(w),
          relaxation_( std::abs( w - 1.0 ) > 1e-15 )
    {
        // BlockMatrix is a Subclass of FieldMatrix that just adds
        // methods. Therefore this cast should be safe.
        init( reinterpret_cast<const Matrix&>(A), 0, milu, redblack,
              reorder_sphere );
    }

    /*!
      \brief Prepare the preconditioner.

      \copydoc Preconditioner::pre(X&,Y&)
    */
    virtual void pre (Domain& x, Range& b)
    {
        DUNE_UNUSED_PARAMETER(x);
        DUNE_UNUSED_PARAMETER(b);
    }

    /*!
      \brief Apply the preconditoner.

      \copydoc Preconditioner::apply(X&,const Y&)
    */
    virtual void apply (Domain& v, const Range& d)
    {
        Range& md = reorderD(d);
        Domain& mv = reorderV(v);
        copyOwnerToAll( md );

        // iterator types
        typedef typename Range ::block_type  dblock;
        typedef typename Domain::block_type  vblock;

        const size_type iEnd = lower_.rows();
        const size_type lastRow = iEnd - 1;
        if( iEnd != upper_.rows() )
        {
            OPM_THROW(std::logic_error,"ILU: number of lower and upper rows must be the same");
        }

        // lower triangular solve
        for( size_type i=0; i<iEnd; ++ i )
        {
          dblock rhs( md[ i ] );
          const size_type rowI     = lower_.rows_[ i ];
          const size_type rowINext = lower_.rows_[ i+1 ];

          for( size_type col = rowI; col < rowINext; ++ col )
          {
            lower_.values_[ col ].mmv( mv[ lower_.cols_[ col ] ], rhs );
          }

          mv[ i ] = rhs;  // Lii = I
        }

        copyOwnerToAll( mv );

        for( size_type i=0; i<iEnd; ++ i )
        {
            vblock& vBlock = mv[ lastRow - i ];
            vblock rhs ( vBlock );
            const size_type rowI     = upper_.rows_[ i ];
            const size_type rowINext = upper_.rows_[ i+1 ];

            for( size_type col = rowI; col < rowINext; ++ col )
            {
                upper_.values_[ col ].mmv( mv[ upper_.cols_[ col ] ], rhs );
            }

            // apply inverse and store result
            inv_[ i ].mv( rhs, vBlock);
        }

        copyOwnerToAll( mv );

        if( relaxation_ ) {
            mv *= w_;
        }
        reorderBack(mv, v);
    }

    template <class V>
    void copyOwnerToAll( V& v ) const
    {
        if( comm_ ) {
            comm_->copyOwnerToAll(v, v);
        }
    }

    /*!
      \brief Clean up.

      \copydoc Preconditioner::post(X&)
    */
    virtual void post (Range& x)
    {
        DUNE_UNUSED_PARAMETER(x);
    }

protected:
    void init( const Matrix& A, const int iluIteration, MILU_VARIANT milu, bool redBlack, bool reorderSpheres )
    {
        // (For older DUNE versions the communicator might be
        // invalid if redistribution in AMG happened on the coarset level.
        // Therefore we check for nonzero size
        if ( comm_ && comm_->communicator().size()<=0 )
        {
            if ( A.N() > 0 )
            {
                OPM_THROW(std::logic_error, "Expected a matrix with zero rows for an invalid communicator.");
            }
            else
            {
                // simply set the communicator to null
                comm_ = nullptr;
            }
        }

        int ilu_setup_successful = 1;
        std::string message;
        const int rank = ( comm_ ) ? comm_->communicator().rank() : 0;

        std::unique_ptr< Matrix > ILU;

        if ( redBlack )
        {
            using Graph = Dune::Amg::MatrixGraph<const Matrix>;
            Graph graph(A);
            auto colorsTuple = colorVerticesWelshPowell(graph);
            const auto& colors = std::get<0>(colorsTuple);
            const auto& verticesPerColor = std::get<2>(colorsTuple);
            auto noColors = std::get<1>(colorsTuple);
            if ( reorderSpheres )
            {
                ordering_ = reorderVerticesSpheres(colors, noColors, verticesPerColor,
                                                   graph, 0);
            }
            else
            {
                ordering_ = reorderVerticesPreserving(colors, noColors, verticesPerColor,
                                                      graph);
            }
        }

        std::vector<std::size_t> inverseOrdering(ordering_.size());
        std::size_t index = 0;
        for( auto newIndex: ordering_)
        {
            inverseOrdering[newIndex] = index++;
        }

        try
        {
            if( iluIteration == 0 ) {
                // create ILU-0 decomposition
                if ( ordering_.empty() )
                {
                    ILU.reset( new Matrix( A ) );
                }
                else
                {
                    ILU.reset( new Matrix(A.N(), A.M(), A.nonzeroes(), Matrix::row_wise));
                    auto& newA = *ILU;
                    // Create sparsity pattern
                    for(auto iter=newA.createbegin(), iend = newA.createend(); iter != iend; ++iter)
                    {
                        const auto& row = A[inverseOrdering[iter.index()]];
                        for(auto col = row.begin(), cend = row.end(); col != cend; ++col)
                        {
                            iter.insert(ordering_[col.index()]);
                        }
                    }
                    // Copy values
                    for(auto iter = A.begin(), iend = A.end(); iter != iend; ++iter)
                    {
                        auto& newRow = newA[ordering_[iter.index()]];
                        for(auto col = iter->begin(), cend = iter->end(); col != cend; ++col)
                        {
                            newRow[ordering_[col.index()]] = *col;
                        }
                    }
                }

                switch ( milu )
                {
                case MILU_VARIANT::MILU_1:
                    detail::milu0_decomposition ( *ILU);
                    break;
                case MILU_VARIANT::MILU_2:
                    detail::milu0_decomposition ( *ILU, detail::IdentityFunctor(),
                                                  detail::SignFunctor() );
                    break;
                case MILU_VARIANT::MILU_3:
                    detail::milu0_decomposition ( *ILU, detail::AbsFunctor(),
                                                  detail::SignFunctor() );
                    break;
                case MILU_VARIANT::MILU_4:
                    detail::milu0_decomposition ( *ILU, detail::IdentityFunctor(),
                                                  detail::IsPositiveFunctor() );
                    break;
                default:
                    bilu0_decomposition( *ILU );
                    break;
                }
            }
            else {
                // create ILU-n decomposition
                ILU.reset( new Matrix( A.N(), A.M(), Matrix::row_wise) );
                std::unique_ptr<detail::Reorderer> reorderer, inverseReorderer;
                if ( ordering_.empty() )
                {
                    reorderer.reset(new detail::NoReorderer());
                    inverseReorderer.reset(new detail::NoReorderer());
                }
                else
                {
                    reorderer.reset(new detail::RealReorderer(ordering_));
                    inverseReorderer.reset(new detail::RealReorderer(inverseOrdering));
                }

                milun_decomposition( A, iluIteration, milu, *ILU, *reorderer, *inverseReorderer );
            }
        }
        catch (const Dune::MatrixBlockError& error)
        {
            message = error.what();
            std::cerr<<"Exception occured on process " << rank << " during " <<
                "setup of ILU0 preconditioner with message: " <<
                message<<std::endl;
            ilu_setup_successful = 0;
        }

        // Check whether there was a problem on some process
        const bool parallel_failure = comm_ && comm_->communicator().min(ilu_setup_successful) == 0;
        const bool local_failure = ilu_setup_successful == 0;
        if ( local_failure || parallel_failure )
        {
            throw Dune::MatrixBlockError();
        }

        // store ILU in simple CRS format
        detail::convertToCRS( *ILU, lower_, upper_, inv_ );
    }

    /// \brief Reorder D if needed and return a reference to it.
    Range& reorderD(const Range& d)
    {
        if ( ordering_.empty())
        {
            // As d is non-const in the apply method of the
            // solver casting away constness in this particular
            // setting is not undefined. It is ugly though but due
            // to the preconditioner interface of dune-istl.
            return const_cast<Range&>(d);
        }
        else
        {
            reorderedD_.resize(d.size());
            std::size_t i = 0;
            for(auto index: ordering_)
            {
                reorderedD_[index]=d[i++];
            }
            return reorderedD_;
        }
    }

    /// \brief Reorder V if needed and return a reference to it.
    Domain& reorderV(Domain& v)
    {
        if ( ordering_.empty())
        {
            return v;
        }
        else
        {
            reorderedV_.resize(v.size());
            std::size_t i = 0;
            for(auto index: ordering_)
            {
                reorderedV_[index]=v[i++];
            }
            return reorderedV_;
        }
    }

    void reorderBack(const Range& reorderedV, Range& v)
    {
        if ( !ordering_.empty() )
        {
            std::size_t i = 0;
            for(auto index: ordering_)
            {
                v[i++] = reorderedV[index];
            }
        }
    }
protected:
    //! \brief The ILU0 decomposition of the matrix.
    CRS lower_;
    CRS upper_;
    std::vector< block_type > inv_;
    //! \brief the reordering of the unknowns
    std::vector< std::size_t > ordering_;
    //! \brief The reordered right hand side
    Range reorderedD_;
    //! \brief The reordered left hand side.
    Domain reorderedV_;

    const ParallelInfo* comm_;
    //! \brief The relaxation factor to use.
    const field_type w_;
    const bool relaxation_;

};

} // end namespace Opm
#endif
