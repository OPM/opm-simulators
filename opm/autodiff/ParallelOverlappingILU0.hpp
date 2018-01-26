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

#include <opm/common/Exceptions.hpp>

#include <dune/common/version.hh>
#include <dune/istl/preconditioner.hh>
#include <dune/istl/paamg/smoother.hh>
#include <dune/istl/paamg/pinfo.hh>

#include <type_traits>

namespace Opm
{

//template<class M, class X, class Y, class C>
//class ParallelOverlappingILU0;
template<class Matrix, class Domain, class Range, class ParallelInfo = Dune::Amg::SequentialInformation>
class ParallelOverlappingILU0;

} // end namespace Opm

namespace Dune
{

namespace Amg
{

/// \brief Tells AMG how to construct the Opm::ParallelOverlappingILU0 smoother
/// \tparam Matrix The type of the Matrix.
/// \tparam Domain The type of the Vector representing the domain.
/// \tparam Range The type of the Vector representing the range.
/// \tparam ParallelInfo The type of the parallel information object
///         used, e.g. Dune::OwnerOverlapCommunication
template<class Matrix, class Domain, class Range, class ParallelInfo>
struct ConstructionTraits<Opm::ParallelOverlappingILU0<Matrix,Domain,Range,ParallelInfo> >
{
    typedef Dune::SeqILU0<Matrix,Domain,Range> T;
    typedef DefaultParallelConstructionArgs<T,ParallelInfo> Arguments;
    typedef ConstructionTraits<T> SeqConstructionTraits;
    static inline Opm::ParallelOverlappingILU0<Matrix,Domain,Range,ParallelInfo>* construct(Arguments& args)
    {
        return new Opm::ParallelOverlappingILU0<Matrix,Domain,Range,ParallelInfo>(args.getMatrix(),
                                                         args.getComm(),
                                                         args.getArgs().relaxationFactor);
    }

    static inline void deconstruct(Opm::ParallelOverlappingILU0<Matrix,Domain,Range,ParallelInfo>* bp)
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
      //! compute ILU decomposition of A. A is overwritten by its decomposition
      template<class M, class CRS, class InvVector>
      void convertToCRS(const M& A, CRS& lower, CRS& upper, InvVector& inv )
      {
        typedef typename M :: size_type size_type;

        lower.resize( A.N() );
        upper.resize( A.N() );
        inv.resize( A.N() );

        lower.reserveAdditional( 2*A.N() );

        // implement left looking variant with stored inverse
        const auto endi = A.end();
        size_type row = 0;
        size_type colcount = 0;
        lower.rows_[ 0 ] = colcount;
        for (auto i=A.begin(); i!=endi; ++i, ++row)
        {
          const size_type iIndex  = i.index();
          lower.reserveAdditional( (*i).size() );

          // eliminate entries left of diagonal; store L factor
          for (auto j=(*i).begin(); j.index() < iIndex; ++j )
          {
            lower.push_back( (*j), j.index() );
            ++colcount;
          }
          lower.rows_[ iIndex+1 ] = colcount;
        }

        const auto rendi = A.beforeBegin();
        row = 0;
        colcount = 0;
        upper.rows_[ 0 ] = colcount ;

        upper.reserveAdditional( lower.nonZeros() + A.N() );

        // NOTE: upper and inv store entries in reverse order, reverse here
        // relative to ILU
        for (auto i=A.beforeEnd(); i!=rendi; --i, ++ row )
        {
          const size_type iIndex = i.index();
          upper.reserveAdditional( (*i).size() );

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
    */
    template<class BlockType, class Alloc>
    ParallelOverlappingILU0 (const Dune::BCRSMatrix<BlockType,Alloc>& A,
                             const int n, const field_type w )
        : lower_(),
          upper_(),
          inv_(),
          comm_(nullptr), w_(w),
          relaxation_( std::abs( w - 1.0 ) > 1e-15 )
    {
        // BlockMatrix is a Subclass of FieldMatrix that just adds
        // methods. Therefore this cast should be safe.
        init( reinterpret_cast<const Matrix&>(A), n );
    }

    /*! \brief Constructor.

      Constructor gets all parameters to operate the prec.
      \param A The matrix to operate on.
      \param w The relaxation factor.
    */
    template<class BlockType, class Alloc>
    ParallelOverlappingILU0 (const Dune::BCRSMatrix<BlockType,Alloc>& A,
                             const field_type w)
        : ParallelOverlappingILU0( A, 0, w )
    {
    }

    /*! \brief Constructor.

      Constructor gets all parameters to operate the prec.
      \param A      The matrix to operate on.
      \param comm   communication object, e.g. Dune::OwnerOverlapCopyCommunication
      \param w      The relaxation factor.
    */
    template<class BlockType, class Alloc>
    ParallelOverlappingILU0 (const Dune::BCRSMatrix<BlockType,Alloc>& A,
                             const ParallelInfo& comm, const field_type w)
        : lower_(),
          upper_(),
          inv_(),
          comm_(&comm), w_(w),
          relaxation_( std::abs( w - 1.0 ) > 1e-15 )
    {
        // BlockMatrix is a Subclass of FieldMatrix that just adds
        // methods. Therefore this cast should be safe.
        init( reinterpret_cast<const Matrix&>(A), 0 );
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
        Range& md = const_cast<Range&>(d);
        copyOwnerToAll( md );

        // iterator types
        typedef typename Range ::block_type  dblock;
        typedef typename Domain::block_type  vblock;

        const size_type iEnd = lower_.rows();
        const size_type lastRow = iEnd - 1;
        if( iEnd != upper_.rows() )
        {
            std::abort();
           // OPM_THROW(std::logic_error,"ILU: lower and upper rows must be the same");
        }

        // lower triangular solve
        for( size_type i=0; i<iEnd; ++ i )
        {
          dblock rhs( d[ i ] );
          const size_type rowI     = lower_.rows_[ i ];
          const size_type rowINext = lower_.rows_[ i+1 ];

          for( size_type col = rowI; col < rowINext; ++ col )
          {
            lower_.values_[ col ].mmv( v[ lower_.cols_[ col ] ], rhs );
          }

          v[ i ] = rhs;  // Lii = I
        }

        copyOwnerToAll( v );

        for( size_type i=0; i<iEnd; ++ i )
        {
            vblock& vBlock = v[ lastRow - i ];
            vblock rhs ( vBlock );
            const size_type rowI     = upper_.rows_[ i ];
            const size_type rowINext = upper_.rows_[ i+1 ];

            for( size_type col = rowI; col < rowINext; ++ col )
            {
                upper_.values_[ col ].mmv( v[ upper_.cols_[ col ] ], rhs );
            }

            // apply inverse and store result
            inv_[ i ].mv( rhs, vBlock);
        }

        copyOwnerToAll( v );

        if( relaxation_ ) {
            v *= w_;
        }
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
    void init( const Matrix& A, const int iluIteration )
    {
        int ilu_setup_successful = 1;
        std::string message;
        const int rank = ( comm_ ) ? comm_->communicator().rank() : 0;

        std::unique_ptr< Matrix > ILU;

        try
        {
            if( iluIteration == 0 ) {
                // create ILU-0 decomposition
                ILU.reset( new Matrix( A ) );
                bilu0_decomposition( *ILU );
            }
            else {
                // create ILU-n decomposition
                ILU.reset( new Matrix( A.N(), A.M(), Matrix::row_wise) );
                bilu_decomposition( A, iluIteration, *ILU );
            }
        }
        catch ( Dune::MatrixBlockError error )
        {
            message = error.what();
            std::cerr<<"Exception occured on process " << rank << " during " <<
                "setup of ILU0 preconditioner with message: " <<
                message<<std::endl;
            ilu_setup_successful = 0;
        }

        // Check whether there was a problem on some process
        if ( comm_ && comm_->communicator().min(ilu_setup_successful) == 0 )
        {
            throw Dune::MatrixBlockError();
        }

        // store ILU in simple CRS format
        detail::convertToCRS( *ILU, lower_, upper_, inv_ );
    }

protected:
    //! \brief The ILU0 decomposition of the matrix.
    CRS lower_;
    CRS upper_;
    std::vector< block_type > inv_;

    const ParallelInfo* comm_;
    //! \brief The relaxation factor to use.
    const field_type w_;
    const bool relaxation_;

};

} // end namespace Opm
#endif
