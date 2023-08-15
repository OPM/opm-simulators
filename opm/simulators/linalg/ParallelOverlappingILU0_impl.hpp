/*
  Copyright 2015, 2022 Dr. Blatt - HPC-Simulation-Software & Services
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

#include <opm/simulators/linalg/ParallelOverlappingILU0.hpp>

#include <dune/common/version.hh>

#include <dune/istl/ilu.hh>
#include <dune/istl/owneroverlapcopy.hh>

#include <opm/common/ErrorMacros.hpp>
#include <opm/common/TimingMacros.hpp>

#include <opm/simulators/linalg/GraphColoring.hpp>
#include <opm/simulators/linalg/matrixblock.hh>

namespace Opm
{
namespace detail
{

//! Compute Blocked ILU0 decomposition, when we know junk ghost rows are located at the end of A
template<class M>
void ghost_last_bilu0_decomposition (M& A, std::size_t interiorSize)
{
    // iterator types
    using rowiterator = typename M::RowIterator;
    using coliterator = typename M::ColIterator;
    using block = typename M::block_type;

    // implement left looking variant with stored inverse
    for (rowiterator i = A.begin(); i.index() < interiorSize; ++i)
    {
        // coliterator is diagonal after the following loop
        coliterator endij=(*i).end();           // end of row i
        coliterator ij;

        // eliminate entries left of diagonal; store L factor
        for (ij=(*i).begin(); ij.index()<i.index(); ++ij)
        {
            // find A_jj which eliminates A_ij
            coliterator jj = A[ij.index()].find(ij.index());

            // compute L_ij = A_jj^-1 * A_ij
            (*ij).rightmultiply(*jj);

            // modify row
            coliterator endjk=A[ij.index()].end();    // end of row j
            coliterator jk=jj; ++jk;
            coliterator ik=ij; ++ik;
            while (ik!=endij && jk!=endjk)
                if (ik.index()==jk.index())
                {
                    block B(*jk);
                    B.leftmultiply(*ij);
                    *ik -= B;
                    ++ik; ++jk;
                }
                else
                {
                    if (ik.index()<jk.index())
                        ++ik;
                    else
                        ++jk;
                }
        }

        // invert pivot and store it in A
        if (ij.index()!=i.index())
            DUNE_THROW(Dune::ISTLError,"diagonal entry missing");
        try {
            (*ij).invert();   // compute inverse of diagonal block
        }
        catch (Dune::FMatrixError & e) {
            DUNE_THROW(Dune::ISTLError,"ILU failed to invert matrix block");
        }
    }
}

//! compute ILU decomposition of A. A is overwritten by its decomposition
template<class M, class CRS, class InvVector>
void convertToCRS(const M& A, CRS& lower, CRS& upper, InvVector& inv)
{
  OPM_TIMEBLOCK(convertToCRS);  
  // No need to do anything for 0 rows. Return to prevent indexing a
  // a zero sized array.
  if ( A.N() == 0 )
  {
    return;
  }

  using size_type = typename M :: size_type;

  lower.clear();
  upper.clear();
  inv.clear();
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


template<class Matrix, class Domain, class Range, class ParallelInfoT>
Dune::SolverCategory::Category
ParallelOverlappingILU0<Matrix,Domain,Range,ParallelInfoT>::category() const
{
    return std::is_same_v<ParallelInfoT, Dune::Amg::SequentialInformation> ?
           Dune::SolverCategory::sequential : Dune::SolverCategory::overlapping;
}

template<class Matrix, class Domain, class Range, class ParallelInfoT>
ParallelOverlappingILU0<Matrix,Domain,Range,ParallelInfoT>::
ParallelOverlappingILU0(const Matrix& A,
                        const int n, const field_type w,
                        MILU_VARIANT milu, bool redblack,
                        bool reorder_sphere)
    : lower_(),
      upper_(),
      inv_(),
      comm_(nullptr), w_(w),
      relaxation_( std::abs( w - 1.0 ) > 1e-15 ),
      A_(&reinterpret_cast<const Matrix&>(A)), iluIteration_(n),
      milu_(milu), redBlack_(redblack), reorderSphere_(reorder_sphere)
{
    interiorSize_ = A.N();
    // BlockMatrix is a Subclass of FieldMatrix that just adds
    // methods. Therefore this cast should be safe.
    update();
}

template<class Matrix, class Domain, class Range, class ParallelInfoT>
ParallelOverlappingILU0<Matrix,Domain,Range,ParallelInfoT>::
ParallelOverlappingILU0(const Matrix& A,
                        const ParallelInfo& comm, const int n, const field_type w,
                        MILU_VARIANT milu, bool redblack,
                        bool reorder_sphere)
    : lower_(),
      upper_(),
      inv_(),
      comm_(&comm), w_(w),
      relaxation_( std::abs( w - 1.0 ) > 1e-15 ),
      A_(&reinterpret_cast<const Matrix&>(A)), iluIteration_(n),
      milu_(milu), redBlack_(redblack), reorderSphere_(reorder_sphere)
{
    interiorSize_ = A.N();
    // BlockMatrix is a Subclass of FieldMatrix that just adds
    // methods. Therefore this cast should be safe.
    update();
}

template<class Matrix, class Domain, class Range, class ParallelInfoT>
ParallelOverlappingILU0<Matrix,Domain,Range,ParallelInfoT>::
ParallelOverlappingILU0(const Matrix& A,
                        const field_type w, MILU_VARIANT milu, bool redblack,
                        bool reorder_sphere)
    : ParallelOverlappingILU0( A, 0, w, milu, redblack, reorder_sphere )
{}

template<class Matrix, class Domain, class Range, class ParallelInfoT>
ParallelOverlappingILU0<Matrix,Domain,Range,ParallelInfoT>::
ParallelOverlappingILU0(const Matrix& A,
                        const ParallelInfo& comm, const field_type w,
                        MILU_VARIANT milu, bool redblack,
                        bool reorder_sphere)
    : lower_(),
      upper_(),
      inv_(),
      comm_(&comm), w_(w),
      relaxation_( std::abs( w - 1.0 ) > 1e-15 ),
      A_(&reinterpret_cast<const Matrix&>(A)), iluIteration_(0),
      milu_(milu), redBlack_(redblack), reorderSphere_(reorder_sphere)
{
    interiorSize_ = A.N();
    // BlockMatrix is a Subclass of FieldMatrix that just adds
    // methods. Therefore this cast should be safe.
    update();
}

template<class Matrix, class Domain, class Range, class ParallelInfoT>
ParallelOverlappingILU0<Matrix,Domain,Range,ParallelInfoT>::
ParallelOverlappingILU0(const Matrix& A,
                        const ParallelInfo& comm,
                        const field_type w, MILU_VARIANT milu,
                        size_type interiorSize, bool redblack,
                        bool reorder_sphere)
    : lower_(),
      upper_(),
      inv_(),
      comm_(&comm), w_(w),
      relaxation_( std::abs( w - 1.0 ) > 1e-15 ),
      interiorSize_(interiorSize),
      A_(&reinterpret_cast<const Matrix&>(A)), iluIteration_(0),
      milu_(milu), redBlack_(redblack), reorderSphere_(reorder_sphere)
{
    // BlockMatrix is a Subclass of FieldMatrix that just adds
    // methods. Therefore this cast should be safe.
    update( );
}

template<class Matrix, class Domain, class Range, class ParallelInfoT>
void ParallelOverlappingILU0<Matrix,Domain,Range,ParallelInfoT>::
apply (Domain& v, const Range& d)
{
    OPM_TIMEBLOCK(apply);
    Range& md = reorderD(d);
    Domain& mv = reorderV(v);

    // iterator types
    using dblock = typename Range ::block_type;
    using vblock = typename Domain::block_type;

    const size_type iEnd = lower_.rows();
    const size_type lastRow = iEnd - 1;
    size_type upperLoopStart = iEnd - interiorSize_;
    size_type lowerLoopEnd = interiorSize_;
    if (iEnd != upper_.rows())
    {
        OPM_THROW(std::logic_error,"ILU: number of lower and upper rows must be the same");
    }

    // lower triangular solve
    for (size_type i = 0; i < lowerLoopEnd; ++i)
    {
        dblock rhs( md[ i ] );
        const size_type rowI     = lower_.rows_[ i ];
        const size_type rowINext = lower_.rows_[ i+1 ];

        for (size_type col = rowI; col < rowINext; ++col)
        {
            lower_.values_[ col ].mmv( mv[ lower_.cols_[ col ] ], rhs );
        }

        mv[ i ] = rhs;  // Lii = I
    }

    for (size_type i = upperLoopStart; i < iEnd; ++i)
    {
        vblock& vBlock = mv[ lastRow - i ];
        vblock rhs ( vBlock );
        const size_type rowI     = upper_.rows_[ i ];
        const size_type rowINext = upper_.rows_[ i+1 ];

        for (size_type col = rowI; col < rowINext; ++col)
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

template<class Matrix, class Domain, class Range, class ParallelInfoT>
template<class V>
void ParallelOverlappingILU0<Matrix,Domain,Range,ParallelInfoT>::
copyOwnerToAll(V& v) const
{
    if( comm_ ) {
        comm_->copyOwnerToAll(v, v);
    }
}

template<class Matrix, class Domain, class Range, class ParallelInfoT>
void ParallelOverlappingILU0<Matrix,Domain,Range,ParallelInfoT>::
update()
{
    OPM_TIMEBLOCK(update);
    // (For older DUNE versions the communicator might be
    // invalid if redistribution in AMG happened on the coarset level.
    // Therefore we check for nonzero size
    if (comm_ && comm_->communicator().size() <= 0)
    {
        if (A_->N() > 0)
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
    const int rank = comm_ ? comm_->communicator().rank() : 0;

    if (redBlack_)
    {
        using Graph = Dune::Amg::MatrixGraph<const Matrix>;
        Graph graph(*A_);
        auto colorsTuple = colorVerticesWelshPowell(graph);
        const auto& colors = std::get<0>(colorsTuple);
        const auto& verticesPerColor = std::get<2>(colorsTuple);
        auto noColors = std::get<1>(colorsTuple);
        if ( reorderSphere_ )
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
    for (const auto newIndex : ordering_)
    {
        inverseOrdering[newIndex] = index++;
    }

    try
    {
        OPM_TIMEBLOCK(iluDecomposition);
        if (iluIteration_ == 0) {
            // create ILU-0 decomposition
            if (ordering_.empty())
            {
                if (ILU_) {
                    OPM_TIMEBLOCK(iluDecompositionMakeMatrix);
                    // The ILU_ matrix is already a copy with the same
                    // sparse structure as A_, but the values of A_ may
                    // have changed, so we must copy all elements.
                    for (std::size_t row = 0; row < A_->N(); ++row) {
                        const auto& Arow = (*A_)[row];
                        auto& ILUrow = (*ILU_)[row];
                        auto Ait = Arow.begin();
                        auto Iit = ILUrow.begin();
                        for (; Ait != Arow.end(); ++Ait, ++Iit) {
                            *Iit = *Ait;
                        }
                    }
                } else {
                    // First call, must duplicate matrix.
                    ILU_ = std::make_unique<Matrix>(*A_);
                }
            }
            else
            {
                ILU_ = std::make_unique<Matrix>(A_->N(), A_->M(),
                                                A_->nonzeroes(), Matrix::row_wise);
                auto& newA = *ILU_;
                // Create sparsity pattern
                for (auto iter = newA.createbegin(), iend = newA.createend(); iter != iend; ++iter)
                {
                    const auto& row = (*A_)[inverseOrdering[iter.index()]];
                    for (auto col = row.begin(), cend = row.end(); col != cend; ++col)
                    {
                        iter.insert(ordering_[col.index()]);
                    }
                }
                // Copy values
                for (auto iter = A_->begin(), iend = A_->end(); iter != iend; ++iter)
                {
                    auto& newRow = newA[ordering_[iter.index()]];
                    for (auto col = iter->begin(), cend = iter->end(); col != cend; ++col)
                    {
                        newRow[ordering_[col.index()]] = *col;
                    }
                }
            }

            switch (milu_)
            {
            case MILU_VARIANT::MILU_1:
                detail::milu0_decomposition ( *ILU_);
                break;
            case MILU_VARIANT::MILU_2:
                detail::milu0_decomposition ( *ILU_, detail::identityFunctor<typename Matrix::field_type>,
                                              detail::signFunctor<typename Matrix::field_type> );
                break;
            case MILU_VARIANT::MILU_3:
                detail::milu0_decomposition ( *ILU_, detail::absFunctor<typename Matrix::field_type>,
                                              detail::signFunctor<typename Matrix::field_type> );
                break;
            case MILU_VARIANT::MILU_4:
                detail::milu0_decomposition ( *ILU_, detail::identityFunctor<typename Matrix::field_type>,
                                              detail::isPositiveFunctor<typename Matrix::field_type> );
                break;
            default:
                if (interiorSize_ == A_->N())
#if DUNE_VERSION_LT(DUNE_GRID, 2, 8)
                    bilu0_decomposition( *ILU_ );
#else
                    Dune::ILU::blockILU0Decomposition( *ILU_ );
#endif
                else
                    detail::ghost_last_bilu0_decomposition(*ILU_, interiorSize_);
                break;
            }
        }
        else {
            // create ILU-n decomposition
            ILU_ = std::make_unique<Matrix>(A_->N(), A_->M(), Matrix::row_wise);
            std::unique_ptr<detail::Reorderer> reorderer, inverseReorderer;
            if (ordering_.empty())
            {
                reorderer.reset(new detail::NoReorderer());
                inverseReorderer.reset(new detail::NoReorderer());
            }
            else
            {
                reorderer.reset(new detail::RealReorderer(ordering_));
                inverseReorderer.reset(new detail::RealReorderer(inverseOrdering));
            }

            milun_decomposition( *A_, iluIteration_, milu_, *ILU_, *reorderer, *inverseReorderer );
        }
    }
    catch (const Dune::MatrixBlockError& error)
    {
        message = error.what();
        std::cerr << "Exception occurred on process " << rank << " during " <<
                     "setup of ILU0 preconditioner with message: "
                  << message<<std::endl;
        ilu_setup_successful = 0;
    }

    // Check whether there was a problem on some process
    const bool parallel_failure = comm_ && comm_->communicator().min(ilu_setup_successful) == 0;
    const bool local_failure = ilu_setup_successful == 0;
    if (local_failure || parallel_failure)
    {
        throw Dune::MatrixBlockError();
    }

    // store ILU in simple CRS format
    detail::convertToCRS(*ILU_, lower_, upper_, inv_);
}

template<class Matrix, class Domain, class Range, class ParallelInfoT>
Range& ParallelOverlappingILU0<Matrix,Domain,Range,ParallelInfoT>::
reorderD(const Range& d)
{
    if (ordering_.empty())
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
        for (const auto index : ordering_)
        {
            reorderedD_[index] = d[i++];
        }
        return reorderedD_;
    }
}

template<class Matrix, class Domain, class Range, class ParallelInfoT>
Domain& ParallelOverlappingILU0<Matrix,Domain,Range,ParallelInfoT>::
reorderV(Domain& v)
{
    if (ordering_.empty())
    {
        return v;
    }
    else
    {
        reorderedV_.resize(v.size());
        std::size_t i = 0;
        for (const auto index : ordering_)
        {
            reorderedV_[index] = v[i++];
        }
        return reorderedV_;
    }
}

template<class Matrix, class Domain, class Range, class ParallelInfoT>
void ParallelOverlappingILU0<Matrix,Domain,Range,ParallelInfoT>::
reorderBack(const Range& reorderedV, Range& v)
{
    if (!ordering_.empty())
    {
        std::size_t i = 0;
        for (const auto index : ordering_)
        {
            v[i++] = reorderedV[index];
        }
    }
}

} // end namespace Opm
