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
#include <opm/common/TimingMacros.hpp>
#include <opm/simulators/linalg/MILU.hpp>
#include <opm/simulators/linalg/PreconditionerWithUpdate.hpp>
#include <dune/istl/paamg/smoother.hh>

#include <cstddef>
#include <vector>
#include <type_traits>

namespace Opm
{

//template<class M, class X, class Y, class C>
//class ParallelOverlappingILU0;
template<class Matrix, class Domain, class Range, class ParallelInfo>
class ParallelOverlappingILU0;

template<class F>
class ParallelOverlappingILU0Args
    : public Dune::Amg::DefaultSmootherArgs<F>
{
 public:
    ParallelOverlappingILU0Args(MILU_VARIANT milu = MILU_VARIANT::ILU )
        : milu_(milu), n_(0)
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
    using T = Opm::ParallelOverlappingILU0<Matrix,Domain,Range,ParallelInfo>;
    using Arguments = DefaultParallelConstructionArgs<T,ParallelInfo>;

    using ParallelOverlappingILU0Pointer = std::shared_ptr<T>;

    static inline ParallelOverlappingILU0Pointer construct(Arguments& args)
    {
        return ParallelOverlappingILU0Pointer(
                new T(args.getMatrix(),
                      args.getComm(),
                      args.getArgs().getN(),
                      args.getArgs().relaxationFactor,
                      args.getArgs().getMilu()) );
    }
};

} // end namespace Amg

} // end namespace Dune

namespace Opm
{
/// \brief A two-step version of an overlapping Schwarz preconditioner using one step ILU0 as
///
/// This preconditioner differs from a ParallelRestrictedOverlappingSchwarz with
/// Dune:SeqILU0 in the following way:
/// During apply we make sure that the current residual is consistent (i.e.
/// each process knows the same value for each index. Then we solve
/// Ly = d for y and make y consistent again. Last we solve Ux = y and
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
    : public Dune::PreconditionerWithUpdate<Domain,Range>
{
    using ParallelInfo = ParallelInfoT;

public:
    //! \brief The matrix type the preconditioner is for.
    using matrix_type = typename std::remove_const<Matrix>::type;
    //! \brief The domain type of the preconditioner.
    using domain_type = Domain;
    //! \brief The range type of the preconditioner.
    using range_type = Range;
    //! \brief The field type of the preconditioner.
    using field_type = typename Domain::field_type;

    using block_type = typename matrix_type::block_type;
    using size_type = typename matrix_type::size_type;

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

      void clear()
      {
          rows_.clear();
          values_.clear();
          cols_.clear();
          nRows_= 0;
      }

      std::vector<size_type> rows_;
      std::vector<block_type> values_;
      std::vector<size_type> cols_;
      size_type nRows_;
    };

public:
    Dune::SolverCategory::Category category() const override;

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
    ParallelOverlappingILU0 (const Matrix& A,
                             const int n, const field_type w,
                             MILU_VARIANT milu, bool redblack = false,
                             bool reorder_sphere = true);

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
    ParallelOverlappingILU0 (const Matrix& A,
                             const ParallelInfo& comm, const int n, const field_type w,
                             MILU_VARIANT milu, bool redblack = false,
                             bool reorder_sphere = true);

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
    ParallelOverlappingILU0 (const Matrix& A,
                             const field_type w, MILU_VARIANT milu,
                             bool redblack = false,
                             bool reorder_sphere = true);

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
    ParallelOverlappingILU0 (const Matrix& A,
                             const ParallelInfo& comm, const field_type w,
                             MILU_VARIANT milu, bool redblack = false,
                             bool reorder_sphere = true);

    /*! \brief Constructor.

      Constructor gets all parameters to operate the prec.
      \param A The matrix to operate on.
      \param n ILU fill in level (for testing). This does not work in parallel.
      \param w The relaxation factor.
      \param milu The modified ILU variant to use. 0 means traditional ILU. \see MILU_VARIANT.
      \param interiorSize The number of interior/owner rows in the matrix.
      \param redblack Whether to use a red-black ordering.
      \param reorder_sphere If true, we start the reordering at a root node.
                            The vertices on each layer aound it (same distance) are
                            ordered consecutivly. If false, we preserver the order of
                            the vertices with the same color.
    */
    ParallelOverlappingILU0 (const Matrix& A,
                             const ParallelInfo& comm,
                             const field_type w, MILU_VARIANT milu,
                             size_type interiorSize, bool redblack = false,
                             bool reorder_sphere = true);

    /*!
      \brief Prepare the preconditioner.

      \copydoc Preconditioner::pre(X&,Y&)
    */
    void pre (Domain&, Range&) override
    {}

    /*!
      \brief Apply the preconditoner.

      \copydoc Preconditioner::apply(X&,const Y&)
    */
    void apply (Domain& v, const Range& d) override;

    template <class V>
    void copyOwnerToAll( V& v ) const;

    /*!
      \brief Clean up.

      \copydoc Preconditioner::post(X&)
    */
    void post (Range&) override
    {}

    void update() override;

protected:
    /// \brief Reorder D if needed and return a reference to it.
    Range& reorderD(const Range& d);

    /// \brief Reorder V if needed and return a reference to it.
    Domain& reorderV(Domain& v);

    void reorderBack(const Range& reorderedV, Range& v);

    //! \brief The ILU0 decomposition of the matrix.
    std::unique_ptr<Matrix> ILU_;
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
    size_type interiorSize_;
    const Matrix* A_;
    int iluIteration_;
    MILU_VARIANT milu_;
    bool redBlack_;
    bool reorderSphere_;
};

} // end namespace Opm
#endif
