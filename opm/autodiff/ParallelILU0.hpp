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
#ifndef OPM_PARALLELILU0_HEADER_INCLUDED
#define OPM_PARALLELILU0_HEADER_INCLUDED

#include <opm/autodiff/SendReceiveCommunicator.hpp>
#include <dune/istl/preconditioner.hh>
#include <dune/istl/paamg/smoother.hh>
#include <dune/common/exceptions.hh>

#if HAVE_MPI
namespace Opm
{

namespace Detail
{

/// \brief Uses the process rank as the label used for global ordering.
class IdProcessMap
{
public:
    template<class ParallelInfo>
    explicit IdProcessMap(ParallelInfo const& comm)
        : rank_(comm.rank())
    {}

    /// \brief Get the label on this process
    std::size_t myMapping() const
    {
        return rank_;
    }

    /// \brief Compare our label with the one of another process.
    /// \tparam Compare Type of the functor for comparison, e.g. std::less<size_t>
    /// \param other_proc The label of the other proess.
    /// \param compare functor for comparison.
    template<class Compare>
    bool compareWithOtherLabel(std::size_t other_proc,
                               Compare const& compare) const
    {
        return compare(other_proc, myMapping());
    }

private:
    const int rank_;
};

class LinearSystemPermutation;

/// \brief A tool for looking up the global indices of a permuted index range.
/// \tparam IndexSet The type of the underlying index set.
template<class IndexSet>
class PermutedGlobalLookupIndexSet
{
public:
    typedef typename IndexSet::GlobalIndex GlobalIndex;
    typedef typename IndexSet::IndexPair IndexPair;
    typedef typename IndexSet::const_iterator const_iterator;

    PermutedGlobalLookupIndexSet(const IndexSet& indexSet,
                                 std::size_t size,
                                 const LinearSystemPermutation& permutation)
        : indexSet_(indexSet), index_pairs_(size, nullptr)
    {
        for( const auto& pair : indexSet_)
        {
            assert ( pair.local() < index_pairs_.size() );
            index_pairs_[ permutation[pair.local()] ] = &pair;
        }
    }

    /**
     * \brief Get the index pair corresponding to a local index.
     * \param The permuted local index
     */
    const IndexPair* pair(const std::size_t& local) const
    {
        return index_pairs_[local];
    }

    const_iterator begin() const
    {
        return indexSet_.begin();
    }

    const_iterator end() const
    {
        return indexSet_.end();
    }
private:
    const IndexSet& indexSet_;
    std::vector<const IndexPair*> index_pairs_;
};

/// \brief A handle to communicate rows of a permuted matrix.
template<class Matrix, class IndexSet>
struct PermutedSparseMatrixHandle
{
    typedef PermutedGlobalLookupIndexSet<IndexSet> GlobalLookupIndexSet;
    typedef typename GlobalLookupIndexSet::GlobalIndex GlobalIndex;
    typedef typename Matrix::block_type block_type;
    typedef std::pair<block_type, GlobalIndex> DataType;

    PermutedSparseMatrixHandle(Matrix& A,
                               PermutedGlobalLookupIndexSet<IndexSet> const& indexset,
                               LinearSystemPermutation const& row_permutation)
        : A_(A), indexset_(indexset),
          row_permutation_(row_permutation)
    {
        MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
    }

    bool fixedsize() const
    {
        return false;
    }

    std::size_t size(std::size_t i) const
    {
        return A_[ i ].size(); // Interface index i is already permuted!
    }

    template<class B>
    void gather(B& buffer, std::size_t i) const
    {
            auto& row = A_[ i ]; // Interface index i is already permuted!
            for( auto entry = row.begin(), end_entry = row.end();
                 entry != end_entry; ++entry)
            {
                // global lookup is also already permuted
                auto index_pair_ptr = indexset_.pair(entry.index());
                // interior nodes do not need a global index.
                typedef typename GlobalLookupIndexSet::GlobalIndex
                    GlobalIndex;
                GlobalIndex index = std::numeric_limits<GlobalIndex>::max();
                if ( index_pair_ptr )
                {
                    index = index_pair_ptr->global();
                }
                buffer.write(std::make_pair(*entry, index));
            }
        }
        template<class B>
        void scatter(B& buffer, std::size_t i, std::size_t n)
        {
            auto& row = A_[ i ]; // Interface index i is already permuted!
            std::vector<DataType> data(n);
            for( auto& datum : data)
            {
                buffer.read(datum);
            }
            // sort by global index. Moved undefined global indices to the end.
            std::sort(data.begin(), data.end(),
                      [](const DataType& d1, const DataType& d2)
                      { return d1.second < d2.second; });
            // skip undefined global indices at the end.
            typedef typename GlobalLookupIndexSet::GlobalIndex GlobalIndex;
            GlobalIndex index = std::numeric_limits<GlobalIndex>::max();
            while ( data.back().second == std::numeric_limits<GlobalIndex>::max() )
            {
                data.pop_back();
            }

            auto hint = indexset_.begin();
            for( auto& datum : data)
            {
                typedef typename GlobalLookupIndexSet::IndexPair IndexPair;
                // \todo This search for every element seems rather expensive.
                // We should try to store the receive values as they are.
                // and later do only one search per column index.
                auto found =
                    std::lower_bound(hint, indexset_.end(), datum,
                                     [](const IndexPair& d1, const DataType& d2)
                                     { return d1.global() < d2.second; });
                if ( found != indexset_.end() )
                {
                    hint = found; // we have found->global() <= datum.second !
                    if( found->global() == datum.second )
                    {
                        // This columns is in our subdomain.
                        // The index set knows nothing about the permutation.
                        // Therefore we need to apply the permutation to the
                        // local index.
                        row[ row_permutation_[found->local()] ] = datum.first;
                    }
                }
            }
        }
private:
    Matrix& A_;
    GlobalLookupIndexSet const& indexset_;
    LinearSystemPermutation const& row_permutation_;
    int rank_;
};

/// \brief A handle for communicating the sparsity of matrix to a map.
template<class Matrix, class GlobalLookupIndexSet>
class SparsityPatternHandle
{
public:
    typedef typename GlobalLookupIndexSet::GlobalIndex GlobalIndex;
    typedef GlobalIndex DataType;
    typedef typename Matrix::size_type size_type;
    typedef typename std::map<size_type,std::set<size_type> > Map;

    SparsityPatternHandle(const Matrix& from, Map& to,
                          GlobalLookupIndexSet const& indexset,
                          LinearSystemPermutation const& row_permutation,
                          std::vector<size_type> const& reverse_permutation)
        : from_(from), to_(to), indexset_(indexset),
          row_permutation_(row_permutation),
          reverse_permutation_(reverse_permutation)
    {
        MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
    }

    bool fixedsize() const
    {
        return false;
    }

    std::size_t size(std::size_t i) const
    {
        // Interface index i is already permuted -> we need to reverse this
        return from_[ reverse_permutation_[ i ] ].size();
    }

    template<class B>
    void gather(B& buffer, std::size_t i) const
    {
        // Interface index i is already permuted -> we need to reverse this
        auto const& row = from_[ reverse_permutation_[ i ] ];

        for( auto entry = row.begin(), end_entry = row.end();
             entry != end_entry; ++entry)
        {
            // global index lookup
            auto index_pair_ptr = indexset_.pair(entry.index());
            // interior nodes do not need a global index.
            GlobalIndex index = std::numeric_limits<GlobalIndex>::max();
            if ( index_pair_ptr )
            {
                index = index_pair_ptr->global();
            }
            buffer.write(index);
        }
    }

    template<class B>
    void scatter(B& buffer, std::size_t i, std::size_t n)
    {
        auto& row = to_[ i ]; // Interface index i is already permuted!
        std::vector<DataType> data(n);
        for( auto& datum : data)
        {
            buffer.read(datum);
        }
        // sort by global index. Moved undefined global indices to the end.
        std::sort(data.begin(), data.end());
        // skip undefined global indices at the end.
        GlobalIndex index = std::numeric_limits<GlobalIndex>::max();
        while ( data.back() == std::numeric_limits<GlobalIndex>::max() )
        {
            data.pop_back();
        }

        auto hint = indexset_.begin();
        for( auto& datum : data)
        {
            typedef typename GlobalLookupIndexSet::IndexPair IndexPair;
            // \todo This search for every element seems rather expensive.
            // We should try to store the receive values as they are.
            // and later do only one search per column index.
            auto found =
                std::lower_bound(hint, indexset_.end(), datum,
                                 [](const IndexPair& d1, const DataType& d2)
                                 { return d1.global() < d2; });
            if ( found != indexset_.end() )
            {
                hint = found; // we have found->global() <= datum.second !
                if( found->global() == datum )
                {
                    // This columns is in our subdomain.
                    // The index set knows nothing about the permutation.
                    // Therefore we need to apply the permutation to the
                    // local index.
                    row.insert(row_permutation_[ found->local() ]);
                }
            }
        }
    }
private:
    Matrix const& from_;
    Map& to_;
    GlobalLookupIndexSet const& indexset_;
    LinearSystemPermutation const& row_permutation_;
    std::vector<size_type> const& reverse_permutation_;
    int rank_;
};

/**
 * \brief A simple Wrapper of a BCRSMatrix that just makes modified row iterators available for subrang of the rows.
 *
 * This is needed for the parallel ILU implementation which needs to neglect
 * rows at the begin and/or end of the matrix.
 */
template<class Matrix>
class SubMatrix
{
public:
    typedef typename Matrix::RowIterator RowIterator;
    typedef typename Matrix::ColIterator ColIterator;
    typedef typename Matrix::block_type block_type;
    typedef typename Matrix::row_type row_type;
    typedef typename Matrix::size_type size_type;

    SubMatrix(Matrix& mat, const RowIterator& begin, const RowIterator& end)
        : mat_(mat), begin_(begin), end_(end)
    {}

    RowIterator begin()
    {
        return begin_;
    }

    RowIterator end()
    {
        return end_;
    }

    row_type& operator[](size_type i)
    {
        return mat_[i];
    }

private:
    Matrix&          mat_;
    RowIterator begin_;
    RowIterator end_;
};

/// \brief Utility class for reordering a linear system.
class LinearSystemPermutation
{
public:
    LinearSystemPermutation(std::size_t size)
        : row_permutation_(size, std::numeric_limits<std::size_t>::max()),
          no_additional_entries_()
    {}

    template<class Y>
    void permutateOrder(const Y& d, Y& permuted_d)
    {
        for ( std::size_t i = 0; i< row_permutation_.size(); ++i)
        {
            permuted_d[ row_permutation_[i] ] = d[ i ];
        }
    }

    template<class X, class T>
    void permutateOrderBackwardsAndScale(X& v, const X& permuted_v,
                                         const T& w)
    {
        for ( std::size_t i = 0; i< row_permutation_.size(); ++i)
        {
            for ( int j = 0; j < X::block_type::dimension; ++j)
            {
                v[ i ] [ j ] = w * permuted_v[ row_permutation_[i] ] [ j ];
            }
        }
    }

    template<class X>
    void permutateOrderBackwards(X& v, const X& permuted_v)
    {
        for ( std::size_t i = 0; i< row_permutation_.size(); ++i)
        {
            v[ i ] = permuted_v[ row_permutation_[i] ];
        }
    }

    template<class B>
    std::unique_ptr<Dune::BCRSMatrix<B> > createPermutedMatrix(Dune::BCRSMatrix<B> const& A)
    {
        return createPermutedMatrix(A, A.N(), A.N(), createInverseRowPermutation());
    }

    std::vector<std::size_t> createInverseRowPermutation()
    {
        std::vector<std::size_t> inverse_row_permutation(row_permutation_.size());

        for(typename std::vector<size_t>::size_type index=0, end=row_permutation_.size();
            index != end; ++index)
        {
            inverse_row_permutation[ row_permutation_[ index ] ] = index;
        }
        return inverse_row_permutation;
    }

    template<class B>
    std::unique_ptr<Dune::BCRSMatrix<B> > createPermutedMatrix(Dune::BCRSMatrix<B> const& A,
                                                               std::size_t interior_end_index,
                                                               std::size_t interface_start_index,
                                                               std::vector<std::size_t>const& inverse_row_permutation)
    {
        std::unique_ptr<Dune::BCRSMatrix<B> >
            permuted_ptr(new Dune::BCRSMatrix<B>(A.N(), A.M(),
                                                 A.nonzeroes() +
                                                 no_additional_entries_,
                                                 Dune::BCRSMatrix<B>::row_wise));
        Dune::BCRSMatrix<B>& permuted = *permuted_ptr;

        auto permuted_row = permuted.createbegin();

        insertPermutedRowIndices(A, permuted_row, 0,  interior_end_index,
                                 inverse_row_permutation);

        // additional_nonzeros in the overlap row might be created by a call to
        // InteriorInterfacePermutation::exchangeMatrixRows
        for(auto overlap_row = additional_nonzeros_.begin(),
            end  = additional_nonzeros_.end(); overlap_row != end;
            ++overlap_row, ++permuted_row)
        {
            assert( permuted_row.index() == overlap_row->first );

            for( auto col : overlap_row->second )
            {
                permuted_row.insert(col);
            }
        }

        assert( permuted_row.index() == interface_start_index );

        insertPermutedRowIndices(A, permuted_row, interface_start_index,  A.N(),
                                 inverse_row_permutation);

        // Copy the values from A
        for(auto row = A.begin(), end_row = A.end(); row != end_row; ++row)
        {
            auto& permuted_row = permuted[ row_permutation_[ row.index() ] ];
            for(auto col = row->begin(), end_col = row->end(); col != end_col; ++col)
            {
                permuted_row[ row_permutation_[col.index()] ] = *col;
            }
        }
        return permuted_ptr;
    }

    std::size_t operator[](std::size_t i) const
    {
        return row_permutation_[i];
    }

    std::size_t size() const
    {
        return row_permutation_.size();
    }

protected:
    /// \brief The permuted indices
    std::vector<std::size_t> row_permutation_;
    /// \brief Additional nonzero entries to be added
    std::map<std::size_t,std::set<std::size_t> > additional_nonzeros_;
    std::size_t no_additional_entries_;

private:
    template<class Block>
    void insertPermutedRowIndices(Dune::BCRSMatrix<Block> const& orig_matrix,
                                  typename Dune::BCRSMatrix<Block>::CreateIterator& permuted_row ,
                                  std::size_t index, std::size_t end_index,
                                  std::vector<std::size_t> const& inverse_permutation)
    {
        for(; index < end_index; ++index, ++permuted_row)
        {
            const auto& original_row = orig_matrix[ inverse_permutation[ index ] ];
            for(auto col = original_row.begin(), endCol = original_row.end();
                col != endCol; ++col)
            {
                permuted_row.insert(row_permutation_[col.index()]);
            }
        }
        assert( end_index == permuted_row.index() );
    }
};

/// \brief Utility class that reorder a linear system as needed for e.g. parallel ILU
///
/// The interior rows, i.e. those that do not have column indices that belong
/// to the overlap come first, then the rows that belong to the overlap but are
/// in the owner region of a process with a lower label, and then the rest.
class InteriorInterfacePermutation
    : public LinearSystemPermutation
{
public:
    template<class M, class C>
    InteriorInterfacePermutation(const M& A, const C& comm)
        : LinearSystemPermutation(A.N()), process_mapping_(comm.communicator())
    {
        computeRowPermutation(A, comm);
    }

    const std::array<std::size_t, 2>& interiorInterval() const
    {
        return interior_interval_;
    }

    const std::array<std::size_t, 2>& interfaceInterval() const
    {
        return interface_interval_;
    }

    const Detail::IdProcessMap& processMapping() const
    {
        return process_mapping_;
    }

    template<class Matrix, class ParallelIndexSet>
    void
    exchangeMatrixRows(Matrix const& A,
                       std::vector<typename Matrix::size_type> const& inverse_permutation,
                       ParallelIndexSet const& indexSet, SendReceiveCommunicator& comm)
    {
        typedef Dune::GlobalLookupIndexSet<ParallelIndexSet> IndexSet;
        typedef Detail::SparsityPatternHandle<Matrix, IndexSet> Handle;

        IndexSet global_indexSet(indexSet, A.N());
        Handle handle(A, additional_nonzeros_, global_indexSet, *this,
                      inverse_permutation);
        no_additional_entries_ = 0;

        comm.receiveData(handle);
        comm.sendData(handle);

        for(auto& row: additional_nonzeros_)
        {
            // \todo investigate why we need to process A, too.
            const auto& original_row =
                A[ inverse_permutation[ row.first ] ];

            for(auto col = original_row.begin(), endCol = original_row.end();
                col != endCol; ++col)
            {
                row.second.insert(row_permutation_[col.index()]);
            }
            no_additional_entries_ += row.second.size() -
                original_row.size();
        }
    }

private:
    template<class M, class C>
    void computeRowPermutation(const M& A, const C& comm)
    {
        std::size_t interior_row_counter = 0;
        std::size_t interface_row_counter = 0;
        Dune::GlobalLookupIndexSet<typename C::ParallelIndexSet>
            lookup_indices(comm.indexSet(), A.N());

        for ( auto row=A.begin(), end_row=A.end(); row != end_row; ++row )
        {
            auto index_pair_ptr = lookup_indices.pair(row.index());

            if ( !index_pair_ptr ||
                 index_pair_ptr->local().attribute() == Dune::OwnerOverlapCopyAttributeSet::owner )
            {
                bool is_interface_row = false;
                for ( auto col = row->begin(), end_col = row->end(); col != end_col;
                      ++col)
                {
                    index_pair_ptr = lookup_indices.pair(col.index());
                    if ( index_pair_ptr &&
                         index_pair_ptr->local().attribute() != Dune::OwnerOverlapCopyAttributeSet::owner )
                    {
                        is_interface_row = true;
                        break;
                    }
                }
                if ( ! is_interface_row )
                {
                    row_permutation_[row.index()] = interior_row_counter++;
                }
                else
                {
                    // interface row numbers are  in the range [2 *A.N(), 2* A.n()+1, ...
                    // This to make sure that the come after interior rows, and
                    // rows in the overlap that are owned by processes with a lower
                    // label than ours.
                    // They will be renumbered afterwards to achieve consecutive
                    // numbers.
                    row_permutation_[row.index()] = 2 * A.N() + interface_row_counter++;
                }
            }
        }
        // Now we need to examine the ownership of the rows associated with the
        // overlap. Those owned by a process with a lower label than ours come
        // between the interior and interface rows, the others after the interface rows.
        // 1. Assign lower numbers to nodes owned by lower procs
        //    range is [A.N(), A.N()+1, ... A.N() + no_lower_overlap_rows)
        //    I.e. the come after the interior but before the interface.
        std::size_t no_lower_overlap_rows =
            computeOverlapRowPermutation(A.N(), comm, std::less<int>());

        // 2. Assign higher numbers to rows owned by higher procs
        //    range is [3*A.N(), 3*A.N()+1, ...
        //    I.e. they are the last rows.
        std::map<std::size_t, std::pair<std::size_t, std::size_t> > backward_sizes;
        computeOverlapRowPermutation(3*A.N(), comm, std::greater<int>());

        // Close the gaps by moving values in [A.N(), A.N()+no_lower_overlap_rows ) to
        // [interior_rows, interior_rows + no_lower_overlap_rows), and
        // [2*A.N(), 2*A.N() + no_interface_rows) to
        // [interior_rows + no_overlap_rows, interior_rows + no_overlap_rows +interface_rows)
        // ansd so on.

        // offset for consecutive numbering
        std::array<std::size_t, 4> offset;
        offset[0] = 0;                                   // interior first
        offset[1] = interior_row_counter;                // then lower overlap
        offset[2] = offset[ 1 ] + no_lower_overlap_rows; // then interface
        offset[3] = offset[ 2 ] + interface_row_counter;        // and finally the rest.

        for(auto& permuted: row_permutation_)
        {
            permuted = offset[ permuted / A.N() ] + permuted % A.N();
        }
        // For convenience we save the intervals that we compute for.
        interior_interval_  =  { offset[0], offset[1] };
        interface_interval_ =  { offset[2], offset[3] };
        if ( interior_interval_[0] == interior_interval_[1] )
        {
            DUNE_THROW(Dune::RangeError, "Subdomain is too small!");
        }
    }

    template<class T, class C>
    std::size_t computeOverlapRowPermutation(std::size_t offset,
                                             const C& comm,
                                             const T& compare)
    {
        std::size_t counter = offset;
        for( const auto& proc_remote_lists : comm.remoteIndices() )
        {
            if ( process_mapping_
                 .compareWithOtherLabel( proc_remote_lists.first, compare ) )
            {
                assert(proc_remote_lists.second.first == proc_remote_lists.second.second );
                std::size_t recv_size = 0;
                std::size_t send_size = 0;

                const auto& rlist = *(proc_remote_lists.second.first);
                for ( const auto& remote_index : rlist )
                {
                    if ( remote_index.attribute() == Dune::OwnerOverlapCopyAttributeSet::owner)
                    {
                        if ( row_permutation_[remote_index.localIndexPair().local()] ==
                             std::numeric_limits<std::size_t>::max() )
                        {
                            row_permutation_[remote_index.localIndexPair().local()] = counter++;
                        }
                        ++recv_size;
                    }
                    else
                    {
                        if( remote_index.localIndexPair().local().attribute() == Dune::OwnerOverlapCopyAttributeSet::owner )
                        {
                            ++send_size;
                        }
                    }
                }
            }
        }
        return counter-offset;
    }
    // The mapping of the process onto an order
    Detail::IdProcessMap process_mapping_;
    /// \brief Interval in which the indices of the interior are.
    std::array<std::size_t, 2> interior_interval_;
    /// \brief Interval in which the indices of the interface are.
    std::array<std::size_t, 2> interface_interval_;
};

/// \brief Utility class describing the communication Interface.
class InterfaceAndCommunicator
{
public:
    typedef Dune::Interface::InformationMap InformationMap;
    typedef InformationMap::mapped_type::first_type IndexList;

    InterfaceAndCommunicator(MPI_Comm comm)
        : comm_(comm)
    {}

    typedef std::map<std::size_t, std::pair<std::size_t, std::size_t> > SizesMap;
    typedef InformationMap::iterator iterator;

    void reserve(const SizesMap& sizes_map)
    {
        for ( const auto& sizes : sizes_map )
        {
            auto& interface_pair = interfaces_[sizes.first];
            interface_pair.first.reserve(std::max(1ul, sizes.second.first));
            interface_pair.second.reserve(std::max(1ul, sizes.second.second));
        }
    }

    ~InterfaceAndCommunicator()
    {
        for( auto& entry : interfaces_)
        {
            entry.second.first.free();
            entry.second.second.free();
        }
    }

    iterator find(std::size_t proc)
    {
        return interfaces_.find(proc);
    }

    MPI_Comm communicator() const
    {
        return comm_;
    }

    const InformationMap& interfaces() const
    {
        return interfaces_;
    }

    InformationMap& interfaces()
    {
        return interfaces_;
    }
private:
    MPI_Comm comm_;
    InformationMap interfaces_;
};

/// \brief Utility class for sending values to higher/lower labelled processors.
template<class Communicator>
class ForwardBackwardLinearSystemCommunicator
{
public:
    ForwardBackwardLinearSystemCommunicator(const Communicator& communicator,
                                            const InteriorInterfacePermutation& permutation)
        : forward_interface_(communicator.communicator()),
          backward_interface_(communicator.communicator())
    {
        calculateInterfaceSizes(communicator.remoteIndices(), permutation);
        calculateInterfaceEntries(communicator.remoteIndices(), permutation);
        forward_communicator_
            .reset(new SendReceiveCommunicator(communicator.communicator(),
                                               forward_interface_.interfaces()));
        backward_communicator_
            .reset(new SendReceiveCommunicator(communicator.communicator(),
                                               backward_interface_.interfaces()));

    }

    SendReceiveCommunicator& getForwardCommunicator()
    {
        return *forward_communicator_;
    }

    SendReceiveCommunicator& getBackwardCommunicator()
    {
        return *backward_communicator_;
    }


private:

    template<class RemoteIndices>
    void calculateInterfaceSizes(const RemoteIndices& remote_indices,
                                 const InteriorInterfacePermutation& permutation)
    {
        typedef typename RemoteIndices::RemoteIndex RemoteIndex;
        std::map<std::size_t, std::pair<std::size_t, std::size_t> > forward_sizes;
        std::map<std::size_t, std::pair<std::size_t, std::size_t> > backward_sizes;
        auto size_adder =
            [](std::size_t& size, const RemoteIndex&)
            {
                ++size;
            };
        processRemoteIndices(remote_indices, size_adder,
                             permutation.processMapping(),
                             forward_sizes, backward_sizes);
        // Reserve interface space
        forward_interface_.reserve(forward_sizes);
        backward_interface_.reserve(backward_sizes);
    }

    template<class RemoteIndices>
    void calculateInterfaceEntries(const RemoteIndices& remote_indices,
                                   const InteriorInterfacePermutation& permutation)
    {
        typedef typename InterfaceAndCommunicator::IndexList IndexList;
        typedef typename RemoteIndices::RemoteIndex RemoteIndex;
        auto index_adder =
            [&permutation] (IndexList& list, const RemoteIndex& index)
            {
                list.add(permutation[index.localIndexPair().local()]);
            };

        processRemoteIndices(remote_indices, index_adder,
                             permutation.processMapping(),
                             forward_interface_.interfaces(),
                             backward_interface_.interfaces());
    }

    template<class RemoteIndices, class IndexProcessor, class ProcessMap, class Interface>
    void processRemoteIndices(const RemoteIndices& remote_indices, const IndexProcessor& index_processor,
                              const ProcessMap& process_mapping,
                              Interface& forward_interface, Interface& backward_interface)
    {
        for( const auto& proc_remote_lists : remote_indices)
        {
            auto& other_proc = proc_remote_lists.first;

            if ( process_mapping
                 .compareWithOtherLabel(other_proc, std::less<std::size_t>()) )
            {
                // other label is smaller than ours
                processRemoteIndexList(*proc_remote_lists.second.first,
                                       forward_interface[other_proc],
                                       backward_interface[other_proc],
                                       index_processor);
            }
            else
            {
                processRemoteIndexList(*proc_remote_lists.second.first,
                                       backward_interface[other_proc],
                                       forward_interface[other_proc],
                                       index_processor);
            }
        }
    }

    template<class RemoteIndices, class InterfacePair, class IndexProcessor>
    void processRemoteIndexList(const RemoteIndices& index_list, InterfacePair& receive_interface_pair,
                                InterfacePair& send_interface_pair, IndexProcessor& index_processor)
    {
        // For each communication there are two lists: the first with indices
        // we send from and the second one with indices we receive data.
        auto& receive_interface = receive_interface_pair.second;
        auto& send_interface    = send_interface_pair.first;

        for ( const auto& remote_index : index_list )
        {
            if ( remote_index.attribute() == Dune::OwnerOverlapCopyAttributeSet::owner )
            {
                index_processor(receive_interface, remote_index);
            }
            else
            {
                if ( remote_index.localIndexPair().local().attribute() ==
                     Dune::OwnerOverlapCopyAttributeSet::owner )
                {
                    index_processor(send_interface, remote_index);
                }
            }
        }
    }

    InterfaceAndCommunicator forward_interface_;
    InterfaceAndCommunicator backward_interface_;
    std::unique_ptr<SendReceiveCommunicator> forward_communicator_;
    std::unique_ptr<SendReceiveCommunicator> backward_communicator_;
};

} // end namespace Detail


template<class Matrix, class Domain, class Range, typename ParallelInfo>
class ParallelILU0 : public Dune::Preconditioner<Domain,Range> {
public:
    //! \brief The matrix type the preconditioner is for.
    typedef typename Dune::remove_const<Matrix>::type matrix_type;
    //! \brief The domain type of the preconditioner.
    typedef Domain domain_type;
    //! \brief The range type of the preconditioner.
    typedef Range range_type;
    //! \brief The field type of the preconditioner.
    typedef typename Domain::field_type field_type;

    // define the category
    enum {
        //! \brief The category the preconditioner is part of.
        category=Dune::SolverCategory::overlapping
    };

    /*! \brief Constructor.

      Constructor gets all parameters to operate the prec.
      \param A The matrix to operate on.
      \param C The parallel information.
      \param w The relaxation factor.
    */
    ParallelILU0 (const Matrix& A, const ParallelInfo& comm, field_type w)
        : ilu_(), comm_(comm),
          row_permutation_(A, comm),
          forward_backward_communicator_(comm,
                                         row_permutation_),
          w_(w)
    {
        std::vector<std::size_t> inverse_row_permutation =
            row_permutation_.createInverseRowPermutation();

        row_permutation_.exchangeMatrixRows(A,
            inverse_row_permutation,
            comm.indexSet(),
            forward_backward_communicator_.getForwardCommunicator());
        ilu_ = row_permutation_.createPermutedMatrix(A,
                                                     row_permutation_.interiorInterval()[1],
                                                     row_permutation_.interfaceInterval()[0],
                                                     inverse_row_permutation);
        int ilu_setup_successful = decompose();
        // Check whether there was a problem on some process
        if ( comm.communicator().min(ilu_setup_successful) == 0 )
        {
            throw Dune::MatrixBlockError();
        }
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
        Range permuted_d(d.size());
        Domain permuted_v(v.size());
        const auto& interior_interval = row_permutation_.interiorInterval();
        const auto& interface_interval = row_permutation_.interfaceInterval();
        row_permutation_.permutateOrder(d, permuted_d);
        forwardSolveRowInterval(permuted_v, permuted_d,
                                interior_interval[0], interior_interval[1]);
        receiveValuesFromLowerProcs(permuted_v);
        forwardSolveRowInterval(permuted_v, permuted_d,
                                interface_interval[0], interface_interval[1]);
        sendValuesToHigherProcs(permuted_v);
        receiveValuesFromHigherProcs(permuted_v);
        backwardSolveRowInterval(permuted_v,
                                 interface_interval[0], interface_interval[1]);
        sendValuesToLowerProcs(permuted_v);
        backwardSolveRowInterval(permuted_v,
                                 interior_interval[0], interior_interval[1]);
        row_permutation_.permutateOrderBackwardsAndScale(v, permuted_v, w_);
        // v is not yet consistent as the overlap/copy region shared owner
        // processes with a lower label is not up to date.
        comm_.copyOwnerToAll(v,v); // We send more than we need here.
        // \todo use asynchronous sends to higher / receives from lower procs.
    }

    /*!
      \brief Clean up.

      \copydoc Preconditioner::post(X&)
    */
    virtual void post (Domain& x)
    {
        DUNE_UNUSED_PARAMETER(x);
    }

private:
    typedef typename ParallelInfo::ParallelIndexSet IndexSet;
    typedef Detail::PermutedGlobalLookupIndexSet<IndexSet>
    PermutedGlobalLookupIndexSet;

    void receiveRowsFromLowerProcs(const PermutedGlobalLookupIndexSet& indexset)
    {
        Detail::PermutedSparseMatrixHandle<Matrix, IndexSet>
            handle(*ilu_, indexset, row_permutation_);
        auto& communicator =
            forward_backward_communicator_.getForwardCommunicator();
        communicator.receiveData(handle);
    }

    void sendRowsToHigherProcs(const PermutedGlobalLookupIndexSet& indexset)
    {
        Detail::PermutedSparseMatrixHandle<Matrix, IndexSet>
            handle(*ilu_, indexset, row_permutation_);
        auto& communicator =
            forward_backward_communicator_.getForwardCommunicator();
        communicator.sendData(handle);
    }

    struct VectorHandle
    {
        typedef typename Domain::block_type DataType;

        VectorHandle(Domain& x)
            : x_(x)
        {}
        bool fixedsize() const
        {
            return true;
        }
        std::size_t size(std::size_t) const
        {
            return 1;
        }
        template<class B>
        void gather(B& buffer, std::size_t i) const
        {
            buffer.write(x_[i]);
        }
        template<class B>
        void scatter(B& buffer, std::size_t i, std::size_t)
        {
            buffer.read(x_[i]);
        }
    private:
        Domain& x_;
    };

    void receiveValuesFromLowerProcs(Domain& v)
    {
        VectorHandle handle(v);
        auto& communicator =
            forward_backward_communicator_.getForwardCommunicator();
        communicator.receiveData(handle);
    }

    void sendValuesToHigherProcs(Domain& v)
    {
        VectorHandle handle(v);
        auto& communicator =
            forward_backward_communicator_.getForwardCommunicator();
        communicator.sendData(handle);
    }

    void receiveValuesFromHigherProcs(Domain& v)
    {
        VectorHandle handle(v);
        auto& communicator =
            forward_backward_communicator_.getBackwardCommunicator();
        communicator.receiveData(handle);
    }

    void sendValuesToLowerProcs(Domain& v)
    {
        VectorHandle handle(v);
        auto& communicator =
            forward_backward_communicator_.getBackwardCommunicator();
        communicator.sendData(handle);
    }


    int decompose()
    {
        // decompose interior rows.
        const auto& interior_interval = row_permutation_.interiorInterval();
        bool success = decomposeRowInterval(interior_interval[0], interior_interval[1]);
        PermutedGlobalLookupIndexSet indexset(comm_.indexSet(),
                                              row_permutation_.size(),
                                              row_permutation_);
        receiveRowsFromLowerProcs(indexset);
        // decompose interface rows
        const auto& interface_interval = row_permutation_.interfaceInterval();
        success = success &&
            decomposeRowInterval(interface_interval[0],
                                 interface_interval[1]);
        sendRowsToHigherProcs(indexset);
        return success;
    }

    bool decomposeRowInterval(std::size_t start, std::size_t end)
    {
        Detail::SubMatrix<Matrix> interior_sub_matrix(*ilu_, ilu_->begin() + start,
                                                 ilu_->begin() + end);
        try
        {
            bilu0_decomposition(interior_sub_matrix);
        }
        catch ( Dune::MatrixBlockError error )
        {
            std::string message = error.what();
            std::cerr<<"Exception occured on process " <<
                comm_.communicator().rank() << " during " <<
                "setup of ILU0 preconditioner with message: " <<
                message<<std::endl;
            return false;
        }
        return true;
    }

    void forwardSolveRowInterval(Domain& v, const Range& d,
                                 std::size_t start, std::size_t end)
    {
        assert(start<=end);
        auto end_row = ilu_->begin() + end;
        for ( auto row = ilu_->begin() + start; row != end_row; ++row )
        {
            auto rhs(d[row.index()]);
            for ( auto col = row->begin(); col.index() < row.index(); ++col )
            {
                col->mmv(v[col.index()],rhs);
            }
            v[row.index()] = rhs;
        }
    }

    void backwardSolveRowInterval(Domain& v,
                                  std::size_t start, std::size_t end)
    {
        assert(start<=end);
        auto end_row = ilu_->beforeEnd() - (ilu_->N() - start);
        for ( auto row = ilu_->beforeEnd() - (ilu_->N() - end);
              row != end_row; --row)
        {
            auto rhs(v[row.index()]);
            auto col = row->beforeEnd();
            for( ; col.index() > row.index(); --col)
                col->mmv(v[col.index()], rhs);
            v[row.index()] = 0;
            col->umv(rhs, v[row.index()]);
        }
    }

    //! \brief The ILU0 decomposition of the matrix.
    std::unique_ptr<matrix_type> ilu_;
    const ParallelInfo& comm_;
    //! \brief The indices of the inner rows.
    Detail::InteriorInterfacePermutation row_permutation_;
    //! \brief The communicator for sending and receiving values
    Detail::ForwardBackwardLinearSystemCommunicator<ParallelInfo> forward_backward_communicator_;
    //! \brief The relaxation factor to use.
    field_type w_;
};
} // end namespace Opm
#endif //HAVE_MPI
#endif
