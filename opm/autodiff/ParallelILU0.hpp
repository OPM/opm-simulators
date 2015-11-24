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
class IdProcessMap
{
public:
    template<class C>
    explicit IdProcessMap(const C& comm)
        : rank_(comm.rank())
    {}
    std::size_t myMapping() const
    {
        return rank_;
    }
    template<class T>
    bool compareWithOtherLabel(std::size_t other_proc,
                               const T& compare) const
    {
        return compare(other_proc, myMapping());
    }
private:
    const int rank_;
};

/**
 * \brief A simple Wrapper of a BCRSMatrix that just makes modified row iterators available.
 *
 * This is needed for the parallel ILU implementation which needs to neglect
 * rows at the begin and/or end of the matrix.
 */
template<class M>
class SubMatrix
{
public:
    typedef typename M::RowIterator RowIterator;
    typedef typename M::ColIterator ColIterator;
    typedef typename M::block_type block_type;
    typedef typename M::row_type row_type;
    typedef typename M::size_type size_type;

    SubMatrix(M& mat, const RowIterator& begin, const RowIterator& end)
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
    M&          mat_;
    RowIterator begin_;
    RowIterator end_;
};

class InteriorInterfacePermutation
{
public:
    template<class M, class C>
    InteriorInterfacePermutation(const M& A, const C& comm)
        : process_mapping_(comm.communicator())
    {
        computeRowPermutation(A, comm);
    }

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

    template<class B>
    std::unique_ptr<Dune::BCRSMatrix<B> > createPermutedMatrix(const Dune::BCRSMatrix<B>& A)
    {
        std::unique_ptr<Dune::BCRSMatrix<B> > permuted_ptr(new Dune::BCRSMatrix<B>());
        Dune::BCRSMatrix<B>& permuted = *permuted_ptr;

        Dune::ImplicitMatrixBuilder<Dune::BCRSMatrix<B> >
            builder(permuted, A.N(), A.M(),
                    std::ceil(static_cast<double>(A.nonzeroes())/A.M()), .1);
        for( std::size_t index = 0; index < A.N(); ++index)
        {
            auto permuted_row = builder[ row_permutation_[index] ];
            const auto& original_row = A[index];
            for(auto col = original_row.begin(), endCol = original_row.end();
                col != endCol; ++col)
            {
                permuted_row[ row_permutation_[col.index()] ] = *col;
            }
        }
        permuted.compress();
        return permuted_ptr;
    }

    std::size_t operator[](std::size_t i) const
    {
        return row_permutation_[i];
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

    std::size_t size() const
    {
        return row_permutation_.size();
    }

private:
    template<class M, class C>
    void computeRowPermutation(const M& A, const C& comm)
    {
        row_permutation_.resize(A.N(), std::numeric_limits<std::size_t>::max());
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
    /// \brief The permuted indices
    std::vector<std::size_t> row_permutation_;
    /// \brief Interval in which the indices of the interior are.
    std::array<std::size_t, 2> interior_interval_;
    /// \brief Interval in which the indices of the interface are.
    std::array<std::size_t, 2> interface_interval_;
};


template<class I>
class PermutedGlobalLookupIndexSet
{
public:
    typedef typename I::GlobalIndex GlobalIndex;
    typedef Dune::IndexPair<typename I::GlobalIndex, typename I::LocalIndex> IndexPair;

    PermutedGlobalLookupIndexSet(const I& indexSet,
                                 const InteriorInterfacePermutation& permutation)
        : indexSet_(indexSet), index_pairs_(permutation.size(), nullptr)
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

    const I& indexSet() const
    {
        return indexSet_;
    }
private:
    const I& indexSet_;
    std::vector<const IndexPair*> index_pairs_;
};



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

template<class C>
class ForwardBackwardLinearSystemCommunicator
{
public:
    ForwardBackwardLinearSystemCommunicator(const C& communicator,
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

    /**
       template<class Container>
       void buildForwardBuffers(const Container& container)
       {
       forward_communicator_.free();
       backward_communicator_.free();
       forward_communicator_.build(container, container, forward_interface_);
       }

       template<class Container>
       void buildBuffers(const Container& container)
       {
       forward_communicator_.build(container, container, forward_interface_);
       backward_communicator_.build(container, container, backward_interface_);
       }

       void freeBuffers()
       {
       forward_communicator_.free();
       backward_communicator_.free();
       }
    */
private:

    template<class R>
    void calculateInterfaceSizes(const R& remote_indices,
                                 const InteriorInterfacePermutation& permutation)
    {
        typedef typename R::RemoteIndex RemoteIndex;
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

    template<class R>
    void calculateInterfaceEntries(const R& remote_indices,
                                   const InteriorInterfacePermutation& permutation)
    {
        typedef typename InterfaceAndCommunicator::IndexList IndexList;
        typedef typename R::RemoteIndex RemoteIndex;
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

    template<class R, class I, class P, class T>
    void processRemoteIndices(const R& remote_indices, const I& index_processor,
                              const P& process_mapping,
                              T& forward_interface, T& backward_interface)
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

    template<class R, class I, class T>
    void processRemoteIndexList(const R& index_list, I& receive_interface_pair,
                                I& send_interface_pair, T& index_processor)
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


template<class M, class X, class Y, typename C>
class ParallelILU0 : public Dune::Preconditioner<X,Y> {
public:
    //! \brief The matrix type the preconditioner is for.
    typedef typename Dune::remove_const<M>::type matrix_type;
    //! \brief The domain type of the preconditioner.
    typedef X domain_type;
    //! \brief The range type of the preconditioner.
    typedef Y range_type;
    //! \brief The field type of the preconditioner.
    typedef typename X::field_type field_type;

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
    ParallelILU0 (const M& A, const C& comm, field_type w)
        : ilu_(), comm_(comm),
          row_permutation_(A, comm),
          forward_backward_communicator_(comm,
                                         row_permutation_),
          w_(w)
    {
        ilu_ = row_permutation_.createPermutedMatrix(A);
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
    virtual void pre (X& x, Y& b)
    {
        DUNE_UNUSED_PARAMETER(x);
        DUNE_UNUSED_PARAMETER(b);
    }

    /*!
      \brief Apply the preconditoner.

      \copydoc Preconditioner::apply(X&,const Y&)
    */
    virtual void apply (X& v, const Y& d)
    {
        Y permuted_d(d.size());
        X permuted_v(v.size());
        const auto& interior_interval = row_permutation_.interiorInterval();
        const auto& interface_interval = row_permutation_.interfaceInterval();
        row_permutation_.permutateOrder(d, permuted_d);
        forwardSolveRowInterval(permuted_v, permuted_d,
                                interior_interval[0], interior_interval[0]);
        receiveValuesFromLowerProcs(v);
        forwardSolveRowInterval(permuted_v, permuted_d,
                                interface_interval[0], interface_interval[1]);
        sendValuesToHigherProcs(permuted_v);
        receiveValuesFromHigherProcs(permuted_v);
        backwardSolveRowInterval(permuted_v,
                                 interface_interval[0], interface_interval[1]);
        sendValuesToLowerProcs(permuted_v);
        backwardSolveRowInterval(permuted_v,
                                 interior_interval[0], interior_interval[0]);
        row_permutation_.permutateOrderBackwardsAndScale(v, permuted_v, w_);
    }

    /*!
      \brief Clean up.

      \copydoc Preconditioner::post(X&)
    */
    virtual void post (X& x)
    {
        DUNE_UNUSED_PARAMETER(x);
    }

private:
    typedef Detail::PermutedGlobalLookupIndexSet<typename C::ParallelIndexSet>
    PermutedGlobalLookupIndexSet;

    struct PermutedSparseMatrixHandle
    {
        typedef typename PermutedGlobalLookupIndexSet::GlobalIndex GlobalIndex;
        typedef typename M::block_type block_type;
        typedef std::pair<block_type, GlobalIndex> DataType;

        PermutedSparseMatrixHandle(M& A, const PermutedGlobalLookupIndexSet& indexset,
                                   const  Detail::InteriorInterfacePermutation& row_permutation)
            : A_(A), indexset_(indexset), row_permutation_(row_permutation)
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
                typedef typename PermutedGlobalLookupIndexSet::GlobalIndex
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
            typedef typename PermutedGlobalLookupIndexSet::GlobalIndex GlobalIndex;
            GlobalIndex index = std::numeric_limits<GlobalIndex>::max();
            while ( data.back().second == std::numeric_limits<GlobalIndex>::max() )
            {
                data.pop_back();
            }

            auto hint = indexset_.indexSet().begin();
            for( auto& datum : data)
            {
                typedef typename C::ParallelIndexSet::IndexPair IndexPair;
                // \todo This search for every element seems rather expensive.
                // We should try to store the receive values as they are.
                // and later do only one search per column index.
                auto found =
                    std::lower_bound(hint, indexset_.indexSet().end(), datum,
                                     [](const IndexPair& d1, const DataType& d2)
                                     { return d1.global() < d2.second; });
                if ( found != indexset_.indexSet().end() )
                {
                    hint = found; // we have found->global() <= datum.second !
                    if( found->global() == datum.second )
                    {
                        // The index set knows nothing about the permutation.
                        // therefore we need to apply the permutation to the
                        // local index.
                        row[ row_permutation_[found->local()] ] = datum.first;
                    }
                    else
                    {
                        // We have discarded all matrix row entries to unknows
                        // that are neither in our owner/interior or our
                        // ghost/copy region. Therefore data contains all
                        // connection/columns that are in the region we store.
                        // If we cannot find a global index than this is an error
                        DUNE_THROW(Dune::RangeError, "Process is missing a"
                                   <<" connection from local index "<<i
                                   <<" to global index "<<datum.second
                                   <<". Please check the matrix or index set.");
                    }
                }
            }
        }
    private:
        M& A_;
        const PermutedGlobalLookupIndexSet& indexset_;
        const Detail::InteriorInterfacePermutation& row_permutation_;
        int rank_;
    };

    void receiveRowsFromLowerProcs(const PermutedGlobalLookupIndexSet& indexset)
    {
        PermutedSparseMatrixHandle handle(*ilu_, indexset, row_permutation_);
        auto& communicator =
            forward_backward_communicator_.getForwardCommunicator();
        communicator.receiveData(handle);
    }

    void sendRowsToHigherProcs(const PermutedGlobalLookupIndexSet& indexset)
    {
        PermutedSparseMatrixHandle handle(*ilu_, indexset, row_permutation_);
        auto& communicator =
            forward_backward_communicator_.getForwardCommunicator();
        communicator.sendData(handle);
    }

    struct VectorHandle
    {
        typedef typename X::block_type DataType;

        VectorHandle(X& x)
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
        X& x_;
    };

    void receiveValuesFromLowerProcs(X& v)
    {
        VectorHandle handle(v);
        auto& communicator =
            forward_backward_communicator_.getForwardCommunicator();
        communicator.receiveData(handle);
    }

    void sendValuesToHigherProcs(X& v)
    {
        VectorHandle handle(v);
        auto& communicator =
            forward_backward_communicator_.getForwardCommunicator();
        communicator.sendData(handle);
    }

    void receiveValuesFromHigherProcs(X& v)
    {
        VectorHandle handle(v);
        auto& communicator =
            forward_backward_communicator_.getBackwardCommunicator();
        communicator.receiveData(handle);
    }

    void sendValuesToLowerProcs(X& v)
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
        Detail::SubMatrix<M> interior_sub_matrix(*ilu_, ilu_->begin() + start,
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
    void forwardSolveRowInterval(X& v, const Y& d,
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
    void backwardSolveRowInterval(X& v,
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
    const C& comm_;
    //! \brief The indices of the inner rows.
    Detail::InteriorInterfacePermutation row_permutation_;

    Detail::ForwardBackwardLinearSystemCommunicator<C> forward_backward_communicator_;
    //! \brief The relaxation factor to use.
    field_type w_;
};
} // end namespace Opm
#endif //HAVE_MPI
#endif 
