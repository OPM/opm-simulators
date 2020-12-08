/*
  Copyright 2020 OPM-OP AS

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
#ifndef OPM_PARALLELWELLINFO_HEADER_INCLUDED

#define OPM_PARALLELWELLINFO_HEADER_INCLUDED
#include <dune/common/version.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/parallel/plocalindex.hh>
#include <dune/istl/owneroverlapcopy.hh>

#include <opm/common/ErrorMacros.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Well/Well.hpp>

#include <memory>
#include <iterator>
#include <numeric>

namespace Opm
{

/// \brief Class to facilitate getting values associated with the above/below perforation
///
class CommunicateAboveBelow
{
    enum Attribute {
      owner=1, overlap=2
    };
    using MPIComm = typename Dune::MPIHelper::MPICommunicator;
#if DUNE_VERSION_NEWER(DUNE_COMMON, 2, 7)
    using Communication = Dune::Communication<MPIComm>;
#else
    using Communication = Dune::CollectiveCommunication<MPIComm>;
#endif
    using LocalIndex = Dune::ParallelLocalIndex<Attribute>;
#if HAVE_MPI
    using IndexSet = Dune::ParallelIndexSet<std::size_t,LocalIndex,50>;
    using RI = Dune::RemoteIndices<IndexSet>;
#endif

public:
    explicit CommunicateAboveBelow(const Communication& comm);
    /// \brief Adds information about original index of the perforations in ECL Schedule.
    ///
    /// \warning Theses indices need to be push in the same order as they
    ///          appear in the ECL well specifiation. Use -1 if there is
    ///          no perforation above.
    /// \param above The ECL index of the next open perforation above.
    /// \param current The ECL index of the current open perforation.
    void pushBackEclIndex(int above, int current, bool owner=true);

    /// \brief Clear all the parallel information
    void clear();

    /// \brief Indicates that we will add the index information
    /// \see pushBackEclIndex
    void beginReset();

    /// \brief Indicates that the index information is complete.
    ///
    /// Sets up the commmunication structures to be used by
    /// communicate()
    void endReset();

    /// \brief Creates an array of values for the perforation above.
    /// \param first_value Value to use for above of the first perforation
    /// \param current C-array of the values at the perforations
    /// \param size The size of the C-array and the returned vector
    /// \return a vector containing the values for the perforation above.
    std::vector<double> communicateAbove(double first_value,
                                         const double* current,
                                         std::size_t size);

    /// \brief Creates an array of values for the perforation below.
    /// \param first_value Value to use for above of the first perforation
    /// \param current C-array of the values at the perforations
    /// \param size The size of the C-array and the returned vector
    /// \return a vector containing the values for the perforation above.
    std::vector<double> communicateBelow(double first_value,
                                         const double* current,
                                         std::size_t size);

    /// \brief Do a (in place) partial sum on values attached to all perforations.
    ///
    /// For distributed wells this may include perforations stored elsewhere.
    /// The result is stored in ther range given as the parameters
    /// \param begin The start of the range
    /// \param ebd The end of the range
    /// \tparam RAIterator The type og random access iterator
    template<class RAIterator>
    void partialSumPerfValues(RAIterator begin, RAIterator end) const
    {
        if (this->comm_.size() < 2)
        {
            std::partial_sum(begin, end, begin);
        }
        else
        {
#if HAVE_MPI
            // The global index used in the index set current_indices
            // is the index of the perforation in ECL Schedule definition.
            // This is assumed to give the topological order that is used
            // when doing the partial sum.
            // allgather the index of the perforation in ECL schedule and the value.
            using Value = typename std::iterator_traits<RAIterator>::value_type;
            std::vector<int> sizes(comm_.size());
            std::vector<int> displ(comm_.size(), 0);
            using GlobalIndex = typename IndexSet::IndexPair::GlobalIndex;
            using Pair = std::pair<GlobalIndex,Value>;
            std::vector<Pair> my_pairs;
            my_pairs.reserve(current_indices_.size());
            for (const auto& pair: current_indices_)
            {
                if (pair.local().attribute() == owner)
                {
                    my_pairs.emplace_back(pair.global(), *(begin+(static_cast<const int&>(pair.local()))));
                }
            }
            int mySize = my_pairs.size();
            comm_.allgather(&mySize, 1, sizes.data());
            std::partial_sum(sizes.begin(), sizes.end(), displ.begin()+1);
            std::vector<Pair> global_pairs(displ.back());
            comm_.allgatherv(my_pairs.data(), my_pairs.size(), global_pairs.data(), sizes.data(), displ.data());
            // sort the complete range to get the correct ordering
            std::sort(global_pairs.begin(), global_pairs.end(),
                      [](const Pair& p1, const Pair& p2){ return p1.first < p2.first; } );
            std::vector<Value> sums(global_pairs.size());
            std::transform(global_pairs.begin(), global_pairs.end(), sums.begin(),
                           [](const Pair& p) { return p.second; });
            std::partial_sum(sums.begin(), sums.end(),sums.begin());
            // assign the values (both ranges are sorted by the ecl index)
            auto global_pair = global_pairs.begin();
            for (const auto& pair: current_indices_)
            {
                global_pair = std::lower_bound(global_pair, global_pairs.end(),
                                               pair.global(),
                                               [](const Pair& val1, const GlobalIndex& val2)
                                               { return val1.first < val2; });
                assert(global_pair != global_pairs.end());
                assert(global_pair->first == pair.global());
                *(begin + static_cast<const int&>(pair.local())) = sums[global_pair - global_pairs.begin()];
            }
#else
            OPM_THROW(std::logic_error, "In a sequential run the size of the communicator should be 1!");
#endif
        }
    }

private:
    Communication comm_;
#if HAVE_MPI
    /// \brief Mapping of the local well index to ecl index
    IndexSet current_indices_;
    /// \brief Mapping of the above well index to ecl index
    IndexSet above_indices_;
    RI remote_indices_;
    Dune::Interface interface_;
    Dune::BufferedCommunicator communicator_;
#endif
    std::size_t count_{};
};


/// \brief Class encapsulating some information about parallel wells
///
/// e.g. It provides a communicator for well information
class ParallelWellInfo
{
public:
    using MPIComm = typename Dune::MPIHelper::MPICommunicator;
#if DUNE_VERSION_NEWER(DUNE_COMMON, 2, 7)
    using Communication = Dune::Communication<MPIComm>;
#else
    using Communication = Dune::CollectiveCommunication<MPIComm>;
#endif


    /// \brief Constructs object using MPI_COMM_SELF
    ParallelWellInfo(const std::string& name = {""},
                     bool hasLocalCells = true);

    /// \brief Constructs object with communication between all rank sharing
    ///        a well
    /// \param well_info Pair of well name and whether local cells might be perforated
    ///        on this rank
    /// \param allComm The communication object with all MPI ranks active in the simulation.
    ///                Default is the one with all ranks available.
    ParallelWellInfo(const std::pair<std::string,bool>& well_info,
                     Communication allComm = Communication());

    const Communication& communication() const
    {
        return *comm_;
    }

    /// \brief Collectively decide which rank has first perforation.
    void communicateFirstPerforation(bool hasFirst);

    template<class T>
    T broadcastFirstPerforationValue(const T& t) const
    {
        T res = t;
#ifndef NDEBUG
        assert(rankWithFirstPerf_ >= 0 && rankWithFirstPerf_ < comm_->size());
        // At least on some OpenMPI version this might broadcast might interfere
        // with other communication if there are bugs
        comm_->barrier();
#endif
        comm_->broadcast(&res, 1, rankWithFirstPerf_);
#ifndef NDEBUG
        comm_->barrier();
#endif
        return res;
    }

    /// \brief Creates an array of values for the perforation above.
    /// \param first_value Value to use for above of the first perforation
    /// \param current C-array of the values at the perforations
    /// \param size The size of the C-array and the returned vector
    /// \return a vector containing the values for the perforation above.
    std::vector<double> communicateAboveValues(double first_value,
                                               const double* current,
                                               std::size_t size) const;

    /// \brief Creates an array of values for the perforation above.
    /// \param first_value Value to use for above of the first perforation
    /// \param current vector of current values
    std::vector<double> communicateAboveValues(double first_value,
                                               const std::vector<double>& current) const;

    /// \brief Creates an array of values for the perforation below.
    /// \param last_value Value to use for below of the last perforation
    /// \param current C-array of the values at the perforations
    /// \param size The size of the C-array and the returned vector
    /// \return a vector containing the values for the perforation above.
    std::vector<double> communicateBelowValues(double last_value,
                                               const double* current,
                                               std::size_t size) const;

    /// \brief Creates an array of values for the perforation above.
    /// \param last_value Value to use for below of the last perforation
    /// \param current vector of current values
    std::vector<double> communicateBelowValues(double last_value,
                                               const std::vector<double>& current) const;

    /// \brief Adds information about the ecl indices of the perforations.
    ///
    /// \warning Theses indices need to be push in the same order as they
    ///          appear in the ECL well specifiation. Use -1 if there is
    ///          no perforation above.
    /// \param above The ECL index of the next open perforation above.
    /// \param current The ECL index of the current open perforation.
    void pushBackEclIndex(int above, int current);

    /// \brief Name of the well.
    const std::string& name() const
    {
        return name_;
    }

    /// \brief Whether local cells are perforated somewhen
    bool hasLocalCells() const
    {
        return hasLocalCells_;
    }
    bool isOwner() const
    {
        return isOwner_;
    }

    /// \brief Inidicate that we will reset the ecl index information
    ///
    /// \see pushBackEclIndex;
    void beginReset();

    /// \brief Inidicate completion of reset of the ecl index information
    void endReset();

    /// \brief Sum all the values of the perforations
    template<typename It>
    typename std::iterator_traits<It>::value_type sumPerfValues(It begin, It end) const
    {
        using V = typename std::iterator_traits<It>::value_type;
        /// \todo cater for overlap later. Currently only owner
        auto local = std::accumulate(begin, end, V());
        return communication().sum(local);
    }

    /// \brief Do a (in place) partial sum on values attached to all perforations.
    ///
    /// For distributed wells this may include perforations stored elsewhere.
    /// The result is stored in ther range given as the parameters
    /// \param begin The start of the range
    /// \param ebd The end of the range
    /// \tparam RAIterator The type og random access iterator
    template<class RAIterator>
    void partialSumPerfValues(RAIterator begin, RAIterator end) const
    {
        commAboveBelow_->partialSumPerfValues(begin, end);
    }

    /// \brief Free data of communication data structures.
    void clearCommunicateAboveBelow();
private:

    /// \brief Deleter that also frees custom MPI communicators
    struct DestroyComm
    {
        void operator()(Communication* comm);
    };


    /// \brief Name of the well.
    std::string name_;
    /// \brief Whether local cells are perforated somewhen
    bool hasLocalCells_;
    /// \brief Whether we own the well and should do reports etc.
    bool isOwner_;
    /// \brief Rank with the first perforation on it.
    int rankWithFirstPerf_;
    /// \brief Communication object for the well
    ///
    /// Contains only ranks where this well will perforate local cells.
    std::unique_ptr<Communication, DestroyComm> comm_;

    /// \brief used to communicate the values for the perforation above.
    std::unique_ptr<CommunicateAboveBelow> commAboveBelow_;
};

/// \brief Class checking that all connections are on active cells
///
/// Works for distributed wells, too
class CheckDistributedWellConnections
{
public:
    CheckDistributedWellConnections(const Well& well,
                                   const ParallelWellInfo& info);

    /// \brief Inidicate that the i-th completion was found
    ///
    /// in the local grid.
    /// \param index The index of the completion in Well::getConnections
    void connectionFound(std::size_t index);

    bool checkAllConnectionsFound();
private:
    std::vector<std::size_t> foundConnections_;
    const Well& well_;
    const ParallelWellInfo& pwinfo_;
};

bool operator<(const ParallelWellInfo& well1, const ParallelWellInfo& well2);

bool operator==(const ParallelWellInfo& well1, const ParallelWellInfo& well2);

bool operator!=(const ParallelWellInfo& well1, const ParallelWellInfo& well2);

bool operator<(const std::pair<std::string, bool>& pair, const ParallelWellInfo& well);

bool operator<( const ParallelWellInfo& well, const std::pair<std::string, bool>& pair);

bool operator==(const std::pair<std::string, bool>& pair, const ParallelWellInfo& well);

bool operator==(const ParallelWellInfo& well, const std::pair<std::string, bool>& pair);

bool operator!=(const std::pair<std::string, bool>& pair, const ParallelWellInfo& well);

bool operator!=(const ParallelWellInfo& well, const std::pair<std::string, bool>& pair);

} // end namespace Opm
#endif //  OPM_PARALLELWELLINFO_HEADER_INCLUDED
