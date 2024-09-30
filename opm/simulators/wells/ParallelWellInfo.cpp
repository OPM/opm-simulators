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
#include <config.h>
#include <opm/simulators/wells/ParallelWellInfo.hpp>

#include <opm/common/ErrorMacros.hpp>
#include <opm/input/eclipse/Schedule/Well/Well.hpp>
#include <opm/input/eclipse/Schedule/Well/WellConnections.hpp>

#include <cassert>
#include <iterator>
#include <numeric>

namespace Dune
{
#if HAVE_MPI
template<>
struct CommPolicy<double*>
{
    using Type = double*;
    using IndexedType = double;
    using IndexedTypeFlag = Dune::SizeOne;
    static const void* getAddress(const double*& v, int index)
    {
        return v + index;
    }
    static int getSize(const double*&, int)
    {
        return 1;
    }
};

template<>
struct CommPolicy<float*>
{
    using Type = float*;
    using IndexedType = float;
    using IndexedTypeFlag = Dune::SizeOne;
    static const void* getAddress(const float*& v, int index)
    {
        return v + index;
    }
    static int getSize(const float*&, int)
    {
        return 1;
    }
};
#endif
}

namespace Opm
{

template<class Scalar>
GlobalPerfContainerFactory<Scalar>::
GlobalPerfContainerFactory(const IndexSet& local_indices,
                           const Parallel::Communication comm,
                           const int num_local_perfs)
    : local_indices_(local_indices), comm_(comm)
{
    if ( comm_.size() > 1 )
    {
        // The global index used in the index set current_indices
        // is the index of the perforation in ECL Schedule definition.
        // This is assumed to give the topological order.
        // allgather the index of the perforation in ECL schedule and the value.
        sizes_.resize(comm_.size());
        displ_.resize(comm_.size() + 1, 0);
        // Send the int for convenience. It will be used to get the place where the
        // data comes from
        using Pair = std::pair<GlobalIndex,int>;
        std::vector<Pair> my_pairs;
        my_pairs.reserve(local_indices_.size());
        for (const auto& pair: local_indices_)
        {
            if (pair.local().attribute() == Attribute::owner)
            {
                my_pairs.emplace_back(pair.global(), -1);
            }
        }
        int mySize = my_pairs.size();
        comm_.allgather(&mySize, 1, sizes_.data());
        std::partial_sum(sizes_.begin(), sizes_.end(), displ_.begin()+1);
        std::vector<Pair> global_pairs(displ_.back());
        comm_.allgatherv(my_pairs.data(), my_pairs.size(), global_pairs.data(), sizes_.data(), displ_.data());
        // Set the the index where we receive
        int count = 0;
        std::for_each(global_pairs.begin(), global_pairs.end(), [&count](Pair& pair){ pair.second = count++;});
        // sort the complete range to get the correct ordering
        std::sort(global_pairs.begin(), global_pairs.end(),
                  [](const Pair& p1, const Pair& p2){ return p1.first < p2.first; } );
        map_received_.resize(global_pairs.size());
        std::transform(global_pairs.begin(), global_pairs.end(), map_received_.begin(),
                       [](const Pair& pair){ return pair.second; });
        perf_ecl_index_.resize(global_pairs.size());
        std::transform(global_pairs.begin(), global_pairs.end(), perf_ecl_index_.begin(),
                       [](const Pair& pair){ return pair.first; });
        num_global_perfs_ = global_pairs.size();
    }
    else
    {
        num_global_perfs_ = num_local_perfs;
    }
}

template<class Scalar>
int GlobalPerfContainerFactory<Scalar>::localToGlobal(std::size_t localIndex) const {
    // todo: maybe Check local global Id with size of communicator!
    // todo: check why this does not work:
    //return local_indices_[localIndex].global();
    if (local_indices_.begin() == local_indices_.end())
        return localIndex;
    auto it = std::find_if(local_indices_.begin(), local_indices_.end(), [localIndex](const auto& index) {return index.local() == localIndex;});
    if (it == local_indices_.end())
        OPM_THROW(std::logic_error, "There is no global index for the localIndex " + std::to_string(localIndex) + ".");
    return it->global();
}

template<class Scalar>
int GlobalPerfContainerFactory<Scalar>::globalToLocal(const int globalIndex) const {
    if (local_indices_.begin() == local_indices_.end())
        return globalIndex;
    auto it = std::find_if(local_indices_.begin(), local_indices_.end(), [globalIndex](const auto& index) {return index.global() == globalIndex;});
    if (it == local_indices_.end())
        return -1;
    return it->local().local();
}

template<class Scalar>
std::vector<Scalar> GlobalPerfContainerFactory<Scalar>::
createGlobal(const std::vector<Scalar>& local_perf_container,
             std::size_t num_components) const
{
    // Could be become templated later.
    using Value = Scalar;

    if (comm_.size() > 1)
    {
        std::vector<Value> global_recv(perf_ecl_index_.size() * num_components);
        if (num_components == 1)
        {
            comm_.allgatherv(local_perf_container.data(), local_perf_container.size(),
                             global_recv.data(), const_cast<int*>(sizes_.data()),
                             const_cast<int*>(displ_.data()));
        }
        else
        {
#if HAVE_MPI
            // Create MPI type for sending num_components entries
            MPI_Datatype data_type;
            MPI_Type_contiguous(num_components, Dune::MPITraits<Value>::getType(), &data_type);
            MPI_Type_commit(&data_type);
            MPI_Allgatherv(local_perf_container.data(),
                           local_perf_container.size()/num_components,
                           data_type, global_recv.data(), sizes_.data(),
                           displ_.data(), data_type, comm_);
            MPI_Type_free(&data_type);
#endif
        }

        // reorder by ascending ecl index.
        std::vector<Value> global_remapped(perf_ecl_index_.size() * num_components);
        auto global = global_remapped.begin();
        for (auto map_entry = map_received_.begin(); map_entry !=  map_received_.end(); ++map_entry)
        {
            auto global_index = *map_entry * num_components;

            for(std::size_t i = 0; i < num_components; ++i)
                *(global++) = global_recv[global_index++];
        }
        assert(global == global_remapped.end());
        return global_remapped;
    }
    else
    {
        return local_perf_container;
    }
}

template<class Scalar>
void GlobalPerfContainerFactory<Scalar>::
copyGlobalToLocal(const std::vector<Scalar>& global,
                  std::vector<Scalar>& local,
                  std::size_t num_components) const
{
    if (global.empty())
    {
        return;
    }

    if (comm_.size() > 1)
    {
        // assign the values (both ranges are sorted by the ecl index)
        auto global_perf = perf_ecl_index_.begin();

        for (const auto& pair: local_indices_)
        {
            global_perf = std::lower_bound(global_perf, perf_ecl_index_.end(), pair.global());
            assert(global_perf != perf_ecl_index_.end());
            assert(*global_perf == pair.global());
            auto local_index = pair.local() * num_components;
            auto global_index = (global_perf - perf_ecl_index_.begin()) * num_components;
            for (std::size_t i = 0; i < num_components; ++i)
                local[local_index++] = global[global_index++];
        }
    }
    else
    {
        std::copy(global.begin(), global.end(), local.begin());
    }
}

template<class Scalar>
int GlobalPerfContainerFactory<Scalar>::numGlobalPerfs() const
{
    return num_global_perfs_;
}

template<class Scalar>
CommunicateAboveBelow<Scalar>::
CommunicateAboveBelow([[maybe_unused]] const Parallel::Communication& comm)
#if HAVE_MPI
    : comm_(comm), interface_(comm_)
#endif
{}

template<class Scalar>
void CommunicateAboveBelow<Scalar>::clear()
{
#if HAVE_MPI
    above_indices_ = {};
    current_indices_ = {};
    interface_.free();
    communicator_.free();
#endif
    num_local_perfs_ = 0;
}

template<class Scalar>
void CommunicateAboveBelow<Scalar>::beginReset()
{
    clear();
#if HAVE_MPI
    if (comm_.size() > 1)
    {
        current_indices_.beginResize();
        above_indices_.beginResize();
    }
#endif
}

template<class Scalar>
int CommunicateAboveBelow<Scalar>::endReset()
{
#if HAVE_MPI
    if (comm_.size() > 1)
    {
        above_indices_.endResize();
        current_indices_.endResize();
        remote_indices_.setIndexSets(current_indices_, above_indices_, comm_);
        // It is mandatory to not set includeSelf to true, as this will skip some entries.
        remote_indices_.template rebuild<true>();
        using FromSet = Dune::EnumItem<Attribute,owner>;
        using ToSet = Dune::AllSet<Attribute>;
        interface_.build(remote_indices_, FromSet(), ToSet());
        communicator_.build<Scalar*>(interface_);
    }
#endif
    return num_local_perfs_;
}

template<class Scalar>
template<class RAIterator>
void CommunicateAboveBelow<Scalar>::partialSumPerfValues(RAIterator begin, RAIterator end) const
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
        std::vector<int> displ(comm_.size() + 1, 0);
        using GlobalIndex = typename IndexSet::IndexPair::GlobalIndex;
        using Pair = std::pair<GlobalIndex,Value>;
        std::vector<Pair> my_pairs;
        my_pairs.reserve(current_indices_.size());
        for (const auto& pair: current_indices_)
        {
            if (pair.local().attribute() == owner)
            {
                my_pairs.emplace_back(pair.global(), begin[pair.local()]);
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
            begin[pair.local()] = sums[global_pair - global_pairs.begin()];
        }
#else
        OPM_THROW(std::logic_error, "In a sequential run the size of the communicator should be 1!");
#endif
    }
}

template<class Scalar>
struct CopyGatherScatter
{
    static const Scalar& gather(const Scalar* a, std::size_t i)
    {
        return a[i];
    }

    static void scatter(Scalar* a, const Scalar& v, std::size_t i)
    {
        a[i] = v;
    }
};

template<class Scalar>
std::vector<Scalar> CommunicateAboveBelow<Scalar>::
communicateAbove(Scalar first_above,
                 const Scalar* current,
                 std::size_t size)
{
    std::vector<Scalar> above(size, first_above);

#if HAVE_MPI
    if (comm_.size() > 1)
    {
        auto aboveData = above.data();
        // Ugly const_cast needed as my compiler says, that
        // passing const double*& and double* as parameter is
        // incompatible with function decl template<Data> forward(const Data&, Data&))
        // That would need the first argument to be double* const&
        communicator_.forward<CopyGatherScatter<Scalar>>(const_cast<Scalar*>(current), aboveData);
    }
    else
#endif
    {
        if (above.size() > 1)
        {
            // No comunication needed, just copy.
            std::copy(current, current + (above.size() - 1), above.begin()+1);
        }
    }
    return above;
}

template<class Scalar>
std::vector<Scalar> CommunicateAboveBelow<Scalar>::
communicateBelow(Scalar last_below,
                 const Scalar* current,
                 std::size_t size)
{
    std::vector<Scalar> below(size, last_below);

#if HAVE_MPI
    if (comm_.size() > 1)
    {
        auto belowData = below.data();
        // Ugly const_cast needed as my compiler says, that
        // passing const double*& and double* as parameter is
        // incompatible with function decl template<Data> backward(Data&, const Data&)
        // That would need the first argument to be double* const&
        communicator_.backward<CopyGatherScatter<Scalar>>(belowData, const_cast<Scalar*>(current));
    }
    else
#endif
    {
        if (below.size() > 1)
        {
            // No comunication needed, just copy.
            std::copy(current+1, current + below.size(), below.begin());
        }
    }
    return below;
}

template<class Scalar>
void CommunicateAboveBelow<Scalar>::
pushBackEclIndex([[maybe_unused]] int above,
                 [[maybe_unused]] int current,
                 [[maybe_unused]] bool isOwner)
{
#if HAVE_MPI
    if (comm_.size() > 1)
    {
        Attribute attr = owner;
        if (!isOwner)
        {
            attr = overlap;
        }
        above_indices_.add(above, {num_local_perfs_, attr, true});
        current_indices_.add(current, {num_local_perfs_, attr, true});
    }
#endif
    ++num_local_perfs_;
}

template<class Scalar>
void ParallelWellInfo<Scalar>::DestroyComm::operator()(Parallel::Communication* comm)
{
#if HAVE_MPI
    // Only delete custom communicators.
    bool del = comm
        && (*comm != Dune::MPIHelper::getLocalCommunicator())
        && (*comm != MPI_COMM_WORLD && *comm != MPI_COMM_NULL);

    if ( del )
    {
        // Not 100% nice but safe as comm is deleted anyway
        // We can only access a copy and no reference.
        MPI_Comm mpi_comm = *comm;
        MPI_Comm_free(&mpi_comm);
    }
#endif
    delete comm;
}

template<class Scalar>
const typename CommunicateAboveBelow<Scalar>::IndexSet&
CommunicateAboveBelow<Scalar>::getIndexSet() const
{
    return current_indices_;
}

template<class Scalar>
int CommunicateAboveBelow<Scalar>::numLocalPerfs() const
{
    return num_local_perfs_;
}

template<class Scalar>
ParallelWellInfo<Scalar>::ParallelWellInfo(const std::string& name,
                                           bool hasLocalCells)
    : name_(name), hasLocalCells_ (hasLocalCells),
      isOwner_(true), rankWithFirstPerf_(-1),
      comm_(new Parallel::Communication(Dune::MPIHelper::getLocalCommunicator())),
      commAboveBelow_(new CommunicateAboveBelow<Scalar>(*comm_))
{}

template<class Scalar>
ParallelWellInfo<Scalar>::ParallelWellInfo(const std::pair<std::string, bool>& well_info,
                                           [[maybe_unused]] Parallel::Communication allComm)
    : name_(well_info.first), hasLocalCells_(well_info.second),
      rankWithFirstPerf_(-1)
{
#if HAVE_MPI
    MPI_Comm newComm;
    int color = hasLocalCells_ ? 1 : MPI_UNDEFINED;
    MPI_Comm_split(allComm, color, allComm.rank(), &newComm);
    comm_.reset(new Parallel::Communication(newComm));
#else
    comm_.reset(new Parallel::Communication(Dune::MPIHelper::getLocalCommunicator()));
#endif
    commAboveBelow_.reset(new CommunicateAboveBelow<Scalar>(*comm_));
    isOwner_ = (comm_->rank() == 0);
}

template<class Scalar>
void ParallelWellInfo<Scalar>::communicateFirstPerforation(bool hasFirst)
{
    int first = hasFirst;
    std::vector<int> firstVec(comm_->size());
    comm_->allgather(&first, 1, firstVec.data());
    auto found = std::find_if(firstVec.begin(), firstVec.end(),
                              [](int i) -> bool{ return i;});
    if (found != firstVec.end())
        rankWithFirstPerf_ = found - firstVec.begin();
}

template<class Scalar>
void ParallelWellInfo<Scalar>::pushBackEclIndex(int above, int current)
{
    commAboveBelow_->pushBackEclIndex(above, current);
}

template<class Scalar>
void ParallelWellInfo<Scalar>::beginReset()
{
    commAboveBelow_->beginReset();
}

template<class Scalar>
void ParallelWellInfo<Scalar>::endReset()
{
    int local_num_perfs = commAboveBelow_->endReset();
    globalPerfCont_
        .reset(new GlobalPerfContainerFactory<Scalar>(commAboveBelow_->getIndexSet(),
                                                      *comm_,
                                                      local_num_perfs));
}

template<class Scalar>
template<typename It>
typename It::value_type
ParallelWellInfo<Scalar>::sumPerfValues(It begin, It end) const
{
    using V = typename It::value_type;
    /// \todo cater for overlap later. Currently only owner
    auto local = std::accumulate(begin, end, V());
    return communication().sum(local);
}

template<class Scalar>
void ParallelWellInfo<Scalar>::clear()
{
    commAboveBelow_->clear();
    globalPerfCont_.reset();
}

template<class Scalar>
template<class T>
T ParallelWellInfo<Scalar>::broadcastFirstPerforationValue(const T& t) const
{
    T res = t;
    if (rankWithFirstPerf_ >= 0) {
#ifndef NDEBUG
        assert(rankWithFirstPerf_ < comm_->size());
        // At least on some OpenMPI version this might broadcast might interfere
        // with other communication if there are bugs
        comm_->barrier();
#endif
        comm_->broadcast(&res, 1, rankWithFirstPerf_);
#ifndef NDEBUG
        comm_->barrier();
#endif
    }
    return res;
}

template<class Scalar>
std::vector<Scalar> ParallelWellInfo<Scalar>::
communicateAboveValues(Scalar zero_value,
                       const Scalar* current_values,
                       std::size_t size) const
{
    return commAboveBelow_->communicateAbove(zero_value, current_values,
                                             size);
}

template<class Scalar>
std::vector<Scalar> ParallelWellInfo<Scalar>::
communicateAboveValues(Scalar zero_value,
                       const std::vector<Scalar>& current_values) const
{
    return commAboveBelow_->communicateAbove(zero_value, current_values.data(),
                                             current_values.size());
}

template<class Scalar>
std::vector<Scalar> ParallelWellInfo<Scalar>::
communicateBelowValues(Scalar last_value,
                       const Scalar* current_values,
                       std::size_t size) const
{
    return commAboveBelow_->communicateBelow(last_value, current_values,
                                             size);
}

template<class Scalar>
std::vector<Scalar> ParallelWellInfo<Scalar>::
communicateBelowValues(Scalar last_value,
                       const std::vector<Scalar>& current_values) const
{
    return commAboveBelow_->communicateBelow(last_value, current_values.data(),
                                             current_values.size());
}

template<class Scalar>
const GlobalPerfContainerFactory<Scalar>&
ParallelWellInfo<Scalar>::getGlobalPerfContainerFactory() const
{
    if(globalPerfCont_)
        return *globalPerfCont_;
    else
        OPM_THROW(std::logic_error,
                  "No ecl indices have been added via beginReset, pushBackEclIndex, endReset");
}

template<class Scalar>
bool operator<(const ParallelWellInfo<Scalar>& well1, const ParallelWellInfo<Scalar>& well2)
{
    return well1.name() < well2.name() || (! (well2.name() < well1.name()) &&
                                              well1.hasLocalCells() < well2.hasLocalCells());
}

template<class Scalar>
bool operator==(const ParallelWellInfo<Scalar>& well1, const ParallelWellInfo<Scalar>& well2)
{
    bool ret = well1.name() == well2.name() && well1.hasLocalCells() == well2.hasLocalCells()
        && well1.isOwner() == well2.isOwner();
#if HAVE_MPI
    using MPIComm = typename Dune::MPIHelper::MPICommunicator;
    ret = ret &&
        static_cast<MPIComm>(well1.communication()) == static_cast<MPIComm>(well2.communication());
#endif
    return ret;
}

template<class Scalar>
bool operator!=(const ParallelWellInfo<Scalar>& well1, const ParallelWellInfo<Scalar>& well2)
{
    return ! (well1 == well2);
}

template<class Scalar>
bool operator<(const std::pair<std::string, bool>& pair, const ParallelWellInfo<Scalar>& well)
{
    return pair.first < well.name() || ( !( well.name() < pair.first ) && pair.second < well.hasLocalCells() );
}

template<class Scalar>
bool operator<( const ParallelWellInfo<Scalar>& well, const std::pair<std::string, bool>& pair)
{
    return well.name() < pair.first || ( !( pair.first < well.name() ) && well.hasLocalCells() < pair.second );
}

template<class Scalar>
bool operator==(const std::pair<std::string, bool>& pair, const ParallelWellInfo<Scalar>& well)
{
    return pair.first == well.name() && pair.second == well.hasLocalCells();
}

template<class Scalar>
bool operator==(const ParallelWellInfo<Scalar>& well, const std::pair<std::string, bool>& pair)
{
    return pair == well;
}

template<class Scalar>
bool operator!=(const std::pair<std::string, bool>& pair, const ParallelWellInfo<Scalar>& well)
{
    return pair.first != well.name() || pair.second != well.hasLocalCells();
}

template<class Scalar>
bool operator!=(const ParallelWellInfo<Scalar>& well, const std::pair<std::string, bool>& pair)
{
    return pair != well;
}

template<class Scalar>
CheckDistributedWellConnections<Scalar>::
CheckDistributedWellConnections(const Well& well,
                                const ParallelWellInfo<Scalar>& info)
    : well_(well), pwinfo_(info)
{
    foundConnections_.resize(well.getConnections().size(), 0);
}

template<class Scalar>
void CheckDistributedWellConnections<Scalar>::
connectionFound(std::size_t i)
{
    foundConnections_[i] = 1;
}

template<class Scalar>
bool CheckDistributedWellConnections<Scalar>::
checkAllConnectionsFound()
{
    // Ecl does not hold any information of remote connections.
    assert(pwinfo_.communication().max(foundConnections_.size()) == foundConnections_.size());
    pwinfo_.communication().sum(foundConnections_.data(),
                                foundConnections_.size());

    std::string msg = std::string("Cells with these i,j,k indices were not found ")
        + "in grid (well = " + pwinfo_.name() + "):";
    bool missingCells = false;
    auto start = foundConnections_.begin();
    for(auto conFound = start; conFound != foundConnections_.end(); ++conFound)
    {
        if (*conFound == 0)
        {
            const auto& completion = well_.getConnections()[conFound - start];
            msg = msg + " " + std::to_string(completion.getI()) + "," +
                std::to_string(completion.getJ()) + ","
                + std::to_string(completion.getK()) + " ";
            missingCells = true;
        }
    }
    if (missingCells && pwinfo_.isOwner())
    {
        OPM_THROW(std::runtime_error, msg);
    }
    return !missingCells;
}

template<class Scalar> using dIter = typename std::vector<Scalar>::iterator;
template<class Scalar> using cdIter = typename std::vector<Scalar>::const_iterator;

#define INSTANTIATE_TYPE(T)                                                         \
    template class CheckDistributedWellConnections<T>;                              \
    template class CommunicateAboveBelow<T>;                                        \
    template class GlobalPerfContainerFactory<T>;                                   \
    template class ParallelWellInfo<T>;                                             \
    template typename cdIter<T>::value_type                                         \
        ParallelWellInfo<T>::sumPerfValues<cdIter<T>>(cdIter<T>,cdIter<T>) const;   \
    template typename dIter<T>::value_type                                          \
        ParallelWellInfo<T>::sumPerfValues<dIter<T>>(dIter<T>,dIter<T>) const;      \
     template int ParallelWellInfo<T>::                                             \
        broadcastFirstPerforationValue<int>(const int&) const;                      \
     template T ParallelWellInfo<T>::                                               \
        broadcastFirstPerforationValue<T>(const T&) const;                          \
    template void CommunicateAboveBelow<T>::                                        \
        partialSumPerfValues<dIter<T>>(dIter<T>, dIter<T>) const;                   \
    template bool operator<(const ParallelWellInfo<T>&,                             \
                            const ParallelWellInfo<T>&);                            \
    template bool operator<(const ParallelWellInfo<T>&,                             \
                            const std::pair<std::string, bool>&);                   \
    template bool operator<(const std::pair<std::string, bool>&,                    \
                            const ParallelWellInfo<T>&);                            \
    template bool operator==(const ParallelWellInfo<T>&,                            \
                             const ParallelWellInfo<T>&);                           \
    template bool operator==(const ParallelWellInfo<T>& well,                       \
                             const std::pair<std::string, bool>&);                  \
    template bool operator==(const std::pair<std::string, bool>&,                   \
                             const ParallelWellInfo<T>&);                           \
    template bool operator!=(const ParallelWellInfo<T>&,                            \
                             const ParallelWellInfo<T>&);                           \
    template bool operator!=(const std::pair<std::string, bool>&,                   \
                             const ParallelWellInfo<T>&);                           \
    template bool operator!=(const ParallelWellInfo<T>&,                            \
                             const std::pair<std::string, bool>&);

INSTANTIATE_TYPE(double)

#if FLOW_INSTANTIATE_FLOAT
INSTANTIATE_TYPE(float)
#endif

} // end namespace Opm
