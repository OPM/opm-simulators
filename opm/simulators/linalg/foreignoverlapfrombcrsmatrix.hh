// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 2 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.

  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the list of
  copyright holders.
*/
/*!
 * \file
 * \copydoc Opm::Linear::ForeignOverlapFromBCRSMatrix
 */
#ifndef EWOMS_FOREIGN_OVERLAP_FROM_BCRS_MATRIX_HH
#define EWOMS_FOREIGN_OVERLAP_FROM_BCRS_MATRIX_HH

#include "overlaptypes.hh"
#include "blacklist.hh"

#include <opm/models/parallel/mpibuffer.hh>

#include <opm/material/common/Unused.hpp>

#include <dune/grid/common/datahandleif.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/scalarproducts.hh>
#include <dune/istl/operators.hh>

#include <algorithm>
#include <iostream>
#include <map>
#include <vector>

#if HAVE_MPI
#include <mpi.h>
#endif // HAVE_MPI

namespace Opm {
namespace Linear {

/*!
 * \brief This class creates and manages the foreign overlap given an
 *        initial list of border indices and a BCRS matrix.
 *
 * The foreign overlap are all (row) indices which overlap with the
 * some of the current process's local indices.
 */
class ForeignOverlapFromBCRSMatrix
{
public:
    // overlaps should never be copied!
    ForeignOverlapFromBCRSMatrix(const ForeignOverlapFromBCRSMatrix&) = delete;

    /*!
     * \brief Constructs the foreign overlap given a BCRS matrix and
     *        an initial list of border indices.
     */
    template <class BCRSMatrix>
    ForeignOverlapFromBCRSMatrix(const BCRSMatrix& A,
                                 const BorderList& borderList,
                                 const BlackList& blackList,
                                 unsigned overlapSize)
        : borderList_(borderList), blackList_(blackList)
    {
        overlapSize_ = overlapSize;

        myRank_ = 0;
#if HAVE_MPI
        {
            int tmp;
            MPI_Comm_rank(MPI_COMM_WORLD, &tmp);
            myRank_ = static_cast<ProcessRank>(tmp);
        }
#endif
        numNative_ = A.N();

        // Computes the local <-> native index maps
        createLocalIndices_();

        // calculate the set of local indices on the border (beware:
        // _not_ the native ones)
        auto it = borderList.begin();
        const auto& endIt = borderList.end();
        for (; it != endIt; ++it) {
            Index localIdx = nativeToLocal(it->localIdx);
            if (localIdx < 0)
                continue;

            localBorderIndices_.insert(localIdx);
        }

        // compute the set of processes which are neighbors of the
        // local process ...
        neighborPeerSet_.update(borderList);
        // ... and the initial set of processes which we will have to
        // communicate with. We must always communicate with our
        // neighbors, but depending on the size of the overlap region,
        // we might have to communicate with additional processes as
        // well (these will be added later).
        peerSet_ = neighborPeerSet_;

        // Create an initial seed list of indices which are in the
        // overlap.
        SeedList initialSeedList;
        initialSeedList.update(borderList);

        // calculate the minimum distance from the border of the
        // initial seed list
        unsigned minBorderDist = overlapSize;
        auto borderIt = borderList.begin();
        const auto& borderEndIt = borderList.end();
        for (; borderIt != borderEndIt; ++borderIt) {
            minBorderDist = std::min(minBorderDist, borderIt->borderDistance);
        }

        // calculate the foreign overlap for the local partition,
        // i.e. find the distance of each row from the seed set.
        foreignOverlapByLocalIndex_.resize(numLocal());
        extendForeignOverlap_(A, initialSeedList, minBorderDist, overlapSize);

        // computes the process with the lowest rank for all local
        // indices.
        computeMasterRanks_();

        // group foreign overlap by peer process rank
        groupForeignOverlapByRank_();
    }

    /*!
     * \brief Returns the size of the overlap region.
     */
    unsigned overlapSize() const
    { return overlapSize_; }

    /*!
     * \brief Returns true iff a local index is a border index.
     */
    bool isBorder(Index localIdx) const
    { return localBorderIndices_.count(localIdx) > 0; }

    /*!
     * \brief Returns true iff a local index is a border index shared with a
     * given peer process.
     */
    bool isBorderWith(Index localIdx, ProcessRank peerRank) const
    {
        const auto& indexOverlap = foreignOverlapByLocalIndex_[static_cast<unsigned>(localIdx)];
        const auto& borderDistIt = indexOverlap.find(peerRank);
        if (borderDistIt == indexOverlap.end())
            return false;

        // border distance of the index needs to be 0
        return borderDistIt->second == 0;
    }

    /*!
     * \brief Return the rank of the master process of an
     *        index.
     */
    ProcessRank masterRank(Index localIdx) const
    { return masterRank_[static_cast<unsigned>(localIdx)]; }

    /*!
     * \brief Return true if the current rank is the "master" of an
     *        index.
     *
     * If the index is at the interior of some process, we define this
     * process as its master, if the index is on the boundary, then
     * the master is defined as the process with the lowest rank.
     */
    bool iAmMasterOf(Index localIdx) const
    { return masterRank_[static_cast<unsigned>(localIdx)] == myRank_; }

    /*!
     * \brief Returns the list of indices which intersect the process
     *        border.
     */
    const BorderList& borderList() const
    { return borderList_; }

    /*!
     * \brief Return the list of (local indices, border distance,
     *        number of processes) triples which are in the overlap of
     *        a given peer rank.
     */
    const OverlapWithPeer& foreignOverlapWithPeer(ProcessRank peerRank) const
    {
        assert(foreignOverlapByRank_.find(peerRank) != foreignOverlapByRank_.end());
        return foreignOverlapByRank_.find(peerRank)->second;
    }

    /*!
     * \brief Return the map of (peer rank, border distance) for a given local
     * index.
     */
    const std::map<ProcessRank, BorderDistance> &
    foreignOverlapByLocalIndex(Index localIdx) const
    {
        assert(isLocal(localIdx));
        return foreignOverlapByLocalIndex_[static_cast<unsigned>(localIdx)];
    }

    /*!
     * \brief Returns true iff a local index is seen by a peer rank.
     */
    bool peerHasIndex(ProcessRank peerRank, Index localIdx) const
    {
        const auto& idxOverlap = foreignOverlapByLocalIndex_[localIdx];
        return idxOverlap.find(peerRank) != idxOverlap.end();
    }

    /*!
     * \brief Returns the number of front indices of a peer process in
     *        the local partition.
     */
    size_t numFront(ProcessRank peerRank) const
    {
        const auto& peerOverlap = foreignOverlapByRank_.find(peerRank)->second;

        size_t n = 0;
        auto it = peerOverlap.begin();
        const auto& endIt = peerOverlap.end();
        for (; it != endIt; ++it) {
            if (it->borderDistance == overlapSize_)
                ++n;
        }
        return n;
    }

    /*!
     * \brief Returns whether a given local index is on the front of a
     *        given peer rank.
     */
    bool isFrontFor(ProcessRank peerRank, Index localIdx) const
    {
        const auto& idxOverlap = foreignOverlapByLocalIndex_[localIdx];

        auto it = idxOverlap.find(peerRank);
        if (it == idxOverlap.end())
            return false; // index is not in overlap

        return it->second == overlapSize_;
    }

    /*!
     * \brief Return the set of process ranks which share an overlap
     *        with the current process.
     */
    const PeerSet& peerSet() const
    { return peerSet_; }

    /*!
     * \brief Return the set of process ranks which share a border index
     *        with the current process.
     */
    const PeerSet& neighborPeerSet() const
    { return neighborPeerSet_; }

    /*!
     * \brief Returns the number of native indices
     */
    size_t numNative() const
    { return numNative_; }

    /*!
     * \brief Returns the number of local indices
     */
    size_t numLocal() const
    { return numLocal_; }

    /*!
     * \brief Returns true iff a domestic index is local
     */
    bool isLocal(Index domesticIdx) const
    { return static_cast<unsigned>(domesticIdx) < numLocal(); }

    /*!
     * \brief Convert a native index to a local one.
     *
     * If a given native index is not in the set of local indices,
     * this method returns -1.
     */
    Index nativeToLocal(Index nativeIdx) const
    { return nativeToLocalIndices_[static_cast<unsigned>(nativeIdx)]; }

    /*!
     * \brief Convert a local index to a native one.
     */
    Index localToNative(Index localIdx) const
    {
        assert(localIdx < static_cast<Index>(localToNativeIndices_.size()));
        return localToNativeIndices_[static_cast<unsigned>(localIdx)];
    }

    /*!
     * \brief Returns the object which represents the black-listed native indices.
     */
    const BlackList& blackList() const
    { return blackList_; }

    /*!
     * \brief Return the number of peer ranks for which a given local
     *        index is visible.
     */
    size_t numPeers(Index localIdx) const
    { return foreignOverlapByLocalIndex_[static_cast<unsigned>(localIdx)].size(); }

    /*!
     * \brief Returns true if a given local index is in the foreign overlap of
     * any rank.
     */
    bool isInOverlap(Index localIdx) const
    { return foreignOverlapByLocalIndex_[static_cast<unsigned>(localIdx)].size() > 0; }

    /*!
     * \brief Print the foreign overlap for debugging purposes.
     */
    void print() const
    {
        auto it = foreignOverlapByRank_.begin();
        const auto& endIt = foreignOverlapByRank_.end();
        for (; it != endIt; ++it) {
            std::cout << "Overlap rows(distance) for rank " << it->first << ": ";

            auto rowIt = it->second.begin();
            const auto& rowEndIt = it->second.end();
            for (; rowIt != rowEndIt; ++rowIt)
                std::cout << rowIt->index << "(" << rowIt->borderDistance << ") ";
            std::cout << "\n" << std::flush;
        }
    }

protected:
    // extend the foreign overlaps by 'overlapSize' levels. this uses
    // a greedy algorithm which extends the region by one level and
    // then calls itself recursively...
    template <class BCRSMatrix>
    void extendForeignOverlap_(const BCRSMatrix& A,
                               SeedList& seedList,
                               BorderDistance borderDistance,
                               BorderDistance overlapSize)
    {
        // communicate the non-neigbor overlap indices
        addNonNeighborOverlapIndices_(A, seedList, borderDistance);

        // add all processes in the seed rows of the current overlap level
        auto seedIt = seedList.begin();
        const auto& seedEndIt = seedList.end();
        for (; seedIt != seedEndIt; ++seedIt) {
            Index localIdx = nativeToLocal(seedIt->index);
            ProcessRank peerRank = seedIt->peerRank;
            unsigned distance = borderDistance;
            if (localIdx < 0)
                continue;
            if (foreignOverlapByLocalIndex_[static_cast<unsigned>(localIdx)].count(peerRank) == 0)
                foreignOverlapByLocalIndex_[static_cast<unsigned>(localIdx)][peerRank] = distance;
        }

        // if we have reached the maximum overlap distance, i.e. we're
        // finished and break the recursion
        if (borderDistance >= overlapSize)
            return;

        // find the seed list for the next overlap level using the
        // seed set for the current level
        SeedList nextSeedList;
        seedIt = seedList.begin();
        for (; seedIt != seedEndIt; ++seedIt) {
            Index nativeRowIdx = seedIt->index;
            if (nativeToLocal(nativeRowIdx) < 0)
                continue; // ignore blacklisted indices
            ProcessRank peerRank = seedIt->peerRank;

            // find all column indices in the row. The indices of the
            // columns are the additional indices of the overlap which
            // we would like to add
            typedef typename BCRSMatrix::ConstColIterator ColIterator;
            ColIterator colIt = A[static_cast<unsigned>(nativeRowIdx)].begin();
            ColIterator colEndIt = A[static_cast<unsigned>(nativeRowIdx)].end();
            for (; colIt != colEndIt; ++colIt) {
                Index nativeColIdx = static_cast<Index>(colIt.index());
                Index localColIdx = nativeToLocal(nativeColIdx);

                // ignore if the native index is not a local one
                if (localColIdx < 0)
                    continue;
                // if the process is already is in the overlap of the
                // column index, ignore this column index!
                else if (foreignOverlapByLocalIndex_[static_cast<unsigned>(localColIdx)].count(peerRank) > 0)
                    continue;

                // check whether the new index is already in the overlap
                bool hasIndex = false;
                typename SeedList::iterator sIt = nextSeedList.begin();
                typename SeedList::iterator sEndIt = nextSeedList.end();
                for (; sIt != sEndIt; ++sIt) {
                    if (sIt->index == nativeColIdx && sIt->peerRank == peerRank) {
                        hasIndex = true;
                        break;
                    }
                }
                if (hasIndex)
                    continue; // we already have this index

                // add the current processes to the seed list for the
                // next overlap level
                IndexRankDist newTuple;
                newTuple.index = nativeColIdx;
                newTuple.peerRank = peerRank;
                newTuple.borderDistance = seedIt->borderDistance + 1;
                nextSeedList.push_back(newTuple);
            }
        }

        // clear the old seed list to save some memory
        seedList.clear();

        // Perform the same excercise for the next overlap distance
        extendForeignOverlap_(A, nextSeedList, borderDistance + 1, overlapSize);
    }

    // Computes the local <-> native index maps
    void createLocalIndices_()
    {
        // create the native <-> local maps
        Index localIdx = 0;
        for (unsigned nativeIdx = 0; nativeIdx < numNative_;) {
            if (!blackList_.hasIndex(static_cast<Index>(nativeIdx))) {
                localToNativeIndices_.push_back(static_cast<Index>(nativeIdx));
                nativeToLocalIndices_.push_back(static_cast<Index>(localIdx));
                ++nativeIdx;
                ++localIdx;
            }
            else {
                nativeToLocalIndices_.push_back(-1);
                ++nativeIdx;
            }
        }

        numLocal_ = localToNativeIndices_.size();
    }

    Index localToPeerIdx_(Index localIdx, ProcessRank peerRank) const
    {
        auto it = borderList_.begin();
        const auto& endIt = borderList_.end();
        for (; it != endIt; ++it) {
            if (it->localIdx == localIdx && it->peerRank == peerRank)
                return it->peerIdx;
        }

        return -1;
    }

    template <class BCRSMatrix>
    void addNonNeighborOverlapIndices_(const BCRSMatrix& A OPM_UNUSED,
                                       SeedList& seedList OPM_UNUSED_NOMPI,
                                       BorderDistance borderDist OPM_UNUSED_NOMPI)
    {
        // TODO: this probably does not work! (the matrix A is unused, but it is needed
        // from a logical POV.)
#if HAVE_MPI
        // first, create the buffers which will contain the number of
        // border indices relevant for a neighbor peer
        std::map<ProcessRank, std::vector<BorderIndex> > borderIndices;

        // get all indices in the border which have borderDist as
        // their distance to the closest border of their local process
        auto it = seedList.begin();
        const auto& endIt = seedList.end();
        for (; it != endIt; ++it) {
            Index localIdx = nativeToLocal(it->index);
            if (!isBorder(localIdx))
                continue;
            BorderIndex borderHandle;
            borderHandle.localIdx = localIdx;
            borderHandle.peerRank = it->peerRank;
            borderHandle.borderDistance = it->borderDistance;

            // add the border index to all the neighboring peers
            auto neighborIt = foreignOverlapByLocalIndex_[static_cast<unsigned>(localIdx)].begin();
            const auto& neighborEndIt = foreignOverlapByLocalIndex_[static_cast<unsigned>(localIdx)].end();
            for (; neighborIt != neighborEndIt; ++neighborIt) {
                if (neighborIt->second != 0)
                    // not a border index for the neighbor
                    continue;
                else if (neighborIt->first == borderHandle.peerRank)
                    // don't communicate the indices which are owned
                    // by the peer to itself
                    continue;

                Index peerIdx = localToPeerIdx_(localIdx, neighborIt->first);
                if (peerIdx < 0)
                    // the index is on the border, but is not on the border
                    // with the considered neighboring process. Ignore it!
                    continue;
                borderHandle.peerIdx = peerIdx;
                borderIndices[neighborIt->first].push_back(borderHandle);
            }
        }

        // now borderIndices contains the lists of indices which we
        // would like to send to each neighbor. Let's create the MPI
        // buffers.
        std::map<ProcessRank, Opm::MpiBuffer<unsigned> > numIndicesSendBufs;
        std::map<ProcessRank, Opm::MpiBuffer<BorderIndex> > indicesSendBufs;
        auto peerIt = neighborPeerSet().begin();
        const auto& peerEndIt = neighborPeerSet().end();
        for (; peerIt != peerEndIt; ++peerIt) {
            ProcessRank peerRank = *peerIt;
            size_t numIndices = borderIndices[peerRank].size();
            numIndicesSendBufs[peerRank].resize(1);
            numIndicesSendBufs[peerRank][0] = static_cast<unsigned>(numIndices);

            const auto& peerBorderIndices = borderIndices[peerRank];
            indicesSendBufs[peerRank].resize(numIndices);

            auto tmpIt = peerBorderIndices.begin();
            const auto& tmpEndIt = peerBorderIndices.end();
            size_t i = 0;
            for (; tmpIt != tmpEndIt; ++tmpIt, ++i) {
                indicesSendBufs[peerRank][i] = *tmpIt;
            }
        }

        // now, send all these nice buffers to our neighbors
        peerIt = neighborPeerSet().begin();
        for (; peerIt != peerEndIt; ++peerIt) {
            ProcessRank neighborPeer = *peerIt;
            numIndicesSendBufs[neighborPeer].send(neighborPeer);
            indicesSendBufs[neighborPeer].send(neighborPeer);
        }

        // receive all data from the neighbors
        std::map<ProcessRank, MpiBuffer<unsigned> > numIndicesRcvBufs;
        std::map<ProcessRank, MpiBuffer<BorderIndex> > indicesRcvBufs;
        peerIt = neighborPeerSet().begin();
        for (; peerIt != peerEndIt; ++peerIt) {
            ProcessRank neighborPeer = *peerIt;
            auto& numIndicesRcvBuf = numIndicesRcvBufs[neighborPeer];
            auto& indicesRcvBuf = indicesRcvBufs[neighborPeer];

            numIndicesRcvBuf.resize(1);
            numIndicesRcvBuf.receive(neighborPeer);
            unsigned numIndices = numIndicesRcvBufs[neighborPeer][0];
            indicesRcvBuf.resize(numIndices);
            indicesRcvBuf.receive(neighborPeer);

            // filter out all indices which are already in the peer
            // processes' overlap and add them to the seed list. also
            // extend the set of peer processes.
            for (unsigned i = 0; i < numIndices; ++i) {
                // swap the local and the peer indices, because they were
                // created with the point view of the sender
                std::swap(indicesRcvBuf[i].localIdx, indicesRcvBuf[i].peerIdx);

                ProcessRank peerRank = indicesRcvBuf[i].peerRank;
                // Index peerIdx = indicesRcvBuf[i].peerIdx;
                Index localIdx = indicesRcvBuf[i].localIdx;

                // check if the index is already in the overlap for
                // the peer
                const auto& distIt = foreignOverlapByLocalIndex_[static_cast<unsigned>(localIdx)].find(peerRank);
                if (distIt != foreignOverlapByLocalIndex_[static_cast<unsigned>(localIdx)].end())
                    continue;

                // make sure the index is not already in the seed list
                bool inSeedList = false;
                auto seedIt = seedList.begin();
                const auto& seedEndIt = seedList.end();
                for (; seedIt != seedEndIt; ++seedIt) {
                    if (seedIt->index == localIdx && seedIt->peerRank == peerRank) {
                        inSeedList = true;
                        break;
                    }
                }
                if (inSeedList)
                    continue;

                IndexRankDist seedEntry;
                seedEntry.index = localIdx;
                seedEntry.peerRank = peerRank;
                seedEntry.borderDistance = borderDist;
                seedList.push_back(seedEntry);

                // update the peer set
                peerSet_.insert(peerRank);
            }
        }

        // make sure all data was send
        peerIt = neighborPeerSet().begin();
        for (; peerIt != peerEndIt; ++peerIt) {
            ProcessRank neighborPeer = *peerIt;
            numIndicesSendBufs[neighborPeer].wait();
            indicesSendBufs[neighborPeer].wait();
        }
#endif // HAVE_MPI
    }

    // given a list of border indices and provided that
    // borderListToSeedList_() was already called, calculate the
    // master process of each local index.
    void computeMasterRanks_()
    {
        // determine the minimum rank for all indices
        masterRank_.resize(numLocal_);
        for (unsigned localIdx = 0; localIdx < numLocal_; ++localIdx) {
            unsigned masterRank = myRank_;
            if (isBorder(static_cast<Index>(localIdx))) {
                // if the local index is a border index, loop over all ranks
                // for which this index is also a border index. the lowest
                // rank wins!
                auto it = foreignOverlapByLocalIndex_[static_cast<unsigned>(localIdx)].begin();
                const auto& endIt = foreignOverlapByLocalIndex_[static_cast<unsigned>(localIdx)].end();
                for (; it != endIt; ++it) {
                    if (it->second == 0) {
                        // if the border distance is zero, the rank with the
                        // minimum
                        masterRank = std::min<ProcessRank>(masterRank, it->first);
                    }
                }
            }
            masterRank_[static_cast<unsigned>(localIdx)] = masterRank;
        }
    }

    // assuming that the foreign overlap has been created for each
    // local index, this method groups the foreign overlap by peer
    // process rank
    void groupForeignOverlapByRank_()
    {
        // loop over all indices which are in the overlap of some
        // process
        size_t numLocal = foreignOverlapByLocalIndex_.size();
        for (unsigned localIdx = 0; localIdx < numLocal; ++localIdx) {
            // loop over the list of processes for the current index
            auto it = foreignOverlapByLocalIndex_[localIdx].begin();
            const auto& endIt = foreignOverlapByLocalIndex_[localIdx].end();
            size_t nRanks = foreignOverlapByLocalIndex_[localIdx].size();
            for (; it != endIt; ++it) {
                IndexDistanceNpeers tmp;
                tmp.index = static_cast<Index>(localIdx);
                tmp.borderDistance = it->second;
                tmp.numPeers = static_cast<unsigned>(nRanks);
                foreignOverlapByRank_[it->first].push_back(tmp);
            }
        }
    }

    // set of processes with which we have to communicate
    PeerSet peerSet_;

    // set of processes which are direct neighbors of us
    PeerSet neighborPeerSet_;

    // the list of indices on the border
    const BorderList& borderList_;

    // the set of indices which should not be considered
    const BlackList& blackList_;

    // local indices are the native indices sans the black listed ones
    std::vector<Index> nativeToLocalIndices_;
    std::vector<Index> localToNativeIndices_;

    // an array which contains the rank of the master process for each
    // index
    std::vector<ProcessRank> masterRank_;

    // set of all local indices which are on the border of some remote
    // process
    std::set<Index> localBorderIndices_;

    // stores the set of process ranks which are in the overlap for a
    // given row index "owned" by the current rank. The second value
    // store the distance from the nearest process border.
    OverlapByIndex foreignOverlapByLocalIndex_;

    // stores a list of foreign overlap indices for each rank
    OverlapByRank foreignOverlapByRank_;

    // size of the overlap region
    unsigned overlapSize_;

    // number of local indices
    size_t numLocal_;

    // number of native indices
    size_t numNative_;

    // the MPI rank of the local process
    ProcessRank myRank_;
};

} // namespace Linear
} // namespace Opm

#endif
