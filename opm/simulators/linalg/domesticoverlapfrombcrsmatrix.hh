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
 * \copydoc Opm::Linear::DomesticOverlapFromBCRSMatrix
 */
#ifndef EWOMS_DOMESTIC_OVERLAP_FROM_BCRS_MATRIX_HH
#define EWOMS_DOMESTIC_OVERLAP_FROM_BCRS_MATRIX_HH

#include "foreignoverlapfrombcrsmatrix.hh"
#include "blacklist.hh"
#include "globalindices.hh"

#include <opm/models/parallel/mpibuffer.hh>

#include <algorithm>
#include <limits>
#include <set>
#include <map>
#include <vector>

namespace Opm {
namespace Linear {

/*!
 * \brief This class creates and manages the foreign overlap given an
 *        initial list of border indices and a BCRS matrix.
 *
 * The foreign overlap are all (row) indices which overlap with the
 * some of the current process's local indices.
 */
class DomesticOverlapFromBCRSMatrix
{
    using ForeignOverlap = Opm::Linear::ForeignOverlapFromBCRSMatrix;
    using GlobalIndices = Opm::Linear::GlobalIndices<ForeignOverlap>;

public:
    // overlaps should never be copied!
    DomesticOverlapFromBCRSMatrix(const DomesticOverlapFromBCRSMatrix&) = delete;

    /*!
     * \brief Constructs the foreign overlap given a BCRS matrix and
     *        an initial list of border indices.
     */
    template <class BCRSMatrix>
    DomesticOverlapFromBCRSMatrix(const BCRSMatrix& A,
                                  const BorderList& borderList,
                                  const BlackList& blackList,
                                  unsigned overlapSize)
        : foreignOverlap_(A, borderList, blackList, overlapSize)
        , blackList_(blackList)
        , globalIndices_(foreignOverlap_)
    {
        myRank_ = 0;
        worldSize_ = 1;

#if HAVE_MPI
        int tmp;
        MPI_Comm_rank(MPI_COMM_WORLD, &tmp);
        myRank_ = static_cast<ProcessRank>(tmp);
        MPI_Comm_size(MPI_COMM_WORLD, &tmp);
        worldSize_ = static_cast<unsigned>(tmp);
#endif // HAVE_MPI

        buildDomesticOverlap_();
        updateMasterRanks_();
        blackList_.updateNativeToDomesticMap(*this);

        setupDebugMapping_();
    }

    void check() const
    {
#ifndef NDEBUG
        // check consistency of global indices
        for (unsigned domIdx = 0; domIdx < numDomestic(); ++domIdx) {
            assert(globalToDomestic(domesticToGlobal(domIdx)) == static_cast<Index>(domIdx));
        }
#endif // NDEBUG

        // send the foreign overlap for which we are master to the
        // peers
        std::map<int, MpiBuffer<unsigned> *> sizeBufferMap;

        auto peerIt = peerSet_.begin();
        const auto& peerEndIt = peerSet_.end();
        for (; peerIt != peerEndIt; ++peerIt) {
            auto& buffer = *(new MpiBuffer<unsigned>(1));
            sizeBufferMap[*peerIt] = &buffer;
            buffer[0] = foreignOverlap_.foreignOverlapWithPeer(*peerIt).size();
            buffer.send(*peerIt);
        }

        peerIt = peerSet_.begin();
        for (; peerIt != peerEndIt; ++peerIt) {
            MpiBuffer<unsigned> rcvBuffer(1);
            rcvBuffer.receive(*peerIt);

            assert(rcvBuffer[0] == domesticOverlapWithPeer_.find(*peerIt)->second.size());
        }

        peerIt = peerSet_.begin();
        for (; peerIt != peerEndIt; ++peerIt) {
            sizeBufferMap[*peerIt]->wait();
            delete sizeBufferMap[*peerIt];
        }
    }

    /*!
     * \brief Returns the rank of the current process.
     */
    ProcessRank myRank() const
    { return myRank_; }

    /*!
     * \brief Returns the number of processes in the global MPI communicator.
     */
    unsigned worldSize() const
    { return worldSize_; }

    /*!
     * \brief Return the set of process ranks which share an overlap
     *        with the current process.
     */
    const PeerSet& peerSet() const
    { return peerSet_; }

    /*!
     * \brief Returns true iff a domestic index is a border index.
     */
    bool isBorder(Index domesticIdx) const
    {
        return isLocal(domesticIdx)
            && foreignOverlap_.isBorder(mapExternalToInternal_(domesticIdx));
    }

    /*!
     * \brief Returns true iff a domestic index is on the border with
     *        a given peer process.
     */
    bool isBorderWith(Index domesticIdx, ProcessRank peerRank) const
    {
        return isLocal(domesticIdx)
            && foreignOverlap_.isBorderWith(mapExternalToInternal_(domesticIdx),
                                            peerRank);
    }

    /*!
     * \brief Returns the number of indices on the front within a given
     *        peer rank's grid partition.
     */
    size_t numFront(ProcessRank peerRank) const
    { return foreignOverlap_.numFront(peerRank); }

    /*!
     * \brief Returns true iff a domestic index is on the front.
     */
    bool isFront(Index domesticIdx) const
    {
        if (isLocal(domesticIdx))
            return false;
        Index internalDomesticIdx = mapExternalToInternal_(domesticIdx);

        // check wether the border distance of the domestic overlap is
        // maximal for the index
        const auto& domOverlap = domesticOverlapByIndex_[internalDomesticIdx];
        return domOverlap.size() > 0
            && domOverlap.begin()->second == foreignOverlap_.overlapSize();
    }

    /*!
     * \brief Returns the object which represents the black-listed native indices.
     */
    const BlackList& blackList() const
    { return blackList_; }

    /*!
     * \brief Returns the number of processes which "see" a given
     *        index.
     */
    size_t numPeers(Index domesticIdx) const
    { return domesticOverlapByIndex_[mapExternalToInternal_(domesticIdx)].size(); }

    /*!
     * \brief Returns the size of the overlap region
     */
    unsigned overlapSize() const
    { return foreignOverlap_.overlapSize(); }

    /*!
     * \brief Returns the number native indices
     *
     * I.e. the number of indices of the "raw" grid partition of the
     * local process (including the indices in ghost and overlap
     * elements).
     */
    size_t numNative() const
    { return foreignOverlap_.numNative(); }

    /*!
     * \brief Returns the number local indices
     *
     * I.e. indices in the interior or on the border of the process'
     * domain.
     */
    size_t numLocal() const
    { return foreignOverlap_.numLocal(); }

    /*!
     * \brief Returns the number domestic indices.
     *
     * The domestic indices are defined as the process' local indices
     * plus its domestic overlap (i.e. indices for which it is not
     * neither master nor are on the process border).
     */
    size_t numDomestic() const
    { return globalIndices_.numDomestic(); }

    /*!
     * \brief Return true if a domestic index is local for the process
     *
     * I.e. the entity for this index is in the interior or on the
     * border of the process' domain.
     */
    bool isLocal(Index domesticIdx) const
    { return mapExternalToInternal_(domesticIdx) < static_cast<Index>(numLocal()); }

    /*!
     * \brief Return true iff the current process is the master of a
     *        given domestic index.
     */
    bool iAmMasterOf(Index domesticIdx) const
    {
        if (!isLocal(domesticIdx))
            return false;
        return foreignOverlap_.iAmMasterOf(mapExternalToInternal_(domesticIdx));
    }

    /*!
     * \brief Return the rank of a master process for a domestic index
     */
    ProcessRank masterRank(Index domesticIdx) const
    { return masterRank_[static_cast<unsigned>(mapExternalToInternal_(domesticIdx))]; }

    /*!
     * \brief Print the foreign overlap for debugging purposes.
     */
    void print() const
    { globalIndices_.print(); }

    /*!
     * \brief Returns a domestic index given a global one
     */
    Index globalToDomestic(Index globalIdx) const
    {
        Index internalIdx = globalIndices_.globalToDomestic(globalIdx);
        if (internalIdx < 0)
            return -1;
        return mapInternalToExternal_(internalIdx);
    }

    /*!
     * \brief Returns a global index given a domestic one
     */
    Index domesticToGlobal(Index domIdx) const
    { return globalIndices_.domesticToGlobal(mapExternalToInternal_(domIdx)); }

    /*!
     * \brief Returns a native index given a domestic one
     */
    Index domesticToNative(Index domIdx) const
    {
        Index internalIdx = mapExternalToInternal_(domIdx);
        if (internalIdx >= static_cast<Index>(numLocal()))
            return -1;
        return foreignOverlap_.localToNative(internalIdx);
    }

    /*!
     * \brief Returns a domestic index given a native one
     */
    Index nativeToDomestic(Index nativeIdx) const
    {
        Index localIdx = foreignOverlap_.nativeToLocal(nativeIdx);
        if (localIdx < 0)
            return localIdx;
        return mapInternalToExternal_(localIdx);
    }

    /*!
     * \brief Returns true if a given domestic index is either in the
     *        foreign or in the domestic overlap.
     */
    bool isInOverlap(Index domesticIdx) const
    {
        return !this->isLocal(domesticIdx)
               || this->foreignOverlap_.isInOverlap(mapExternalToInternal_(domesticIdx));
    }

    /*!
     * \brief Returns true if a given domestic index is a front index
     *        for a peer rank.
     */
    bool isFrontFor(ProcessRank peerRank, Index domesticIdx) const
    {
        Index internalIdx = mapExternalToInternal_(domesticIdx);
        return this->foreignOverlap_.isFrontFor(peerRank, internalIdx);
    }

    /*!
     * \brief Returns true iff a domestic index is seen by a peer rank.
     */
    bool peerHasIndex(int peerRank, Index domesticIdx) const
    {
        return foreignOverlap_.peerHasIndex(peerRank,
                                            mapExternalToInternal_(domesticIdx));
    }

    /*!
     * \brief Returns number of indices which are contained in the
     *        foreign overlap with a peer.
     */
    size_t foreignOverlapSize(ProcessRank peerRank) const
    { return foreignOverlap_.foreignOverlapWithPeer(peerRank).size(); }

    /*!
     * \brief Returns the domestic index given an offset in the
     *        foreign overlap of a peer process with the local
     *        process.
     */
    Index foreignOverlapOffsetToDomesticIdx(ProcessRank peerRank, unsigned overlapOffset) const
    {
        Index internalIdx =
            foreignOverlap_.foreignOverlapWithPeer(peerRank)[overlapOffset].index;
        return mapInternalToExternal_(internalIdx);
    }

    /*!
     * \brief Returns number of indices which are contained in the
     *        domestic overlap with a peer.
     */
    size_t domesticOverlapSize(ProcessRank peerRank) const
    { return domesticOverlapWithPeer_.at(peerRank).size(); }

    /*!
     * \brief Returns the domestic index given an offset in the
     *        domestic overlap of a peer process with the local
     *        process.
     */
    Index domesticOverlapOffsetToDomesticIdx(ProcessRank peerRank, Index overlapOffset) const
    {
        Index internalIdx = domesticOverlapWithPeer_.at(peerRank)[overlapOffset];
        return mapInternalToExternal_(internalIdx);
    }

protected:
    void buildDomesticOverlap_()
    {
        // copy the set of peers from the foreign overlap
        peerSet_ = foreignOverlap_.peerSet();

        // resize the array which stores the number of peers for
        // each entry.
        domesticOverlapByIndex_.resize(numLocal());
        borderDistance_.resize(numLocal(), 0);

        PeerSet::const_iterator peerIt;
        PeerSet::const_iterator peerEndIt = peerSet_.end();

        // send the overlap indices to all peer processes
        peerIt = peerSet_.begin();
        for (; peerIt != peerEndIt; ++peerIt) {
            ProcessRank peerRank = *peerIt;
            sendIndicesToPeer_(peerRank);
        }

        // receive our overlap from the processes to all peer processes
        peerIt = peerSet_.begin();
        for (; peerIt != peerEndIt; ++peerIt) {
            ProcessRank peerRank = *peerIt;
            receiveIndicesFromPeer_(peerRank);
        }

        // wait until all send operations complete
        peerIt = peerSet_.begin();
        for (; peerIt != peerEndIt; ++peerIt) {
            ProcessRank peerRank = *peerIt;
            waitSendIndices_(peerRank);
        }
    }

    void updateMasterRanks_()
    {
        size_t nLocal = numLocal();
        size_t nDomestic = numDomestic();
        masterRank_.resize(nDomestic);

        // take the master ranks for the local indices from the
        // foreign overlap
        for (unsigned i = 0; i < nLocal; ++i) {
            masterRank_[i] = foreignOverlap_.masterRank(static_cast<Index>(i));
        }

        // for non-local indices, initially use INT_MAX as their master
        // rank
        for (size_t i = nLocal; i < nDomestic; ++i)
            masterRank_[i] = std::numeric_limits<ProcessRank>::max();

        // for the non-local indices, take the peer process for which
        // a given local index is in the interior
        auto peerIt = peerSet_.begin();
        const auto& peerEndIt = peerSet_.end();
        for (; peerIt != peerEndIt; ++peerIt) {
            const auto& overlapWithPeer = domesticOverlapWithPeer_.find(*peerIt)->second;

            auto idxIt = overlapWithPeer.begin();
            const auto& idxEndIt = overlapWithPeer.end();
            for (; idxIt != idxEndIt; ++idxIt) {
                if (*idxIt >= 0 && foreignOverlap_.isLocal(*idxIt))
                    continue; // ignore border indices

                masterRank_[static_cast<unsigned>(*idxIt)] = std::min(masterRank_[static_cast<unsigned>(*idxIt)], *peerIt);
            }
        }
    }

    void sendIndicesToPeer_([[maybe_unused]] ProcessRank peerRank)
    {
#if HAVE_MPI
        const auto& foreignOverlap = foreignOverlap_.foreignOverlapWithPeer(peerRank);

        // first, send a message containing the number of additional
        // indices stemming from the overlap (i.e. without the border
        // indices)
        size_t numIndices = foreignOverlap.size();
        numIndicesSendBuffer_[peerRank] = new MpiBuffer<size_t>(1);
        (*numIndicesSendBuffer_[peerRank])[0] = numIndices;
        numIndicesSendBuffer_[peerRank]->send(peerRank);

        // create MPI buffers
        indicesSendBuffer_[peerRank] = new MpiBuffer<IndexDistanceNpeers>(numIndices);

        // then send the additional indices themselfs
        auto overlapIt = foreignOverlap.begin();
        const auto& overlapEndIt = foreignOverlap.end();
        for (unsigned i = 0; overlapIt != overlapEndIt; ++overlapIt, ++i) {
            Index localIdx = overlapIt->index;
            BorderDistance borderDistance = overlapIt->borderDistance;
            size_t numPeers = foreignOverlap_.foreignOverlapByLocalIndex(localIdx).size();

            IndexDistanceNpeers tmp;
            tmp.index = globalIndices_.domesticToGlobal(localIdx);
            tmp.borderDistance = borderDistance;
            tmp.numPeers = static_cast<unsigned>(numPeers);

            (*indicesSendBuffer_[peerRank])[i] = tmp;
        }

        indicesSendBuffer_[peerRank]->send(peerRank);
#endif // HAVE_MPI
    }

    void waitSendIndices_(ProcessRank peerRank)
    {
        numIndicesSendBuffer_[peerRank]->wait();
        delete numIndicesSendBuffer_[peerRank];

        indicesSendBuffer_[peerRank]->wait();
        delete indicesSendBuffer_[peerRank];
    }

    void receiveIndicesFromPeer_([[maybe_unused]] ProcessRank peerRank)
    {
#if HAVE_MPI
        // receive the number of additional indices
        int numIndices = -1;
        MpiBuffer<size_t> numIndicesRecvBuff(1);
        numIndicesRecvBuff.receive(peerRank);
        numIndices = static_cast<int>(numIndicesRecvBuff[0]);

        // receive the additional indices themselfs
        MpiBuffer<IndexDistanceNpeers> recvBuff(static_cast<size_t>(numIndices));
        recvBuff.receive(peerRank);
        for (unsigned i = 0; i < static_cast<unsigned>(numIndices); ++i) {
            Index globalIdx = recvBuff[i].index;
            BorderDistance borderDistance = recvBuff[i].borderDistance;

            // if the index is not already known, add it to the
            // domestic indices
            if (!globalIndices_.hasGlobalIndex(globalIdx)) {
                Index newDomesticIdx = static_cast<Index>(globalIndices_.numDomestic());
                globalIndices_.addIndex(newDomesticIdx, globalIdx);

                size_t newSize = globalIndices_.numDomestic();
                borderDistance_.resize(newSize, std::numeric_limits<int>::max());
                domesticOverlapByIndex_.resize(newSize);
            }

            // convert the global index into a domestic one
            Index domesticIdx = globalIndices_.globalToDomestic(globalIdx);

            // extend the domestic overlap
            domesticOverlapByIndex_[static_cast<unsigned>(domesticIdx)][static_cast<unsigned>(peerRank)] = borderDistance;
            domesticOverlapWithPeer_[static_cast<unsigned>(peerRank)].push_back(domesticIdx);

            //assert(borderDistance >= 0);
            assert(globalIdx >= 0);
            assert(domesticIdx >= 0);
            assert(!(borderDistance == 0 && !foreignOverlap_.isLocal(domesticIdx)));
            assert(!(borderDistance > 0 && foreignOverlap_.isLocal(domesticIdx)));

            borderDistance_[static_cast<unsigned>(domesticIdx)] = std::min(borderDistance, borderDistance_[static_cast<unsigned>(domesticIdx)]);
        }
#endif // HAVE_MPI
    }

    // this method is intended to set up the code mapping code for
    // mapping domestic indices to the same ones used by a sequential
    // grid. this requires detailed knowledge about how a grid
    // distributes the degrees of freedom over multiple processes, but
    // it can simplify debugging considerably because the indices can
    // be made identical for the parallel and the sequential
    // computations.
    //
    // by default, this method does nothing
    void setupDebugMapping_()
    {}

    // this method is intended to map domestic indices to the ones
    // used by a sequential grid.
    //
    // by default, this method does nothing
    Index mapInternalToExternal_(Index internalIdx) const
    { return internalIdx; }

    // this method is intended to map the indices used by a sequential
    // to grid domestic indices ones.
    //
    // by default, this method does nothing
    Index mapExternalToInternal_(Index externalIdx) const
    { return externalIdx; }

    ProcessRank myRank_;
    unsigned worldSize_;
    ForeignOverlap foreignOverlap_;

    BlackList blackList_;

    DomesticOverlapByRank domesticOverlapWithPeer_;
    OverlapByIndex domesticOverlapByIndex_;
    std::vector<BorderDistance> borderDistance_;
    std::vector<ProcessRank> masterRank_;

    std::map<ProcessRank, MpiBuffer<size_t> *> numIndicesSendBuffer_;
    std::map<ProcessRank, MpiBuffer<IndexDistanceNpeers> *> indicesSendBuffer_;
    GlobalIndices globalIndices_;
    PeerSet peerSet_;
};

} // namespace Linear
} // namespace Opm

#endif
