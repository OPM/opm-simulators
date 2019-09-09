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
 * \copydoc Opm::Linear::GlobalIndices
 */
#ifndef EWOMS_GLOBAL_INDICES_HH
#define EWOMS_GLOBAL_INDICES_HH

#include <dune/grid/common/datahandleif.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/scalarproducts.hh>
#include <dune/istl/operators.hh>

#include <algorithm>
#include <set>
#include <map>
#include <iostream>
#include <tuple>

#if HAVE_MPI
#include <mpi.h>
#endif

#include "overlaptypes.hh"

namespace Opm {
namespace Linear {
/*!
 * \brief This class maps domestic row indices to and from "global"
 *        indices which is used to construct an algebraic overlap
 *        for the parallel linear solvers.
 */
template <class ForeignOverlap>
class GlobalIndices
{
    GlobalIndices(const GlobalIndices& ) = delete;

    typedef std::map<Index, Index> GlobalToDomesticMap;
    typedef std::map<Index, Index> DomesticToGlobalMap;

public:
    GlobalIndices(const ForeignOverlap& foreignOverlap)
        : foreignOverlap_(foreignOverlap)
    {
        myRank_ = 0;
        mpiSize_ = 1;

#if HAVE_MPI
        {
            int tmp;
            MPI_Comm_rank(MPI_COMM_WORLD, &tmp);
            myRank_ = static_cast<ProcessRank>(tmp);
            MPI_Comm_size(MPI_COMM_WORLD, &tmp);
            mpiSize_ = static_cast<size_t>(tmp);
        }
#endif

        // calculate the domestic overlap (i.e. all overlap indices in
        // foreign processes which the current process overlaps.)
        // This requires communication via MPI.
        buildGlobalIndices_();
    }

    /*!
     * \brief Converts a domestic index to a global one.
     */
    Index domesticToGlobal(Index domesticIdx) const
    {
        assert(domesticToGlobal_.find(domesticIdx) != domesticToGlobal_.end());

        return domesticToGlobal_.find(domesticIdx)->second;
    }

    /*!
     * \brief Converts a global index to a domestic one.
     */
    Index globalToDomestic(Index globalIdx) const
    {
        const auto& tmp = globalToDomestic_.find(globalIdx);

        if (tmp == globalToDomestic_.end())
            return -1;

        return tmp->second;
    }

    /*!
     * \brief Returns the number of indices which are in the interior or
     *        on the border of the current rank.
     */
    size_t numLocal() const
    { return foreignOverlap_.numLocal(); }

    /*!
     * \brief Returns the number domestic indices.
     *
     * The domestic indices are defined as the process' local indices
     * plus its copies of indices in the overlap regions
     */
    size_t numDomestic() const
    { return numDomestic_; }

    /*!
     * \brief Add an index to the domestic<->global mapping.
     */
    void addIndex(Index domesticIdx, Index globalIdx)
    {
        domesticToGlobal_[domesticIdx] = globalIdx;
        globalToDomestic_[globalIdx] = domesticIdx;
        numDomestic_ = domesticToGlobal_.size();

        assert(domesticToGlobal_.size() == globalToDomestic_.size());
    }

    /*!
     * \brief Send a border index to a remote process.
     */
    void sendBorderIndex(ProcessRank peerRank OPM_UNUSED_NOMPI,
                         Index domesticIdx OPM_UNUSED_NOMPI,
                         Index peerLocalIdx OPM_UNUSED_NOMPI)
    {
#if HAVE_MPI
        PeerIndexGlobalIndex sendBuf;
        sendBuf.peerIdx = peerLocalIdx;
        sendBuf.globalIdx = domesticToGlobal(domesticIdx);
        MPI_Send(&sendBuf,                     // buff
                 sizeof(PeerIndexGlobalIndex), // count
                 MPI_BYTE,                     // data type
                 static_cast<int>(peerRank),   // peer process
                 0,                            // tag
                 MPI_COMM_WORLD);              // communicator
#endif
    }

    /*!
     * \brief Receive an index on the border from a remote
     *        process and add it the translation maps.
     */
    void receiveBorderIndex(ProcessRank peerRank OPM_UNUSED_NOMPI)
    {
#if HAVE_MPI
        PeerIndexGlobalIndex recvBuf;
        MPI_Recv(&recvBuf,                     // buff
                 sizeof(PeerIndexGlobalIndex), // count
                 MPI_BYTE,                     // data type
                 static_cast<int>(peerRank),   // peer process
                 0,                            // tag
                 MPI_COMM_WORLD,               // communicator
                 MPI_STATUS_IGNORE);           // status

        Index domesticIdx = foreignOverlap_.nativeToLocal(recvBuf.peerIdx);
        if (domesticIdx >= 0) {
            Index globalIdx = recvBuf.globalIdx;
            addIndex(domesticIdx, globalIdx);
        }
#endif // HAVE_MPI
    }

    /*!
     * \brief Return true iff a given global index already exists
     */
    bool hasGlobalIndex(Index globalIdx) const
    { return globalToDomestic_.find(globalIdx) != globalToDomestic_.end(); }

    /*!
     * \brief Prints the global indices of all domestic indices
     *        for debugging purposes.
     */
    void print() const
    {
        std::cout << "(domestic index, global index, domestic->global->domestic)"
                  << " list for rank " << myRank_ << "\n";

        for (size_t domIdx = 0; domIdx < domesticToGlobal_.size(); ++domIdx)
            std::cout << "(" << domIdx << ", " << domesticToGlobal(domIdx)
                      << ", " << globalToDomestic(domesticToGlobal(domIdx)) << ") ";
        std::cout << "\n" << std::flush;
    }

protected:
    // retrieve the offset for the indices where we are master in the
    // global index list
    void buildGlobalIndices_()
    {
#if HAVE_MPI
        numDomestic_ = 0;
#else
        numDomestic_ = foreignOverlap_.numLocal();
#endif

#if HAVE_MPI
        if (myRank_ == 0) {
            // the first rank starts at index zero
            domesticOffset_ = 0;
        }
        else {
            // all other ranks retrieve their offset from the next
            // lower rank
            MPI_Recv(&domesticOffset_, // buffer
                     1,                // count
                     MPI_INT,          // data type
                     static_cast<int>(myRank_ - 1), // peer rank
                     0,                // tag
                     MPI_COMM_WORLD,   // communicator
                     MPI_STATUS_IGNORE);
        }

        // create maps for all indices for which the current process
        // is the master
        int numMaster = 0;
        for (unsigned i = 0; i < foreignOverlap_.numLocal(); ++i) {
            if (!foreignOverlap_.iAmMasterOf(static_cast<Index>(i)))
                continue;

            addIndex(static_cast<Index>(i),
                     static_cast<Index>(domesticOffset_ + numMaster));
            ++numMaster;
        }

        if (myRank_ < mpiSize_ - 1) {
            // send the domestic offset plus the number of master
            // indices to the process which is one rank higher
            int tmp = domesticOffset_ + numMaster;
            MPI_Send(&tmp,            // buff
                     1,               // count
                     MPI_INT,         // data type
                     static_cast<int>(myRank_ + 1), // peer rank
                     0,               // tag
                     MPI_COMM_WORLD); // communicator
        }

        typename PeerSet::const_iterator peerIt;
        typename PeerSet::const_iterator peerEndIt = peerSet_().end();
        // receive the border indices from the lower ranks
        peerIt = peerSet_().begin();
        for (; peerIt != peerEndIt; ++peerIt) {
            if (*peerIt < myRank_)
                receiveBorderFrom_(*peerIt);
        }

        // send the border indices to the higher ranks
        peerIt = peerSet_().begin();
        for (; peerIt != peerEndIt; ++peerIt) {
            if (*peerIt > myRank_)
                sendBorderTo_(*peerIt);
        }

        // receive the border indices from the higher ranks
        peerIt = peerSet_().begin();
        for (; peerIt != peerEndIt; ++peerIt) {
            if (*peerIt > myRank_)
                receiveBorderFrom_(*peerIt);
        }

        // send the border indices to the lower ranks
        peerIt = peerSet_().begin();
        for (; peerIt != peerEndIt; ++peerIt) {
            if (*peerIt < myRank_)
                sendBorderTo_(*peerIt);
        }
#endif // HAVE_MPI
    }

    void sendBorderTo_(ProcessRank peerRank OPM_UNUSED_NOMPI)
    {
#if HAVE_MPI
        // send (local index on myRank, global index) pairs to the
        // peers
        BorderList::const_iterator borderIt = borderList_().begin();
        BorderList::const_iterator borderEndIt = borderList_().end();
        for (; borderIt != borderEndIt; ++borderIt) {
            ProcessRank borderPeer = borderIt->peerRank;
            BorderDistance borderDistance = borderIt->borderDistance;
            if (borderPeer != peerRank || borderDistance != 0)
                continue;

            Index localIdx = foreignOverlap_.nativeToLocal(borderIt->localIdx);
            Index peerIdx = borderIt->peerIdx;
            assert(localIdx >= 0);
            if (foreignOverlap_.iAmMasterOf(localIdx)) {
                sendBorderIndex(borderPeer, localIdx, peerIdx);
            }
        }
#endif // HAVE_MPI
    }

    void receiveBorderFrom_(ProcessRank peerRank OPM_UNUSED_NOMPI)
    {
#if HAVE_MPI
        // retrieve the global indices for which we are not master
        // from the processes with lower rank
        BorderList::const_iterator borderIt = borderList_().begin();
        BorderList::const_iterator borderEndIt = borderList_().end();
        for (; borderIt != borderEndIt; ++borderIt) {
            ProcessRank borderPeer = borderIt->peerRank;
            BorderDistance borderDistance = borderIt->borderDistance;
            if (borderPeer != peerRank || borderDistance != 0)
                continue;

            Index nativeIdx = borderIt->localIdx;
            Index localIdx = foreignOverlap_.nativeToLocal(nativeIdx);
            if (localIdx >= 0 && foreignOverlap_.masterRank(localIdx) == borderPeer)
                receiveBorderIndex(borderPeer);
        }
#endif // HAVE_MPI
    }

    const PeerSet& peerSet_() const
    { return foreignOverlap_.peerSet(); }

    const BorderList& borderList_() const
    { return foreignOverlap_.borderList(); }

    ProcessRank myRank_;
    size_t mpiSize_;

    int domesticOffset_;
    size_t numDomestic_;
    const ForeignOverlap& foreignOverlap_;

    GlobalToDomesticMap globalToDomestic_;
    DomesticToGlobalMap domesticToGlobal_;
};

} // namespace Linear
} // namespace Opm

#endif
