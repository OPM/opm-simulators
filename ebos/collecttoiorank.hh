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
#ifndef EWOMS_COLLECT_TO_IO_RANK_HH
#define EWOMS_COLLECT_TO_IO_RANK_HH

#include <opm/output/data/Cells.hpp>
#include <opm/output/data/Solution.hpp>
#include <opm/output/data/Wells.hpp>

#include <opm/grid/common/p2pcommunicator.hh>
#include <dune/grid/utility/persistentcontainer.hh>
#include <dune/grid/common/gridenums.hh>

#include <opm/parser/eclipse/EclipseState/EclipseState.hpp>
#include <opm/parser/eclipse/EclipseState/Grid/EclipseGrid.hpp>

#include <opm/material/common/Exceptions.hpp>
#include <opm/material/common/Unused.hpp>

#include <dune/grid/common/mcmgmapper.hh>

#include <stdexcept>

namespace Ewoms {

template <class Vanguard>
class CollectDataToIORank
{
public:
    typedef typename Vanguard::Grid  Grid;
    typedef typename Grid::CollectiveCommunication  CollectiveCommunication;

    // global id
    class GlobalCellIndex
    {
        int globalId_;
        int localIndex_;
        bool isInterior_;

    public:
        GlobalCellIndex()
            : globalId_(-1)
            , localIndex_(-1)
            , isInterior_(true)
        {}
        void setGhost()
        { isInterior_ = false; }

        void setId(int globalId)
        { globalId_ = globalId; }
        void setIndex(int localIndex)
        { localIndex_ = localIndex; }

        int localIndex () const
        { return localIndex_; }
        int id () const
        { return globalId_; }
        bool isInterior() const
        { return isInterior_; }
    };

    typedef typename Dune::PersistentContainer<Grid, GlobalCellIndex> GlobalIndexContainer;

    static const int dimension = Grid::dimension;

    typedef typename Grid::LeafGridView GridView;
    typedef GridView AllGridView;

    typedef Dune::Point2PointCommunicator<Dune::SimpleMessageBuffer> P2PCommunicatorType;
    typedef typename P2PCommunicatorType::MessageBufferType MessageBufferType;

    typedef std::vector<GlobalCellIndex> LocalIndexMapType;

    typedef std::vector<int> IndexMapType;
    typedef std::vector<IndexMapType> IndexMapStorageType;

    class DistributeIndexMapping : public P2PCommunicatorType::DataHandleInterface
    {
    protected:
        const std::vector<int>& distributedGlobalIndex_;
        IndexMapType& localIndexMap_;
        IndexMapStorageType& indexMaps_;
        std::map<int, int> globalPosition_;
        std::vector<int>& ranks_;

    public:
        DistributeIndexMapping(const std::vector<int>& globalIndex,
                               const std::vector<int>& distributedGlobalIndex,
                               IndexMapType& localIndexMap,
                               IndexMapStorageType& indexMaps,
                               std::vector<int>& ranks)
            : distributedGlobalIndex_(distributedGlobalIndex)
            , localIndexMap_(localIndexMap)
            , indexMaps_(indexMaps)
            , globalPosition_()
            , ranks_(ranks)
        {
            size_t size = globalIndex.size();
            // create mapping globalIndex --> localIndex
            for (size_t index = 0; index < size; ++index)
                globalPosition_.insert(std::make_pair(globalIndex[index], index));

            // we need to create a mapping from local to global
            if (!indexMaps_.empty()) {
                ranks_.resize(size, -1);
                // for the ioRank create a localIndex to index in global state map
                IndexMapType& indexMap = indexMaps_.back();
                size_t localSize = localIndexMap_.size();
                indexMap.resize(localSize);
                for (size_t i=0; i<localSize; ++i)
                {
                    int id = distributedGlobalIndex_[localIndexMap_[i]];
                    indexMap[i] = globalPosition_[id];
                    ranks_[indexMap[i]] = ioRank;
                }
            }
        }

        void pack(int link, MessageBufferType& buffer)
        {
            // we should only get one link
            if (link != 0)
                throw std::logic_error("link in method pack is not 0 as execpted");

            // pack all interior global cell id's
            int size = localIndexMap_.size();
            buffer.write(size);

            for (int index = 0; index < size; ++index) {
                int globalIdx = distributedGlobalIndex_[localIndexMap_[index]];
                buffer.write(globalIdx);
            }
        }

        void unpack(int link, MessageBufferType& buffer)
        {
            // get index map for current link
            IndexMapType& indexMap = indexMaps_[link];
            assert(!globalPosition_.empty());

            // unpack all interior global cell id's
            int numCells = 0;
            buffer.read(numCells);
            indexMap.resize(numCells);
            for (int index = 0; index < numCells; ++index) {
                int globalId = -1;
                buffer.read(globalId);
                assert(globalPosition_.find(globalId) != globalPosition_.end());
                indexMap[index] = globalPosition_[globalId];
                ranks_[indexMap[index]] = link + 1;
            }
        }
    };

    enum { ioRank = 0 };

    static const bool needsReordering =
        !std::is_same<typename Vanguard::Grid, typename Vanguard::EquilGrid>::value;

    CollectDataToIORank(const Vanguard& vanguard)
        : toIORankComm_()
    {
        // index maps only have to be build when reordering is needed
        if (!needsReordering && !isParallel())
            return;

        const CollectiveCommunication& comm = vanguard.grid().comm();

        {
            std::set<int> send, recv;
            typedef typename Vanguard::EquilGrid::LeafGridView EquilGridView;
            const EquilGridView equilGridView = vanguard.equilGrid().leafGridView();

#if DUNE_VERSION_NEWER(DUNE_GRID, 2,6)
            typedef Dune::MultipleCodimMultipleGeomTypeMapper<EquilGridView> EquilElementMapper;
            EquilElementMapper equilElemMapper(equilGridView, Dune::mcmgElementLayout());
#else
            typedef Dune::MultipleCodimMultipleGeomTypeMapper<EquilGridView, Dune::MCMGElementLayout> EquilElementMapper;
            EquilElementMapper equilElemMapper(equilGridView);
#endif

            // We need a mapping from local to global grid, here we
            // use equilGrid which represents a view on the global grid
            const size_t globalSize = vanguard.equilGrid().leafGridView().size(0);
            // reserve memory
            globalCartesianIndex_.resize(globalSize, -1);

            // loop over all elements (global grid) and store Cartesian index
            auto elemIt = vanguard.equilGrid().leafGridView().template begin<0>();
            const auto& elemEndIt = vanguard.equilGrid().leafGridView().template end<0>();
            for (; elemIt != elemEndIt; ++elemIt) {
                int elemIdx = equilElemMapper.index(*elemIt);
                int cartElemIdx = vanguard.equilCartesianIndexMapper().cartesianIndex(elemIdx);
                globalCartesianIndex_[elemIdx] = cartElemIdx;
            }

            // the I/O rank receives from all other ranks
            if (isIORank()) {
                for (int i = 0; i < comm.size(); ++i) {
                    if (i != ioRank)
                        recv.insert(i);
                }
            }
            else // all other simply send to the I/O rank
                send.insert(ioRank);

            localIndexMap_.clear();
            const size_t gridSize = vanguard.grid().size(0);
            localIndexMap_.reserve(gridSize);

            // store the local Cartesian index
            IndexMapType distributedCartesianIndex;
            distributedCartesianIndex.resize(gridSize, -1);

            typedef typename Vanguard::GridView LocalGridView;
            const LocalGridView localGridView = vanguard.gridView();

#if DUNE_VERSION_NEWER(DUNE_GRID, 2,6)
            typedef Dune::MultipleCodimMultipleGeomTypeMapper<LocalGridView> ElementMapper;
            ElementMapper elemMapper(localGridView, Dune::mcmgElementLayout());
#else
            typedef Dune::MultipleCodimMultipleGeomTypeMapper<LocalGridView, Dune::MCMGElementLayout> ElementMapper;
            ElementMapper elemMapper(localGridView);
#endif

            // A mapping for the whole grid (including the ghosts) is needed for restarts
            auto eIt = localGridView.template begin<0>();
            const auto& eEndIt = localGridView.template end<0>();
            for (; eIt != eEndIt; ++eIt) {
                const auto element = *eIt;
                int elemIdx = elemMapper.index(element);
                distributedCartesianIndex[elemIdx] = vanguard.cartesianIndex(elemIdx);

                // only store interior element for collection
                //assert(element.partitionType() == Dune::InteriorEntity);

                localIndexMap_.push_back(elemIdx);
            }

            // insert send and recv linkage to communicator
            toIORankComm_.insertRequest(send, recv);

            // need an index map for each rank
            indexMaps_.clear();
            indexMaps_.resize(comm.size());

            // distribute global id's to io rank for later association of dof's
            DistributeIndexMapping distIndexMapping(globalCartesianIndex_,
                                                    distributedCartesianIndex,
                                                    localIndexMap_,
                                                    indexMaps_,
                                                    globalRanks_);
            toIORankComm_.exchange(distIndexMapping);
        }
    }

    class PackUnPackCellData : public P2PCommunicatorType::DataHandleInterface
    {
        const Opm::data::Solution& localCellData_;
        Opm::data::Solution& globalCellData_;

        const IndexMapType& localIndexMap_;
        const IndexMapStorageType& indexMaps_;

    public:
        PackUnPackCellData(const Opm::data::Solution& localCellData,
                           Opm::data::Solution& globalCellData,
                           const IndexMapType& localIndexMap,
                           const IndexMapStorageType& indexMaps,
                           size_t globalSize,
                           bool isIORank)
            : localCellData_(localCellData)
            , globalCellData_(globalCellData)
            , localIndexMap_(localIndexMap)
            , indexMaps_(indexMaps)
        {
            if (isIORank) {
                // add missing data to global cell data
                for (const auto& pair : localCellData_) {
                    const std::string& key = pair.first;
                    std::size_t containerSize = globalSize;
                    auto OPM_OPTIM_UNUSED ret = globalCellData_.insert(key, pair.second.dim,
                                                                       std::vector<double>(containerSize),
                                                                       pair.second.target);
                    assert(ret.second);
                }

                MessageBufferType buffer;
                pack(0, buffer);

                // the last index map is the local one
                doUnpack(indexMaps.back(), buffer);
            }
        }

        // pack all data associated with link
        void pack(int link, MessageBufferType& buffer)
        {
            // we should only get one link
            if (link != 0)
                throw std::logic_error("link in method pack is not 0 as expected");

            // write all cell data registered in local state
            for (const auto& pair : localCellData_) {
                const auto& data = pair.second.data;

                // write all data from local data to buffer
                write(buffer, localIndexMap_, data);
            }
        }

        void doUnpack(const IndexMapType& indexMap, MessageBufferType& buffer)
        {
            // we loop over the data as
            // its order governs the order the data got received.
            for (auto& pair : localCellData_) {
                const std::string& key = pair.first;
                auto& data = globalCellData_.data(key);

                //write all data from local cell data to buffer
                read(buffer, indexMap, data);
            }
        }

        // unpack all data associated with link
        void unpack(int link, MessageBufferType& buffer)
        { doUnpack(indexMaps_[link], buffer); }

    protected:
        template <class Vector>
        void write(MessageBufferType& buffer,
                   const IndexMapType& localIndexMap,
                   const Vector& vector,
                   unsigned int offset = 0,
                   unsigned int stride = 1) const
        {
            unsigned int size = localIndexMap.size();
            buffer.write(size);
            assert(vector.size() >= stride * size);
            for (unsigned int i=0; i<size; ++i)
            {
                unsigned int index = localIndexMap[i] * stride + offset;
                assert(index < vector.size());
                buffer.write(vector[index]);
            }
        }

        template <class Vector>
        void read(MessageBufferType& buffer,
                  const IndexMapType& indexMap,
                  Vector& vector,
                  unsigned int offset = 0,
                  unsigned int stride = 1) const
        {
            unsigned int size = 0;
            buffer.read(size);
            assert(size == indexMap.size());
            for (unsigned int i=0; i<size; ++i) {
                unsigned int index = indexMap[i] * stride + offset;
                assert(index < vector.size());
                buffer.read(vector[index]);
            }
        }
    };

    class PackUnPackWellData : public P2PCommunicatorType::DataHandleInterface
    {
        const Opm::data::Wells& localWellData_;
        Opm::data::Wells& globalWellData_;

    public:
        PackUnPackWellData(const Opm::data::Wells& localWellData,
                           Opm::data::Wells& globalWellData,
                           bool isIORank)
            : localWellData_(localWellData)
            , globalWellData_(globalWellData)
        {
            if (isIORank) {
                MessageBufferType buffer;
                pack(0, buffer);

                // pass a dummy_link to satisfy virtual class
                int dummyLink = -1;
                unpack(dummyLink, buffer);
            }
        }

        // pack all data associated with link
        void pack(int link, MessageBufferType& buffer)
        {
            // we should only get one link
            if (link != 0)
                throw std::logic_error("link in method pack is not 0 as expected");

            localWellData_.write(buffer);
        }

        // unpack all data associated with link
        void unpack(int /*link*/, MessageBufferType& buffer)
        { globalWellData_.read(buffer); }

    };

    class PackUnPackBlockData : public P2PCommunicatorType::DataHandleInterface
    {
        const std::map<std::pair<std::string, int>, double>& localBlockData_;
        std::map<std::pair<std::string, int>, double>& globalBlockValues_;

    public:
        PackUnPackBlockData(const std::map<std::pair<std::string, int>, double>& localBlockData,
                            std::map<std::pair<std::string, int>, double>& globalBlockValues,
                            bool isIORank)
            : localBlockData_(localBlockData)
            , globalBlockValues_(globalBlockValues)
        {
            if (isIORank) {
                MessageBufferType buffer;
                pack(0, buffer);

                // pass a dummyLink to satisfy virtual class
                int dummyLink = -1;
                unpack(dummyLink, buffer);
            }
        }

        // pack all data associated with link
        void pack(int link, MessageBufferType& buffer)
        {
            // we should only get one link
            if (link != 0)
                throw std::logic_error("link in method pack is not 0 as expected");

            // write all block data
            unsigned int size = localBlockData_.size();
            buffer.write(size);
            for (const auto& map : localBlockData_) {
                buffer.write(map.first.first);
                buffer.write(map.first.second);
                buffer.write(map.second);
            }
        }

        // unpack all data associated with link
        void unpack(int /*link*/, MessageBufferType& buffer)
        {
            // read all block data
            unsigned int size = 0;
            buffer.read(size);
            for (size_t i = 0; i < size; ++i) {
                std::string name;
                int idx;
                double data;
                buffer.read(name);
                buffer.read(idx);
                buffer.read(data);
                globalBlockValues_[std::make_pair(name, idx)] = data;
            }
        }

    };

    // gather solution to rank 0 for EclipseWriter
    void collect(const Opm::data::Solution& localCellData,
                 const std::map<std::pair<std::string, int>, double>& localBlockData,
                 const Opm::data::Wells& localWellData)
    {
        globalCellData_ = {};
        globalBlockData_.clear();
        globalWellData_.clear();

        // index maps only have to be build when reordering is needed
        if(!needsReordering && !isParallel())
            return;

        // this also packs and unpacks the local buffers one ioRank
        PackUnPackCellData
            packUnpackCellData(localCellData,
                                globalCellData_,
                                localIndexMap_,
                                indexMaps_,
                                numCells(),
                                isIORank());

        if (!isParallel())
            // no need to collect anything.
            return;

        PackUnPackWellData
            packUnpackWellData(localWellData,
                               globalWellData_,
                               isIORank());

        PackUnPackBlockData
            packUnpackBlockData(localBlockData,
                                globalBlockData_,
                                isIORank());

        toIORankComm_.exchange(packUnpackCellData);
        toIORankComm_.exchange(packUnpackWellData);
        toIORankComm_.exchange(packUnpackBlockData);



#ifndef NDEBUG
        // mkae sure every process is on the same page
        toIORankComm_.barrier();
#endif
    }

    const std::map<std::pair<std::string, int>, double>& globalBlockData() const
    { return globalBlockData_; }

    const Opm::data::Solution& globalCellData() const
    { return globalCellData_; }

    const Opm::data::Wells& globalWellData() const
    { return globalWellData_; }

    bool isIORank() const
    { return toIORankComm_.rank() == ioRank; }

    bool isParallel() const
    { return toIORankComm_.size() > 1; }

    int localIdxToGlobalIdx(unsigned localIdx) const
    {
        if (!isParallel())
            return localIdx;

        // the last indexMap is the local one
        const IndexMapType& indexMap = indexMaps_.back();
        if (indexMap.empty())
            throw std::logic_error("index map is not created on this rank");

        if (localIdx > indexMap.size())
            throw std::logic_error("local index is outside map range");

        return indexMap[localIdx];
    }

    size_t numCells () const
    { return globalCartesianIndex_.size(); }

    const std::vector<int>& globalRanks() const
    { return globalRanks_; }

    bool isGlobalIdxOnThisRank(unsigned globalIdx) const
    {
        if (!isParallel())
            return true;

        // the last indexMap is the local one
        const IndexMapType& indexMap = indexMaps_.back();
        if (indexMap.empty())
            throw std::logic_error("index map is not created on this rank");

        return std::find(indexMap.begin(), indexMap.end(), globalIdx) != indexMap.end();
    }

protected:
    P2PCommunicatorType toIORankComm_;
    IndexMapType globalCartesianIndex_;
    IndexMapType localIndexMap_;
    IndexMapStorageType indexMaps_;
    std::vector<int> globalRanks_;
    Opm::data::Solution globalCellData_;
    std::map<std::pair<std::string, int>, double> globalBlockData_;
    Opm::data::Wells globalWellData_;
};

} // end namespace Ewoms

#endif
