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

#include <config.h>
#include <ebos/collecttoiorank.hh>

#include <opm/grid/common/CartesianIndexMapper.hpp>
#include <opm/grid/CpGrid.hpp>
#include <opm/grid/polyhedralgrid.hh>

#include <dune/grid/common/gridenums.hh>
#include <dune/grid/common/mcmgmapper.hh>

#include <cassert>
#include <stdexcept>
#include <string>

namespace Opm {

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

using IndexMapType = std::vector<int>;
using IndexMapStorageType = std::vector<IndexMapType>;
using P2PCommunicatorType = Dune::Point2PointCommunicator<Dune::SimpleMessageBuffer>;
using MessageBufferType = typename P2PCommunicatorType::MessageBufferType;
class DistributeIndexMapping : public P2PCommunicatorType::DataHandleInterface
{
protected:
    const std::vector<int>& distributedGlobalIndex_;
    IndexMapType& localIndexMap_;
    IndexMapStorageType& indexMaps_;
    std::map<int, int> globalPosition_;
    std::set<int>& recv_;
    std::vector<int>& ranks_;

public:
    DistributeIndexMapping(const std::vector<int>& globalIndex,
                           const std::vector<int>& distributedGlobalIndex,
                           IndexMapType& localIndexMap,
                           IndexMapStorageType& indexMaps,
                           std::vector<int>& ranks,
                           std::set<int>& recv,
                           bool isIORank)
    : distributedGlobalIndex_(distributedGlobalIndex)
    , localIndexMap_(localIndexMap)
    , indexMaps_(indexMaps)
    , globalPosition_()
    , recv_(recv)
    , ranks_(ranks)
    {
        size_t size = globalIndex.size();
        // create mapping globalIndex --> localIndex
        if ( isIORank ) // ioRank
            for (size_t index = 0; index < size; ++index)
                globalPosition_.insert(std::make_pair(globalIndex[index], index));

        // we need to create a mapping from local to global
        if (!indexMaps_.empty()) {
            if (isIORank)
                ranks_.resize(size, -1);
            // for the ioRank create a localIndex to index in global state map
            IndexMapType& indexMap = indexMaps_.back();
            size_t localSize = localIndexMap_.size();
            indexMap.resize(localSize);
            for (size_t i=0; i<localSize; ++i)
            {
                int id = distributedGlobalIndex_[localIndexMap_[i]];
                indexMap[i] = id;
            }
        }
    }

    ~DistributeIndexMapping()
    {
        if (!indexMaps_.size())
            return;

        if ( ranks_.size() )
        {
            auto rankIt = recv_.begin();
            // translate index maps from global cartesian to index
            for (auto& indexMap: indexMaps_)
            {
                int rank = 0;
                if (rankIt != recv_.end())
                    rank = *rankIt;

                for (auto&& entry: indexMap)
                {
                    auto candidate = globalPosition_.find(entry);
                    assert(candidate != globalPosition_.end());
                    entry = candidate->second;
                    // Using max should be backwards compatible
                    ranks_[entry] = std::max(ranks_[entry], rank);
                }
                if (rankIt != recv_.end())
                    ++rankIt;
            }
    #ifndef NDEBUG
            for (const auto& rank: ranks_)
                assert(rank>=0);
    #endif
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
            buffer.read(indexMap[index]);
        }
    }
};


/// \brief Communication handle to scatter the global index
template<class EquilMapper, class Mapper>
class ElementIndexScatterHandle
{
public:
    ElementIndexScatterHandle(const EquilMapper& sendMapper, const Mapper& recvMapper, std::vector<int>& elementIndices)
        : sendMapper_(sendMapper), recvMapper_(recvMapper), elementIndices_(elementIndices)
    {}
    using DataType = int;
    bool fixedsize(int /*dim*/, int /*codim*/)
    {
        return true;
    }

    template<class T>
    std::size_t size(const T&)
    {
        return 1;
    }
    template<class B, class T>
    void gather(B& buffer, const T& t)
    {
        buffer.write(sendMapper_.index(t));
    }
    template<class B, class T>
    void scatter(B& buffer, const T& t, std::size_t)
    {
        buffer.read(elementIndices_[recvMapper_.index(t)]);
    }

    bool contains(int dim, int codim)
    {
        return dim==3 && codim==0;
    }
private:
    const EquilMapper& sendMapper_;
    const Mapper& recvMapper_;
    std::vector<int>& elementIndices_;
};

/// \brief Communication handle to scatter the global index
template<class Mapper>
class ElementIndexHandle
{
public:
    ElementIndexHandle(const Mapper& mapper, std::vector<int>& elementIndices)
        : mapper_(mapper), elementIndices_(elementIndices)
    {}
    using DataType = int;
    bool fixedsize(int /*dim*/, int /*codim*/)
    {
        return true;
    }

    template<class T>
    std::size_t size(const T&)
    {
        return 1;
    }
    template<class B, class T>
    void gather(B& buffer, const T& t)
    {
        buffer.write(elementIndices_[mapper_.index(t)]);
    }
    template<class B, class T>
    void scatter(B& buffer, const T& t, std::size_t)
    {
        buffer.read(elementIndices_[mapper_.index(t)]);
    }

    bool contains(int dim, int codim)
    {
        return dim==3 && codim==0;
    }
private:
    const Mapper& mapper_;
    std::vector<int>& elementIndices_;
};

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
                [[maybe_unused]] auto ret = globalCellData_.insert(key, pair.second.dim,
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

class PackUnPackGroupAndNetworkValues : public P2PCommunicatorType::DataHandleInterface
{
    const Opm::data::GroupAndNetworkValues& localGroupAndNetworkData_;
    Opm::data::GroupAndNetworkValues&       globalGroupAndNetworkData_;

public:
    PackUnPackGroupAndNetworkValues(const Opm::data::GroupAndNetworkValues& localGroupAndNetworkData,
                                    Opm::data::GroupAndNetworkValues&       globalGroupAndNetworkData,
                                    const bool                              isIORank)
        : localGroupAndNetworkData_ (localGroupAndNetworkData)
        , globalGroupAndNetworkData_(globalGroupAndNetworkData)
    {
        if (! isIORank) { return; }

        MessageBufferType buffer;
        this->pack(0, buffer);

        // pass a dummy_link to satisfy virtual class
        const int dummyLink = -1;
        this->unpack(dummyLink, buffer);
    }

    // pack all data associated with link
    void pack(int link, MessageBufferType& buffer)
    {
        // we should only get one link
        if (link != 0) {
            throw std::logic_error {
                "link in method pack is not 0 as expected"
            };
        }

        // write all group and network (node/branch) data
        this->localGroupAndNetworkData_.write(buffer);
    }

    // unpack all data associated with link
    void unpack(int /*link*/, MessageBufferType& buffer)
    { this->globalGroupAndNetworkData_.read(buffer); }
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

class PackUnPackWBPData : public P2PCommunicatorType::DataHandleInterface
{
    const std::map<std::size_t, double>& localWBPData_;
    std::map<std::size_t, double>& globalWBPValues_;

public:
    PackUnPackWBPData(const std::map<std::size_t, double>& localWBPData,
                      std::map<std::size_t, double>& globalWBPValues,
                      bool isIORank)
        : localWBPData_(localWBPData)
        , globalWBPValues_(globalWBPValues)
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
        unsigned int size = localWBPData_.size();
        buffer.write(size);
        for (const auto& [global_index, wbp_value] : localWBPData_) {
            buffer.write(global_index);
            buffer.write(wbp_value);
        }
    }

    // unpack all data associated with link
    void unpack(int /*link*/, MessageBufferType& buffer)
    {
        // read all block data
        unsigned int size = 0;
        buffer.read(size);
        for (size_t i = 0; i < size; ++i) {
            std::size_t idx;
            double data;
            buffer.read(idx);
            buffer.read(data);
            globalWBPValues_[idx] = data;
        }
    }
};

class PackUnPackAquiferData : public P2PCommunicatorType::DataHandleInterface
{
    const Opm::data::Aquifers& localAquiferData_;
    data::Aquifers& globalAquiferData_;

public:
    PackUnPackAquiferData(const Opm::data::Aquifers& localAquiferData,
                          Opm::data::Aquifers& globalAquiferData,
                          bool isIORank)
        : localAquiferData_(localAquiferData)
        , globalAquiferData_(globalAquiferData)
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

        int size = localAquiferData_.size();
        buffer.write(size);
        for (const auto& [key, data] : localAquiferData_) {
            buffer.write(key);
            data.write(buffer);
        }
    }

    // unpack all data associated with link
    void unpack(int /*link*/, MessageBufferType& buffer)
    {
        int size;
        buffer.read(size);
        for (int i = 0; i < size; ++i) {
            int key;
            buffer.read(key);
            data::AquiferData data;
            data.read(buffer);

            auto& aq = this->globalAquiferData_[key];

            this->unpackCommon(data, aq);

            if (data.aquFet != nullptr) {
                this->unpackAquFet(std::move(data.aquFet), aq);
            }

            if (data.aquCT != nullptr) {
                this->unpackCarterTracy(std::move(data.aquCT), aq);
            }
        }
    }

private:
    void unpackCommon(const data::AquiferData& data, data::AquiferData& aq)
    {
        aq.aquiferID = std::max(aq.aquiferID, data.aquiferID);
        aq.pressure = std::max(aq.pressure, data.pressure);
        aq.initPressure = std::max(aq.initPressure, data.initPressure);
        aq.datumDepth = std::max(aq.datumDepth, data.datumDepth);
        aq.fluxRate += data.fluxRate;
        aq.volume += data.volume;
    }

    void unpackAquFet(std::shared_ptr<data::FetkovichData> aquFet, data::AquiferData& aq)
    {
        if (aq.aquFet == nullptr) {
            aq.aquFet = std::move(aquFet);
        }
        else {
            const auto& src = *aquFet;
            auto&       dst = *aq.aquFet;

            dst.initVolume   = std::max(dst.initVolume  , src.initVolume);
            dst.prodIndex    = std::max(dst.prodIndex   , src.prodIndex);
            dst.timeConstant = std::max(dst.timeConstant, src.timeConstant);
        }
    }

    void unpackCarterTracy(std::shared_ptr<data::CarterTracyData> aquCT, data::AquiferData& aq)
    {
        if (aq.aquCT == nullptr) {
            aq.aquCT = std::move(aquCT);
        }
        else {
            const auto& src = *aquCT;
            auto&       dst = *aq.aquCT;

            dst.dimensionless_time     = std::max(dst.dimensionless_time    , src.dimensionless_time);
            dst.dimensionless_pressure = std::max(dst.dimensionless_pressure, src.dimensionless_pressure);
        }
    }
};


template <class Grid, class EquilGrid, class GridView>
CollectDataToIORank<Grid,EquilGrid,GridView>::
CollectDataToIORank(const Grid& grid, const EquilGrid& equilGrid,
                    const GridView& localGridView,
                    const Dune::CartesianIndexMapper<Grid>& cartMapper,
                    const Dune::CartesianIndexMapper<EquilGrid>& equilCartMapper)
    : toIORankComm_()
{
    // index maps only have to be build when reordering is needed
    if (!needsReordering && !isParallel())
        return;

    const CollectiveCommunication& comm = grid.comm();

    {
        std::set<int> send, recv;
        using EquilGridView = typename EquilGrid::LeafGridView;

        typedef Dune::MultipleCodimMultipleGeomTypeMapper<GridView> ElementMapper;
        ElementMapper elemMapper(localGridView, Dune::mcmgElementLayout());
        sortedCartesianIdx_.reserve(localGridView.size(0));

        for(const auto& elem: elements(localGridView))
        {
            auto idx = elemMapper.index(elem);
            sortedCartesianIdx_.push_back(cartMapper.cartesianIndex(idx));
        }

        std::sort(sortedCartesianIdx_.begin(), sortedCartesianIdx_.end());
        localIdxToGlobalIdx_.resize(localGridView.size(0), -1);

        // the I/O rank receives from all other ranks
        if (isIORank()) {
            // We need a mapping from local to global grid, here we
            // use equilGrid which represents a view on the global grid
            // reserve memory
            const size_t globalSize = equilGrid.leafGridView().size(0);
            globalCartesianIndex_.resize(globalSize, -1);
            const EquilGridView equilGridView = equilGrid.leafGridView();

            using EquilElementMapper = Dune::MultipleCodimMultipleGeomTypeMapper<EquilGridView>;
            EquilElementMapper equilElemMapper(equilGridView, Dune::mcmgElementLayout());

            // Scatter the global index to local index for lookup during restart
            ElementIndexScatterHandle<EquilElementMapper,ElementMapper> handle(equilElemMapper, elemMapper, localIdxToGlobalIdx_);
            grid.scatterData(handle);

            // loop over all elements (global grid) and store Cartesian index
            auto elemIt = equilGrid.leafGridView().template begin<0>();
            const auto& elemEndIt = equilGrid.leafGridView().template end<0>();
            for (; elemIt != elemEndIt; ++elemIt) {
                int elemIdx = equilElemMapper.index(*elemIt);
                int cartElemIdx = equilCartMapper.cartesianIndex(elemIdx);
                globalCartesianIndex_[elemIdx] = cartElemIdx;
            }

            for (int i = 0; i < comm.size(); ++i) {
                if (i != ioRank)
                    recv.insert(i);
            }
        }
        else
        {
            // all other simply send to the I/O rank
            send.insert(ioRank);

            // Scatter the global index to local index for lookup during restart
            // This is a bit hacky since the type differs from the iorank.
            // But should work since we only receive, i.e. use the second parameter.
            ElementIndexScatterHandle<ElementMapper, ElementMapper> handle(elemMapper, elemMapper,
                                                                           localIdxToGlobalIdx_);
            grid.scatterData(handle);
        }

        // Sync the global element indices
        ElementIndexHandle<ElementMapper> handle(elemMapper, localIdxToGlobalIdx_);
        grid.communicate(handle, Dune::InteriorBorder_All_Interface,
                           Dune::ForwardCommunication);
        localIndexMap_.clear();
        const size_t gridSize = grid.size(0);
        localIndexMap_.reserve(gridSize);

        // store the local Cartesian index
        IndexMapType distributedCartesianIndex;
        distributedCartesianIndex.resize(gridSize, -1);

        // A mapping for the whole grid (including the ghosts) is needed for restarts
        auto eIt = localGridView.template begin<0>();
        const auto& eEndIt = localGridView.template end<0>();
        for (; eIt != eEndIt; ++eIt) {
            const auto element = *eIt;
            int elemIdx = elemMapper.index(element);
            distributedCartesianIndex[elemIdx] = cartMapper.cartesianIndex(elemIdx);

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
                                                globalRanks_,
                                                recv,
                                                isIORank());
        toIORankComm_.exchange(distIndexMapping);
    }
}

template <class Grid, class EquilGrid, class GridView>
void CollectDataToIORank<Grid,EquilGrid,GridView>::
collect(const Opm::data::Solution& localCellData,
        const std::map<std::pair<std::string, int>, double>& localBlockData,
        const std::map<std::size_t, double>& localWBPData,
        const Opm::data::Wells& localWellData,
        const Opm::data::GroupAndNetworkValues& localGroupAndNetworkData,
        const Opm::data::Aquifers& localAquiferData)
{
    globalCellData_ = {};
    globalBlockData_.clear();
    globalWBPData_.clear();
    globalWellData_.clear();
    globalGroupAndNetworkData_.clear();
    globalAquiferData_.clear();

    // index maps only have to be build when reordering is needed
    if(!needsReordering && !isParallel())
        return;

    // this also linearises the local buffers on ioRank
    PackUnPackCellData packUnpackCellData {
        localCellData,
        this->globalCellData_,
        this->localIndexMap_,
        this->indexMaps_,
        this->numCells(),
        this->isIORank()
    };

    if (! isParallel()) {
        // no need to collect anything.
        return;
    }

    PackUnPackWellData packUnpackWellData {
        localWellData,
                this->globalWellData_,
                this->isIORank()
    };

    PackUnPackGroupAndNetworkValues packUnpackGroupAndNetworkData {
        localGroupAndNetworkData,
                this->globalGroupAndNetworkData_,
                this->isIORank()
    };

    PackUnPackBlockData packUnpackBlockData {
        localBlockData,
                this->globalBlockData_,
                this->isIORank()
    };

    PackUnPackWBPData packUnpackWBPData {
        localWBPData,
                this->globalWBPData_,
                this->isIORank()
    };

    PackUnPackAquiferData packUnpackAquiferData {
        localAquiferData,
                this->globalAquiferData_,
                this->isIORank()
    };

    toIORankComm_.exchange(packUnpackCellData);
    toIORankComm_.exchange(packUnpackWellData);
    toIORankComm_.exchange(packUnpackGroupAndNetworkData);
    toIORankComm_.exchange(packUnpackBlockData);
    toIORankComm_.exchange(packUnpackWBPData);
    toIORankComm_.exchange(packUnpackAquiferData);

#ifndef NDEBUG
    // make sure every process is on the same page
    toIORankComm_.barrier();
#endif
}

template <class Grid, class EquilGrid, class GridView>
int CollectDataToIORank<Grid,EquilGrid,GridView>::
localIdxToGlobalIdx(unsigned localIdx) const
{
    if (!isParallel())
        return localIdx;

    if (localIdxToGlobalIdx_.empty())
        throw std::logic_error("index map is not created on this rank");

    if (localIdx > localIdxToGlobalIdx_.size())
        throw std::logic_error("local index is outside map range");

    return localIdxToGlobalIdx_[localIdx];
}

template <class Grid, class EquilGrid, class GridView>
bool CollectDataToIORank<Grid,EquilGrid,GridView>::
isCartIdxOnThisRank(int cartIdx) const
{
    if (!isParallel())
        return true;

    assert(!needsReordering);
    auto candidate = std::lower_bound(sortedCartesianIdx_.begin(), sortedCartesianIdx_.end(), cartIdx);
    return (candidate != sortedCartesianIdx_.end() && *candidate == cartIdx);
}



template class CollectDataToIORank<Dune::CpGrid,
                                   Dune::CpGrid,
                                   Dune::GridView<Dune::DefaultLeafGridViewTraits<Dune::CpGrid>>>;
template class CollectDataToIORank<Dune::PolyhedralGrid<3,3,double>,
                                   Dune::PolyhedralGrid<3,3,double>,
                                   Dune::GridView<Dune::PolyhedralGridViewTraits<3,3,double,Dune::PartitionIteratorType(4)>>>;

} // end namespace Opm
