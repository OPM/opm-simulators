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

#ifndef EWOMS_COLLECT_TO_IO_RANK_IMPL_HH
#define EWOMS_COLLECT_TO_IO_RANK_IMPL_HH

#include <ebos/collecttoiorank.hh>

#include <opm/grid/common/CartesianIndexMapper.hpp>

#include <dune/common/version.hh>
#include <dune/grid/common/gridenums.hh>
#include <dune/grid/common/mcmgmapper.hh>

#include <algorithm>
#include <cassert>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
    std::vector<std::string> toVector(const std::set<std::string>& string_set)
    {
        return { string_set.begin(), string_set.end() };
    }
}

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
        std::size_t size = globalIndex.size();
        // create mapping globalIndex --> localIndex
        if ( isIORank ) // ioRank
            for (std::size_t index = 0; index < size; ++index)
                globalPosition_.insert(std::make_pair(globalIndex[index], index));

        // we need to create a mapping from local to global
        if (!indexMaps_.empty()) {
            if (isIORank)
                ranks_.resize(size, -1);
            // for the ioRank create a localIndex to index in global state map
            IndexMapType& indexMap = indexMaps_.back();
            std::size_t localSize = localIndexMap_.size();
            indexMap.resize(localSize);
            for (std::size_t i = 0; i < localSize; ++i)
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
    bool fixedSize(int /*dim*/, int /*codim*/)
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
    bool fixedSize(int /*dim*/, int /*codim*/)
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
    const data::Solution& localCellData_;
    data::Solution& globalCellData_;

    const IndexMapType& localIndexMap_;
    const IndexMapStorageType& indexMaps_;

public:
    PackUnPackCellData(const data::Solution& localCellData,
                       data::Solution& globalCellData,
                       const IndexMapType& localIndexMap,
                       const IndexMapStorageType& indexMaps,
                       std::size_t globalSize,
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
            const auto& data = pair.second.data<double>();

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
            auto& data = globalCellData_.data<double>(key);

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
    const data::Wells& localWellData_;
    data::Wells& globalWellData_;

public:
    PackUnPackWellData(const data::Wells& localWellData,
                       data::Wells& globalWellData,
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
    const data::GroupAndNetworkValues& localGroupAndNetworkData_;
    data::GroupAndNetworkValues&       globalGroupAndNetworkData_;

public:
    PackUnPackGroupAndNetworkValues(const data::GroupAndNetworkValues& localGroupAndNetworkData,
                                    data::GroupAndNetworkValues&       globalGroupAndNetworkData,
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
        for (std::size_t i = 0; i < size; ++i) {
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
    const data::WellBlockAveragePressures& localWBPData_;
    data::WellBlockAveragePressures& globalWBPValues_;

public:
    PackUnPackWBPData(const data::WellBlockAveragePressures& localWBPData,
                      data::WellBlockAveragePressures&       globalWBPValues,
                      const bool                             isIORank)
        : localWBPData_   (localWBPData)
        , globalWBPValues_(globalWBPValues)
    {
        if (! isIORank) {
            return;
        }

        MessageBufferType buffer;
        pack(0, buffer);

        // Pass a dummy link to satisfy base class API requirement
        const int dummyLink = -1;
        unpack(dummyLink, buffer);
    }

    // Pack all data associated with link
    void pack(int link, MessageBufferType& buffer)
    {
        // We should only get one link
        if (link != 0) {
            throw std::logic_error {
                "link (" + std::to_string(link) +
                ") in method pack() is not 0 as expected"
            };
        }

        // Write all local, per-well, WBP data
        this->localWBPData_.write(buffer);
    }

    // Unpack all data associated with link
    void unpack([[maybe_unused]] const int link,
                MessageBufferType&         buffer)
    {
        this->globalWBPValues_.read(buffer);
    }
};

class PackUnPackWellTestState : public P2PCommunicatorType::DataHandleInterface
{
public:
    PackUnPackWellTestState(const WellTestState& localWTestState,
                            WellTestState& globalWTestState,
                            bool isIORank)
        : local_(localWTestState)
        , global_(globalWTestState)
    {
        if (isIORank) {
            MessageBufferType buffer;
            pack(0, buffer);

            // pass a dummyLink to satisfy virtual class
            int dummyLink = -1;
            unpack(dummyLink, buffer);
        }
    }

    void pack(int link, MessageBufferType& buffer) {
        if (link != 0)
            throw std::logic_error("link in method pack is not 0 as expected");
        this->local_.pack(buffer);
    }

    void unpack(int, MessageBufferType& buffer) {
        this->global_.unpack(buffer);
    }

private:
    const WellTestState& local_;
    WellTestState& global_;
};

class PackUnPackAquiferData : public P2PCommunicatorType::DataHandleInterface
{
    const data::Aquifers& localAquiferData_;
    data::Aquifers& globalAquiferData_;

public:
    PackUnPackAquiferData(const data::Aquifers& localAquiferData,
                          data::Aquifers& globalAquiferData,
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

            if (auto const* aquFet = data.typeData.get<data::AquiferType::Fetkovich>();
                aquFet != nullptr)
            {
                this->unpackAquFet(*aquFet, aq);
            }
            else if (auto const* aquCT = data.typeData.get<data::AquiferType::CarterTracy>();
                     aquCT != nullptr)
            {
                this->unpackCarterTracy(*aquCT, aq);
            }
            else if (auto const* aquNum = data.typeData.get<data::AquiferType::Numerical>();
                     aquNum != nullptr)
            {
                this->unpackNumericAquifer(*aquNum, aq);
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

    void unpackAquFet(const data::FetkovichData& aquFet, data::AquiferData& aq)
    {
        if (! aq.typeData.is<data::AquiferType::Fetkovich>()) {
            auto* tData = aq.typeData.create<data::AquiferType::Fetkovich>();
            *tData = aquFet;
        }
        else {
            const auto& src = aquFet;
            auto&       dst = *aq.typeData.getMutable<data::AquiferType::Fetkovich>();

            dst.initVolume   = std::max(dst.initVolume  , src.initVolume);
            dst.prodIndex    = std::max(dst.prodIndex   , src.prodIndex);
            dst.timeConstant = std::max(dst.timeConstant, src.timeConstant);
        }
    }

    void unpackCarterTracy(const data::CarterTracyData& aquCT, data::AquiferData& aq)
    {
        if (! aq.typeData.is<data::AquiferType::CarterTracy>()) {
            auto* tData = aq.typeData.create<data::AquiferType::CarterTracy>();
            *tData = aquCT;
        }
        else {
            const auto& src = aquCT;
            auto&       dst = *aq.typeData.getMutable<data::AquiferType::CarterTracy>();

            dst.timeConstant   = std::max(dst.timeConstant  , src.timeConstant);
            dst.influxConstant = std::max(dst.influxConstant, src.influxConstant);
            dst.waterDensity   = std::max(dst.waterDensity  , src.waterDensity);
            dst.waterViscosity = std::max(dst.waterViscosity, src.waterViscosity);

            dst.dimensionless_time     = std::max(dst.dimensionless_time    , src.dimensionless_time);
            dst.dimensionless_pressure = std::max(dst.dimensionless_pressure, src.dimensionless_pressure);
        }
    }

    void unpackNumericAquifer(const data::NumericAquiferData& aquNum, data::AquiferData& aq)
    {
        if (! aq.typeData.is<data::AquiferType::Numerical>()) {
            auto* tData = aq.typeData.create<data::AquiferType::Numerical>();
            *tData = aquNum;
        }
        else {
            const auto& src = aquNum;
            auto&       dst = *aq.typeData.getMutable<data::AquiferType::Numerical>();

            std::transform(src.initPressure.begin(),
                           src.initPressure.end(),
                           dst.initPressure.begin(),
                           dst.initPressure.begin(),
                           [](const double p0_1, const double p0_2)
                           {
                               return std::max(p0_1, p0_2);
                           });
        }
    }
};

class PackUnpackInterRegFlows : public P2PCommunicatorType::DataHandleInterface
{
    const EclInterRegFlowMap& localInterRegFlows_;
    EclInterRegFlowMap&       globalInterRegFlows_;

public:
    PackUnpackInterRegFlows(const EclInterRegFlowMap& localInterRegFlows,
                            EclInterRegFlowMap&       globalInterRegFlows,
                            const bool                isIORank)
        : localInterRegFlows_(localInterRegFlows)
        , globalInterRegFlows_(globalInterRegFlows)
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

        // write all inter-region flow data
        this->localInterRegFlows_.write(buffer);
    }

    // unpack all data associated with link
    void unpack(int /*link*/, MessageBufferType& buffer)
    { this->globalInterRegFlows_.read(buffer); }
};

class PackUnpackFlows : public P2PCommunicatorType::DataHandleInterface
{
    const std::array<std::pair<std::string, std::pair<std::vector<int>, std::vector<double>>>, 3>& localFlows_;
    std::array<std::pair<std::string, std::pair<std::vector<int>, std::vector<double>>>, 3>& globalFlows_;

public:
    PackUnpackFlows(const std::array<std::pair<std::string, std::pair<std::vector<int>, std::vector<double>>>, 3> & localFlows,
                     std::array<std::pair<std::string, std::pair<std::vector<int>, std::vector<double>>>, 3>& globalFlows,
                     const bool isIORank)
        : localFlows_(localFlows)
        , globalFlows_(globalFlows)
    {
        if (! isIORank) { return; }

        MessageBufferType buffer;
        this->pack(0, buffer);

        // pass a dummy_link to satisfy virtual class
        const int dummyLink = -1;
        this->unpack(dummyLink, buffer);
    }

    void pack(int link, MessageBufferType& buffer)
    {
        if (link != 0)
            throw std::logic_error("link in method pack is not 0 as expected");
        for (int i = 0; i < 3; ++i) {
            const auto& name = localFlows_[i].first;
            buffer.write(name);
            unsigned int size = localFlows_[i].second.first.size();
            buffer.write(size);
            for (unsigned int j = 0; j < size; ++j) {
                const auto& nncIdx = localFlows_[i].second.first[j];
                buffer.write(nncIdx);
                const auto& flows = localFlows_[i].second.second[j];
                buffer.write(flows);
            }
        }
    }

    void unpack(int /*link*/, MessageBufferType& buffer)
    {
        for (int i = 0; i < 3; ++i) {
            std::string name;
            buffer.read(name);
            globalFlows_[i].first = name;
            unsigned int size = 0;
            buffer.read(size);
            for (unsigned int j = 0; j < size; ++j) {
                int nncIdx;
                double data;
                buffer.read(nncIdx);
                buffer.read(data);
                if (nncIdx < 0)
                    continue;
                // This array is initialized in the collect(...) method below
                globalFlows_[i].second.second[nncIdx] = data;
            }
        }
    }
};

template <class Grid, class EquilGrid, class GridView>
CollectDataToIORank<Grid,EquilGrid,GridView>::
CollectDataToIORank(const Grid& grid, const EquilGrid* equilGrid,
                    const GridView& localGridView,
                    const Dune::CartesianIndexMapper<Grid>& cartMapper,
                    const Dune::CartesianIndexMapper<EquilGrid>* equilCartMapper,
                    const std::set<std::string>& fipRegionsInterregFlow)
    : toIORankComm_(grid.comm())
    , globalInterRegFlows_(EclInterRegFlowMap::createMapFromNames(toVector(fipRegionsInterregFlow)))
{
    // index maps only have to be build when reordering is needed
    if (!needsReordering && !isParallel())
        return;

    const CollectiveCommunication& comm = grid.comm();

    {
        std::set<int> send, recv;
        using EquilGridView = typename EquilGrid::LeafGridView;
        typename std::is_same<Grid, EquilGrid>::type isSameGrid;

        typedef Dune::MultipleCodimMultipleGeomTypeMapper<GridView> ElementMapper;
        ElementMapper elemMapper(localGridView, Dune::mcmgElementLayout());
        sortedCartesianIdx_.reserve(localGridView.size(0));

        for (const auto& elem : elements(localGridView))
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
            const std::size_t globalSize = equilGrid->leafGridView().size(0);
            globalCartesianIndex_.resize(globalSize, -1);
            const EquilGridView equilGridView = equilGrid->leafGridView();

            using EquilElementMapper = Dune::MultipleCodimMultipleGeomTypeMapper<EquilGridView>;
            EquilElementMapper equilElemMapper(equilGridView, Dune::mcmgElementLayout());

           // Scatter the global index to local index for lookup during restart
           if constexpr (isSameGrid) {
             ElementIndexScatterHandle<EquilElementMapper,ElementMapper> handle(equilElemMapper, elemMapper, localIdxToGlobalIdx_);
             grid.scatterData(handle);
           }

            // loop over all elements (global grid) and store Cartesian index
            for (const auto& elem : elements(equilGrid->leafGridView())) {
                int elemIdx = equilElemMapper.index(elem);
                int cartElemIdx = equilCartMapper->cartesianIndex(elemIdx);
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
            if constexpr (isSameGrid) {
              ElementIndexScatterHandle<ElementMapper, ElementMapper> handle(elemMapper, elemMapper, localIdxToGlobalIdx_);
              grid.scatterData(handle);
            }
        }

        // Sync the global element indices
        if constexpr (isSameGrid) {
          ElementIndexHandle<ElementMapper> handle(elemMapper, localIdxToGlobalIdx_);
          grid.communicate(handle, Dune::InteriorBorder_All_Interface,
                           Dune::ForwardCommunication);
        }

        localIndexMap_.clear();
        const std::size_t gridSize = grid.size(0);
        localIndexMap_.reserve(gridSize);

        // store the local Cartesian index
        IndexMapType distributedCartesianIndex;
        distributedCartesianIndex.resize(gridSize, -1);

        // A mapping for the whole grid (including the ghosts) is needed for restarts
        for (const auto& elem : elements(localGridView)) {
            int elemIdx = elemMapper.index(elem);
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
collect(const data::Solution&                                localCellData,
        const std::map<std::pair<std::string, int>, double>& localBlockData,
        const data::Wells&                                   localWellData,
        const data::WellBlockAveragePressures&               localWBPData,
        const data::GroupAndNetworkValues&                   localGroupAndNetworkData,
        const data::Aquifers&                                localAquiferData,
        const WellTestState&                                 localWellTestState,
        const EclInterRegFlowMap&                            localInterRegFlows,
        const std::array<std::pair<std::string, std::pair<std::vector<int>, std::vector<double>>>, 3>& localFlowsn,
        const std::array<std::pair<std::string, std::pair<std::vector<int>, std::vector<double>>>, 3>& localFloresn)
{
    globalCellData_ = {};
    globalBlockData_.clear();
    globalWellData_.clear();
    globalWBPData_.values.clear();
    globalGroupAndNetworkData_.clear();
    globalAquiferData_.clear();
    globalWellTestState_.clear();
    this->globalInterRegFlows_.clear();
    globalFlowsn_ = {};
    globalFloresn_ = {};

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

    // Set the right sizes for Flowsn and Floresn
    for (int i = 0; i < 3; ++i) {
        unsigned int sizeFlr = localFloresn[i].second.first.size();
        globalFloresn_[i].second.first.resize(sizeFlr, 0);
        globalFloresn_[i].second.second.resize(sizeFlr, 0.0);
        unsigned int sizeFlo = localFlowsn[i].second.first.size();
        globalFlowsn_[i].second.first.resize(sizeFlo, 0);
        globalFlowsn_[i].second.second.resize(sizeFlo, 0.0);
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

    PackUnPackWellTestState packUnpackWellTestState {
        localWellTestState,
        this->globalWellTestState_,
        this->isIORank()
    };

    PackUnpackInterRegFlows packUnpackInterRegFlows {
        localInterRegFlows,
        this->globalInterRegFlows_,
        this->isIORank()
    };

    PackUnpackFlows packUnpackFlowsn {
        localFlowsn,
        this->globalFlowsn_,
        this->isIORank()
    };

    PackUnpackFlows packUnpackFloresn {
        localFloresn,
        this->globalFloresn_,
        this->isIORank()
    };

    toIORankComm_.exchange(packUnpackCellData);
    toIORankComm_.exchange(packUnpackWellData);
    toIORankComm_.exchange(packUnpackGroupAndNetworkData);
    toIORankComm_.exchange(packUnpackBlockData);
    toIORankComm_.exchange(packUnpackWBPData);
    toIORankComm_.exchange(packUnpackAquiferData);
    toIORankComm_.exchange(packUnpackWellTestState);
    toIORankComm_.exchange(packUnpackInterRegFlows);
    toIORankComm_.exchange(packUnpackFlowsn);
    toIORankComm_.exchange(packUnpackFloresn);

#ifndef NDEBUG
    // make sure every process is on the same page
    toIORankComm_.barrier();
#endif
}

template <class Grid, class EquilGrid, class GridView>
int CollectDataToIORank<Grid,EquilGrid,GridView>::
localIdxToGlobalIdx(unsigned localIdx) const
{
    if (!isParallel()) {
        return localIdx;
    }

    if (this->localIdxToGlobalIdx_.empty()) {
        throw std::logic_error("index map is not created on this rank");
    }

    if (localIdx > this->localIdxToGlobalIdx_.size()) {
        throw std::logic_error("local index is outside map range");
    }

    return this->localIdxToGlobalIdx_[localIdx];
}

template <class Grid, class EquilGrid, class GridView>
bool CollectDataToIORank<Grid,EquilGrid,GridView>::
isCartIdxOnThisRank(int cartIdx) const
{
    if (! this->isParallel()) {
        return true;
    }

    assert (! needsReordering);

    auto candidate = std::lower_bound(this->sortedCartesianIdx_.begin(),
                                      this->sortedCartesianIdx_.end(),
                                      cartIdx);

    return (candidate != sortedCartesianIdx_.end())
        && (*candidate == cartIdx);
}

} // end namespace Opm
#endif
