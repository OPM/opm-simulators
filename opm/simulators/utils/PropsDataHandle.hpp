/*
  Copyright 2020, 2023 Equinor AS.
  Copyright 2026 SINTEF Digital

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
/*!
 * \file PropsDatahandle.hpp
 * \brief File containing a data handle for communicating the FieldProperties
 *
 * \author Markus Blatt, OPM-OP AS
 */

#ifndef PROPS_DATAHANDLE_HPP
#define PROPS_DATAHANDLE_HPP

#if HAVE_MPI

#include <opm/input/eclipse/EclipseState/Grid/FieldData.hpp>

#include <opm/simulators/utils/MPISerializer.hpp>
#include <opm/simulators/utils/ParallelEclipseState.hpp>
#include <opm/simulators/utils/ParallelRestart.hpp>
#include <opm/grid/LookUpData.hh>

#include <dune/grid/common/datahandleif.hh>
#include <dune/grid/common/mcmgmapper.hh>
#include <dune/grid/common/partitionset.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <algorithm>
#include <numeric>
#include <stdexcept>
#include <unordered_map>
#include <iostream>

namespace Opm
{

/*!
 * \brief A Data handle to communicate the field properties during load balance.
 * \tparam Grid The type of grid where the load balancing is happening.
 * \todo Maybe specialize this for CpGrid to save some space, later.
 */
template<class Grid>
class PropsDataHandle
    : public Dune::CommDataHandleIF< PropsDataHandle<Grid>, double>
{
public:
    //! \brief the data type we send (ints are converted to double)
    using DataType = std::pair<double, unsigned char>;

    //! \brief Constructor
    //! \param grid The grid where the loadbalancing is happening.
    //! \param globalProps The field properties of the global grid
    //! \param distributedProps The distributed field properties
    PropsDataHandle(const Grid& grid, ParallelEclipseState& eclState)
        : m_grid(grid),
          m_distributed_fieldProps(eclState.m_fieldProps)
    {
        // Scatter the keys
        const Parallel::Communication comm = m_grid.comm();
        if (comm.rank() == 0)
        {
            const FieldPropsManager& globalProps = eclState.globalFieldProps();
            m_intKeys = globalProps.keys<int>();
            m_doubleKeys = globalProps.keys<double>();
            // PORV is derived on demand and else read through the parent's index on every rank.
            if (m_grid.maxLevel() > 0 &&
                std::find(m_doubleKeys.begin(), m_doubleKeys.end(), "PORV") == m_doubleKeys.end()) {
                static_cast<void>(globalProps.get_double("PORV"));
                m_doubleKeys.push_back("PORV");
            }
            m_distributed_fieldProps.copyTran(globalProps);

            // Multi-valued fields such as ZMF use component-major storage:
            // the cell index varies fastest.
            m_doubleMult.reserve(m_doubleKeys.size());
            for (const auto& doubleKey : m_doubleKeys)
            {
                const auto& fieldData = globalProps.get_double_field_data(doubleKey,
                                                                          /* allow_unsupported = */ true);
                m_doubleMult.push_back(fieldData.numValuePerCell());
            }
        }

        Parallel::MpiSerializer ser(comm);
        ser.broadcast(Parallel::RootRank{0}, *this);

        int hasMultiValuedIntField = 0;
        if (comm.rank() == 0) {
            const auto& globalProps = eclState.globalFieldProps();
            hasMultiValuedIntField = std::any_of(
                m_intKeys.begin(), m_intKeys.end(),
                [&globalProps](const std::string& key) {
                    return globalProps.get_int_field_data(key).numValuePerCell() > 1;
                });
        }
        comm.broadcast(&hasMultiValuedIntField, 1, 0);
        if (hasMultiValuedIntField) {
            throw std::runtime_error {
                "Distributing multi-valued integer field properties is not supported"
            };
        }

        int hasLgr = (comm.rank() == 0) && (eclState.getLgrs().size() > 0);
        comm.broadcast(&hasLgr, 1, 0);
        const bool hasMultiValuedField =
            std::any_of(m_doubleMult.begin(), m_doubleMult.end(),
                        [](const std::size_t multiplicity) { return multiplicity > 1; });
        if (hasLgr && hasMultiValuedField) {
            throw std::runtime_error {
                "Distributing multi-valued field properties with LGRs is not supported"
            };
        }

        m_no_data = m_intKeys.size() +
                    std::accumulate(m_doubleMult.begin(), m_doubleMult.end(), std::size_t{0});

        if (comm.rank() == 0) {
            const FieldPropsManager& globalProps = eclState.globalFieldProps();
            const auto& idSet = m_grid.localIdSet();
            const std::size_t numCells = m_grid.levelGridView(0).size(0);

            // Record a single source cell's properties into elementData_, keyed
            // by the cell's local id; propIndex indexes the (level-zero/Cartesian
            // sized) global field properties.
            // A grid refined before load balancing ships its leaf flat: give refined cells the
            // values LookUpData would (PORV by volume share, an LGR's own arrays), since the
            // distributed cells have no father to derive them from.
            const bool refinedLeaf = m_grid.maxLevel() > 0;
            using LeafView = typename Grid::LeafGridView;
            const LeafView leafView = m_grid.leafGridView();  // LookUpData keeps a reference
            const LookUpData<Grid, LeafView> lookup(leafView);
            const auto leafPorv = (refinedLeaf && globalProps.has_double("PORV"))
                ? lookup.assignFieldPropsDoubleOnLeaf(globalProps, "PORV") : std::vector<double>{};

            auto record = [&](const auto& element, const std::size_t propIndex)
            {
                const bool refined = refinedLeaf && element.hasFather();
                const auto& id = idSet.id(element);
                auto& data = elementData_[id];
                data.reserve(m_no_data);

                for (const auto& intKey : m_intKeys)
                {
                    const auto& fieldData = globalProps.get_int_field_data(intKey);
                    const double value = refined && globalProps.has_int(intKey)
                        ? lookup.fieldPropInt(globalProps, intKey, element) : fieldData.data[propIndex];
                    data.emplace_back(value,
                                      static_cast<unsigned char>(fieldData.value_status[propIndex]));
                }

                for (std::size_t keyIdx = 0; keyIdx < m_doubleKeys.size(); ++keyIdx)
                {
                    // We need to allow unsupported keywords to get the data
                    // for TranCalculator, too.
                    const auto& fieldData = globalProps.get_double_field_data(m_doubleKeys[keyIdx],
                                                                              /* allow_unsupported = */ true);
                    for (std::size_t comp = 0; comp < m_doubleMult[keyIdx]; ++comp)
                    {
                        const auto dataIdx = comp * numCells + propIndex;
                        double value = fieldData.data[dataIdx];
                        if (refined && m_doubleKeys[keyIdx] == "PORV" && !leafPorv.empty()) {
                            value = leafPorv[leafView.indexSet().index(element)];
                        }
                        else if (refined && globalProps.has_double(m_doubleKeys[keyIdx])) {
                            value = lookup.fieldPropDouble(globalProps, m_doubleKeys[keyIdx], element);
                        }
                        data.emplace_back(value,
                                          static_cast<unsigned char>(fieldData.value_status[dataIdx]));
                    }
                }
            };

            // The view being scattered is the leaf: level zero, or the refined leaf
            // in refine-before-redistribute. The serial grid is all on rank 0 here,
            // and each cell takes the properties of its level-zero origin.
            for (const auto& element : elements(m_grid.leafGridView(), Dune::Partitions::all))
            {
                record(element, element.getOrigin().index());
            }
        }
    }

    ~PropsDataHandle()
    {
        // distributed grid is now correctly set up.
        const auto& gridView = m_grid.levelGridView(0);
        const std::size_t numCells = gridView.size(0);
        for (const auto& intKey : m_intKeys)
        {
            m_distributed_fieldProps.m_intProps[intKey].data.resize(numCells);
            m_distributed_fieldProps.m_intProps[intKey].value_status.resize(numCells);
        }

        for (std::size_t keyIdx = 0; keyIdx < m_doubleKeys.size(); ++keyIdx)
        {
            auto& props = m_distributed_fieldProps.m_doubleProps[m_doubleKeys[keyIdx]];
            props.data.resize(m_doubleMult[keyIdx] * numCells);
            props.value_status.resize(m_doubleMult[keyIdx] * numCells);
            // Preserve the number of values per cell for numCells() and compress().
            props.kw_info.num_value_per_cell(m_doubleMult[keyIdx]);
        }

        // copy data for the persistent mao to the field properties
        const auto& idSet = m_grid.localIdSet();
        using ElementMapper =
            Dune::MultipleCodimMultipleGeomTypeMapper<typename Grid::LevelGridView>;
        ElementMapper elemMapper(gridView, Dune::mcmgElementLayout());

        for (const auto &element : elements( gridView, Dune::Partitions::all))
        {
            std::size_t counter{};
            const auto& id = idSet.id(element);
            auto index = elemMapper.index(element);
            auto data = elementData_.find(id);
            assert(data != elementData_.end());

            for (const auto& intKey : m_intKeys)
            {
                const auto& pair = data->second[counter++];
                m_distributed_fieldProps.m_intProps[intKey].data[index] = static_cast<int>(pair.first);
                m_distributed_fieldProps.m_intProps[intKey].value_status[index] = static_cast<value::status>(pair.second);
            }

            for (std::size_t keyIdx = 0; keyIdx < m_doubleKeys.size(); ++keyIdx)
            {
                auto& props = m_distributed_fieldProps.m_doubleProps[m_doubleKeys[keyIdx]];
                for (std::size_t comp = 0; comp < m_doubleMult[keyIdx]; ++comp)
                {
                    const auto& pair = data->second[counter++];
                    props.data[comp * numCells + index] = pair.first;
                    props.value_status[comp * numCells + index] = static_cast<value::status>(pair.second);
                }
            }
        }
    }

    bool contains(int /* dim */, int codim)
    {
        return codim == 0;
    }

    bool fixedsize(int /* dim */, int /* codim */)
    {
        return true;
    }
    bool fixedSize(int /* dim */, int /* codim */)
    {
        return true;
    }

    template<class EntityType>
    std::size_t size(const EntityType /* entity */)
    {
        return m_no_data;
    }

    template<class BufferType, class EntityType>
    void gather(BufferType& buffer, const EntityType& e) const
    {
        auto iter = elementData_.find(m_grid.localIdSet().id(e));
        assert(iter != elementData_.end());
        for (const auto& data : iter->second)
        {
            buffer.write(data);
        }
    }

    template<class BufferType, class EntityType>
    void scatter(BufferType& buffer, const EntityType& e, std::size_t n)
    {
        assert(n == m_no_data);
        auto& array = elementData_[m_grid.localIdSet().id(e)];
        array.resize(n);
        for (auto& data : array)
        {
            buffer.read(data);
        }
    }

    template<class Serializer>
    void serializeOp(Serializer& serializer)
    {
        serializer(m_intKeys);
        serializer(m_doubleKeys);
        serializer(m_doubleMult);
        m_distributed_fieldProps.serializeOp(serializer);
    }

private:
    using LocalIdSet = typename Grid::LocalIdSet;
    const Grid& m_grid;
    //! \brief The distributed field properties for receiving
    ParallelFieldPropsManager& m_distributed_fieldProps;
    //! \brief The names of the keys of the integer fields.
    std::vector<std::string> m_intKeys;
    //! \brief The names of the keys of the double fields.
    std::vector<std::string> m_doubleKeys;
    //! \brief The number of values per cell of each double field.
    std::vector<std::size_t> m_doubleMult;
    /// \brief The data per element as a vector mapped from the local id.
    ///
    /// each entry is a pair of data and value_status.
    std::unordered_map<typename LocalIdSet::IdType, std::vector<std::pair<double,unsigned char> > > elementData_;
    /// \brief The amount of data to send for each element
    std::size_t m_no_data;
};

} // end namespace Opm

#endif // HAVE_MPI
#endif // PROPS_DATAHANDLE_HPP
