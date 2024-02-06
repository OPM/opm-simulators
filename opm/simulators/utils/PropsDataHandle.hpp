/*
  Copyright 2020, 2023 Equinor AS.

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
#include <dune/grid/common/datahandleif.hh>
#include <dune/grid/common/mcmgmapper.hh>
#include <dune/grid/common/partitionset.hh>
#include <dune/common/parallel/mpihelper.hh>
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
            m_distributed_fieldProps.copyTran(globalProps);
        }

        Parallel::MpiSerializer ser(comm);
        ser.broadcast(*this);

        m_no_data = m_intKeys.size() + m_doubleKeys.size();

        if (comm.rank() == 0) {
            const FieldPropsManager& globalProps = eclState.globalFieldProps();
            const auto& idSet = m_grid.localIdSet();
            const auto& gridView = m_grid.levelGridView(0);
            using ElementMapper =
                Dune::MultipleCodimMultipleGeomTypeMapper<typename Grid::LevelGridView>;
            ElementMapper elemMapper(gridView, Dune::mcmgElementLayout());

            for (const auto &element : elements(gridView, Dune::Partitions::interiorBorder))
            {
                const auto& id = idSet.id(element);
                auto index = elemMapper.index(element);
                auto& data = elementData_[id];
                data.reserve(m_no_data);

                for (const auto& intKey : m_intKeys)
                {
                    const auto& fieldData = globalProps.get_int_field_data(intKey);
                    data.emplace_back(fieldData.data[index],
                                      static_cast<unsigned char>(fieldData.value_status[index]));
                }

                for (const auto& doubleKey : m_doubleKeys)
                {
                    // We need to allow unsupported keywords to get the data
                    // for TranCalculator, too.
                    const auto& fieldData = globalProps.get_double_field_data(doubleKey,
                                                                              /* allow_unsupported = */ true);
                    data.emplace_back(fieldData.data[index],
                                      static_cast<unsigned char>(fieldData.value_status[index]));
                }
            }
        }
    }

    ~PropsDataHandle()
    {
        // distributed grid is now correctly set up.
        for (const auto& intKey : m_intKeys)
        {
            m_distributed_fieldProps.m_intProps[intKey].data.resize(m_grid.size(0));
            m_distributed_fieldProps.m_intProps[intKey].value_status.resize(m_grid.size(0));
        }

        for (const auto& doubleKey : m_doubleKeys)
        {
            m_distributed_fieldProps.m_doubleProps[doubleKey].data.resize(m_grid.size(0));
            m_distributed_fieldProps.m_doubleProps[doubleKey].value_status.resize(m_grid.size(0));
        }

        // copy data for the persistent mao to the field properties
        const auto& idSet = m_grid.localIdSet();
        const auto& gridView = m_grid.levelGridView(0);
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

            for (const auto& doubleKey : m_doubleKeys)
            {
                const auto& pair = data->second[counter++];
                m_distributed_fieldProps.m_doubleProps[doubleKey].data[index] = pair.first;
                m_distributed_fieldProps.m_doubleProps[doubleKey].value_status[index] = static_cast<value::status>(pair.second);
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
