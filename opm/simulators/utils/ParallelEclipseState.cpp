/*
  Copyright 2019 Equinor AS.

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

#include "ParallelEclipseState.hpp"
#include "ParallelRestart.hpp"
#include <ebos/eclmpiserializer.hh>

namespace Opm {


ParallelFieldPropsManager::ParallelFieldPropsManager(FieldPropsManager& manager)
    : m_manager(manager)
    , m_comm(Dune::MPIHelper::getCollectiveCommunication())
{
}


std::vector<int> ParallelFieldPropsManager::actnum() const
{
    if (m_comm.rank() == 0)
        return m_manager.actnum();

    return{};
}


void ParallelFieldPropsManager::reset_actnum(const std::vector<int>& actnum)
{
    if (m_comm.rank() != 0)
        OPM_THROW(std::runtime_error, "reset_actnum should only be called on root process.");
    m_manager.reset_actnum(actnum);
}


std::vector<double> ParallelFieldPropsManager::porv(bool global) const
{
    std::vector<double> result;
    if (m_comm.rank() == 0)
        result = m_manager.porv(global);
    size_t size = result.size();
    m_comm.broadcast(&size, 1, 0);
    result.resize(size);
    m_comm.broadcast(result.data(), size, 0);
    return result;
}


const std::vector<int>& ParallelFieldPropsManager::get_int(const std::string& keyword) const
{
    auto it = m_intProps.find(keyword);
    if (it == m_intProps.end())
        OPM_THROW(std::runtime_error, "No integer property field: " + keyword);

    return it->second;
}

std::vector<int> ParallelFieldPropsManager::get_global_int(const std::string& keyword) const
{
    std::vector<int> result;
    if (m_comm.rank() == 0)
        result = m_manager.get_global_int(keyword);
    size_t size = result.size();
    m_comm.broadcast(&size, 1, 0);
    result.resize(size);
    m_comm.broadcast(result.data(), size, 0);

    return result;
}


const std::vector<double>& ParallelFieldPropsManager::get_double(const std::string& keyword) const
{
    auto it = m_doubleProps.find(keyword);
    if (it == m_doubleProps.end())
        OPM_THROW(std::runtime_error, "No double property field: " + keyword);

    return it->second;
}


std::vector<double> ParallelFieldPropsManager::get_global_double(const std::string& keyword) const
{
    std::vector<double> result;
    if (m_comm.rank() == 0)
        result = m_manager.get_global_double(keyword);
    size_t size = result.size();
    m_comm.broadcast(&size, 1, 0);
    result.resize(size);
    m_comm.broadcast(result.data(), size, 0);

    return result;
}


bool ParallelFieldPropsManager::has_int(const std::string& keyword) const
{
    auto it = m_intProps.find(keyword);
    return it != m_intProps.end();
}


bool ParallelFieldPropsManager::has_double(const std::string& keyword) const
{
    auto it = m_doubleProps.find(keyword);
    return it != m_doubleProps.end();
}


ParallelEclipseState::ParallelEclipseState()
    : m_fieldProps(field_props)
{
}


ParallelEclipseState::ParallelEclipseState(const Deck& deck)
    : EclipseState(deck)
    , m_fieldProps(field_props)
{
}


std::size_t ParallelEclipseState::packSize(EclMpiSerializer& serializer) const
{
    return serializer.packSize(m_tables) +
           serializer.packSize(m_runspec) +
           serializer.packSize(m_eclipseConfig) +
           serializer.packSize(m_deckUnitSystem) +
           serializer.packSize(m_inputNnc) +
           serializer.packSize(m_inputEditNnc) +
           serializer.packSize(m_gridDims) +
           serializer.packSize(m_simulationConfig) +
           serializer.packSize(m_transMult) +
           serializer.packSize(m_faults) +
           serializer.packSize(m_title);

}


void ParallelEclipseState::pack(std::vector<char>& buffer, int& position,
                                EclMpiSerializer& serializer) const
{
    serializer.pack(m_tables, buffer, position);
    serializer.pack(m_runspec, buffer, position);
    serializer.pack(m_eclipseConfig, buffer, position);
    serializer.pack(m_deckUnitSystem, buffer, position);
    serializer.pack(m_inputNnc, buffer, position);
    serializer.pack(m_inputEditNnc, buffer, position);
    serializer.pack(m_gridDims, buffer, position);
    serializer.pack(m_simulationConfig, buffer, position);
    serializer.pack(m_transMult, buffer, position);
    serializer.pack(m_faults, buffer, position);
    serializer.pack(m_title, buffer, position);
}


void ParallelEclipseState::unpack(std::vector<char>& buffer, int& position,
                                  EclMpiSerializer& serializer)
{
    serializer.unpack(m_tables, buffer, position);
    serializer.unpack(m_runspec, buffer, position);
    serializer.unpack(m_eclipseConfig, buffer, position);
    serializer.unpack(m_deckUnitSystem, buffer, position);
    serializer.unpack(m_inputNnc, buffer, position);
    serializer.unpack(m_inputEditNnc, buffer, position);
    serializer.unpack(m_gridDims, buffer, position);
    serializer.unpack(m_simulationConfig, buffer, position);
    serializer.unpack(m_transMult, buffer, position);
    serializer.unpack(m_faults, buffer, position);
    serializer.unpack(m_title, buffer, position);
}


const FieldPropsManager& ParallelEclipseState::fieldProps() const
{
    if (!m_parProps && Dune::MPIHelper::getCollectiveCommunication().rank() != 0)
        OPM_THROW(std::runtime_error, "Attempt to access field properties on no-root process before switch to parallel properties");

    if (!m_parProps || Dune::MPIHelper::getCollectiveCommunication().size() == 1)
        return this->EclipseState::fieldProps();

    return m_fieldProps;
}


const FieldPropsManager& ParallelEclipseState::globalFieldProps() const
{
    if (Dune::MPIHelper::getCollectiveCommunication().rank() != 0)
        OPM_THROW(std::runtime_error, "Attempt to access global field properties on non-root process");
    return this->EclipseState::globalFieldProps();
}


const EclipseGrid& ParallelEclipseState::getInputGrid() const
{
    if (Dune::MPIHelper::getCollectiveCommunication().rank() != 0)
        OPM_THROW(std::runtime_error, "Attempt to access eclipse grid on non-root process");
    return this->EclipseState::getInputGrid();
}


void ParallelEclipseState::switchToGlobalProps()
{
    m_parProps = false;
}


void ParallelEclipseState::switchToDistributedProps()
{
    const auto& comm = Dune::MPIHelper::getCollectiveCommunication();
    if (comm.size() == 1) // No need for the parallel frontend
        return;

    m_parProps = true;
}


namespace {


template<class T>
struct GetField {
  GetField(const FieldPropsManager& propMan) : props(propMan) {}
  std::vector<T> getField(const std::string& key) const;
  const FieldPropsManager& props;
};


template<>
std::vector<int> GetField<int>::getField(const std::string& key) const {
  return props.get_global_int(key);
}


template<>
std::vector<double> GetField<double>::getField(const std::string& key) const {
  return props.get_global_double(key);
}


template<class T>
void extractRootProps(const std::vector<int>& localToGlobal,
                    const std::vector<std::string>& keys,
                    const GetField<T>& getter,
                    std::map<std::string,std::vector<T>>& localMap)
{
    for (const std::string& key : keys) {
        auto prop = getter.getField(key);
        std::vector<T>& local = localMap[key];
        local.reserve(localToGlobal.size());
        for (int cell : localToGlobal) {
            local.push_back(prop[cell]);
        }
    }
}


template<class T>
void packProps(const std::vector<int>& l2gCell,
             const std::vector<std::string>& keys,
             const GetField<T>& getter,
             std::vector<char>& buffer, int& position)
{
    const auto& comm = Dune::MPIHelper::getCollectiveCommunication();
    std::vector<T> sendData(l2gCell.size());
    for (const std::string& key : keys) {
        auto prop = getter.getField(key);
        size_t idx = 0;
        for (int cell : l2gCell)
            sendData[idx++] = prop[cell];
        Mpi::pack(sendData, buffer, position, comm);
    }
}


}


void ParallelEclipseState::setupLocalProps(const std::vector<int>& localToGlobal)
{
#if HAVE_MPI
    const auto& comm = Dune::MPIHelper::getCollectiveCommunication();
    if (comm.rank() == 0) {
        extractRootProps(localToGlobal, this->globalFieldProps().keys<int>(),
                         GetField<int>(this->globalFieldProps()),
                         m_fieldProps.m_intProps);
        extractRootProps(localToGlobal, this->globalFieldProps().keys<double>(),
                         GetField<double>(this->globalFieldProps()),
                         m_fieldProps.m_doubleProps);
        for (int i = 1; i < comm.size(); ++i) {
            std::vector<int> l2gCell;
            size_t size;
            MPI_Recv(&size, 1, Dune::MPITraits<size_t>::getType(), i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            l2gCell.resize(size);
            MPI_Recv(l2gCell.data(), size, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            size_t cells = l2gCell.size();
            const auto& intKeys = this->globalFieldProps().keys<int>();
            const auto& dblKeys = this->globalFieldProps().keys<double>();
            size = Mpi::packSize(intKeys, comm) +
                   Mpi::packSize(dblKeys,comm) +
                   intKeys.size() * Mpi::packSize(std::vector<int>(cells), comm) +
                   dblKeys.size() * Mpi::packSize(std::vector<double>(cells), comm);

            std::vector<char> buffer(size);
            int position = 0;
            Mpi::pack(intKeys, buffer, position, comm);
            Mpi::pack(dblKeys, buffer, position, comm);
            packProps(l2gCell, intKeys, GetField<int>(this->globalFieldProps()),
                      buffer, position);
            packProps(l2gCell, dblKeys, GetField<double>(this->globalFieldProps()),
                      buffer, position);
            MPI_Send(&position, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
            MPI_Send(buffer.data(), position, MPI_CHAR, i, 0, MPI_COMM_WORLD);
        }
    } else {
        size_t l2gSize = localToGlobal.size();
        MPI_Send(&l2gSize, 1, Dune::MPITraits<size_t>::getType(), 0, 0, MPI_COMM_WORLD);
        MPI_Send(localToGlobal.data(), localToGlobal.size(), MPI_INT, 0, 0, MPI_COMM_WORLD);
        int size;
        MPI_Recv(&size, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        std::vector<char> buffer(size);
        MPI_Recv(buffer.data(), size, MPI_CHAR, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        std::vector<std::string> intKeys, dblKeys;
        int position = 0;
        Mpi::unpack(intKeys, buffer, position, comm);
        Mpi::unpack(dblKeys, buffer, position, comm);
        for (const std::string& key : intKeys) {
            Mpi::unpack(m_fieldProps.m_intProps[key], buffer, position, comm);
        }
        for (const std::string& key : dblKeys) {
            Mpi::unpack(m_fieldProps.m_doubleProps[key], buffer, position, comm);
        }
    }
#endif
}


} // end namespace Opm
