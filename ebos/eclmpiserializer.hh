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
#ifndef ECL_MPI_SERIALIZER_HH
#define ECL_MPI_SERIALIZER_HH

#include <opm/simulators/utils/ParallelRestart.hpp>

namespace Opm {

class EclMpiSerializer {
public:
    enum class Operation {
        PACKSIZE,
        PACK,
        UNPACK
    };

    explicit EclMpiSerializer(Dune::CollectiveCommunication<Dune::MPIHelper::MPICommunicator> comm) :
        m_comm(comm)
    {}

    template<class T>
    void operator()(const T& data)
    {
        if constexpr (is_shared_ptr<T>::value) {
            shared_ptr(data);
        } else {
          if (m_op == Operation::PACKSIZE)
              m_packSize += Mpi::packSize(data, m_comm);
          else if (m_op == Operation::PACK)
              Mpi::pack(data, m_buffer, m_position, m_comm);
          else if (m_op == Operation::UNPACK)
              Mpi::unpack(const_cast<T&>(data), m_buffer, m_position, m_comm);
        }
    }

    template<class T>
    void vector(std::vector<T>& data)
    {
        auto handle = [&](auto& d)
        {
            for (auto& it : d) {
              if constexpr (is_pair<T>::value)
                  pair(it);
              else if constexpr (is_shared_ptr<T>::value)
                  shared_ptr(it);
              else
                  it.serializeOp(*this);
            }
        };

        if (m_op == Operation::PACKSIZE) {
            m_packSize += Mpi::packSize(data.size(), m_comm);
            handle(data);
        } else if (m_op == Operation::PACK) {
            Mpi::pack(data.size(), m_buffer, m_position, m_comm);
            handle(data);
        } else if (m_op == Operation::UNPACK) {
            size_t size;
            Mpi::unpack(size, m_buffer, m_position, m_comm);
            data.resize(size);
            handle(data);
        }
    }

    template<template<class Key, class Data> class Map, class Key, class Data>
    void map(Map<Key, Data>& data)
    {
        auto handle = [&](auto& d)
        {
            if constexpr (is_vector<Data>::value)
                vector(d);
            else if constexpr (is_shared_ptr<Data>::value)
                shared_ptr(d);
            else
                d.serializeOp(*this);
        };

        if (m_op == Operation::PACKSIZE) {
            m_packSize += Mpi::packSize(data.size(), m_comm);
            for (auto& it : data) {
                m_packSize += Mpi::packSize(it.first, m_comm);
                handle(it.second);
            }
        } else if (m_op == Operation::PACK) {
            Mpi::pack(data.size(), m_buffer, m_position, m_comm);
            for (auto& it : data) {
                Mpi::pack(it.first, m_buffer, m_position, m_comm);
                handle(it.second);
            }
        } else if (m_op == Operation::UNPACK) {
            size_t size;
            Mpi::unpack(size, m_buffer, m_position, m_comm);
            for (size_t i = 0; i < size; ++i) {
                Key key;
                Mpi::unpack(key, m_buffer, m_position, m_comm);
                Data entry;
                handle(entry);
                data.insert(std::make_pair(key, entry));
            }
        }
    }

    template<class T>
    void pack(T& data)
    {
        m_op = Operation::PACKSIZE;
        m_packSize = 0;
        data.serializeOp(*this);
        m_position = 0;
        m_buffer.resize(m_packSize);
        m_op = Operation::PACK;
        data.serializeOp(*this);
    }

    template<class T>
    void unpack(T& data)
    {
        m_position = 0;
        m_op = Operation::UNPACK;
        data.serializeOp(*this);
    }

    template<class T>
    void broadcast(T& data)
    {
        if (m_comm.size() == 1)
            return;

#if HAVE_MPI
        if (m_comm.rank() == 0) {
            pack(data);
            m_comm.broadcast(&m_position, 1, 0);
            m_comm.broadcast(m_buffer.data(), m_position, 0);
        } else {
            m_comm.broadcast(&m_packSize, 1, 0);
            m_buffer.resize(m_packSize);
            m_comm.broadcast(m_buffer.data(), m_packSize, 0);
            unpack(data);
        }
#endif
    }

    size_t position() const
    {
        return m_position;
    }

protected:
    template<class T>
    struct is_pair {
        constexpr static bool value = false;
    };

    template<class T1, class T2>
    struct is_pair<std::pair<T1,T2>> {
        constexpr static bool value = true;
    };

    template<class T>
    struct is_vector {
        constexpr static bool value = false;
    };

    template<class T1>
    struct is_vector<std::vector<T1>> {
        constexpr static bool value = true;
    };

    template<class T>
    struct is_shared_ptr {
        constexpr static bool value = false;
    };

    template<class T1>
    struct is_shared_ptr<std::shared_ptr<T1>> {
        constexpr static bool value = true;
    };

    template<class T1, class T2>
    void pair(const std::pair<T1,T2>& data)
    {
        if constexpr (std::is_pod<T1>::value || std::is_same<T1,std::string>::value)
            (*this)(data.first);
        else
            data.first.serializeOp(*this);

        if constexpr (std::is_pod<T2>::value || std::is_same<T2,std::string>::value)
            (*this)(data.second);
        else
            const_cast<T2&>(data.second).serializeOp(*this);
    }

    template<class T1>
    void shared_ptr(const std::shared_ptr<T1>& data)
    {
        bool value = data ? true : false;
        (*this)(value);
        if (m_op == Operation::UNPACK && value) {
            const_cast<std::shared_ptr<T1>&>(data) = std::make_shared<T1>();
        }
        if (data)
            data->serializeOp(*this);
    }

    Dune::CollectiveCommunication<Dune::MPIHelper::MPICommunicator> m_comm;

    Operation m_op = Operation::PACKSIZE;
    size_t m_packSize = 0;
    int m_position = 0;
    std::vector<char> m_buffer;
};

}

#endif
