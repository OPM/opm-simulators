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
#ifndef MPI_SERIALIZER_HPP
#define MPI_SERIALIZER_HPP

#include <opm/common/utility/Serializer.hpp>
#include <opm/simulators/utils/MPIPacker.hpp>
#include <opm/simulators/utils/ParallelCommunication.hpp>

namespace Opm::Parallel {

//! \brief Class for serializing and broadcasting data using MPI.
class MpiSerializer : public Serializer<Mpi::Packer> {
public:
    MpiSerializer(Parallel::Communication comm)
        : Serializer<Mpi::Packer>(m_packer)
        , m_packer(comm)
        , m_comm(comm)
    {}

    //! \brief Serialize and broadcast on root process, de-serialize on
    //! others.
    //!
    //! \tparam T Type of class to broadcast
    //! \param data Class to broadcast
    //! \param root Process to broadcast from
    template<class T>
    void broadcast(T& data, int root = 0)
    {
        if (m_comm.size() == 1)
            return;

        if (m_comm.rank() == root) {
            try {
                this->pack(data);
                m_comm.broadcast(&m_packSize, 1, root);
                m_comm.broadcast(m_buffer.data(), m_packSize, root);
            } catch (...) {
                m_packSize = std::numeric_limits<size_t>::max();
                m_comm.broadcast(&m_packSize, 1, root);
                throw;
            }
        } else {
            m_comm.broadcast(&m_packSize, 1, root);
            if (m_packSize == std::numeric_limits<size_t>::max()) {
                throw std::runtime_error("Error detected in parallel serialization");
            }

            m_buffer.resize(m_packSize);
            m_comm.broadcast(m_buffer.data(), m_packSize, root);
            this->unpack(data);
        }
    }

    template<typename... Args>
    void broadcast(int root, Args&&... args)
    {
        if (m_comm.size() == 1)
            return;

        if (m_comm.rank() == root) {
            try {
                this->pack(std::forward<Args>(args)...);
                m_comm.broadcast(&m_packSize, 1, root);
                m_comm.broadcast(m_buffer.data(), m_packSize, root);
            } catch (...) {
                m_packSize = std::numeric_limits<size_t>::max();
                m_comm.broadcast(&m_packSize, 1, root);
                throw;
            }
        } else {
            m_comm.broadcast(&m_packSize, 1, root);
            if (m_packSize == std::numeric_limits<size_t>::max()) {
                throw std::runtime_error("Error detected in parallel serialization");
            }
            m_buffer.resize(m_packSize);
            m_comm.broadcast(m_buffer.data(), m_packSize, root);
            this->unpack(std::forward<Args>(args)...);
        }
    }

    //! \brief Serialize and broadcast on root process, de-serialize and append on
    //! others.
    //!
    //! \tparam T Type of class to broadcast
    //! \param data Class to broadcast
    //! \param root Process to broadcast from
    template<class T>
    void append(T& data, int root = 0)
    {
        if (m_comm.size() == 1)
            return;

        T tmp;
        T& bcast = m_comm.rank() == root ? data : tmp;
        broadcast(bcast, root);

        if (m_comm.rank() != root)
            data.append(tmp);
    }

private:
    const Mpi::Packer m_packer; //!< Packer instance
    Parallel::Communication m_comm; //!< Communicator to use
};

}

#endif
