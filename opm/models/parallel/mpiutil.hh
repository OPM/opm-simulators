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
 * \copydoc Opm::MpiBuffer
 */
#ifndef OPM_MATERIAL_MPIUTIL_HH
#define OPM_MATERIAL_MPIUTIL_HH

#include <dune/common/parallel/mpitraits.hh>

#include <cassert>
#include <numeric>
#include <string>
#include <vector>


#if HAVE_MPI

#include <mpi.h>



namespace mpiutil_details
{

    template <typename T>
    int packSize()
    {
        int pack_size;
        MPI_Pack_size(1, Dune::MPITraits<T>::getType(), MPI_COMM_WORLD, &pack_size);
        return pack_size;
    }

    // --------  Packer --------
    template <typename T>
    struct Packer
    {
        static int size(const T&)
        {
            return packSize<T>();
        }

        static void pack(const T& content, std::vector<char>& buf, int& offset)
        {
            MPI_Pack(&content, 1, Dune::MPITraits<T>::getType(), buf.data(), buf.size(), &offset, MPI_COMM_WORLD);
        }

        static T unpack(const std::vector<char>& recv_buffer, int& offset)
        {
            T content;
            auto* data = const_cast<char*>(recv_buffer.data());
            MPI_Unpack(data, recv_buffer.size(), &offset, &content, 1, Dune::MPITraits<T>::getType(), MPI_COMM_WORLD);
            return content;
        }
    };

    // --------  Packer, string specialization --------
    template <>
    struct Packer<std::string>
    {
        static int size(const std::string& content)
        {
            return packSize<unsigned int>() + content.size()*packSize<char>();
        }

        static void pack(const std::string& content, std::vector<char>& buf, int& offset)
        {
            unsigned int size = content.size();
            Packer<unsigned int>::pack(size, buf, offset);
            if (size > 0) {
                MPI_Pack(const_cast<char*>(content.c_str()), size, MPI_CHAR, buf.data(), buf.size(), &offset, MPI_COMM_WORLD);
            }
        }

        static std::string unpack(const std::vector<char>& recv_buffer, int& offset)
        {
            unsigned int size = Packer<unsigned int>::unpack(recv_buffer, offset);
            std::string text;
            if (size > 0) {
                auto* data = const_cast<char*>(recv_buffer.data());
                std::vector<char> chars(size);
                MPI_Unpack(data, recv_buffer.size(), &offset, chars.data(), size, MPI_CHAR, MPI_COMM_WORLD);
                text = std::string(chars.data(), size);
            }
            return text;
        }
    };

    // --------  Packer, vector partial specialization --------
    template <typename T>
    struct Packer<std::vector<T>>
    {
        static int size(const std::string& content)
        {
            int sz = 0;
            sz += packSize<unsigned int>();
            for (const T& elem : content) {
                sz += Packer<T>::size(elem);
            }
            return sz;
        }

        static void pack(const std::vector<T>& content, std::vector<char>& buf, int& offset)
        {
            unsigned int size = content.size();
            Packer<unsigned int>::pack(size, buf, offset);
            for (const T& elem : content) {
                Packer<T>::pack(elem);
            }
        }

        static std::vector<T> unpack(const std::vector<char>& recv_buffer, int& offset)
        {
            unsigned int size = Packer<T>::unpack(recv_buffer, offset);
            std::vector<T> content;
            content.reserve(size);
            for (unsigned int i = 0; i < size; ++i) {
                content.push_back(Packer<T>::unpack(recv_buffer, offset));
            }
            return content;
        }
    };


} // anonymous namespace


namespace Opm
{

    /// From each rank, gather its string (if not empty) into a vector.
    inline std::vector<std::string> gatherStrings(const std::string& local_string)
    {
        using StringPacker = mpiutil_details::Packer<std::string>;

        // Pack local messages.
        const int message_size = StringPacker::size(local_string);
        std::vector<char> buffer(message_size);
        int offset = 0;
        StringPacker::pack(local_string, buffer, offset);
        assert(offset == message_size);

        // Get message sizes and create offset/displacement array for gathering.
        int num_processes = -1;
        MPI_Comm_size(MPI_COMM_WORLD, &num_processes);
        std::vector<int> message_sizes(num_processes);
        MPI_Allgather(&message_size, 1, MPI_INT, message_sizes.data(), 1, MPI_INT, MPI_COMM_WORLD);
        std::vector<int> displ(num_processes + 1, 0);
        std::partial_sum(message_sizes.begin(), message_sizes.end(), displ.begin() + 1);

        // Gather.
        std::vector<char> recv_buffer(displ.back());
        MPI_Allgatherv(buffer.data(), buffer.size(), MPI_PACKED,
                       const_cast<char*>(recv_buffer.data()), message_sizes.data(),
                       displ.data(), MPI_PACKED,
                       MPI_COMM_WORLD);

        // Unpack and return.
        std::vector<std::string> ret;
        for (int process = 0; process < num_processes; ++process) {
            offset = displ[process];
            std::string s = StringPacker::unpack(recv_buffer, offset);
            if (!s.empty()) {
                ret.push_back(s);
            }
            assert(offset == displ[process + 1]);
        }
        return ret;
    }

} // namespace Opm

#else // HAVE_MPI

namespace Opm
{
    inline std::vector<std::string> gatherStrings(const std::string& local_string)
    {
        if (local_string.empty()) {
            return {};
        } else {
            return { local_string };
        }
    }
} // namespace Opm

#endif // HAVE_MPI

#endif // OPM_MATERIAL_MPIUTIL_HH

