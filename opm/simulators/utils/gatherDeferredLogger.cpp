/*
  Copyright 2019 SINTEF Digital, Mathematics and Cybernetics.
  Copyright 2019 Equinor.

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

#include "config.h"

#include <opm/simulators/utils/gatherDeferredLogger.hpp>

#if HAVE_MPI

#include <cassert>
#include <cstdint>
#include <numeric>
#include <mpi.h>

namespace
{

    void packMessages(const std::vector<Opm::DeferredLogger::Message>& local_messages, std::vector<char>& buf, int& offset, const Opm::Parallel::Communication mpi_communicator)
    {

        int messagesize = local_messages.size();
        MPI_Pack(&messagesize, 1, MPI_UNSIGNED, buf.data(), buf.size(), &offset, mpi_communicator);

        for (const auto& lm : local_messages) {
            MPI_Pack(static_cast<void*>(const_cast<std::int64_t*>(&lm.flag)), 1, MPI_INT64_T, buf.data(), buf.size(), &offset, mpi_communicator);
            int tagsize = lm.tag.size();
            MPI_Pack(&tagsize, 1, MPI_UNSIGNED, buf.data(), buf.size(), &offset, mpi_communicator);
            if (tagsize>0) {
                MPI_Pack(const_cast<char*>(lm.tag.c_str()), lm.tag.size(), MPI_CHAR, buf.data(), buf.size(), &offset, mpi_communicator);
            }
            int textsize = lm.text.size();
            MPI_Pack(&textsize, 1, MPI_UNSIGNED, buf.data(), buf.size(), &offset, mpi_communicator);
            if (textsize>0) {
                MPI_Pack(const_cast<char*>(lm.text.c_str()), lm.text.size(), MPI_CHAR, buf.data(), buf.size(), &offset, mpi_communicator);
            }
        }
    }

    Opm::DeferredLogger::Message unpackSingleMessage(const std::vector<char>& recv_buffer, int& offset, const MPI_Comm mpi_communicator)
    {
        int64_t flag;
        auto* data = const_cast<char*>(recv_buffer.data());
        MPI_Unpack(data, recv_buffer.size(), &offset, &flag, 1, MPI_INT64_T, mpi_communicator);

        // unpack tag
        unsigned int tagsize;
        MPI_Unpack(data, recv_buffer.size(), &offset, &tagsize, 1, MPI_UNSIGNED, mpi_communicator);
        std::string tag;
        if (tagsize>0) {
            std::vector<char> tagchars(tagsize);
            MPI_Unpack(data, recv_buffer.size(), &offset, tagchars.data(), tagsize, MPI_CHAR, mpi_communicator);
            tag = std::string(tagchars.data(), tagsize);
        }
        // unpack text
        unsigned int textsize;
        MPI_Unpack(data, recv_buffer.size(), &offset, &textsize, 1, MPI_UNSIGNED, mpi_communicator);
        std::string text;
        if (textsize>0) {
            std::vector<char> textchars(textsize);
            MPI_Unpack(data, recv_buffer.size(), &offset, textchars.data(), textsize, MPI_CHAR, mpi_communicator);
            text = std::string (textchars.data(), textsize);
        }
        return Opm::DeferredLogger::Message({flag, tag, text});
    }

    std::vector<Opm::DeferredLogger::Message> unpackMessages(const std::vector<char>& recv_buffer, const std::vector<int>& displ, const MPI_Comm mpi_communicator)
    {
        std::vector<Opm::DeferredLogger::Message> messages;
        const int num_processes = displ.size() - 1;
        auto* data = const_cast<char*>(recv_buffer.data());
        for (int process = 0; process < num_processes; ++process) {
            int offset = displ[process];
            // unpack number of messages
            unsigned int messagesize;
            MPI_Unpack(data, recv_buffer.size(), &offset, &messagesize, 1, MPI_UNSIGNED, mpi_communicator);
            for (unsigned int i=0; i<messagesize; i++) {
                messages.push_back(unpackSingleMessage(recv_buffer, offset, mpi_communicator));
            }
            assert(offset == displ[process + 1]);
        }
        return messages;
    }

} // anonymous namespace


namespace Opm
{

    /// combine (per-process) messages
    Opm::DeferredLogger gatherDeferredLogger(const Opm::DeferredLogger& local_deferredlogger,
                                             Opm::Parallel::Communication mpi_communicator)
    {

        int num_messages = local_deferredlogger.messages_.size();

        int int64_mpi_pack_size;
        MPI_Pack_size(1, MPI_INT64_T, mpi_communicator, &int64_mpi_pack_size);
        int unsigned_int_mpi_pack_size;
        MPI_Pack_size(1, MPI_UNSIGNED, mpi_communicator, &unsigned_int_mpi_pack_size);

        // store number of messages;
        int message_size = unsigned_int_mpi_pack_size;
        // store 1 int64 per message for flag
        message_size += num_messages*int64_mpi_pack_size;
        // store 2 unsigned ints per message for length of tag and length of text
        message_size += num_messages*2*unsigned_int_mpi_pack_size;

        for (const auto& lm : local_deferredlogger.messages_) {
            int string_mpi_pack_size;
            MPI_Pack_size(lm.tag.size(), MPI_CHAR, mpi_communicator, &string_mpi_pack_size);
            message_size += string_mpi_pack_size;
            MPI_Pack_size(lm.text.size(), MPI_CHAR, mpi_communicator, &string_mpi_pack_size);
            message_size += string_mpi_pack_size;
        }

        // Pack local messages.
        std::vector<char> buffer(message_size);

        int offset = 0;
        packMessages(local_deferredlogger.messages_, buffer, offset, mpi_communicator);
        assert(offset == message_size);

        // Get message sizes and create offset/displacement array for gathering.
        int num_processes = -1;
        MPI_Comm_size(mpi_communicator, &num_processes);
        std::vector<int> message_sizes(num_processes);
        MPI_Allgather(&message_size, 1, MPI_INT, message_sizes.data(), 1, MPI_INT, mpi_communicator);
        std::vector<int> displ(num_processes + 1, 0);
        std::partial_sum(message_sizes.begin(), message_sizes.end(), displ.begin() + 1);

        // Gather.
        std::vector<char> recv_buffer(displ.back());
        MPI_Allgatherv(buffer.data(), buffer.size(), MPI_PACKED,
                       const_cast<char*>(recv_buffer.data()), message_sizes.data(),
                       displ.data(), MPI_PACKED,
                       mpi_communicator);

        // Unpack.
        Opm::DeferredLogger global_deferredlogger;
        global_deferredlogger.messages_ = unpackMessages(recv_buffer, displ, mpi_communicator);
        return global_deferredlogger;
    }

} // namespace Opm

#else // HAVE_MPI

namespace Opm
{
    Opm::DeferredLogger gatherDeferredLogger(const Opm::DeferredLogger& local_deferredlogger,
                                             Opm::Parallel::Communication /* dummy communicator */)
    {
        return local_deferredlogger;
    }
} // namespace Opm

#endif // HAVE_MPI
