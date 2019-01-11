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

#include <opm/simulators/gatherDeferredLogger.hpp>

#if HAVE_MPI

#include <cassert>
#include <numeric>
#include <mpi.h>

namespace
{

    void packMessages(const std::vector<Opm::DeferredLogger::Message>& local_messages, std::vector<char>& buf, int& offset)
    {

        for (const auto lm : local_messages) {
            MPI_Pack(&lm.flag, 1, MPI_INT64_T, buf.data(), buf.size(), &offset, MPI_COMM_WORLD);
            int tagsize = lm.tag.size();
            MPI_Pack(&tagsize, 1, MPI_UNSIGNED, buf.data(), buf.size(), &offset, MPI_COMM_WORLD);
            if (tagsize>0) {
                MPI_Pack(lm.tag.c_str(), lm.tag.size(), MPI_CHAR, buf.data(), buf.size(), &offset, MPI_COMM_WORLD);
            }
            int textsize = lm.text.size();
            MPI_Pack(&textsize, 1, MPI_UNSIGNED, buf.data(), buf.size(), &offset, MPI_COMM_WORLD);
            if (textsize>0) {
                MPI_Pack(lm.text.c_str(), lm.text.size(), MPI_CHAR, buf.data(), buf.size(), &offset, MPI_COMM_WORLD);
            }
        }
    }

    Opm::DeferredLogger::Message unpackSingleMessage(const std::vector<char>& recv_buffer, int& offset)
    {
        int64_t flag;
        auto* data = const_cast<char*>(recv_buffer.data());
        MPI_Unpack(data, recv_buffer.size(), &offset, &flag, 1, MPI_INT64_T, MPI_COMM_WORLD);

        // unpack tag
        unsigned int tagsize;
        MPI_Unpack(data, recv_buffer.size(), &offset, &tagsize, 1, MPI_UNSIGNED, MPI_COMM_WORLD);
        std::string tag;
        if (tagsize>0) {
            std::vector<char> tagchars(tagsize);
            MPI_Unpack(data, recv_buffer.size(), &offset, tagchars.data(), tagsize, MPI_CHAR, MPI_COMM_WORLD);
            tag = std::string(tagchars.data(), tagsize);
        }
        // unpack text
        unsigned int textsize;
        MPI_Unpack(data, recv_buffer.size(), &offset, &textsize, 1, MPI_UNSIGNED, MPI_COMM_WORLD);
        std::string text;
        if (textsize>0) {
            std::vector<char> textchars(textsize);
            MPI_Unpack(data, recv_buffer.size(), &offset, textchars.data(), textsize, MPI_CHAR, MPI_COMM_WORLD);
            text = std::string (textchars.data(), textsize);
        }
        return Opm::DeferredLogger::Message({flag, tag, text});
    }

    std::vector<Opm::DeferredLogger::Message> unpackMessages(const std::vector<char>& recv_buffer, const std::vector<int>& displ)
    {
        std::vector<Opm::DeferredLogger::Message> messages;
        const int num_processes = displ.size() - 1;
        for (int process = 0; process < num_processes; ++process) {
            int offset = displ[process];
            messages.push_back(unpackSingleMessage(recv_buffer, offset));
            assert(offset == displ[process + 1]);
        }
        return messages;
    }

} // anonymous namespace


namespace Opm
{

    /// combine (per-process) messages
    Opm::DeferredLogger gatherDeferredLogger(const Opm::DeferredLogger& local_deferredlogger)
    {
        // Pack local messages.
        int message_size = 0;
        for (const auto lm : local_deferredlogger.messages_) {
            message_size += sizeof(lm.flag);
            message_size += sizeof(unsigned int);// to store the length of tag
            message_size += lm.tag.size();
            message_size += sizeof(unsigned int);// to store the length of text
            message_size += lm.text.size();
        }

        //////int message_size = local_messages.size()
        std::vector<char> buffer(message_size);
        int offset = 0;
        packMessages(local_deferredlogger.messages_, buffer, offset);
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

        // Unpack.
        std::vector<Opm::DeferredLogger::Message> m = unpackMessages(recv_buffer, displ);
        Opm::DeferredLogger global_deferredlogger;
        global_deferredlogger.messages_ = m;
        return global_deferredlogger;
    }

} // namespace Opm

#else // HAVE_MPI

namespace Opm
{
    Opm::DeferredLogger gatherDeferredLogger(const Opm::DeferredLogger& local_deferredlogger)
    {
        return local_deferredlogger;
    }
} // namespace Opm

#endif // HAVE_MPI
