/*
  Copyright 2016 Dr. Blatt - HPC-Simulation-Software & Services
  Copyright 2016 Statoil AS

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

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif // HAVE_CONFIG_H

#include "WellSwitchingLogger.hpp"
#include <numeric>

namespace Opm
{
namespace wellhelpers
{

#if HAVE_MPI
int WellSwitchingLogger::calculateMessageSize(std::vector<int>& well_name_lengths)
{

    // Each process will send a message to the root process with
    // the following data:
    // total number of switches, for each switch the length of the
    // well name, for each switch the well name and the two controls.
    well_name_lengths.reserve(switchMap_.size());

    for(const auto& switchEntry : switchMap_)
    {
        int length = switchEntry.first.size() +1;  //we write an additional \0
        well_name_lengths.push_back(length);
    }

    // compute the message size
    int message_size = 0;
    int increment = 0;
    // number of switches
    MPI_Pack_size(1, MPI_INT, MPI_COMM_WORLD, &message_size);
    // const char* length include delimiter for each switch
    MPI_Pack_size(switchMap_.size(), MPI_INT, MPI_COMM_WORLD, &increment);
    message_size += increment;

    // for each well the name + two controls in one write
    for(const auto& length : well_name_lengths)
    {
        // well name
        MPI_Pack_size(length, MPI_CHAR, MPI_COMM_WORLD, &increment);
        message_size += increment;
        // controls
        MPI_Pack_size(2, MPI_CHAR, MPI_COMM_WORLD, &increment);
        message_size += increment;
    }
    return message_size;
}

void WellSwitchingLogger::packData(std::vector<int>& well_name_lengths,
                                   std::vector<char>& buffer)
{
    // Pack the data
    // number of switches
    int offset = 0;
    int no_switches = switchMap_.size();
    MPI_Pack(&no_switches, 1, MPI_INT, buffer.data(), buffer.size(),
             &offset, MPI_COMM_WORLD);
    MPI_Pack(well_name_lengths.data(), well_name_lengths.size(),
             MPI_INT, buffer.data(), buffer.size(),
             &offset, MPI_COMM_WORLD);

    for(const auto& switchEntry : switchMap_)
    {
        // well name
        auto& well_name = switchEntry.first;
        MPI_Pack(const_cast<char*>(well_name.c_str()), well_name.size()+1,
                 MPI_CHAR, buffer.data(), buffer.size(),
                 &offset, MPI_COMM_WORLD);

        // controls
        MPI_Pack(const_cast<char*>(switchEntry.second.data()), 2 , MPI_CHAR,
                 buffer.data(), buffer.size(), &offset, MPI_COMM_WORLD);
    }
}

void WellSwitchingLogger::unpackDataAndLog(std::vector<char>& recv_buffer,
                                           const std::vector<int>& displ)
{
    for(int p=1; p < cc_.size(); ++p)
    {
        int offset = displ[p];
        int no_switches = 0;
        MPI_Unpack(recv_buffer.data(), recv_buffer.size(), &offset,
                   &no_switches, 1, MPI_INT, MPI_COMM_WORLD);

        if ( no_switches == 0 )
        {
            continue;
        }

        std::vector<int> well_name_lengths(no_switches);

        MPI_Unpack(recv_buffer.data(), recv_buffer.size(), &offset,
                   well_name_lengths.data(), well_name_lengths.size(),
                   MPI_INT, MPI_COMM_WORLD);

        std::vector<char> well_name;
        for ( int i = 0; i < no_switches; ++i )
        {
            well_name.resize(well_name_lengths[i]);
            MPI_Unpack(recv_buffer.data(), recv_buffer.size(), &offset,
                       well_name.data(), well_name_lengths[i], MPI_CHAR,
                       MPI_COMM_WORLD);

            std::array<char,2> fromto{{}};
            MPI_Unpack(recv_buffer.data(), recv_buffer.size(), &offset,
                       fromto.data(), 2, MPI_CHAR, MPI_COMM_WORLD);

            logSwitch(well_name.data(), fromto, p);
        }
    }
}

void WellSwitchingLogger::logSwitch(const char* name, std::array<char,2> fromto,
                                    int rank)
{
            std::ostringstream ss;
            ss << "    Switching control mode for well " << name
               << " from " << modestring[WellControlType(fromto[0])]
               << " to " <<  modestring[WellControlType(fromto[1])]
               << " on rank " << rank;
            OpmLog::info(ss.str());
}
#endif

void WellSwitchingLogger::gatherDataAndLog()
{

#if HAVE_MPI
    if(cc_.size() == 1)
    {
        return;
    }

    std::vector<int> message_sizes;
    std::vector<int> well_name_lengths;
    int message_size = calculateMessageSize(well_name_lengths);

    if ( cc_.rank() == 0 ){
        for(const auto& entry : switchMap_)
        {
            logSwitch(entry.first.c_str(), entry.second,0);
        }

        message_sizes.resize(cc_.size());
    }

    MPI_Gather(&message_size, 1, MPI_INT, message_sizes.data(),
               1, MPI_INT, 0, MPI_COMM_WORLD);

    std::vector<char> buffer(message_size);
    packData(well_name_lengths, buffer);

    std::vector<int> displ;

    if ( cc_.rank() == 0){
        // last entry will be total size of
        displ.resize(cc_.size() + 1, 0);
        std::partial_sum(message_sizes.begin(), message_sizes.end(),
                         displ.begin()+1);
    }
    std::vector<char> recv_buffer;
    if ( cc_.rank() == 0 ){
        recv_buffer.resize(displ[cc_.size()]);
    }

    MPI_Gatherv(buffer.data(), buffer.size(), MPI_PACKED,
                recv_buffer.data(), message_sizes.data(),
                displ.data(), MPI_PACKED, 0, MPI_COMM_WORLD);
    if ( cc_.rank() == 0 )
    {
        unpackDataAndLog(recv_buffer, displ);
    }
#endif
}


WellSwitchingLogger::~WellSwitchingLogger()
{
    gatherDataAndLog();
}
} // end namespace wellhelpers
} // end namespace Opm
