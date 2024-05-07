/*
  Copyright 2018, 2022 Equinor ASA.
  Copyright 2018 SINTEF Digital, Mathematics and Cybernetics.

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

#include <opm/simulators/timestepping/gatherConvergenceReport.hpp>

#include <algorithm>
#include <tuple>
#include <utility>
#include <vector>

#if HAVE_MPI

#include <mpi.h>

namespace
{

    using Opm::ConvergenceReport;

    void packReservoirFailure(const ConvergenceReport::ReservoirFailure& f,
                              std::vector<char>& buf, int& offset,
                              MPI_Comm mpi_communicator)
    {
        auto pack = [&buf, &offset, mpi_communicator]
            (const auto* ptr, const int size, const auto type)
        {
            MPI_Pack(ptr, size, type,
                     buf.data(), buf.size(), &offset,
                     mpi_communicator);
        };

        const int type = static_cast<int>(f.type());
        const int severity = static_cast<int>(f.severity());
        const int phase = f.phase();

        pack(&type, 1, MPI_INT);
        pack(&severity, 1, MPI_INT);
        pack(&phase, 1, MPI_INT);
    }

    void packReservoirConvergenceMetric(const ConvergenceReport::ReservoirConvergenceMetric& m,
                                        std::vector<char>& buf, int& offset,
                                        MPI_Comm mpi_communicator)
    {
        auto pack = [&buf, &offset, mpi_communicator]
            (const auto* ptr, const int size, const auto type)
        {
            MPI_Pack(ptr, size, type,
                     buf.data(), buf.size(), &offset,
                     mpi_communicator);
        };

        const int type = static_cast<int>(m.type());
        const int phase = m.phase();
        const double value = m.value();

        pack(&type, 1, MPI_INT);
        pack(&phase, 1, MPI_INT);
        pack(&value, 1, MPI_DOUBLE);
    }

    void packWellFailure(const ConvergenceReport::WellFailure& f,
                         std::vector<char>& buf, int& offset,
                         MPI_Comm mpi_communicator)
    {
        auto pack = [&buf, &offset, mpi_communicator]
            (const auto* ptr, const int size, const auto type)
        {
            MPI_Pack(ptr, size, type,
                     buf.data(), buf.size(), &offset,
                     mpi_communicator);
        };

        const int type = static_cast<int>(f.type());
        const int severity = static_cast<int>(f.severity());
        const int phase = f.phase();
        pack(&type, 1, MPI_INT);
        pack(&severity, 1, MPI_INT);
        pack(&phase, 1, MPI_INT);

        // Add one for the null terminator.
        const int name_length = f.wellName().size() + 1;
        pack(&name_length, 1, MPI_INT);
        pack(f.wellName().c_str(), name_length, MPI_CHAR);
    }

    void packConvergenceReport(const ConvergenceReport& local_report,
                               std::vector<char>& buf, int& offset,
                               MPI_Comm mpi_communicator)
    {
        auto pack = [&buf, &offset, mpi_communicator]
            (const auto* ptr, const int size, const auto type)
        {
            MPI_Pack(ptr, size, type,
                     buf.data(), buf.size(), &offset,
                     mpi_communicator);
        };

        // Pack the data.
        // Status will not be packed, it is possible to deduce from the other data.
        const auto reportTime = local_report.reportTime();
        pack(&reportTime, 1, MPI_DOUBLE);

        // CNV pore-volume split
        {
            const auto& cnvPvSplit = local_report.cnvPvSplit();
            const auto  eligiblePoreVolume = local_report.eligiblePoreVolume();
            const auto  num_cnv_pv = static_cast<int>(cnvPvSplit.first.size());

            pack(&eligiblePoreVolume     , 1         , MPI_DOUBLE);
            pack(&num_cnv_pv             , 1         , MPI_INT);
            pack(cnvPvSplit.first .data(), num_cnv_pv, MPI_DOUBLE);
            pack(cnvPvSplit.second.data(), num_cnv_pv, MPI_INT);
        }

        // Reservoir failures.
        {
            const auto rf = local_report.reservoirFailures();
            const int num_rf = rf.size();

            pack(&num_rf, 1, MPI_INT);

            for (const auto& f : rf) {
                packReservoirFailure(f, buf, offset, mpi_communicator);
            }
        }

        // Reservoir convergence metrics.
        {
            const auto rm = local_report.reservoirConvergence();
            const int num_rm = rm.size();

            pack(&num_rm, 1, MPI_INT);

            for (const auto& m : rm) {
                packReservoirConvergenceMetric(m, buf, offset, mpi_communicator);
            }
        }

        // Well failures.
        {
            const auto wf = local_report.wellFailures();
            const int num_wf = wf.size();

            pack(&num_wf, 1, MPI_INT);

            for (const auto& f : wf) {
                packWellFailure(f, buf, offset, mpi_communicator);
            }
        }
    }

    int messageSize(const ConvergenceReport& local_report, MPI_Comm mpi_communicator)
    {
        int int_pack_size = 0;
        MPI_Pack_size(1, MPI_INT, mpi_communicator, &int_pack_size);

        int double_pack_size = 0;
        MPI_Pack_size(1, MPI_DOUBLE, mpi_communicator, &double_pack_size);

        int char_pack_size = 0;
        MPI_Pack_size(1, MPI_CHAR, mpi_communicator, &char_pack_size);

        const int num_cnv_pv = local_report.cnvPvSplit().first.size();
        const int num_rf = local_report.reservoirFailures().size();
        const int num_rm = local_report.reservoirConvergence().size();
        const int num_wf = local_report.wellFailures().size();

        int wellnames_length = 0;
        for (const auto& f : local_report.wellFailures()) {
            // Add one for the null terminator.
            wellnames_length += f.wellName().size() + 1;
        }

        return (3 + 1 + num_cnv_pv + 3*num_rf + 2*num_rm + 4*num_wf)*int_pack_size
            +  (1 + 1 + num_cnv_pv + 1*num_rm)*double_pack_size
            +  wellnames_length*char_pack_size;
    }

    ConvergenceReport::ReservoirFailure
    unpackReservoirFailure(const std::vector<char>& recv_buffer, int& offset, MPI_Comm mpi_communicator)
    {
        auto unpackInt = [data = recv_buffer.data(),
                          size = static_cast<int>(recv_buffer.size()),
                          &offset, mpi_communicator]()
        {
            auto x = -1;
            MPI_Unpack(data, size, &offset, &x, 1, MPI_INT, mpi_communicator);

            return x;
        };

        const auto type     = unpackInt();
        const auto severity = unpackInt();
        const auto phase    = unpackInt();

        return {
            static_cast<ConvergenceReport::ReservoirFailure::Type>(type),
            static_cast<ConvergenceReport::Severity>(severity),
            phase
        };
    }

    ConvergenceReport::ReservoirConvergenceMetric
    unpackReservoirConvergenceMetric(const std::vector<char>& recv_buffer, int& offset, MPI_Comm mpi_communicator)
    {
        auto unpack = [data = recv_buffer.data(),
                       size = static_cast<int>(recv_buffer.size()),
                       &offset, mpi_communicator](const auto type, auto x)
        {
            MPI_Unpack(data, size, &offset, &x, 1, type, mpi_communicator);

            return x;
        };

        const auto type  = unpack(MPI_INT, -1);
        const auto phase = unpack(MPI_INT, -1);
        const auto value = unpack(MPI_DOUBLE, -1.0);

        return { static_cast<ConvergenceReport::ReservoirFailure::Type>(type), phase, value };
    }

    ConvergenceReport::WellFailure
    unpackWellFailure(const std::vector<char>& recv_buffer, int& offset, MPI_Comm mpi_communicator)
    {
        auto unpackInt = [data = recv_buffer.data(),
                          size = static_cast<int>(recv_buffer.size()),
                          &offset, mpi_communicator]()
        {
            auto x = -1;
            MPI_Unpack(data, size, &offset, &x, 1, MPI_INT, mpi_communicator);

            return x;
        };

        const auto type     = unpackInt();
        const auto severity = unpackInt();
        const auto phase    = unpackInt();

        const auto name_length = unpackInt();
        std::vector<char> namechars(name_length);
        MPI_Unpack(recv_buffer.data(), recv_buffer.size(), &offset,
                   namechars.data(), name_length,
                   MPI_CHAR, mpi_communicator);

        return {
            static_cast<ConvergenceReport::WellFailure::Type>(type),
            static_cast<ConvergenceReport::Severity>(severity),
            phase,
            const_cast<const std::vector<char>&>(namechars).data()
        };
    }

    ConvergenceReport
    unpackSingleConvergenceReport(const std::vector<char>& recv_buffer, int& offset, MPI_Comm mpi_communicator)
    {
        auto unpack = [data = recv_buffer.data(),
                       size = static_cast<int>(recv_buffer.size()),
                       &offset, mpi_communicator](const auto type, auto x)
        {
            MPI_Unpack(data, size, &offset, &x, 1, type, mpi_communicator);

            return x;
        };

        auto cr = ConvergenceReport { unpack(MPI_DOUBLE, 0.0) };

        {
            const auto eligiblePoreVolume = unpack(MPI_DOUBLE, 0.0);
            const auto num_cnv_pv = unpack(MPI_INT, -1);

            auto cnvPvSplit = ConvergenceReport::CnvPvSplit {
                std::piecewise_construct,
                std::forward_as_tuple(num_cnv_pv),
                std::forward_as_tuple(num_cnv_pv)
            };

            const auto* data = recv_buffer.data();

            MPI_Unpack(data, recv_buffer.size(), &offset,
                       cnvPvSplit.first.data(), num_cnv_pv,
                       MPI_DOUBLE, mpi_communicator);

            MPI_Unpack(data, recv_buffer.size(), &offset,
                       cnvPvSplit.second.data(), num_cnv_pv,
                       MPI_DOUBLE, mpi_communicator);

            cr.setCnvPoreVolSplit(cnvPvSplit, eligiblePoreVolume);
        }

        {
            const auto num_rf = unpack(MPI_INT, -1);

            for (int rf = 0; rf < num_rf; ++rf) {
                cr.setReservoirFailed(unpackReservoirFailure(recv_buffer, offset, mpi_communicator));
            }
        }

        {
            const auto num_rm = unpack(MPI_INT, -1);

            for (int rm = 0; rm < num_rm; ++rm) {
                cr.setReservoirConvergenceMetric(unpackReservoirConvergenceMetric(recv_buffer, offset, mpi_communicator));
            }
        }

        {
            const auto num_wf = unpack(MPI_INT, -1);

            for (int wf = 0; wf < num_wf; ++wf) {
                cr.setWellFailed(unpackWellFailure(recv_buffer, offset, mpi_communicator));
            }
        }

        return cr;
    }

    ConvergenceReport
    unpackConvergenceReports(const std::vector<char>& recv_buffer,
                             const std::vector<int>& displ,
                             MPI_Comm mpi_communicator)
    {
        ConvergenceReport cr;

        const int num_processes = displ.size() - 1;
        for (int process = 0; process < num_processes; ++process) {
            int offset = displ[process];
            cr += unpackSingleConvergenceReport(recv_buffer, offset, mpi_communicator);
            assert(offset == displ[process + 1]);
        }

        return cr;
    }

} // anonymous namespace


namespace Opm
{

    /// Create a global convergence report combining local
    /// (per-process) reports.
    ConvergenceReport gatherConvergenceReport(const ConvergenceReport& local_report,
                                              Parallel::Communication mpi_communicator)
    {
        // Pack local report.
        const int message_size = messageSize(local_report, mpi_communicator);
        std::vector<char> buffer(message_size);
        {
            int offset = 0;
            packConvergenceReport(local_report, buffer, offset, mpi_communicator);
            assert(offset == message_size);
        }

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
                       recv_buffer.data(), message_sizes.data(),
                       displ.data(), MPI_PACKED,
                       mpi_communicator);

        // Unpack.
        return unpackConvergenceReports(recv_buffer, displ, mpi_communicator);
    }

} // namespace Opm

#else // HAVE_MPI

namespace Opm
{
    ConvergenceReport gatherConvergenceReport(const ConvergenceReport& local_report,
                                              [[maybe_unused]] Parallel::Communication mpi_communicator)
    {
        return local_report;
    }
} // namespace Opm

#endif // HAVE_MPI
