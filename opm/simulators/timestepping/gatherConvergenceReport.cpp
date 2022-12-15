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

#if HAVE_MPI

#include <mpi.h>

namespace
{

    using Opm::ConvergenceReport;

    void packReservoirFailure(const ConvergenceReport::ReservoirFailure& f,
                              std::vector<char>& buf,
                              int& offset, MPI_Comm mpi_communicator)
    {
        int type = static_cast<int>(f.type());
        int severity = static_cast<int>(f.severity());
        int phase = f.phase();
        MPI_Pack(&type, 1, MPI_INT, buf.data(), buf.size(), &offset, mpi_communicator);
        MPI_Pack(&severity, 1, MPI_INT, buf.data(), buf.size(), &offset, mpi_communicator);
        MPI_Pack(&phase, 1, MPI_INT, buf.data(), buf.size(), &offset, mpi_communicator);
    }

    void packReservoirConvergenceMetric(const ConvergenceReport::ReservoirConvergenceMetric& m,
                                        std::vector<char>& buf,
                                        int& offset, MPI_Comm mpi_communicator)
    {
        int type = static_cast<int>(m.type());
        int phase = m.phase();
        double value = m.value();
        MPI_Pack(&type, 1, MPI_INT, buf.data(), buf.size(), &offset, mpi_communicator);
        MPI_Pack(&phase, 1, MPI_INT, buf.data(), buf.size(), &offset, mpi_communicator);
        MPI_Pack(&value, 1, MPI_DOUBLE, buf.data(), buf.size(), &offset, mpi_communicator);
    }

    void packWellFailure(const ConvergenceReport::WellFailure& f,
                         std::vector<char>& buf,
                         int& offset, MPI_Comm mpi_communicator)
    {
        int type = static_cast<int>(f.type());
        int severity = static_cast<int>(f.severity());
        int phase = f.phase();
        MPI_Pack(&type, 1, MPI_INT, buf.data(), buf.size(), &offset, mpi_communicator);
        MPI_Pack(&severity, 1, MPI_INT, buf.data(), buf.size(), &offset, mpi_communicator);
        MPI_Pack(&phase, 1, MPI_INT, buf.data(), buf.size(), &offset, mpi_communicator);
        int name_length = f.wellName().size() + 1; // Adding 1 for the null terminator.
        MPI_Pack(&name_length, 1, MPI_INT, buf.data(), buf.size(), &offset, mpi_communicator);
        MPI_Pack(const_cast<char*>(f.wellName().c_str()), name_length, MPI_CHAR, buf.data(), buf.size(), &offset, mpi_communicator);
    }

    void packConvergenceReport(const ConvergenceReport& local_report,
                               std::vector<char>& buf,
                               int& offset, MPI_Comm mpi_communicator)
    {
        // Pack the data.
        // Status will not be packed, it is possible to deduce from the other data.
        // Reservoir failures.
        double reportTime = local_report.reportTime();
        MPI_Pack(&reportTime, 1, MPI_DOUBLE, buf.data(), buf.size(), &offset, mpi_communicator);
        const auto rf = local_report.reservoirFailures();
        int num_rf = rf.size();
        MPI_Pack(&num_rf, 1, MPI_INT, buf.data(), buf.size(), &offset, mpi_communicator);
        for (const auto& f : rf) {
            packReservoirFailure(f, buf, offset, mpi_communicator);
        }
        // Reservoir convergence metrics.
        const auto rm = local_report.reservoirConvergence();
        int num_rm = rm.size();
        MPI_Pack(&num_rm, 1, MPI_INT, buf.data(), buf.size(), &offset, mpi_communicator);
        for (const auto& m : rm) {
            packReservoirConvergenceMetric(m, buf, offset, mpi_communicator);
        }
        // Well failures.
        const auto wf = local_report.wellFailures();
        int num_wf = wf.size();
        MPI_Pack(&num_wf, 1, MPI_INT, buf.data(), buf.size(), &offset, mpi_communicator);
        for (const auto& f : wf) {
            packWellFailure(f, buf, offset, mpi_communicator);
        }
    }

    int messageSize(const ConvergenceReport& local_report, MPI_Comm mpi_communicator)
    {
        int int_pack_size = 0;
        MPI_Pack_size(1, MPI_INT, mpi_communicator, &int_pack_size);
        int double_pack_size = 0;
        MPI_Pack_size(1, MPI_DOUBLE, mpi_communicator, &double_pack_size);
        const int num_rf = local_report.reservoirFailures().size();
        const int num_rm = local_report.reservoirConvergence().size();
        const int num_wf = local_report.wellFailures().size();
        int wellnames_length = 0;
        for (const auto& f : local_report.wellFailures()) {
            wellnames_length += (f.wellName().size() + 1);
        }
        return (3 + 3*num_rf + 2*num_rm + 4*num_wf)*int_pack_size + (1 + 1*num_rm)*double_pack_size + wellnames_length;
    }

    ConvergenceReport::ReservoirFailure unpackReservoirFailure(const std::vector<char>& recv_buffer, int& offset, MPI_Comm mpi_communicator)
    {
        int type = -1;
        int severity = -1;
        int phase = -1;
        auto* data = const_cast<char*>(recv_buffer.data());
        MPI_Unpack(data, recv_buffer.size(), &offset, &type, 1, MPI_INT, mpi_communicator);
        MPI_Unpack(data, recv_buffer.size(), &offset, &severity, 1, MPI_INT, mpi_communicator);
        MPI_Unpack(data, recv_buffer.size(), &offset, &phase, 1, MPI_INT, mpi_communicator);
        return ConvergenceReport::ReservoirFailure(static_cast<ConvergenceReport::ReservoirFailure::Type>(type),
                                                   static_cast<ConvergenceReport::Severity>(severity),
                                                   phase);
    }

    ConvergenceReport::ReservoirConvergenceMetric
    unpackReservoirConvergenceMetric(const std::vector<char>& recv_buffer, int& offset, MPI_Comm mpi_communicator)
    {
        int type = -1;
        int phase = -1;
        double value = -1.0;
        auto* data = const_cast<char*>(recv_buffer.data());
        MPI_Unpack(data, recv_buffer.size(), &offset, &type, 1, MPI_INT, mpi_communicator);
        MPI_Unpack(data, recv_buffer.size(), &offset, &phase, 1, MPI_INT, mpi_communicator);
        MPI_Unpack(data, recv_buffer.size(), &offset, &value, 1, MPI_DOUBLE, mpi_communicator);
        return { static_cast<ConvergenceReport::ReservoirFailure::Type>(type), phase, value };
    }

    ConvergenceReport::WellFailure unpackWellFailure(const std::vector<char>& recv_buffer, int& offset, MPI_Comm mpi_communicator)
    {
        int type = -1;
        int severity = -1;
        int phase = -1;
        auto* data = const_cast<char*>(recv_buffer.data());
        MPI_Unpack(data, recv_buffer.size(), &offset, &type, 1, MPI_INT, mpi_communicator);
        MPI_Unpack(data, recv_buffer.size(), &offset, &severity, 1, MPI_INT, mpi_communicator);
        MPI_Unpack(data, recv_buffer.size(), &offset, &phase, 1, MPI_INT, mpi_communicator);
        int name_length = -1;
        MPI_Unpack(data, recv_buffer.size(), &offset, &name_length, 1, MPI_INT, mpi_communicator);
        std::vector<char> namechars(name_length);
        MPI_Unpack(data, recv_buffer.size(), &offset, namechars.data(), name_length, MPI_CHAR, mpi_communicator);
        std::string name(namechars.data());
        return ConvergenceReport::WellFailure(static_cast<ConvergenceReport::WellFailure::Type>(type),
                                              static_cast<ConvergenceReport::Severity>(severity),
                                              phase,
                                              name);
    }

    ConvergenceReport unpackSingleConvergenceReport(const std::vector<char>& recv_buffer, int& offset, MPI_Comm mpi_communicator)
    {
        auto* data = const_cast<char*>(recv_buffer.data());
        double reportTime{0.0};
        MPI_Unpack(data, recv_buffer.size(), &offset, &reportTime, 1, MPI_DOUBLE, mpi_communicator);
        ConvergenceReport cr{reportTime};
        int num_rf = -1;
        MPI_Unpack(data, recv_buffer.size(), &offset, &num_rf, 1, MPI_INT, mpi_communicator);
        for (int rf = 0; rf < num_rf; ++rf) {
            ConvergenceReport::ReservoirFailure f = unpackReservoirFailure(recv_buffer, offset, mpi_communicator);
            cr.setReservoirFailed(f);
        }
        int num_rm = -1;
        MPI_Unpack(data, recv_buffer.size(), &offset, &num_rm, 1, MPI_INT, mpi_communicator);
        for (int rm = 0; rm < num_rm; ++rm) {
            cr.setReservoirConvergenceMetric(unpackReservoirConvergenceMetric(recv_buffer, offset, mpi_communicator));
        }
        int num_wf = -1;
        MPI_Unpack(data, recv_buffer.size(), &offset, &num_wf, 1, MPI_INT, mpi_communicator);
        for (int wf = 0; wf < num_wf; ++wf) {
            ConvergenceReport::WellFailure f = unpackWellFailure(recv_buffer, offset, mpi_communicator);
            cr.setWellFailed(f);
        }
        return cr;
    }

    ConvergenceReport unpackConvergenceReports(const std::vector<char>& recv_buffer,
                                               const std::vector<int>& displ, MPI_Comm mpi_communicator)
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
    ConvergenceReport gatherConvergenceReport(const ConvergenceReport& local_report, Parallel::Communication mpi_communicator)
    {
        // Pack local report.
        int message_size = messageSize(local_report, mpi_communicator);
        std::vector<char> buffer(message_size);
        int offset = 0;
        packConvergenceReport(local_report, buffer, offset,mpi_communicator);
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
        ConvergenceReport global_report = unpackConvergenceReports(recv_buffer, displ, mpi_communicator);
        return global_report;
    }

} // namespace Opm

#else // HAVE_MPI

namespace Opm
{
    ConvergenceReport gatherConvergenceReport(const ConvergenceReport& local_report,
                                              Parallel::Communication mpi_communicator [[maybe_unused]])
    {
        return local_report;
    }
} // namespace Opm

#endif // HAVE_MPI
