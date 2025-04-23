/*
  Copyright 2018, 2022, 2024 Equinor ASA.
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

#include <opm/simulators/utils/ParallelCommunication.hpp>

#if HAVE_MPI

#include <opm/common/utility/Serializer.hpp>

#include <opm/simulators/utils/MPIPacker.hpp>

#include <opm/grid/common/CommunicationUtils.hpp>

#include <algorithm>
#include <tuple>
#include <vector>

namespace {
    /// Special purpose utility to collect each rank's local convergence
    /// report object and distribute those to all ranks.
    ///
    /// This particular feature could arguably be a part of the
    /// Opm::MpiSerializer class, but we create it here first as a narrowly
    /// scoped testing ground.
    ///
    /// Inherits from Opm::Serializer<> in order to use that type's pack()
    /// and unpack() member functions along with its internal 'm_buffer'
    /// data member.
    class CollectConvReports : private Opm::Serializer<Opm::Mpi::Packer>
    {
    public:
        /// Constructor.
        ///
        /// \param[in] packer Packing object to convert in-memory objects to
        ///   a byte representation.  Lifetime must exceed that of the
        ///   CollectConvReports object.
        ///
        /// \param[in] comm MPI communicator.  Should be the same as the
        ///   internal communicator in the \p packer object.
        explicit CollectConvReports(const Opm::Mpi::Packer&            packer,
                                    const Opm::Parallel::Communication comm)
            : Opm::Serializer<Opm::Mpi::Packer> { packer }
            , comm_                             { comm }
        {}

        /// Collect local convergence reports from each MPI rank.
        ///
        /// \param[in] report Local convergence report from the current MPI
        /// rank.
        ///
        /// \return Collection of local convergence reports from all MPI
        /// ranks in the current communicator.  Report objects returned in
        /// rank order.
        std::vector<Opm::ConvergenceReport>
        operator()(const Opm::ConvergenceReport& report);

    private:
        /// MPI communicator.
        Opm::Parallel::Communication comm_;

        /// Byte representation of local convergence report objects from all
        /// ranks.
        ///
        /// Only valid during a call to operator().
        std::vector<char> rankBuffers_{};

        /// Start pointers into rankBuffers_ for each rank's local
        /// convergence report object.
        std::vector<int> rankStart_{};

        /// Reconstitute a convergence report from byte stream representation.
        ///
        /// \param[in] rank MPI rank for which to reconstitute the local
        ///   convergence report object.
        ///
        /// \param[in,out] report \p rank's local convergence report object.
        ///   On entry, a default constructed object usable as a destination
        ///   for deserialisation.  On exit, a fully populated convergence
        ///   report object filled from the byte stream of \p rank.
        void deserialise(const std::vector<int>::size_type rank,
                         Opm::ConvergenceReport&           report);
    };

    std::vector<Opm::ConvergenceReport>
    CollectConvReports::operator()(const Opm::ConvergenceReport& report)
    {
        auto rankReports = std::vector<Opm::ConvergenceReport>(this->comm_.size());

        this->pack(report);

        std::tie(this->rankBuffers_, this->rankStart_) =
            Opm::allGatherv(this->m_buffer, this->comm_);

        auto rank = std::vector<int>::size_type{0};
        for (auto& rankReport : rankReports) {
            this->deserialise(rank++, rankReport);
        }

        return rankReports;
    }

    void CollectConvReports::deserialise(const std::vector<int>::size_type rank,
                                         Opm::ConvergenceReport&           report)
    {
        auto begin = this->rankBuffers_.begin() + this->rankStart_[rank + 0];
        auto end   = this->rankBuffers_.begin() + this->rankStart_[rank + 1];

        this->m_buffer.assign(begin, end);
        this->m_packSize = std::distance(begin, end);

        this->unpack(report);
    }
} // Anonymous namespace

namespace Opm
{
    /// Create a global convergence report combining local (per-process)
    /// reports.
    ConvergenceReport
    gatherConvergenceReport(const ConvergenceReport& local_report,
                            Parallel::Communication  mpi_communicator)
    {
        if (mpi_communicator.size() == 1) {
            // Sequential run, no communication needed.
            return local_report;
        }

        // Multi-process run (common case).  Need object distribution.
        auto combinedReport = ConvergenceReport {};

        const auto packer = Mpi::Packer { mpi_communicator };

        const auto rankReports =
            CollectConvReports { packer, mpi_communicator }(local_report);

        std::for_each(rankReports.begin(), rankReports.end(),
                      [&combinedReport](const auto& rankReport)
                      { combinedReport += rankReport; });

        return combinedReport;
    }

} // namespace Opm

#else // !HAVE_MPI

namespace Opm
{
    ConvergenceReport
    gatherConvergenceReport(const ConvergenceReport& local_report,
                            [[maybe_unused]] Parallel::Communication mpi_communicator)
    {
        return local_report;
    }
} // namespace Opm

#endif // HAVE_MPI
