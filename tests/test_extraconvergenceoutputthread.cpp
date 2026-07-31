/*
  Copyright 2026 Equinor ASA.

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

#define BOOST_TEST_MODULE TestExtraConvergenceOutputThread

#include <boost/test/unit_test.hpp>

#include <opm/simulators/flow/ConvergenceOutputConfiguration.hpp>
#include <opm/simulators/flow/ExtraConvergenceOutputThread.hpp>

#include <opm/simulators/timestepping/ConvergenceReport.hpp>

#include <filesystem>
#include <fstream>
#include <iterator>
#include <string>
#include <string_view>
#include <thread>
#include <vector>

namespace {

    std::string_view phaseName(const int comp)
    {
        return (comp == 0) ? "Water" : "Oil";
    }

    class OutputDir
    {
    public:
        explicit OutputDir(const std::string& name)
            : dir_ { std::filesystem::temp_directory_path() / name }
        {
            std::filesystem::remove_all(this->dir_);
            std::filesystem::create_directories(this->dir_);
        }

        ~OutputDir()
        {
            std::filesystem::remove_all(this->dir_);
        }

        std::string path() const { return this->dir_.string(); }

        std::vector<std::string> infoIterLines(const std::string& baseName) const
        {
            auto is = std::ifstream { this->dir_ / (baseName + ".INFOITER") };
            auto lines = std::vector<std::string>{};

            for (auto line = std::string{}; std::getline(is, line); ) {
                lines.push_back(line);
            }

            return lines;
        }

    private:
        std::filesystem::path dir_;
    };

    // Run the output thread against 'requests', then signal end of
    // production and join.
    void runOutputThread(const OutputDir&                                          dir,
                         const std::string&                                        baseName,
                         std::vector<Opm::ConvergenceReportQueue::OutputRequest>&&  requests)
    {
        auto queue = Opm::ConvergenceReportQueue{};
        auto writer = Opm::ConvergenceOutputThread {
            dir.path(), baseName, &phaseName,
            [](const double t) { return t; },
            Opm::ConvergenceOutputConfiguration{"steps,iterations"},
            queue
        };

        auto thread = std::thread {
            &Opm::ConvergenceOutputThread::writeASynchronous, &writer
        };

        if (! requests.empty()) {
            queue.enqueue(std::move(requests));
        }

        queue.signalLastOutputRequest();
        thread.join();
    }

    Opm::ConvergenceReport makeReport()
    {
        auto report = Opm::ConvergenceReport { 0.0 };

        report.setReservoirConvergenceMetric
            (Opm::ConvergenceReport::ReservoirFailure::Type::MassBalance,
             0, 1.0e-6, 1.0e-5);

        report.setCnvPoreVolSplit({{1.0, 0.0, 0.0}, {1, 0, 0}}, 1.0);

        return report;
    }

} // Anonymous namespace

// The end-of-production sentinel is an OutputRequest with no convergence
// reports.  A run that fails before completing a single non-linear iteration
// sends it as the very first request, so the output thread must not try to
// form the column header from it.
BOOST_AUTO_TEST_CASE(SentinelOnlyDoesNotWriteHeader)
{
    const auto dir = OutputDir { "test_extraconvergenceoutputthread_sentinel" };

    runOutputThread(dir, "CASE", {});

    BOOST_CHECK(dir.infoIterLines("CASE").empty());
}

BOOST_AUTO_TEST_CASE(HeaderAndOneIteration)
{
    const auto dir = OutputDir { "test_extraconvergenceoutputthread_normal" };

    auto requests = std::vector<Opm::ConvergenceReportQueue::OutputRequest>{};
    requests.push_back({ 0, 0, { makeReport() } });

    runOutputThread(dir, "CASE", std::move(requests));

    const auto lines = dir.infoIterLines("CASE");

    BOOST_REQUIRE_EQUAL(lines.size(), std::size_t{2});
    BOOST_CHECK(lines[0].find("ReportStep") != std::string::npos);
    BOOST_CHECK(lines[1].find("CONV") != std::string::npos);
}
