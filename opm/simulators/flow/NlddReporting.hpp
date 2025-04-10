/*
  Copyright 2025, SINTEF Digital

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

#ifndef OPM_NLDD_REPORTING_HEADER_INCLUDED
#define OPM_NLDD_REPORTING_HEADER_INCLUDED

#include <dune/grid/common/gridenums.hh>
#include <dune/grid/common/partitionset.hh>

#include <opm/common/ErrorMacros.hpp>
#include <opm/common/OpmLog/OpmLog.hpp>
#include <opm/simulators/timestepping/SimulatorReport.hpp>
#include <opm/simulators/utils/DeferredLogger.hpp>
#include <opm/simulators/utils/ParallelCommunication.hpp>

#include <fmt/format.h>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <limits>
#include <numeric>
#include <sstream>
#include <vector>

namespace Opm {

/**
 * Reports NLDD statistics after simulation.
 *
 * @param domain_reports The accumulated reports per domain
 * @param local_report The accumulated reports per rank
 * @param output_cout Whether to output to cout
 * @param comm The communication object for parallel runs
 */
void reportNlddStatistics(const std::vector<SimulatorReport>& domain_reports,
                          const SimulatorReport& local_report,
                          const bool output_cout,
                          const Parallel::Communication& comm);

/**
 * Writes the number of nonlinear iterations per cell to a file in ResInsight compatible format
 *
 * @param odir The output directory
 * @param domains The subdomains
 * @param domain_reports The accumulated reports per domain
 * @param grid The grid
 * @param elementMapper The element mapper
 * @param cartMapper The cartesian index mapper
 */
template <class Grid, class Domain, class ElementMapper, class CartMapper>
void writeNonlinearIterationsPerCell(
    const std::filesystem::path& odir,
    const std::vector<Domain>& domains,
    const std::vector<SimulatorReport>& domain_reports,
    const Grid& grid,
    const ElementMapper& elementMapper,
    const CartMapper& cartMapper)
{
    const auto& dims = cartMapper.cartesianDimensions();
    const auto total_size = dims[0] * dims[1] * dims[2];
    const auto& comm = grid.comm();
    const int rank = comm.rank();

    // Create a cell-to-iterations mapping for this process
    std::vector<int> cell_iterations(grid.size(0), 0);

    // Populate the mapping with iteration counts for each domain
    const auto ds = domains.size();
    for (auto domain_idx = 0*ds; domain_idx < ds; ++domain_idx) {
        const auto& domain = domains[domain_idx];
        const auto& report = domain_reports[domain_idx];
        const int iterations = report.success.total_newton_iterations + 
                              report.failure.total_newton_iterations;

        for (const int cell_idx : domain.cells) {
            cell_iterations[cell_idx] = iterations;
        }
    }

    // Create a full-sized vector initialized with zeros (indicating inactive cells)
    auto full_iterations = std::vector<int>(total_size, 0);

    // Convert local cell indices to cartesian indices
    const auto& gridView = grid.leafGridView();
    for (const auto& cell : elements(gridView, Dune::Partitions::interior)) {
        const int cell_idx = elementMapper.index(cell);
        const int cart_idx = cartMapper.cartesianIndex(cell_idx);
        full_iterations[cart_idx] = cell_iterations[cell_idx];
    }

    // Gather all iteration data using max operation
    comm.max(full_iterations.data(), full_iterations.size());

    // Only rank 0 writes the file
    if (rank == 0) {
        auto fname = odir / "ResInsight_nonlinear_iterations.txt";
        std::ofstream resInsightFile { fname };
        // Write header
        resInsightFile << "NLDD_ITER" << '\n';

        // Write all cells, including inactive ones
        for (const auto& val : full_iterations) {
            resInsightFile << val << '\n';
        }
        resInsightFile << "/" << '\n';
    }
}

/**
 * Writes the partition vector to a file in ResInsight compatible format and a partition file for each rank
 *
 * @param odir The output directory
 * @param domains Vector of domains
 * @param grid The grid
 * @param elementMapper The element mapper
 * @param cartMapper The cartesian index mapper
 */
template <class Grid, class Domain, class ElementMapper, class CartMapper>
void writePartitions(
    const std::filesystem::path& odir,
    const std::vector<Domain>& domains,
    const Grid& grid,
    const ElementMapper& elementMapper,
    const CartMapper& cartMapper)
{
    const auto& dims = cartMapper.cartesianDimensions();
    const auto total_size = dims[0] * dims[1] * dims[2];
    const auto& comm = grid.comm();
    const int rank = comm.rank();

    const auto& partition_vector = reconstitutePartitionVector(domains, grid);

    // Create a full-sized vector initialized with -1 (indicating inactive cells)
    auto full_partition = std::vector<int>(total_size, -1);

    // Fill active cell values for this rank
    auto i = 0;
    for (const auto& cell : elements(grid.leafGridView(), Dune::Partitions::interior)) {
        full_partition[cartMapper.cartesianIndex(elementMapper.index(cell))] = partition_vector[i++];
    }

    // Gather all partitions using max operation
    comm.max(full_partition.data(), full_partition.size());

    // Only rank 0 writes the file
    if (rank == 0) {
        auto fname = odir / "ResInsight_compatible_partition.txt";
        std::ofstream resInsightFile { fname };

        // Write header
        resInsightFile << "NLDD_DOM" << '\n';

        // Write all cells, including inactive ones
        for (const auto& val : full_partition) {
            resInsightFile << val << '\n';
        }
        resInsightFile << "/" << '\n';
    }

    const auto nDigit = 1 + static_cast<int>(std::floor(std::log10(comm.size())));
    auto partition_fname = odir / fmt::format("{1:0>{0}}", nDigit, rank);
    std::ofstream pfile { partition_fname };

    auto cell_index = 0;
    for (const auto& cell : elements(grid.leafGridView(), Dune::Partitions::interior)) {
        pfile << rank << ' '
              << cartMapper.cartesianIndex(elementMapper.index(cell)) << ' '
              << partition_vector[cell_index++] << '\n';
    }
}

/**
 * Prints a summary of domain distribution across ranks
 *
 * @param partition_vector The partition vector
 * @param domains The subdomains
 * @param local_reports_accumulated The accumulated reports per rank
 * @param domain_reports_accumulated The accumulated reports per domain
 * @param grid The grid
 * @param num_wells The number of wells
 */
template <class Grid, class Domain>
void printDomainDistributionSummary(
    const std::vector<int>& partition_vector,
    const std::vector<Domain>& domains,
    SimulatorReport& local_reports_accumulated,
    std::vector<SimulatorReport>& domain_reports_accumulated,
    const Grid& grid,
    int num_wells)
{
    const auto& gridView = grid.leafGridView();
    const auto& comm = grid.comm();
    const int rank = comm.rank();

    const int num_domains = domains.size();
    const int owned_cells = partition_vector.size();

    // Count overlap cells using grid view iteration
    int overlap_cells = std::count_if(elements(gridView).begin(), elements(gridView).end(),
                                      [](const auto& cell) { return cell.partitionType() == Dune::OverlapEntity; });

    // Store data for summary output
    local_reports_accumulated.success.num_wells = num_wells;
    local_reports_accumulated.success.num_domains = num_domains;
    local_reports_accumulated.success.num_overlap_cells = overlap_cells;
    local_reports_accumulated.success.num_owned_cells = owned_cells;

    // Set statistics for each domain report
    const auto dr_size = domain_reports_accumulated.size();
    for (auto i = 0*dr_size; i < dr_size; ++i) {
        auto& domain_report = domain_reports_accumulated[i];
        domain_report.success.num_domains = 1;
        domain_report.success.num_overlap_cells = 0;
        domain_report.success.num_owned_cells = domains[i].cells.size();
    }

    // Gather data from all ranks
    std::vector<int> all_owned(comm.size());
    std::vector<int> all_overlap(comm.size());
    std::vector<int> all_wells(comm.size());
    std::vector<int> all_domains(comm.size());

    comm.gather(&owned_cells, all_owned.data(), 1, 0);
    comm.gather(&overlap_cells, all_overlap.data(), 1, 0);
    comm.gather(&num_wells, all_wells.data(), 1, 0);
    comm.gather(&num_domains, all_domains.data(), 1, 0);

    if (rank == 0) {
        std::ostringstream ss;
        ss << "\nNLDD domain distribution summary:\n"
        << "  rank   owned cells   overlap cells   total cells   wells   domains\n"
        << "--------------------------------------------------------------------\n";

        int total_owned = 0;
        int total_overlap = 0;
        int total_wells = 0;
        int total_domains = 0;

        for (int r = 0; r < comm.size(); ++r) {
            ss << std::setw(6) << r
            << std::setw(13) << all_owned[r]
            << std::setw(15) << all_overlap[r]
            << std::setw(14) << (all_owned[r] + all_overlap[r])
            << std::setw(8) << all_wells[r]
            << std::setw(9) << all_domains[r] << '\n';

            total_owned += all_owned[r];
            total_overlap += all_overlap[r];
            total_wells += all_wells[r];
            total_domains += all_domains[r];
        }

        ss << "--------------------------------------------------------------------\n"
        << "   sum"
        << std::setw(13) << total_owned
        << std::setw(15) << total_overlap
        << std::setw(14) << (total_owned + total_overlap)
        << std::setw(8) << total_wells
        << std::setw(9) << total_domains << '\n';

        OpmLog::info(ss.str());
    }
}

/**
 * Reconstructs the partition vector that maps each grid cell to its corresponding domain ID,
 * accounting for domain distribution across MPI ranks.
 *
 * @param domains Vector of domains
 * @param grid The grid
 * @return The reconstructed partition vector
 */
template <class Grid, class Domain>
std::vector<int> reconstitutePartitionVector(
    const std::vector<Domain>& domains,
    const Grid& grid)
{
    const auto& comm = grid.comm();
    const int rank = comm.rank();

    auto numD = std::vector<int>(comm.size() + 1, 0);
    numD[rank + 1] = static_cast<int>(domains.size());
    comm.sum(numD.data(), numD.size());
    std::partial_sum(numD.begin(), numD.end(), numD.begin());

    auto p = std::vector<int>(grid.size(0));
    auto maxCellIdx = std::numeric_limits<int>::min();

    auto d = numD[rank];
    for (const auto& domain : domains) {
        for (const auto& cell : domain.cells) {
            p[cell] = d;
            if (cell > maxCellIdx) {
                maxCellIdx = cell;
            }
        }

        ++d;
    }

    p.erase(p.begin() + maxCellIdx + 1, p.end());
    return p;
}

} // namespace Opm

#endif // OPM_NLDD_REPORTING_HEADER_INCLUDED
