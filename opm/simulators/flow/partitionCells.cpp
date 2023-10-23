/*
  Copyright 2021 Total SE

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

#include <opm/simulators/flow/partitionCells.hpp>

#include <opm/simulators/flow/countGlobalCells.hpp>

#include <opm/simulators/utils/compressPartition.hpp>

#if HAVE_MPI && HAVE_ZOLTAN
#include <opm/simulators/utils/DeferredLogger.hpp>
#include <opm/simulators/utils/DeferredLoggingErrorHelpers.hpp>
#include <opm/simulators/utils/gatherDeferredLogger.hpp>
#include <opm/simulators/utils/ParallelNLDDPartitioningZoltan.hpp>
#include <opm/simulators/utils/SetupZoltanParams.hpp>

#include <opm/input/eclipse/Schedule/Well/Connection.hpp>
#include <opm/input/eclipse/Schedule/Well/Well.hpp>
#include <opm/input/eclipse/Schedule/Well/WellConnections.hpp>
#endif // HAVE_MPI && HAVE_ZOLTAN

#include <opm/grid/CpGrid.hpp>
#include <opm/grid/polyhedralgrid.hh>

#if HAVE_DUNE_ALUGRID
#include <dune/alugrid/grid.hh>
#endif // HAVE_DUNE_ALUGRID

#include <opm/common/utility/CSRGraphFromCoordinates.hpp>

#include <fmt/format.h>

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <numeric>
#include <stdexcept>
#include <string>
#include <string_view>
#include <tuple>
#include <type_traits>
#include <unordered_map>
#include <utility>
#include <vector>

// ===========================================================================
// Zoltan-based partitioning if available in this build configuration.
// ---------------------------------------------------------------------------

#if HAVE_MPI && HAVE_ZOLTAN

namespace {

/// Partition rank's interior cells into non-overlapping domains using the
/// Zoltan graph partitioning software package.
class ZoltanPartitioner
{
public:
    /// Constructor.
    ///
    /// \tparam GridView DUNE grid view type
    ///
    /// \param[in] grid_view Current rank's reachable cells.
    ///
    /// \param[in] local_to_global Callback for mapping rank's local cells
    ///   to globally unique cell IDs.  May for instance return the
    ///   Cartesian cell ID of each local cell.
    template <class GridView>
    explicit ZoltanPartitioner(const GridView&                                          grid_view,
                               const Opm::ParallelNLDDPartitioningZoltan::GlobalCellID& local_to_global);

    /// Form internal connectivity graph between rank's reachable, interior
    /// cells.
    ///
    /// \tparam GridView DUNE grid view type
    ///
    /// \tparam Element Grid view's representation of a cell.
    ///
    /// \param[in] grid_view Current rank's reachable cells.
    ///
    /// \param[in] wells Collection of wells to which this rank's partition
    ///   vector must be adapted.  Wells must not intersect multiple local
    ///   domains/blocks.
    ///
    /// \param[in] zoltan_ctrl Control parameters for on-rank subdomain
    ///   partitioning.
    template <class GridView, class Element>
    void buildLocalGraph(const GridView&                                grid_view,
                         const std::vector<Opm::Well>&                  wells,
                         const Opm::ZoltanPartitioningControl<Element>& zoltan_ctrl);

    /// Partition rank's interior cells into non-overlapping domains using
    /// the Zoltan graph partitioning software package.
    ///
    /// \tparam GridView DUNE grid view type
    ///
    /// \tparam Element Grid view's representation of a cell.
    ///
    /// \param[in] num_local_domains Target number of domains into which the
    ///   full model will be partitioned.  Note: This number is treated as a
    ///   global parameter across all MPI ranks and is not used as a number
    ///   of domains per rank.  The final number of model subdomains may
    ///   differ from this number when adapting to well constraints and if,
    ///   for whatever reason, Zoltan happens to generate a global
    ///   partitioning for which one domain/block happens to cross process
    ///   boundaries.
    ///
    /// \param[in] grid_view Current rank's reachable cells.
    ///
    /// \param[in] zoltan_ctrl Control parameters for on-rank subdomain
    ///   partitioning.
    ///
    /// \return Partition vector for rank's interior cells and number of
    ///   unique domains across those interior cells.
    template <class GridView, class Element>
    std::pair<std::vector<int>, int>
    partition(const int                                      num_local_domains,
              const GridView&                                grid_view,
              const Opm::ZoltanPartitioningControl<Element>& zoltan_ctrl) const;

private:
    /// Zoltan partitioning backend.
    Opm::ParallelNLDDPartitioningZoltan partitioner_;

    /// Form internal connectivity graph between rank's reachable, interior
    /// cells.
    ///
    /// \tparam GridView DUNE grid view type
    ///
    /// \tparam Element Grid view's representation of a cell.
    ///
    /// \param[in] grid_view Current rank's reachable cells.
    ///
    /// \param[in] zoltan_ctrl Control parameters for on-rank subdomain
    ///   partitioning.
    ///
    /// \return Mapping from globally unique cell IDs to local, on-rank
    ///   active cell IDs.  Effectively the cached result of \code
    ///   zoltan_ctrl.index() \endcode for each globally unique cell ID.
    template <class GridView, class Element>
    std::unordered_map<int, int>
    connectElements(const GridView&                                grid_view,
                    const Opm::ZoltanPartitioningControl<Element>& zoltan_ctrl);

    /// Adapt internal connectivity graph to account for connections induced
    /// by well reservoir connections.
    ///
    /// \tparam Comm MPI communication type.
    ///
    /// \param[in] comm MPI communication object.  Needed to get a global
    ///   view of any exceptions thrown in the implementation.
    ///
    /// \param[in] wells Collection of wells to which this rank's partition
    ///   vector must be adapted.  Wells must not intersect multiple local
    ///   domains/blocks.
    ///
    /// \param[in] g2l Mapping from globally unique cell IDs to local,
    ///   on-rank active cell IDs.  Return value from \c connectElements().
    template <typename Comm>
    void connectWells(const Comm                          comm,
                      const std::vector<Opm::Well>&       wells,
                      const std::unordered_map<int, int>& g2l);
};

// Note: "grid_view.size(0)" is intentional here.  It is not an error.  The
// ParallelNLDDPartitioningZoltan type will internally enumerate interior
// cells in a consecutive sequence, regardless of the ordering in the
// grid_view object.
template <class GridView>
ZoltanPartitioner::ZoltanPartitioner(const GridView&                                          grid_view,
                                     const Opm::ParallelNLDDPartitioningZoltan::GlobalCellID& local_to_global)
    : partitioner_(grid_view.comm(), grid_view.size(0), local_to_global)
{}

template <class GridView, class Element>
void ZoltanPartitioner::buildLocalGraph(const GridView&                                grid_view,
                                        const std::vector<Opm::Well>&                  wells,
                                        const Opm::ZoltanPartitioningControl<Element>& zoltan_ctrl)
{
    this->connectWells(grid_view.comm(), wells, this->connectElements(grid_view, zoltan_ctrl));
}

template <class GridView, class Element>
std::pair<std::vector<int>, int>
ZoltanPartitioner::partition(const int                                      num_local_domains,
                             const GridView&                                grid_view,
                             const Opm::ZoltanPartitioningControl<Element>& zoltan_ctrl) const
{
    auto partition = std::pair<std::vector<int>, int>{};

    auto param = Opm::setupZoltanParams("graph");
    param.insert_or_assign("NUM_GLOBAL_PARTS", fmt::format("{}", num_local_domains));
    param.insert_or_assign("IMBALANCE_TOL"   , fmt::format("{}", zoltan_ctrl.domain_imbalance));

    const auto domainsAll = this->partitioner_.partitionElements(param);

    auto& [domains, num_domains] = partition;
    domains.reserve(domainsAll.size());
    num_domains = 0;

    for (const auto& cell : elements(grid_view, Dune::Partitions::interior)) {
        const auto domainId = domains.emplace_back(domainsAll[zoltan_ctrl.index(cell)]);

        if (domainId >= num_domains) {
            num_domains = domainId + 1;
        }
    }

    // Fix-up for an extreme case: Interior cells whose neighbours are all
    // on a different rank--i.e., interior cells which do not have on-rank
    // neighbours.  Make each such cell be a separate domain on this rank.
    // Don't remove this code unless you've taken remedial action elsewhere
    // to ensure the situation doesn't happen.  For what it's worth, we've
    // seen this case occur in practice when testing on field cases.
    for (auto& domainId : domains) {
        if (domainId < 0) {
            domainId = num_domains++;
        }
    }

    return partition;
}

template <class GridView, class Element>
std::unordered_map<int, int>
ZoltanPartitioner::connectElements(const GridView&                                grid_view,
                                   const Opm::ZoltanPartitioningControl<Element>& zoltan_ctrl)
{
    auto g2l = std::unordered_map<int, int>{};

    for (const auto& element : elements(grid_view, Dune::Partitions::interior)) {
        {
            const auto c = zoltan_ctrl.index(element);
            g2l.insert_or_assign(zoltan_ctrl.local_to_global(c), c);
        }

        for (const auto& is : intersections(grid_view, element)) {
            if (! is.neighbor()) {
                continue;
            }

            const auto& in  = is.inside();
            const auto& out = is.outside();

            if ((in .partitionType() != Dune::InteriorEntity) ||
                (out.partitionType() != Dune::InteriorEntity))
            {
                // Connect cells in interior only.
                continue;
            }

            const auto c1 = zoltan_ctrl.index(in);
            const auto c2 = zoltan_ctrl.index(out);

            this->partitioner_.registerConnection(c1, c2);
            this->partitioner_.registerConnection(c2, c1);
        }
    }

    return g2l;
}

template <typename Comm>
void ZoltanPartitioner::connectWells(const Comm                          comm,
                                     const std::vector<Opm::Well>&       wells,
                                     const std::unordered_map<int, int>& g2l)
{
    auto distributedWells = 0;

    for (const auto& well : wells) {
        auto cellIx = std::vector<int>{};
        auto otherProc = 0;

        for (const auto& conn : well.getConnections()) {
            auto locPos = g2l.find(conn.global_index());
            if (locPos == g2l.end()) {
                ++otherProc;
                continue;
            }

            cellIx.push_back(locPos->second);
        }

        if ((otherProc > 0) && !cellIx.empty()) {
            ++distributedWells;
            continue;
        }

        const auto nc = cellIx.size();
        for (auto ic1 = 0*nc; ic1 < nc; ++ic1) {
            for (auto ic2 = ic1 + 1; ic2 < nc; ++ic2) {
                this->partitioner_.registerConnection(cellIx[ic1], cellIx[ic2]);
                this->partitioner_.registerConnection(cellIx[ic2], cellIx[ic1]);
            }
        }

        if (! cellIx.empty()) {
            // All cells intersected by a single well go in the same domain.
            this->partitioner_.forceSameDomain(std::move(cellIx));
        }
    }

    OPM_BEGIN_PARALLEL_TRY_CATCH()

    if (distributedWells > 0) {
        const auto* pl = (distributedWells == 1) ? "" : "s";

        throw std::invalid_argument {
            fmt::format(R"({} distributed well{} detected on rank {}.
This is not supported in the current NLDD implementation.)",
                        distributedWells, pl, comm.rank())
        };
    }

    OPM_END_PARALLEL_TRY_CATCH(std::string { "ZoltanPartitioner::connectWells()" }, comm)
}

template <class GridView, class Element>
std::pair<std::vector<int>, int>
partitionCellsZoltan(const int                                      num_domains,
                     const GridView&                                grid_view,
                     const std::vector<Opm::Well>&                  wells,
                     const Opm::ZoltanPartitioningControl<Element>& zoltan_ctrl)
{
    if (num_domains <= 1) {     // No partitioning => every cell in domain zero.
        const auto num_interior_cells =
            Opm::detail::countLocalInteriorCellsGridView(grid_view);

        return {
            std::piecewise_construct,
            std::forward_as_tuple(num_interior_cells, 0),
            std::forward_as_tuple(1)
        };
    }

    auto partitioner = ZoltanPartitioner { grid_view, zoltan_ctrl.local_to_global };
    partitioner.buildLocalGraph(grid_view, wells, zoltan_ctrl);

    return partitioner.partition(num_domains, grid_view, zoltan_ctrl);
}

} // Anonymous namespace

#endif // HAVE_MPI && HAVE_ZOLTAN

namespace {
    std::vector<int> integerVectorFromFile(const std::filesystem::path& fileName)
    {
        std::ifstream is(fileName);

        return { std::istream_iterator<int>(is), std::istream_iterator<int>{} };
    }
} // Anonymous namespace

#if HAVE_MPI

namespace {

    template <class Communicator>
    std::vector<int>
    groupInputByProcessID(const std::vector<int>& rawInput,
                          const Communicator&     comm)
    {
        auto groupedCells = std::vector<int>{};
        auto startPtr = std::vector<int>(comm.size() + 1, 0);

        if (comm.rank() == 0) {
            // VertexID = std::vector<int>::size_type
            // TrackCompressedIdx = false
            // PermitSelfConnections = true (e.g., process 0 has row 0).
            auto cellGroup = Opm::utility::CSRGraphFromCoordinates<
                std::vector<int>::size_type, false, true
            >{};

            const auto numCells = rawInput.size() / 3;
            for (auto cell = 0*numCells; cell < numCells; ++cell) {
                cellGroup.addConnection(rawInput[3*cell + 0], cell);
            }
            cellGroup.compress(comm.size());

            const auto& inputIx = cellGroup.columnIndices();
            groupedCells.reserve(2 * inputIx.size());

            for (const auto& cell : inputIx) {
                const auto* cellVal = &rawInput[3*cell + 0];
                groupedCells.push_back(cellVal[1]); // Cartesian index
                groupedCells.push_back(cellVal[2]); // Block ID
            }

            const auto& start = cellGroup.startPointers();
            std::copy(start.begin(), start.end(), startPtr.begin());
        }

        comm.broadcast(startPtr.data(), comm.size() + 1, 0);

        // We're scattering two ints per cell
        for (auto& startIx : startPtr) {
            startIx *= 2;
        }

        auto sendLength = std::vector<int>(comm.size());
        std::adjacent_difference(startPtr.begin() + 1,
                                 startPtr.end(),
                                 sendLength.begin());

        auto perProcessInput = std::vector<int>(sendLength[comm.rank()], 0);
        comm.scatterv(groupedCells.data(),
                      sendLength.data(),
                      startPtr.data(),
                      perProcessInput.data(),
                      static_cast<int>(perProcessInput.size()),
                      0);

        return perProcessInput;
    }

    template <class GridView, class Element>
    std::unordered_map<int,int>
    interiorLocalToGlobalIndexMap(const GridView&                                gridView,
                                  const std::size_t                              numInterior,
                                  const Opm::ZoltanPartitioningControl<Element>& zoltanCtrl)
    {
        auto g2l = std::unordered_map<int,int>{};

        auto localCell = 0;
        for (const auto& cell : elements(gridView, Dune::Partitions::interior)) {
            const auto globIx = zoltanCtrl.local_to_global(zoltanCtrl.index(cell));
            g2l.insert_or_assign(globIx, localCell++);
        }

        OPM_BEGIN_PARALLEL_TRY_CATCH()

        if (g2l.size() != numInterior) {
            throw std::invalid_argument {
                fmt::format("Local-to-global mapping is not bijective on rank {}.",
                            gridView.comm().rank())
            };
        }

        OPM_END_PARALLEL_TRY_CATCH(std::string { "scatterPartition()/UniqueL2G" },
                                   gridView.comm())

        return g2l;
    }

    template <class Communicator>
    std::vector<int>
    interiorPartitionInput(const std::vector<int>& rawInput,
                           const Communicator&     comm,
                           const std::size_t       numInterior)
    {
        auto myPartitionInput = groupInputByProcessID(rawInput, comm);

        OPM_BEGIN_PARALLEL_TRY_CATCH()

        if (myPartitionInput.size() != 2*numInterior) {
            throw std::out_of_range {
                fmt::format("Input partition of size {} does not match "
                            "the number of interior cells ({}) on rank {}.",
                            myPartitionInput.size(),
                            numInterior,
                            comm.rank())
            };
        }

        OPM_END_PARALLEL_TRY_CATCH(std::string { "scatterPartition()/InputSize" },
                                   comm)

        return myPartitionInput;
    }

    template <class GridView, class Element>
    std::vector<int>
    scatterPartition(const std::vector<int>&                        rawInput,
                     const GridView&                                gridView,
                     const std::size_t                              numInterior,
                     const Opm::ZoltanPartitioningControl<Element>& zoltanCtrl)
    {
        auto partition = std::vector<int>(numInterior, -1);

        const auto g2l = interiorLocalToGlobalIndexMap(gridView, numInterior, zoltanCtrl);
        const auto myPartitionInput =
            interiorPartitionInput(rawInput, gridView.comm(), numInterior);

        for (auto i = 0*numInterior; i < numInterior; ++i) {
            auto cellPos = g2l.find(myPartitionInput[2*i + 0]);
            if (cellPos != g2l.end()) {
                partition[cellPos->second] = myPartitionInput[2*i + 1];
            }
        }

        OPM_BEGIN_PARALLEL_TRY_CATCH()

        const auto allCovered =
            std::all_of(partition.begin(), partition.end(),
                        [](const int d) { return d >= 0; });

        if (! allCovered) {
            throw std::out_of_range {
                fmt::format("Input partition does not cover "
                            "all {} interior cells on rank {}.",
                            numInterior, gridView.comm().rank())
            };
        }

        OPM_END_PARALLEL_TRY_CATCH(std::string { "scatterPartition()/CellCoverage" },
                                   gridView.comm())

        return partition;
    }

    template <class GridView, class Element>
    std::pair<std::vector<int>, int>
    partitionCellsFromFileMPI(const std::filesystem::path&                   fileName,
                              const GridView&                                grid_view,
                              const Opm::ZoltanPartitioningControl<Element>& zoltan_ctrl)
    {
        const auto input = (grid_view.comm().rank() == 0)
            ? integerVectorFromFile(fileName)
            : std::vector<int>{};

        const auto num_interior = Opm::detail::
            countLocalInteriorCellsGridView(grid_view);

        const auto num_cells = grid_view.comm().sum(num_interior);

        OPM_BEGIN_PARALLEL_TRY_CATCH()

        if ((grid_view.comm().rank() == 0) &&
            (input.size() != 3 * num_cells))
        {
            throw std::out_of_range {
                fmt::format("Number of partition file elements {} does not "
                            "match model's number of active cells ({})",
                            input.size(), num_cells)
            };
        }

        OPM_END_PARALLEL_TRY_CATCH(std::string { "partitionCellsFromFileMPI()" },
                                   grid_view.comm())

        return Opm::util::compressAndCountPartitionIDs
            (scatterPartition(input, grid_view, num_interior, zoltan_ctrl));
    }
} // Anonymous namespace

#endif // HAVE_MPI

// ===========================================================================

template <class GridView, class Element>
std::pair<std::vector<int>, int>
Opm::partitionCells(const std::string& method,
                    const int          num_local_domains,
                    const GridView&    grid_view,
                    [[maybe_unused]] const std::vector<Well>&                  wells,
                    [[maybe_unused]] const ZoltanPartitioningControl<Element>& zoltan_ctrl)
{
    if (method == "zoltan") {
#if HAVE_MPI && HAVE_ZOLTAN

        return partitionCellsZoltan(num_local_domains, grid_view, wells, zoltan_ctrl);

#else // !HAVE_MPI || !HAVE_ZOLTAN

        OPM_THROW(std::runtime_error, "Zoltan requested for local domain partitioning, "
                  "but is not available in the current build configuration.");

#endif // HAVE_MPI && HAVE_ZOLTAN
    }
    else if (method == "simple") {
        const int num_cells = detail::countLocalInteriorCellsGridView(grid_view);
        return partitionCellsSimple(num_cells, num_local_domains);
    }
    else if (const auto ext = std::string_view { ".partition" };
             method.rfind(ext) == method.size() - ext.size())
    {
#if HAVE_MPI
        if (grid_view.comm().size() > 1) {
            return partitionCellsFromFileMPI(method, grid_view, zoltan_ctrl);
        }
#endif  // HAVE_MPI

        // Method string ends with ".partition", interpret as filename for partitioning.
        const int num_cells = detail::countLocalInteriorCellsGridView(grid_view);
        return partitionCellsFromFile(method, num_cells);
    }
    else {
        OPM_THROW(std::runtime_error, "Unknown local domain partitioning method requested: " + method);
    }
}

// Read from file, containing one or three numbers per cell.
std::pair<std::vector<int>, int>
Opm::partitionCellsFromFile(const std::string& filename, const int num_cells)
{
    // Read file into single vector.
    auto cellpart = integerVectorFromFile(filename);
    if (cellpart.size() == static_cast<std::size_t>(num_cells)) {
        return util::compressAndCountPartitionIDs(std::move(cellpart));
    }

    if (cellpart.size() == 3 * static_cast<std::size_t>(num_cells)) {
        auto p = std::vector<int>(num_cells);
        for (auto c = 0*num_cells; c < num_cells; ++c) {
            p[c] = cellpart[3*c + 2];
        }

        return util::compressAndCountPartitionIDs(std::move(p));
    }

    throw std::runtime_error {
        fmt::format("Partition file contains {} entries, "
                    "but there are {} cells.",
                    cellpart.size(), num_cells)
    };
}


// Trivially simple partitioner
std::pair<std::vector<int>, int>
Opm::partitionCellsSimple(const int num_cells, const int num_domains)
{
    // Build the partitions.
    const int dom_sz = num_cells / num_domains;
    std::vector<int> bounds(num_domains + 1, dom_sz);
    bounds[0] = 0;
    for (int i = 0; i < num_cells % num_domains; ++i) {
        ++bounds[i + 1];
    }
    std::partial_sum(bounds.begin(), bounds.end(), bounds.begin());
    assert(bounds[num_domains] == num_cells);
    std::vector<int> part(num_cells);
    for (int i = 0; i < num_domains; ++i) {
        std::fill(part.begin() + bounds[i], part.begin() + bounds[i + 1], i);
    }
    return { part, num_domains };
}

// ===========================================================================

// ---------------------------------------------------------------------------
// Explicit specialisations/instantiations of partitionCells template.
// Deliberately placed at end of file.  No other code beyond this separator.
// ---------------------------------------------------------------------------

#define InstantiatePartitionCells(Grid)                                 \
    template std::pair<std::vector<int>, int>                           \
    Opm::partitionCells(const std::string&,                             \
                        const int,                                      \
                        const std::remove_cv_t<std::remove_reference_t< \
                        decltype(std::declval<Grid>().leafGridView())>>&, \
                        const std::vector<Opm::Well>&,                  \
                        const Opm::ZoltanPartitioningControl<           \
                        typename std::remove_cv_t<std::remove_reference_t< \
                        decltype(std::declval<Grid>().leafGridView())>>::template Codim<0>::Entity>&)

// ---------------------------------------------------------------------------
// Grid types built into Flow.
// ---------------------------------------------------------------------------
InstantiatePartitionCells(Dune::CpGrid);

using PolyGrid = Dune::PolyhedralGrid<3, 3, double>;
InstantiatePartitionCells(PolyGrid);

// ---------------------------------------------------------------------------
// ALUGrid, if available.
// ---------------------------------------------------------------------------
#if HAVE_DUNE_ALUGRID
#if HAVE_MPI
using ALUGridComm = Dune::ALUGridMPIComm;
#else
using ALUGridComm = Dune::ALUGridNoComm;
#endif // HAVE_MPI

// ---------------------------------------------------------------------------
// Base case.
// ---------------------------------------------------------------------------
using BaseALUGrid = Dune::ALUGrid<3, 3, Dune::cube, Dune::nonconforming, ALUGridComm>;
InstantiatePartitionCells(BaseALUGrid);
#endif // HAVE_DUNE_ALUGRID

// ===========================================================================

#undef InstantiatePartitionCells

// ---------------------------------------------------------------------------
// End of file.  No statements below this separator.
// ---------------------------------------------------------------------------
