/*
  Copyright 2021,2025 SINTEF Digital, Mathematics and Cybernetics.
  Copyright 2025 Equinor ASA

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

#define BOOST_TEST_MODULE OPM_test_partitionCells
#define BOOST_TEST_NO_MAIN
#include <boost/test/unit_test.hpp>

#include <dune/grid/common/gridview.hh>

#include <opm/grid/common/WellConnections.hpp>
#include <opm/grid/CpGrid.hpp>
#include <opm/grid/utility/OpmWellType.hpp>

#include <opm/input/eclipse/Schedule/Well/Connection.hpp>
#include <opm/input/eclipse/Schedule/Well/Well.hpp>
#include <opm/input/eclipse/Schedule/Well/WellConnections.hpp>
#include <opm/input/eclipse/Parser/Parser.hpp>
#include <opm/input/eclipse/Deck/Deck.hpp>
#include <opm/input/eclipse/EclipseState/Grid/EclipseGrid.hpp>

#include <opm/simulators/flow/partitionCells.hpp>

// Helper functions
namespace {
    Opm::Connection createConnection(int i, int j, int k, int global_index) {
        return Opm::Connection(i, j, k,                           // i,j,k
                             global_index,                        // global_index
                             0,                                   // complnum
                             Opm::Connection::State::OPEN,        // state
                             Opm::Connection::Direction::Z,       // direction
                             Opm::Connection::CTFKind::DeckValue, // kind
                             0,                                   // satTableId
                             0.0,                                 // depth
                             Opm::Connection::CTFProperties(),    // properties
                             0,                                   // sort_value
                             false);                             // defaultSatTabId
    }

    Dune::cpgrid::OpmWellType createWell(const std::string& name) {
        using namespace Opm;
        return Dune::cpgrid::OpmWellType(name, name, 0, 0, 0, 0, 0.0, WellType(),
                                        Well::ProducerCMode(), Connection::Order::TRACK,
                                        UnitSystem::newMETRIC(),
                                        0.0, 0.0, false, false, 0, Well::GasInflowEquation());
    }

    std::vector<Dune::cpgrid::OpmWellType> createWellsWithConnections(
        const std::vector<std::pair<std::string, std::vector<int>>>& well_specs) {
        std::vector<Dune::cpgrid::OpmWellType> wells;
        for (const auto& [well_name, cell_indices] : well_specs) {
            auto well_conn = std::make_shared<Opm::WellConnections>();
            for (int idx : cell_indices) {
                well_conn->add(createConnection(idx, 0, 0, idx));
            }
            auto well = createWell(well_name);
            well.updateConnections(well_conn, true);
            wells.push_back(well);
        }
        return wells;
    }

    Dune::CpGrid createTestGrid(const std::array<int, 3>& dims,
                               const std::array<double, 3>& size) {
        Dune::CpGrid grid;
        grid.createCartesian(dims, size);
        return grid;
    }

    template<typename Entity>
    Opm::ZoltanPartitioningControl<Entity> createZoltanControl(const Dune::CpGrid& grid) {
        Opm::ZoltanPartitioningControl<Entity> zoltan_ctrl;
        zoltan_ctrl.domain_imbalance = 1.1;
        zoltan_ctrl.index = [](const auto& element) {
            return element.index();
        };
        zoltan_ctrl.local_to_global = [&grid](const int local_idx) {
            return grid.globalCell()[local_idx];
        };
        return zoltan_ctrl;
    }
}

BOOST_AUTO_TEST_CASE(FileBased)
{
    auto [part, num_part] = Opm::partitionCellsFromFile("test10.partition", 10);
    BOOST_CHECK_EQUAL(num_part, 3);
    std::vector<int> expected = { 0, 0, 1, 1, 2, 2, 1, 1, 0, 0 };
    BOOST_CHECK_EQUAL_COLLECTIONS(expected.begin(), expected.end(), part.begin(), part.end());
}

BOOST_AUTO_TEST_CASE(FileBasedWrongNumberOfCells)
{
    auto func = []() { auto [part, num_part] = Opm::partitionCellsFromFile("test10.partition", 11); };
    BOOST_CHECK_THROW(func(), std::runtime_error);
}

BOOST_AUTO_TEST_CASE(Simple1)
{
    auto [part, num_part]  = Opm::partitionCellsSimple(10, 3);
    BOOST_CHECK_EQUAL(num_part, 3);
    std::vector<int> expected = { 0, 0, 0, 0, 1, 1, 1, 2, 2, 2 };
    BOOST_CHECK_EQUAL_COLLECTIONS(expected.begin(), expected.end(), part.begin(), part.end());

}

BOOST_AUTO_TEST_CASE(Simple2)
{
    auto [part, num_part]  = Opm::partitionCellsSimple(10, 7);
    BOOST_CHECK_EQUAL(num_part, 7);
    std::vector<int> expected = { 0, 0, 1, 1, 2, 2, 3, 4, 5, 6 };
    BOOST_CHECK_EQUAL_COLLECTIONS(expected.begin(), expected.end(), part.begin(), part.end());
}


BOOST_AUTO_TEST_CASE(PartitionCellsTest)
{
    auto grid = createTestGrid({3, 4, 1}, {3.0, 4.0, 1.0});
    using Entity = typename Dune::CpGrid::LeafGridView::template Codim<0>::Entity;
    auto zoltan_ctrl = createZoltanControl<Entity>(grid);

    // Test simple partitioning
    {
        auto [part, num_part] = Opm::partitionCells("simple", 3, grid.leafGridView(),
                                                   std::vector<Opm::Well>{},
                                                   std::unordered_map<std::string, std::set<int>>{},
                                                   zoltan_ctrl,
                                                   0);
        BOOST_CHECK_EQUAL(num_part, 3);
        BOOST_CHECK_EQUAL(part.size(), 12);

        // Check that all partition numbers are valid
        for (const auto& p : part) {
            BOOST_CHECK(p >= 0 && p < 3);
        }
    }

#if HAVE_MPI && HAVE_ZOLTAN
    // Test Zoltan partitioning
    {
        auto [part, num_part] = Opm::partitionCells("zoltan", 3, grid.leafGridView(),
                                                   std::vector<Opm::Well>{},
                                                   std::unordered_map<std::string, std::set<int>>{},
                                                   zoltan_ctrl,
                                                   0);
        BOOST_CHECK_EQUAL(num_part, 3);
        BOOST_CHECK_EQUAL(part.size(), 12);

        // Check that all partition numbers are valid
        for (const auto& p : part) {
            BOOST_CHECK(p >= 0 && p < 3);
        }
    }
#endif
}

BOOST_AUTO_TEST_CASE(PartitionCellsWithWellMergeTest)
{
    // Create a 1D grid with 7 cells
    // Visual representation of the grid and wells:
    // Cell indices:    0    1    2    3    4    5    6
    // Well layout:     W1---W1   |    W2---W2   |    |
    // Expected:        Cell 0 and 1 should be merged
    //                  Cell 3 and 4 should be merged
    //                  Cells 2, 5, and 6 can be in any partition
    auto grid = createTestGrid({7, 1, 1}, {7.0, 1.0, 1.0});

    // Setup Zoltan control parameters
    using Entity = typename Dune::CpGrid::LeafGridView::template Codim<0>::Entity;
    auto zoltan_ctrl = createZoltanControl<Entity>(grid);

    auto wells = createWellsWithConnections({
        {"TESTW1", {0, 1}},    // Well 1 connects cells 0,1
        {"TESTW2", {3, 4}}     // Well 2 connects cells 3,4
    });

    // Create well connections object
    Dune::cpgrid::WellConnections wellConnections(wells, std::unordered_map<std::string, std::set<int>>{}, grid);

#if HAVE_MPI && HAVE_ZOLTAN
    // Test Zoltan partitioning with well cells
    auto [part, num_part] = Opm::partitionCells("zoltan", 4, grid.leafGridView(),
                                               wells,
                                               std::unordered_map<std::string, std::set<int>>{},
                                               zoltan_ctrl,
                                               0);

    // Verify number of partitions
    BOOST_CHECK_EQUAL(num_part, 4);
    BOOST_CHECK_EQUAL(part.size(), 7);

    // The key tests:
    // 1. Cells 0 and 1 (well1) should be in the same partition
    BOOST_CHECK_EQUAL(part[0], part[1]);

    // 2. Cells 3 and 4 (well2) should be in the same partition
    BOOST_CHECK_EQUAL(part[3], part[4]);

    // 3. Cells 2, 5, and 6 can be in any partition, but well1 and well2 cells
    // don't need to be in the same partition
#endif
}

BOOST_AUTO_TEST_CASE(PartitionCellsWithOverlappingWellsTest)
{
    // Create a 1D grid with 7 cells
    // Visual representation of the grid and wells:
    // Cell indices:    0    1    2    3    4    5    6
    // Well layout:          W1---W1---W1   |    |    |
    //                                W2---W2
    // Expected:        Cells 1,2,3 should be in same partition (Well 1)
    //                  Cells 3,4 should be in same partition (Well 2)
    //                  Therefore cells 1,2,3,4 should all end up in same partition

    auto grid = createTestGrid({7, 1, 1}, {7.0, 1.0, 1.0});
    using Entity = typename Dune::CpGrid::LeafGridView::template Codim<0>::Entity;
    auto zoltan_ctrl = createZoltanControl<Entity>(grid);

    auto wells = createWellsWithConnections({
        {"TESTW1", {1, 2, 3}}, // Well 1 connects cells 1,2,3
        {"TESTW2", {3, 4}}     // Well 2 connects cells 3,4
    });

#if HAVE_MPI && HAVE_ZOLTAN
    // Test Zoltan partitioning with overlapping well cells
    auto [part, num_part] = Opm::partitionCells("zoltan", 4, grid.leafGridView(),
                                               wells,
                                               std::unordered_map<std::string, std::set<int>>{},
                                               zoltan_ctrl,
                                               0);

    // Verify number of partitions
    BOOST_CHECK_EQUAL(num_part, 4);
    BOOST_CHECK_EQUAL(part.size(), 7);

    // The key tests:
    // 1. All cells connected by Well 1 should be in the same partition
    BOOST_CHECK_EQUAL(part[1], part[2]);
    BOOST_CHECK_EQUAL(part[2], part[3]);

    // 2. All cells connected by Well 2 should be in the same partition
    BOOST_CHECK_EQUAL(part[3], part[4]);

    // 3. Due to the overlap at cell 3, all well-connected cells should be in the same partition
    BOOST_CHECK_EQUAL(part[1], part[4]);

    // 4. Cells 0, 5, and 6 can be in any partition
#endif
}


BOOST_AUTO_TEST_CASE(PartitionCellsWithOverlappingWells3DTest)
{
    // Create a 3x3x3 grid with three wells:
    // - Well1: Vertical well through (1,1,0) -> (1,1,1) -> (1,1,2)
    // - Well2: Diagonal well through (0,0,0) -> (1,1,1) -> (2,2,2)  [overlaps with Well1 at (1,1,1)]
    // - Well3: Horizontal well through (0,2,0) -> (1,2,0) -> (2,2,0) [no overlap]

    auto grid = createTestGrid({3, 3, 3}, {3.0, 3.0, 3.0});
    using Entity = typename Dune::CpGrid::LeafGridView::template Codim<0>::Entity;
    auto zoltan_ctrl = createZoltanControl<Entity>(grid);

    // Helper to convert i,j,k coordinates to global index for 3x3x3 grid
    auto ijkToGlobal = [](int i, int j, int k) { return i + (3 * j) + (9 * k); };

    // Create wells with proper 3D connections
    auto wells = std::vector<Dune::cpgrid::OpmWellType>();

    // Well 1 - vertical well
    {
        auto well_conn = std::make_shared<Opm::WellConnections>();
        well_conn->add(createConnection(1, 1, 0, ijkToGlobal(1, 1, 0)));  // Global index: 4
        well_conn->add(createConnection(1, 1, 1, ijkToGlobal(1, 1, 1)));  // Global index: 13
        well_conn->add(createConnection(1, 1, 2, ijkToGlobal(1, 1, 2)));  // Global index: 22
        auto well = createWell("VERTICAL");
        well.updateConnections(well_conn, true);
        wells.push_back(well);
    }

    // Well 2 - diagonal well (overlaps with Well 1 at (1,1,1))
    {
        auto well_conn = std::make_shared<Opm::WellConnections>();
        well_conn->add(createConnection(0, 0, 0, ijkToGlobal(0, 0, 0)));  // Global index: 0
        well_conn->add(createConnection(1, 1, 1, ijkToGlobal(1, 1, 1)));  // Global index: 13
        well_conn->add(createConnection(2, 2, 2, ijkToGlobal(2, 2, 2)));  // Global index: 26
        auto well = createWell("DIAGONAL");
        well.updateConnections(well_conn, true);
        wells.push_back(well);
    }

    // Well 3 - horizontal well (no overlap)
    {
        auto well_conn = std::make_shared<Opm::WellConnections>();
        well_conn->add(createConnection(0, 2, 0, ijkToGlobal(0, 2, 0)));  // Global index: 6
        well_conn->add(createConnection(1, 2, 0, ijkToGlobal(1, 2, 0)));  // Global index: 7
        well_conn->add(createConnection(2, 2, 0, ijkToGlobal(2, 2, 0)));  // Global index: 8
        auto well = createWell("HORIZONTAL");
        well.updateConnections(well_conn, true);
        wells.push_back(well);
    }

    // Create well connections object
    Dune::cpgrid::WellConnections wellConnections(wells, std::unordered_map<std::string, std::set<int>>{}, grid);

#if HAVE_MPI && HAVE_ZOLTAN
    // Test Zoltan partitioning with overlapping well cells
    auto [part, num_part] = Opm::partitionCells("zoltan", 10, grid.leafGridView(),
                                               wells,
                                               std::unordered_map<std::string, std::set<int>>{},
                                               zoltan_ctrl,
                                               0);

    // Verify number of partitions
    BOOST_CHECK_EQUAL(num_part, 10);
    BOOST_CHECK_EQUAL(part.size(), 27);  // 3x3x3 grid

    // The key tests:
    // 1. All cells in Well 1 (vertical) should be in the same partition
    BOOST_CHECK_EQUAL(part[ijkToGlobal(1,1,0)], part[ijkToGlobal(1,1,1)]);
    BOOST_CHECK_EQUAL(part[ijkToGlobal(1,1,1)], part[ijkToGlobal(1,1,2)]);

    // 2. All cells in Well 2 (diagonal) should be in the same partition
    BOOST_CHECK_EQUAL(part[ijkToGlobal(0,0,0)], part[ijkToGlobal(1,1,1)]);
    BOOST_CHECK_EQUAL(part[ijkToGlobal(1,1,1)], part[ijkToGlobal(2,2,2)]);

    // 3. Due to overlap at (1,1,1), all cells from Well 1 and Well 2 should be in same partition
    BOOST_CHECK_EQUAL(part[ijkToGlobal(1,1,0)], part[ijkToGlobal(0,0,0)]);
    BOOST_CHECK_EQUAL(part[ijkToGlobal(1,1,2)], part[ijkToGlobal(2,2,2)]);

    // 4. All cells in Well 3 (horizontal) should be in the same partition
    BOOST_CHECK_EQUAL(part[ijkToGlobal(0,2,0)], part[ijkToGlobal(1,2,0)]);
    BOOST_CHECK_EQUAL(part[ijkToGlobal(1,2,0)], part[ijkToGlobal(2,2,0)]);

    // 5. Well 3 can be in a different partition from Well 1 and 2 since it doesn't overlap
    // No test needed as it can be in any partition
#endif
}

BOOST_AUTO_TEST_CASE(PartitionCellsComplexWellNetworkTest)
{
    // Create a 5x5x4 grid to accommodate a complex well network
    auto grid = createTestGrid({5, 5, 4}, {5.0, 5.0, 4.0});
    using Entity = typename Dune::CpGrid::LeafGridView::template Codim<0>::Entity;
    auto zoltan_ctrl = createZoltanControl<Entity>(grid);

    // Helper to convert i,j,k coordinates to global index for 5x5x4 grid
    auto ijkToGlobal = [](int i, int j, int k) { return i + (5 * j) + (25 * k); };

    // Create wells with complex connections
    auto wells = std::vector<Dune::cpgrid::OpmWellType>();

    // Well 1 - Main vertical producer with two branches
    {
        auto well_conn = std::make_shared<Opm::WellConnections>();
        // Main wellbore - vertical
        well_conn->add(createConnection(2, 2, 0, ijkToGlobal(2, 2, 0)));  // Bottom
        well_conn->add(createConnection(2, 2, 1, ijkToGlobal(2, 2, 1)));
        well_conn->add(createConnection(2, 2, 2, ijkToGlobal(2, 2, 2)));
        well_conn->add(createConnection(2, 2, 3, ijkToGlobal(2, 2, 3)));  // Top
        // Branch 1 - horizontal in x direction
        well_conn->add(createConnection(3, 2, 2, ijkToGlobal(3, 2, 2)));
        well_conn->add(createConnection(4, 2, 2, ijkToGlobal(4, 2, 2)));
        // Branch 2 - horizontal in y direction
        well_conn->add(createConnection(2, 3, 1, ijkToGlobal(2, 3, 1)));
        well_conn->add(createConnection(2, 4, 1, ijkToGlobal(2, 4, 1)));
        auto well = createWell("PRODUCER1");
        well.updateConnections(well_conn, true);
        wells.push_back(well);
    }

    // Well 2 - Diagonal injector crossing the producer
    {
        auto well_conn = std::make_shared<Opm::WellConnections>();
        well_conn->add(createConnection(0, 0, 0, ijkToGlobal(0, 0, 0)));
        well_conn->add(createConnection(1, 1, 1, ijkToGlobal(1, 1, 1)));
        well_conn->add(createConnection(2, 2, 2, ijkToGlobal(2, 2, 2)));  // Intersects with Well 1
        well_conn->add(createConnection(3, 3, 3, ijkToGlobal(3, 3, 3)));
        auto well = createWell("INJECTOR1");
        well.updateConnections(well_conn, true);
        wells.push_back(well);
    }

    // Well 3 - Horizontal well that crosses both previous wells
    {
        auto well_conn = std::make_shared<Opm::WellConnections>();
        well_conn->add(createConnection(0, 2, 2, ijkToGlobal(0, 2, 2)));
        well_conn->add(createConnection(1, 2, 2, ijkToGlobal(1, 2, 2)));
        well_conn->add(createConnection(2, 2, 2, ijkToGlobal(2, 2, 2)));  // Intersects with Well 1 and 2
        well_conn->add(createConnection(3, 2, 2, ijkToGlobal(3, 2, 2)));  // Intersects with Well 1's branch
        well_conn->add(createConnection(4, 2, 2, ijkToGlobal(4, 2, 2)));
        auto well = createWell("PRODUCER2");
        well.updateConnections(well_conn, true);
        wells.push_back(well);
    }

    // Well 4 - L-shaped well connecting to the network
    {
        auto well_conn = std::make_shared<Opm::WellConnections>();
        well_conn->add(createConnection(2, 4, 1, ijkToGlobal(2, 4, 1)));  // Connects to Well 1's branch
        well_conn->add(createConnection(2, 4, 2, ijkToGlobal(2, 4, 2)));
        well_conn->add(createConnection(2, 4, 3, ijkToGlobal(2, 4, 3)));
        well_conn->add(createConnection(3, 4, 3, ijkToGlobal(3, 4, 3)));
        well_conn->add(createConnection(4, 4, 3, ijkToGlobal(4, 4, 3)));
        auto well = createWell("PRODUCER3");
        well.updateConnections(well_conn, true);
        wells.push_back(well);
    }

    // Well 5 - Isolated well (not connected to the network)
    {
        auto well_conn = std::make_shared<Opm::WellConnections>();
        well_conn->add(createConnection(0, 4, 0, ijkToGlobal(0, 4, 0)));
        well_conn->add(createConnection(1, 4, 0, ijkToGlobal(1, 4, 0)));
        auto well = createWell("ISOLATED");
        well.updateConnections(well_conn, true);
        wells.push_back(well);
    }

#if HAVE_MPI && HAVE_ZOLTAN
    // Test Zoltan partitioning with complex well network
    auto [part, num_part] = Opm::partitionCells("zoltan", 15, grid.leafGridView(),
                                               wells,
                                               std::unordered_map<std::string, std::set<int>>{},
                                               zoltan_ctrl,
                                               0);

    // Verify number of partitions and grid size
    BOOST_CHECK_EQUAL(num_part, 15);
    BOOST_CHECK_EQUAL(part.size(), 100);  // 5x5x4 grid

    // Helper to check if two cells are in the same partition
    auto inSamePartition = [&ppart = part](int idx1, int idx2) {
        return ppart[idx1] == ppart[idx2];
    };

    // Test 1: All cells in Well 1 (main bore and branches) should be in the same partition
    // Main wellbore
    BOOST_CHECK(inSamePartition(ijkToGlobal(2,2,0), ijkToGlobal(2,2,1)));
    BOOST_CHECK(inSamePartition(ijkToGlobal(2,2,1), ijkToGlobal(2,2,2)));
    BOOST_CHECK(inSamePartition(ijkToGlobal(2,2,2), ijkToGlobal(2,2,3)));
    // Branch 1
    BOOST_CHECK(inSamePartition(ijkToGlobal(2,2,2), ijkToGlobal(3,2,2)));
    BOOST_CHECK(inSamePartition(ijkToGlobal(3,2,2), ijkToGlobal(4,2,2)));
    // Branch 2
    BOOST_CHECK(inSamePartition(ijkToGlobal(2,2,1), ijkToGlobal(2,3,1)));
    BOOST_CHECK(inSamePartition(ijkToGlobal(2,3,1), ijkToGlobal(2,4,1)));

    // Test 2: Well 2 (diagonal) should be in same partition as Well 1 due to intersection
    BOOST_CHECK(inSamePartition(ijkToGlobal(0,0,0), ijkToGlobal(2,2,2)));
    BOOST_CHECK(inSamePartition(ijkToGlobal(3,3,3), ijkToGlobal(2,2,2)));

    // Test 3: Well 3 (horizontal) should be in same partition as Well 1 and 2
    BOOST_CHECK(inSamePartition(ijkToGlobal(0,2,2), ijkToGlobal(2,2,2)));
    BOOST_CHECK(inSamePartition(ijkToGlobal(4,2,2), ijkToGlobal(2,2,2)));

    // Test 4: Well 4 (L-shaped) should be in same partition as Well 1 due to connection
    BOOST_CHECK(inSamePartition(ijkToGlobal(2,4,1), ijkToGlobal(2,2,1)));
    BOOST_CHECK(inSamePartition(ijkToGlobal(4,4,3), ijkToGlobal(2,2,1)));

    // Test 5: Well 5 (isolated) should maintain its own cells in same partition
    BOOST_CHECK(inSamePartition(ijkToGlobal(0,4,0), ijkToGlobal(1,4,0)));
    // But should be in different partition from the connected network
    BOOST_CHECK(part[ijkToGlobal(0,4,0)] != part[ijkToGlobal(2,2,2)]);

    // Test 6: Verify the connected network (Wells 1-4) forms a single partition
    std::set<int> network_partition_ids;
    network_partition_ids.insert(part[ijkToGlobal(2,2,2)]);  // Well 1
    network_partition_ids.insert(part[ijkToGlobal(1,1,1)]);  // Well 2
    network_partition_ids.insert(part[ijkToGlobal(2,2,2)]);  // Well 3
    network_partition_ids.insert(part[ijkToGlobal(2,4,1)]);  // Well 4
    BOOST_CHECK_EQUAL(network_partition_ids.size(), 1);
#endif
}

BOOST_AUTO_TEST_CASE(PartitionCellsWithNonReachableCellsTest)
{
    // Create a 3x3x3 grid where:
    // - cell (0,0,0) is isolated (surrounded by inactive cells)
    // - cell (0,2,0) is isolated (surrounded by inactive cells)
    // - cells needed for wells are active
    // Grid structure (top view):
    //   Layer k=0:      Layer k=1:       Layer k=2:     W  = well cells
    //   [I][X][X ]      [X][  ][  ]      [X][  ][  ]    X  = inactive cells
    //   [X][X][W2]      [X][W1][W2]      [X][W1][W2]    I  = isolated active cell
    //   [I][X][  ]      [X][  ][  ]      [X][  ][  ]   [ ] = active cells
    const std::string deckString =
    R"(RUNSPEC

    DIMENS
    3 3 3 /

    GRID
    COORD
    0 0 0  0 0 1
    1 0 0  1 0 1
    2 0 0  2 0 1
    3 0 0  3 0 1
    0 1 0  0 1 1
    1 1 0  1 1 1
    2 1 0  2 1 1
    3 1 0  3 1 1
    0 2 0  0 2 1
    1 2 0  1 2 1
    2 2 0  2 2 1
    3 2 0  3 2 1
    0 3 0  0 3 1
    1 3 0  1 3 1
    2 3 0  2 3 1
    3 3 0  3 3 1
    /

    ZCORN
    36*0
    36*1
    36*1
    36*2
    36*2
    36*3
    /

    ACTNUM
    -- First layer (3x3)
    1 0 0  --  cell (0,0,0) is isolated
    0 0 1
    1 0 1  --  cell (1,1,1) is isolated
    -- Middle layer (3x3)
    0 1 1
    0 1 1
    0 1 1
    -- Top layer (3x3)
    0 1 1
    0 1 1
    0 1 1
    /

    END
    )";

    // Parse and create grid
    Opm::Parser parser;
    const auto deck = parser.parseString(deckString);
    Opm::EclipseGrid ecl_grid(deck);
    Dune::CpGrid grid;
    grid.processEclipseFormat(&ecl_grid, nullptr, false, false, false);

    // Verify grid properties
    const auto& grid_view = grid.leafGridView();
    const std::size_t num_active_cells = grid_view.size(0);
    BOOST_CHECK_EQUAL(num_active_cells, 16);  // 27 total cells - 11 inactive = 16 active

    // Helper functions
    using Entity = typename Dune::CpGrid::LeafGridView::template Codim<0>::Entity;
    auto zoltan_ctrl = createZoltanControl<Entity>(grid);
    auto ijkToGlobal = [](int i, int j, int k) { return i + (3 * j) + (9 * k); };

    auto wells = std::vector<Dune::cpgrid::OpmWellType>();
    // Well 1 - Middle column well
    {
        auto well_conn = std::make_shared<Opm::WellConnections>();
        well_conn->add(createConnection(1, 1, 1, ijkToGlobal(1, 1, 1)));
        well_conn->add(createConnection(1, 1, 2, ijkToGlobal(1, 1, 2)));
        auto well = createWell("WELL1");
        well.updateConnections(well_conn, true);
        wells.push_back(well);
    }

    // Well 2 - Right column well
    {
        auto well_conn = std::make_shared<Opm::WellConnections>();
        well_conn->add(createConnection(2, 1, 0, ijkToGlobal(2, 1, 0)));
        well_conn->add(createConnection(2, 1, 1, ijkToGlobal(2, 1, 1)));
        well_conn->add(createConnection(2, 1, 2, ijkToGlobal(2, 1, 2)));
        auto well = createWell("WELL2");
        well.updateConnections(well_conn, true);
        wells.push_back(well);
    }

#if HAVE_MPI && HAVE_ZOLTAN
    auto [part, num_part] = Opm::partitionCells("zoltan", 5, grid.leafGridView(),
                                               wells,
                                               std::unordered_map<std::string, std::set<int>>{},
                                               zoltan_ctrl,
                                               0);

    BOOST_CHECK_EQUAL(num_part, 5);
    BOOST_CHECK_EQUAL(part.size(), 16);

    // For this test global indices != local indices, so we need to
    // create a map from global index to local index
    std::map<int, int> g2l;
    for (const auto& element : elements(grid.leafGridView(), Dune::Partitions::interior)) {
        const auto globalIndex = grid.globalCell()[zoltan_ctrl.index(element)];
        g2l[globalIndex] = zoltan_ctrl.index(element);
    }

    // Verify that the two isolated cells are assigned to the -1 partition
    BOOST_CHECK_EQUAL(part[g2l[ijkToGlobal(0,0,0)]], -1);
    BOOST_CHECK_EQUAL(part[g2l[ijkToGlobal(0,2,0)]], -1);

    // Well 1: All cells should be in the same partition
    BOOST_CHECK(part[g2l[ijkToGlobal(1,1,1)]] == part[g2l[ijkToGlobal(1,1,2)]]);

    // Well 2: All cells should be in the same partition
    BOOST_CHECK(part[g2l[ijkToGlobal(2,1,0)]] == part[g2l[ijkToGlobal(2,1,1)]]);
    BOOST_CHECK(part[g2l[ijkToGlobal(2,1,1)]] == part[g2l[ijkToGlobal(2,1,2)]]);

    // Wells can be in different partitions as they don't connect
    BOOST_CHECK(part[g2l[ijkToGlobal(1,1,1)]] != part[g2l[ijkToGlobal(2,1,2)]]);
#endif
}

BOOST_AUTO_TEST_CASE(PartitionCellsWithNeighborConnectivityTest)
{
    // Create a 1D grid with 7 cells and two wells that are not directly connected
    // but will be connected through their neighboring cells when using neighbor_levels=1
    // Visual representation:
    // Cell indices:    0    1    2    3    4    5    6
    // Well layout:          W1        W2
    // Connected:       *    *    *    *    *
    // The wells will be connected because cell 2 is a neighbor of W1 and W2
    auto grid = createTestGrid({7, 1, 1}, {7.0, 1.0, 1.0});
    using Entity = typename Dune::CpGrid::LeafGridView::template Codim<0>::Entity;
    auto zoltan_ctrl = createZoltanControl<Entity>(grid);

    // Create two wells that are not directly connected
    auto wells = createWellsWithConnections({
        {"WELL1", {1}},    // Well 1 in cell 1
        {"WELL2", {3}}     // Well 2 in cell 3 (not adjacent to Well 1)
    });

#if HAVE_MPI && HAVE_ZOLTAN
    // Test Zoltan partitioning with neighbor connectivity
    auto [part, num_part] = Opm::partitionCells("zoltan", 3, grid.leafGridView(),
                                               wells,
                                               std::unordered_map<std::string, std::set<int>>{},
                                               zoltan_ctrl,
                                               1);

    // Verify number of partitions
    BOOST_CHECK_EQUAL(num_part, 3);
    BOOST_CHECK_EQUAL(part.size(), 7);  // 7x1x1 grid

    // Check that the wells are connected through their neighbors
    BOOST_CHECK(part[0] == part[1]);
    BOOST_CHECK(part[1] == part[2]);
    BOOST_CHECK(part[2] == part[3]);
    BOOST_CHECK(part[3] == part[4]);

    // Check that the remaining cells are in different partitions
    BOOST_CHECK(part[0] != part[5]);
    BOOST_CHECK(part[5] != part[6]);
#endif
}

bool
init_unit_test_func()
{
    return true;
}

int main(int argc, char** argv)
{
    Dune::MPIHelper::instance(argc, argv);
    boost::unit_test::unit_test_main(&init_unit_test_func,
                                     argc, argv);
}
