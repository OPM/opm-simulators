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

#include <opm/simulators/flow/partitionCells.hpp>

// Helper functions
namespace {
    Opm::Connection createConnection(int i, int j, int k) {
        // Calculate global index from i,j,k coordinates
        // For this 1D test grid, global_index is just i since j=k=0
        std::size_t global_index = i;
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
                well_conn->add(createConnection(idx, 0, 0));
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
                                                   zoltan_ctrl);
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
                                                   zoltan_ctrl);
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
                                               zoltan_ctrl);

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
                                               zoltan_ctrl);

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
