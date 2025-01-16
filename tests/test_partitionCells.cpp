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

#include <opm/simulators/flow/partitionCells.hpp>
#include <dune/grid/common/gridview.hh>
#include <opm/grid/CpGrid.hpp>

#include <opm/input/eclipse/Schedule/Well/Well.hpp>

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

BOOST_AUTO_TEST_CASE(PartitionCellsTestSimple)
{
    // Create a 3x2x2 grid (12 cells total)
    Dune::CpGrid grid;
    std::array<int, 3> dims{{3, 4, 1}};
    std::array<double, 3> size{{3.0, 4.0, 1.0}};
    grid.createCartesian(dims, size);

    // Setup Zoltan control parameters
    using Entity = typename Dune::CpGrid::LeafGridView::template Codim<0>::Entity;
    Opm::ZoltanPartitioningControl<Entity> zoltan_ctrl;
    zoltan_ctrl.domain_imbalance = 1.1;  // Allow 10% imbalance

    // Initialize index function to return entity index
    std::function<int(const Entity&)> index_func = [](const Entity& e) { return e.index(); };
    zoltan_ctrl.index = index_func;

    // Initialize local_to_global function to use index as global ID for this test
    std::function<int(int)> l2g_func = [](int idx) { return idx; };
    zoltan_ctrl.local_to_global = l2g_func;

    // Test simple partitioning first
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
}

BOOST_AUTO_TEST_CASE(PartitionCellsTestZoltan)
{
    // Create a 3x2x2 grid (12 cells total)
    Dune::CpGrid grid;
    std::array<int, 3> dims{{3, 4, 1}};
    std::array<double, 3> size{{3.0, 4.0, 1.0}};
    grid.createCartesian(dims, size);

    // Setup Zoltan control parameters
    using Entity = typename Dune::CpGrid::LeafGridView::template Codim<0>::Entity;
    Opm::ZoltanPartitioningControl<Entity> zoltan_ctrl;
    zoltan_ctrl.domain_imbalance = 1.1;  // Allow 10% imbalance

    // Initialize index function to return entity index
    std::function<int(const Entity&)> index_func = [](const Entity& e) { return e.index(); };
    zoltan_ctrl.index = index_func;

    // Initialize local_to_global function to use index as global ID for this test
    std::function<int(int)> l2g_func = [](int idx) { return idx; };
    zoltan_ctrl.local_to_global = l2g_func;

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
