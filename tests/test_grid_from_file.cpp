// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:

#include "config.h"

#define BOOST_TEST_MODULE GridFromFileSimulator

#include <boost/test/unit_test.hpp>

#include "SimulatorFixture.hpp"

#include <opm/grid/UnstructuredGrid.h>
#include <opm/grid/CpGrid.hpp>
#include <opm/grid/cart_grid.h>
#include <opm/grid/cpgrid/GridHelpers.hpp>
#include <opm/grid/polyhedralgrid.hh>

#include <opm/simulators/flow/FlowProblemBlackoil.hpp>

#include <opm/models/blackoil/blackoillocalresidualtpfa.hh>
#include <opm/models/discretization/common/tpfalinearizer.hh>

#include <opm/simulators/aquifers/SupportsFaceTag.hpp>
#include <opm/simulators/flow/PolyhedralGridVanguard.hpp>

// These templates are not explicitly instantiated in the library for polyhedral grids.
#include <opm/simulators/flow/CollectDataOnIORank_impl.hpp>
#include <opm/simulators/flow/EclGenericWriter_impl.hpp>
#include <opm/simulators/flow/FlowGenericProblem_impl.hpp>
#include <opm/simulators/flow/GenericThresholdPressure_impl.hpp>
#include <opm/simulators/flow/GenericTracerModel_impl.hpp>
#include <opm/simulators/flow/Transmissibility_impl.hpp>
#include <opm/simulators/flow/equil/InitStateEquil_impl.hpp>
#include <opm/simulators/utils/GridDataOutput_impl.hpp>

#include <cstdio>
#include <fstream>
#include <memory>
#include <type_traits>

namespace Opm::Properties {

namespace TTag {

struct TestPolyhedralTypeTag {
    using InheritsFrom = std::tuple<TestTypeTag>;
};

} // namespace TTag

template<class TypeTag>
struct Linearizer<TypeTag, TTag::TestPolyhedralTypeTag>
{
    using type = TpfaLinearizer<TypeTag>;
};

template<class TypeTag>
struct LocalResidual<TypeTag, TTag::TestPolyhedralTypeTag>
{
    using type = BlackOilLocalResidualTPFA<TypeTag>;
};

template<class TypeTag>
struct EnableDiffusion<TypeTag, TTag::TestPolyhedralTypeTag>
{
    static constexpr bool value = false;
};

template<class TypeTag>
struct Grid<TypeTag, TTag::TestPolyhedralTypeTag>
{
    using type = Dune::PolyhedralGrid<3, 3>;
};

template<class TypeTag>
struct EquilGrid<TypeTag, TTag::TestPolyhedralTypeTag>
{
    using type = GetPropType<TypeTag, Properties::Grid>;
};

template<class TypeTag>
struct Vanguard<TypeTag, TTag::TestPolyhedralTypeTag>
{
    using type = Opm::PolyhedralGridVanguard<TypeTag>;
};

} // namespace Opm::Properties

namespace Opm {

template<>
class SupportsFaceTag<Dune::PolyhedralGrid<3, 3>>
    : public std::bool_constant<true>
{};

} // namespace Opm

namespace {

class ScopedFile
{
public:
    explicit ScopedFile(std::string path)
        : path_(std::move(path))
    {}

    ~ScopedFile()
    {
        std::remove(path_.c_str());
    }

    const std::string& path() const
    {
        return path_;
    }

private:
    std::string path_;
};

template <class T>
void writeValues(std::ostream& os, const T* values, std::size_t count)
{
    for (std::size_t i = 0; i < count; ++i) {
        if (i > 0) {
            os << ' ';
        }
        os << values[i];
    }
    os << '\n';
}

void writeTrivialDeckFile(const std::string& deckFileName, int numCells)
{
    std::ofstream deck(deckFileName);
    BOOST_REQUIRE(deck);

    deck << "RUNSPEC\n\n"
         << "WATER\n"
         << "GAS\n"
         << "OIL\n\n"
         << "METRIC\n\n"
         << "DIMENS\n"
         << "    " << numCells << " 1 1 /\n\n"
         << "GRID\n\n"
         << "DX\n"
         << "    " << numCells << "*1 /\n"
         << "DY\n"
         << "    " << numCells << "*1 /\n"
         << "DZ\n"
         << "    " << numCells << "*1 /\n\n"
         << "TOPS\n"
         << "    " << numCells << "*0 /\n\n"
         << "PORO\n"
         << "    " << numCells << "*0.3 /\n\n"
         << "PERMX\n"
         << "    " << numCells << "*500 /\n\n"
         << "PROPS\n\n"
         << "PVTW\n"
         << "    4017.55 1.038 3.22E-6 0.318 0.0 /\n\n"
         << "ROCK\n"
         << "    14.7 3E-6 /\n\n"
         << "SWOF\n"
         << "0.12 0 1 0\n"
         << "0.30 0.1 0.9 0\n"
         << "1.00 1 0 0 /\n\n"
         << "SGOF\n"
         << "0.00 0 1 0\n"
         << "0.30 0.1 0.9 0\n"
         << "0.88 1 0 0 /\n\n"
         << "DENSITY\n"
         << "    53.66 64.49 0.0533 /\n\n"
         << "PVDG\n"
         << "14.700 166.666 0.008000\n"
         << "9014.7 0.38600 0.047000 /\n\n"
         << "PVTO\n"
         << "0.0010 14.7 1.0620 1.0400 /\n"
         << "1.2700 4014.7 1.6950 0.5100\n"
         << "    9014.7 1.5790 0.7400 /\n"
         << "/\n\n"
         << "SOLUTION\n\n"
         << "SWAT\n"
         << "    " << numCells << "*0.12 /\n\n"
         << "SGAS\n"
         << "    " << numCells << "*0 /\n\n"
         << "PRESSURE\n"
         << "    " << numCells << "*300 /\n\n"
         << "SUMMARY\n\n"
         << "SCHEDULE\n\n"
         << "TSTEP\n"
         << "1 /\n";
}

void writeCartesianUnstructuredGridFile(const std::string& gridFileName, int numCells)
{
    std::unique_ptr<UnstructuredGrid, decltype(&destroy_grid)> ugPtr(create_grid_cart3d(numCells, 1, 1), &destroy_grid);
    BOOST_REQUIRE(ugPtr);

    const auto& grid = *ugPtr;

    std::ofstream gridFile(gridFileName);
    BOOST_REQUIRE(gridFile);

    const auto numFaceNodes = static_cast<std::size_t>(grid.face_nodepos[grid.number_of_faces]);
    const auto numCellFaces = static_cast<std::size_t>(grid.cell_facepos[grid.number_of_cells]);
    const bool hasFaceTags = grid.cell_facetag != nullptr;
    const bool hasIndexMap = grid.global_cell != nullptr;

    gridFile << grid.dimensions << ' '
             << grid.number_of_cells << ' '
             << grid.number_of_faces << ' '
             << grid.number_of_nodes << ' '
             << numFaceNodes << ' '
             << numCellFaces << ' '
             << hasFaceTags << ' ' << hasIndexMap << '\n';

    writeValues(gridFile, grid.cartdims, static_cast<std::size_t>(grid.dimensions));
    writeValues(gridFile, grid.node_coordinates, static_cast<std::size_t>(grid.dimensions * grid.number_of_nodes));
    writeValues(gridFile, grid.face_nodepos, static_cast<std::size_t>(grid.number_of_faces + 1));
    writeValues(gridFile, grid.face_nodes, numFaceNodes);
    writeValues(gridFile, grid.face_cells, static_cast<std::size_t>(2 * grid.number_of_faces));
    writeValues(gridFile, grid.face_areas, static_cast<std::size_t>(grid.number_of_faces));
    writeValues(gridFile, grid.face_centroids, static_cast<std::size_t>(grid.dimensions * grid.number_of_faces));
    writeValues(gridFile, grid.face_normals, static_cast<std::size_t>(grid.dimensions * grid.number_of_faces));
    writeValues(gridFile, grid.cell_facepos, static_cast<std::size_t>(grid.number_of_cells + 1));

    if (hasFaceTags) {
        for (std::size_t i = 0; i < numCellFaces; ++i) {
            if (i > 0) {
                gridFile << ' ';
            }
            gridFile << grid.cell_faces[i] << ' ' << grid.cell_facetag[i];
        }
        gridFile << '\n';
    }
    else {
        writeValues(gridFile, grid.cell_faces, numCellFaces);
    }

    if (hasIndexMap) {
        writeValues(gridFile, grid.global_cell, static_cast<std::size_t>(grid.number_of_cells));
    }

    writeValues(gridFile, grid.cell_volumes, static_cast<std::size_t>(grid.number_of_cells));
    writeValues(gridFile, grid.cell_centroids, static_cast<std::size_t>(grid.dimensions * grid.number_of_cells));
}

} // namespace

using SimulatorFixture = Opm::SimulatorFixture;
BOOST_GLOBAL_FIXTURE(SimulatorFixture);

BOOST_AUTO_TEST_CASE(init_polyhedral_simulator_from_unstructured_grid_file)
{
    using TypeTag = Opm::Properties::TTag::TestPolyhedralTypeTag;

    constexpr int numCells = 3;
    const ScopedFile deckFile{"test_polyhedral_unstructured.DATA"};
    const ScopedFile gridFile{"test_polyhedral_unstructured.grid"};

    writeTrivialDeckFile(deckFile.path(), numCells);
    writeCartesianUnstructuredGridFile(gridFile.path(), numCells);

    auto simulator = Opm::initSimulator<TypeTag>(deckFile.path().c_str(),
                                                 "test_polyhedral_unstructured",
                                                 /*threads_per_process=*/1,
                                                 {std::string{"--unstructured-grid-file-name="} + gridFile.path()});
    BOOST_REQUIRE(simulator);

    const auto& grid = simulator->vanguard().grid();
    const auto& eclGrid = simulator->vanguard().eclState().getInputGrid();

    BOOST_CHECK_EQUAL(eclGrid.getNumActive(), numCells);
    BOOST_CHECK_EQUAL(static_cast<std::size_t>(grid.size(0)), eclGrid.getNumActive());
    BOOST_CHECK_EQUAL(simulator->vanguard().globalCell().size(), grid.size(0));
}

BOOST_AUTO_TEST_CASE(init_cpgrid_simulator_from_unstructured_grid_file)
{
    using TypeTag = Opm::Properties::TTag::TestTypeTag;

    constexpr int numCells = 3;
    const ScopedFile deckFile{"test_cpgrid_unstructured.DATA"};
    const ScopedFile gridFile{"test_cpgrid_unstructured.grid"};

    writeTrivialDeckFile(deckFile.path(), numCells);
    writeCartesianUnstructuredGridFile(gridFile.path(), numCells);

    auto simulator = Opm::initSimulator<TypeTag>(deckFile.path().c_str(),
                                                 "test_cpgrid_unstructured",
                                                 /*threads_per_process=*/1,
                                                 {std::string{"--unstructured-grid-file-name="} + gridFile.path()});
    BOOST_REQUIRE(simulator);

    const auto& grid = simulator->vanguard().grid();
    const auto& gridView = simulator->vanguard().gridView();
    const auto& eclGrid = simulator->vanguard().eclState().getInputGrid();

    BOOST_CHECK_EQUAL(eclGrid.getNumActive(), numCells);
    BOOST_CHECK_EQUAL(grid.numCells(), numCells);
    BOOST_CHECK_EQUAL(static_cast<std::size_t>(grid.numCells()), gridView.size(0));
    BOOST_CHECK_EQUAL(grid.globalCell().size(), static_cast<std::size_t>(grid.numCells()));
}

BOOST_AUTO_TEST_CASE(init_polyhedral_simulator_from_deck)
{
    using TypeTag = Opm::Properties::TTag::TestPolyhedralTypeTag;
    using Grid = Opm::GetPropType<TypeTag, Opm::Properties::Grid>;
    using Vanguard = Opm::GetPropType<TypeTag, Opm::Properties::Vanguard>;

    static_assert(std::is_same_v<Grid, Dune::PolyhedralGrid<3, 3>>);
    static_assert(std::is_same_v<Vanguard, Opm::PolyhedralGridVanguard<TypeTag>>);

    auto simulator = Opm::initSimulator<TypeTag>("equil_base.DATA", "test_polyhedral");
    BOOST_REQUIRE(simulator);

    const auto& grid = simulator->vanguard().grid();
    const auto& gridView = simulator->vanguard().gridView();
    const auto& eclGrid = simulator->vanguard().eclState().getInputGrid();

    BOOST_CHECK_EQUAL(static_cast<std::size_t>(grid.size(0)), eclGrid.getNumActive());
    BOOST_CHECK_GT(gridView.size(0), 0u);
    BOOST_CHECK_GT(gridView.size(1), 0u);
    BOOST_CHECK_GT(gridView.size(3), 0u);
    BOOST_CHECK_EQUAL(simulator->vanguard().globalCell().size(), gridView.size(0));
}

BOOST_AUTO_TEST_CASE(init_cpgrid_simulator_from_deck)
{
    using TypeTag = Opm::Properties::TTag::TestTypeTag;

    auto simulator = Opm::initSimulator<TypeTag>("equil_base.DATA", "test_cpgrid");
    BOOST_REQUIRE(simulator);

    const auto& grid = simulator->vanguard().grid();
    const auto& gridView = simulator->vanguard().gridView();
    const auto& eclGrid = simulator->vanguard().eclState().getInputGrid();

    BOOST_CHECK_EQUAL(static_cast<std::size_t>(grid.numCells()), eclGrid.getNumActive());
    BOOST_CHECK_GT(gridView.size(0), 0u);
    BOOST_CHECK_EQUAL(grid.globalCell().size(), gridView.size(0));
}
