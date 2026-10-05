// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
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

  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the list of
  copyright holders.
*/
#include "config.h"

#define BOOST_TEST_MODULE ParallelTransmissibility

#include <opm/models/utils/propertysystem.hh>

#include <opm/simulators/flow/FlowGenericVanguard.hpp>
#include <opm/simulators/flow/Main.hpp>
#include <opm/simulators/flow/TTagFlowProblemTPFA.hpp>

#include <opm/grid/common/CommunicationUtils.hpp>

#if HAVE_DUNE_FEM
#include <dune/fem/misc/mpimanager.hh>
#else
#include <dune/common/parallel/mpihelper.hh>
#endif

#include <dune/grid/common/partitionset.hh>
#include <dune/grid/common/rangegenerators.hh>

#include <boost/test/unit_test.hpp>
#include <boost/version.hpp>
#if BOOST_VERSION / 100000 == 1 && BOOST_VERSION / 100 % 1000 < 71
#include <boost/test/floating_point_comparison.hpp>
#else
#include <boost/test/tools/floating_point_comparison.hpp>
#endif

#include <algorithm>
#include <cstddef>
#include <cstdlib>
#include <map>
#include <memory>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace Opm::Properties::TTag {

struct TestTransmissibilityTypeTag
{ using InheritsFrom = std::tuple<FlowProblemTPFA>; };

}

namespace Opm::Properties {

template<class TypeTag>
struct EnableConvectiveMixing<TypeTag, TTag::TestTransmissibilityTypeTag>
{ static constexpr bool value = false; };

template<class TypeTag>
struct EnableDiffusion<TypeTag, TTag::TestTransmissibilityTypeTag>
{ static constexpr bool value = false; };

}

namespace {

using TypeTag = Opm::Properties::TTag::TestTransmissibilityTypeTag;
using Simulator = Opm::GetPropType<TypeTag, Opm::Properties::Simulator>;

// Runs Flow's set-up, which distributes the grid and computes the
// transmissibilities, but no time steps.
class InitOnlyMain : public Opm::Main
{
public:
    InitOnlyMain(int argc, char** argv)
        // MPI is initialised once for the whole test process.
        : Main(argc, argv, /*ownMPI=*/false)
    {
        auto exitCode = EXIT_SUCCESS;
        if (this->initialize_<Opm::Properties::TTag::FlowEarlyBird>(exitCode, /*keep_keywords=*/false)) {
            this->setupVanguard();
            flowMain_ = std::make_unique<Opm::FlowMain<TypeTag>>(this->argc_, this->argv_,
                                                                  this->outputCout_, this->outputFiles_);
            if (flowMain_->executeInitStep() == EXIT_SUCCESS) {
                simulator_ = flowMain_->getSimulatorPtr();
            }
        }
    }

    const Simulator* simulator() const
    { return simulator_; }

private:
    std::unique_ptr<Opm::FlowMain<TypeTag>> flowMain_;
    Simulator* simulator_ = nullptr;
};

class FlowSetup
{
public:
    explicit FlowSetup(std::vector<std::string> args)
        : args_(std::move(args))
    {
        for (auto& arg : args_) {
            argv_.push_back(arg.data());
        }
        // Main expects argv[argc] to be a null pointer.
        argv_.push_back(nullptr);
        main_ = std::make_unique<InitOnlyMain>(static_cast<int>(args_.size()), argv_.data());
    }

    const Simulator& simulator() const
    {
        BOOST_REQUIRE(main_->simulator() != nullptr);
        return *main_->simulator();
    }

private:
    std::vector<std::string> args_;
    std::vector<char*> argv_;
    std::unique_ptr<InitOnlyMain> main_;
};

struct GlobalFixture
{
    GlobalFixture()
    {
        auto argc = boost::unit_test::framework::master_test_suite().argc;
        auto argv = boost::unit_test::framework::master_test_suite().argv;
#if HAVE_DUNE_FEM
        Dune::Fem::MPIManager::initialize(argc, argv);
#else
        Dune::MPIHelper::instance(argc, argv);
#endif
        Opm::FlowGenericVanguard::setCommunication(std::make_unique<Opm::Parallel::Communication>());
    }
};

} // Anonymous namespace

BOOST_GLOBAL_FIXTURE(GlobalFixture);

BOOST_AUTO_TEST_CASE(OwnersAgreeOnProcessBoundaryFaces)
{
    const auto numRanks = Opm::FlowGenericVanguard::comm().size();
    const FlowSetup flow {{"test_parallel_transmissibility",
                           "--edge-conformal=true",
                           "--output-dir=parallel_transmissibility_np" + std::to_string(numRanks),
                           "parallel_transmissibility.DATA"}};

    const auto& simulator = flow.simulator();
    const auto& vanguard = simulator.vanguard();
    const auto& gridView = vanguard.gridView();
    const auto& elemMapper = simulator.model().elementMapper();

    // Every face between an owned cell and a cell owned by another rank, seen
    // from the owned cell: Cartesian indices of both cells, transmissibility.
    std::vector<double> faces;
    for (const auto& elem : elements(gridView, Dune::Partitions::interior)) {
        const unsigned inside = elemMapper.index(elem);
        for (const auto& intersection : intersections(gridView, elem)) {
            if (!intersection.neighbor() ||
                intersection.outside().partitionType() == Dune::InteriorEntity)
            {
                continue;
            }
            const unsigned outside = elemMapper.index(intersection.outside());
            faces.insert(faces.end(), {static_cast<double>(vanguard.cartesianIndex(inside)),
                                       static_cast<double>(vanguard.cartesianIndex(outside)),
                                       simulator.problem().transmissibility(inside, outside)});
        }
    }

    const auto [allFaces, offsets] = Opm::allGatherv(faces, gridView.comm());

    std::map<std::pair<int, int>, std::vector<double>> transByFace;
    for (std::size_t i = 0; i < allFaces.size(); i += 3) {
        const auto cell1 = static_cast<int>(allFaces[i]);
        const auto cell2 = static_cast<int>(allFaces[i + 1]);
        transByFace[std::minmax(cell1, cell2)].push_back(allFaces[i + 2]);
    }

    // All active cells are merged ones, so any face between ranks exercises them.
    if (numRanks > 1) {
        BOOST_REQUIRE(!transByFace.empty());
    }

    for (const auto& [face, trans] : transByFace) {
        BOOST_TEST_CONTEXT("Cells " << face.first << " and " << face.second) {
            BOOST_REQUIRE_EQUAL(trans.size(), 2u);
            BOOST_CHECK_CLOSE(trans[0], trans[1], 1.0e-10);
        }
    }
}
