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

#define BOOST_TEST_MODULE EditDepth

#include <opm/models/utils/propertysystem.hh>

#include <opm/simulators/flow/FlowGenericVanguard.hpp>
#include <opm/simulators/flow/Main.hpp>
#include <opm/simulators/flow/TTagFlowProblemTPFA.hpp>

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

#include <cstdlib>
#include <map>
#include <memory>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace Opm::Properties::TTag {

struct TestEditDepthTypeTag
{ using InheritsFrom = std::tuple<FlowProblemTPFA>; };

}

namespace Opm::Properties {

template<class TypeTag>
struct EnableConvectiveMixing<TypeTag, TTag::TestEditDepthTypeTag>
{ static constexpr bool value = false; };

template<class TypeTag>
struct EnableDiffusion<TypeTag, TTag::TestEditDepthTypeTag>
{ static constexpr bool value = false; };

}

namespace {

using TypeTag = Opm::Properties::TTag::TestEditDepthTypeTag;
using Simulator = Opm::GetPropType<TypeTag, Opm::Properties::Simulator>;

// Runs Flow's set-up, which distributes the grid and computes the cell
// depths, but no time steps.
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

BOOST_AUTO_TEST_CASE(OnlyEditedCellsTakeTheirDepth)
{
    const auto numRanks = Opm::FlowGenericVanguard::comm().size();
    const FlowSetup flow {{"test_edit_depth",
                           "--edge-conformal=true",
                           "--output-dir=edit_depth_np" + std::to_string(numRanks),
                           "edit_depth.DATA"}};

    const auto& simulator = flow.simulator();
    const auto& vanguard = simulator.vanguard();
    const auto& gridView = vanguard.gridView();
    const auto& elemMapper = simulator.model().elementMapper();

    // Depth of each active cell by Cartesian index.  Cells 4 and 5 have
    // absorbed the cells above them and span 1000-1002 m, while the input
    // grid has them between 1001 and 1002 m.
    const std::map<unsigned, double> expected {
        {2, 1005.0}, // edited
        {3, 1000.5},
        {4, 1001.0}, {5, 1001.0},
        {6, 1001.5}, {7, 1001.5},
    };

    int numChecked = 0;
    for (const auto& elem : elements(gridView, Dune::Partitions::interior)) {
        const unsigned elemIdx = elemMapper.index(elem);
        const auto cartIdx = vanguard.cartesianIndex(elemIdx);
        BOOST_TEST_CONTEXT("Cell " << cartIdx) {
            BOOST_REQUIRE(expected.count(cartIdx) == 1);
            BOOST_CHECK_CLOSE(vanguard.cellCenterDepth(elemIdx), expected.at(cartIdx), 1.0e-10);
            ++numChecked;
        }
    }

    BOOST_CHECK_EQUAL(gridView.comm().sum(numChecked), static_cast<int>(expected.size()));
}
