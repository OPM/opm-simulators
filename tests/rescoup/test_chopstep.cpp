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

#define BOOST_TEST_MODULE ResCoup_ChopStep
#include <opm/models/utils/propertysystem.hh>
#include <opm/models/utils/parametersystem.hpp>
#include <opm/models/utils/start.hh>

#include <opm/input/eclipse/EclipseState/SummaryConfig/SummaryConfig.hpp>
#include <opm/input/eclipse/Schedule/Schedule.hpp>
#include <opm/input/eclipse/Schedule/ResCoup/ReservoirCouplingInfo.hpp>
#include <opm/input/eclipse/Schedule/ResCoup/MasterGroup.hpp>
#include <opm/input/eclipse/Schedule/ResCoup/Slaves.hpp>
#include <opm/simulators/flow/FlowGenericVanguard.hpp>
#include <opm/simulators/flow/Main.hpp>
#include <opm/simulators/flow/TTagFlowProblemTPFA.hpp>
#include <opm/simulators/flow/ReservoirCoupling.hpp>
#include <opm/simulators/flow/ReservoirCouplingMaster.hpp>
#include <opm/simulators/flow/ReservoirCouplingSlave.hpp>
#include <opm/simulators/utils/DeferredLogger.hpp>
#include <opm/simulators/flow/BlackoilModel.hpp>
#include <opm/simulators/flow/FlowProblemBlackoil.hpp>

#if HAVE_DUNE_FEM
#include <dune/fem/misc/mpimanager.hh>
#else
#include <dune/common/parallel/mpihelper.hh>
#endif

#include <boost/test/unit_test.hpp>
#include <boost/version.hpp>
#if BOOST_VERSION / 100000 == 1 && BOOST_VERSION / 100 % 1000 < 71
#include <boost/test/floating_point_comparison.hpp>
#else
#include <boost/test/tools/floating_point_comparison.hpp>
#endif

namespace Opm {

class MainTestWrapper : public Main
{
public:
    using TypeTag = Properties::TTag::FlowProblemTPFA;
    using FlowMainPtr = std::unique_ptr<FlowMain<TypeTag>>;
    using Simulator = GetPropType<TypeTag, Properties::Simulator>;

    MainTestWrapper(int argc, char** argv) :
    // NOTE: passing ownMPI=false to prevent each test from initializing MPI (which is illegal)
          Main{argc, argv, /*ownMPI=*/false}
        , exit_code_{EXIT_SUCCESS}
        , flow_main_{nullptr}
        , simulator_{nullptr}
    {
        this->exit_code_ = EXIT_SUCCESS;
        if (initialize_<Properties::TTag::FlowEarlyBird>(this->exit_code_, /*keep_keywords=*/false)) {
            this->setupVanguard();
            this->flow_main_ = std::make_unique<FlowMain<TypeTag>>(
                this->argc_, this->argv_, this->outputCout_, this->outputFiles_
            );
            this->exit_code_ = this->flow_main_->executeInitStep();
            if (this->exit_code_ == EXIT_SUCCESS) {
                this->simulator_ = this->flow_main_->getSimulatorPtr();
            }
        }
    }
    Simulator* getSimulator() { return simulator_; }

private:
    int exit_code_;
    FlowMainPtr flow_main_;
    Simulator *simulator_;
};

} // namespace Opm

namespace {
class OpmSimulatorTestCase
{
public:
    using MainTestWrapper = Opm::MainTestWrapper;
    using Simulator = Opm::MainTestWrapper::Simulator;

    OpmSimulatorTestCase(
        const std::string& prog_name,
        const std::string& deck_filename,
        const std::vector<std::string>& args
    )
    {
        this->save_argv_.push_back(prog_name);
        this->save_argv_.insert(save_argv_.end(), args.begin(), args.end());
        this->save_argv_.push_back(deck_filename);
        this->argv_.reserve(this->save_argv_.size() + 1);  // +1 for nullptr
        for (const auto& arg : this->save_argv_) {
            this->argv_.push_back(const_cast<char*>(arg.c_str()));
        }
        // Make sure the last element is nullptr, this is required by MPI_Init()
        this->argv_.push_back(nullptr);
        char **argv = this->argv_.data();
        int argc = this->argv_.size() - 1;
        this->main_ = std::make_unique<MainTestWrapper>(argc, argv);
    }

    MainTestWrapper& getMain()
    {
        return *main_;
    }

    Simulator* getSimulatorPtr()
    {
        return main_->getSimulator();
    }

    private:
        std::unique_ptr<MainTestWrapper> main_;
        std::vector<std::string> save_argv_;
        std::vector<char *> argv_;
};

struct SimulatorFixture
{
    using Simulator = Opm::MainTestWrapper::Simulator;
    SimulatorFixture() :
        simulator_wrapper{
            /*prog_name=*/"flow_test_chopstep",
            /*deck_filename=*/"RC-01_MAST_PRED.DATA",
            /*args=*/{"--parsing-strictness=low"}  // Need this until PR #5643 is merged
        },
        simulator{simulator_wrapper.getSimulatorPtr()},
        schedule{simulator->vanguard().schedule()},
        rc_master{Opm::FlowGenericVanguard::comm(), schedule, 0, nullptr},
        start_date{static_cast<double>(schedule.getStartTime())}
    {
        rc_master.addSlaveName("RES-1");
        rc_master.addSlaveName("RES-2");
        rc_master.resizeSlaveStartDates(2);
        rc_master.resizeNextReportDates(2);
    }

    void setSlaveStartDateAndNextReportTimeOffset(int slave_number, double slave_start_date, double report_time_step_size)
    {
        rc_master.setSlaveStartDate(slave_number, slave_start_date);
        rc_master.setSlaveNextReportTimeOffset(slave_number, report_time_step_size);
    }
    OpmSimulatorTestCase simulator_wrapper;
    Simulator *simulator;
    Opm::Schedule& schedule;
    Opm::ReservoirCouplingMaster rc_master;
    double start_date;
};

struct GlobalTestFixture
{
    // NOTE: We can only call MPI_Init() and MPI_Finalize() once per process, so we prevent
    // Opm::Main() from initializing MPI by setting ownMPI=false in the constructor, and instead do
    // the initialization once for the whole test case process here.
    GlobalTestFixture()
    {
        int argc = boost::unit_test::framework::master_test_suite().argc;
        char** argv = boost::unit_test::framework::master_test_suite().argv;
    #if HAVE_DUNE_FEM
        Dune::Fem::MPIManager::initialize(argc, argv);
    #else
        Dune::MPIHelper::instance(argc, argv);
#endif
        Opm::FlowGenericVanguard::setCommunication(std::make_unique<Opm::Parallel::Communication>());
    }
};

} // Anonymous namespace


BOOST_GLOBAL_FIXTURE(GlobalTestFixture);

BOOST_FIXTURE_TEST_SUITE(ResCoupTestSuite, SimulatorFixture)

BOOST_AUTO_TEST_CASE(NoChop)
{
    BOOST_CHECK_EQUAL(start_date, 1538352000);  // 1538352000 = 2018-10-01 00:00:00
    double time_step = 60*60*24;  // 1 day
    double elapsed_time = 0.0;
    auto tol = Opm::ReservoirCoupling::Seconds::reltol / 2;

    // Check that the time step is not chopped when slave processes have identical start dates
    setSlaveStartDateAndNextReportTimeOffset(/*slave_number=*/0, /*start_date=*/start_date, /*report_time_step_size=*/time_step);
    setSlaveStartDateAndNextReportTimeOffset(/*slave_number=*/1, /*start_date=*/start_date, /*report_time_step_size=*/time_step);
    auto new_time_step = rc_master.maybeChopSubStep(time_step, elapsed_time);
    BOOST_CHECK_CLOSE(new_time_step, time_step, 1e-16);

    // Check that the time step is not chopped when slave processes start at the same time as
    // the master but their report steps end after the master process
    setSlaveStartDateAndNextReportTimeOffset(/*slave_number=*/0, /*start_date=*/start_date, /*report_time_step_size=*/time_step+1);
    setSlaveStartDateAndNextReportTimeOffset(/*slave_number=*/1, /*start_date=*/start_date, /*report_time_step_size=*/time_step+2);
    new_time_step = rc_master.maybeChopSubStep(time_step, elapsed_time);
    BOOST_CHECK_CLOSE(new_time_step, time_step, 1e-16);

    // Check that the time step is not chopped when slave processes start at the same time as
    // the master but their report steps end before the master process but within the tolerance
    setSlaveStartDateAndNextReportTimeOffset(/*slave_number=*/0, /*start_date=*/start_date, /*report_time_step_size=*/time_step-tol);
    setSlaveStartDateAndNextReportTimeOffset(/*slave_number=*/1, /*start_date=*/start_date, /*report_time_step_size=*/time_step+2);
    new_time_step = rc_master.maybeChopSubStep(time_step, elapsed_time);
    BOOST_CHECK_CLOSE(new_time_step, time_step, 1e-16);

    // Check that the time step is not chopped when slave processes start before
    // the master but their report steps end after the master process
    setSlaveStartDateAndNextReportTimeOffset(/*slave_number=*/0, /*start_date=*/start_date-1, /*report_time_step_size=*/time_step+1+tol);
    setSlaveStartDateAndNextReportTimeOffset(/*slave_number=*/1, /*start_date=*/start_date-2, /*report_time_step_size=*/time_step+2+tol);
    new_time_step = rc_master.maybeChopSubStep(time_step, elapsed_time);
    BOOST_CHECK_CLOSE(new_time_step, time_step, 1e-16);

    // Check that the time step is not chopped when slave processes starts after the master's
    // report step has ended
    setSlaveStartDateAndNextReportTimeOffset(/*slave_number=*/0,
                /*start_date=*/start_date+time_step+1,
                /*report_time_step_size=*/time_step);
    setSlaveStartDateAndNextReportTimeOffset(/*slave_number=*/1,
                 /*start_date=*/start_date+time_step+tol,
                 /*report_time_step_size=*/time_step);
    new_time_step = rc_master.maybeChopSubStep(time_step, elapsed_time);
    BOOST_CHECK_CLOSE(new_time_step, time_step, 1e-16);
}

BOOST_AUTO_TEST_CASE(Chop)
{
    BOOST_CHECK_EQUAL(start_date, 1538352000);  // 1538352000 = 2018-10-01 00:00:00
    double time_step = 60*60*24;  // 1 day
    double elapsed_time = 0.0;
    double tol;
    double new_time_step;

    // Check that the time step is chopped when slave processes start in the middle of the master's
    // report step
    setSlaveStartDateAndNextReportTimeOffset(/*slave_number=*/0,
                /*start_date=*/start_date+time_step/2,
                /*report_time_step_size=*/time_step);
    setSlaveStartDateAndNextReportTimeOffset(/*slave_number=*/1,
                 /*start_date=*/start_date+time_step/2,
                 /*report_time_step_size=*/time_step);
    new_time_step = rc_master.maybeChopSubStep(time_step, elapsed_time);
    BOOST_CHECK_CLOSE(new_time_step, time_step/2, 1e-16);

    // Check that the time step is chopped when slave processes start at the beginning of the master's
    // report step, but ends within but just outside the tolerance of the master's report step
    // NOTE: microseconds are approximately the smallest time units that can be represented accurately
    // for epoch values in this century (2000-2100), so first check this
    tol = 1e-8;  // This is then too small, and should not cause the time step to be chopped
    setSlaveStartDateAndNextReportTimeOffset(/*slave_number=*/0,
                /*start_date=*/start_date,
                /*report_time_step_size=*/time_step - tol);
    setSlaveStartDateAndNextReportTimeOffset(/*slave_number=*/1,
                 /*start_date=*/start_date+time_step,
                 /*report_time_step_size=*/time_step);
    new_time_step = rc_master.maybeChopSubStep(time_step, elapsed_time);
    BOOST_CHECK(Opm::ReservoirCoupling::Seconds::compare_eq(new_time_step, time_step));
    tol = 1e-5;  // This is greater than 1e-6 small, and should cause the time step to be chopped
    setSlaveStartDateAndNextReportTimeOffset(/*slave_number=*/0,
                /*start_date=*/start_date,
                /*report_time_step_size=*/time_step - tol);
    setSlaveStartDateAndNextReportTimeOffset(/*slave_number=*/1,
                 /*start_date=*/start_date+time_step,
                 /*report_time_step_size=*/time_step);
    new_time_step = rc_master.maybeChopSubStep(time_step, elapsed_time);
    BOOST_CHECK(!Opm::ReservoirCoupling::Seconds::compare_gt_or_eq(new_time_step, time_step));

}

BOOST_AUTO_TEST_SUITE_END() 
