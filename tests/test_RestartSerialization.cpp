/*
  Copyright 2019 SINTEF Digital, Mathematics and Cybernetics.

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

#include <ebos/ebos.hh>
#include <ebos/eclgenericvanguard.hh>
#include <ebos/femcpgridcompat.hh>

#include <opm/common/utility/Serializer.hpp>

#include <opm/input/eclipse/Schedule/Action/State.hpp>
#include <opm/input/eclipse/EclipseState/SummaryConfig/SummaryConfig.hpp>
#include <opm/input/eclipse/Schedule/SummaryState.hpp>
#include <opm/input/eclipse/Schedule/UDQ/UDQState.hpp>

#include <opm/output/eclipse/Inplace.hpp>

#include <opm/input/eclipse/EclipseState/WagHysteresisConfig.hpp>
#include <opm/material/fluidmatrixinteractions/EclHysteresisTwoPhaseLawParams.hpp>

#include <opm/models/blackoil/blackoilprimaryvariables.hh>

#include <opm/simulators/timestepping/AdaptiveTimeStepping.hpp>
#include <opm/simulators/timestepping/SimulatorReport.hpp>
#include <opm/simulators/timestepping/SimulatorTimer.hpp>
#include <opm/simulators/timestepping/TimeStepControl.hpp>
#include <opm/simulators/utils/SerializationPackers.hpp>
#include <opm/simulators/wells/ALQState.hpp>
#include <opm/simulators/wells/BlackoilWellModelGeneric.hpp>
#include <opm/simulators/wells/GroupState.hpp>
#include <opm/simulators/wells/PerfData.hpp>
#include <opm/simulators/wells/SegmentState.hpp>
#include <opm/simulators/wells/SingleWellState.hpp>

#define BOOST_TEST_MODULE TestRestartSerialization
#define BOOST_TEST_NO_MAIN

#include <boost/date_time/gregorian/gregorian.hpp>
#include <boost/test/unit_test.hpp>

#include <algorithm>

namespace Opm::Properties {
    namespace TTag {
    struct TestRestartTypeTag {
            using InheritsFrom = std::tuple<EbosTypeTag, FlowTimeSteppingParameters>;
        };
    }

    template<>
    struct LinearSolverBackend<TTag::TestRestartTypeTag, TTag::FlowIstlSolverParams> {
        using type = ISTLSolver<TTag::TestRestartTypeTag>;
    };
}

template<class T>
std::tuple<T,int,int> PackUnpack(T& in)
{
    Opm::Serialization::MemPacker packer;
    Opm::Serializer ser(packer);
    ser.pack(in);
    const size_t pos1 = ser.position();
    T out{};
    ser.unpack(out);
    const size_t pos2 = ser.position();

    return std::make_tuple(std::move(out), pos1, pos2);
}

#define TEST_FOR_TYPE_NAMED_OBJ(TYPE, NAME, OBJ) \
BOOST_AUTO_TEST_CASE(NAME) \
{ \
    auto val1 = Opm::TYPE::OBJ(); \
    auto val2 = PackUnpack(val1); \
    BOOST_CHECK_MESSAGE(std::get<1>(val2) == std::get<2>(val2), "Packed size differ from unpack size for " #TYPE); \
    BOOST_CHECK_MESSAGE(val1 == std::get<0>(val2), "Deserialized " #TYPE " differ"); \
}

#define TEST_FOR_TYPE_NAMED(TYPE, NAME) \
    TEST_FOR_TYPE_NAMED_OBJ(TYPE, NAME, serializationTestObject)

#define TEST_FOR_TYPE(TYPE) \
    TEST_FOR_TYPE_NAMED(TYPE, TYPE)

TEST_FOR_TYPE(ALQState)
TEST_FOR_TYPE(GroupState)
TEST_FOR_TYPE(HardcodedTimeStepControl)
TEST_FOR_TYPE(Inplace)
TEST_FOR_TYPE(PerfData)
TEST_FOR_TYPE(PIDAndIterationCountTimeStepControl)
TEST_FOR_TYPE(PIDTimeStepControl)
TEST_FOR_TYPE(SegmentState)
TEST_FOR_TYPE(SimpleIterationCountTimeStepControl)
TEST_FOR_TYPE(SimulatorReport)
TEST_FOR_TYPE(SimulatorReportSingle)
TEST_FOR_TYPE(SimulatorTimer)

namespace Opm { using ATS = AdaptiveTimeStepping<Properties::TTag::EbosTypeTag>; }
TEST_FOR_TYPE_NAMED_OBJ(ATS, AdaptiveTimeSteppingHardcoded, serializationTestObjectHardcoded)
TEST_FOR_TYPE_NAMED_OBJ(ATS, AdaptiveTimeSteppingPID, serializationTestObjectPID)
TEST_FOR_TYPE_NAMED_OBJ(ATS, AdaptiveTimeSteppingPIDIt, serializationTestObjectPIDIt)
TEST_FOR_TYPE_NAMED_OBJ(ATS, AdaptiveTimeSteppingSimple, serializationTestObjectSimple)

namespace Opm { using BPV = BlackOilPrimaryVariables<Properties::TTag::EbosTypeTag>; }
TEST_FOR_TYPE_NAMED(BPV, BlackoilPrimaryVariables)

namespace Opm {
    struct DummyMaterial {
        struct Params {
            struct Traits {
                using Scalar = double;
            };
        };
    };

    using HystParam = EclHysteresisTwoPhaseLawParams<DummyMaterial>;
}

TEST_FOR_TYPE_NAMED(HystParam, EclHysteresisTwoPhaseLawParams)

namespace Opm {
    using Disc = Opm::FvBaseDiscretization<Opm::Properties::TTag::EbosTypeTag>;
    using BVec = typename Disc::BlockVectorWrapper;
}
TEST_FOR_TYPE_NAMED(BVec, BlockVectorWrapper)

BOOST_AUTO_TEST_CASE(SingleWellState)
{
    Opm::ParallelWellInfo dummy;
    auto data_out = Opm::SingleWellState::serializationTestObject(dummy);
    Opm::Serialization::MemPacker packer;
    Opm::Serializer ser(packer);
    ser.pack(data_out);
    const size_t pos1 = ser.position();
    decltype(data_out) data_in("", dummy, false, 0.0, {}, Opm::PhaseUsage{}, 0.0);
    ser.unpack(data_in);
    const size_t pos2 = ser.position();
    BOOST_CHECK_MESSAGE(pos1 == pos2, "Packed size differ from unpack size for SingleWellState");
    BOOST_CHECK_MESSAGE(data_out == data_in, "Deserialized SingleWellState differ");
}

BOOST_AUTO_TEST_CASE(WellContainer)
{
    auto data_out = Opm::WellContainer<double>::serializationTestObject(1.0);
    Opm::Serialization::MemPacker packer;
    Opm::Serializer ser(packer);
    ser.pack(data_out);
    const size_t pos1 = ser.position();
    decltype(data_out) data_in;
    ser.unpack(data_in);
    const size_t pos2 = ser.position();
    BOOST_CHECK_MESSAGE(pos1 == pos2, "Packed size differ from unpack size for WellContainer");
    BOOST_CHECK_MESSAGE(data_out == data_in, "Deserialized WellContainer differ");
}

BOOST_AUTO_TEST_CASE(WellState)
{
    Opm::ParallelWellInfo dummy;
    auto data_out = Opm::WellState::serializationTestObject(dummy);
    Opm::Serialization::MemPacker packer;
    Opm::Serializer ser(packer);
    ser.pack(data_out);
    const size_t pos1 = ser.position();
    decltype(data_out) data_in(dummy);
    ser.unpack(data_in);
    const size_t pos2 = ser.position();
    BOOST_CHECK_MESSAGE(pos1 == pos2, "Packed size differ from unpack size for WellState");
    BOOST_CHECK_MESSAGE(data_out == data_in, "Deserialized WellState differ");
}

BOOST_AUTO_TEST_CASE(WGState)
{
    Opm::ParallelWellInfo dummy;
    auto data_out = Opm::WGState::serializationTestObject(dummy);
    Opm::Serialization::MemPacker packer;
    Opm::Serializer ser(packer);
    ser.pack(data_out);
    const size_t pos1 = ser.position();
    decltype(data_out) data_in(Opm::PhaseUsage{});
    data_in.well_state = Opm::WellState(dummy);
    ser.unpack(data_in);
    const size_t pos2 = ser.position();
    BOOST_CHECK_MESSAGE(pos1 == pos2, "Packed size differ from unpack size for WGState");
    BOOST_CHECK_MESSAGE(data_out == data_in, "Deserialized WGState differ");
}

BOOST_AUTO_TEST_CASE(EclGenericVanguard)
{
    auto in_params = Opm::EclGenericVanguard::serializationTestParams();
    Opm::EclGenericVanguard val1(std::move(in_params));
    Opm::Serialization::MemPacker packer;
    Opm::Serializer ser(packer);
    ser.pack(val1);
    const size_t pos1 = ser.position();
    Opm::EclGenericVanguard::SimulationModelParams out_params;
    out_params.setupTime_ = 0.0;
    out_params.actionState_ = std::make_unique<Opm::Action::State>();
    out_params.udqState_ = std::make_unique<Opm::UDQState>();
    out_params.eclSchedule_ = std::make_shared<Opm::Schedule>();
    out_params.summaryState_ = std::make_unique<Opm::SummaryState>();
    Opm::EclGenericVanguard val2(std::move(out_params));
    ser.unpack(val2);
    const size_t pos2 = ser.position();

    BOOST_CHECK_MESSAGE(pos1 == pos2, "Packed size differ from unpack size for EclGenericVanguard");
    BOOST_CHECK_MESSAGE(val1 == val2, "Deserialized EclGenericVanguard differ");
}

BOOST_AUTO_TEST_CASE(EclGenericProblem)
{
    Opm::EclipseState eclState;
    Opm::Schedule schedule;
    Dune::CpGrid grid;
#if HAVE_DUNE_FEM
    using GridPart = Dune::Fem::AdaptiveLeafGridPart<Dune::CpGrid, Dune::PartitionIteratorType(4), false>;
    using GridView = Dune::GridView<Dune::Fem::GridPart2GridViewTraits<GridPart>>;
    auto gridPart = GridPart(grid);
    auto gridView = GridView(static_cast<GridView>(gridPart));
#else
    using GridView = Dune::GridView<Dune::DefaultLeafGridViewTraits<Dune::CpGrid>>;
    auto gridView = grid.leafGridView();
#endif // HAVE_DUNE_FEM
    auto data_out
        = Opm::EclGenericProblem<GridView, Opm::BlackOilFluidSystem<double, Opm::BlackOilDefaultIndexTraits>, double>::
            serializationTestObject(eclState, schedule, gridView);
    Opm::Serialization::MemPacker packer;
    Opm::Serializer ser(packer);
    ser.pack(data_out);
    const size_t pos1 = ser.position();
    decltype(data_out) data_in(eclState, schedule, gridView);
    ser.unpack(data_in);
    const size_t pos2 = ser.position();
    BOOST_CHECK_MESSAGE(pos1 == pos2, "Packed size differ from unpack size for EclGenericProblem");
    BOOST_CHECK_MESSAGE(data_out == data_in, "Deserialized EclGenericProblem differ");
}

namespace Opm {

class BlackoilWellModelGenericTest : public BlackoilWellModelGeneric
{
public:
    BlackoilWellModelGenericTest(Schedule& schedule,
                                 const SummaryState& summaryState,
                                 const EclipseState& eclState,
                                 const PhaseUsage& phase_usage,
                                 const Parallel::Communication& comm,
                                 bool deserialize)
        : BlackoilWellModelGeneric(schedule, summaryState,
                                   eclState, phase_usage, comm)
    {
        if (deserialize) {
            active_wgstate_.well_state = WellState(dummy);
            last_valid_wgstate_.well_state = WellState(dummy);
            nupcol_wgstate_.well_state = WellState(dummy);
        }
    }

    void setSerializationTestData()
    {
        initial_step_ =  true;
        report_step_starts_ = true;
        last_run_wellpi_ = 1;
        local_shut_wells_ = {2, 3};
        closed_this_step_ = {"test1", "test2"};
        guideRate_.setSerializationTestData();
        node_pressures_ = {{"test3", 4.0}};
        active_wgstate_ = WGState::serializationTestObject(dummy);
        last_valid_wgstate_ = WGState::serializationTestObject(dummy);
        nupcol_wgstate_ = WGState::serializationTestObject(dummy);
        last_glift_opt_time_ = 5.0;
        switched_prod_groups_ = {{"test4", "test5"}};
        switched_inj_groups_ = {{{"test4", Phase::SOLVENT}, "test5"}};
    }

    void calcRates(const int, const int, const std::vector<double>&, std::vector<double>&) override
    {}

    void calcInjRates(const int, const int, std::vector<double>&) override
    {}

    void computePotentials(const std::size_t,
                           const WellState&,
                           std::string&,
                           ExceptionType::ExcEnum&,
                           DeferredLogger&) override
    {}

    void createWellContainer(const int) override
    {}

    void initWellContainer(const int) override
    {}

    void calculateProductivityIndexValuesShutWells(const int,
                                                   DeferredLogger&) override
    {}

    void calculateProductivityIndexValues(DeferredLogger&) override
    {}


    int compressedIndexForInterior(int) const override
    {
        return 0;
    }

private:
    ParallelWellInfo dummy;
};

}

BOOST_AUTO_TEST_CASE(BlackoilWellModelGeneric)
{
    Opm::Schedule schedule{};
    Opm::SummaryState summaryState{};
    Opm::EclipseState eclState{};
    Opm::PhaseUsage phase_usage{};
    Opm::Parallel::Communication comm{};
    Opm::BlackoilWellModelGenericTest data_out(schedule, summaryState,
                                               eclState, phase_usage, comm, false);
    data_out.setSerializationTestData();
    Opm::Serialization::MemPacker packer;
    Opm::Serializer ser(packer);
    ser.pack(data_out);
    const size_t pos1 = ser.position();
    Opm::BlackoilWellModelGenericTest data_in(schedule, summaryState,
                                              eclState, phase_usage, comm, true);
    ser.unpack(data_in);
    const size_t pos2 = ser.position();
    BOOST_CHECK_MESSAGE(pos1 == pos2, "Packed size differ from unpack size for BlackoilWellModelGeneric");
    BOOST_CHECK_MESSAGE(data_out == data_in, "Deserialized BlackoilWellModelGeneric differ");
}

template<class Grid, class GridView, class DofMapper, class Stencil, class Scalar>
class EclGenericTracerModelTest : public Opm::EclGenericTracerModel<Grid,GridView,DofMapper,Stencil,Scalar> {
    using Base = Opm::EclGenericTracerModel<Grid,GridView,DofMapper,Stencil,Scalar>;
public:
    EclGenericTracerModelTest(const GridView& gridView,
                              const Opm::EclipseState& eclState,
                              const Dune::CartesianIndexMapper<Grid>& cartMapper,
                              const DofMapper& dofMapper,
                              const std::function<std::array<double,Grid::dimensionworld>(int)> centroids) :
        Base(gridView, eclState, cartMapper, dofMapper, centroids)
    {}

    static EclGenericTracerModelTest
    serializationTestObject(const GridView& gridView,
                            const Opm::EclipseState& eclState,
                            const Dune::CartesianIndexMapper<Grid>& cartMapper,
                            const DofMapper& dofMapper,
                            const std::function<std::array<double,Grid::dimensionworld>(int)> centroids)
    {
        EclGenericTracerModelTest result(gridView, eclState, cartMapper, dofMapper, centroids);
        result.tracerConcentration_ = {{1.0}, {2.0}, {3.0}};
        result.wellTracerRate_.insert({{"foo", "bar"}, 4.0});

        return result;
    }

    bool operator==(const EclGenericTracerModelTest& rhs) const
    {
        if (this->tracerConcentration_.size() != rhs.tracerConcentration_.size()) {
            return false;
        }
        for (size_t i = 0; i < this->tracerConcentration_.size(); ++i) {
            if (!std::equal(this->tracerConcentration_[i].begin(),
                            this->tracerConcentration_[i].end(),
                            rhs.tracerConcentration_[i].begin(),
                            rhs.tracerConcentration_[i].end())) {
                return false;
            }
        }
        return this->wellTracerRate_ == rhs.wellTracerRate_;
    }
};

BOOST_AUTO_TEST_CASE(EclGenericTracerModel)
{
    Dune::CpGrid grid;
    Opm::EclipseState eclState;
    Dune::CartesianIndexMapper<Dune::CpGrid> mapper(grid);
    auto centroids = [](int) { return std::array<double,Dune::CpGrid::dimensionworld>{}; };
#if HAVE_DUNE_FEM
    using GridPart = Dune::Fem::AdaptiveLeafGridPart<Dune::CpGrid, Dune::PartitionIteratorType(4), false>;
    using GridView = Dune::GridView<Dune::Fem::GridPart2GridViewTraits<GridPart>>;
    auto gridPart = GridPart(grid);
    auto gridView = GridView(static_cast<GridView>(gridPart));
#else
    using GridView = Dune::GridView<Dune::DefaultLeafGridViewTraits<Dune::CpGrid>>;
    auto gridView = grid.leafGridView();
#endif // HAVE_DUNE_FEM
    Dune::MultipleCodimMultipleGeomTypeMapper<GridView> dofMapper(gridView, Dune::mcmgElementLayout());
    auto data_out = EclGenericTracerModelTest<Dune::CpGrid,
                                              GridView,
                                              Dune::MultipleCodimMultipleGeomTypeMapper<GridView>,
                                              Opm::EcfvStencil<double, GridView, false, false>,
                                              double>
        ::serializationTestObject(gridView, eclState, mapper, dofMapper, centroids);
    Opm::Serialization::MemPacker packer;
    Opm::Serializer ser(packer);
    ser.pack(data_out);
    const size_t pos1 = ser.position();
    decltype(data_out) data_in(gridView, eclState, mapper, dofMapper, centroids);
    ser.unpack(data_in);
    const size_t pos2 = ser.position();
    BOOST_CHECK_MESSAGE(pos1 == pos2, "Packed size differ from unpack size for EclGenericTracerModel");
    BOOST_CHECK_MESSAGE(data_out == data_in, "Deserialized EclGenericTracerModel differ");
}

namespace Opm {

class TBatchExport : public EclTracerModel<Properties::TTag::EbosTypeTag> {
public:
    using TBatch = TracerBatch<double>;
};

}

TEST_FOR_TYPE_NAMED(TBatchExport::TBatch, TracerBatch)

namespace {

struct AquiferFixture {
    AquiferFixture() {
        using TT = Opm::Properties::TTag::TestRestartTypeTag;
        const char* argv[] = {
            "test_RestartSerialization",
            "--ecl-deck-file-name=GLIFT1.DATA"
        };
        Opm::AdaptiveTimeStepping<TT>::registerParameters();
        Opm::setupParameters_<TT>(2, argv, /*registerParams=*/true);
        Opm::EclGenericVanguard::setCommunication(std::make_unique<Opm::Parallel::Communication>());
    }
};

}

BOOST_GLOBAL_FIXTURE(AquiferFixture);

#define TEST_FOR_AQUIFER(TYPE) \
BOOST_AUTO_TEST_CASE(TYPE) \
{ \
    using TT = Opm::Properties::TTag::TestRestartTypeTag; \
    Opm::EclGenericVanguard::readDeck("GLIFT1.DATA"); \
    using Simulator = Opm::GetPropType<TT, Opm::Properties::Simulator>; \
    Simulator sim; \
    auto data_out = Opm::TYPE<TT>::serializationTestObject(sim); \
    Opm::Serialization::MemPacker packer; \
    Opm::Serializer ser(packer); \
    ser.pack(data_out); \
    const size_t pos1 = ser.position(); \
    decltype(data_out) data_in({}, sim, {}); \
    ser.unpack(data_in); \
    const size_t pos2 = ser.position(); \
    BOOST_CHECK_MESSAGE(pos1 == pos2, "Packed size differ from unpack size for " #TYPE); \
    BOOST_CHECK_MESSAGE(data_out == data_in, "Deserialized " #TYPE " differ"); \
}

TEST_FOR_AQUIFER(AquiferCarterTracy)
TEST_FOR_AQUIFER(AquiferFetkovich)

BOOST_AUTO_TEST_CASE(AquiferNumerical)
{
    using TT = Opm::Properties::TTag::TestRestartTypeTag;
    Opm::EclGenericVanguard::readDeck("GLIFT1.DATA");
    using Simulator = Opm::GetPropType<TT, Opm::Properties::Simulator>;
    Simulator sim;
    auto data_out = Opm::AquiferNumerical<TT>::serializationTestObject(sim);
    Opm::Serialization::MemPacker packer;
    Opm::Serializer ser(packer);
    ser.pack(data_out);
    const size_t pos1 = ser.position();
    decltype(data_out) data_in({}, sim);
    ser.unpack(data_in);
    const size_t pos2 = ser.position();
    BOOST_CHECK_MESSAGE(pos1 == pos2, "Packed size differ from unpack size for AquiferNumerical");
    BOOST_CHECK_MESSAGE(data_out == data_in, "Deserialized AquiferNumerical differ");
}

BOOST_AUTO_TEST_CASE(AquiferConstantFlux)
{
    using TT = Opm::Properties::TTag::TestRestartTypeTag;
    Opm::EclGenericVanguard::readDeck("GLIFT1.DATA");
    using Simulator = Opm::GetPropType<TT, Opm::Properties::Simulator>;
    Simulator sim;
    auto data_out = Opm::AquiferConstantFlux<TT>::serializationTestObject(sim);
    Opm::Serialization::MemPacker packer;
    Opm::Serializer ser(packer);
    ser.pack(data_out);
    const size_t pos1 = ser.position();
    decltype(data_out) data_in({}, sim, {});
    ser.unpack(data_in);
    const size_t pos2 = ser.position();
    BOOST_CHECK_MESSAGE(pos1 == pos2, "Packed size differ from unpack size for AquiferConstantFlux");
    BOOST_CHECK_MESSAGE(data_out == data_in, "Deserialized AquiferConstantFlux differ");
}

bool init_unit_test_func()
{
    return true;
}

int main(int argc, char** argv)
{
    // MPI setup.
    int argcDummy = 1;
    const char *tmp[] = {"test_RestartSerialization"};
    char **argvDummy = const_cast<char**>(tmp);
#if HAVE_DUNE_FEM
    Dune::Fem::MPIManager::initialize(argcDummy, argvDummy);
#else
    Dune::MPIHelper::instance(argcDummy, argvDummy);
#endif

    return boost::unit_test::unit_test_main(&init_unit_test_func, argc, argv);
}
