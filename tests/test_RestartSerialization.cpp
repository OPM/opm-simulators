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

#include <opm/common/utility/Serializer.hpp>

#include <opm/input/eclipse/Schedule/Action/State.hpp>
#include <opm/input/eclipse/EclipseState/SummaryConfig/SummaryConfig.hpp>
#include <opm/input/eclipse/Schedule/SummaryState.hpp>
#include <opm/input/eclipse/Schedule/UDQ/UDQState.hpp>

#include <opm/models/blackoil/blackoilprimaryvariables.hh>

#include <opm/simulators/timestepping/AdaptiveTimeSteppingEbos.hpp>
#include <opm/simulators/timestepping/SimulatorTimer.hpp>
#include <opm/simulators/timestepping/TimeStepControl.hpp>
#include <opm/simulators/utils/SerializationPackers.hpp>

#define BOOST_TEST_MODULE TestRestartSerialization
#define BOOST_TEST_NO_MAIN

#include <boost/date_time/gregorian/gregorian.hpp>
#include <boost/test/unit_test.hpp>

template<class T>
std::tuple<T,int,int> PackUnpack(T& in)
{
    Opm::Serialization::MemPacker packer;
    Opm::Serializer ser(packer);
    ser.pack(in);
    size_t pos1 = ser.position();
    T out{};
    ser.unpack(out);
    size_t pos2 = ser.position();

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

TEST_FOR_TYPE(HardcodedTimeStepControl)
TEST_FOR_TYPE(PIDAndIterationCountTimeStepControl)
TEST_FOR_TYPE(PIDTimeStepControl)
TEST_FOR_TYPE(SimpleIterationCountTimeStepControl)
TEST_FOR_TYPE(SimulatorTimer)

namespace Opm { using ATE = AdaptiveTimeSteppingEbos<Properties::TTag::EbosTypeTag>; }
TEST_FOR_TYPE_NAMED_OBJ(ATE, AdaptiveTimeSteppingEbosHardcoded, serializationTestObjectHardcoded)
TEST_FOR_TYPE_NAMED_OBJ(ATE, AdaptiveTimeSteppingEbosPID, serializationTestObjectPID)
TEST_FOR_TYPE_NAMED_OBJ(ATE, AdaptiveTimeSteppingEbosPIDIt, serializationTestObjectPIDIt)
TEST_FOR_TYPE_NAMED_OBJ(ATE, AdaptiveTimeSteppingEbosSimple, serializationTestObjectSimple)

namespace Opm { using BPV = BlackOilPrimaryVariables<Properties::TTag::EbosTypeTag>; }
TEST_FOR_TYPE_NAMED(BPV, BlackoilPrimaryVariables)

namespace Opm {
    using Disc = Opm::FvBaseDiscretization<Opm::Properties::TTag::EbosTypeTag>;
    using BVec = typename Disc::BlockVectorWrapper;
}
TEST_FOR_TYPE_NAMED(BVec, BlockVectorWrapper)

BOOST_AUTO_TEST_CASE(EclGenericVanguard)
{
    auto in_params = Opm::EclGenericVanguard::serializationTestParams();
    Opm::EclGenericVanguard val1(in_params);
    Opm::Serialization::MemPacker packer;
    Opm::Serializer ser(packer);
    ser.pack(val1);
    size_t pos1 = ser.position();
    Opm::EclGenericVanguard::SetupParams out_params;
    out_params.setupTime_ = 0.0;
    out_params.actionState_ = std::make_unique<Opm::Action::State>();
    out_params.udqState_ = std::make_unique<Opm::UDQState>();
    out_params.eclSchedule_ = std::make_shared<Opm::Schedule>();
    out_params.summaryState_ = std::make_unique<Opm::SummaryState>();
    Opm::EclGenericVanguard val2(out_params);
    ser.unpack(val2);
    size_t pos2 = ser.position();

    BOOST_CHECK_MESSAGE(pos1 == pos2, "Packed size differ from unpack size for EclGenericVanguard");
    BOOST_CHECK_MESSAGE(val1 == val2, "Deserialized EclGenericVanguard differ");
}

BOOST_AUTO_TEST_CASE(EclGenericProblem)
{
    Opm::EclipseState eclState;
    Opm::Schedule schedule;
    Dune::CpGrid grid;
    auto data_out = Opm::EclGenericProblem<Dune::GridView<Dune::DefaultLeafGridViewTraits<Dune::CpGrid>>,
                                           Opm::BlackOilFluidSystem<double,Opm::BlackOilDefaultIndexTraits>,
                                           double>::
        serializationTestObject(eclState, schedule, grid.leafGridView());
    Opm::Serialization::MemPacker packer;
    Opm::Serializer ser(packer);
    ser.pack(data_out);
    decltype(data_out) data_in(eclState, schedule, grid.leafGridView());
    ser.unpack(data_in);
    BOOST_CHECK_MESSAGE(data_out == data_in, "Deserialized EclGenericProblem differ");
}

namespace {

struct AquiferFixture {
    AquiferFixture() {
        using TT = Opm::Properties::TTag::EbosTypeTag;
        const char* argv[] = {
            "test_RestartSerialization",
            "--ecl-deck-file-name=GLIFT1.DATA"
        };
        Opm::setupParameters_<TT>(2, argv, /*registerParams=*/true);
        Opm::EclGenericVanguard::setCommunication(std::make_unique<Opm::Parallel::Communication>());
    }
};

}

BOOST_GLOBAL_FIXTURE(AquiferFixture);

BOOST_AUTO_TEST_CASE(AquiferCarterTracy)
{
    using TT = Opm::Properties::TTag::EbosTypeTag;
    Opm::EclGenericVanguard::readDeck("GLIFT1.DATA");
    using Simulator = Opm::GetPropType<TT, Opm::Properties::Simulator>;
    Simulator sim;
    auto data_out = Opm::AquiferCarterTracy<TT>::serializationTestObject(sim);
    Opm::Serialization::MemPacker packer;
    Opm::Serializer ser(packer);
    ser.pack(data_out);
    size_t pos1 = ser.position();
    decltype(data_out) data_in({}, sim, {});
    ser.unpack(data_in);
    size_t pos2 = ser.position();
    BOOST_CHECK_MESSAGE(pos1 == pos2, "Packed size differ from unpack size for AquiferCarterTracy");
    BOOST_CHECK_MESSAGE(data_out == data_in, "Deserialized AquiferCarterTracy differ");
}

BOOST_AUTO_TEST_CASE(AquiferFetkovich)
{
    using TT = Opm::Properties::TTag::EbosTypeTag;
    Opm::EclGenericVanguard::readDeck("GLIFT1.DATA");
    using Simulator = Opm::GetPropType<TT, Opm::Properties::Simulator>;
    Simulator sim;
    auto data_out = Opm::AquiferFetkovich<TT>::serializationTestObject(sim);
    Opm::Serialization::MemPacker packer;
    Opm::Serializer ser(packer);
    ser.pack(data_out);
    size_t pos1 = ser.position();
    decltype(data_out) data_in({}, sim, {});
    ser.unpack(data_in);
    size_t pos2 = ser.position();
    BOOST_CHECK_MESSAGE(pos1 == pos2, "Packed size differ from unpack size for AquiferFetkovich");
    BOOST_CHECK_MESSAGE(data_out == data_in, "Deserialized AquiferFetkovich differ");
}

BOOST_AUTO_TEST_CASE(AquiferNumerical)
{
    using TT = Opm::Properties::TTag::EbosTypeTag;
    Opm::EclGenericVanguard::readDeck("GLIFT1.DATA");
    using Simulator = Opm::GetPropType<TT, Opm::Properties::Simulator>;
    Simulator sim;
    auto data_out = Opm::AquiferNumerical<TT>::serializationTestObject(sim);
    Opm::Serialization::MemPacker packer;
    Opm::Serializer ser(packer);
    ser.pack(data_out);
    size_t pos1 = ser.position();
    decltype(data_out) data_in({}, sim);
    ser.unpack(data_in);
    size_t pos2 = ser.position();
    BOOST_CHECK_MESSAGE(pos1 == pos2, "Packed size differ from unpack size for AquiferFetkovich");
    BOOST_CHECK_MESSAGE(data_out == data_in, "Deserialized AquiferFetkovich differ");
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
        result.tracerResidual_ = {{1.0}};
        result.wellTracerRate_.insert({{"foo", "bar"}, 4.0});

        return result;
    }

    bool operator==(const EclGenericTracerModelTest& rhs) const
    {
        if (this->tracerResidual_.size() != rhs.tracerResidual_.size()) {
            return false;
        }
        for (size_t i = 0; i < this->tracerResidual_.size(); ++i) {
            if (this->tracerResidual_[i] != rhs.tracerResidual_[i]) {
                return false;
            }
        }
        if (this->tracerConcentration_.size() != rhs.tracerConcentration_.size()) {
            return false;
        }
        for (size_t i = 0; i < this->tracerConcentration_.size(); ++i) {
            if (this->tracerConcentration_[i].size() != rhs.tracerConcentration_[i].size()) {
                return false;
            }
            for (size_t j = 0; j < this->tracerConcentration_[i].size(); ++j) {
                if (this->tracerConcentration_[i][j] != rhs.tracerConcentration_[i][j]) {
                    return false;
                }
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
    Dune::MultipleCodimMultipleGeomTypeMapper<Dune::CpGrid::LeafGridView> dofMapper(grid.leafGridView(), Dune::mcmgElementLayout());
    auto centroids = [](int) { return std::array<double,Dune::CpGrid::dimensionworld>{}; };
    auto data_out = EclGenericTracerModelTest<Dune::CpGrid,
                                              Dune::GridView<Dune::DefaultLeafGridViewTraits<Dune::CpGrid>>,
                                              Dune::MultipleCodimMultipleGeomTypeMapper<Dune::GridView<Dune::DefaultLeafGridViewTraits<Dune::CpGrid>>>,
                                              Opm::EcfvStencil<double,Dune::GridView<Dune::DefaultLeafGridViewTraits<Dune::CpGrid>>,false,false>,
                                              double>::serializationTestObject(grid.leafGridView(), eclState, mapper, dofMapper, centroids);
    Opm::Serialization::MemPacker packer;
    Opm::Serializer ser(packer);
    ser.pack(data_out);
    size_t pos1 = ser.position();
    decltype(data_out) data_in(grid.leafGridView(), eclState, mapper, dofMapper, centroids);
    ser.unpack(data_in);
    size_t pos2 = ser.position();
    BOOST_CHECK_MESSAGE(pos1 == pos2, "Packed size differ from unpack size for EclGenericTracerModel");
    BOOST_CHECK_MESSAGE(data_out == data_in, "Deserialized EclGenericTracerModel differ");
}

bool init_unit_test_func()
{
    return true;
}

int main(int argc, char** argv)
{
    Dune::MPIHelper::instance(argc, argv);
    return boost::unit_test::unit_test_main(&init_unit_test_func, argc, argv);
}
