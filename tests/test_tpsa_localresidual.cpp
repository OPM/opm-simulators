// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  Copyright 2025, NORCE AS

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

#define BOOST_TEST_MODULE TpsaLocalResidualTests

#include <boost/test/unit_test.hpp>

#include <dune/common/fvector.hh>

#include <opm/input/eclipse/Schedule/BCProp.hpp>

#include <opm/material/materialstates/MaterialStateTPSA.hpp>

#include <opm/models/blackoil/blackoilonephaseindices.hh>
#include <opm/models/utils/parametersystem.hpp>
#include <opm/models/utils/propertysystem.hh>
#include <opm/models/utils/start.hh>

#include <opm/simulators/flow/BlackoilModelProperties.hpp>
#include <opm/simulators/flow/FlowProblemTPSA.hpp>
#include <opm/simulators/flow/TTagFlowProblemTPSA.hpp>

#if HAVE_DUNE_FEM
#include <dune/fem/misc/mpimanager.hh>
#else
#include <dune/common/parallel/mpihelper.hh>
#endif

#include <tuple>
#include <vector>

namespace Opm::Properties {
    namespace TTag {
        struct TpsaTestTypeTag {
            using InheritsFrom = std::tuple<FlowProblem, FlowProblemTpsa>;
        };
    }

    template <class TypeTag>
    struct Problem<TypeTag, TTag::TpsaTestTypeTag>
    { using type = FlowProblemTPSA<TypeTag>; };

    template <class TypeTag>
    struct EnableMech<TypeTag, TTag::TpsaTestTypeTag>
    { static constexpr bool value = true; };

    template<class TypeTag>
    struct Indices<TypeTag, TTag::TpsaTestTypeTag>
    {
    private:
        using BaseTypeTag = TTag::FlowProblem;
        using FluidSystem = GetPropType<BaseTypeTag, Properties::FluidSystem>;
        static constexpr EnergyModules energyModuleType = getPropValue<TypeTag, Properties::EnergyModuleType>();
        static constexpr int numEnergyVars = energyModuleType == EnergyModules::FullyImplicitThermal;
        static constexpr bool enableSeqImpEnergy = energyModuleType == EnergyModules::SequentialImplicitThermal;

    public:
        using type = BlackOilOnePhaseIndices<getPropValue<TypeTag, Properties::EnableSolvent>(),
                                             getPropValue<TypeTag, Properties::EnableExtbo>(),
                                             getPropValue<TypeTag, Properties::EnablePolymer>(),
                                             numEnergyVars,
                                             enableSeqImpEnergy,
                                             getPropValue<TypeTag, Properties::EnableFoam>(),
                                             getPropValue<TypeTag, Properties::EnableBrine>(),
                                             /*PVOffset=*/0,
                                             /*enabledCompIdx=*/FluidSystem::waterCompIdx,
                                             getPropValue<TypeTag, Properties::EnableBioeffects>()>;
    };
}

namespace {

template <class TypeTag>
std::unique_ptr<Opm::GetPropType<TypeTag, Opm::Properties::Simulator>>
initSimulator(const char *filename)
{
    using Simulator = Opm::GetPropType<TypeTag, Opm::Properties::Simulator>;
    const auto filenameArg = std::string {"--ecl-deck-file-name="} + filename;
    const char* argv[] = {
        "test_tpsa",
        filenameArg.c_str()
    };

    Opm::Parameters::reset();
    Opm::registerAllParameters_<TypeTag>(false);
    Opm::BlackoilModelParameters<double>::registerParameters();
    Opm::registerEclTimeSteppingParameters<double>();
    Opm::Parameters::Register<Opm::Parameters::EnableTerminalOutput>("Do *NOT* use!");
    Opm::Parameters::SetDefault<Opm::Parameters::ThreadsPerProcess>(1);
    Opm::Parameters::endRegistration();
    Opm::setupParameters_<TypeTag>(/*argc=*/sizeof(argv) / sizeof(argv[0]),
                                   argv,
                                   /*registerParams=*/false,
                                   /*allowUnused*/false,
                                   /*handleHelp*/true,
                                   /*myRank*/0);

    Opm::FlowGenericVanguard::readDeck(filename);

    return std::make_unique<Simulator>();
}

struct TpsaLocalResidualFixture
{
    // Constructor
    TpsaLocalResidualFixture()
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

struct BoundaryConditionData
{
    Opm::BCMECHType type;
    std::vector<double> displacement;
    unsigned boundaryFaceIndex;
    double faceArea;
};

}  // namespace Anonymous

BOOST_GLOBAL_FIXTURE(TpsaLocalResidualFixture);

BOOST_AUTO_TEST_CASE(TestElasticityResidual) {
    using TypeTag = Opm::Properties::TTag::TpsaTestTypeTag;
    using LocalResidual = Opm::ElasticityLocalResidual<TypeTag>;
    using Evaluation = Opm::GetPropType<TypeTag, Opm::Properties::EvaluationTPSA>;
    using MaterialState = Opm::MaterialStateTPSA<Evaluation>;
    using FluidSystem = Opm::GetPropType<TypeTag, Opm::Properties::FluidSystem>;
    using FVector = Dune::FieldVector<Evaluation, 7>;

    enum { waterPhaseIdx = FluidSystem::waterPhaseIdx };

    // Setup
    auto simulator = initSimulator<TypeTag>("tpsa_ex.data");
    simulator->model().applyInitialSolution();
    // simulator->model().invalidateAndUpdateIntensiveQuantities(/*timeIdx=*/0);
    auto& problem = simulator->problem();
    MaterialState materialStateIn;
    MaterialState materialStateEx;
    unsigned globI = 0;
    unsigned globJ = 1;

    // Set TPSA quantities
    materialStateIn.setDisplacement(0, 1.0);
    materialStateIn.setDisplacement(1, 1.0);
    materialStateIn.setDisplacement(2, 1.0);
    materialStateIn.setRotation(0, 1.0);
    materialStateIn.setRotation(1, 1.0);
    materialStateIn.setRotation(2, 1.0);
    materialStateIn.setSolidPressure(1.0);

    materialStateEx.setDisplacement(0, 2.0);
    materialStateEx.setDisplacement(1, 2.0);
    materialStateEx.setDisplacement(2, 2.0);
    materialStateEx.setRotation(0, 2.0);
    materialStateEx.setRotation(1, 2.0);
    materialStateEx.setRotation(2, 2.0);
    materialStateEx.setSolidPressure(2.0);

    // Set pressure s.t. dP in TPSA source is nonzero
    auto& sol = simulator->model().solution(/*timeIdx=*/0)[globI];
    sol[0] = 12.5e5;
    simulator->model().invalidateAndUpdateIntensiveQuantities(/*timeIdx=*/0);

    // Boundary info
    // TODO: FIXED and FREE tests
    BoundaryConditionData bcdata { Opm::BCMECHType::NONE, std::vector<double> {0.0}, 0, 0.0 };

    //
    // Volume term
    //
    FVector volTerm;
    LocalResidual::computeVolumeTerm(volTerm,
                                     materialStateIn,
                                     problem,
                                     globI);

    for (std::size_t i = 0; i < volTerm.size(); ++i) {
        if (i < 3) {
            BOOST_CHECK_CLOSE(volTerm[i].value(), 0.0, 1.0e-3);
        }
        else if (i >= 3 && i < volTerm.size() - 1) {
            // = rotation / SMODULUS
            BOOST_CHECK_CLOSE(volTerm[i].value(), 1.0 / 3.0e9, 1.0e-3);
        }
        else {
            // = solidPressure / LAME
            BOOST_CHECK_CLOSE(volTerm[i].value(), 1.0 / 2.0e9, 1.0e-3);
        }
    }

    //
    // Face term
    //
    FVector faceTerm;
    LocalResidual::computeFaceTerm(faceTerm,
                                   materialStateIn,
                                   materialStateEx,
                                   problem,
                                   globI,
                                   globJ);

    for (std::size_t i = 0; i < faceTerm.size(); ++i) {
        if (i < 2) {
            BOOST_CHECK_CLOSE(faceTerm[i].value(), -40000001.66666, 1.0e-3);
        }
        else if (i == 2) {
            BOOST_CHECK_CLOSE(faceTerm[i].value(), -39999998.33333, 1.0e-3);
        }
        else if (i == 3) {
            BOOST_CHECK_CLOSE(faceTerm[i].value(), 0.0, 1.0e-3);
        }
        else if (i == 4 || i == 6) {
            BOOST_CHECK_CLOSE(faceTerm[i].value(), -1.33333, 1.0e-3);
        }
        else {
            BOOST_CHECK_CLOSE(faceTerm[i].value(), 1.33333, 1.0e-3);
        }
    }

    //
    // Source term
    //
    FVector srcTerm;
    LocalResidual::computeSourceTerm(srcTerm,
                                     problem,
                                     globI,
                                     /*timeIdx=*/0);

    for (std::size_t i = 0; i < srcTerm.size(); ++i) {
        if (i == 6) {
            BOOST_CHECK_CLOSE(srcTerm[i].value(), -0.000125, 1.0e-2);
        }
        else {
            BOOST_CHECK_CLOSE(srcTerm[i].value(), 0.0, 1.0e-3);
        }
    }

    //
    // Boundary term
    //
    FVector bndryTerm;
    LocalResidual::computeBoundaryTerm(bndryTerm,
                                       materialStateIn,
                                       bcdata,
                                       problem,
                                       globI);
    for (std::size_t i = 0; i < bndryTerm.size(); ++i) {
        if (i < 2) {
            BOOST_CHECK_CLOSE(bndryTerm[i].value(), 120000001.0, 1.0e-6);
        }
        else if (i == 2) {
            BOOST_CHECK_CLOSE(bndryTerm[i].value(), 119999999.0, 1.0e-6);
        }
        else {
            BOOST_CHECK_CLOSE(bndryTerm[i].value(), 0.0, 1.0e-6);
        }
    }
}
