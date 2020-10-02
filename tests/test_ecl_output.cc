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

#include <ebos/equil/equilibrationhelpers.hh>
#include <ebos/eclproblem.hh>
#include <opm/models/utils/start.hh>

#include <opm/grid/UnstructuredGrid.h>
#include <opm/grid/GridManager.hpp>

#include <opm/parser/eclipse/Units/Units.hpp>

#include <opm/io/eclipse/ESmry.hpp>

#include <opm/output/eclipse/Summary.hpp>
#include <ebos/collecttoiorank.hh>
#include <ebos/ecloutputblackoilmodule.hh>
#include <ebos/eclwriter.hh>
#include <opm/parser/eclipse/EclipseState/Schedule/Action/State.hpp>

#if HAVE_DUNE_FEM
#include <dune/fem/misc/mpimanager.hh>
#else
#include <dune/common/parallel/mpihelper.hh>
#endif

#include <array>
#include <iostream>
#include <limits>
#include <memory>
#include <numeric>
#include <sstream>
#include <string>
#include <vector>
#include <string.h>

#define CHECK(value, expected)             \
    {                                      \
         if ((value) != (expected)) {      \
             std::cerr << "Test failure: ";     \
             std::cerr << "expected value " << expected << " != " << value << std::endl; \
             throw std::runtime_error("Test failed"); \
         }                                 \
    }

#define CHECK_CLOSE(value, expected, reltol)                            \
    {                                                                   \
        if (std::fabs((expected) - (value)) > 1e-14 &&                  \
            std::fabs(((expected) - (value))/((expected) + (value))) > reltol) \
            { \
            std::cerr << "Test failure: "; \
            std::cerr << "expected value " << expected << " is not close to value " << value << std::endl; \
            throw std::runtime_error("Test failed");                    \
        } \
    }


namespace Opm::Properties {
namespace TTag {
struct TestEclOutputTypeTag {
    using InheritsFrom = std::tuple<EclBaseProblem, BlackOilModel>;
};
}
template<class TypeTag>
struct EnableGravity<TypeTag, TTag::TestEclOutputTypeTag> {
    static constexpr bool value = false;
};
template<class TypeTag>
struct EnableAsyncEclOutput<TypeTag, TTag::TestEclOutputTypeTag> {
    static constexpr bool value = false;
};

} // namespace Opm::Properties

namespace {
std::unique_ptr<Opm::EclIO::ESmry> readsum(const std::string& base)
{
    return std::make_unique<Opm::EclIO::ESmry>(base);
}

double ecl_sum_get_field_var(const Opm::EclIO::ESmry* smry,
                             const int                timeIdx,
                             const std::string&       var)
{
    return smry->get(var)[timeIdx];
}

double ecl_sum_get_general_var(const Opm::EclIO::ESmry* smry,
                               const int                timeIdx,
                               const std::string&       var)
{
    return smry->get(var)[timeIdx];
}

template <class TypeTag>
std::unique_ptr<Opm::GetPropType<TypeTag, Opm::Properties::Simulator>>
initSimulator(const char *filename)
{
    using Simulator = Opm::GetPropType<TypeTag, Opm::Properties::Simulator>;

    std::string filenameArg = "--ecl-deck-file-name=";
    filenameArg += filename;

    const char* argv[] = {
        "test_equil",
        filenameArg.c_str()
    };

    Opm::setupParameters_<TypeTag>(/*argc=*/sizeof(argv)/sizeof(argv[0]), argv, /*registerParams=*/false);

    return std::unique_ptr<Simulator>(new Simulator);
}

void test_summary()
{
    using TypeTag = Opm::Properties::TTag::TestEclOutputTypeTag;
    const std::string filename = "SUMMARY_DECK_NON_CONSTANT_POROSITY.DATA";
    const std::string casename = "SUMMARY_DECK_NON_CONSTANT_POROSITY";

    auto simulator = initSimulator<TypeTag>(filename.data());
    using Vanguard = Opm::GetPropType<TypeTag, Opm::Properties::Vanguard>;
    typedef Opm::CollectDataToIORank< Vanguard > CollectDataToIORankType;
    CollectDataToIORankType collectToIORank(simulator->vanguard());
    Opm::EclOutputBlackOilModule<TypeTag> eclOutputModule(*simulator, collectToIORank);

    typedef Opm::EclWriter<TypeTag> EclWriterType;
    // create the actual ECL writer
    std::unique_ptr<EclWriterType> eclWriter = std::unique_ptr<EclWriterType>(new EclWriterType(*simulator));

    simulator->model().applyInitialSolution();
    Opm::data::Wells dw;
    bool substep = false;
    simulator->startNextEpisode(0.0, 1e30);

    simulator->setEpisodeIndex(0);
    eclWriter->evalSummaryState(substep);
    eclWriter->writeOutput(substep);

    simulator->setEpisodeIndex(1);
    eclWriter->evalSummaryState(substep);
    eclWriter->writeOutput(substep);

    simulator->setEpisodeIndex(2);
    eclWriter->evalSummaryState(substep);
    eclWriter->writeOutput(substep);

    auto res = readsum( casename );
    const auto* resp = res.get();

    // fpr = sum_ (p * hcpv ) / hcpv, hcpv = pv * (1 - sw)
    const double fpr =  ( (3 * 0.1 + 8 * 0.2) * 500 * (1 - 0.2) ) / ( (500*0.1 + 500*0.2) * (1 - 0.2));
    CHECK_CLOSE( fpr, ecl_sum_get_field_var( resp, 1, "FPR" ) , 1e-5 );

    // foip = sum_ (b * s * pv), rs == 0;
    const double foip = ( (0.3 * 0.1 + 0.8 * 0.2) * 500 * (1 - 0.2) );
    CHECK_CLOSE(foip, ecl_sum_get_field_var( resp, 1, "FOIP" ), 1e-3 );

    // fgip = sum_ (b * pv * s), sg == 0;
    const double fgip = 0.0;
    CHECK_CLOSE(fgip, ecl_sum_get_field_var( resp, 1, "FGIP" ), 1e-3 );

    // fgip = sum_ (b * pv * s),
    const double fwip = 1.0/1000 * ( 0.1 + 0.2) * 500 * 0.2;
    CHECK_CLOSE(fwip, ecl_sum_get_field_var( resp, 1, "FWIP" ), 1e-3 );

    // region 1
    // rpr = sum_ (p * hcpv ) / hcpv, hcpv = pv * (1 - sw)
    const double rpr1 =  ( 2.5 * 0.1 * 400 * (1 - 0.2) ) / (400*0.1 * (1 - 0.2));
    CHECK_CLOSE( rpr1, ecl_sum_get_general_var( resp, 1, "RPR:1" ) , 1e-5 );
    // roip = sum_ (b * s * pv) // rs == 0;
    const double roip1 = ( 0.25 * 0.1 * 400 * (1 - 0.2) );
    CHECK_CLOSE(roip1, ecl_sum_get_general_var( resp, 1, "ROIP:1" ), 1e-3 );


    // region 2
    // rpr = sum_ (p * hcpv ) / hcpv, hcpv = pv * (1 - sw)
    const double rpr2 =  ( (5 * 0.1 * 100 + 6 * 0.2 * 100) * (1 - 0.2) ) / ( (100*0.1 + 100*0.2) * (1 - 0.2));
    CHECK_CLOSE( rpr2, ecl_sum_get_general_var( resp, 1, "RPR:2" ) , 1e-5 );
    // roip = sum_ (b * s * pv) // rs == 0;
    const double roip2 = ( (0.5 * 0.1 * 100 + 0.6 * 0.2 * 100) * (1 - 0.2) );
    CHECK_CLOSE(roip2, ecl_sum_get_general_var( resp, 1, "ROIP:2" ), 1e-3 );
}

void test_readWriteWells()
{
    using opt = Opm::data::Rates::opt;

    Opm::data::Rates r1, r2, rc1, rc2, rc3;
    r1.set( opt::wat, 5.67 );
    r1.set( opt::oil, 6.78 );
    r1.set( opt::gas, 7.89 );

    r2.set( opt::wat, 8.90 );
    r2.set( opt::oil, 9.01 );
    r2.set( opt::gas, 10.12 );

    rc1.set( opt::wat, 20.41 );
    rc1.set( opt::oil, 21.19 );
    rc1.set( opt::gas, 22.41 );

    rc2.set( opt::wat, 23.19 );
    rc2.set( opt::oil, 24.41 );
    rc2.set( opt::gas, 25.19 );

    rc3.set( opt::wat, 26.41 );
    rc3.set( opt::oil, 27.19 );
    rc3.set( opt::gas, 28.41 );

    Opm::data::Well w1, w2;
    w1.rates = r1;
    w1.bhp = 1.23;
    w1.temperature = 3.45;
    w1.control = 1;
    //w1.injectionControl = 1;
    //w1.productionControl = 1;


    /*
     *  the connection keys (active indices) and well names correspond to the
     *  input deck. All other entries in the well structures are arbitrary.
     */
    w1.connections.push_back( { 88, rc1, 30.45, 123.45, 0.0, 0.0, 0.0, 0.0, 123.456 } );
    w1.connections.push_back( { 288, rc2, 33.19, 67.89, 0.0, 0.0, 0.0, 0.0, 123.456 } );

    w2.rates = r2;
    w2.bhp = 2.34;
    w2.temperature = 4.56;
    w2.control = 1;
    //w1.injectionControl = 2;
    //w1.productionControl = 2;
    w2.connections.push_back( { 188, rc3, 36.22, 19.28, 0.0, 0.0, 0.0, 0.0, 123.456 } );

    Opm::data::Wells wellRates;

    wellRates["OP_1"] = w1;
    wellRates["OP_2"] = w2;

    typedef Dune :: Point2PointCommunicator< Dune :: SimpleMessageBuffer > P2PCommunicatorType;
    typedef typename P2PCommunicatorType :: MessageBufferType MessageBufferType;
    MessageBufferType buffer;

    wellRates.write(buffer);

    Opm::data::Wells wellRatesCopy;
    wellRatesCopy.read(buffer);

    CHECK( wellRatesCopy.get( "OP_1" , opt::wat) , wellRates.get( "OP_1" , opt::wat));
    CHECK( wellRatesCopy.get( "OP_2" , 188 , opt::wat) , wellRates.get( "OP_2" , 188 , opt::wat));
}
} // Anonymous namespace


int main(int argc, char** argv)
{
#if HAVE_DUNE_FEM
    Dune::Fem::MPIManager::initialize(argc, argv);
#else
    Dune::MPIHelper::instance(argc, argv);
#endif

    using TypeTag = Opm::Properties::TTag::TestEclOutputTypeTag;
    Opm::registerAllParameters_<TypeTag>();
    test_summary();
    test_readWriteWells();

    return 0;
}


