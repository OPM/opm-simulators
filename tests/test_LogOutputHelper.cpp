/*
  Copyright 2018 SINTEF Digital, Mathematics and Cybernetics.
  Copyright 2018 Equinor.

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

#define BOOST_TEST_MODULE TestLogOutputHelper

#include <boost/test/unit_test.hpp>

#include <opm/common/OpmLog/OpmLog.hpp>
#include <opm/common/OpmLog/StreamLog.hpp>
#include <opm/common/utility/String.hpp>

#include <opm/input/eclipse/Deck/Deck.hpp>
#include <opm/input/eclipse/EclipseState/EclipseState.hpp>
#include <opm/input/eclipse/Parser/Parser.hpp>
#include <opm/input/eclipse/Python/Python.hpp>
#include <opm/input/eclipse/Schedule/Schedule.hpp>
#include <opm/input/eclipse/Schedule/SummaryState.hpp>

#include <opm/simulators/flow/LogOutputHelper.hpp>

#include <memory>
#include <string>

namespace {

const std::string input = R"(
RUNSPEC
FIELD
DIMENS
  10 10 3 /
DX
    300*1000 /
DY
  300*1000 /
DZ
  100*20 100*30 100*50 /
TOPS
  100*8325 /
START             -- 0
31 AUG 1993 /
SCHEDULE
WELSPECS
-- Item #: 1   2  3 4 5  6
  'PROD'  'G1'  10  10  8400  'OIL' /
  'INJ' 'G1'  1 1 8335  'GAS' /
/
WCONPROD
-- Item #:1 2      3     4     5  9
  'PROD' 'OPEN' 'ORAT' 20000 4* 1000 /
/
WCONINJE
-- Item #:1  2   3   4  5      6  7
  'INJ' 'GAS' 'OPEN'  'RATE'  100000 1* 9014 /
/)";

std::string trimStream(std::stringstream& str)
{
    char buffer[1024];
    std::string data;
    do {
        str.getline(buffer, 1024, '\n');
        std::string tmp(buffer);
        if (!tmp.empty()) {
            tmp = Opm::trim_copy(tmp);
            data += tmp;
            data += '\n';
        }
    } while (!str.eof());

    return data;
}

}

BOOST_AUTO_TEST_CASE(Cumulative)
{
    const std::string reference = R"(=================================================== CUMULATIVE PRODUCTION/INJECTION REPORT =========================================
:  WELL  :  LOCATION :  WELL  :CTRL:    OIL    :   WATER   :    GAS    :   Prod    :    OIL    :   WATER   :    GAS    :   INJ     :
:  NAME  :  (I,J,K)  :  TYPE  :MODE:    PROD   :   PROD    :    PROD   :  RES.VOL. :    INJ    :   INJ     :    INJ    :  RES.VOL. :
:        :           :        :    :    MSTB   :   MSTB    :    MMSCF  :   MRB     :    MSTB   :   MSTB    :    MMSCF  :   MRB     :
====================================================================================================================================
:   FIELD:           :        :    :        1.0:        2.0:        3.0:        4.0:        5.0:        6.0:        7.0:        8.0:
:--------:-----------:--------:----:-----------:-----------:-----------:-----------:------------:----------:-----------:-----------:
:      G1:           :        :    :        9.0:       10.0:       11.0:       12.0:       13.0:       14.0:       15.0:       15.0:
:--------:-----------:--------:----:-----------:-----------:-----------:-----------:------------:----------:-----------:-----------:
:    PROD:   10,   10:    PROD:ORAT:       16.0:       17.0:       18.0:       19.0:       20.0:       21.0:       22.0:       23.0:
:--------:-----------:--------:----:-----------:-----------:-----------:-----------:------------:----------:-----------:-----------:
:     INJ:    1,    1:     INJ:GRAT:       24.0:       25.0:       26.0:       27.0:       28.0:       29.0:       30.0:       31.0:
:--------:-----------:--------:----:-----------:-----------:-----------:-----------:------------:----------:-----------:-----------:
)";

    std::stringstream str;
    Opm::OpmLog::addBackend("stream",
                            std::make_shared<Opm::StreamLog>(str, Opm::Log::MessageType::Note));


    Opm::Parser parser;
    auto python = std::make_shared<Opm::Python>();
    auto deck = parser.parseString(input);
    Opm::EclipseGrid grid(10,10,3);
    Opm::TableManager table ( deck );
    Opm::FieldPropsManager fp( deck, Opm::Phases{true, true, true}, grid, table);
    Opm::Runspec runspec (deck );
    Opm::Schedule schedule(deck,  grid, fp, runspec, python);

    Opm::EclipseState eclState(deck);
    Opm::SummaryState st;
    constexpr auto fields = std::array {
        std::pair{"FOPT", 1.0},
        std::pair{"FWPT", 2.0},
        std::pair{"FGPT", 3.0},
        std::pair{"FVPT", 4.0},
        std::pair{"FOIT", 5.0},
        std::pair{"FWIT", 6.0},
        std::pair{"FGIT", 7.0},
        std::pair{"FVIT", 8.0},
        std::pair{"GOPT:G1", 9.0},
        std::pair{"GWPT:G1", 10.0},
        std::pair{"GGPT:G1", 11.0},
        std::pair{"GVPT:G1", 12.0},
        std::pair{"GOIT:G1", 13.0},
        std::pair{"GWIT:G1", 14.0},
        std::pair{"GGIT:G1", 15.0},
        std::pair{"GVIT:G1", 15.0},
        std::pair{"WOPT:PROD", 16.0},
        std::pair{"WWPT:PROD", 17.0},
        std::pair{"WGPT:PROD", 18.0},
        std::pair{"WVPT:PROD", 19.0},
        std::pair{"WOIT:PROD", 20.0},
        std::pair{"WWIT:PROD", 21.0},
        std::pair{"WGIT:PROD", 22.0},
        std::pair{"WVIT:PROD", 23.0},
        std::pair{"WOPT:INJ", 24.0},
        std::pair{"WWPT:INJ", 25.0},
        std::pair{"WGPT:INJ", 26.0},
        std::pair{"WVPT:INJ", 27.0},
        std::pair{"WOIT:INJ", 28.0},
        std::pair{"WWIT:INJ", 29.0},
        std::pair{"WGIT:INJ", 30.0},
        std::pair{"WVIT:INJ", 31.0},
    };
    for (const auto& p : fields) {
        st.set(p.first, p.second * 1e3);
    }

    Opm::LogOutputHelper<double> helper(eclState, schedule, st);
    helper.cumulative(0, [](const std::string&) { return false; });
    std::string data = trimStream(str);
    BOOST_CHECK_EQUAL(data, reference);
}

BOOST_AUTO_TEST_CASE(Error)
{
    const std::string reference = R"(Finding the bubble point pressure failed for 3 cells [(2,1,1), (1,3,1), (1,4,1)]
Finding the dew point pressure failed for 3 cells [(5,1,1), (6,1,1), (7,1,1)]
)";

    std::stringstream str;
    Opm::OpmLog::addBackend("stream",
                            std::make_shared<Opm::StreamLog>(str, Opm::Log::MessageType::Warning));

    Opm::Parser parser;
    auto python = std::make_shared<Opm::Python>();
    auto deck = parser.parseString(input);
    Opm::EclipseGrid grid(10,10,3);
    Opm::TableManager table ( deck );
    Opm::FieldPropsManager fp( deck, Opm::Phases{true, true, true}, grid, table);
    Opm::Runspec runspec (deck );
    Opm::Schedule schedule(deck,  grid, fp, runspec, python);
    Opm::EclipseState eclState(deck);

    Opm::SummaryState st;
    Opm::LogOutputHelper<double> helper(eclState, schedule, st);

    str.str(""); // clear out parser errors
    helper.error({1,20,30}, {4,5,6});
    std::string data = trimStream(str);
    BOOST_CHECK_EQUAL(data, reference);
}

BOOST_AUTO_TEST_CASE(Fip)
{
    const std::string reference = R"(Field total pressure dependent pore volume = 50 RB
===================================================
:                   Field Totals                  :
:      PAV  =             0  PSIA                 :
:      PORV =           157   RB                  :
: Pressure is weighted by hydrocarbon pore volume :
: Pore volumes are taken at reference conditions  :
:--------------- Oil    STB ---------------:-- Wat    STB --:--------------- Gas   MSCF ---------------:
:      Liquid        Vapour        Total   :      Total     :      Free        Dissolved       Total   :
:------------------------:------------------------------------------:----------------:------------------------------------------:
:Currently   in place    :           132           138           120:          113   :             1             1             1:
:------------------------:------------------------------------------:----------------:------------------------------------------:
:Originally  in place    :            25            31            13:            6   :             0             0             0:
:========================:==========================================:================:==========================================:
FIPNUM report region 1 pressure dependent pore volume = 50 RB
===================================================
:        FIPNUM report region   1                 :
:      PAV  =             0  PSIA                 :
:      PORV =           371   RB                  :
:--------------- Oil    STB ---------------:-- Wat    STB --:--------------- Gas   MSCF ---------------:
:      Liquid        Vapour        Total   :      Total     :      Free        Dissolved       Total   :
:------------------------:------------------------------------------:----------------:------------------------------------------:
:Currently   in place    :           346           352           333:          327   :             2             2             2:
:------------------------:------------------------------------------:----------------:------------------------------------------:
:Originally  in place    :           239           245           226:          220   :             1             1             1:
:========================:==========================================:================:==========================================:
)";

    std::stringstream str;
    Opm::OpmLog::addBackend("stream",
                            std::make_shared<Opm::StreamLog>(str, Opm::Log::MessageType::Note));

    Opm::Parser parser;
    auto python = std::make_shared<Opm::Python>();
    auto deck = parser.parseString(input);
    Opm::EclipseGrid grid(10,10,3);
    Opm::TableManager table ( deck );
    Opm::FieldPropsManager fp( deck, Opm::Phases{true, true, true}, grid, table);
    Opm::Runspec runspec (deck );
    Opm::Schedule schedule(deck,  grid, fp, runspec, python);

    Opm::EclipseState eclState(deck);
    Opm::SummaryState st;

    Opm::LogOutputHelper<double> helper(eclState, schedule, st);
    Opm::Inplace initial, current;
    const auto& phases = current.phases();
    double j = 1.0;
    for (const auto& phase : phases) {
        initial.add(phase, j);
        initial.add("FIPNUM", phase, 0, j + 2*phases.size());
        initial.add("FIPNUM", phase, 1, j + 2*phases.size());
        current.add(phase, j + phases.size());
        current.add("FIPNUM", phase, 0, j + 3*phases.size());
        current.add("FIPNUM", phase, 1, j + 3*phases.size());
        ++j;
    }

    initial.add(Opm::Inplace::Phase::PressureHydroCarbonPV, 1.0);
    initial.add(Opm::Inplace::Phase::HydroCarbonPV, 2.0);
    initial.add(Opm::Inplace::Phase::PressurePV, 3.0);
    initial.add(Opm::Inplace::Phase::DynamicPoreVolume, 4.0);

    current.add(Opm::Inplace::Phase::PressureHydroCarbonPV, 2.0);
    current.add(Opm::Inplace::Phase::HydroCarbonPV, 4.0);
    current.add(Opm::Inplace::Phase::PressurePV, 6.0);
    current.add(Opm::Inplace::Phase::DynamicPoreVolume, 8.0);
    current.add("FIPNUM", Opm::Inplace::Phase::PressureHydroCarbonPV, 1, 2.0);
    current.add("FIPNUM", Opm::Inplace::Phase::HydroCarbonPV, 1, 4.0);
    current.add("FIPNUM", Opm::Inplace::Phase::PressurePV, 1, 6.0);
    current.add("FIPNUM", Opm::Inplace::Phase::DynamicPoreVolume, 1, 8.0);

    helper.fip(current, initial, "");
    helper.fip(current, initial, "FIPNUM");
    std::string data = trimStream(str);
    BOOST_CHECK_EQUAL(data, reference);
}

BOOST_AUTO_TEST_CASE(FipResv)
{
    const std::string reference = R"(===================================
:  RESERVOIR VOLUMES      RB      :
:---------:---------------:---------------:---------------:---------------:---------------:
: REGION  :  TOTAL PORE   :  PORE VOLUME  :  PORE VOLUME  : PORE VOLUME   :  PORE VOLUME  :
:         :   VOLUME      :  CONTAINING   :  CONTAINING   : CONTAINING    :  CONTAINING   :
:         :               :     OIL       :    WATER      :    GAS        :  HYDRO-CARBON :
:---------:---------------:---------------:---------------:---------------:---------------
:        1:            176:            170:            164:            176:            346:
:---------:---------------:---------------:---------------:---------------:---------------:
)";

    std::stringstream str;
    Opm::OpmLog::addBackend("stream",
                            std::make_shared<Opm::StreamLog>(str, Opm::Log::MessageType::Note));

    Opm::Parser parser;
    auto python = std::make_shared<Opm::Python>();
    auto deck = parser.parseString(input);
    Opm::EclipseGrid grid(10,10,3);
    Opm::TableManager table ( deck );
    Opm::FieldPropsManager fp( deck, Opm::Phases{true, true, true}, grid, table);
    Opm::Runspec runspec (deck );
    Opm::Schedule schedule(deck,  grid, fp, runspec, python);

    Opm::EclipseState eclState(deck);
    Opm::SummaryState st;

    Opm::LogOutputHelper<double> helper(eclState, schedule, st);
    Opm::Inplace current;
    const auto& phases = current.phases();
    double j = 1.0;
    for (const auto& phase : phases) {
        current.add(phase, phases.size());
        current.add("FIPNUM", phase, 1, j + phases.size());
        ++j;
    }

    current.add(Opm::Inplace::Phase::DynamicPoreVolume, 1.0);
    current.add(Opm::Inplace::Phase::OilResVolume, 2.0);
    current.add(Opm::Inplace::Phase::WaterResVolume, 3.0);
    current.add(Opm::Inplace::Phase::GasResVolume, 4.0);
    current.add("FIPNUM", Opm::Inplace::Phase::DynamicPoreVolume, 1, 11.0 + phases.size());

    helper.fipResv(current);
    std::string data = trimStream(str);
    BOOST_CHECK_EQUAL(data, reference);
}

BOOST_AUTO_TEST_CASE(Injection)
{
    const std::string reference = R"(=================================================== INJECTION REPORT ========================================
:  WELL  :  LOCATION : CTRL : CTRL : CTRL :    OIL    :   WATER   :    GAS    :   FLUID   : BHP OR : THP OR :
:  NAME  :  (I,J,K)  : MODE : MODE : MODE :    RATE   :   RATE    :    RATE   :  RES.VOL. : CON.PR.: BLK.PR.:
:        :           : OIL  : WAT  : GAS  :  STB/DAY  :  STB/DAY  :  MSCF/DAY :  RB/DAY   :  PSIA  :  PSIA  :
==============================================================================================================
:   FIELD:           :      :      :      :        1.0:        2.0:        3.0:        4.0:        :        :
:--------:-----------:------:------:------:-----------:-----------:-----------:-----------:--------:--------:
:      G1:           :      :      :      :        5.0:        6.0:        7.0:        8.0:        :        :
:--------:-----------:------:------:------:-----------:-----------:-----------:-----------:--------:--------:
:     INJ:    1,    1:      :      :  GRAT:        9.0:       10.0:       11.0:       12.0:    13.0:    14.0:
:--------:-----------:------:------:------:-----------:-----------:-----------:-----------:--------:--------:
)";

    std::stringstream str;
    Opm::OpmLog::addBackend("stream",
                            std::make_shared<Opm::StreamLog>(str, Opm::Log::MessageType::Note));

    Opm::Parser parser;
    auto python = std::make_shared<Opm::Python>();
    auto deck = parser.parseString(input);
    Opm::EclipseGrid grid(10,10,3);
    Opm::TableManager table ( deck );
    Opm::FieldPropsManager fp( deck, Opm::Phases{true, true, true}, grid, table);
    Opm::Runspec runspec (deck );
    Opm::Schedule schedule(deck,  grid, fp, runspec, python);

    Opm::EclipseState eclState(deck);
    Opm::SummaryState st;
    constexpr auto fields = std::array {
        std::pair{"FOIR", 1.0},
        std::pair{"FWIR", 2.0},
        std::pair{"FGIR", 3.0},
        std::pair{"FVIR", 4.0},
        std::pair{"GOIR:G1", 5.0},
        std::pair{"GWIR:G1", 6.0},
        std::pair{"GGIR:G1", 7.0},
        std::pair{"GVIR:G1", 8.0},
        std::pair{"WOIR:INJ", 9.0},
        std::pair{"WWIR:INJ", 10.0},
        std::pair{"WGIR:INJ", 11.0},
        std::pair{"WVIR:INJ", 12.0},
        std::pair{"WBHP:INJ", 13.0},
        std::pair{"WTHP:INJ", 14.0},
    };
    for (const auto& p : fields) {
        st.set(p.first, p.second);
    }

    Opm::LogOutputHelper<double> helper(eclState, schedule, st);
    helper.injection(0, [](const std::string&) { return false; });
    std::string data = trimStream(str);
    BOOST_CHECK_EQUAL(data, reference);
}

BOOST_AUTO_TEST_CASE(Production)
{
    const std::string reference = R"(======================================================= PRODUCTION REPORT =======================================================
:  WELL  :  LOCATION :CTRL:    OIL    :   WATER   :    GAS    :   FLUID   :   WATER   : GAS/OIL  :  WAT/GAS   : BHP OR : THP OR :
:  NAME  :  (I,J,K)  :MODE:    RATE   :   RATE    :    RATE   :  RES.VOL. :    CUT    :  RATIO   :   RATIO    : CON.PR.: BLK.PR.:
:        :           :    :  STB/DAY  :  STB/DAY  :  MSCF/DAY :  RB/DAY   :           : MSCF/STB :  STB/MSCF  :  PSIA  :  PSIA  :
=================================================================================================================================
:   FIELD:           :    :        1.0:        2.0:        3.0:        4.0:      5.000:      6.00:      0.6667:        :        :
:--------:-----------:----:-----------:-----------:-----------:-----------:-----------:----------:------------:--------:--------:
:      G1:           :    :        7.0:        8.0:        9.0:       10.0:     11.000:     12.00:      0.8889:        :        :
:--------:-----------:----:-----------:-----------:-----------:-----------:-----------:----------:------------:--------:--------:
:    PROD:   10,   10:ORAT:       13.0:       14.0:       15.0:       16.0:     17.000:     18.00:      0.9333:    19.0:    20.0:
:--------:-----------:----:-----------:-----------:-----------:-----------:-----------:----------:------------:--------:--------:
)";

    std::stringstream str;
    Opm::OpmLog::addBackend("stream",
                            std::make_shared<Opm::StreamLog>(str, Opm::Log::MessageType::Note));

    Opm::Parser parser;
    auto python = std::make_shared<Opm::Python>();
    auto deck = parser.parseString(input);
    Opm::EclipseGrid grid(10,10,3);
    Opm::TableManager table ( deck );
    Opm::FieldPropsManager fp( deck, Opm::Phases{true, true, true}, grid, table);
    Opm::Runspec runspec (deck );
    Opm::Schedule schedule(deck,  grid, fp, runspec, python);

    Opm::EclipseState eclState(deck);
    Opm::SummaryState st;
    constexpr auto fields = std::array {
        std::pair{"FOPR", 1.0},
        std::pair{"FWPR", 2.0},
        std::pair{"FGPR", 3.0},
        std::pair{"FVPR", 4.0},
        std::pair{"FWCT", 5.0},
        std::pair{"FGOR", 6.0},
        std::pair{"GOPR:G1", 7.0},
        std::pair{"GWPR:G1", 8.0},
        std::pair{"GGPR:G1", 9.0},
        std::pair{"GVPR:G1", 10.0},
        std::pair{"GWCT:G1", 11.0},
        std::pair{"GGOR:G1", 12.0},
        std::pair{"WOPR:PROD", 13.0},
        std::pair{"WWPR:PROD", 14.0},
        std::pair{"WGPR:PROD", 15.0},
        std::pair{"WVPR:PROD", 16.0},
        std::pair{"WWCT:PROD", 17.0},
        std::pair{"WGOR:PROD", 18.0},
        std::pair{"WBHP:PROD", 19.0},
        std::pair{"WTHP:PROD", 20.0},
    };
    for (const auto& p : fields) {
        st.set(p.first, p.second);
    }

    Opm::LogOutputHelper<double> helper(eclState, schedule, st);
    helper.production(0, [](const std::string&) { return false; });
    std::string data = trimStream(str);
    BOOST_CHECK_EQUAL(data, reference);
}
