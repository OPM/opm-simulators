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

#include <opm/simulators/flow/LogOutputHelper.hpp>

#include <opm/common/OpmLog/OpmLog.hpp>
#include <opm/common/OpmLog/StreamLog.hpp>

#include <opm/common/utility/String.hpp>

#include <opm/input/eclipse/EclipseState/EclipseState.hpp>
#include <opm/input/eclipse/Schedule/Schedule.hpp>
#include <opm/input/eclipse/Schedule/SummaryState.hpp>

#include <opm/input/eclipse/Deck/Deck.hpp>

#include <opm/input/eclipse/Parser/Parser.hpp>

#include <array>
#include <memory>
#include <sstream>
#include <string>

namespace {

    Opm::Deck deck()
    {
        return Opm::Parser{}.parseString(R"(RUNSPEC
OIL
GAS
WATER
FIELD
DIMENS
  10 10 3 /
DXV
  10*1000.0 /
DYV
  10*1000.0 /
DZV
  20.0 30.0 50.0 /
TOPS
  100*8325.0 /
GRID
PORO
  300*0.3 /
PERMX
  300*0.3 /
PERMY
  300*0.3 /
PERMZ
  300*0.3 /
START        -- 0
31 AUG 1993 /
SCHEDULE
WELSPECS
-- Item #: 1   2  3 4 5  6
  'PROD'  'G1'  10  10  8400  'OIL' /
  'INJ' 'G1'  1 1 8335  'GAS' /
/
COMPDAT
   'INJ'   1   1   1    1    OPEN     1*    1.172656E+2   0.21600   1*   0.00000   1*   'Z' /
   'PROD'  2   2   1    1    OPEN     1*    1.172656E+2   0.21600   1*   0.00000   1*   'Z' /
/
WCONPROD
-- Item #:1 2      3     4     5  9
  'PROD' 'OPEN' 'ORAT' 20000 4* 1000 /
/
WCONINJE
-- Item #:1  2   3   4  5      6  7
  'INJ' 'GAS' 'OPEN'  'RATE'  100000 1* 9014 /
/
END
)");
    }

    Opm::Deck mswDeck()
    {
        return Opm::Parser{}.parseString(R"(RUNSPEC
OIL
GAS
WATER
FIELD
DIMENS
  10 10 3 /
DXV
  10*1000.0 /
DYV
  10*1000.0 /
DZV
  20.0 30.0 50.0 /
TOPS
  100*8325.0 /
GRID
PORO
  300*0.3 /
PERMX
  300*0.3 /
PERMY
  300*0.3 /
PERMZ
  300*0.3 /
START        -- 0
31 AUG 1993 /
SCHEDULE
WELSPECS
  'PROD'  'G1'  10  10  8400  'OIL' /
  'INJ' 'G1'  1 1 8335  'GAS' /
/
WELSEGS
PROD          123.4          0.0000       1*        INC         'HF-' /
   2             2            1              1              49.36537      3.56589          0.15200     0.00001 /
/
COMPSEGS
   PROD /
   2     2     1     2             0.00000          0.10000         /
/
END
)");
    }

    std::string trimStream(std::stringstream& str)
    {
        std::string data;

        std::array<char, 1024> buffer{};
        do {
            str.getline(buffer.data(), buffer.size(), '\n');

            const auto tmp = std::string { buffer.data() };
            if (!tmp.empty()) {
                data += Opm::trim_copy(tmp);
                data += '\n';
            }
        } while (!str.eof());

        return data;
    }

    template <int type, bool msw = false>
    struct LogFixture
    {
        LogFixture() : LogFixture { msw ? mswDeck() : deck() } {}

        explicit LogFixture(const Opm::Deck& deck)
            : eclState { deck }
            , schedule { deck, eclState }
        {
            Opm::OpmLog::addBackend("stream", std::make_shared<Opm::StreamLog>(str, type));
        }

        ~LogFixture()
        {
            Opm::OpmLog::removeBackend("stream");
        }

        Opm::EclipseState eclState{};
        Opm::Schedule schedule{};
        Opm::SummaryState st{};

        std::stringstream str;
    };

    using LogNoteFixture = LogFixture<Opm::Log::MessageType::Note>;
    using LogNoteFixtureMSW = LogFixture<Opm::Log::MessageType::Note, true>;
    using LogWarningFixture = LogFixture<Opm::Log::MessageType::Warning>;

} // Anonymous namespace

BOOST_FIXTURE_TEST_CASE(Cumulative, LogNoteFixture)
{
    const auto reference = std::string {
        R"(============================================== CUMULATIVE PRODUCTION/INJECTION TOTALS ==============================================
:  WELL  : LOCATION  :  WELL  :CTRL:    OIL    :   WATER   :    GAS    :   Prod    :    OIL    :   WATER   :    GAS    :    INJ    :
:  NAME  :  (I,J,K)  :  TYPE  :MODE:   PROD    :   PROD    :   PROD    : RES.VOL.  :    INJ    :    INJ    :    INJ    : RES.VOL.  :
:        :           :        :    :   MSTB    :   MSTB    :   MMSCF   :    MRB    :   MSTB    :   MSTB    :   MMSCF   :    MRB    :
====================================================================================================================================
:FIELD   :           :        :    :        1.0:        2.0:        3.0:        4.0:        5.0:        6.0:        7.0:        8.0:
:G1      :           :        :    :        9.0:       10.0:       11.0:       12.0:       13.0:       14.0:       15.0:       15.0:
:PROD    : 10, 10    :    PROD:ORAT:       16.0:       17.0:       18.0:       19.0:       20.0:       21.0:       22.0:       23.0:
:INJ     :  1,  1    :     INJ:GRAT:       24.0:       25.0:       26.0:       27.0:       28.0:       29.0:       30.0:       31.0:
:--------:-----------:--------:----:-----------:-----------:-----------:-----------:-----------:-----------:-----------:-----------:
)"
    };

    // Note: Cumulative gas values--e.g., FGPT--multiplied by an additional
    // factor of 1000, for a total multiplicative factor of one million, in
    // order to produce the expected balance sheet output in MM* units.
    st.set("FOPT", 1.0e3);
    st.set("FWPT", 2.0e3);
    st.set("FGPT", 3.0e6);
    st.set("FVPT", 4.0e3);
    st.set("FOIT", 5.0e3);
    st.set("FWIT", 6.0e3);
    st.set("FGIT", 7.0e6);
    st.set("FVIT", 8.0e3);

    st.update_group_var("G1", "GOPT",  9.0e3);
    st.update_group_var("G1", "GWPT", 10.0e3);
    st.update_group_var("G1", "GGPT", 11.0e6);
    st.update_group_var("G1", "GVPT", 12.0e3);
    st.update_group_var("G1", "GOIT", 13.0e3);
    st.update_group_var("G1", "GWIT", 14.0e3);
    st.update_group_var("G1", "GGIT", 15.0e6);
    st.update_group_var("G1", "GVIT", 15.0e3);

    st.update_well_var("PROD", "WOPT", 16.0e3);
    st.update_well_var("PROD", "WWPT", 17.0e3);
    st.update_well_var("PROD", "WGPT", 18.0e6);
    st.update_well_var("PROD", "WVPT", 19.0e3);
    st.update_well_var("PROD", "WOIT", 20.0e3);
    st.update_well_var("PROD", "WWIT", 21.0e3);
    st.update_well_var("PROD", "WGIT", 22.0e6);
    st.update_well_var("PROD", "WVIT", 23.0e3);

    st.update_well_var("INJ", "WOPT", 24.0e3);
    st.update_well_var("INJ", "WWPT", 25.0e3);
    st.update_well_var("INJ", "WGPT", 26.0e6);
    st.update_well_var("INJ", "WVPT", 27.0e3);
    st.update_well_var("INJ", "WOIT", 28.0e3);
    st.update_well_var("INJ", "WWIT", 29.0e3);
    st.update_well_var("INJ", "WGIT", 30.0e6);
    st.update_well_var("INJ", "WVIT", 31.0e3);

    Opm::LogOutputHelper<double> helper(eclState, schedule, st, "dummy version");
    helper.cumulative(0, false);

    const auto data = trimStream(str);
    BOOST_CHECK_EQUAL(data, reference);
}

BOOST_FIXTURE_TEST_CASE(CumulativeW2, LogNoteFixture)
{
    const auto reference = std::string {
        R"(============================================== CUMULATIVE PRODUCTION/INJECTION TOTALS ==============================================
:  WELL  : LOCATION  :  WELL  :CTRL:    OIL    :   WATER   :    GAS    :   Prod    :    OIL    :   WATER   :    GAS    :    INJ    :
:  NAME  :  (I,J,K)  :  TYPE  :MODE:   PROD    :   PROD    :   PROD    : RES.VOL.  :    INJ    :    INJ    :    INJ    : RES.VOL.  :
:        :           :        :    :   MSTB    :   MSTB    :   MMSCF   :    MRB    :   MSTB    :   MSTB    :   MMSCF   :    MRB    :
====================================================================================================================================
:FIELD   :           :        :    :        1.0:        2.0:        3.0:        4.0:        5.0:        6.0:        7.0:        8.0:
:G1      :           :        :    :        9.0:       10.0:       11.0:       12.0:       13.0:       14.0:       15.0:       15.0:
:PROD    : 10, 10    :    PROD:ORAT:       16.0:       17.0:       18.0:       19.0:       20.0:       21.0:       22.0:       23.0:
:  BLOCK :  2,  2,  1:        :    :       40.0:       41.0:       42.0:       43.0:       44.0:       45.0:       46.0:       47.0:
:INJ     :  1,  1    :     INJ:GRAT:       24.0:       25.0:       26.0:       27.0:       28.0:       29.0:       30.0:       31.0:
:  BLOCK :  1,  1,  1:        :    :       32.0:       33.0:       34.0:       35.0:       36.0:       37.0:       38.0:       39.0:
:--------:-----------:--------:----:-----------:-----------:-----------:-----------:-----------:-----------:-----------:-----------:
)"
    };

    // Note: Cumulative gas values--e.g., FGPT--multiplied by an additional
    // factor of 1000, for a total multiplicative factor of one million, in
    // order to produce the expected balance sheet output in MM* units.
    st.set("FOPT", 1.0e3);
    st.set("FWPT", 2.0e3);
    st.set("FGPT", 3.0e6);
    st.set("FVPT", 4.0e3);
    st.set("FOIT", 5.0e3);
    st.set("FWIT", 6.0e3);
    st.set("FGIT", 7.0e6);
    st.set("FVIT", 8.0e3);

    st.update_group_var("G1", "GOPT",  9.0e3);
    st.update_group_var("G1", "GWPT", 10.0e3);
    st.update_group_var("G1", "GGPT", 11.0e6);
    st.update_group_var("G1", "GVPT", 12.0e3);
    st.update_group_var("G1", "GOIT", 13.0e3);
    st.update_group_var("G1", "GWIT", 14.0e3);
    st.update_group_var("G1", "GGIT", 15.0e6);
    st.update_group_var("G1", "GVIT", 15.0e3);

    st.update_well_var("PROD", "WOPT", 16.0e3);
    st.update_well_var("PROD", "WWPT", 17.0e3);
    st.update_well_var("PROD", "WGPT", 18.0e6);
    st.update_well_var("PROD", "WVPT", 19.0e3);
    st.update_well_var("PROD", "WOIT", 20.0e3);
    st.update_well_var("PROD", "WWIT", 21.0e3);
    st.update_well_var("PROD", "WGIT", 22.0e6);
    st.update_well_var("PROD", "WVIT", 23.0e3);

    st.update_well_var("INJ", "WOPT", 24.0e3);
    st.update_well_var("INJ", "WWPT", 25.0e3);
    st.update_well_var("INJ", "WGPT", 26.0e6);
    st.update_well_var("INJ", "WVPT", 27.0e3);
    st.update_well_var("INJ", "WOIT", 28.0e3);
    st.update_well_var("INJ", "WWIT", 29.0e3);
    st.update_well_var("INJ", "WGIT", 30.0e6);
    st.update_well_var("INJ", "WVIT", 31.0e3);

    st.update_conn_var("INJ", "COPT", 1, 32.0e3);
    st.update_conn_var("INJ", "CWPT", 1, 33.0e3);
    st.update_conn_var("INJ", "CGPT", 1, 34.0e6);
    st.update_conn_var("INJ", "CVPT", 1, 35.0e3);
    st.update_conn_var("INJ", "COIT", 1, 36.0e3);
    st.update_conn_var("INJ", "CWIT", 1, 37.0e3);
    st.update_conn_var("INJ", "CGIT", 1, 38.0e6);
    st.update_conn_var("INJ", "CVIT", 1, 39.0e3);

    st.update_conn_var("PROD", "COPT", 12, 40.0e3);
    st.update_conn_var("PROD", "CWPT", 12, 41.0e3);
    st.update_conn_var("PROD", "CGPT", 12, 42.0e6);
    st.update_conn_var("PROD", "CVPT", 12, 43.0e3);
    st.update_conn_var("PROD", "COIT", 12, 44.0e3);
    st.update_conn_var("PROD", "CWIT", 12, 45.0e3);
    st.update_conn_var("PROD", "CGIT", 12, 46.0e6);
    st.update_conn_var("PROD", "CVIT", 12, 47.0e3);

    Opm::LogOutputHelper<double> helper(eclState, schedule, st, "dummy version");
    helper.cumulative(0, true);

    const auto data = trimStream(str);
    BOOST_CHECK_EQUAL(data, reference);
}

BOOST_FIXTURE_TEST_CASE(Error, LogWarningFixture)
{
    const auto reference = std::string {
        R"(Finding the bubble point pressure failed for 3 cells [(2,1,1), (1,3,1), (1,4,1)]
Finding the dew point pressure failed for 3 cells [(5,1,1), (6,1,1), (7,1,1)]
)"
    };

    Opm::LogOutputHelper<double> helper(eclState, schedule, st, "dummy version");

    str.str(""); // clear out parser errors
    helper.error({1,20,30}, {4,5,6});

    const auto data = trimStream(str);
    BOOST_CHECK_EQUAL(data, reference);
}

BOOST_FIXTURE_TEST_CASE(Fip, LogNoteFixture)
{
    const auto reference = std::string {
        R"(
                                                     ==================================================
                                                     :                  FIELD TOTALS                  :
                                                     :        PAV  =             0  PSIA              :
                                                     :        PORV =           157  RB                :
                                                     : Pressure is weighted by hydrocarbon pore volume:
                                                     : Pore volumes are taken at reference conditions :
                           :--------------- OIL  STB ------------------:-- WAT    STB --:--------------  GAS    MSCF ---------------:
                           :    LIQUID         VAPOUR        TOTAL     :     TOTAL      :     FREE        DISSOLVED       TOTAL     :
 :-------------------------:-------------------------------------------:----------------:-------------------------------------------:
 :CURRENTLY IN PLACE       :      132           138           120      :      113       :       1             1             1       :
 :-------------------------:-------------------------------------------:----------------:-------------------------------------------:
 :ORIGINALLY IN PLACE      :      25             31            13      :       6        :       0             0             0       :
 ====================================================================================================================================



                                                     ==================================================
                                                     :              FIPNUM REPORT REGION  1           :
                                                     :        PAV  =             0  PSIA              :
                                                     :        PORV =           371  RB                :
                           :--------------- OIL  STB ------------------:-- WAT    STB --:--------------  GAS    MSCF ---------------:
                           :    LIQUID         VAPOUR        TOTAL     :     TOTAL      :     FREE        DISSOLVED       TOTAL     :
 :-------------------------:-------------------------------------------:----------------:-------------------------------------------:
 :CURRENTLY IN PLACE       :      346           352           333      :      327       :       2             2             2       :
 :-------------------------:-------------------------------------------:----------------:-------------------------------------------:
 :ORIGINALLY IN PLACE      :      239           245           226      :      220       :       1             1             1       :
 :-------------------------:-------------------------------------------:----------------:-------------------------------------------:
 ====================================================================================================================================


)"
    };

    Opm::Inplace initial, current;
    {
        const auto offset = 17.0;

        double j = 1.0;
        for (const auto& phase : current.phases()) {
            initial.add(phase, j);
            initial.add("FIPNUM", phase, 0, j + 2*offset);
            initial.add("FIPNUM", phase, 1, j + 2*offset);

            current.add(phase, j + offset);
            current.add("FIPNUM", phase, 0, j + 3*offset);
            current.add("FIPNUM", phase, 1, j + 3*offset);

            ++j;
        }
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

    Opm::LogOutputHelper<double> helper(eclState, schedule, st, "dummy version");
    helper.fip(current, initial, "");
    helper.fip(current, initial, "FIPNUM");

    BOOST_CHECK_EQUAL(str.str(), reference);
}

BOOST_FIXTURE_TEST_CASE(FipResv, LogNoteFixture)
{
    const auto reference = std::string {
        R"(
                                                     ===================================
                                                     :  RESERVOIR VOLUMES      RB      :
 :---------:---------------:---------------:---------------:---------------:---------------:
 : REGION  :  TOTAL PORE   :  PORE VOLUME  :  PORE VOLUME  :  PORE VOLUME  :  PORE VOLUME  :
 :         :    VOLUME     :  CONTAINING   :  CONTAINING   :  CONTAINING   :  CONTAINING   :
 :         :               :      OIL      :     WATER     :      GAS      : HYDRO-CARBON  :
 :---------:---------------:---------------:---------------:---------------:---------------:
 :FIELD    :            176:             13:             19:             25:             38:
 :1        :            176:            170:            164:            176:            346:
 ===========================================================================================
)"
    };

    const auto offset = 17.0;

    Opm::Inplace current;
    {
        double j = 1.0;
        for (const auto& phase : current.phases()) {
            current.add(phase, offset);
            current.add("FIPNUM", phase, 1, j + offset);
            ++j;
        }
    }

    current.add(Opm::Inplace::Phase::DynamicPoreVolume, 1.0);
    current.add(Opm::Inplace::Phase::OilResVolume, 2.0);
    current.add(Opm::Inplace::Phase::WaterResVolume, 3.0);
    current.add(Opm::Inplace::Phase::GasResVolume, 4.0);
    current.add("FIPNUM", Opm::Inplace::Phase::DynamicPoreVolume, 1, 11.0 + offset);

    Opm::LogOutputHelper<double> helper(eclState, schedule, st, "dummy version");
    helper.fipResv(current, "FIPNUM");

    BOOST_CHECK_EQUAL(str.str(), reference);
}

BOOST_FIXTURE_TEST_CASE(Injection, LogNoteFixture)
{
    const auto reference = std::string {
        R"(============================================= INJECTION REPORT ==============================================
:  WELL  : LOCATION  : CTRL : CTRL : CTRL :    OIL    :   WATER   :    GAS    :   FLUID   : BHP OR : THP OR :
:  NAME  :  (I,J,K)  : MODE : MODE : MODE :   RATE    :   RATE    :   RATE    : RES.VOL.  :CON.PR. :BLK.PR. :
:        :           : OIL  : WAT  : GAS  :  STB/DAY  :  STB/DAY  : MSCF/DAY  :  RB/DAY   :  PSIA  :  PSIA  :
=============================================================================================================
:FIELD   :           :      :      :      :        1.0:        2.0:        3.0:        4.0:        :        :
:G1      :           :      :      :      :        5.0:        6.0:        7.0:        8.0:        :        :
:INJ     :  1,  1    :      :      :  GRAT:        9.0:       10.0:       11.0:       12.0:    13.0:    14.0:
:--------:-----------:------:------:------:-----------:-----------:-----------:-----------:--------:--------:
)"
    };

    st.set("FOIR", 1.0);
    st.set("FWIR", 2.0);
    st.set("FGIR", 3.0);
    st.set("FVIR", 4.0);

    st.update_group_var("G1", "GOIR", 5.0);
    st.update_group_var("G1", "GWIR", 6.0);
    st.update_group_var("G1", "GGIR", 7.0);
    st.update_group_var("G1", "GVIR", 8.0);

    st.update_well_var("INJ", "WOIR",  9.0);
    st.update_well_var("INJ", "WWIR", 10.0);
    st.update_well_var("INJ", "WGIR", 11.0);
    st.update_well_var("INJ", "WVIR", 12.0);
    st.update_well_var("INJ", "WBHP", 13.0);
    st.update_well_var("INJ", "WTHP", 14.0);

    Opm::LogOutputHelper<double> helper(eclState, schedule, st, "dummy version");
    helper.injection(0, {});

    const auto data = trimStream(str);
    BOOST_CHECK_EQUAL(data, reference);
}

BOOST_FIXTURE_TEST_CASE(InjectionW2, LogNoteFixture)
{
    const auto reference = std::string {
        R"(============================================= INJECTION REPORT ==============================================
:  WELL  : LOCATION  : CTRL : CTRL : CTRL :    OIL    :   WATER   :    GAS    :   FLUID   : BHP OR : THP OR :
:  NAME  :  (I,J,K)  : MODE : MODE : MODE :   RATE    :   RATE    :   RATE    : RES.VOL.  :CON.PR. :BLK.PR. :
:        :           : OIL  : WAT  : GAS  :  STB/DAY  :  STB/DAY  : MSCF/DAY  :  RB/DAY   :  PSIA  :  PSIA  :
=============================================================================================================
:FIELD   :           :      :      :      :        1.0:        2.0:        3.0:        4.0:        :        :
:G1      :           :      :      :      :        5.0:        6.0:        7.0:        8.0:        :        :
:INJ     :  1,  1    :      :      :  GRAT:        9.0:       10.0:       11.0:       12.0:    13.0:    14.0:
:  BLOCK :  1,  1,  1:      :      :      :       15.0:       16.0:       17.0:       18.0:    19.0:    20.0:
:--------:-----------:------:------:------:-----------:-----------:-----------:-----------:--------:--------:
)"
    };

    st.set("FOIR", 1.0);
    st.set("FWIR", 2.0);
    st.set("FGIR", 3.0);
    st.set("FVIR", 4.0);

    st.update_group_var("G1", "GOIR", 5.0);
    st.update_group_var("G1", "GWIR", 6.0);
    st.update_group_var("G1", "GGIR", 7.0);
    st.update_group_var("G1", "GVIR", 8.0);

    st.update_well_var("INJ", "WOIR",  9.0);
    st.update_well_var("INJ", "WWIR", 10.0);
    st.update_well_var("INJ", "WGIR", 11.0);
    st.update_well_var("INJ", "WVIR", 12.0);
    st.update_well_var("INJ", "WBHP", 13.0);
    st.update_well_var("INJ", "WTHP", 14.0);

    st.update_conn_var("INJ", "COIR", 1, 15.0);
    st.update_conn_var("INJ", "CWIR", 1, 16.0);
    st.update_conn_var("INJ", "CGIR", 1, 17.0);
    st.update_conn_var("INJ", "CVIR", 1, 18.0);
    st.update_conn_var("INJ", "CPR", 1, 19.0);

    const auto bprs = std::map<std::pair<std::string,int>,double> {
        {{"BPR", 1}, eclState.getUnits().to_si(Opm::UnitSystem::measure::pressure, 20.0)},
    };

    Opm::LogOutputHelper<double> helper(eclState, schedule, st, "dummy version");
    helper.injection(0, bprs);

    const auto data = trimStream(str);
    BOOST_CHECK_EQUAL(data, reference);
}

BOOST_FIXTURE_TEST_CASE(MultiSegment, LogNoteFixtureMSW)
{
    const auto reference = std::string {
        R"(================================================ MULTI-SEGMENT WELL REPORT =================================================
:   WELL   : BRN : SEG :    OIL    :   WATER   :    GAS    : MIXTURE :HOLDUP FRACTION :PRESSURE :   PRESSURE HEAD LOSSES   :
:   NAME   : NO. : NO. :   FLOW    :   FLOW    :   FLOW    :VELOCITY : OIL  WAT  GAS  :         :H-STATIC:FRICTION:ACCELRTN:
:          :     :     :  STB/DAY  :  STB/DAY  : MSCF/DAY  : FT/SEC  :                :  PSIA   :  PSI   :  PSI   :  PSI   :
============================================================================================================================
: PROD     :  1  :  1  :        1.0:        2.0:        3.0:   12.200: 0.40 0.50 0.60 :     10.0:  11.000:  12.000:  13.000:
:          :     :  2  :       14.0:       15.0:       16.0:   11.360: 0.17 0.18 0.19 :     23.0:  24.000:  25.000:  27.000:
:----------:-----:-----:-----------:-----------:-----------:---------:----------------:---------:--------:--------:--------:
)"
    };

    st.update_segment_var("PROD", "SOFR", 1, 1.0);
    st.update_segment_var("PROD", "SWFR", 1, 2.0);
    st.update_segment_var("PROD", "SGFR", 1, 3.0);
    st.update_segment_var("PROD", "SOHF", 1, 0.4);
    st.update_segment_var("PROD", "SWHF", 1, 0.5);
    st.update_segment_var("PROD", "SGHF", 1, 0.6);
    st.update_segment_var("PROD", "SOFV", 1, 7.0);
    st.update_segment_var("PROD", "SWFV", 1, 8.0);
    st.update_segment_var("PROD", "SGFV", 1, 9.0);
    st.update_segment_var("PROD", "SPR", 1, 10.0);
    st.update_segment_var("PROD", "SPRDH", 1, 11.0);
    st.update_segment_var("PROD", "SPRDF", 1, 12.0);
    st.update_segment_var("PROD", "SPRDA", 1, 13.0);

    st.update_segment_var("PROD", "SOFR", 2, 14.0);
    st.update_segment_var("PROD", "SWFR", 2, 15.0);
    st.update_segment_var("PROD", "SGFR", 2, 16.0);
    st.update_segment_var("PROD", "SOHF", 2, 0.17);
    st.update_segment_var("PROD", "SWHF", 2, 0.18);
    st.update_segment_var("PROD", "SGHF", 2, 0.19);
    st.update_segment_var("PROD", "SOFV", 2, 20.0);
    st.update_segment_var("PROD", "SWFV", 2, 21.0);
    st.update_segment_var("PROD", "SGFV", 2, 22.0);
    st.update_segment_var("PROD", "SPR", 2, 23.0);
    st.update_segment_var("PROD", "SPRDH", 2, 24.0);
    st.update_segment_var("PROD", "SPRDF", 2, 25.0);
    st.update_segment_var("PROD", "SPRDA", 2, 26.0);
    st.update_segment_var("PROD", "SPRDA", 2, 27.0);

    Opm::LogOutputHelper<double> helper(eclState, schedule, st, "dummy version");
    helper.msw(0);
    const auto data = trimStream(str);
    BOOST_CHECK_EQUAL(data, reference);
}


BOOST_FIXTURE_TEST_CASE(Production, LogNoteFixture)
{
    const auto reference = std::string {
        R"(======================================================= PRODUCTION REPORT =======================================================
:  WELL  : LOCATION  :CTRL:    OIL    :   WATER   :    GAS    :   FLUID   :   WATER   : GAS/OIL  :  WAT/GAS   : BHP OR : THP OR :
:  NAME  :  (I,J,K)  :MODE:   RATE    :   RATE    :   RATE    : RES.VOL.  :    CUT    :  RATIO   :   RATIO    :CON.PR. :BLK.PR. :
:        :           :    :  STB/DAY  :  STB/DAY  : MSCF/DAY  :  RB/DAY   :           : MSCF/STB :  STB/MSCF  :  PSIA  :  PSIA  :
=================================================================================================================================
:FIELD   :           :    :        1.0:        2.0:        3.0:        4.0:      5.000:      6.00:      0.6667:        :        :
:G1      :           :    :        7.0:        8.0:        9.0:       10.0:     11.000:     12.00:      0.8889:        :        :
:PROD    : 10, 10    :ORAT:       13.0:       14.0:       15.0:       16.0:     17.000:     18.00:      0.9333:    19.0:    20.0:
:--------:-----------:----:-----------:-----------:-----------:-----------:-----------:----------:------------:--------:--------:
)"
    };

    st.set("FOPR", 1.0);
    st.set("FWPR", 2.0);
    st.set("FGPR", 3.0);
    st.set("FVPR", 4.0);
    st.set("FWCT", 5.0);
    st.set("FGOR", 6.0);

    st.update_group_var("G1", "GOPR",  7.0);
    st.update_group_var("G1", "GWPR",  8.0);
    st.update_group_var("G1", "GGPR",  9.0);
    st.update_group_var("G1", "GVPR", 10.0);
    st.update_group_var("G1", "GWCT", 11.0);
    st.update_group_var("G1", "GGOR", 12.0);

    st.update_well_var("PROD", "WOPR", 13.0);
    st.update_well_var("PROD", "WWPR", 14.0);
    st.update_well_var("PROD", "WGPR", 15.0);
    st.update_well_var("PROD", "WVPR", 16.0);
    st.update_well_var("PROD", "WWCT", 17.0);
    st.update_well_var("PROD", "WGOR", 18.0);
    st.update_well_var("PROD", "WBHP", 19.0);
    st.update_well_var("PROD", "WTHP", 20.0);

    Opm::LogOutputHelper<double> helper(eclState, schedule, st, "dummy version");
    helper.production(0, {});

    const auto data = trimStream(str);
    BOOST_CHECK_EQUAL(data, reference);
}

BOOST_FIXTURE_TEST_CASE(ProductionW2, LogNoteFixture)
{
    const auto reference = std::string {
        R"(======================================================= PRODUCTION REPORT =======================================================
:  WELL  : LOCATION  :CTRL:    OIL    :   WATER   :    GAS    :   FLUID   :   WATER   : GAS/OIL  :  WAT/GAS   : BHP OR : THP OR :
:  NAME  :  (I,J,K)  :MODE:   RATE    :   RATE    :   RATE    : RES.VOL.  :    CUT    :  RATIO   :   RATIO    :CON.PR. :BLK.PR. :
:        :           :    :  STB/DAY  :  STB/DAY  : MSCF/DAY  :  RB/DAY   :           : MSCF/STB :  STB/MSCF  :  PSIA  :  PSIA  :
=================================================================================================================================
:FIELD   :           :    :        1.0:        2.0:        3.0:        4.0:      5.000:      6.00:      0.6667:        :        :
:G1      :           :    :        7.0:        8.0:        9.0:       10.0:     11.000:     12.00:      0.8889:        :        :
:PROD    : 10, 10    :ORAT:       13.0:       14.0:       15.0:       16.0:     17.000:     18.00:      0.9333:    19.0:    20.0:
:  BLOCK :  2,  2,  1:    :       21.0:       22.0:       23.0:       24.0:     25.000:      0.91:      0.9565:    26.0:    27.0:
:--------:-----------:----:-----------:-----------:-----------:-----------:-----------:----------:------------:--------:--------:
)"
    };

    st.set("FOPR", 1.0);
    st.set("FWPR", 2.0);
    st.set("FGPR", 3.0);
    st.set("FVPR", 4.0);
    st.set("FWCT", 5.0);
    st.set("FGOR", 6.0);

    st.update_group_var("G1", "GOPR",  7.0);
    st.update_group_var("G1", "GWPR",  8.0);
    st.update_group_var("G1", "GGPR",  9.0);
    st.update_group_var("G1", "GVPR", 10.0);
    st.update_group_var("G1", "GWCT", 11.0);
    st.update_group_var("G1", "GGOR", 12.0);

    st.update_well_var("PROD", "WOPR", 13.0);
    st.update_well_var("PROD", "WWPR", 14.0);
    st.update_well_var("PROD", "WGPR", 15.0);
    st.update_well_var("PROD", "WVPR", 16.0);
    st.update_well_var("PROD", "WWCT", 17.0);
    st.update_well_var("PROD", "WGOR", 18.0);
    st.update_well_var("PROD", "WBHP", 19.0);
    st.update_well_var("PROD", "WTHP", 20.0);

    st.update_conn_var("PROD", "COPR", 12, 21.0);
    st.update_conn_var("PROD", "CWPR", 12, 22.0);
    st.update_conn_var("PROD", "CGPR", 12, 23.0);
    st.update_conn_var("PROD", "CVPR", 12, 24.0);
    st.update_conn_var("PROD", "CWCT", 12, 25.0);
    st.update_conn_var("PROD", "CGOR", 12, 21.0 / 23.0);
    st.update_conn_var("PROD", "CPR", 12, 26.0);

    const auto bprs = std::map<std::pair<std::string,int>,double> {
        {{"BPR", 12}, eclState.getUnits().to_si(Opm::UnitSystem::measure::pressure, 27.0)},
    };

    Opm::LogOutputHelper<double> helper(eclState, schedule, st, "dummy version");
    helper.production(0, bprs);

    const auto data = trimStream(str);
    BOOST_CHECK_EQUAL(data, reference);
}
