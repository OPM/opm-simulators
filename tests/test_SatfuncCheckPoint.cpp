/*
  Copyright 2024 Equinor AS

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

#define BOOST_TEST_MODULE TestSatfuncCheckPoint

#ifndef HAVE_MPI
// Suppress GCC diagnostics of the form
//
//   warning: "HAVE_MPI" is not defined, evaluates to 0
//
// when compiling with "-Wundef".
#define HAVE_MPI 0
#endif  // HAVE_MPI

#include <boost/test/unit_test.hpp>

#include <opm/simulators/utils/satfunc/ScaledSatfuncCheckPoint.hpp>
#include <opm/simulators/utils/satfunc/UnscaledSatfuncCheckPoint.hpp>

#include <opm/input/eclipse/EclipseState/Grid/SatfuncPropertyInitializers.hpp>
#include <opm/input/eclipse/EclipseState/EclipseState.hpp>

#include <opm/material/fluidmatrixinteractions/EclEpsGridProperties.hpp>
#include <opm/material/fluidmatrixinteractions/EclEpsScalingPoints.hpp>

#include <opm/input/eclipse/Deck/Deck.hpp>

#include <opm/input/eclipse/Parser/Parser.hpp>

#include <cstddef>
#include <initializer_list>
#include <vector>

// ###########################################################################

namespace {
    Opm::satfunc::RawTableEndPoints exampleRTep()
    {
        auto rtep = Opm::satfunc::RawTableEndPoints{};

        rtep.connate.gas  .assign({ 0.0, 0.0, 0.0, 0.0 , });
        rtep.connate.water.assign({ 0.2, 0.1, 0.0, 0.25, });

        rtep.critical.oil_in_gas  .assign({ 0.1 , 0.125, 0.02, 0.35, });
        rtep.critical.oil_in_water.assign({ 0.2 , 0.225, 0.07, 0.5 , });
        rtep.critical.gas         .assign({ 0.15, 0.257, 0.12, 0.04, });
        rtep.critical.water       .assign({ 0.27, 0.325, 0.25, 0.42, });

        rtep.maximum.gas  .assign({ 0.8, 0.9, 1.0, 0.75, });
        rtep.maximum.water.assign({ 1.0, 1.0, 1.0, 0.8 , });

        return rtep;
    }

    Opm::satfunc::RawFunctionValues exampleRFunc()
    {
        auto rfunc = Opm::satfunc::RawFunctionValues{};

        rfunc.kro.max.assign({ 0.75, 0.8, 0.85, 0.7, });
        rfunc.kro.rg .assign({ 0.62, 0.6, 0.78, 0.7, });
        rfunc.kro.rw .assign({ 0.53, 0.4, 0.82, 0.5, });

        rfunc.krg.max.assign({ 1.0, 0.99, 0.98, 0.95, });
        rfunc.krg.r  .assign({ 0.8, 0.79, 0.87, 0.75, });

        rfunc.krw.max.assign({ 0.7, 0.69, 0.68, 0.65, });
        rfunc.krw.r  .assign({ 0.5, 0.49, 0.48, 0.45, });

        rfunc.pc.g.assign({  0.0, 1.0, 2.0, 3.0, });
        rfunc.pc.w.assign({ 10.0, 8.5, 7.0, 5.5, });

        return rfunc;
    }
} // Anonymous namespace

// ===========================================================================

BOOST_AUTO_TEST_SUITE(Unscaled_Curve)

BOOST_AUTO_TEST_CASE(All_Regions)
{
    using CurvePt = Opm::Satfunc::PhaseChecks::UnscaledSatfuncCheckPoint<double>;

    const auto satnum = std::vector { 1, 2, 3, 4, };
    const auto rtep   = exampleRTep();
    const auto rfunc  = exampleRFunc();

    const auto curvePt = CurvePt {
        &satnum, 1, CurvePt::UnscaledEndPoints { &rtep, &rfunc }
    };

    {
        const auto id = curvePt.pointID(0);
        BOOST_REQUIRE_MESSAGE(id.has_value(),
                              "pointID(0) must return a value on first call");

        BOOST_CHECK_EQUAL(*id, std::size_t{1});

        BOOST_CHECK_MESSAGE(! curvePt.pointID(0).has_value(),
                            "pointID(0) must NOT return a value on second call");
    }

    {
        const auto id = curvePt.pointID(1);
        BOOST_REQUIRE_MESSAGE(id.has_value(),
                              "pointID(1) must return a value on first call");

        BOOST_CHECK_EQUAL(*id, std::size_t{2});

        BOOST_CHECK_MESSAGE(! curvePt.pointID(1).has_value(),
                            "pointID(1) must NOT return a value on second call");
    }

    {
        const auto id = curvePt.pointID(2);
        BOOST_REQUIRE_MESSAGE(id.has_value(),
                              "pointID(2) must return a value on first call");

        BOOST_CHECK_EQUAL(*id, std::size_t{3});

        BOOST_CHECK_MESSAGE(! curvePt.pointID(0).has_value(),
                            "pointID(2) must NOT return a value on second call");
    }

    {
        const auto id = curvePt.pointID(3);
        BOOST_REQUIRE_MESSAGE(id.has_value(),
                              "pointID(3) must return a value on first call");

        BOOST_CHECK_EQUAL(*id, std::size_t{4});

        BOOST_CHECK_MESSAGE(! curvePt.pointID(0).has_value(),
                            "pointID(3) must NOT return a value on second call");
    }

    // Unscaled saturation function end-points in cell 0 (SATNUM = 1).
    {
        auto endPoints = Opm::EclEpsScalingPointsInfo<double>{};
        curvePt.populateCheckPoint(0, endPoints);

        BOOST_CHECK_CLOSE(endPoints.Swl, 0.2, 1.0e-8);
        BOOST_CHECK_CLOSE(endPoints.Sgl, 0.0, 1.0e-8);

        BOOST_CHECK_CLOSE(endPoints.Swcr,  0.27, 1.0e-8);
        BOOST_CHECK_CLOSE(endPoints.Sgcr,  0.15, 1.0e-8);
        BOOST_CHECK_CLOSE(endPoints.Sowcr, 0.2 , 1.0e-8);
        BOOST_CHECK_CLOSE(endPoints.Sogcr, 0.1 , 1.0e-8);

        BOOST_CHECK_CLOSE(endPoints.Swu, 1.0, 1.0e-8);
        BOOST_CHECK_CLOSE(endPoints.Sgu, 0.8, 1.0e-8);

        BOOST_CHECK_CLOSE(endPoints.maxPcow, 10.0, 1.0e-8);
        BOOST_CHECK_CLOSE(endPoints.maxPcgo,  0.0, 1.0e-8);

        BOOST_CHECK_CLOSE(endPoints.pcowLeverettFactor, 1.0, 1.0e-8);
        BOOST_CHECK_CLOSE(endPoints.pcgoLeverettFactor, 1.0, 1.0e-8);

        BOOST_CHECK_CLOSE(endPoints.Krwr,  0.5,  1.0e-8);
        BOOST_CHECK_CLOSE(endPoints.Krgr,  0.8,  1.0e-8);
        BOOST_CHECK_CLOSE(endPoints.Krorw, 0.53, 1.0e-8);
        BOOST_CHECK_CLOSE(endPoints.Krorg, 0.62, 1.0e-8);

        BOOST_CHECK_CLOSE(endPoints.maxKrw,  0.7,  1.0e-8);
        BOOST_CHECK_CLOSE(endPoints.maxKrow, 0.75, 1.0e-8);
        BOOST_CHECK_CLOSE(endPoints.maxKrog, 0.75, 1.0e-8);
        BOOST_CHECK_CLOSE(endPoints.maxKrg,  1.0,  1.0e-8);
    }

    // Unscaled saturation function end-points in cell 0.  Intentionally
    // same as previous.  Member function populateCheckPoint() should not
    // care that we've called the function with the same arguments before.
    {
        auto endPoints = Opm::EclEpsScalingPointsInfo<double>{};
        curvePt.populateCheckPoint(0, endPoints);

        BOOST_CHECK_CLOSE(endPoints.Swl, 0.2, 1.0e-8);
        BOOST_CHECK_CLOSE(endPoints.Sgl, 0.0, 1.0e-8);

        BOOST_CHECK_CLOSE(endPoints.Swcr,  0.27, 1.0e-8);
        BOOST_CHECK_CLOSE(endPoints.Sgcr,  0.15, 1.0e-8);
        BOOST_CHECK_CLOSE(endPoints.Sowcr, 0.2 , 1.0e-8);
        BOOST_CHECK_CLOSE(endPoints.Sogcr, 0.1 , 1.0e-8);

        BOOST_CHECK_CLOSE(endPoints.Swu, 1.0, 1.0e-8);
        BOOST_CHECK_CLOSE(endPoints.Sgu, 0.8, 1.0e-8);

        BOOST_CHECK_CLOSE(endPoints.maxPcow, 10.0, 1.0e-8);
        BOOST_CHECK_CLOSE(endPoints.maxPcgo,  0.0, 1.0e-8);

        BOOST_CHECK_CLOSE(endPoints.pcowLeverettFactor, 1.0, 1.0e-8);
        BOOST_CHECK_CLOSE(endPoints.pcgoLeverettFactor, 1.0, 1.0e-8);

        BOOST_CHECK_CLOSE(endPoints.Krwr,  0.5,  1.0e-8);
        BOOST_CHECK_CLOSE(endPoints.Krgr,  0.8,  1.0e-8);
        BOOST_CHECK_CLOSE(endPoints.Krorw, 0.53, 1.0e-8);
        BOOST_CHECK_CLOSE(endPoints.Krorg, 0.62, 1.0e-8);

        BOOST_CHECK_CLOSE(endPoints.maxKrw,  0.7,  1.0e-8);
        BOOST_CHECK_CLOSE(endPoints.maxKrow, 0.75, 1.0e-8);
        BOOST_CHECK_CLOSE(endPoints.maxKrog, 0.75, 1.0e-8);
        BOOST_CHECK_CLOSE(endPoints.maxKrg,  1.0,  1.0e-8);
    }

    // Unscaled saturation function end-points in cell 1 (SATNUM = 2).
    {
        auto endPoints = Opm::EclEpsScalingPointsInfo<double>{};
        curvePt.populateCheckPoint(1, endPoints);

        BOOST_CHECK_CLOSE(endPoints.Swl, 0.1, 1.0e-8);
        BOOST_CHECK_CLOSE(endPoints.Sgl, 0.0, 1.0e-8);

        BOOST_CHECK_CLOSE(endPoints.Swcr,  0.325, 1.0e-8);
        BOOST_CHECK_CLOSE(endPoints.Sgcr,  0.257, 1.0e-8);
        BOOST_CHECK_CLOSE(endPoints.Sowcr, 0.225, 1.0e-8);
        BOOST_CHECK_CLOSE(endPoints.Sogcr, 0.125, 1.0e-8);

        BOOST_CHECK_CLOSE(endPoints.Swu, 1.0, 1.0e-8);
        BOOST_CHECK_CLOSE(endPoints.Sgu, 0.9, 1.0e-8);

        BOOST_CHECK_CLOSE(endPoints.maxPcow, 8.5, 1.0e-8);
        BOOST_CHECK_CLOSE(endPoints.maxPcgo, 1.0, 1.0e-8);

        BOOST_CHECK_CLOSE(endPoints.pcowLeverettFactor, 1.0, 1.0e-8);
        BOOST_CHECK_CLOSE(endPoints.pcgoLeverettFactor, 1.0, 1.0e-8);

        BOOST_CHECK_CLOSE(endPoints.Krwr,  0.49, 1.0e-8);
        BOOST_CHECK_CLOSE(endPoints.Krgr,  0.79, 1.0e-8);
        BOOST_CHECK_CLOSE(endPoints.Krorw, 0.4,  1.0e-8);
        BOOST_CHECK_CLOSE(endPoints.Krorg, 0.6,  1.0e-8);

        BOOST_CHECK_CLOSE(endPoints.maxKrw,  0.69, 1.0e-8);
        BOOST_CHECK_CLOSE(endPoints.maxKrow, 0.8,  1.0e-8);
        BOOST_CHECK_CLOSE(endPoints.maxKrog, 0.8,  1.0e-8);
        BOOST_CHECK_CLOSE(endPoints.maxKrg,  0.99, 1.0e-8);
    }

    // Unscaled saturation function end-points in cell 2 (SATNUM = 3).
    {
        auto endPoints = Opm::EclEpsScalingPointsInfo<double>{};
        curvePt.populateCheckPoint(2, endPoints);

        BOOST_CHECK_CLOSE(endPoints.Swl, 0.0, 1.0e-8);
        BOOST_CHECK_CLOSE(endPoints.Sgl, 0.0, 1.0e-8);

        BOOST_CHECK_CLOSE(endPoints.Swcr,  0.25, 1.0e-8);
        BOOST_CHECK_CLOSE(endPoints.Sgcr,  0.12, 1.0e-8);
        BOOST_CHECK_CLOSE(endPoints.Sowcr, 0.07, 1.0e-8);
        BOOST_CHECK_CLOSE(endPoints.Sogcr, 0.02, 1.0e-8);

        BOOST_CHECK_CLOSE(endPoints.Swu, 1.0, 1.0e-8);
        BOOST_CHECK_CLOSE(endPoints.Sgu, 1.0, 1.0e-8);

        BOOST_CHECK_CLOSE(endPoints.maxPcow, 7.0, 1.0e-8);
        BOOST_CHECK_CLOSE(endPoints.maxPcgo, 2.0, 1.0e-8);

        BOOST_CHECK_CLOSE(endPoints.pcowLeverettFactor, 1.0, 1.0e-8);
        BOOST_CHECK_CLOSE(endPoints.pcgoLeverettFactor, 1.0, 1.0e-8);

        BOOST_CHECK_CLOSE(endPoints.Krwr,  0.48, 1.0e-8);
        BOOST_CHECK_CLOSE(endPoints.Krgr,  0.87, 1.0e-8);
        BOOST_CHECK_CLOSE(endPoints.Krorw, 0.82, 1.0e-8);
        BOOST_CHECK_CLOSE(endPoints.Krorg, 0.78, 1.0e-8);

        BOOST_CHECK_CLOSE(endPoints.maxKrw,  0.68, 1.0e-8);
        BOOST_CHECK_CLOSE(endPoints.maxKrow, 0.85, 1.0e-8);
        BOOST_CHECK_CLOSE(endPoints.maxKrog, 0.85, 1.0e-8);
        BOOST_CHECK_CLOSE(endPoints.maxKrg,  0.98, 1.0e-8);
    }

    // Unscaled saturation function end-points in cell 3 (SATNUM = 4).
    {
        auto endPoints = Opm::EclEpsScalingPointsInfo<double>{};
        curvePt.populateCheckPoint(3, endPoints);

        BOOST_CHECK_CLOSE(endPoints.Swl, 0.25, 1.0e-8);
        BOOST_CHECK_CLOSE(endPoints.Sgl, 0.0,  1.0e-8);

        BOOST_CHECK_CLOSE(endPoints.Swcr,  0.42, 1.0e-8);
        BOOST_CHECK_CLOSE(endPoints.Sgcr,  0.04, 1.0e-8);
        BOOST_CHECK_CLOSE(endPoints.Sowcr, 0.5,  1.0e-8);
        BOOST_CHECK_CLOSE(endPoints.Sogcr, 0.35, 1.0e-8);

        BOOST_CHECK_CLOSE(endPoints.Swu, 0.8,  1.0e-8);
        BOOST_CHECK_CLOSE(endPoints.Sgu, 0.75, 1.0e-8);

        BOOST_CHECK_CLOSE(endPoints.maxPcow, 5.5, 1.0e-8);
        BOOST_CHECK_CLOSE(endPoints.maxPcgo, 3.0, 1.0e-8);

        BOOST_CHECK_CLOSE(endPoints.pcowLeverettFactor, 1.0, 1.0e-8);
        BOOST_CHECK_CLOSE(endPoints.pcgoLeverettFactor, 1.0, 1.0e-8);

        BOOST_CHECK_CLOSE(endPoints.Krwr,  0.45, 1.0e-8);
        BOOST_CHECK_CLOSE(endPoints.Krgr,  0.75, 1.0e-8);
        BOOST_CHECK_CLOSE(endPoints.Krorw, 0.5,  1.0e-8);
        BOOST_CHECK_CLOSE(endPoints.Krorg, 0.7,  1.0e-8);

        BOOST_CHECK_CLOSE(endPoints.maxKrw,  0.65, 1.0e-8);
        BOOST_CHECK_CLOSE(endPoints.maxKrow, 0.7,  1.0e-8);
        BOOST_CHECK_CLOSE(endPoints.maxKrog, 0.7,  1.0e-8);
        BOOST_CHECK_CLOSE(endPoints.maxKrg,  0.95, 1.0e-8);
    }
}

BOOST_AUTO_TEST_CASE(Same_Region)
{
    using CurvePt = Opm::Satfunc::PhaseChecks::UnscaledSatfuncCheckPoint<double>;

    const auto satnum = std::vector { 2, 2, 2, 2, };
    const auto rtep   = exampleRTep();
    const auto rfunc  = exampleRFunc();

    const auto curvePt = CurvePt {
        &satnum, 1, CurvePt::UnscaledEndPoints { &rtep, &rfunc }
    };

    {
        const auto id = curvePt.pointID(0);
        BOOST_REQUIRE_MESSAGE(id.has_value(),
                              "pointID(0) must return a value on first call");

        BOOST_CHECK_EQUAL(*id, std::size_t{2});

        BOOST_CHECK_MESSAGE(! curvePt.pointID(0).has_value(),
                            "pointID(0) must NOT return a value on second call");
    }

    BOOST_CHECK_MESSAGE(! curvePt.pointID(1).has_value(),
                        "pointID(1) must NOT return a value on first call");

    BOOST_CHECK_MESSAGE(! curvePt.pointID(2).has_value(),
                        "pointID(2) must NOT return a value on first call");

    BOOST_CHECK_MESSAGE(! curvePt.pointID(3).has_value(),
                        "pointID(3) must NOT return a value on first call");

    // Unscaled saturation function end-points in cell 0 (SATNUM = 2).
    {
        auto endPoints = Opm::EclEpsScalingPointsInfo<double>{};
        curvePt.populateCheckPoint(0, endPoints);

        BOOST_CHECK_CLOSE(endPoints.Swl, 0.1, 1.0e-8);
        BOOST_CHECK_CLOSE(endPoints.Sgl, 0.0, 1.0e-8);

        BOOST_CHECK_CLOSE(endPoints.Swcr,  0.325, 1.0e-8);
        BOOST_CHECK_CLOSE(endPoints.Sgcr,  0.257, 1.0e-8);
        BOOST_CHECK_CLOSE(endPoints.Sowcr, 0.225, 1.0e-8);
        BOOST_CHECK_CLOSE(endPoints.Sogcr, 0.125, 1.0e-8);

        BOOST_CHECK_CLOSE(endPoints.Swu, 1.0, 1.0e-8);
        BOOST_CHECK_CLOSE(endPoints.Sgu, 0.9, 1.0e-8);

        BOOST_CHECK_CLOSE(endPoints.maxPcow, 8.5, 1.0e-8);
        BOOST_CHECK_CLOSE(endPoints.maxPcgo, 1.0, 1.0e-8);

        BOOST_CHECK_CLOSE(endPoints.pcowLeverettFactor, 1.0, 1.0e-8);
        BOOST_CHECK_CLOSE(endPoints.pcgoLeverettFactor, 1.0, 1.0e-8);

        BOOST_CHECK_CLOSE(endPoints.Krwr,  0.49, 1.0e-8);
        BOOST_CHECK_CLOSE(endPoints.Krgr,  0.79, 1.0e-8);
        BOOST_CHECK_CLOSE(endPoints.Krorw, 0.4,  1.0e-8);
        BOOST_CHECK_CLOSE(endPoints.Krorg, 0.6,  1.0e-8);

        BOOST_CHECK_CLOSE(endPoints.maxKrw,  0.69, 1.0e-8);
        BOOST_CHECK_CLOSE(endPoints.maxKrow, 0.8,  1.0e-8);
        BOOST_CHECK_CLOSE(endPoints.maxKrog, 0.8,  1.0e-8);
        BOOST_CHECK_CLOSE(endPoints.maxKrg,  0.99, 1.0e-8);
    }

    // Unscaled saturation function end-points in cell 1 (SATNUM = 2).
    {
        auto endPoints = Opm::EclEpsScalingPointsInfo<double>{};
        curvePt.populateCheckPoint(1, endPoints);

        BOOST_CHECK_CLOSE(endPoints.Swl, 0.1, 1.0e-8);
        BOOST_CHECK_CLOSE(endPoints.Sgl, 0.0, 1.0e-8);

        BOOST_CHECK_CLOSE(endPoints.Swcr,  0.325, 1.0e-8);
        BOOST_CHECK_CLOSE(endPoints.Sgcr,  0.257, 1.0e-8);
        BOOST_CHECK_CLOSE(endPoints.Sowcr, 0.225, 1.0e-8);
        BOOST_CHECK_CLOSE(endPoints.Sogcr, 0.125, 1.0e-8);

        BOOST_CHECK_CLOSE(endPoints.Swu, 1.0, 1.0e-8);
        BOOST_CHECK_CLOSE(endPoints.Sgu, 0.9, 1.0e-8);

        BOOST_CHECK_CLOSE(endPoints.maxPcow, 8.5, 1.0e-8);
        BOOST_CHECK_CLOSE(endPoints.maxPcgo, 1.0, 1.0e-8);

        BOOST_CHECK_CLOSE(endPoints.pcowLeverettFactor, 1.0, 1.0e-8);
        BOOST_CHECK_CLOSE(endPoints.pcgoLeverettFactor, 1.0, 1.0e-8);

        BOOST_CHECK_CLOSE(endPoints.Krwr,  0.49, 1.0e-8);
        BOOST_CHECK_CLOSE(endPoints.Krgr,  0.79, 1.0e-8);
        BOOST_CHECK_CLOSE(endPoints.Krorw, 0.4,  1.0e-8);
        BOOST_CHECK_CLOSE(endPoints.Krorg, 0.6,  1.0e-8);

        BOOST_CHECK_CLOSE(endPoints.maxKrw,  0.69, 1.0e-8);
        BOOST_CHECK_CLOSE(endPoints.maxKrow, 0.8,  1.0e-8);
        BOOST_CHECK_CLOSE(endPoints.maxKrog, 0.8,  1.0e-8);
        BOOST_CHECK_CLOSE(endPoints.maxKrg,  0.99, 1.0e-8);
    }

    // Unscaled saturation function end-points in cell 2 (SATNUM = 2).
    {
        auto endPoints = Opm::EclEpsScalingPointsInfo<double>{};
        curvePt.populateCheckPoint(2, endPoints);

        BOOST_CHECK_CLOSE(endPoints.Swl, 0.1, 1.0e-8);
        BOOST_CHECK_CLOSE(endPoints.Sgl, 0.0, 1.0e-8);

        BOOST_CHECK_CLOSE(endPoints.Swcr,  0.325, 1.0e-8);
        BOOST_CHECK_CLOSE(endPoints.Sgcr,  0.257, 1.0e-8);
        BOOST_CHECK_CLOSE(endPoints.Sowcr, 0.225, 1.0e-8);
        BOOST_CHECK_CLOSE(endPoints.Sogcr, 0.125, 1.0e-8);

        BOOST_CHECK_CLOSE(endPoints.Swu, 1.0, 1.0e-8);
        BOOST_CHECK_CLOSE(endPoints.Sgu, 0.9, 1.0e-8);

        BOOST_CHECK_CLOSE(endPoints.maxPcow, 8.5, 1.0e-8);
        BOOST_CHECK_CLOSE(endPoints.maxPcgo, 1.0, 1.0e-8);

        BOOST_CHECK_CLOSE(endPoints.pcowLeverettFactor, 1.0, 1.0e-8);
        BOOST_CHECK_CLOSE(endPoints.pcgoLeverettFactor, 1.0, 1.0e-8);

        BOOST_CHECK_CLOSE(endPoints.Krwr,  0.49, 1.0e-8);
        BOOST_CHECK_CLOSE(endPoints.Krgr,  0.79, 1.0e-8);
        BOOST_CHECK_CLOSE(endPoints.Krorw, 0.4,  1.0e-8);
        BOOST_CHECK_CLOSE(endPoints.Krorg, 0.6,  1.0e-8);

        BOOST_CHECK_CLOSE(endPoints.maxKrw,  0.69, 1.0e-8);
        BOOST_CHECK_CLOSE(endPoints.maxKrow, 0.8,  1.0e-8);
        BOOST_CHECK_CLOSE(endPoints.maxKrog, 0.8,  1.0e-8);
        BOOST_CHECK_CLOSE(endPoints.maxKrg,  0.99, 1.0e-8);
    }

    // Unscaled saturation function end-points in cell 3 (SATNUM = 2).
    {
        auto endPoints = Opm::EclEpsScalingPointsInfo<double>{};
        curvePt.populateCheckPoint(3, endPoints);

        BOOST_CHECK_CLOSE(endPoints.Swl, 0.1, 1.0e-8);
        BOOST_CHECK_CLOSE(endPoints.Sgl, 0.0, 1.0e-8);

        BOOST_CHECK_CLOSE(endPoints.Swcr,  0.325, 1.0e-8);
        BOOST_CHECK_CLOSE(endPoints.Sgcr,  0.257, 1.0e-8);
        BOOST_CHECK_CLOSE(endPoints.Sowcr, 0.225, 1.0e-8);
        BOOST_CHECK_CLOSE(endPoints.Sogcr, 0.125, 1.0e-8);

        BOOST_CHECK_CLOSE(endPoints.Swu, 1.0, 1.0e-8);
        BOOST_CHECK_CLOSE(endPoints.Sgu, 0.9, 1.0e-8);

        BOOST_CHECK_CLOSE(endPoints.maxPcow, 8.5, 1.0e-8);
        BOOST_CHECK_CLOSE(endPoints.maxPcgo, 1.0, 1.0e-8);

        BOOST_CHECK_CLOSE(endPoints.pcowLeverettFactor, 1.0, 1.0e-8);
        BOOST_CHECK_CLOSE(endPoints.pcgoLeverettFactor, 1.0, 1.0e-8);

        BOOST_CHECK_CLOSE(endPoints.Krwr,  0.49, 1.0e-8);
        BOOST_CHECK_CLOSE(endPoints.Krgr,  0.79, 1.0e-8);
        BOOST_CHECK_CLOSE(endPoints.Krorw, 0.4,  1.0e-8);
        BOOST_CHECK_CLOSE(endPoints.Krorg, 0.6,  1.0e-8);

        BOOST_CHECK_CLOSE(endPoints.maxKrw,  0.69, 1.0e-8);
        BOOST_CHECK_CLOSE(endPoints.maxKrow, 0.8,  1.0e-8);
        BOOST_CHECK_CLOSE(endPoints.maxKrog, 0.8,  1.0e-8);
        BOOST_CHECK_CLOSE(endPoints.maxKrg,  0.99, 1.0e-8);
    }
}

BOOST_AUTO_TEST_SUITE_END() // Unscaled_Curve

// ===========================================================================

BOOST_AUTO_TEST_SUITE(Scaled_Curve)

namespace {
    Opm::EclipseState scaledEndpointProps()
    {
        return Opm::EclipseState {
            Opm::Parser{}.parseString(R"(RUNSPEC
DIMENS
1 5 1 /
ENDSCALE
/
TABDIMS
 3 /
OIL
GAS
WATER
GRID
DXV
  100.0 /
DYV
  5*20.0 /
DZV
  5.0 /
DEPTHZ
  12*2000.0 /
PORO
  5*0.3 /
PROPS
SCALECRS
  'YES' /
SWOF
 0.0 0.0 1.0 0.0
 1.0 1.0 0.0 0.0
/
/
/
SGOF
 0.0 0.0 1.0 0.0
 1.0 1.0 0.0 0.0
/
/
/
SWL
  0.05 0.1 0.15 0.2 0.27 /
SWCR
  0.06 0.11 0.16 0.21 0.3 /
SWU
  0.98 0.97 0.96 0.95 0.94 /
SGL
  0.01 0.02 0.03 0.04 0.05 /
SGCR
  0.1 0.11 0.12 0.13 0.14 /
SGU
  0.89 0.9 0.91 0.92 0.93 /
SOWCR
  0.21 0.23 0.25 0.27 0.29 /
SOGCR
  0.15 0.16 0.17 0.18 0.19 /
KRW
  0.72 0.77 0.82 0.87 0.92 /
KRWR
  0.6 0.63 0.66 0.71 0.65 /
KRG
  0.85 0.84 0.86 0.83 0.87 /
KRGR
  5*0.775 /
KRORW
  0.55 0.65 0.65 0.7 0.6 /
-- No KRORG => use unscaled function values
REGIONS
SATNUM
 1 3 2 2 1 /
END
)")
        };
    }
} // Anonymous namespace

BOOST_AUTO_TEST_CASE(All_Cells)
{
    using UnscaledCurvePt = Opm::Satfunc::PhaseChecks::UnscaledSatfuncCheckPoint<float>;
    using ScaledCurvePt = Opm::Satfunc::PhaseChecks::ScaledSatfuncCheckPoint<float>;

    const auto es = scaledEndpointProps();

    const auto* satnum = &es.fieldProps().get_int("SATNUM"); // {1, 3, 2, 2, 1}
    const auto rtep  = exampleRTep();
    const auto rfunc = exampleRFunc();

    const auto unscaledCurvePt = UnscaledCurvePt {
        satnum, 1, UnscaledCurvePt::UnscaledEndPoints { &rtep, &rfunc }
    };

    const auto epsGridProps = Opm::EclEpsGridProperties { es, /* useImbibition = */ false };

    const auto scaledCurvePt = ScaledCurvePt {
        unscaledCurvePt, &es, &epsGridProps,
        [](const auto& i) { return static_cast<std::size_t>(i); }
    };

    {
        const auto id = scaledCurvePt.pointID(0);
        BOOST_REQUIRE_MESSAGE(id.has_value(),
                              "pointID(0) must return a value on first call");

        BOOST_CHECK_EQUAL(*id, std::size_t{0});

        BOOST_CHECK_MESSAGE(scaledCurvePt.pointID(0).has_value(),
                            "pointID(0) must return a value on second call");
    }

    // Scaled saturation function end-points in cell 0 (SATNUM = 1).
    {
        auto endPoints = Opm::EclEpsScalingPointsInfo<float>{};
        scaledCurvePt.populateCheckPoint(0, endPoints);

        BOOST_CHECK_CLOSE(endPoints.Swl, 0.05f, 1.0e-6f);
        BOOST_CHECK_CLOSE(endPoints.Sgl, 0.01f, 1.0e-6f);

        BOOST_CHECK_CLOSE(endPoints.Swcr,  0.06f, 1.0e-6f);
        BOOST_CHECK_CLOSE(endPoints.Sgcr,  0.1f,  1.0e-6f);
        BOOST_CHECK_CLOSE(endPoints.Sowcr, 0.21f, 1.0e-6f);
        BOOST_CHECK_CLOSE(endPoints.Sogcr, 0.15f, 1.0e-6f);

        BOOST_CHECK_CLOSE(endPoints.Swu, 0.98f, 1.0e-6f);
        BOOST_CHECK_CLOSE(endPoints.Sgu, 0.89f, 1.0e-6f);

        BOOST_CHECK_CLOSE(endPoints.maxPcow, 10.0f, 1.0e-6f);
        BOOST_CHECK_CLOSE(endPoints.maxPcgo,  0.0f, 1.0e-6f);

        BOOST_CHECK_CLOSE(endPoints.pcowLeverettFactor, 1.0f, 1.0e-6f);
        BOOST_CHECK_CLOSE(endPoints.pcgoLeverettFactor, 1.0f, 1.0e-6f);

        BOOST_CHECK_CLOSE(endPoints.Krwr,  0.6f,   1.0e-6f);
        BOOST_CHECK_CLOSE(endPoints.Krgr,  0.775f, 1.0e-6f);
        BOOST_CHECK_CLOSE(endPoints.Krorw, 0.55f,  1.0e-6f);
        BOOST_CHECK_CLOSE(endPoints.Krorg, 0.62f,  1.0e-6f);

        BOOST_CHECK_CLOSE(endPoints.maxKrw,  0.72f, 1.0e-6f);
        BOOST_CHECK_CLOSE(endPoints.maxKrow, 0.75f, 1.0e-6f);
        BOOST_CHECK_CLOSE(endPoints.maxKrog, 0.75f, 1.0e-6f);
        BOOST_CHECK_CLOSE(endPoints.maxKrg,  0.85f, 1.0e-6f);
    }

    // Scaled saturation function end-points in cell 1 (SATNUM = 3).
    {
        auto endPoints = Opm::EclEpsScalingPointsInfo<float>{};
        scaledCurvePt.populateCheckPoint(1, endPoints);

        BOOST_CHECK_CLOSE(endPoints.Swl, 0.1f, 1.0e-6f);
        BOOST_CHECK_CLOSE(endPoints.Sgl, 0.02f, 1.0e-6f);

        BOOST_CHECK_CLOSE(endPoints.Swcr,  0.11f, 1.0e-6f);
        BOOST_CHECK_CLOSE(endPoints.Sgcr,  0.11f, 1.0e-6f);
        BOOST_CHECK_CLOSE(endPoints.Sowcr, 0.23f, 1.0e-6f);
        BOOST_CHECK_CLOSE(endPoints.Sogcr, 0.16f, 1.0e-6f);

        BOOST_CHECK_CLOSE(endPoints.Swu, 0.97f, 1.0e-6f);
        BOOST_CHECK_CLOSE(endPoints.Sgu, 0.9f,  1.0e-6f);

        BOOST_CHECK_CLOSE(endPoints.maxPcow, 7.0f, 1.0e-6f);
        BOOST_CHECK_CLOSE(endPoints.maxPcgo, 2.0f, 1.0e-6f);

        BOOST_CHECK_CLOSE(endPoints.pcowLeverettFactor, 1.0f, 1.0e-6f);
        BOOST_CHECK_CLOSE(endPoints.pcgoLeverettFactor, 1.0f, 1.0e-6f);

        BOOST_CHECK_CLOSE(endPoints.Krwr,  0.63f,  1.0e-6f);
        BOOST_CHECK_CLOSE(endPoints.Krgr,  0.775f, 1.0e-6f);
        BOOST_CHECK_CLOSE(endPoints.Krorw, 0.65f,  1.0e-6f);
        BOOST_CHECK_CLOSE(endPoints.Krorg, 0.78f,  1.0e-6f);

        BOOST_CHECK_CLOSE(endPoints.maxKrw,  0.77f, 1.0e-6f);
        BOOST_CHECK_CLOSE(endPoints.maxKrow, 0.85f, 1.0e-6f);
        BOOST_CHECK_CLOSE(endPoints.maxKrog, 0.85f, 1.0e-6f);
        BOOST_CHECK_CLOSE(endPoints.maxKrg,  0.84f, 1.0e-6f);
    }

    // Scaled saturation function end-points in cell 2 (SATNUM = 2).
    {
        auto endPoints = Opm::EclEpsScalingPointsInfo<float>{};
        scaledCurvePt.populateCheckPoint(2, endPoints);

        BOOST_CHECK_CLOSE(endPoints.Swl, 0.15f, 1.0e-6f);
        BOOST_CHECK_CLOSE(endPoints.Sgl, 0.03f, 1.0e-6f);

        BOOST_CHECK_CLOSE(endPoints.Swcr,  0.16f, 1.0e-6f);
        BOOST_CHECK_CLOSE(endPoints.Sgcr,  0.12f, 1.0e-6f);
        BOOST_CHECK_CLOSE(endPoints.Sowcr, 0.25f, 1.0e-6f);
        BOOST_CHECK_CLOSE(endPoints.Sogcr, 0.17f, 1.0e-6f);

        BOOST_CHECK_CLOSE(endPoints.Swu, 0.96f, 1.0e-6f);
        BOOST_CHECK_CLOSE(endPoints.Sgu, 0.91f, 1.0e-6f);

        BOOST_CHECK_CLOSE(endPoints.maxPcow, 8.5f, 1.0e-6f);
        BOOST_CHECK_CLOSE(endPoints.maxPcgo, 1.0f, 1.0e-6f);

        BOOST_CHECK_CLOSE(endPoints.pcowLeverettFactor, 1.0f, 1.0e-6f);
        BOOST_CHECK_CLOSE(endPoints.pcgoLeverettFactor, 1.0f, 1.0e-6f);

        BOOST_CHECK_CLOSE(endPoints.Krwr,  0.66f,  1.0e-6f);
        BOOST_CHECK_CLOSE(endPoints.Krgr,  0.775f, 1.0e-6f);
        BOOST_CHECK_CLOSE(endPoints.Krorw, 0.65f,  1.0e-6f);
        BOOST_CHECK_CLOSE(endPoints.Krorg, 0.6f,   1.0e-6f);

        BOOST_CHECK_CLOSE(endPoints.maxKrw,  0.82f, 1.0e-6f);
        BOOST_CHECK_CLOSE(endPoints.maxKrow, 0.8f,  1.0e-6f);
        BOOST_CHECK_CLOSE(endPoints.maxKrog, 0.8f,  1.0e-6f);
        BOOST_CHECK_CLOSE(endPoints.maxKrg,  0.86f, 1.0e-6f);
    }

    // Scaled saturation function end-points in cell 3 (SATNUM = 2).
    {
        auto endPoints = Opm::EclEpsScalingPointsInfo<float>{};
        scaledCurvePt.populateCheckPoint(3, endPoints);

        BOOST_CHECK_CLOSE(endPoints.Swl, 0.2f,  1.0e-6f);
        BOOST_CHECK_CLOSE(endPoints.Sgl, 0.04f, 1.0e-6f);

        BOOST_CHECK_CLOSE(endPoints.Swcr,  0.21f, 1.0e-6f);
        BOOST_CHECK_CLOSE(endPoints.Sgcr,  0.13f, 1.0e-6f);
        BOOST_CHECK_CLOSE(endPoints.Sowcr, 0.27f, 1.0e-6f);
        BOOST_CHECK_CLOSE(endPoints.Sogcr, 0.18f, 1.0e-6f);

        BOOST_CHECK_CLOSE(endPoints.Swu, 0.95f, 1.0e-6f);
        BOOST_CHECK_CLOSE(endPoints.Sgu, 0.92f, 1.0e-6f);

        BOOST_CHECK_CLOSE(endPoints.maxPcow, 8.5f, 1.0e-6f);
        BOOST_CHECK_CLOSE(endPoints.maxPcgo, 1.0f, 1.0e-6f);

        BOOST_CHECK_CLOSE(endPoints.pcowLeverettFactor, 1.0f, 1.0e-6f);
        BOOST_CHECK_CLOSE(endPoints.pcgoLeverettFactor, 1.0f, 1.0e-6f);

        BOOST_CHECK_CLOSE(endPoints.Krwr,  0.71f,  1.0e-6f);
        BOOST_CHECK_CLOSE(endPoints.Krgr,  0.775f, 1.0e-6f);
        BOOST_CHECK_CLOSE(endPoints.Krorw, 0.7f,   1.0e-6f);
        BOOST_CHECK_CLOSE(endPoints.Krorg, 0.6f,   1.0e-6f);

        BOOST_CHECK_CLOSE(endPoints.maxKrw,  0.87f, 1.0e-6f);
        BOOST_CHECK_CLOSE(endPoints.maxKrow, 0.8f,  1.0e-6f);
        BOOST_CHECK_CLOSE(endPoints.maxKrog, 0.8f,  1.0e-6f);
        BOOST_CHECK_CLOSE(endPoints.maxKrg,  0.83f, 1.0e-6f);
    }

    // Scaled saturation function end-points in cell 4 (SATNUM = 1).
    {
        auto endPoints = Opm::EclEpsScalingPointsInfo<float>{};
        scaledCurvePt.populateCheckPoint(4, endPoints);

        BOOST_CHECK_CLOSE(endPoints.Swl, 0.27f, 1.0e-6f);
        BOOST_CHECK_CLOSE(endPoints.Sgl, 0.05f, 1.0e-6f);

        BOOST_CHECK_CLOSE(endPoints.Swcr,  0.3f,  1.0e-6f);
        BOOST_CHECK_CLOSE(endPoints.Sgcr,  0.14f, 1.0e-6f);
        BOOST_CHECK_CLOSE(endPoints.Sowcr, 0.29f, 1.0e-6f);
        BOOST_CHECK_CLOSE(endPoints.Sogcr, 0.19f, 1.0e-6f);

        BOOST_CHECK_CLOSE(endPoints.Swu, 0.94f, 1.0e-6f);
        BOOST_CHECK_CLOSE(endPoints.Sgu, 0.93f, 1.0e-6f);

        BOOST_CHECK_CLOSE(endPoints.maxPcow, 10.0f, 1.0e-6f);
        BOOST_CHECK_CLOSE(endPoints.maxPcgo,  0.0f, 1.0e-6f);

        BOOST_CHECK_CLOSE(endPoints.pcowLeverettFactor, 1.0f, 1.0e-6f);
        BOOST_CHECK_CLOSE(endPoints.pcgoLeverettFactor, 1.0f, 1.0e-6f);

        BOOST_CHECK_CLOSE(endPoints.Krwr,  0.65f,  1.0e-6f);
        BOOST_CHECK_CLOSE(endPoints.Krgr,  0.775f, 1.0e-6f);
        BOOST_CHECK_CLOSE(endPoints.Krorw, 0.6f,   1.0e-6f);
        BOOST_CHECK_CLOSE(endPoints.Krorg, 0.62f,  1.0e-6f);

        BOOST_CHECK_CLOSE(endPoints.maxKrw,  0.92f, 1.0e-6f);
        BOOST_CHECK_CLOSE(endPoints.maxKrow, 0.75f, 1.0e-6f);
        BOOST_CHECK_CLOSE(endPoints.maxKrog, 0.75f, 1.0e-6f);
        BOOST_CHECK_CLOSE(endPoints.maxKrg,  0.87f, 1.0e-6f);
    }
}

BOOST_AUTO_TEST_SUITE_END() // Scaled_Curve
