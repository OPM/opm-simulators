/*
  Copyright 2023 Equinor.

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

#define BOOST_TEST_MODULE Running_Statistics

#include <boost/test/unit_test.hpp>

#include <opm/simulators/wells/RunningStatistics.hpp>

#include <boost/mpl/list.hpp>

#include <cstddef>
#include <initializer_list>
#include <limits>
#include <optional>

namespace {

#if FLOW_INSTANTIATE_FLOAT
    using Types = boost::mpl::list<float, double>;
#else
    using Types = boost::mpl::list<double>;
#endif  // FLOW_INSTANTIATE_FLOAT

} // Anonymous namespace

BOOST_AUTO_TEST_CASE_TEMPLATE(Empty, Scalar, Types)
{
    const auto rstat = Opm::RunningStatistics<Scalar>{};

    BOOST_CHECK_EQUAL(rstat.sampleSize(), std::size_t{0});
    BOOST_CHECK_MESSAGE(! (rstat.min() < std::numeric_limits<Scalar>::max()),
                        "Minimum value must not be less than Scalar::max() in empty sample");
    BOOST_CHECK_MESSAGE(! (rstat.max() > std::numeric_limits<Scalar>::min()),
                        "Maximum value must not be greater than Scalar::min() in empty sample");
    BOOST_CHECK_MESSAGE(! rstat.stdev().has_value(),
                        "Standard deviation must not exist in empty sample");
}

BOOST_AUTO_TEST_CASE_TEMPLATE(SinglePoint, Scalar, Types)
{
    auto rstat = Opm::RunningStatistics<Scalar>{};

    const auto samplePoint = static_cast<Scalar>(17.29);
    rstat.addSamplePoint(samplePoint);

    BOOST_CHECK_EQUAL(rstat.sampleSize(), std::size_t{1});
    BOOST_CHECK_CLOSE(rstat.min(), samplePoint, static_cast<Scalar>(1.0e-8));
    BOOST_CHECK_CLOSE(rstat.max(), samplePoint, static_cast<Scalar>(1.0e-8));
    BOOST_CHECK_CLOSE(rstat.mean(), samplePoint, static_cast<Scalar>(1.0e-8));
    BOOST_CHECK_MESSAGE(! rstat.stdev().has_value(),
                        "Standard deviation must not exist in single point sample");
}

BOOST_AUTO_TEST_CASE_TEMPLATE(Constant_Sampled_Twice, Scalar, Types)
{
    auto rstat = Opm::RunningStatistics<Scalar>{};

    const auto samplePoint = static_cast<Scalar>(17.29);
    rstat.addSamplePoint(samplePoint);
    rstat.addSamplePoint(samplePoint);

    BOOST_CHECK_EQUAL(rstat.sampleSize(), std::size_t{2});
    BOOST_CHECK_CLOSE(rstat.min(), samplePoint, static_cast<Scalar>(1.0e-8));
    BOOST_CHECK_CLOSE(rstat.max(), samplePoint, static_cast<Scalar>(1.0e-8));
    BOOST_CHECK_CLOSE(rstat.mean(), samplePoint, static_cast<Scalar>(1.0e-8));
    BOOST_CHECK_MESSAGE(rstat.stdev().has_value(),
                        "Standard deviation must exist in two-point sample");
    BOOST_CHECK_CLOSE(*rstat.stdev(), static_cast<Scalar>(0), static_cast<Scalar>(1.0e-8));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(Constant_Sampled_Twice_Reset, Scalar, Types)
{
    auto rstat = Opm::RunningStatistics<Scalar>{};

    const auto samplePoint = static_cast<Scalar>(17.29);
    rstat.addSamplePoint(samplePoint);
    rstat.addSamplePoint(samplePoint);

    BOOST_CHECK_EQUAL(rstat.sampleSize(), std::size_t{2});

    rstat.reset();

    BOOST_CHECK_EQUAL(rstat.sampleSize(), std::size_t{0});
    BOOST_CHECK_MESSAGE(! (rstat.min() < std::numeric_limits<Scalar>::max()),
                        "Minimum value must not be less than Scalar::max() in reset() sample");
    BOOST_CHECK_MESSAGE(! (rstat.max() > std::numeric_limits<Scalar>::min()),
                        "Maximum value must not be greater than Scalar::min() in reset() sample");
    BOOST_CHECK_MESSAGE(! rstat.stdev().has_value(),
                        "Standard deviation must not exist in reset() sample");
}

BOOST_AUTO_TEST_CASE_TEMPLATE(Two_Samples_Zero_Avg, Scalar, Types)
{
    auto rstat = Opm::RunningStatistics<Scalar>{};

    const auto samplePoint = static_cast<Scalar>(17.29);
    rstat.addSamplePoint(  samplePoint);
    rstat.addSamplePoint(- samplePoint);

    BOOST_CHECK_EQUAL(rstat.sampleSize(), std::size_t{2});
    BOOST_CHECK_CLOSE(rstat.min(), - samplePoint, static_cast<Scalar>(1.0e-8));
    BOOST_CHECK_CLOSE(rstat.max(),   samplePoint, static_cast<Scalar>(1.0e-8));
    BOOST_CHECK_CLOSE(rstat.mean(), static_cast<Scalar>(0), static_cast<Scalar>(1.0e-8));
    BOOST_CHECK_CLOSE(*rstat.stdev(), static_cast<Scalar>(24.451752493), static_cast<Scalar>(1.0e-5));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(One_Two_Five_Repetitions, Scalar, Types)
{
    auto rstat = Opm::RunningStatistics<Scalar>{};

    const auto rpt = 5;
    const auto one = static_cast<Scalar>(1);
    const auto two = static_cast<Scalar>(2);

    for (auto i = 0*rpt; i < rpt; ++i) { rstat.addSamplePoint(one); }
    for (auto i = 0*rpt; i < rpt; ++i) { rstat.addSamplePoint(two); }

    BOOST_CHECK_EQUAL(rstat.sampleSize(), 2 * static_cast<std::size_t>(rpt));
    BOOST_CHECK_CLOSE(rstat.min(), one, static_cast<Scalar>(1.0e-8));
    BOOST_CHECK_CLOSE(rstat.max(), two, static_cast<Scalar>(1.0e-8));
    BOOST_CHECK_CLOSE(rstat.mean(), static_cast<Scalar>(1.5), static_cast<Scalar>(1.0e-8));
    BOOST_CHECK_CLOSE(*rstat.stdev(), static_cast<Scalar>(5.27046276695e-01), static_cast<Scalar>(1.0e-8));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(One_Two_Five_Repetitions_Interleaved, Scalar, Types)
{
    auto rstat = Opm::RunningStatistics<Scalar>{};

    const auto rpt = 5;
    const auto one = static_cast<Scalar>(1);
    const auto two = static_cast<Scalar>(2);

    for (auto i = 0*rpt; i < rpt; ++i) {
        rstat.addSamplePoint(one);
        rstat.addSamplePoint(two);
    }

    BOOST_CHECK_EQUAL(rstat.sampleSize(), 2 * static_cast<std::size_t>(rpt));
    BOOST_CHECK_CLOSE(rstat.min(), one, static_cast<Scalar>(1.0e-8));
    BOOST_CHECK_CLOSE(rstat.max(), two, static_cast<Scalar>(1.0e-8));
    BOOST_CHECK_CLOSE(rstat.mean(), static_cast<Scalar>(1.5), static_cast<Scalar>(1.0e-8));
    BOOST_CHECK_CLOSE(*rstat.stdev(), static_cast<Scalar>(5.27046276695e-01), static_cast<Scalar>(1.0e-8));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(One_Two_Twenty_Repetitions, Scalar, Types)
{
    auto rstat = Opm::RunningStatistics<Scalar>{};

    const auto rpt = 20;
    const auto one = static_cast<Scalar>(1);
    const auto two = static_cast<Scalar>(2);

    for (auto i = 0*rpt; i < rpt; ++i) { rstat.addSamplePoint(one); }
    for (auto i = 0*rpt; i < rpt; ++i) { rstat.addSamplePoint(two); }

    BOOST_CHECK_EQUAL(rstat.sampleSize(), 2 * static_cast<std::size_t>(rpt));
    BOOST_CHECK_CLOSE(rstat.min(), one, static_cast<Scalar>(1.0e-8));
    BOOST_CHECK_CLOSE(rstat.max(), two, static_cast<Scalar>(1.0e-8));
    BOOST_CHECK_CLOSE(rstat.mean(), static_cast<Scalar>(1.5), static_cast<Scalar>(1.0e-5));
    BOOST_CHECK_CLOSE(*rstat.stdev(), static_cast<Scalar>(5.06369683542e-01), static_cast<Scalar>(5.0e-5));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(RandN_50, Scalar, Types)
{
    auto rstat = Opm::RunningStatistics<Scalar>{};

    for (const auto& samplePoint : {
            // Values from Octave (RANDN([50, 1]))

            -5.983051973110e-01, -2.129185497807e-01, -5.160782188139e-01,
            -5.405962754204e-02, 1.081538561527e+00, -5.444361014868e-01,
            6.802216201756e-01, -1.732112159825e+00, -1.681240838370e+00,
            9.026897835500e-01, -2.510481192847e-01, 2.503216953163e+00,
            2.096538154630e+00, -1.157169480540e+00, 2.777014458586e-01,
            3.206654292492e-01, 2.720894169117e+00, 3.394916550830e-01,
            1.603500564897e+00, 1.305786976012e+00, 3.358587787445e-01,
            1.688835457016e+00, 6.464543963139e-01, -8.880888352071e-01,
            1.785948404587e+00, 7.344602418137e-01, 1.272049108856e+00,
            3.618201220834e-02, -1.254183439007e+00, 8.551411128509e-01,
            2.002540536438e+00, -8.442308733039e-01, -5.880385774749e-01,
            5.134590162252e-01, 2.242601346140e-01, -1.632624091153e+00,
            -2.041052498197e-01, 8.535062014928e-01, 1.218406596883e-01,
            4.866545018493e-01, -1.249277350665e+00, -5.014606488004e-01,
            2.795291286626e-01, -6.459643961814e-01, -1.061751877930e+00,
            -1.422156588627e+00, 1.026939772058e-01, -3.551330895391e-01,
            8.945851907053e-01, -2.250951241491e-01,
        })
    {
        rstat.addSamplePoint(static_cast<Scalar>(samplePoint));
    }

    BOOST_CHECK_EQUAL(rstat.sampleSize(), std::size_t{50});
    BOOST_CHECK_CLOSE(rstat.min(), static_cast<Scalar>(-1.732112159824814), static_cast<Scalar>(1.0e-6));
    BOOST_CHECK_CLOSE(rstat.max(), static_cast<Scalar>(2.720894169117450), static_cast<Scalar>(1.0e-6));
    BOOST_CHECK_CLOSE(rstat.mean(), static_cast<Scalar>(0.1809353147544384), static_cast<Scalar>(5.0e-5));
    BOOST_CHECK_CLOSE(*rstat.stdev(), static_cast<Scalar>(1.097155256132325), static_cast<Scalar>(1.0e-6));
}
