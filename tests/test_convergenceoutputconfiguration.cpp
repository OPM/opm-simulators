/*
  Copyright 2022 Equinor.

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

#define BOOST_TEST_MODULE TestConvergenceOutputConfiguration

#include <boost/test/unit_test.hpp>

#include <opm/simulators/flow/ConvergenceOutputConfiguration.hpp>

#include <cstddef>
#include <stdexcept>
#include <string_view>
#include <string>

BOOST_AUTO_TEST_SUITE(Common_Operations)

BOOST_AUTO_TEST_CASE(None)
{
    const auto config = Opm::ConvergenceOutputConfiguration{"none"};
    BOOST_CHECK_MESSAGE(! config.any(),
                        "Configuration object with option value "
                        "\"none\" must NOT activate output");
}

BOOST_AUTO_TEST_CASE(Steps)
{
    const auto config = Opm::ConvergenceOutputConfiguration{"steps"};
    BOOST_CHECK_MESSAGE(config.any(),
                        "Configuration object with supported "
                        "option value must activate option");
    BOOST_CHECK_MESSAGE(config.want(Opm::ConvergenceOutputConfiguration::Option::Steps),
                        "Configuration object with \"steps\" "
                        "option value must activate Steps option");
}

BOOST_AUTO_TEST_CASE(Steps_Alias)
{
    const auto config = Opm::ConvergenceOutputConfiguration{"step"};
    BOOST_CHECK_MESSAGE(config.any(),
                        "Configuration object with supported "
                        "option value must activate option");
    BOOST_CHECK_MESSAGE(config.want(Opm::ConvergenceOutputConfiguration::Option::Steps),
                        "Configuration object with \"step\" "
                        "option value must activate Steps option");
}

BOOST_AUTO_TEST_CASE(Iterations)
{
    const auto config = Opm::ConvergenceOutputConfiguration{"iterations"};
    BOOST_CHECK_MESSAGE(config.any(),
                        "Configuration object with supported "
                        "option value must activate option");
    BOOST_CHECK_MESSAGE(config.want(Opm::ConvergenceOutputConfiguration::Option::Iterations),
                        "Configuration object with \"iterations\" "
                        "option value must activate Steps option");
}

BOOST_AUTO_TEST_CASE(Iterations_Alias)
{
    const auto config = Opm::ConvergenceOutputConfiguration{"iteration"};
    BOOST_CHECK_MESSAGE(config.any(),
                        "Configuration object with supported "
                        "option value must activate option");
    BOOST_CHECK_MESSAGE(config.want(Opm::ConvergenceOutputConfiguration::Option::Iterations),
                        "Configuration object with \"iterations\" "
                        "option value must activate Steps option");
}

BOOST_AUTO_TEST_CASE(Combinations)
{
    const auto steps_iter = Opm::ConvergenceOutputConfiguration{"steps,iterations"};
    BOOST_CHECK_MESSAGE(steps_iter.any(),
                        "Configuration object with supported "
                        "option value must activate option");
    BOOST_CHECK_MESSAGE(steps_iter.want(Opm::ConvergenceOutputConfiguration::Option::Steps),
                        "Configuration object with \"steps\" "
                        "option value must activate Steps option");
    BOOST_CHECK_MESSAGE(steps_iter.want(Opm::ConvergenceOutputConfiguration::Option::Iterations),
                        "Configuration object with \"iterations\" "
                        "option value must activate Steps option");

    const auto iter_steps = Opm::ConvergenceOutputConfiguration{"iterations,steps"};
    BOOST_CHECK_MESSAGE(iter_steps.any(),
                        "Configuration object with supported "
                        "option value must activate option");
    BOOST_CHECK_MESSAGE(iter_steps.want(Opm::ConvergenceOutputConfiguration::Option::Steps),
                        "Configuration object with \"steps\" "
                        "option value must activate Steps option");
    BOOST_CHECK_MESSAGE(iter_steps.want(Opm::ConvergenceOutputConfiguration::Option::Iterations),
                        "Configuration object with \"iterations\" "
                        "option value must activate Steps option");

    const auto none_iter_steps = Opm::ConvergenceOutputConfiguration{"none,iterations,steps"};
    BOOST_CHECK_MESSAGE(! none_iter_steps.any(),
                        "Configuration object with any option "
                        "value \"none\" must NOT activate output");

    const auto iter_none_steps = Opm::ConvergenceOutputConfiguration{"iterations,none,steps"};
    BOOST_CHECK_MESSAGE(! iter_none_steps.any(),
                        "Configuration object with any option "
                        "value \"none\" must NOT activate output");

    const auto steps_iter_none = Opm::ConvergenceOutputConfiguration{"steps,iterations,   none"};
    BOOST_CHECK_MESSAGE(! steps_iter_none.any(),
                        "Configuration object with any option "
                        "value \"none\" must NOT activate output");
}

BOOST_AUTO_TEST_SUITE_END() // Common_Operations

// ---------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE(Failed_Construction)

BOOST_AUTO_TEST_CASE(Misprint)
{
    BOOST_CHECK_THROW(Opm::ConvergenceOutputConfiguration{"nonce"},
                      std::invalid_argument);

    BOOST_CHECK_THROW(Opm::ConvergenceOutputConfiguration("nonce", "X"),
                      std::invalid_argument);

    BOOST_CHECK_THROW(Opm::ConvergenceOutputConfiguration{"stepS"},
                      std::invalid_argument);

    BOOST_CHECK_THROW(Opm::ConvergenceOutputConfiguration{"steps, iter"},
                      std::invalid_argument);

    BOOST_CHECK_THROW(Opm::ConvergenceOutputConfiguration{"steps, iterations, non"},
                      std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(Unknown)
{
    BOOST_CHECK_THROW(Opm::ConvergenceOutputConfiguration{"Hello"},
                      std::invalid_argument);
    BOOST_CHECK_THROW(Opm::ConvergenceOutputConfiguration{"meow"},
                      std::invalid_argument);
    BOOST_CHECK_THROW(Opm::ConvergenceOutputConfiguration{""},
                      std::invalid_argument);
    BOOST_CHECK_THROW(Opm::ConvergenceOutputConfiguration{"xyz,zy;;;"},
                      std::invalid_argument);
}

BOOST_AUTO_TEST_SUITE_END()
