/* Copyright (c) 2013 Uni Research AS.
   This file is licensed under the GNU General Public License v3.0 or later. */
#include <config.h>

#if HAVE_DYNAMIC_BOOST_TEST
#define BOOST_TEST_DYN_LINK
#endif
#define NVERBOSE // to suppress our messages when throwing

#define BOOST_TEST_MODULE ShadowTest
#include <boost/test/unit_test.hpp>

#include <opm/core/utility/parameters/ParameterGroup.hpp>
#include <opm/core/props/IncompPropertiesBasic.hpp>
#include <opm/core/props/IncompPropertiesShadow.hpp>

using namespace Opm;

BOOST_AUTO_TEST_CASE(shadowPorosity)
{
    const double defaultPorosity = 1.0;
    const double newPorosity = 0.5;

    parameter::ParameterGroup param;
    IncompPropertiesBasic basic (param, 2, 1);
    IncompPropertiesShadow shadow (basic);
    BOOST_CHECK_CLOSE (*(shadow.porosity()), defaultPorosity, 0.001);
    shadow.usePorosity (&newPorosity);
    BOOST_CHECK_CLOSE (*(shadow.porosity()), newPorosity, 0.001);
    shadow.usePorosity (basic);
    BOOST_CHECK_CLOSE (*(shadow.porosity()), defaultPorosity, 0.001);
}
