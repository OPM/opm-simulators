/*
  Copyright 2014 SINTEF ICT, Applied Mathematics.
  Copyright 2014 Statoil ASA.

  This file is part of the Open Porous Media Project (OPM).

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

#if HAVE_DYNAMIC_BOOST_TEST
#define BOOST_TEST_DYN_LINK
#endif

#define BOOST_TEST_MODULE VoidageRateConversionTest

#include <opm/autodiff/RateConverter.hpp>

#include <opm/autodiff/BlackoilPropsAdFromDeck.hpp>

#include <boost/test/unit_test.hpp>

#include <opm/core/grid/GridHelpers.hpp>
#include <opm/core/grid/GridManager.hpp>
#include <opm/core/props/BlackoilPropertiesFromDeck.hpp>
#include <opm/parser/eclipse/Units/Units.hpp>
#include <opm/common/utility/parameters/ParameterGroup.hpp>
#include <opm/core/simulator/BlackoilState.hpp>

#include <opm/parser/eclipse/Parser/Parser.hpp>
#include <opm/parser/eclipse/Parser/ParseContext.hpp>
#include <opm/parser/eclipse/EclipseState/EclipseState.hpp>
#include <opm/parser/eclipse/Deck/Deck.hpp>


struct SetupSimple {
    SetupSimple() :
        deck( Opm::Parser{}.parseFile( "fluid.data" ) ),
        eclState( deck, Opm::ParseContext() )
    {
        param.disableOutput();
        param.insertParameter("init_rock"       , "false" );
        param.insertParameter("threephase_model", "simple");
        param.insertParameter("pvt_tab_size"    , "0"     );
        param.insertParameter("sat_tab_size"    , "0"     );
    }

    Opm::ParameterGroup  param;
    Opm::Deck                       deck;
    Opm::EclipseState               eclState;
};


template <class Setup>
struct TestFixture : public Setup
{
    TestFixture()
        : Setup()
        , grid (eclState.getInputGrid())
        , ad_props(deck, eclState, *grid.c_grid(), param.getDefault("init_rock", false))
    {
    }

    using Setup::param;
    using Setup::deck;
    using Setup::eclState;

    Opm::GridManager             grid;
    Opm::BlackoilPropsAdFromDeck ad_props;
};


BOOST_FIXTURE_TEST_CASE(Construction, TestFixture<SetupSimple>)
{
    typedef std::vector<int>                     Region;
    typedef Opm::BlackoilPropsAdFromDeck         Props;
    typedef Opm::RateConverter::
        SurfaceToReservoirVoidage<Props::FluidSystem, Region> RCvrt;

    Region reg{ 0 };
    RCvrt  cvrt(ad_props.phaseUsage(), reg);
}


BOOST_FIXTURE_TEST_CASE(ThreePhase, TestFixture<SetupSimple>)
{
    // Immiscible and incompressible two-phase fluid
    typedef std::vector<int>                     Region;
    typedef Opm::BlackoilPropsAdFromDeck         Props;
    typedef Opm::RateConverter::
        SurfaceToReservoirVoidage<Props::FluidSystem, Region> RCvrt;

    Region reg{ 0 };
    int numCells = Opm::UgGridHelpers::numCells(*grid.c_grid());
    RCvrt  cvrt(ad_props.phaseUsage(), reg);

    Opm::BlackoilState x(numCells, Opm::UgGridHelpers::numFaces( *grid.c_grid()) , 3);

    cvrt.defineState(x);

    std::vector<double> coeff(3, 0.0);

    // Immiscible and incompressible: All coefficients are one (1),
    // irrespective of actual surface rates.
    cvrt.calcCoeff(0, 0, coeff);
    BOOST_CHECK_CLOSE(coeff[0], 1.0, 1.0e-6);
    BOOST_CHECK_CLOSE(coeff[1], 1.0, 1.0e-6);
    BOOST_CHECK_CLOSE(coeff[2], 1.0, 1.0e-6);
}
