// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 2 of the License, or
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
/*!
 * \file
 *
 * \brief This test makes sure that mandated API is adhered to by all component classes
 */
#include "config.h"

#include <opm/material/densead/Evaluation.hpp>
#include <opm/material/densead/Math.hpp>

#include "checkComponent.hpp"

// include all components shipped with opm-material
#include <opm/material/components/Unit.hpp>
#include <opm/material/components/NullComponent.hpp>
#include <opm/material/components/Component.hpp>
#include <opm/material/components/Dnapl.hpp>
#include <opm/material/components/SimpleH2O.hpp>
#include <opm/material/components/Lnapl.hpp>
#include <opm/material/components/iapws/Region2.hpp>
#include <opm/material/components/iapws/Region1.hpp>
#include <opm/material/components/iapws/Common.hpp>
#include <opm/material/components/iapws/Region4.hpp>
#include <opm/material/components/H2O.hpp>
#include <opm/material/components/SimpleHuDuanH2O.hpp>
#include <opm/material/components/CO2.hpp>
#include <opm/material/components/Mesitylene.hpp>
#include <opm/material/components/TabulatedComponent.hpp>
#include <opm/material/components/Brine.hpp>
#include <opm/material/components/N2.hpp>
#include <opm/material/components/Xylene.hpp>
#include <opm/material/components/Air.hpp>
#include <opm/material/components/SimpleCO2.hpp>

#include <opm/material/common/UniformTabulated2DFunction.hpp>

namespace Opm {
namespace ComponentsTest {
#include <opm/material/components/co2tables.inc>
}}

#include <dune/common/parallel/mpihelper.hh>

template <class Scalar, class Evaluation>
void testSimpleH2O()
{
    typedef Opm::H2O<Scalar> H2O;
    typedef Opm::SimpleHuDuanH2O<Scalar> SimpleHuDuanH2O;
    typedef Opm::MathToolbox<Evaluation> EvalToolbox;

    int numT = 67;
    int numP = 45;
    Evaluation T = 280;
    Evaluation p = 1e6;

    for (int iT = 0; iT < numT; ++iT) {
        p = 1e6;
        T += 5;
        for (int iP = 0; iP < numP; ++iP) {
            p *= 1.1;
            if (!EvalToolbox::isSame(H2O::liquidDensity(T,p), SimpleHuDuanH2O::liquidDensity(T,p), /*tolerance=*/1e-3*H2O::liquidDensity(T,p).value()))
                throw std::logic_error("oops: the water density based on Hu-Duan has more then 1e-3 deviation from IAPWS'97");

            if (T >= 570) // for temperature larger then 570 the viscosity based on HuDuan is too far from IAPWS.
                continue;

            if (!EvalToolbox::isSame(H2O::liquidViscosity(T,p), SimpleHuDuanH2O::liquidViscosity(T,p), /*tolerance=*/5.e-2*H2O::liquidViscosity(T,p).value())){
                throw std::logic_error("oops: the water viscosity based on Hu-Duan has more then 5e-2 deviation from IAPWS'97");
            }
        }
    }
}

template <class Scalar, class Evaluation>
void testAllComponents()
{
    typedef Opm::H2O<Scalar> H2O;

    checkComponent<Opm::Air<Scalar>, Evaluation>();
    checkComponent<Opm::Brine<Scalar, H2O>, Evaluation>();
    checkComponent<Opm::CO2<Scalar, Opm::ComponentsTest::CO2Tables>, Evaluation>();
    checkComponent<Opm::DNAPL<Scalar>, Evaluation>();
    checkComponent<Opm::H2O<Scalar>, Evaluation>();
    checkComponent<Opm::LNAPL<Scalar>, Evaluation>();
    checkComponent<Opm::Mesitylene<Scalar>, Evaluation>();
    checkComponent<Opm::N2<Scalar>, Evaluation>();
    checkComponent<Opm::NullComponent<Scalar>, Evaluation>();
    checkComponent<Opm::SimpleCO2<Scalar>, Evaluation>();
    checkComponent<Opm::SimpleH2O<Scalar>, Evaluation>();
    checkComponent<Opm::TabulatedComponent<Scalar, H2O>, Evaluation>();
    checkComponent<Opm::Unit<Scalar>, Evaluation>();
    checkComponent<Opm::Xylene<Scalar>, Evaluation>();
}

template <class Scalar>
inline void testAll()
{
    typedef Opm::DenseAd::Evaluation<Scalar, 3> Evaluation;

    // ensure that all components are API-compliant
    testAllComponents<Scalar, Scalar>();
    testAllComponents<Scalar, Evaluation>();
    testSimpleH2O<Scalar, Evaluation>();

}


int main(int argc, char **argv)
{
    Dune::MPIHelper::instance(argc, argv);

    testAll<double>();
    testAll<float>();

    return 0;
}
