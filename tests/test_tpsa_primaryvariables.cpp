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

#define BOOST_TEST_MODULE TpsaPrimaryVariablesTests

#include <boost/test/unit_test.hpp>

#include <opm/material/materialstates/MaterialStateTPSA.hpp>

#include <opm/models/tpsa/elasticityprimaryvariables.hpp>
#include <opm/models/utils/propertysystem.hh>

#include <opm/simulators/flow/BlackoilModelProperties.hpp>
#include <opm/simulators/flow/TTagFlowProblemTPSA.hpp>

#include <tuple>


namespace Opm::Properties::TTag {
    struct TpsaTestTypeTag {
        using InheritsFrom = std::tuple<FlowProblem, FlowProblemTpsa>;
    };
}

BOOST_AUTO_TEST_CASE(ElasticityPrimVarTest) {
    using TypeTag = Opm::Properties::TTag::TpsaTestTypeTag;
    using Scalar = Opm::GetPropType<TypeTag, Opm::Properties::Scalar>;

    // Check straightforward assignment
    Opm::ElasticityPrimaryVariables<TypeTag> priVars;
    for (std::size_t i = 0; i < priVars.size(); ++i) {
        priVars[i] = static_cast<double>(i);
    }
    BOOST_CHECK_EQUAL(priVars.size(), 7U);
    for (std::size_t i = 0; i < priVars.size(); ++i) {
        const auto eval_var = priVars.makeEvaluation(i, 0);
        BOOST_CHECK_CLOSE(eval_var.value(), static_cast<double>(i), 1.0e-6);
        BOOST_CHECK_CLOSE(eval_var.derivative(i), 1.0, 1.0e-6);
    }

    // Assign from material state container
    Opm::MaterialStateTPSA<Scalar> ms;
    for (int i = 0; i < 3; ++i) {
        ms.setDisplacement(i, static_cast<double>(i) + 10.0);
        ms.setRotation(i, static_cast<double>(i) + 20.0);
    }
    ms.setSolidPressure(30.0);
    Opm::ElasticityPrimaryVariables<TypeTag> priVarsFromMs;
    priVarsFromMs.assignNaive(ms);
    for (std::size_t i = 0; i < 3; ++i) {
        const auto eval_disp = priVarsFromMs.makeEvaluation(i, 0);
        BOOST_CHECK_CLOSE(eval_disp.value(), ms.displacement(i), 1.0e-6);
        BOOST_CHECK_CLOSE(eval_disp.derivative(i), 1.0, 1.0e-6);

        const auto eval_rot = priVarsFromMs.makeEvaluation(i + 3, 0);
        BOOST_CHECK_CLOSE(eval_rot.value(), ms.rotation(i), 1.0e-6);
        BOOST_CHECK_CLOSE(eval_rot.derivative(i + 3), 1.0, 1.0e-6);
    }
    const auto eval_sp = priVarsFromMs.makeEvaluation(6, 0);
    BOOST_CHECK_CLOSE(eval_sp.value(), ms.solidPressure(), 1.0e-6);
    BOOST_CHECK_CLOSE(eval_sp.derivative(6), 1.0, 1.0e-6);
}
