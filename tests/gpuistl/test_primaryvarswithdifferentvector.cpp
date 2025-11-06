/*
  Copyright 2025 Equinor ASA
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
#define HAVE_ECL_INPUT 1
#include <config.h>


#include <stdexcept>

#define BOOST_TEST_MODULE TestPrimaryVariablesWithDifferentVector

#include <boost/test/unit_test.hpp>
#include <opm/common/ErrorMacros.hpp>
#include <opm/models/blackoil/blackoilmodel.hh>
#include <opm/models/blackoil/blackoilprimaryvariables.hh>
#include <opm/models/discretization/common/fvbaseprimaryvariables.hh>
#include <opm/models/discretization/common/tpfalinearizer.hh>
#include <opm/models/io/dgfvanguard.hh>
#include <opm/models/utils/start.hh>
#include <opm/simulators/flow/Main.hpp>
#include <opm/simulators/linalg/gpuistl/MiniVector.hpp>
#include <opm/simulators/linalg/gpuistl/detail/gpu_safe_call.hpp>
#include <opm/simulators/linalg/gpuistl/gpu_smart_pointer.hpp>

BOOST_AUTO_TEST_CASE(TestPrimaryVariablesCreationWithFieldVector)
{
    // Just make sure we can instantiate
    using TypeTag = Opm::Properties::TTag::FlowProblem;
    Opm::BlackOilPrimaryVariables<TypeTag> primaryVariables;
    for (std::size_t i = 0; i < primaryVariables.size(); ++i) {
        primaryVariables[i] = static_cast<double>(i);
    }

    primaryVariables.setPrimaryVarsMeaningBrine(Opm::BlackOil::BrineMeaning::Sp);
    primaryVariables.setPrimaryVarsMeaningGas(Opm::BlackOil::GasMeaning::Rv);
    primaryVariables.setPrimaryVarsMeaningWater(Opm::BlackOil::WaterMeaning::Rvw);
    primaryVariables.setPrimaryVarsMeaningPressure(Opm::BlackOil::PressureMeaning::Pg);
    primaryVariables.setPrimaryVarsMeaningSolvent(Opm::BlackOil::SolventMeaning::Rsolw);

    // Now we instantiate with MiniVector
    Opm::BlackOilPrimaryVariables<TypeTag, Opm::gpuistl::MiniVector> primaryVariablesMiniVector(primaryVariables);

    // Check results
    for (std::size_t i = 0; i < primaryVariablesMiniVector.size(); ++i) {
        BOOST_CHECK_EQUAL(primaryVariablesMiniVector[i], static_cast<double>(i));
    }
    BOOST_CHECK(primaryVariablesMiniVector.primaryVarsMeaningBrine() == Opm::BlackOil::BrineMeaning::Sp);
    BOOST_CHECK(primaryVariablesMiniVector.primaryVarsMeaningGas() == Opm::BlackOil::GasMeaning::Rv);
    BOOST_CHECK(primaryVariablesMiniVector.primaryVarsMeaningWater() == Opm::BlackOil::WaterMeaning::Rvw);
    BOOST_CHECK(primaryVariablesMiniVector.primaryVarsMeaningPressure() == Opm::BlackOil::PressureMeaning::Pg);
    BOOST_CHECK(primaryVariablesMiniVector.primaryVarsMeaningSolvent() == Opm::BlackOil::SolventMeaning::Rsolw);
}
