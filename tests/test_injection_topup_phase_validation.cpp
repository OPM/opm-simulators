// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
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

  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the list of
  copyright holders.
*/
/*!
 * \file
 * \brief Test GCONINJE top-up phase validation across the group hierarchy
 *
 * Validates that conflicting RESV/VREP injection top-up phases are rejected
 * when they appear on the same group or on the same ancestor-descendant path,
 * while disconnected branches are allowed to use different top-up phases.
 *
 * Group hierarchy (from INJECTION_TOPUP_PHASE_VALIDATION.DATA):
 * \code
 *   FIELD
 *   ├── SAME (INJ2)
 *   ├── PARENT
 *   │   └── CHILD (INJ1)
 *   ├── DISCON-G (INJ3)
 *   └── DISCON-W (PROD1)
 * \endcode
 *
 * Test cases:
 *   - SameGroupConflictThrows: report step 1 — SAME has both WATER RESV
 *     and GAS RESV, expect rejection.
 *   - AncestorConflictThrows: report step 2 — PARENT has WATER VREP while
 *     CHILD has GAS RESV on the same path, expect rejection.
 *   - DisconnectedBranchesAreAllowed: report step 3 — DISCON-G has GAS RESV
 *     and DISCON-W has WATER VREP on separate branches, expect no error.
 */

#define BOOST_TEST_MODULE InjectionTopupPhaseValidation

#include "SimulatorFixture.hpp"

#include <opm/simulators/flow/FlowProblemBlackoil.hpp>
#include <opm/simulators/wells/BlackoilWellModel.hpp>

#include <boost/test/unit_test.hpp>

#include <memory>
#include <string_view>

namespace Opm::Properties::TTag {
    struct TestInjectionTopupPhaseTypeTag {
        using InheritsFrom = std::tuple<TestTypeTag>;
    };
}

using SimulatorFixture = Opm::SimulatorFixture;
BOOST_GLOBAL_FIXTURE(SimulatorFixture);

namespace {

template <class TypeTag>
std::unique_ptr<Opm::Simulator<TypeTag>>
initTopupPhaseSimulator(const char* filename)
{
    return Opm::initSimulator<TypeTag>(filename, "test_injection_topup_phase_validation");
}

bool hasTopupPhaseError(const std::runtime_error& error)
{
    return std::string_view{error.what()}.find("Invalid GCONINJE top-up phase configuration")
        != std::string_view::npos;
}

template <class TypeTag>
void validateReportStepTopupConfiguration(Opm::Simulator<TypeTag>& simulator, const int report_step)
{
    auto& well_model = simulator.problem().wellModel();
    auto& helper = well_model.groupStateHelper();
    helper.setReportStep(report_step);
    const auto& field_group = simulator.vanguard().schedule().getGroup("FIELD", report_step);
    helper.setCmodeGroup(field_group);
}

} // namespace

BOOST_AUTO_TEST_CASE(SameGroupConflictThrows)
{
    using TypeTag = Opm::Properties::TTag::TestInjectionTopupPhaseTypeTag;

    auto simulator = initTopupPhaseSimulator<TypeTag>("INJECTION_TOPUP_PHASE_VALIDATION.DATA");

    BOOST_CHECK_EXCEPTION(
        validateReportStepTopupConfiguration(*simulator, 1),
        std::runtime_error,
        hasTopupPhaseError
    );
}

BOOST_AUTO_TEST_CASE(AncestorConflictThrows)
{
    using TypeTag = Opm::Properties::TTag::TestInjectionTopupPhaseTypeTag;

    auto simulator = initTopupPhaseSimulator<TypeTag>("INJECTION_TOPUP_PHASE_VALIDATION.DATA");

    BOOST_CHECK_EXCEPTION(
        validateReportStepTopupConfiguration(*simulator, 2),
        std::runtime_error,
        hasTopupPhaseError
    );
}

BOOST_AUTO_TEST_CASE(DisconnectedBranchesAreAllowed)
{
    using TypeTag = Opm::Properties::TTag::TestInjectionTopupPhaseTypeTag;

    auto simulator = initTopupPhaseSimulator<TypeTag>("INJECTION_TOPUP_PHASE_VALIDATION.DATA");

    BOOST_CHECK_NO_THROW(validateReportStepTopupConfiguration(*simulator, 3));
}
