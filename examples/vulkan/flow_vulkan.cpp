// SPDX-License-Identifier: GPL-3.0-or-later
#include <opm/grid/CpGrid.hpp>
#include <opm/models/blackoil/blackoillocalresidualtpfa.hh>
#include <opm/models/discretization/common/tpfalinearizer.hh>
#include <opm/simulators/flow/Main.hpp>
#include <opm/simulators/flow/SimulatorFullyImplicitBlackoil.hpp>
#include <opm/simulators/flow/TTagFlowProblemTPFA.hpp>

// The OPM solver interface expects the Flow property and Dune declarations above.
#include "VulkanISTLSolver.hpp"

namespace Opm::Properties
{
namespace TTag
{
    struct ReslabVulkanFlow {
        using InheritsFrom = std::tuple<FlowProblemTPFA>;
    };
} // namespace TTag
template <class TypeTag>
struct Linearizer<TypeTag, TTag::ReslabVulkanFlow> {
    using type = TpfaLinearizer<TypeTag>;
};
template <class TypeTag>
struct LocalResidual<TypeTag, TTag::ReslabVulkanFlow> {
    using type = BlackOilLocalResidualTPFA<TypeTag>;
};
template <class TypeTag>
struct EnableDiffusion<TypeTag, TTag::ReslabVulkanFlow> {
    static constexpr bool value = false;
};
template <class TypeTag>
struct AvoidElementContext<TypeTag, TTag::ReslabVulkanFlow> {
    static constexpr bool value = true;
};
#ifndef RESLAB_CPU_REFERENCE
template <class TypeTag>
struct LinearSolverBackend<TypeTag, TTag::ReslabVulkanFlow> {
    using type = Reslab::VulkanISTLSolver<TypeTag>;
};
#endif
} // namespace Opm::Properties
int
main(int argc, char** argv)
{
    auto simulator = std::make_unique<Opm::Main>(argc, argv);
    const int result = simulator->runStatic<Opm::Properties::TTag::ReslabVulkanFlow>();
    simulator.reset();
    return result;
}
