// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  Copyright 2026 Equinor ASA

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
#ifndef OPM_GPU_BLACKOIL_NOHYST_INTENSIVE_QUANTITIES_HPP
#define OPM_GPU_BLACKOIL_NOHYST_INTENSIVE_QUANTITIES_HPP

#if HAVE_CUDA

#include <opm/models/utils/propertysystem.hh>

#include <cstddef>
#include <memory>
#include <span>

namespace Opm::Properties::TTag {
    struct FlowProblemNohyst;
}

namespace Opm {
template <class TypeTag> class TpfaLinearizerGpuParams;
}

namespace Opm::gpuistl {

namespace dispatcher_detail {
template <class CpuTypeTag, class PvtApproaches>
struct ImplFor;
}

/// Compile-time predicate: is the GPU intensive-quantities dispatcher
/// explicitly instantiated for this CPU \c TypeTag? Each supported TypeTag
/// has a list of configurations in \c GpuBlackoilNohystIntensiveQuantities.cpp.
/// Whether a given deck is supported is checked at run time, when the
/// dispatcher is initialized.
template <class CpuTypeTag>
struct GpuBlackoilNohystIntensiveQuantitiesSupport {
    static constexpr bool value = false;
};

/// Enable the dispatcher for the TypeTag that \c flow uses for black-oil decks
/// without hysteresis. Three-phase decks with live oil (PVTO), dry gas
/// (PVDG) or wet gas (PVTG with VAPOIL), and PVTW water are accepted at run time,
/// with or without endpoint scaling and SWATINIT.
template <>
struct GpuBlackoilNohystIntensiveQuantitiesSupport<
    Opm::Properties::TTag::FlowProblemNohyst> {
    static constexpr bool value = true;
};

/// Runs the supported BlackOil intensive-quantities update on the GPU for all
/// grid degrees of freedom. Each dispatcher instance owns its GPU-side state
/// and lazily constructs the \c GpuFlowProblem from the supplied CPU problem
/// on the first call. On that call it also chooses the configuration, from the
/// PVT approaches of the CPU fluid system, and validates the deck against it.
///
/// On every call, primary variables for the requested DoFs are uploaded to
/// the device, the per-cell update kernel is launched (one thread per DoF),
/// and the resulting intensive quantities are read back to host memory.
/// The supported configurations compute the complete intensive-quantity
/// state needed by this dispatcher, including mobility.
///
/// The class is a template on the CPU \c TypeTag and a facade: the GPU state
/// and the kernel launch of a configuration live behind the abstract \c Impl.
/// A configuration is a (CPU \c TypeTag, PVT approaches) pair. Each has its own
/// \c .cu translation unit with its \c GpuDispatcherTraits specialization and
/// the explicit instantiation of \c makeImpl, which builds its \c Impl (see
/// \c GpuBlackoilNohystIntensiveQuantities_impl.hpp). The facade itself is
/// host code, in \c GpuBlackoilNohystIntensiveQuantities.cpp.
template <class CpuTypeTag>
class GpuBlackoilNohystIntensiveQuantities
{
public:
    using Scalar = Opm::GetPropType<CpuTypeTag, Opm::Properties::Scalar>;
    using Problem = Opm::GetPropType<CpuTypeTag, Opm::Properties::Problem>;
    using SolutionVector = Opm::GetPropType<CpuTypeTag, Opm::Properties::SolutionVector>;
    using IntensiveQuantities = Opm::GetPropType<CpuTypeTag, Opm::Properties::IntensiveQuantities>;

    GpuBlackoilNohystIntensiveQuantities();
    ~GpuBlackoilNohystIntensiveQuantities();

    GpuBlackoilNohystIntensiveQuantities(const GpuBlackoilNohystIntensiveQuantities&) = delete;
    GpuBlackoilNohystIntensiveQuantities&
    operator=(const GpuBlackoilNohystIntensiveQuantities&) = delete;

    /// Update the complete CPU cache. The current solution is uploaded once
    /// and the evaluated IQ state is downloaded once for CPU assembly.
    void update(const Problem& cpuProblem,
                const SolutionVector& solution,
                std::span<IntensiveQuantities> outIQ,
                unsigned timeIdx);

    /// Assemble the reservoir from the dispatcher-owned current and previous IQs.
    void assemble(TpfaLinearizerGpuParams<CpuTypeTag>& buffers,
                  const Problem& cpuProblem, Scalar dt) const;

    /// Refresh the previous slot after the host recomputes its history, for example on restart.
    void invalidatePrevious();

    /// Publish the accepted current IQs as the previous timestep's state.
    void advanceTimeLevel();

    /// The name of the chosen configuration, or \c nullptr before the first
    /// call to \c update().
    const char* configurationName() const;

private:
    /// The GPU state and kernel launch of one configuration.
    struct Impl {
        virtual ~Impl() = default;
        virtual void update(const Problem& cpuProblem,
                            const SolutionVector& solution,
                            std::span<IntensiveQuantities> outIQ,
                            unsigned timeIdx) = 0;
        virtual void assemble(TpfaLinearizerGpuParams<CpuTypeTag>& buffers,
                              const Problem& cpuProblem, Scalar dt) = 0;
        virtual void invalidatePrevious() = 0;
        virtual void advanceTimeLevel() = 0;
        virtual const char* name() const = 0;
    };

    template <class, class>
    friend struct dispatcher_detail::ImplFor;

    /// Validate the deck against the configuration of \p PvtApproaches and
    /// build its \c Impl. Explicitly instantiated in that configuration's
    /// \c .cu file.
    template <class PvtApproaches>
    static std::unique_ptr<Impl> makeImpl(const Problem& cpuProblem);

    std::unique_ptr<Impl> impl_;
};

} // namespace Opm::gpuistl

#endif // HAVE_CUDA

#endif // OPM_GPU_BLACKOIL_NOHYST_INTENSIVE_QUANTITIES_HPP
