/*
  Copyright 2026 Equinor ASA
  This file is part of OPM, distributed under the GNU General Public License,
  version 3 or (at your option) any later version.
*/
#ifndef OPM_SCOPED_GPU_MEMORY_ACCOUNTING_HPP
#define OPM_SCOPED_GPU_MEMORY_ACCOUNTING_HPP

#include <cstddef>
#include <cstdint>

namespace Opm::gpuistl::detail {

struct GpuMemoryCounts {
    std::uint64_t allocations{0};
    std::uint64_t allocationBytes{0};
    std::uint64_t hostToDeviceCalls{0};
    std::uint64_t hostToDeviceBytes{0};
};

// Accounting is opt-in at an owner-construction boundary. Nested scopes
// count both the whole owner and a selected static-data subtree. Unrelated
// threads and ordinary solver allocations have no active accounting scope.
class ScopedGpuMemoryAccounting
{
public:
    explicit ScopedGpuMemoryAccounting(GpuMemoryCounts& counts)
        : counts_(counts), parent_(active_)
    { active_ = this; }

    ~ScopedGpuMemoryAccounting() { active_ = parent_; }
    ScopedGpuMemoryAccounting(const ScopedGpuMemoryAccounting&) = delete;
    ScopedGpuMemoryAccounting& operator=(const ScopedGpuMemoryAccounting&) = delete;

    static void allocation(std::size_t bytes)
    {
        for (auto* scope = active_; scope; scope = scope->parent_) {
            ++scope->counts_.allocations;
            scope->counts_.allocationBytes += bytes;
        }
    }

    static void hostToDevice(std::size_t bytes)
    {
        for (auto* scope = active_; scope; scope = scope->parent_) {
            ++scope->counts_.hostToDeviceCalls;
            scope->counts_.hostToDeviceBytes += bytes;
        }
    }

private:
    GpuMemoryCounts& counts_;
    ScopedGpuMemoryAccounting* parent_;
    inline static thread_local ScopedGpuMemoryAccounting* active_{nullptr};
};

} // namespace Opm::gpuistl::detail
#endif
