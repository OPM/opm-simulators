// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

namespace Reslab
{
struct SolveResult {
    bool converged = false;
    int iterations = 0;
    double reduction = 1;
    double seconds = 0;
};

// Sequential FP64 solver. Host assembly/setup, resident GPU iteration, true residual check.
class VulkanSolver
{
public:
    explicit VulkanSolver(const std::string& deviceMatch,
                          bool useDilu = false,
                          bool parallelSweeps = false,
                          bool useCpr = false);
    ~VulkanSolver();
    VulkanSolver(const VulkanSolver&) = delete;
    VulkanSolver& operator=(const VulkanSolver&) = delete;
    void prepare(std::vector<std::uint32_t> rows,
                 std::vector<std::uint32_t> columns,
                 std::vector<double> values,
                 const std::vector<double>& pressureWeights = {});
    SolveResult solve(const std::vector<double>& rhs,
                      std::vector<double>& x,
                      double tolerance,
                      int maxIterations);
    const std::string& deviceName() const;

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};
} // namespace Reslab
