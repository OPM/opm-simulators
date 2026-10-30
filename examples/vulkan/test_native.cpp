// SPDX-License-Identifier: GPL-3.0-or-later
#include "VulkanSolver.hpp"
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <stdexcept>

void
require(bool value, const char* message)
{
    if (!value)
        throw std::runtime_error(message);
}
int
main(int argc, char** argv)
{
    const char* selected = std::getenv("OPM_VULKAN_TEST_DEVICE");
    const std::string device = selected ? selected : "";
    const bool cpr = argc > 1 && std::string(argv[1]) == "--cpr";
    const bool sweeps = cpr || (argc > 1 && std::string(argv[1]) == "--dilu-sweeps");
    const bool dilu = sweeps || (argc > 1 && std::string(argv[1]) == "--dilu");
    Reslab::VulkanSolver solver(device, dilu, sweeps, cpr);
    std::cout << "Hardware device: " << solver.deviceName() << '\n';
    std::vector<uint32_t> rows {0, 2, 5, 7}, cols {0, 1, 0, 1, 2, 1, 2};
    std::vector<double> values {4, 1, 1, 3, 1, 1, 2}, x;
    solver.prepare(rows, cols, values);
    auto result = solver.solve({5, 5, 3}, x, 1e-10, 1);
    require(result.converged && result.iterations == 1, "Exact block solve failed");
    for (double value : x)
        require(std::abs(value - 1) < 1e-10, "Wrong solution");
    for (auto& value : values)
        value *= 2;
    solver.prepare(rows, cols, values);
    result = solver.solve({10, 10, 6}, x, 1e-10, 7);
    require(result.converged, "Matrix value update failed");
    result = solver.solve({0, 0, 0}, x, 1e-10, 200);
    require(result.converged && result.iterations == 0, "Zero RHS failed");
    rows = {0};
    cols.clear();
    values.clear();
    constexpr uint32_t n = 129;
    std::vector<double> rhs(n);
    for (uint32_t i = 0; i < n; ++i) {
        if (i >= 4) {
            cols.push_back(i - 4);
            values.push_back(-0.2);
        }
        if (i) {
            cols.push_back(i - 1);
            values.push_back(-1);
        }
        cols.push_back(i);
        values.push_back(4);
        if (i + 1 < n) {
            cols.push_back(i + 1);
            values.push_back(-0.5);
        }
        if (i + 4 < n) {
            cols.push_back(i + 4);
            values.push_back(-0.1);
        }
        rows.push_back(cols.size());
        rhs[i]
            = 4 - (i ? 1 : 0) - (i + 1 < n ? 0.5 : 0) - (i >= 4 ? 0.2 : 0) - (i + 4 < n ? 0.1 : 0);
    }
    solver.prepare(rows, cols, values);
    result = solver.solve(rhs, x, 1e-10, 200);
    require(result.converged && result.reduction <= 1e-10, "Nonsymmetric solve failed");
    result = solver.solve(rhs, x, 1e-14, 1);
    require(!result.converged && result.iterations == 1, "Iteration limit ignored");
    for (uint32_t i = 0; i < n; ++i) {
        const double scale = std::pow(10.0, int(i % 3) * 12 - 12);
        rhs[i] *= scale;
        for (uint32_t j = rows[i]; j < rows[i + 1]; ++j)
            values[j] *= scale;
    }
    solver.prepare(rows, cols, values);
    result = solver.solve(rhs, x, 1e-10, 200);
    require(result.converged, "Row equilibration failed");
    for (double value : x)
        require(std::abs(value - 1) < 1e-8, "Scaled system solution incorrect");
    bool invalid = false;
    try {
        solver.prepare({0, 5, 4, 3}, {0, 1, 2}, {1, 1, 1});
    } catch (const std::invalid_argument&) {
        invalid = true;
    }
    require(invalid, "Malformed CSR offsets accepted");
    invalid = false;
    try {
        solver.solve(rhs, x, 1e-10, 200);
    } catch (const std::invalid_argument&) {
        invalid = true;
    }
    require(invalid, "Solver used stale state after failed preparation");
    solver.prepare(rows, cols, values);
    require(solver.solve(rhs, x, 1e-10, 200).converged, "Recovery after failed preparation failed");
    // 129 independent cell pairs: each triangular level spans three workgroups,
    // including a partial group. Use coupled 3x3 blocks and a nonconstant solution.
    rows = {0};
    cols.clear();
    values.clear();
    constexpr uint32_t cells = 258;
    std::vector<double> exact(cells * 3);
    for (uint32_t i = 0; i < exact.size(); ++i)
        exact[i] = 0.5 + std::sin(i * 0.23);
    for (uint32_t cell = 0; cell < cells; ++cell)
        for (uint32_t r = 0; r < 3; ++r) {
            const uint32_t low = cell % 129, high = low + 129;
            for (uint32_t neighbor : {low, high})
                for (uint32_t c = 0; c < 3; ++c) {
                    cols.push_back(3 * neighbor + c);
                    values.push_back(neighbor == cell ? (r == c ? 4.0 : 0.1)
                                                      : (r == c ? -0.4 : 0.0));
                }
            rows.push_back(cols.size());
        }
    rhs.assign(exact.size(), 0);
    for (uint32_t i = 0; i < rhs.size(); ++i)
        for (uint32_t j = rows[i]; j < rows[i + 1]; ++j)
            rhs[i] += values[j] * exact[cols[j]];
    solver.prepare(rows, cols, values);
    result = solver.solve(rhs, x, 1e-10, 200);
    require(result.converged, "Parallel level solve failed");
    for (size_t i = 0; i < x.size(); ++i)
        require(std::abs(x[i] - exact[i]) < 1e-8, "Parallel level solution incorrect");
    // A nine-color clique exercises genuinely truncated triangular sweeps.
    constexpr uint32_t denseN = 27;
    rows = {0};
    cols.clear();
    values.clear();
    rhs.assign(denseN, 0);
    exact.resize(denseN);
    for (uint32_t i = 0; i < denseN; ++i)
        exact[i] = std::sin(0.37 * i) + 0.4;
    for (uint32_t i = 0; i < denseN; ++i) {
        for (uint32_t j = 0; j < denseN; ++j) {
            double a = i == j ? 10.0 : (i / 3 == j / 3 ? 0.3 : (i % 3 == j % 3 ? -0.2 : 0.01));
            cols.push_back(j);
            values.push_back(a);
            rhs[i] += a * exact[j];
        }
        rows.push_back(cols.size());
    }
    solver.prepare(rows, cols, values);
    require(solver.solve(rhs, x, 1e-10, 200).converged, "Nine-color solve failed");
    for (size_t i = 0; i < x.size(); ++i)
        require(std::abs(x[i] - exact[i]) < 1e-8, "Nine-color solution incorrect");
    if (cpr) {
        // Connected 3D pressure problem: multiple aggregation levels, external
        // pressure restriction weights, nonconstant solution and weak storage.
        constexpr uint32_t side = 8, count = side * side * side;
        rows = {0};
        cols.clear();
        values.clear();
        rhs.assign(3 * count, 0);
        exact.resize(3 * count);
        std::vector<double> weights(3 * count, 0);
        for (uint32_t i = 0; i < exact.size(); ++i)
            exact[i] = 0.5 + std::sin(0.17 * i);
        for (uint32_t cell = 0; cell < count; ++cell) {
            weights[3 * cell + 1] = 1;
            std::vector<uint32_t> neighbors;
            for (uint32_t axis = 0, stride = 1; axis < 3; ++axis, stride *= side) {
                uint32_t coordinate = (cell / stride) % side;
                if (coordinate > 0)
                    neighbors.push_back(cell - stride);
                if (coordinate + 1 < side)
                    neighbors.push_back(cell + stride);
            }
            for (uint32_t r = 0; r < 3; ++r) {
                std::vector<std::pair<uint32_t, double>> entries {
                    {3 * cell + r, r == 1 ? double(neighbors.size()) + 0.02 : 3.0}};
                if (r == 1)
                    for (auto neighbor : neighbors)
                        entries.emplace_back(3 * neighbor + 1, -1.0);
                std::sort(entries.begin(), entries.end());
                for (auto [column, value] : entries) {
                    cols.push_back(column);
                    values.push_back(value);
                    rhs[3 * cell + r] += value * exact[column];
                }
                rows.push_back(cols.size());
            }
        }
        solver.prepare(rows, cols, values, weights);
        result = solver.solve(rhs, x, 1e-10, 500);
        require(result.converged, "Multilevel 3D pressure solve failed");
        for (size_t i = 0; i < x.size(); ++i)
            require(std::abs(x[i] - exact[i]) < 1e-7, "Multilevel pressure solution incorrect");
        // Same sparsity, different matrix and restriction: cached coarse values
        // and the dense coarse inverse must both be updated.
        for (auto& a : values)
            a *= 1.7;
        for (auto& b : rhs)
            b *= 1.7;
        for (auto& a : weights)
            a *= 0.3;
        solver.prepare(rows, cols, values, weights);
        require(solver.solve(rhs, x, 1e-10, 500).converged,
                "Cached pressure hierarchy update failed");
        for (size_t i = 0; i < x.size(); ++i)
            require(std::abs(x[i] - exact[i]) < 1e-7, "Cached pressure solution incorrect");
        // Anisotropy and a weak interface exercise strength-aware aggregation,
        // rather than just the isotropic graph above.
        rhs.assign(exact.size(), 0);
        for (uint32_t i = 0; i < exact.size(); ++i) {
            double diagonal = 0.01;
            uint32_t diagonalEntry = UINT32_MAX;
            for (uint32_t k = rows[i]; k < rows[i + 1]; ++k) {
                if (cols[k] == i) {
                    diagonalEntry = k;
                    continue;
                }
                uint32_t a = i / 3, b = cols[k] / 3, difference = a > b ? a - b : b - a;
                double coupling = difference == 1 ? 1.0 : (difference == side ? 0.1 : 0.001);
                if ((a % side < side / 2) != (b % side < side / 2))
                    coupling *= 1e-4;
                values[k] = -coupling;
                diagonal += coupling;
            }
            if (i % 3 == 1)
                values[diagonalEntry] = diagonal;
            for (uint32_t k = rows[i]; k < rows[i + 1]; ++k)
                rhs[i] += values[k] * exact[cols[k]];
        }
        Reslab::VulkanSolver heterogeneous(device, true, true, true);
        heterogeneous.prepare(rows, cols, values, weights);
        require(heterogeneous.solve(rhs, x, 1e-12, 500).converged,
                "Anisotropic pressure solve failed");
        for (size_t i = 0; i < x.size(); ++i)
            require(std::abs(x[i] - exact[i]) < 1e-7, "Anisotropic pressure solution incorrect");
        // The symbolic matrix stays dense while numerical edges are activated.
        // All affected DILU and AMG coefficients must update through the cache.
        constexpr uint32_t activationN = 12;
        std::vector<uint32_t> ar {0}, ac;
        std::vector<double> av, ab(activationN), ae(activationN), aw(activationN, 0);
        for (uint32_t i = 0; i < activationN; ++i) {
            ae[i] = 0.5 + std::sin(0.29 * i);
            aw[i] = i % 3 == 1 ? 1 : 0;
            for (uint32_t j = 0; j < activationN; ++j) {
                ac.push_back(j);
                av.push_back(i == j ? 4.0 : 0.0);
            }
            ar.push_back(ac.size());
            ab[i] = 4 * ae[i];
        }
        heterogeneous.prepare(ar, ac, av, aw);
        require(heterogeneous.solve(ab, x, 1e-11, 100).converged, "Zero-edge initial solve failed");
        for (uint32_t i = 0; i < activationN; ++i) {
            ab[i] = 0;
            for (uint32_t j = 0; j < activationN; ++j) {
                if (i + 3 == j || j + 3 == i)
                    av[i * activationN + j] = -0.4;
                ab[i] += av[i * activationN + j] * ae[j];
            }
        }
        heterogeneous.prepare(ar, ac, av, aw);
        require(heterogeneous.solve(ab, x, 1e-11, 100).converged,
                "Activated symbolic edges failed");
        for (size_t i = 0; i < x.size(); ++i)
            require(std::abs(x[i] - ae[i]) < 1e-8, "Stale symbolic edge coefficients");
        weights[0] = NAN;
        bool invalidWeights = false;
        try {
            solver.prepare(rows, cols, values, weights);
        } catch (const std::invalid_argument&) {
            invalidWeights = true;
        }
        require(invalidWeights, "Nonfinite pressure weights accepted");
    }
    bool rejected = false;
    try {
        Reslab::VulkanSolver software("llvmpipe");
    } catch (const std::runtime_error&) {
        rejected = true;
    }
    require(rejected, "Software renderer accepted");
    std::cout << "Native Vulkan tests passed\n";
}
