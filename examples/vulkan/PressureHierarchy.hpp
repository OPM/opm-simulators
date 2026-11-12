// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

namespace Reslab
{
// Unsmoothed aggregation, piecewise-constant prolongation, Galerkin coarse
// operators. All preparation is CPU; the packed hierarchy is applied on Vulkan.
struct PressureHierarchy {
    double strength = 0;
    PressureHierarchy()
    {
        if (const char* value = std::getenv("RESLAB_VULKAN_AGGREGATION_STRENGTH")) {
            std::size_t end = 0;
            strength = std::stod(value, &end);
            if (end != std::string(value).size() || !std::isfinite(strength) || strength < 0
                || strength > 1)
                throw std::invalid_argument("Aggregation strength must be in [0,1]");
        }
    }
    struct Level {
        std::vector<uint32_t> rows, cols, aggregate, members, memberRows;
        std::vector<double> values;
        std::vector<uint32_t> coarseScatter;
        uint32_t meta = 0;
        uint32_t size() const
        {
            return rows.size() - 1;
        }
    };
    std::vector<Level> levels;
    std::vector<uint32_t> indices;
    std::vector<double> data;
    uint32_t weightsOffset = 0;
    bool rebuilt = true;
    bool valid = false;
    std::vector<uint32_t> sourceRows, sourceCols, scalarScatter;
    static std::vector<double> invert(std::vector<double> a, uint32_t n)
    {
        std::vector<double> inverse(size_t(n) * n, 0);
        for (uint32_t i = 0; i < n; ++i)
            inverse[size_t(i) * n + i] = 1;
        for (uint32_t k = 0; k < n; ++k) {
            uint32_t pivot = k;
            for (uint32_t i = k + 1; i < n; ++i)
                if (std::abs(a[size_t(i) * n + k]) > std::abs(a[size_t(pivot) * n + k]))
                    pivot = i;
            double d = a[size_t(pivot) * n + k];
            if (d == 0 || !std::isfinite(d))
                throw std::runtime_error("Singular pressure coarse matrix");
            for (uint32_t j = 0; j < n; ++j) {
                std::swap(a[size_t(k) * n + j], a[size_t(pivot) * n + j]);
                std::swap(inverse[size_t(k) * n + j], inverse[size_t(pivot) * n + j]);
            }
            for (uint32_t j = 0; j < n; ++j) {
                a[size_t(k) * n + j] /= d;
                inverse[size_t(k) * n + j] /= d;
            }
            for (uint32_t i = 0; i < n; ++i)
                if (i != k) {
                    double f = a[size_t(i) * n + k];
                    for (uint32_t j = 0; j < n; ++j) {
                        a[size_t(i) * n + j] -= f * a[size_t(k) * n + j];
                        inverse[size_t(i) * n + j] -= f * inverse[size_t(k) * n + j];
                    }
                }
        }
        return inverse;
    }
    void build(const std::vector<uint32_t>& rows,
               const std::vector<uint32_t>& cols,
               const std::vector<double>& values,
               const std::vector<double>& suppliedWeights,
               const std::vector<uint32_t>& symbolicRows,
               const std::vector<uint32_t>& symbolicCols)
    {
        rebuilt = true;
        // Saturation coefficients may appear/disappear without changing the
        // pressure graph. Rebuild only the scalar scatter in that case.
        if (valid && !suppliedWeights.empty() && rows.size() - 1 == 3 * levels.front().size()
            && (rows != sourceRows || cols != sourceCols)) {
            const auto& f = levels.front();
            std::vector<uint32_t> scatter(cols.size(), UINT32_MAX);
            bool fits = true;
            for (uint32_t i = 0; i + 1 < rows.size() && fits; ++i)
                for (uint32_t k = rows[i]; k < rows[i + 1]; ++k)
                    if (cols[k] % 3 == 1) {
                        auto first = f.cols.begin() + f.rows[i / 3],
                             last = f.cols.begin() + f.rows[i / 3 + 1];
                        auto found = std::lower_bound(first, last, cols[k] / 3);
                        if (found == last || *found != cols[k] / 3) {
                            fits = false;
                            break;
                        }
                        scatter[k] = found - f.cols.begin();
                    }
            if (fits) {
                sourceRows = rows;
                sourceCols = cols;
                scalarScatter = std::move(scatter);
            }
        }
        bool reuse = valid && !suppliedWeights.empty() && rows == sourceRows && cols == sourceCols
            && !levels.empty();
        valid = false;
        if (reuse) {
            rebuilt = false;
            auto& fine = levels.front();
            std::fill(fine.values.begin(), fine.values.end(), 0);
            for (uint32_t i = 0; i + 1 < rows.size(); ++i)
                for (uint32_t k = rows[i]; k < rows[i + 1]; ++k)
                    if (scalarScatter[k] != UINT32_MAX)
                        fine.values[scalarScatter[k]] += suppliedWeights[i] * values[k];
            for (size_t l = 0; l + 1 < levels.size(); ++l) {
                auto& f = levels[l];
                auto& c = levels[l + 1];
                std::fill(c.values.begin(), c.values.end(), 0);
                for (size_t k = 0; k < f.values.size(); ++k)
                    c.values[f.coarseScatter[k]] += f.values[k];
            }
            std::copy(suppliedWeights.begin(), suppliedWeights.end(), data.begin() + weightsOffset);
            for (auto& f : levels) {
                uint32_t m = f.meta, n = f.size();
                std::copy(f.values.begin(), f.values.end(), data.begin() + indices[m + 3]);
                for (uint32_t i = 0; i < n; ++i) {
                    double diagonal = 0;
                    for (uint32_t j = f.rows[i]; j < f.rows[i + 1]; ++j)
                        if (f.cols[j] == i)
                            diagonal = f.values[j];
                    if (diagonal == 0 || !std::isfinite(diagonal))
                        throw std::runtime_error("Invalid updated pressure diagonal");
                    data[indices[m + 4] + i] = 1 / diagonal;
                }
            }
            auto& c = levels.back();
            uint32_t n = c.size();
            std::vector<double> dense(size_t(n) * n, 0);
            for (uint32_t i = 0; i < n; ++i)
                for (uint32_t j = c.rows[i]; j < c.rows[i + 1]; ++j)
                    dense[size_t(i) * n + c.cols[j]] = c.values[j];
            auto inverse = invert(dense, n);
            std::copy(inverse.begin(), inverse.end(), data.begin() + indices[c.meta + 12]);
            valid = true;
            return;
        }
        levels.clear();
        indices.clear();
        data.clear();
        uint32_t nb = (rows.size() - 1) / 3;
        std::vector<double> weights(3 * nb);
        Level fine;
        fine.rows.push_back(0);
        for (uint32_t cell = 0; cell < nb; ++cell) {
            if (suppliedWeights.empty()) {
                std::vector<double> diagonal(9, 0);
                for (uint32_t r = 0; r < 3; ++r) {
                    for (uint32_t j = rows[3 * cell + r]; j < rows[3 * cell + r + 1]; ++j) {
                        if (cols[j] / 3 == cell)
                            diagonal[3 * r + cols[j] % 3] = values[j];
                    }
                }
                auto inv = invert(diagonal, 3);
                for (uint32_t r = 0; r < 3; ++r)
                    weights[3 * cell + r] = inv[3 + r];
            } else {
                // Preserve OPM's physical restriction. Normalizing every row by
                // its diagonal changes Galerkin aggregation weights and can
                // severely weaken coarse pressure corrections on heterogeneous grids.
                for (uint32_t r = 0; r < 3; ++r)
                    weights[3 * cell + r] = suppliedWeights[3 * cell + r];
            }
            std::map<uint32_t, double> entries;
            // Keep the symbolic graph even when all pressure coefficients of a
            // connection are temporarily zero. Later activation updates values
            // instead of rebuilding every level of AMG.
            for (uint32_t k = symbolicRows[cell]; k < symbolicRows[cell + 1]; ++k)
                entries[symbolicCols[k]] = 0;
            for (uint32_t r = 0; r < 3; ++r)
                for (uint32_t j = rows[3 * cell + r]; j < rows[3 * cell + r + 1]; ++j)
                    if (cols[j] % 3 == 1)
                        entries[cols[j] / 3] += weights[3 * cell + r] * values[j];
            for (auto [c, v] : entries) {
                fine.cols.push_back(c);
                fine.values.push_back(v);
            }
            fine.rows.push_back(fine.cols.size());
        }
        levels.push_back(std::move(fine));
        while (levels.back().size() > 48) {
            auto& f = levels.back();
            uint32_t n = f.size(), nc = 0;
            std::vector<double> rowStrength(n, 0);
            for (uint32_t i = 0; i < n; ++i)
                for (uint32_t j = f.rows[i]; j < f.rows[i + 1]; ++j)
                    if (f.cols[j] != i)
                        rowStrength[i] = std::max(rowStrength[i], std::abs(f.values[j]));
            f.aggregate.assign(n, UINT32_MAX);
            for (uint32_t i = 0; i < n; ++i)
                if (f.aggregate[i] == UINT32_MAX) {
                    std::vector<std::pair<double, uint32_t>> neighbors;
                    for (uint32_t j = f.rows[i]; j < f.rows[i + 1]; ++j)
                        if (f.cols[j] != i && f.values[j] != 0
                            && std::abs(f.values[j]) >= strength * rowStrength[i]
                            && f.aggregate[f.cols[j]] == UINT32_MAX)
                            neighbors.emplace_back(std::abs(f.values[j]), f.cols[j]);
                    if (neighbors.empty())
                        continue;
                    f.aggregate[i] = nc;
                    std::sort(neighbors.rbegin(), neighbors.rend());
                    unsigned added = 0;
                    for (auto [strength, c] : neighbors)
                        if (strength > 0 && f.aggregate[c] == UINT32_MAX && added < 3) {
                            f.aggregate[c] = nc;
                            ++added;
                        }
                    ++nc;
                }
            // Attach leftovers to an existing strong neighbor. Creating a new
            // singleton for every leaf makes well/star graphs coarsen by only
            // a few vertices per level and produces dozens of GPU V-cycle levels.
            for (uint32_t i = 0; i < n; ++i)
                if (f.aggregate[i] == UINT32_MAX) {
                    double strongest = 0;
                    uint32_t target = UINT32_MAX;
                    for (uint32_t j = f.rows[i]; j < f.rows[i + 1]; ++j)
                        if (f.cols[j] != i && f.aggregate[f.cols[j]] != UINT32_MAX
                            && std::abs(f.values[j]) >= strength * rowStrength[i]
                            && std::abs(f.values[j]) > strongest) {
                            strongest = std::abs(f.values[j]);
                            target = f.aggregate[f.cols[j]];
                        }
                    f.aggregate[i] = target == UINT32_MAX ? nc++ : target;
                }
            if (nc >= n) {
                if (n <= 256)
                    break;
                throw std::runtime_error("Pressure aggregation cannot coarsen disconnected system");
            }
            f.memberRows.assign(nc + 1, 0);
            for (auto c : f.aggregate)
                ++f.memberRows[c + 1];
            for (uint32_t c = 0; c < nc; ++c)
                f.memberRows[c + 1] += f.memberRows[c];
            auto cursor = f.memberRows;
            f.members.resize(n);
            for (uint32_t i = 0; i < n; ++i)
                f.members[cursor[f.aggregate[i]]++] = i;
            std::vector<std::map<uint32_t, double>> entries(nc);
            for (uint32_t i = 0; i < n; ++i)
                for (uint32_t j = f.rows[i]; j < f.rows[i + 1]; ++j)
                    entries[f.aggregate[i]][f.aggregate[f.cols[j]]] += f.values[j];
            Level coarse;
            coarse.rows.push_back(0);
            for (auto& row : entries) {
                for (auto [c, v] : row) {
                    coarse.cols.push_back(c);
                    coarse.values.push_back(v);
                }
                coarse.rows.push_back(coarse.cols.size());
            }
            levels.push_back(std::move(coarse));
        }
        // Per-level metadata: n,row,col,val,invdiag,b,x,temp,aggregate,
        // memberRows,members,coarseMeta,denseInverse,weights. Offsets are elements.
        indices.resize(levels.size() * 14, 0);
        auto putI = [&](const auto& v) {
            uint32_t o = indices.size();
            indices.insert(indices.end(), v.begin(), v.end());
            return o;
        };
        auto putD = [&](const auto& v) {
            uint32_t o = data.size();
            data.insert(data.end(), v.begin(), v.end());
            return o;
        };
        weightsOffset = putD(weights);
        for (uint32_t l = 0; l < levels.size(); ++l) {
            auto& f = levels[l];
            uint32_t m = l * 14, n = f.size();
            f.meta = m;
            indices[m] = n;
            indices[m + 1] = putI(f.rows);
            indices[m + 2] = putI(f.cols);
            indices[m + 3] = putD(f.values);
            std::vector<double> diagonal(n);
            for (uint32_t i = 0; i < n; ++i) {
                for (uint32_t j = f.rows[i]; j < f.rows[i + 1]; ++j)
                    if (f.cols[j] == i)
                        diagonal[i] = f.values[j];
                if (diagonal[i] == 0 || !std::isfinite(diagonal[i]))
                    throw std::runtime_error("Invalid pressure diagonal");
                diagonal[i] = 1 / diagonal[i];
            }
            indices[m + 4] = putD(diagonal);
            for (uint32_t k = 5; k <= 7; ++k)
                indices[m + k] = putD(std::vector<double>(n, 0));
            indices[m + 8] = putI(f.aggregate);
            indices[m + 9] = putI(f.memberRows);
            indices[m + 10] = putI(f.members);
            indices[m + 11] = (l + 1) * 14;
            indices[m + 13] = weightsOffset;
            if (l + 1 == levels.size()) {
                std::vector<double> dense(size_t(n) * n, 0);
                for (uint32_t i = 0; i < n; ++i)
                    for (uint32_t j = f.rows[i]; j < f.rows[i + 1]; ++j)
                        dense[size_t(i) * n + f.cols[j]] = f.values[j];
                indices[m + 12] = putD(invert(dense, n));
            }
        }
        sourceRows = rows;
        sourceCols = cols;
        scalarScatter.assign(cols.size(), UINT32_MAX);
        const auto& fineLevel = levels.front();
        for (uint32_t i = 0; i + 1 < rows.size(); ++i)
            for (uint32_t k = rows[i]; k < rows[i + 1]; ++k)
                if (cols[k] % 3 == 1) {
                    auto first = fineLevel.cols.begin() + fineLevel.rows[i / 3],
                         last = fineLevel.cols.begin() + fineLevel.rows[i / 3 + 1];
                    scalarScatter[k]
                        = std::lower_bound(first, last, cols[k] / 3) - fineLevel.cols.begin();
                }
        for (size_t l = 0; l + 1 < levels.size(); ++l) {
            auto& f = levels[l];
            auto& c = levels[l + 1];
            f.coarseScatter.resize(f.values.size());
            for (uint32_t i = 0; i < f.size(); ++i)
                for (uint32_t k = f.rows[i]; k < f.rows[i + 1]; ++k) {
                    uint32_t cr = f.aggregate[i], cc = f.aggregate[f.cols[k]];
                    f.coarseScatter[k] = std::lower_bound(c.cols.begin() + c.rows[cr],
                                                          c.cols.begin() + c.rows[cr + 1],
                                                          cc)
                        - c.cols.begin();
                }
        }
        valid = true;
    }
};
} // namespace Reslab
