/*
  Copyright 2023 Equinor ASA.

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
#include <opm/simulators/flow/RSTConv.hpp>

#include <dune/common/fvector.hh>
#include <dune/istl/bvector.hh>

#include <opm/input/eclipse/Schedule/RSTConfig.hpp>
#include <opm/input/eclipse/Schedule/Schedule.hpp>

#include <opm/output/data/Solution.hpp>

#include <algorithm>
#include <cmath>
#include <numeric>

namespace Opm {

void RSTConv::init(const std::size_t numCells,
                   const RSTConfig& rst_config,
                   const std::array<int,6>& compIdx)
{
    const auto kw = rst_config.keywords.find("CONV");
    if (kw == rst_config.keywords.end()) {
        N_ = 0;
        cnv_X_.clear();
        conv_new_.clear();
        return;
    }

    N_ = kw->second;
    compIdx_ = compIdx;

    cnv_X_.resize(6);
    for (std::size_t i = 0; i < 6; ++i) {
        if (compIdx_[i] > -1) {
            cnv_X_[i].resize(numCells);
        }
    }
    conv_new_.resize(numCells, 1);
}

void RSTConv::outputRestart(data::Solution& sol)
{
    if (!this->cnv_X_.empty()) {
        constexpr const std::array names{"CNV_OIL", "CNV_GAS", "CNV_WAT",
                                         "CNV_PLY", "CNV_SAL", "CNV_SOL"};
        std::ranges::for_each(this->cnv_X_,
                              [i = 0, &names, &sol](auto& cnv) mutable
                              {
                                  if (!cnv.empty()) {
                                      sol.insert(names[i], std::move(cnv), data::TargetType::RESTART_SOLUTION);
                                  }
                                  ++i;
                              });
        sol.insert("CONV_NEW", std::move(conv_new_), data::TargetType::RESTART_SOLUTION);
    }
}

bool RSTConv::hasConv() const
{
    return !conv_new_.empty();
}

void RSTConv::prepareConv()
{
    std::fill(conv_new_.begin(), conv_new_.end(), 1);
}

void RSTConv::updateNewton(const std::vector<int>& convNewt)
{
    const std::size_t numGloCells = conv_new_.size();

    if (comm_.size() == 1) {
        for (std::size_t n = 0; n < numGloCells; ++n) {
            conv_new_[n] += convNewt[n];
        }
        return;
    }

    std::vector<int> values(numGloCells, 0);
    const std::size_t numLocCells = convNewt.size();

    for (std::size_t i = 0; i < numLocCells; ++i) {
        values[globalCell_(i)] = convNewt[i];
    }

    comm_.sum(values.data(), numGloCells);

    for (std::size_t n = 0; n < numGloCells; ++n) {
        conv_new_[n] += values[n];
    }
}

template<class ResidualVector>
void RSTConv::update(const ResidualVector& resid)
{
    if (N_ == 0) {
        return;
    }

    std::vector<int> cix(resid.size());
    std::iota(cix.begin(), cix.end(), 0);
    for (std::size_t i = 0; i < cnv_X_.size(); ++i) {
        if (cnv_X_[i].empty()) {
            continue;
        }
        std::partial_sort(cix.begin(), cix.begin() + N_, cix.end(),
                          [&resid,c = compIdx_[i]](const int a, const int b)
                          { return std::abs(resid[a][c]) > std::abs(resid[b][c]); });

        this->gatherAndAccumulate(cix, resid, i);
    }
}

template<class ResidualVector>
void RSTConv::gatherAndAccumulate(const std::vector<int>& lIdx,
                                  const ResidualVector& resid, int comp)
{
    if (comm_.size() == 1) {
        for (int n = 0; n < N_; ++n) {
            ++cnv_X_[comp][lIdx[n]];
        }
        return;
    }

    std::vector<double> send_values(N_);
    std::vector<double> values(comm_.rank() == 0 ? comm_.size() * N_ : 0);
    std::vector<int> send_idx(N_);
    std::vector<int> gIdx(comm_.rank() == 0 ? comm_.size() * N_ : 0);
    for (int i = 0; i < N_; ++i) {
        send_values[i] = std::abs(resid[lIdx[i]][compIdx_[comp]]);
        send_idx[i] = globalCell_(lIdx[i]);
    }

    comm_.gather(send_idx.data(), gIdx.data(), N_, 0);
    comm_.gather(send_values.data(), values.data(), N_, 0);
    if (comm_.rank() != 0) {
        return;
    }

    std::vector<int> valIdx(values.size());
    std::iota(valIdx.begin(), valIdx.end(), 0);

    std::partial_sort(valIdx.begin(), valIdx.begin() + N_, valIdx.end(),
                      [&values](const int i1, const int i2)
                      { return values[i1] > values[i2]; });

    for (int n = 0; n < N_; ++n) {
        ++cnv_X_[comp][gIdx[valIdx[n]]];
    }
}

template<class Scalar, std::size_t Size>
using BFV = Dune::BlockVector<Dune::FieldVector<Scalar,Size>>;

#define INSTANTIATE(T,SIZE)                                              \
    template void RSTConv::update(const BFV<T,SIZE>&);                   \
    template void RSTConv::gatherAndAccumulate(const std::vector<int>&,  \
                                               const BFV<T,SIZE>&, int);

#define INSTANTIATE_TYPE(T) \
    INSTANTIATE(T,1)        \
    INSTANTIATE(T,2)        \
    INSTANTIATE(T,3)        \
    INSTANTIATE(T,4)        \
    INSTANTIATE(T,5)        \
    INSTANTIATE(T,6)

INSTANTIATE_TYPE(double)

#if FLOW_INSTANTIATE_FLOAT
INSTANTIATE_TYPE(float)
#endif

} // namespace Opm
