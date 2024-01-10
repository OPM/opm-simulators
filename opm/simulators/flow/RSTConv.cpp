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

#include <opm/input/eclipse/Schedule/Schedule.hpp>
#include <opm/input/eclipse/Schedule/RSTConfig.hpp>

#include <dune/common/fvector.hh>
#include <dune/istl/bvector.hh>

#include <algorithm>
#include <cmath>
#include <numeric>

namespace Opm {

void RSTConv::init(const std::size_t numCells,
                   const RSTConfig& rst_config,
                   const std::array<int,3>& compIdx)
{
    const auto kw = rst_config.keywords.find("CONV");
    if (kw == rst_config.keywords.end()) {
        N_ = 0;
        cnv_X_.clear();
        return;
    }

    N_ = kw->second;
    compIdx_ = compIdx;

    cnv_X_.resize(3);
    for (std::size_t i = 0; i < 3; ++i) {
        if (compIdx_[i] != -1) {
            cnv_X_[i].resize(numCells);
        }
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

    std::vector<double> values(comm_.rank() == 0 ? comm_.size() * N_ : N_);
    std::vector<int> gIdx(comm_.rank() == 0 ? comm_.size() * N_ : N_);
    for (int i = 0; i < N_; ++i) {
        values[i] = std::abs(resid[lIdx[i]][compIdx_[comp]]);
        gIdx[i] = globalCell_[lIdx[i]];
    }

    comm_.gather(gIdx.data(), gIdx.data(), N_, 0);
    comm_.gather(values.data(), values.data(), N_, 0);
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

template<std::size_t Size>
using BFV = Dune::BlockVector<Dune::FieldVector<double,Size>>;

#define INSTANCE(SIZE) \
    template void RSTConv::update<BFV<SIZE>>(const BFV<SIZE>&); \
    template void RSTConv::gatherAndAccumulate<BFV<SIZE>>(const std::vector<int>&, \
                                                          const BFV<SIZE>&, int);

INSTANCE(1)
INSTANCE(2)
INSTANCE(3)
INSTANCE(4)
INSTANCE(5)
INSTANCE(6)

} // namespace Opm
