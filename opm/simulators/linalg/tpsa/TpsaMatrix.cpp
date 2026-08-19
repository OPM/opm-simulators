// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  Copyright 2025 NORCE AS

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

#include <opm/simulators/linalg/tpsa/TpsaMatrix.hpp>

#include <opm/common/ErrorMacros.hpp>

#include <algorithm>

namespace Opm::Linear
{

//
// TpsaBlockRef implementation
//
template <class Scalar>
TpsaBlockRef<Scalar>&
TpsaBlockRef<Scalar>::operator+=(const MatrixBlock& b)
{
    apply_([](Scalar& stored, const Scalar& dense) {
               stored += dense;
           },
           b);

    return *this;
}

template <class Scalar>
TpsaBlockRef<Scalar>&
TpsaBlockRef<Scalar>::operator=(const MatrixBlock& b)
{
    apply_([](Scalar& stored, const Scalar& dense) {
               stored = dense;
           },
           b);

    return *this;
}

template <class Scalar>
TpsaBlockRef<Scalar>&
TpsaBlockRef<Scalar>::operator=(Scalar value)
{
    const MatrixBlock b(value);

    return *this = b;
}

template <class Scalar>
void
TpsaBlockRef<Scalar>::gather(MatrixBlock& b) const
{
    b = Scalar(0.0);
    apply_([](Scalar& stored, Scalar& dense) {
               dense = stored;
           },
           b);
}

//
// TpsaMatrix implementation
//
template <class Scalar>
void
TpsaMatrix<Scalar>::clear()
{
    // Return early, if called before reserve()
    if (nnz_ == 0) {
        return;
    }

    // Note: base_ points to start entry of sub-matrix thus suitable to use with std::fill_n
    forEachSubMatrixWithIndex_([&](auto& subMatrix, std::size_t subIdx) {
        std::fill_n(base_[subIdx], nnz_ * blockScalars_(subMatrix), Scalar(0.0));
    });
}

template <class Scalar>
void
TpsaMatrix<Scalar>::clearRow(const std::size_t row, const Scalar diag)
{
    const MatrixBlock zeroBlock(Scalar(0.0));

    MatrixBlock diagBlock(Scalar(0.0));
    for (int eqIdx = 0; eqIdx < numTpsaEq; ++eqIdx) {
        diagBlock[eqIdx][eqIdx] = diag;
    }

    for (std::size_t k = rowStart_[row]; k < rowStart_[row + 1]; ++k) {
        BlockAddress(*this, k) = (colIdx_[k] == row) ? diagBlock : zeroBlock;
    }
}

template <class Scalar>
void
TpsaMatrix<Scalar>::makeOverlapRowsInvalid(const std::vector<int>& overlapRows)
{
    for (const int row : overlapRows) {
        clearRow(static_cast<std::size_t>(row), Scalar(1.0));
    }
}

template <class Scalar>
void
TpsaMatrix<Scalar>::scaleFields(const std::array<Scalar, numTpsaFields>& rowFac,
                                const std::array<Scalar, numTpsaFields>& colFac)
{
    // Return early, if called before reserve()
    if (nnz_ == 0) {
        return;
    }

    forEachSubMatrixWithIndex_([&](auto& subMatrix, std::size_t subIdx) {
        const auto [rowField, colField] = subMatrixFields_[subIdx];
        const Scalar factor = rowFac[rowField] * colFac[colField];
        if (factor == Scalar(1.0)) {
            return;
        }

        Scalar* values = base_[subIdx];
        const std::size_t numValues = nnz_ * blockScalars_(subMatrix);
        for (std::size_t k = 0; k < numValues; ++k) {
            values[k] *= factor;
        }
    });
}

template <class Scalar>
void
TpsaMatrix<Scalar>::cacheValueArrays_()
{
    forEachSubMatrixWithIndex_([&](auto& subMatrix, std::size_t subIdx) {
        const std::size_t stride = blockScalars_(subMatrix);
        Scalar* base = nullptr;
        std::size_t k = 0;

        for (auto rowIt = subMatrix.begin(); rowIt != subMatrix.end(); ++rowIt) {
            for (auto colIt = rowIt->begin(); colIt != rowIt->end(); ++colIt, ++k) {
                Scalar* entry = &(*colIt)[0][0];
                if (k == 0) {
                    base = entry;
                }

                if (entry != base + k * stride || colIt.index() != colIdx_[k]) {
                    OPM_THROW(std::logic_error,
                              "TPSA: sub-matrix storage is not the contiguous, pattern-ordered "
                              "layout needed in TpsaMatrix");
                }
            }
        }

        if (k != nnz_) {
            OPM_THROW(std::logic_error,
                      "TPSA: sub-matrix has an unexpected number of nonzeroes");
        }

        base_[subIdx] = base;
    });
}

template <class Scalar>
void
TpsaMatrix<Scalar>::setMatrixView_()
{
    view_.M11_00 = &dd00_;
    view_.M12_00 = &dr0_;
    view_.M13_00 = &dsp0_;

    view_.M11_11 = &dd11_;
    view_.M12_10 = &dr1_;
    view_.M13_10 = &dsp1_;

    view_.M11_22 = &dd22_;
    view_.M12_20 = &dr2_;
    view_.M13_20 = &dsp2_;

    view_.M21_00 = &rd0_;
    view_.M21_01 = &rd1_;
    view_.M21_02 = &rd2_;
    view_.M22 = &rr_;
    view_.M23 = &rsp_;

    view_.M31_00 = &spd0_;
    view_.M31_01 = &spd1_;
    view_.M31_02 = &spd2_;
    view_.M32 = &spr_;
    view_.M33 = &spsp_;
}

template <class Scalar>
std::size_t
TpsaMatrix<Scalar>::flatIndex_(std::size_t rowIdx, std::size_t colIdx) const
{
    const auto begin = colIdx_.begin() + rowStart_[rowIdx];
    const auto end = colIdx_.begin() + rowStart_[rowIdx + 1];
    const auto it = std::lower_bound(begin, end, static_cast<unsigned>(colIdx));
    if (it == end || *it != static_cast<unsigned>(colIdx)) {
        OPM_THROW(std::logic_error,
                  "TPSA: requested a matrix block outside the sparsity pattern");
    }

    return rowStart_[rowIdx] + static_cast<std::size_t>(it - begin);
}

template class TpsaBlockRef<double>;
template class TpsaMatrix<double>;

#if FLOW_INSTANTIATE_FLOAT
template class TpsaBlockRef<float>;
template class TpsaMatrix<float>;
#endif

} // namespace Opm::Linear
