/*
  Copyright 2024 SINTEF AS

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
#ifndef OPM_ILU0_KERNELS_HPP
#define OPM_ILU0_KERNELS_HPP
#include <cstddef>
#include <vector>
namespace Opm::cuistl::detail::ILU0
{
template <class T, int blocksize>
void ILUUpperSolveLevelSet(T* reorderedMat,
                               int* rowIndices,
                               int* colIndices,
                               int* indexConversion,
                               int startIdx,
                               int rowsInLevelSet,
                               T* v,
                               int threadBlockSize);


template <class T, int blocksize>
void ILULowerSolveLevelSet(T* reorderedMat,
                               int* rowIndices,
                               int* colIndices,
                               int* indexConversion,
                               int startIdx,
                               int rowsInLevelSet,
                               const T* d,
                               T* v,
                               int threadBlockSize);

template <class T, int blocksize>
void ILUUpperSolveLevelSetSplit(T* reorderedMat,
                               int* rowIndices,
                               int* colIndices,
                               int* indexConversion,
                               int startIdx,
                               int rowsInLevelSet,
                               const T* dInv,
                               T* v,
                               int threadBlockSize);


template <class T, int blocksize>
void ILULowerSolveLevelSetSplit(T* reorderedUpperMat,
                               int* rowIndices,
                               int* colIndices,
                               int* indexConversion,
                               int startIdx,
                               int rowsInLevelSet,
                               const T* dInv,
                               const T* d,
                               T* v,
                               int threadBlockSize);

template <class T, int blocksize>
void LUFactorization(
    T* srcMatrix, int* srcRowIndices, int* srcColumnIndices, int* naturalToReordered, int* reorderedToNatual, size_t rowsInLevelSet, int startIdx, int threadBlockSize);

template <class T, int blocksize>
void LUFactorizationSplit(T* reorderedLowerMat,
                         int* lowerRowIndices,
                         int* lowerColIndices,
                         T* reorderedUpperMat,
                         int* upperRowIndices,
                         int* upperColIndices,
                         T* diagonal,
                         int* reorderedToNatural,
                         int* naturalToReordered,
                         int startIdx,
                         int rowsInLevelSet,
                         int threadBlockSize);

} // namespace Opm::cuistl::detail
#endif
