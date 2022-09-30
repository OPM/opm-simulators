/*
  Copyright SINTEF AS 2022

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
#ifndef OPM_CUISTL_VECTOR_OPERATIONS_HEADER
#define OPM_CUISTL_VECTOR_OPERATIONS_HEADER
namespace Opm::cuistl::impl
{

template <class T>
void setVectorValue(T* deviceData, size_t numberOfElements, const T& value);

template <class T>
void setZeroAtIndexSet(T* deviceData, size_t numberOfElements, const int* indices);

template <class T>
T innerProductAtIndices(const T* deviceA, const T* deviceB, T* buffer, size_t numberOfElements, const int* indices);
} // namespace Opm::cuistl::impl
#endif
