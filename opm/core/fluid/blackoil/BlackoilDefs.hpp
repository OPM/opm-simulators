/*
  Copyright 2010 SINTEF ICT, Applied Mathematics.

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

#ifndef OPM_BLACKOILDEFS_HEADER_INCLUDED
#define OPM_BLACKOILDEFS_HEADER_INCLUDED


#include <tr1/array>
#include <iostream>
#include <boost/static_assert.hpp>

namespace Opm
{

    class BlackoilDefs
    {
    public:
        enum { numComponents = 3 };
        enum { numPhases = 3 };

        enum ComponentIndex { Water = 0, Oil = 1, Gas = 2 };
        enum PhaseIndex { Aqua = 0, Liquid = 1, Vapour = 2 };

        // We need types with operator= and constructor taking scalars
        // for the small vectors and matrices, to save us from having to
        // rewrite a large amount of code.
        template <typename T, int N>
        class SmallVec
        {
        public:
            SmallVec()
            {}
            SmallVec(const T& elem)
            { data_.assign(elem); } // In C++11, assign -> fill
            SmallVec& operator=(const T& elem)
            { data_.assign(elem); return *this; }
            const T& operator[](int i) const
            { return data_[i]; }
            T& operator[](int i)
            { return data_[i]; }
            template <typename U>
            void assign(const U& elem)
            {
                for (int i = 0; i < N; ++i) {
                    data_[i] = elem;
                }
            }
        private:
            std::tr1::array<T, N> data_;
        };
        template <typename T, int Rows, int Cols>
        class SmallMat
        {
        public:
            SmallMat()
            {}
            SmallMat(const T& elem)
            { data_.assign(elem); } // In C++11, assign -> fill
            SmallMat& operator=(const T& elem)
            { data_.assign(elem); return *this; }
            typedef SmallVec<T, Cols> RowType;
            const RowType& operator[](int i) const
            { return data_[i]; }
            RowType& operator[](int i)
            { return data_[i]; }
        private:
            SmallVec<RowType, Rows> data_;
        };

        typedef double Scalar;
        typedef SmallVec<Scalar, numComponents> CompVec;
        typedef SmallVec<Scalar, numPhases> PhaseVec;
        BOOST_STATIC_ASSERT(int(numComponents) == int(numPhases));
        typedef SmallMat<Scalar, numComponents, numPhases> PhaseToCompMatrix;
        typedef SmallMat<Scalar, numPhases, numPhases> PhaseJacobian;
        // Attempting to guard against alignment issues.
        BOOST_STATIC_ASSERT(sizeof(CompVec) == numComponents*sizeof(Scalar));
        BOOST_STATIC_ASSERT(sizeof(PhaseVec) == numPhases*sizeof(Scalar));
        BOOST_STATIC_ASSERT(sizeof(PhaseToCompMatrix) == numComponents*numPhases*sizeof(Scalar));
        BOOST_STATIC_ASSERT(sizeof(PhaseJacobian) == numPhases*numPhases*sizeof(Scalar));
    };

} // namespace Opm

#endif // OPM_BLACKOILDEFS_HEADER_INCLUDED
