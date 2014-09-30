/*
  Copyright 2014 SINTEF ICT, Applied Mathematics.

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

#ifndef OPM_DUNEMATRIX_HEADER_INCLUDED
#define OPM_DUNEMATRIX_HEADER_INCLUDED

#ifdef DUNE_BCRSMATRIX_HH
#error This header must be included before any bcrsmatrix.hh is included (directly or indirectly)
#endif

#include <opm/core/utility/platform_dependent/disable_warnings.h>

#include <dune/common/fmatrix.hh>
#include <dune/common/version.hh>

// Include matrix header with hackery to make it possible to inherit.
#define private protected
#include <dune/istl/bcrsmatrix.hh>
#undef private

#include <opm/core/utility/platform_dependent/reenable_warnings.h>

namespace Opm
{

    template <class B>
    class DuneMatrix : public Dune::BCRSMatrix<B>
    {
    public:
        DuneMatrix(const int rows, const int cols, const int* ia, const int* ja, const double* sa)
            // : build_mode(unknown), ready(built), n(rows), m(cols), nnz(ia[rows]),
            //   allocationSize(nnz), r(0), a(0),
            //   avg(0), overflowsize(-1.0)
        {
            typedef Dune::BCRSMatrix<B> Super;
            this->build_mode = Super::unknown;
            this->ready = Super::built;
            this->n = rows;
            this->m = cols;
            this->nnz = ia[rows];

#if DUNE_VERSION_NEWER(DUNE_ISTL, 2, 3)
            this->allocationSize = this->nnz;
            this->avg = 0;
            this->overflowsize = -1.0;
#endif

            this->a = new B[this->nnz];
            static_assert(sizeof(B) == sizeof(double), "This constructor requires a block type that is the same as a double.");
            std::copy(sa, sa + this->nnz, reinterpret_cast<double*>(this->a));
            this->j.reset(new typename Super::size_type[this->nnz]);
            std::copy(ja, ja +this-> nnz, this->j.get());
            this->r = new typename Super::row_type[rows];
            for (int row = 0; row < rows; ++row) {
                this->r[row].set(ia[row+1] - ia[row], this->a + ia[row], this->j.get() + ia[row]);
            }
        }
    };

} // namespace Opm

#endif // OPM_DUNEMATRIX_HEADER_INCLUDED
