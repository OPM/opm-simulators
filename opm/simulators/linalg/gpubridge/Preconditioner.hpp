/*
  Copyright 2024 Equinor ASA

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

#ifndef OPM_PRECONDITIONER_HEADER_INCLUDED
#define OPM_PRECONDITIONER_HEADER_INCLUDED

#if HAVE_OPENCL
#include <opm/simulators/linalg/gpubridge/opencl/opencl.hpp>
#endif

#include <opm/simulators/linalg/gpubridge/WellContributions.hpp>

#include <memory>

namespace Opm::Accelerator {

enum PreconditionerType {
    BILU0,
    CPR,
    BISAI
};

template<class Scalar> class BlockedMatrix;

template<class Scalar, unsigned int block_size, class ApplyScalar = Scalar>
class Preconditioner
{
protected:
    int N = 0;       // number of rows of the matrix
    int Nb = 0;      // number of blockrows of the matrix
    int nnz = 0;     // number of nonzeroes of the matrix (scalar)
    int nnzb = 0;    // number of blocks of the matrix
    int verbosity = 0;

    Preconditioner(int verbosity_)
        : verbosity(verbosity_)
    {}

public:

    virtual ~Preconditioner() = default;

    static std::unique_ptr<Preconditioner> create(PreconditionerType type,
                                                  bool opencl_ilu_parallel,
                                                  int verbosity);

    // apply preconditioner, x = prec(y)
    virtual void apply(const ApplyScalar& y,
                       ApplyScalar& x,
                       WellContributions<Scalar>& wellContribs) = 0;

    // analyze matrix, e.g. the sparsity pattern
    // probably only called once
    // the version with two params can be overloaded, if not, it will default to using the one param version
    virtual bool analyze_matrix(BlockedMatrix<Scalar>* mat) = 0;
    virtual bool analyze_matrix(BlockedMatrix<Scalar>* mat,
                                BlockedMatrix<Scalar>* jacMat) = 0;

    // create/update preconditioner, probably used every linear solve
    // the version with two params can be overloaded, if not, it will default to using the one param version
    virtual bool create_preconditioner(BlockedMatrix<Scalar>* mat) = 0;
    virtual bool create_preconditioner(BlockedMatrix<Scalar>* mat,
                                       BlockedMatrix<Scalar>* jacMat) = 0;
};

} // namespace Opm::Accelerator

#endif
