/*
  Copyright 2025 Equinor ASA

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

#ifndef OPM_ABSTRACTISTLSOLVER_HEADER_INCLUDED
#define OPM_ABSTRACTISTLSOLVER_HEADER_INCLUDED


namespace Opm
{
template <class TypeTag>
class AbstractISTLSolver
{
public:

#if HAVE_MPI
    using CommunicationType = Dune::OwnerOverlapCopyCommunication<int,int>;
#else
    using CommunicationType = Dune::Communication<int>;
#endif

    using SparseMatrixAdapter = GetPropType<TypeTag, Properties::SparseMatrixAdapter>;
    using Vector = GetPropType<TypeTag, Properties::GlobalEqVector>;
    using Matrix = typename SparseMatrixAdapter::IstlMatrix;

    virtual ~AbstractISTLSolver() = default;
    
    virtual void eraseMatrix() = 0;

    virtual void setActiveSolver(int num) = 0;

    virtual int numAvailableSolvers() const = 0;

    virtual void prepare(const Matrix& M, Vector& b) = 0;

    virtual void prepare(const SparseMatrixAdapter& M, Vector& b) = 0;

    virtual void setResidual(Vector& b) = 0;

    virtual void getResidual(Vector& b) const = 0;

    virtual void setMatrix(const SparseMatrixAdapter& M) = 0;


    virtual bool solve(Vector& x) = 0;

    virtual int iterations() const = 0;

    virtual const CommunicationType* comm() const = 0;

    virtual int getSolveCount() const = 0;


};
} // namespace Opm

#endif
