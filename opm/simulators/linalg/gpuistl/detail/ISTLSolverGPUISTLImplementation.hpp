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

#ifndef OPM_ISTLSOLVERGPUISTLIMPLEMENTATION_HEADER_INCLUDED
#define OPM_ISTLSOLVERGPUISTLIMPLEMENTATION_HEADER_INCLUDED

#include "dune/istl/solver.hh"
#include <memory>

#include <dune/istl/operators.hh>

#include <opm/simulators/linalg/FlowLinearSolverParameters.hpp>
#include <opm/simulators/linalg/PropertyTree.hpp>
#include <opm/simulators/linalg/gpuistl/detail/FlexibleSolverWrapper.hpp>


namespace Opm::gpuistl
{
// We need to forward declare these classes, but not include the headers,
// to be able to include this header in main files.
//
// TODO: Fix hipification so that we can include the headers here.
//       Once hipification is fixed, this class can be integrated into ISTLSolverGPUISTL.
template <class T>
class GpuSparseMatrix;

template <class T>
class GpuVector;
} // namespace Opm::gpuistl

namespace Opm::gpuistl::detail
{

/**
 * This is an implementation detail class for the ISTLSolverGPUISTL class.
 * It is not intended to be used directly.
 *
 * In short, this class just exists until we can make the hipification
 * a bit more robust.
 */
template <class CpuMatrix, class CpuVector, class Comm>
class ISTLSolverGPUISTLImplementation
{
public:
    using real_type = typename CpuVector::field_type;
    using GPUMatrix = Opm::gpuistl::GpuSparseMatrix<real_type>;
    using GPUVector = Opm::gpuistl::GpuVector<real_type>;

    using SolverType = Opm::gpuistl::detail::FlexibleSolverWrapper<GPUMatrix, GPUVector, Comm>;

    ISTLSolverGPUISTLImplementation(const FlowLinearSolverParameters& parameters,
                                    const PropertyTree& propertyTree,
                                    bool forceSerial,
                                    std::size_t pressureIndex);

    ISTLSolverGPUISTLImplementation() = delete;
    ISTLSolverGPUISTLImplementation(const ISTLSolverGPUISTLImplementation&) = delete;
    ISTLSolverGPUISTLImplementation& operator=(const ISTLSolverGPUISTLImplementation&) = delete;

    ~ISTLSolverGPUISTLImplementation();

    void prepare(const CpuMatrix& M, CpuVector& b);

    void apply(CpuVector& x, Dune::InverseOperatorResult& result);

    void getResidual(CpuVector& b) const;

private:
    void updateMatrix(const CpuMatrix& M);
    void updateRhs(const CpuVector& b);

    std::unique_ptr<GPUMatrix> m_matrix;

    std::unique_ptr<SolverType> m_gpuSolver;

    std::unique_ptr<GPUVector> m_rhs;
    std::unique_ptr<GPUVector> m_x;

    const FlowLinearSolverParameters& m_parameters;
    const PropertyTree& m_propertyTree;
    const bool m_forceSerial;
    const std::size_t m_pressureIndex;
};

} // namespace Opm::gpuistl::detail

#endif // OPM_ISTLSOLVERGPUISTLIMPLEMENTATION_HEADER_INCLUDED
