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

#include "config.h"
#include "dune/common/fvector.hh"
#include "dune/istl/bvector.hh"
#include <dune/common/fmatrix.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <opm/simulators/linalg/gpuistl/GpuSparseMatrix.hpp>
#include <opm/simulators/linalg/gpuistl/GpuVector.hpp>
#include <opm/simulators/linalg/gpuistl/detail/ISTLSolverGPUISTLImplementation.hpp>
#include <opm/simulators/linalg/matrixblock.hh>



namespace Opm::gpuistl::detail
{
template <class CpuMatrix, class CpuVector, class Comm>
ISTLSolverGPUISTLImplementation<CpuMatrix, CpuVector, Comm>::ISTLSolverGPUISTLImplementation(
    const FlowLinearSolverParameters& parameters,
    const PropertyTree& propertyTree,
    bool forceSerial,
    std::size_t pressureIndex)
    : m_parameters(parameters)
    , m_propertyTree(propertyTree)
    , m_forceSerial(forceSerial)
    , m_pressureIndex(pressureIndex)
{
}

template <class CpuMatrix, class CpuVector, class Comm>
void
ISTLSolverGPUISTLImplementation<CpuMatrix, CpuVector, Comm>::prepare(const CpuMatrix& M, CpuVector& b)
{
    updateMatrix(M);
    updateRhs(b);
}

template <class CpuMatrix, class CpuVector, class Comm>
void
ISTLSolverGPUISTLImplementation<CpuMatrix, CpuVector, Comm>::apply(CpuVector& x, Dune::InverseOperatorResult& result)
{
    if (!m_matrix) {
        OPM_THROW(std::runtime_error, "m_matrix not initialized, prepare(matrix, rhs); needs to be called");
    }
    if (!m_rhs) {
        OPM_THROW(std::runtime_error, "m_rhs not initialized, prepare(matrix, rhs); needs to be called");
    }
    if (!m_gpuSolver) {
        OPM_THROW(std::runtime_error, "m_gpuFlexibleSolver not initialized, prepare(matrix, rhs); needs to be called");
    }

    if (!m_x) {
        m_x = std::make_unique<GpuVector<real_type>>(x);
    } else {
        m_x->copyFromHost(x);
    }
    m_gpuSolver->apply(*m_x, *m_rhs, result);

    m_x->copyToHost(x);

    // The rest is handled "upstream"
}

template <class CpuMatrix, class CpuVector, class Comm>
void
ISTLSolverGPUISTLImplementation<CpuMatrix, CpuVector, Comm>::getResidual(CpuVector& b) const
{
    if (!m_rhs) {
        OPM_THROW(std::runtime_error, "m_rhs not initialized, prepare(matrix, rhs); needs to be called");
    }
    m_rhs->copyToHost(b);
}

template <class CpuMatrix, class CpuVector, class Comm>
void
ISTLSolverGPUISTLImplementation<CpuMatrix, CpuVector, Comm>::updateMatrix(const CpuMatrix& M)
{
    if (!m_matrix) {

        m_matrix.reset(new auto(GPUMatrix::fromMatrix(M)));
        std::function<GPUVector()> weightsCalculator = {};
        const bool parallel = false;
        m_gpuSolver = std::make_unique<SolverType>(
            *m_matrix, parallel, m_propertyTree, m_pressureIndex, weightsCalculator, m_forceSerial, nullptr);
    } else {
        m_matrix->updateNonzeroValues(M);
    }

    m_gpuSolver->update();
}

template <class CpuMatrix, class CpuVector, class Comm>
void
ISTLSolverGPUISTLImplementation<CpuMatrix, CpuVector, Comm>::updateRhs(const CpuVector& b)
{
    if (!m_rhs) {
        m_rhs = std::make_unique<GPUVector>(b);
    } else {
        m_rhs->copyFromHost(b);
    }
}

template <class CpuMatrix, class CpuVector, class Comm>
ISTLSolverGPUISTLImplementation<CpuMatrix, CpuVector, Comm>::~ISTLSolverGPUISTLImplementation() = default;

} // namespace Opm::gpuistl::detail


#if HAVE_MPI
// Parallel communicator and operators.
using Comm = Dune::OwnerOverlapCopyCommunication<int, int>;
#else
using Comm = Dune::Communication<int>;
#endif


#define INSTANTIATE_ISTLSOLVER_GPUISTL_DUNE(realtype, blockdim)                                                        \
    template class ::Opm::gpuistl::detail::ISTLSolverGPUISTLImplementation<                                            \
        ::Dune::BCRSMatrix<Dune::FieldMatrix<realtype, blockdim, blockdim>>,                                           \
        ::Dune::BlockVector<::Dune::FieldVector<realtype, blockdim>>,                                                  \
        Comm>;                                                                                                         \
    template class ::Opm::gpuistl::detail::ISTLSolverGPUISTLImplementation<                                            \
        ::Dune::BCRSMatrix<Opm::MatrixBlock<realtype, blockdim, blockdim>>,                                            \
        ::Dune::BlockVector<::Dune::FieldVector<realtype, blockdim>>,                                                  \
        Comm>


INSTANTIATE_ISTLSOLVER_GPUISTL_DUNE(double, 1);
INSTANTIATE_ISTLSOLVER_GPUISTL_DUNE(double, 2);
INSTANTIATE_ISTLSOLVER_GPUISTL_DUNE(double, 3);
INSTANTIATE_ISTLSOLVER_GPUISTL_DUNE(double, 4);
INSTANTIATE_ISTLSOLVER_GPUISTL_DUNE(double, 5);
INSTANTIATE_ISTLSOLVER_GPUISTL_DUNE(double, 6);

#if FLOW_INSTANTIATE_FLOAT
INSTANTIATE_ISTLSOLVER_GPUISTL_DUNE(float, 1);
INSTANTIATE_ISTLSOLVER_GPUISTL_DUNE(float, 2);
INSTANTIATE_ISTLSOLVER_GPUISTL_DUNE(float, 3);
INSTANTIATE_ISTLSOLVER_GPUISTL_DUNE(float, 4);
INSTANTIATE_ISTLSOLVER_GPUISTL_DUNE(float, 5);
INSTANTIATE_ISTLSOLVER_GPUISTL_DUNE(float, 6);
#endif
