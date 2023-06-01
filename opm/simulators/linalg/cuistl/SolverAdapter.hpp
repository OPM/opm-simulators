/*
  Copyright 2022-2023 SINTEF AS

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
#ifndef OPM_SOLVERADAPTER_HPP
#define OPM_SOLVERADAPTER_HPP

#include <memory>

#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/operators.hh>
#include <dune/istl/schwarz.hh>
#include <dune/istl/solver.hh>
#include <opm/common/ErrorMacros.hpp>
#include <opm/simulators/linalg/cuistl/CuBlockPreconditioner.hpp>
#include <opm/simulators/linalg/cuistl/CuOwnerOverlapCopy.hpp>
#include <opm/simulators/linalg/cuistl/CuSparseMatrix.hpp>
#include <opm/simulators/linalg/cuistl/CuVector.hpp>
#include <opm/simulators/linalg/cuistl/PreconditionerAdapter.hpp>
#include <opm/simulators/linalg/cuistl/detail/has_function.hpp>



namespace Opm::cuistl
{
//! @brief Wraps a CUDA solver to work with CPU data.
//!
//! @tparam Operator the Dune::LinearOperator to use
//! @tparam UnderlyingSolver a Dune solver like class, eg Dune::BiCGSTABSolver
//! @tparam X the outer type to use (eg. Dune::BlockVector<Dune::FieldVector<...>>)
template <class Operator, template <class> class UnderlyingSolver, class X>
class SolverAdapter : public Dune::IterativeSolver<X, X>
{
public:
    using typename Dune::IterativeSolver<X, X>::domain_type;
    using typename Dune::IterativeSolver<X, X>::range_type;
    using typename Dune::IterativeSolver<X, X>::field_type;
    using typename Dune::IterativeSolver<X, X>::real_type;
    using typename Dune::IterativeSolver<X, X>::scalar_real_type;
    static constexpr auto block_size = domain_type::block_type::dimension;
    using XGPU = Opm::cuistl::CuVector<real_type>;

    // TODO: Use a std::forward
    SolverAdapter(Operator& op,
                  Dune::ScalarProduct<X>& sp,
                  std::shared_ptr<Dune::Preconditioner<X, X>> prec,
                  scalar_real_type reduction,
                  int maxit,
                  int verbose)
        : Dune::IterativeSolver<X, X>(op, sp, *prec, reduction, maxit, verbose)
        , m_opOnCPUWithMatrix(op)
        , m_matrix(CuSparseMatrix<real_type>::fromMatrix(op.getmat()))
        , m_underlyingSolver(constructSolver(prec, reduction, maxit, verbose))
    {
    }

    virtual void apply(X& x, X& b, double reduction, Dune::InverseOperatorResult& res) override
    {
        // TODO: Can we do this without reimplementing the other function?
        // TODO: [perf] Do we need to update the matrix every time? Probably yes
        m_matrix.updateNonzeroValues(m_opOnCPUWithMatrix.getmat());

        if (!m_inputBuffer) {
            m_inputBuffer.reset(new XGPU(b.dim()));
            m_outputBuffer.reset(new XGPU(x.dim()));
        }

        m_inputBuffer->copyFromHost(b);
        // TODO: [perf] do we need to copy x here?
        m_outputBuffer->copyFromHost(x);

        m_underlyingSolver.apply(*m_outputBuffer, *m_inputBuffer, reduction, res);

        // TODO: [perf] do we need to copy b here?
        m_inputBuffer->copyToHost(b);
        m_outputBuffer->copyToHost(x);
    }
    virtual void apply(X& x, X& b, Dune::InverseOperatorResult& res) override
    {
        // TODO: [perf] Do we need to update the matrix every time? Probably yes
        m_matrix.updateNonzeroValues(m_opOnCPUWithMatrix.getmat());

        if (!m_inputBuffer) {
            m_inputBuffer.reset(new XGPU(b.dim()));
            m_outputBuffer.reset(new XGPU(x.dim()));
        }

        m_inputBuffer->copyFromHost(b);
        // TODO: [perf] do we need to copy x here?
        m_outputBuffer->copyFromHost(x);

        m_underlyingSolver.apply(*m_outputBuffer, *m_inputBuffer, res);

        // TODO: [perf] do we need to copy b here?
        m_inputBuffer->copyToHost(b);
        m_outputBuffer->copyToHost(x);
    }

private:
    Operator& m_opOnCPUWithMatrix;
    CuSparseMatrix<real_type> m_matrix;

    UnderlyingSolver<XGPU> m_underlyingSolver;


    // TODO: Use a std::forward
    UnderlyingSolver<XGPU> constructSolver(std::shared_ptr<Dune::Preconditioner<X, X>> prec,
                                           scalar_real_type reduction,
                                           int maxit,
                                           int verbose)
    {

        OPM_ERROR_IF(
            detail::is_a_well_operator<Operator>::value,
            "Currently we only support operators of type MatrixAdapter in the CUDA solver. "
            "Use --matrix-add-well-contributions=true. "
            "Using WellModelMatrixAdapter with SolverAdapter is not well-defined so it will not work well, or at all.");





        if constexpr (detail::has_communication<Operator>::value) {
            // TODO: See the below TODO over the definition of precHolder in the other branch
            // TODO: We are currently double wrapping preconditioners in the preconditioner factory to be extra
            // compatible
            //       with CPU. Probably a cleaner way eventually would be to do more modifications to the flexible
            //       solver to accomodate the pure GPU better.
            auto precAsHolder = std::dynamic_pointer_cast<PreconditionerHolder<X, X>>(prec);
            if (!precAsHolder) {
                OPM_THROW(std::invalid_argument,
                          "The preconditioner needs to be a CUDA preconditioner (eg. CuILU0) wrapped in a "
                          "Opm::cuistl::PreconditionerAdapter wrapped in a "
                          "Opm::cuistl::CuBlockPreconditioner. If you are unsure what this means, set "
                          "preconditioner to 'CUILU0'"); // TODO: Suggest a better preconditioner
            }

            auto preconditionerAdapter = precAsHolder->getUnderlyingPreconditioner();
            auto preconditionerAdapterAsHolder
                = std::dynamic_pointer_cast<PreconditionerHolder<XGPU, XGPU>>(preconditionerAdapter);
            if (!preconditionerAdapterAsHolder) {
                OPM_THROW(std::invalid_argument,
                          "The preconditioner needs to be a CUDA preconditioner (eg. CuILU0) wrapped in a "
                          "Opm::cuistl::PreconditionerAdapter wrapped in a "
                          "Opm::cuistl::CuBlockPreconditioner. If you are unsure what this means, set "
                          "preconditioner to 'CUILU0'"); // TODO: Suggest a better preconditioner
            }
            // We need to get the underlying preconditioner:
            auto preconditionerReallyOnGPU = preconditionerAdapterAsHolder->getUnderlyingPreconditioner();
            const auto& communication = m_opOnCPUWithMatrix.getCommunication();

            using CudaCommunication = CuOwnerOverlapCopy<real_type, block_size, typename Operator::communication_type>;
            using SchwarzOperator
                = Dune::OverlappingSchwarzOperator<CuSparseMatrix<real_type>, XGPU, XGPU, CudaCommunication>;
            auto cudaCommunication = std::make_shared<CudaCommunication>(communication);

            auto mpiPreconditioner = std::make_shared<CuBlockPreconditioner<XGPU, XGPU, CudaCommunication>>(
                preconditionerReallyOnGPU, cudaCommunication);

            auto scalarProduct = std::make_shared<Dune::ParallelScalarProduct<XGPU, CudaCommunication>>(
                cudaCommunication, m_opOnCPUWithMatrix.category());


            // NOTE: Ownsership of cudaCommunication is handled by mpiPreconditioner. However, just to make sure we
            // remember
            //       this, we add this check. So remember that we hold one count in this scope, and one in the
            //       CuBlockPreconditioner instance. We accomedate for the fact that it could be passed around in
            //       CuBlockPreconditioner, hence we do not test for != 2
            OPM_ERROR_IF(cudaCommunication.use_count() < 2, "Internal error. Shared pointer not owned properly.");
            auto overlappingCudaOperator = std::make_shared<SchwarzOperator>(m_matrix, *cudaCommunication);

            return UnderlyingSolver<XGPU>(
                overlappingCudaOperator, scalarProduct, mpiPreconditioner, reduction, maxit, verbose);
        } else {
            // TODO: Fix the reliance on casting here. This is a code smell to a certain degree, and assumes
            //       a certain setup beforehand. The only reason we do it this way is to minimize edits to the
            //       flexible solver. We could design it differently, but keep this for the time being until
            //       we figure out how we want to GPU-ify the rest of the system.
            auto precAsHolder = std::dynamic_pointer_cast<PreconditionerHolder<XGPU, XGPU>>(prec);
            if (!precAsHolder) {
                OPM_THROW(std::invalid_argument,
                          "The preconditioner needs to be a CUDA preconditioner wrapped in a "
                          "Opm::cuistl::PreconditionerHolder (eg. CuILU0).");
            }
            auto preconditionerOnGPU = precAsHolder->getUnderlyingPreconditioner();

            auto matrixOperator
                = std::make_shared<Dune::MatrixAdapter<CuSparseMatrix<real_type>, XGPU, XGPU>>(m_matrix);
            auto scalarProduct = std::make_shared<Dune::SeqScalarProduct<XGPU>>();
            return UnderlyingSolver<XGPU>(
                matrixOperator, scalarProduct, preconditionerOnGPU, reduction, maxit, verbose);
        }
    }

    std::unique_ptr<XGPU> m_inputBuffer;
    std::unique_ptr<XGPU> m_outputBuffer;
};
} // namespace Opm::cuistl

#endif
