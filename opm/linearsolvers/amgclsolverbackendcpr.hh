// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 2 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.

  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the list of
  copyright holders.
*/
/*!
 * \file
 * \copydoc Ewoms::Linear::SuperLUBackend
 */
#ifndef EWOMS_AMGCLSOLVER_BACKEND_CPR_HH
#define EWOMS_AMGCLSOLVER_BACKEND_CPR_HH

#include <opm/linearsolvers/amgclsolverbackend.hh>
#include <opm/core/linalg/LinearSolverAmgclCpr.hpp>
BEGIN_PROPERTIES

END_PROPERTIES

namespace Ewoms {
    namespace Linear {
        template <class TypeTag>
        class AMGCLSolverBackendCpr :
            public AMGCLSolverBackend<TypeTag>
        {
            typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
            typedef typename GET_PROP_TYPE(TypeTag, FluidSystem)       FluidSystem;
            typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
            typedef typename GET_PROP_TYPE(TypeTag, GlobalEqVector) Vector;
            typedef typename GET_PROP_TYPE(TypeTag, SparseMatrixAdapter) SparseMatrixAdapter;
            typedef typename GET_PROP_TYPE(TypeTag, EclWellModel) WellModel;
            typedef typename GET_PROP_TYPE(TypeTag, Simulator) Simulator;
            typedef typename SparseMatrixAdapter::MatrixBlock MatrixBlockType;
            typedef typename SparseMatrixAdapter::IstlMatrix Matrix;
            typedef typename Vector::block_type BlockVector;
            typedef typename GET_PROP_TYPE(TypeTag, Evaluation) Evaluation;
            typedef typename GET_PROP_TYPE(TypeTag, ThreadManager) ThreadManager;
            //typedef typename SparseMatrixAdapter::block_type Matrix;
            //    static_assert(std::is_same<SparseMatrixAdapter, IstlSparseMatrixAdapter<MatrixBlock>::value,
            //              "The AMGCLSolver linear solver backend requires the IstlSparseMatrixAdapter");
            typedef typename GridView::template Codim<0>::Entity Element;
            typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
            typedef AMGCLSolverBackend<TypeTag>  SuperClass;
        public:
            typedef Dune::AssembledLinearOperator< Matrix, Vector, Vector > AssembledLinearOperatorType;
            
            AMGCLSolverBackendCpr(Simulator& simulator):
                SuperClass(simulator)
            {
            }
            static void registerParameters()
            {
                SuperClass::registerParameters();
            }
            void prepare(const Matrix& M, Vector& b)
            {
                SuperClass::prepare(M, b);
                solver_.reset(new Opm::LinearSolverAmgclCpr(this->prm_));
                solver_->init(this->sz_, this->M_.ptr, this->M_.col, this->M_.val);
            }

            bool solve(Vector& sol)
            {
                std::vector<double> x(this->sz_, 0.0);
                const int nnz = this->M_.ptr[this->sz_];
                double error;
                solver_->solve(this->b_, x, this->iters_, error);
                for (int cell = 0; cell < this->n_; ++cell) {
                    for (int phase = 0; phase < this->np_; ++phase) {
                        sol[cell][phase] = x[this->np_*cell + phase];
                    }
                }
                this->multBlocksVector(sol, this->rightTrans_);
                return true;
            }

            int iterations () const { return this->iters_; }
        private:
            std::unique_ptr<Opm::LinearSolverAmgclCpr> solver_;
                
            
        };


        // BEGIN_PROPERTIES

        // //SET_INT_PROP(AMGCLSolver, LinearSolverVerbosity, 0);
        // SET_TYPE_PROP(AMGCLSolver, LinearSolverBackend,
        //               Ewoms::Linear::AMGCLSolverBackend<TypeTag>);

        // END_PROPERTIES
    }
}
#endif 


