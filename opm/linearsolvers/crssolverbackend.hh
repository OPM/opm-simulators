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
#ifndef EWOMS_CRSSOLVER_BACKEND_HH
#define EWOMS_CRSSOLVER_BACKEND_HH

//#include <ewoms/linear/istlsparsematrixbackend.hh>
#include <ewoms/common/parametersystem.hh>
#include <opm/material/common/Unused.hpp>
#include <dune/common/fmatrix.hh>
#include <dune/common/version.hh>
#include <opm/core/linalg/LinearSolverUmfpack.hpp>
//#include <ewoms/linear/matrixmarket_ewoms.hh>
BEGIN_PROPERTIES

// forward declaration of the required property tags
NEW_PROP_TAG(Scalar);
NEW_PROP_TAG(NumEq);
NEW_PROP_TAG(SparseMatrixAdapter);
NEW_PROP_TAG(GlobalEqVector);
NEW_PROP_TAG(EclWellModel);
NEW_PROP_TAG(Simulator);
NEW_PROP_TAG(LinearSolverVerbosity);
NEW_PROP_TAG(LinearSolverBackend);
NEW_TYPE_TAG(CRSSolver);
NEW_PROP_TAG(MatrixAddWellContributions);

END_PROPERTIES

namespace Ewoms {
     namespace Linear {
         struct CRSMatrixHelper
        {
            CRSMatrixHelper()
            {
            }
            CRSMatrixHelper(const int sz, const int nnz)
                : ptr(sz + 1, -1)
                , col(nnz, -1)
                , val(nnz, 0.0)
            {
            }
            void print(){
                {
                    std::ofstream file("ptr.txt");
                    for(auto x: ptr){
                        file << x << std::endl;
                    }
                }
                {
                    std::ofstream file("col.txt");
                    for(auto x: col){
                        file << x << std::endl;
                    }
                }
                {
                    std::ofstream file("val.txt");
                    for(auto x: val){
                        file << x << std::endl;
                    }
                }
            }
            std::vector<int>    ptr;
            std::vector<int>    col;
            std::vector<double> val;
        };

        
         
        /*!
         * \ingroup Linear
         * \brief A linear solver backend for the CRS based matrixes uses copy of the BCRS matrix using umfpack.
         */
        template <class TypeTag>
        class CRSSolverBackend
        {
            typedef typename GET_PROP_TYPE(TypeTag, FluidSystem)       FluidSystem;
            typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
            typedef typename GET_PROP_TYPE(TypeTag, GlobalEqVector) Vector;
            typedef typename GET_PROP_TYPE(TypeTag, SparseMatrixAdapter) SparseMatrixAdapter;
            typedef typename GET_PROP_TYPE(TypeTag, EclWellModel) WellModel;
            typedef typename GET_PROP_TYPE(TypeTag, Simulator) Simulator;
            typedef typename SparseMatrixAdapter::IstlMatrix Matrix;
            //typedef typename SparseMatrixAdapter::block_type Matrix;
            //    static_assert(std::is_same<SparseMatrixAdapter, IstlSparseMatrixAdapter<MatrixBlock>::value,
            //              "The CRSSolver linear solver backend requires the IstlSparseMatrixAdapter");

        public:
            typedef Dune::AssembledLinearOperator< Matrix, Vector, Vector > AssembledLinearOperatorType;
            CRSSolverBackend(Simulator& simulator OPM_UNUSED)
            {
                matrixAddWellContribution_ = EWOMS_GET_PARAM(TypeTag, bool, MatrixAddWellContributions);
            }

            static void registerParameters()
            {
                EWOMS_REGISTER_PARAM(TypeTag, bool, MatrixAddWellContributions, "Explicitly specify the influences of wells between cells in the Jacobian and preconditioner matrices");
                //FlowLinearSolverParameters::registerParameters<TypeTag>();
                //EWOMS_REGISTER_PARAM(TypeTag, int, LinearSolverVerbosity,
                //                     "The verbosity level of the linear solver");
            }

            /*!
             * \brief Causes the solve() method to discared the structure of the linear system of
             *        equations the next time it is called.
             *
             * Since the SuperLU backend does not create any internal matrices, this is a no-op.
             */
            void eraseMatrix()
            { }

            void prepare(const Matrix& M, Vector& b)
            {
                if (not(matrixAddWellContribution_)) {
                    std::cout << "CRSSolverBackend:: only partyally solving the system" << std::endl;
                    //OPM_THROW(std::runtime_error,
                    //          "Cannot run with amgcl or UMFPACK without also using 'matrix_add_well_contributions=true'.");
                }  
                if(printMatrixSystem_){
                    // used for testing matixes outside flow
                    /*
                    std::ofstream filem("matrix.txt");
                    Dune::writeMatrixMarket(M, filem);
                    std::ofstream fileb("rhs.txt");
                    Dune::writeMatrixMarket(b, fileb);
                    */
                }
                
                // should this be done with move??  
                
                np_ = FluidSystem::numPhases;
                n_ = M.N();
                sz_ = np_*n_;
                //bool do_transpose = false;
                M_ = buildCRSMatrixNoBlocks(M);//, do_transpose);
                b_.resize(sz_, 0.0);
                
                for (int cell = 0; cell < n_; ++cell) {
                    for (int phase = 0; phase < np_; ++phase) {
                        b_[np_*cell + phase] = b[cell][phase];
                    }
                }
            }

            bool solve(Vector& sol)
            {
                std::vector<double> x(sz_, 0.0);
                const int nnz = M_.ptr[sz_];
                Opm::LinearSolverUmfpack solver;
                solver.solve(sz_, nnz, M_.ptr.data(), M_.col.data(), M_.val.data(), b_.data(), x.data());
                for (int cell = 0; cell < n_; ++cell) {
                    for (int phase = 0; phase < np_; ++phase) {
                         sol[cell][phase] = x[np_*cell + phase];
                    }
                }
                return true;
                /*CRSMatrix A  = {
                  (size_t) sz,
                  (size_t) nnz,
                  const_cast<int*>(M_.prt.data()),
                  const_cast<int*>(M_.col.data()),
                  const_cast<double*>(M_.val.data())
                  };
                  call_UMFPACK(&A, b_.data(), sol.data());
                */
            }

           int iterations () const { return 0; }
        private:
            static CRSMatrixHelper buildCRSMatrixNoBlocks(const Matrix& ebosJac) 
            {
                /*
                  const auto& ebosJacOrg = ebosSimulator_.model().linearizer().matrix();
                  auto ebosJac = ebosJacOrg;
                  if(do_transpose){
                  Dune::MatrixVector::transpose(ebosJacOrg, ebosJac);
              }
                */
                const int n = ebosJac.N();
                const int np = FluidSystem::numPhases;
                const int sz = n * np;
                const int nnz = ebosJac.nonzeroes() * np * np;
                CRSMatrixHelper A(sz, nnz);
                A.ptr[0] = 0;
                int index = 0;
                for (int row_index = 0; row_index < n; ++row_index) {
                    const auto& row = ebosJac[row_index];
                    const auto* dataptr = row.getptr();
                    const auto* indexptr = row.getindexptr();
                    //for (int brow = np-1; brow > -1; --brow) {
                    for (int brow = 0; brow < np; ++brow) {
                        A.ptr[row_index*np + brow + 1] = A.ptr[row_index*np + brow] + row.N() * np;
                        for (int elem = 0; elem < row.N(); ++elem) {
                            const int istlcol = indexptr[elem];
                            const auto& block = dataptr[elem];
                            for (int bcol = 0; bcol < np; ++bcol) {                            
                                A.val[index] = block[brow][bcol];
                                //A.val[index] = block[np - brow - 1][bcol];
                                A.col[index] = np*istlcol + bcol;
                                ++index;
                            }
                        }
                    }
                }
                return A;
            }

            bool matrixAddWellContribution_;
            bool printMatrixSystem_;
            CRSMatrixHelper M_;
            std::vector<double> b_;
            int sz_;
            int np_;
            int n_;
            
        };


        // BEGIN_PROPERTIES

        // //SET_INT_PROP(CRSSolver, LinearSolverVerbosity, 0);
        // SET_TYPE_PROP(CRSSolver, LinearSolverBackend,
        //               Ewoms::Linear::CRSSolverBackend<TypeTag>);

        // END_PROPERTIES
     }
}
#endif 


