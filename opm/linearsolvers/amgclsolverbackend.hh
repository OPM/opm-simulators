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
#ifndef EWOMS_AMGCLSOLVER_BACKEND_HH
#define EWOMS_AMGCLSOLVER_BACKEND_HH

//#include <ewoms/linear/istlsparsematrixbackend.hh>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/optional.hpp>
#include <iostream>
#include <ewoms/common/parametersystem.hh>
#include <opm/material/common/Unused.hpp>
#include <dune/common/fmatrix.hh>
#include <dune/common/version.hh>
#include <opm/core/linalg/LinearSolverAmgcl.hpp>
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
NEW_TYPE_TAG(AMGCLSolver);
NEW_PROP_TAG(MatrixAddWellContributions);

END_PROPERTIES

namespace Ewoms {
     namespace Linear {
         struct AMGCLMatrixHelper
        {
            AMGCLMatrixHelper()
            {
            }
            AMGCLMatrixHelper(const int sz, const int nnz)
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
         * \brief A linear solver backend for the AMGCL based matrixes uses copy of the BAMGCL matrix using umfpack.
         */
        template <class TypeTag>
        class AMGCLSolverBackend
        {
            typedef typename GET_PROP_TYPE(TypeTag, FluidSystem)       FluidSystem;
            typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
            typedef typename GET_PROP_TYPE(TypeTag, GlobalEqVector) Vector;
            typedef typename GET_PROP_TYPE(TypeTag, SparseMatrixAdapter) SparseMatrixAdapter;
            typedef typename GET_PROP_TYPE(TypeTag, EclWellModel) WellModel;
            typedef typename GET_PROP_TYPE(TypeTag, Simulator) Simulator;
            typedef typename SparseMatrixAdapter::MatrixBlock MatrixBlockType;
            typedef typename SparseMatrixAdapter::IstlMatrix Matrix;
            //typedef typename SparseMatrixAdapter::block_type Matrix;
            //    static_assert(std::is_same<SparseMatrixAdapter, IstlSparseMatrixAdapter<MatrixBlock>::value,
            //              "The AMGCLSolver linear solver backend requires the IstlSparseMatrixAdapter");

        public:
            typedef Dune::AssembledLinearOperator< Matrix, Vector, Vector > AssembledLinearOperatorType;
            AMGCLSolverBackend(Simulator& simulator OPM_UNUSED)
            {
                matrixAddWellContribution_ = EWOMS_GET_PARAM(TypeTag, bool, MatrixAddWellContributions);
                std::ifstream file("amgcl_setup.json");
                boost::property_tree::json_parser::read_json(file, prm_);
                prm_.put("solver_type","cpr_drs");
                np_ = FluidSystem::numPhases;
                prm_.put("block_size", np_);

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
                    std::cout << "AMGCLSolverBackend:: only partyally solving the system" << std::endl;
                    //OPM_THROW(std::runtime_error,
                    //          "Cannot run with amgcl or UMFPACK without also using 'matrix_add_well_contributions=true'.");
                }
                MatrixBlockType leftTrans(0.0);
                bool use_drs = (prm_.get<std::string>("solver_type") == "amgcl_drs");
                if(use_drs){
                    leftTrans=getBlockTransform(2);
                }else{
                    auto eqChange=getBlockTransform(2);
                    auto cprTrans = getBlockTransform(1);
                    leftTrans= cprTrans.rightmultiply(eqChange);
                }
                auto M_cp = M;// = ebosSimulator_.model().linearizer().jacobian().istlMatrix();
                auto b_cp = b;
                multBlocksInMatrix(M_cp, leftTrans, true);
                multBlocksVector(b_cp, leftTrans);
                
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
                
                n_ = M.N();
                sz_ = np_*n_;
                //bool do_transpose = false;
                M_ = buildAMGCLMatrixNoBlocks(M_cp);//, do_transpose);
                b_.resize(sz_, 0.0);
                
                for (int cell = 0; cell < n_; ++cell) {
                    for (int phase = 0; phase < np_; ++phase) {
                        b_[np_*cell + phase] = b_cp[cell][phase];
                    }
                }
            }

            bool solve(Vector& sol)
            {
                std::vector<double> x(sz_, 0.0);
                const int nnz = M_.ptr[sz_];
                Opm::LinearSolverAmgcl solver(prm_);
                double error;
                //solver.solve(sz_, M_.ptr.data(), M_.col.data(), M_.val.data(), b_.data(), x.data(), iters_, error);
                solver.solve(sz_, M_.ptr, M_.col, M_.val, b_, x, iters_, error);
                for (int cell = 0; cell < n_; ++cell) {
                    for (int phase = 0; phase < np_; ++phase) {
                         sol[cell][phase] = x[np_*cell + phase];
                    }
                }
                return true;
             }

            int iterations () const { return iters_; }
        private:
            static MatrixBlockType getBlockTransform(int meth_trans){
                int np = FluidSystem::numPhases;
                MatrixBlockType leftTrans(0.0);
                switch(meth_trans)
                    {
                    case 1 :
                        //cpr
                        for (int row = 0; row < np; ++row) {
                            for (int col = 0; col < np; ++col) {
                                if(row==0){
                                    leftTrans[row][col]=1.0;
                                }else{
                                    if(row==col){
                                        leftTrans[row][col]=1.0;
                                    }
                                }
                            }
                        }
                        break;
                        //permute equations
                    case 2 :
                        for (int row = 0; row < 2; ++row) {
                            for (int col = 0; col < 2; ++col) {
                                if(row!=col){
                                    leftTrans[row][col]=1.0;
                                }
                            }
                        }
                        if(np==3){
                            leftTrans[2][2]=1.0;
                        }
                        break;
                    default:
                        OPM_THROW(std::logic_error,"return zero tranformation matrix");
                    }
                return leftTrans;
            }
            
            static void multBlocksInMatrix(Matrix& ebosJac,const MatrixBlockType& trans,bool left=true){
                const int n = ebosJac.N();
                const int np = FluidSystem::numPhases;
                for (int row_index = 0; row_index < n; ++row_index) {
                    auto& row = ebosJac[row_index];
                    auto* dataptr = row.getptr();
                    //auto* indexptr = row.getindexptr();
                    for (int elem = 0; elem < row.N(); ++elem) {
                        auto& block = dataptr[elem];
                        if(left){
                            block = block.leftmultiply(trans);
                        }else{
                            block = block.rightmultiply(trans);
                        }
                    }
                }
            }
            static void multBlocksVector(Vector& ebosResid_cp,const MatrixBlockType& leftTrans){
                for( auto& bvec: ebosResid_cp){
                    auto bvec_new=bvec;
                    leftTrans.mv(bvec, bvec_new);
                    bvec=bvec_new;
                }
            }
            static AMGCLMatrixHelper buildAMGCLMatrixNoBlocks(const Matrix& ebosJac) 
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
                AMGCLMatrixHelper A(sz, nnz);
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
            AMGCLMatrixHelper M_;
            std::vector<double> b_;
            int sz_;
            int np_;
            int n_;
            int iters_;
            boost::property_tree::ptree prm_;
            std::string solver_type_;
            
        };


        // BEGIN_PROPERTIES

        // //SET_INT_PROP(AMGCLSolver, LinearSolverVerbosity, 0);
        // SET_TYPE_PROP(AMGCLSolver, LinearSolverBackend,
        //               Ewoms::Linear::AMGCLSolverBackend<TypeTag>);

        // END_PROPERTIES
     }
}
#endif 


