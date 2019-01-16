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
NEW_PROP_TAG(AmgclDrs);
NEW_PROP_TAG(AmgclSetupFile);
NEW_PROP_TAG(LinearSolverBackend);
NEW_TYPE_TAG(AMGCLSolver);
NEW_PROP_TAG(MatrixAddWellContributions);
NEW_PROP_TAG(LinearSolverReduction);
NEW_PROP_TAG(LinearSolverMaxIter);
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
                std::string fileName = EWOMS_GET_PARAM(TypeTag, std::string, AmgclSetupFile);
                //std::string file_name("amgcl_setup.json");
                std::ifstream file(fileName);
                if (file.is_open()) {
                    boost::property_tree::json_parser::read_json(file, prm_);
                }else {
                    // show message:
                    std::cout << "Error opening file " <<  fileName <<std::endl;
                    OPM_THROW(std::runtime_error,"Error opening file");
                }
                
                bool use_drs = EWOMS_GET_PARAM(TypeTag, bool, AmgclDrs);
                if(use_drs){
                    prm_.put("solver_type","cpr_drs");
                    
                }else{
                   prm_.put("solver_type","cpr");
                }
                np_ = FluidSystem::numPhases;
                prm_.put("block_size", np_);
                //prm_.put("precond.block_size", np_);
                //prm_.put("precond.active_rows", 0);
                double tolerance = EWOMS_GET_PARAM(TypeTag, double, LinearSolverReduction);
                int maxiter = EWOMS_GET_PARAM(TypeTag, int, LinearSolverMaxIter);
                int verb_int = EWOMS_GET_PARAM(TypeTag, int, LinearSolverVerbosity);
                if(verb_int>0){
                  prm_.put("verbose", true);
                }else{
                  prm_.put("verbose", false);
                }
                prm_.put("solver.tol", tolerance);
                prm_.put("solver.maxiter", maxiter);

            }

            static void registerParameters()
            {
                EWOMS_REGISTER_PARAM(TypeTag, bool, MatrixAddWellContributions, "Explicitly specify the influences of wells between cells in the Jacobian and preconditioner matrices");
                //FlowLinearSolverParameters::registerParameters<TypeTag>();
                EWOMS_REGISTER_PARAM(TypeTag, int, LinearSolverVerbosity,
                                     "The verbosity level of the linear solver");
                EWOMS_REGISTER_PARAM(TypeTag, bool, AmgclDrs,
                                     "Use amgcl_drs solver");
                EWOMS_REGISTER_PARAM(TypeTag, std::string, AmgclSetupFile,
                                     "Setup jeson file for amgcl");
                EWOMS_REGISTER_PARAM(TypeTag, double, LinearSolverReduction, "The minimum reduction of the residual which the linear solver must achieve");
                EWOMS_REGISTER_PARAM(TypeTag, int, LinearSolverMaxIter, "The maximum number of iterations of the linear solver");
                
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
                pressure_scale_ = 1;
                MatrixBlockType rightTrans=getPressureTransform(pressure_scale_);                
                if(use_drs){
                    leftTrans=getBlockTransform(2);
                    MatrixBlockType scale_eq =getBlockTransform(3);
                    //leftTrans = leftTrans.leftmultiply(scale_eq);
                }else{
                     auto eqChange=getBlockTransform(2);
                     auto cprTrans = getBlockTransform(1);
                     leftTrans= cprTrans.rightmultiply(eqChange);
                }
                auto M_cp = M;// = ebosSimulator_.model().linearizer().jacobian().istlMatrix();
                auto b_cp = b;
                // right transforme to scale pressure
                multBlocksInMatrix(M_cp, rightTrans, false);
                // make system by manipulating equations
                scaleCPRSystem(M_cp, b_cp, leftTrans);
                //multBlocksInMatrix(M_cp, leftTrans, true);
                //multBlocksVector(b_cp, leftTrans);
                
                
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
                // if(use_drs){
                //     if(prm.get<std::string>("strategy") == "quasimpes"){
                //         std::vector<double> drs_weights = getQuasiImpesWeights(M_cp);
                //         prm.put("precond.weights", drs_weights.data());
                //         prm.put("precond.weights_size", drs_weights.size());
                //     }
                // }
                n_ = M.N();
                sz_ = np_*n_;
                // set active rows
                prm_.put("precond.active_rows", sz_);
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
                        if(phase==0){
                            // have transformed pressure variable
                            sol[cell][phase] = x[np_*cell + phase]*pressure_scale_;
                        }else{
                            sol[cell][phase] = x[np_*cell + phase];
                        }
                    }
                }
                return true;
             }

            int iterations () const { return iters_; }
        private:
            static MatrixBlockType getPressureTransform(double pressure_scale){
                int np = FluidSystem::numPhases;
                MatrixBlockType leftTrans(0.0);
                for (int row = 0; row < np; ++row) {
                    for (int col = 0; col < np; ++col) {
                        if(row==col){
                            if(row==0){
                                leftTrans[row][col]=pressure_scale;                                
                            }else{
                                leftTrans[row][col]=1.0;
                            }
                        }
                    }
                }
                return leftTrans;
            }
                
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
                    case 3 :
                        //cpr
                        for (int row = 0; row < np; ++row) {
                            for (int col = 0; col < np; ++col) {
                                if(row==col){
                                    if(row==3){
                                        leftTrans[row][col]=1.0/100;
                                    }else{
                                        leftTrans[row][col]=1.0;
                                    }
                                }
                            }
                        }    
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
            static void scaleCPRSystem(Matrix& M_cp,Vector& b_cp,const MatrixBlockType& leftTrans){
                multBlocksInMatrix(M_cp, leftTrans, true);
                multBlocksVector(b_cp, leftTrans);
            }
            /*
            static std::vector<double> getQuasiImpesWeights(const Matrix &M){               
                std::vector<double> weights(np_*sz_);
                BlockVector rhs(0.0);
                rhs[0] = 1;
                int index = 0;
                const auto endi = A.end();
                for (const auto i=A.begin(); i!=endi; ++i){
                    const auto endj = (*i).end();
                    MatrixBlockType diag_block;
                    for (const auto j=(*i).begin(); j!=endj; ++j){
                        if(i.index() == j.index()){
                            diag_block = (*j);
                            break;
                        }                        
                    }
                    auto bweights = diag_block.inv()*rhs;
                    for(int bind=0; bind < np_;++bind){
                        weights[index] =bweights[bind];
                        ++index;
                    }
                }
                return weights;
            }

            static std::vector<double> getSimpleWeights(){               
                std::vector<double> weights(np_*sz_);
                for(int i=0; i< sz_; ++i){
                    for(int j=0; j< np_; ++j){
                        if(j==3){
                            weights[i*np_+j] = 1;
                        }else{
                            weights[i*np_+j] = 1/100;
                        }
                    }
                }
                return weights;
            }
            
            static std::vector<double> getStorageWeights(){
                BlockVector rhs(0.0);
                rhs[0] = 1.0;
                    
                auto model& simulator_.model();
                std::vector<double> weights(np_*sz_);
                int index = 0;
                // iterate the element context ie. cells
                {    
                    Dune::FieldVector<LhsEval, numEq>& storage;
                    model.computeStorage(storage,elemCtx, dofIdx, timeIdx);
                    BlockMatrix block;
                    int offset = 0;
                    for(int ii=0; ii< np_; ++ii){
                        for(int jj=0; jj< np_; ++jj){
                            const auto& vec = storage[ii].derivative();
                            block[ii][jj] = vec[jj];
                        }
                    }
                    auto bweights = diag_block.inv()*rhs;
                    for(int bind=0; bind < np_;++bind){
                        weights[index] =bweights[bind];
                        ++index;
                    }
                }
                return weights;
     
                
                
            }
            
            static void multBlocksInMatrix(Matrix& A,Vector& b){
                BlockVector rhs(0.0);
                rhs[0] = 1;
                const auto endi = A.end();
                for (const auto i=A.begin(); i!=endi; ++i){
                    const auto endj = (*i).end();
                    
                    BlockVector weights;
                    switch (strategy){
                    case "quasiimpes" :
                        MatrixBlockType diag_block;
                        for (const auto j=(*i).begin(); j!=endj; ++j){
                            if(i.index() == j.index()){
                                diag_block = (*j);
                                break;
                            }                        
                        }
                        weights = diag_block.inv()*rhs;
                    case "impes" :
                        
                        
                    default :
                        for(int ii=0; ii< np_; ++ii){
                            if(ii==3){//gas
                                weights[ii] = 1/100;
                            }else{
                                weights[ii] = 1;
                            }
                            
                        }                        
                    }
                        
                    for (auto j=(*i).begin(); j!=endj; ++j){
                        auto block& = (*j);
                        for(int jj=0;jj < np_;++jj){
                            for(int ii=0;ii < np_;++ii){
                                if(jj==0){
                                    block[0][ii] *=weights[jj];
                                    b[0] *=weights[jj];
                                }else{
                                    block[0][ii] += block[jj][ii]*weights[jj];
                                    b[0] += b[jj]*weights[jj];
                                }
                            }
                        }
                    }
                    
                    }
                }

            */
            
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
            double pressure_scale_;
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


