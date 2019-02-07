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
#include <ewoms/linear/matrixmarket_ewoms.hh>
BEGIN_PROPERTIES

// forward declaration of the required property tags
NEW_PROP_TAG(Scalar);
NEW_PROP_TAG(Evaluation);
NEW_PROP_TAG(NumEq);
NEW_PROP_TAG(SparseMatrixAdapter);
NEW_PROP_TAG(GlobalEqVector);
NEW_PROP_TAG(EclWellModel);
NEW_PROP_TAG(Simulator);
NEW_PROP_TAG(LinearSolverVerbosity);
NEW_PROP_TAG(AmgclSolverStrategy);
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
            
        public:
            typedef Dune::AssembledLinearOperator< Matrix, Vector, Vector > AssembledLinearOperatorType;
            AMGCLSolverBackend(Simulator& simulator):
                simulator_(simulator),
                solver_type_("amgcl_drs"),
                np_(FluidSystem::numPhases)
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
                
                solver_type_ = EWOMS_GET_PARAM(TypeTag, std::string, AmgclSolverStrategy);
                prm_.put("solver_type",solver_type_);
                np_ = FluidSystem::numPhases;
                prm_.put("block_size", np_);
                
                double tolerance = EWOMS_GET_PARAM(TypeTag, double, LinearSolverReduction);
                int maxiter = EWOMS_GET_PARAM(TypeTag, int, LinearSolverMaxIter);
                int verb_int = EWOMS_GET_PARAM(TypeTag, int, LinearSolverVerbosity);
                printMatrixSystem_ = (verb_int>5);
                if(verb_int>0){
                    prm_.put("verbose", true);
                    std::cout << "Solve system using :" << solver_type_ << std::endl;
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
                EWOMS_REGISTER_PARAM(TypeTag, std::string, AmgclSolverStrategy,
                                     "Solver Strategy for amgcl");
                EWOMS_REGISTER_PARAM(TypeTag, std::string, AmgclSetupFile,
                                     "Setup jeson file for amgcl");
                EWOMS_REGISTER_PARAM(TypeTag, double, LinearSolverReduction, "The minimum reduction of the residual which the linear solver must achieve");
                EWOMS_REGISTER_PARAM(TypeTag, int, LinearSolverMaxIter, "The maximum number of iterations of the linear solver");
                
            }

            /*!
             * \brief Causes the solve() method to discared the structure of the linear system of
             *        equations the next time it is called.
             *
             * Since the SuperLU backend does not create any internal matrices, this is a no-op.linearsolvers
             */
            void eraseMatrix()
            { }

            void prepare(const Matrix& M, Vector& b)
            {
                n_ = M.N();
                sz_ = np_*n_;// assumes block size and number of phases is equal
                if (not(matrixAddWellContribution_)) {
                    std::cout << "AMGCLSolverBackend:: only partyally solving the system" << std::endl;
                    //OPM_THROW(std::runtime_error,
                    //          "Cannot run with amgcl or UMFPACK without also using 'matrix_add_well_contributions=true'.");
                }
                auto M_cp = M;// = ebosSimulator_.model().linearizer().jacobian().istlMatrix();
                auto b_cp = b;
                
                double pressure_scale_ = 50e5;
                rightTrans_ =getPressureTransform(pressure_scale_);
                multBlocksInMatrix(M_cp, rightTrans_, false);
                MatrixBlockType leftTrans(0.0);
                if( solver_type_ == "amgcl_quasiimpes"){
                    //leftTrans=getBlockTransform(3);
                    leftTrans=getBlockTransform(4);
                    scaleCPRSystem(M_cp, b_cp, leftTrans);
                    prm_.put<bool>("use_drs",true);                    
                }
                else if ( solver_type_ == "working_hack_drs"){
                    // leftTrans=getBlockTransform(1);
                    // MatrixBlockType eqChange=getBlockTransform(2);
                    // MatrixBlockType scaleEq =getBlockTransform(3);
                    // leftTrans = leftTrans.rightmultiply(eqChange);
                    // leftTrans = leftTrans.leftmultiply(scaleEq);
                    leftTrans=getBlockTransform(4);
                    scaleCPRSystem(M_cp, b_cp, leftTrans);
                    prm_.put<bool>("use_drs",true);
                }
                else if ( solver_type_ == "working_hack"){
                    // leftTrans=getBlockTransform(1);
                    // MatrixBlockType eqChange=getBlockTransform(2);
                    // MatrixBlockType scaleEq =getBlockTransform(3);
                    // leftTrans = leftTrans.rightmultiply(eqChange);
                    // leftTrans = leftTrans.leftmultiply(scaleEq);
                    leftTrans=getBlockTransform(4);
                    scaleCPRSystem(M_cp, b_cp, leftTrans);
                    prm_.put<bool>("use_drs",false);
                }                
                else if (solver_type_ == "new_amgcl_quasiimpes"){
                    //leftTrans=getBlockTransform(3);
                    leftTrans=getBlockTransform(4);
                    prm_.put<bool>("use_drs",true);
                    scaleCPRSystem(M_cp, b_cp, leftTrans);                    
                    // set up quasi impes weights
                    // avoid local scope for vertor
                    weights_ = getQuasiImpesWeights(M_cp);
                    prm_.put("precond.weights", weights_.data());
                    prm_.put("precond.weights_size", weights_.size());
                    prm_.put("precond.eps_dd", -1e8);
                    prm_.put("precond.eps_ps", -1e8);
                }
                else if (solver_type_ == "amgcl_trueimpes") {
                    //leftTrans=getBlockTransform(2);
                    prm_.put<bool>("use_drs",true);
                    weights_ = getStorageWeights();
                    leftTrans=getBlockTransform(4);
                    scaleCPRSystem(M_cp, b_cp, leftTrans);
                    auto leftTinv = leftTrans.transpose();
                    leftTinv.invert();
                    auto bweights = toVector(weights_);
                    multBlocksVector(bweights, leftTinv);
                    weights_ = toStdVector(bweights);    
                    prm_.put("precond.weights", weights_.data());
                    prm_.put("precond.weights_size", weights_.size());
                    prm_.put("precond.eps_dd", -1e8);
                    prm_.put("precond.eps_ps", -1e8);
                }
                //else if (solver_type_ == "amgcl_trueimpes_pressure"){
                //     prm_.put<bool>("use_drs",true);                    
                // trueImpesBlocksInMatrix(M_cp,b_cp);
                // multBlocksInMatrix(M_cp, rightTrans_, false);
                //     weights_ = std::vector<double>(sz_,0.0);
                //     std::vector<double>& fak_weights = weights_;
                //     for(int i=0; i < n_; ++i){
                //         fak_weights[i*np_] = 1;
                //     }
                //     prm_.put("precond.weights", fak_weights.data());
                //     prm_.put("precond.weights_size", fak_weights.size());
                //     prm_.put("precond.eps_dd", -1e8);
                //     prm_.put("precond.eps_ps", -1e8); 
                // }
                // else if (solver_type_ == "amgcl_quasiimpes_pressure"){
                //     // form quasi impes pressure equation
                //     // fake drs to get correct behavoir
                //     prm_.put<bool>("use_drs",true);
                //     quasiImpesBlocksInMatrix(M_cp, b_cp);
                //     //std::vector<double> fak_weights(sz_,0.0);
                //     weights_ = std::vector<double>(sz_,0.0);
                //     std::vector<double>& fak_weights = weights_;
                //     for(int i=0; i < n_; ++i){
                //         fak_weights[i*np_] = 1;
                //     }
                //     prm_.put("precond.weights", fak_weights.data());
                //     prm_.put("precond.weights_size", fak_weights.size());
                //     prm_.put("precond.eps_dd", -1e8);
                //     prm_.put("precond.eps_ps", -1e8);   
                // }
                else{
                    std::cout << "Solver type set to :" << solver_type_ << std::endl;
                    std::cout << "Use amgcl on original system quasi impes" << std::endl;
                    prm_.put("use_drs",false);
                }
                                       
                if(printMatrixSystem_){
                    std::cout << "Matrix rightTrans" << std::endl;
                    Dune::writeMatrixMarket(rightTrans_, std::cout);
                    std::cout << "Matrix leftTrans" << std::endl;
                    Dune::writeMatrixMarket(leftTrans, std::cout);
                    {
                        std::ofstream file("matrix.txt");
                        Dune::writeMatrixMarket(M, file);
                    }
                    {
                        std::ofstream file("rhs.txt");
                        Dune::writeMatrixMarket(b, file);
                    }
                    {
                        std::ofstream file("weights.txt");
                        for(const auto& w: weights_){
                            file << w << std::endl;
                        }
                    }
                    {
                        const auto& sol = simulator_.model().solution(0);
                        std::ofstream file("solution.txt");
                        // Dune::writeMatrixMarket(sol, file);
                        for(int i=0; i < n_; ++i){
                            const auto& priVars = sol[i];
                            for (unsigned eqIdx = 0; eqIdx < np_; ++eqIdx){
                                file << priVars[eqIdx] << " ";
                            }
                            file << priVars.primaryVarsMeaning() << std::endl;
                        }
                    }
                }                
                
                // set active rows
                prm_.put("precond.block_size", np_);
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
                solver.solve(sz_, M_.ptr, M_.col, M_.val, b_, x, iters_, error);
                for (int cell = 0; cell < n_; ++cell) {
                    for (int phase = 0; phase < np_; ++phase) {
                        sol[cell][phase] = x[np_*cell + phase];
                    }
                }
                multBlocksVector(sol, rightTrans_);
                return true;
            }

            int iterations () const { return iters_; }
        protected:
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
                    case 1 : {
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
                    } break;
                        //permute equations
                    case 2 : {
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
                    } break;
                    case 3 :{
                        //cpr
                        for (int row = 0; row < np; ++row) {
                            for (int col = 0; col < np; ++col) {
                                if(row==col){
                                    if(row==2){
                                        leftTrans[row][col]=1.0/100.0;
                                    }else{
                                        leftTrans[row][col]=1.0;
                                    }
                                }
                            }
                        }
                    } break;
                    case 4 :{ // hack which seems to avoid pivot
                        leftTrans=getBlockTransform(1);                                                                            
                        MatrixBlockType eqChange=getBlockTransform(2);                                                             
                        MatrixBlockType scaleEq =getBlockTransform(3);                                                             
                        //leftTrans = leftTrans.rightmultiply(eqChange);// correct if water is first equation                                                             
                        leftTrans = leftTrans.leftmultiply(scaleEq);
                    } break;
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

            Vector toVector(const std::vector<double>& weights){
                Vector bweights(sz_);
                int ind = 0;
                for(int i=0; i < n_; ++i){
                    for(int j=0; j < np_; ++j){
                        bweights[i][j] = weights[ind];
                        ++ind;
                    }                   
                }
                return bweights;
            }
            
            std::vector<double> toStdVector(const Vector& bweights){
                std::vector<double> weights(sz_,0.0);
                int ind = 0;
                for(int i=0; i < n_; ++i){
                    for(int j=0; j < np_; ++j){
                        weights[ind] = bweights[i][j]; 
                        ++ind;
                    }                   
                }
                return weights;
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
           
            std::vector<double> getQuasiImpesWeights(const Matrix &M){               
                std::vector<double> weights(sz_);
                BlockVector rhs(0.0);
                rhs[0] = 1;
                int index = 0;
                const auto endi = M.end();
                for (auto i=M.begin(); i!=endi; ++i){
                    const auto endj = (*i).end();
                    MatrixBlockType diag_block;
                    for (auto j=(*i).begin(); j!=endj; ++j){
                        if(i.index() == j.index()){
                            diag_block = (*j);
                            break;
                        }                        
                    }
                    //auto bweights = rhs;//diag_block.invert()*rhs;
                    BlockVector bweights;
                    auto diag_block_transpose = diag_block.transpose();
                    diag_block_transpose.solve(bweights, rhs);
                    for(int bind=0; bind < np_;++bind){
                        weights[index] = bweights[bind];
                        ++index;
                    }
                }
                return weights;
            }

            std::vector<double> getSimpleWeights(){               
                std::vector<double> weights(sz_);
                for(int i=0; i< n_; ++i){
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
            
            std::vector<double> getStorageWeights(){
                BlockVector rhs(0.0);
                rhs[0] = 1.0;
                //auto model& simulator_.model();
                std::vector<double> weights(sz_);
                int index = 0;
                ElementContext elemCtx(simulator_);
                const auto& vanguard = simulator_.vanguard();
                auto elemIt = vanguard.gridView().template begin</*codim=*/0>();
                const auto& elemEndIt = vanguard.gridView().template end</*codim=*/0>();
                for (; elemIt != elemEndIt; ++elemIt) {
                    const Element& elem = *elemIt;
                    elemCtx.updatePrimaryStencil(elem);
                    elemCtx.updatePrimaryIntensiveQuantities(/*timeIdx=*/0);
                    Dune::FieldVector<Evaluation, FluidSystem::numPhases> storage;
                    //for(int cell=0; cell < n_;++cell){
                    //EqVector& cachedStorage(cell, /*timeIdx=*/0);
                    unsigned threadId = ThreadManager::threadId();
                    simulator_.model().localLinearizer(threadId).localResidual().computeStorage(storage,elemCtx,/*spaceIdx=*/0, /*timeIdx=*/0);
                    Scalar extrusionFactor =
                        elemCtx.intensiveQuantities(0, /*timeIdx=*/0).extrusionFactor();
                    Scalar scvVolume =
                        elemCtx.stencil(/*timeIdx=*/0).subControlVolume(0).volume() * extrusionFactor;
                    Scalar storage_scale = scvVolume / elemCtx.simulator().timeStepSize();
                    MatrixBlockType block;
                    int offset = 0;
                    for(int ii=0; ii< np_; ++ii){
                        for(int jj=0; jj< np_; ++jj){
                            //const auto& vec = storage[ii].derivative(jj);
                            block[ii][jj] = storage[ii].derivative(jj)/storage_scale;
                            if(jj==0){
                                block[ii][jj] *=pressure_scale_;
                            }
                        }
                    }
                    //auto bweights = block.invert()*rhs;
                    //auto bweights = rhs;//diag_block.invert()*rhs;
                    BlockVector bweights;
                    MatrixBlockType block_transpose = block.transpose();
                    block_transpose.solve(bweights, rhs);
                    for(int bind=0; bind < np_;++bind){
                        weights[index] =bweights[bind];
                        ++index;
                    }
                }
                return weights;
            }
            std::vector<double> quasiImpesWeights(Matrix& A){
                std::vector<double> weights(sz_);
                BlockVector rhs(0.0);
                rhs[0] = 1;
                const auto endi = A.end();
                int index = 0;
                for (auto i=A.begin(); i!=endi; ++i){
                    const auto endj = (*i).end();                                        
                    MatrixBlockType diag_block(0.0);
                    for (auto j=(*i).begin(); j!=endj; ++j){
                        if(i.index() == j.index()){
                            diag_block = (*j);
                            break;
                        }                        
                    }                    
                    //bweights = diag_block.invert()*rhs;
                    //auto bweights = rhs;//diag_block.invert()*rhs;
                    BlockVector bweights;
                    auto diag_block_transpose = diag_block.transpose();
                    diag_block_transpose.solve(bweights, rhs);
                    for(int bind=0; bind < np_;++bind){
                        weights[index] =bweights[bind];
                        ++index;
                    }
                }
                return weights;
            }
            
            void quasiImpesBlocksInMatrix(Matrix& A,Vector& b){
                BlockVector rhs(0.0);
                rhs[0] = 1;
                const auto endi = A.end();
                for (auto i=A.begin(); i!=endi; ++i){
                    const auto endj = (*i).end();
                    MatrixBlockType diag_block;
                    for (auto j=(*i).begin(); j!=endj; ++j){
                        if(i.index() == j.index()){
                            diag_block = (*j);
                            break;
                        }                        
                    }
                    //weights = diag_block.invert()*rhs;
                    //auto weights = rhs;//diag_block.invert()*rhs;
                    BlockVector weights;
                    auto diag_block_transpose = diag_block.transpose();
                    diag_block_transpose.solve(weights, rhs);
                    
                    for (auto j=(*i).begin(); j!=endj; ++j){
                        auto & block = (*j);
                        for(int jj=0;jj < np_;++jj){
                            for(int ii=0;ii < np_;++ii){
                                if(jj==0){
                                    block[0][ii] *=weights[jj];
                                }else{
                                    block[0][ii] += block[jj][ii]*weights[jj];
                                }
                            }                           
                        }
                    }
                    auto & bb = b[i.index()];
                    for(int jj=0;jj < np_;++jj){
                        if(jj==0){
                            bb[jj] *= weights[jj];
                        }
                        else{
                            bb[00] += bb[jj]*weights[jj];  
                        }
                    }
                }
            }

            // void trueImpesBlocksInMatrix(Matrix& M,Vector& b){
            //     BlockVector rhs(0.0);
            //     rhs[0] = 1.0;                    
            //     //auto model& simulator_.model();
            //     std::vector<double> weights(sz_);
            //     int index = 0;
            //     ElementContext elemCtx(simulator_);
            //     const auto& vanguard = simulator_.vanguard();
            //     auto elemIt = vanguard.gridView().template begin</*codim=*/0>();
            //     const auto& elemEndIt = vanguard.gridView().template end</*codim=*/0>();
            //     int row_index = 0;
            //     for (; elemIt != elemEndIt; ++elemIt) {
            //         const Element& elem = *elemIt;
            //         elemCtx.updatePrimaryStencil(elem);
            //         elemCtx.updatePrimaryIntensiveQuantities(/*timeIdx=*/0);
            //         Dune::FieldVector<Evaluation, FluidSystem::NumPhases>& storage;
            //         simulator_.model().computeStorage(storage,elemCtx,/*spaceIdx=*/0, /*timeIdx=*/0);
            //         MatrixBlockType block;
            //         int offset = 0;
            //         for(int ii=0; ii< np_; ++ii){
            //             for(int jj=0; jj< np_; ++jj){
            //                 const auto& vec = storage[ii].derivative();
            //                 block[ii][jj] = vec[jj];
            //             }
            //         }
            //         //auto bweights = block.invert()*rhs;
            //         //auto bweights = rhs;//diag_block.invert()*rhs;
            //         BlockVector bweights;
            //         block.solve(bweights, rhs);
            //         // assume iteration and row indexes correspond
            //         auto& row = M[row_index];
            //         auto& bv = b[row_index];
            //         auto* dataptr = row.getptr();
            //         for (int elem = 0; elem < row.N(); ++elem) {
            //             auto& block = dataptr[elem];
            //             for(int j=0; j< np_; ++j){
            //                 for(int i=1; i< np_; ++i){
            //                     block[0][j] += block[i][j]; 
            //                 }
            //             }                        
            //         }
            //         for(int i = 1 ; i < np_ ; ++i){
            //             bv[0] += bv[i]; 
            //         }
            //         ++row_index;
            //     }
            // }
        
            
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
        
            Simulator& simulator_;
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
            std::vector<double> weights_;
            MatrixBlockType rightTrans_;
            
        };


        // BEGIN_PROPERTIES

        // //SET_INT_PROP(AMGCLSolver, LinearSolverVerbosity, 0);
        // SET_TYPE_PROP(AMGCLSolver, LinearSolverBackend,
        //               Ewoms::Linear::AMGCLSolverBackend<TypeTag>);

        // END_PROPERTIES
    }
}
#endif 


