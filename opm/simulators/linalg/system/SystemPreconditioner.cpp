#include <config.h>
#include <opm/simulators/linalg/system/SystemPreconditioner.hpp>
#include <opm/simulators/linalg/FlexibleSolver_impl.hpp>



namespace Opm {
    SystemPreconditioner::SystemPreconditioner(const SystemMatrix& S,const std::function<ResVector()> &weightsCalculator,int pressureIndex, 
    const Opm::PropertyTree& prm) :
        S_(S), prm_(prm) {
            auto rop = std::make_unique<ResOperator>(S[_0][_0]);
            auto wop = std::make_unique<WellOperator>(S[_1][_1]);
            auto resprm = prm_.get_child("reservoir_solver");
            auto resprmsmoother = prm_.get_child("reservoir_smoother");
            auto wellprm = prm_.get_child("well_solver");
            //std::function<ResVector()> weightsCalculatorRes;
            auto rsol = std::make_unique<ResFlexibleSolverType>(*rop, resprm,
                                                             weightsCalculator,
                                                             pressureIndex);

            auto rsmoother = std::make_unique<ResFlexibleSolverType>(*rop, resprmsmoother,
                                                             weightsCalculator,
                                                             pressureIndex);

            std::function<WellVector()> weightsCalculatorWell;
            auto wsol = std::make_unique<WellFlexibleSolverType>(*wop, wellprm,
                                                             weightsCalculatorWell,
                                                             pressureIndex); 
            this->rop_ = std::move(rop);
            this->wop_ = std::move(wop);
            this->resSolver_ = std::move(rsol);
            this->wellSolver_ = std::move(wsol);
            this->resSmoother_ = std::move(rsmoother);
        }
    
    void 
    SystemPreconditioner::apply(SystemVector& v, const SystemVector& d) {
            // change order?
            
            double well_tol = prm_.get<double>("well_solver.tol");
            double res_tol = prm_.get<double>("reservoir_solver.tol");
            Dune::InverseOperatorResult well_result;
            // references to system
            const auto& A = S_[_0][_0];
            const auto& C = S_[_1][_0];
            const auto& B = S_[_0][_1];
            const auto& D = S_[_1][_1];
            const auto& r_r = d[_0];
            const auto& r_w = d[_1];

            auto wRes = r_w;
            auto resRes = r_r;
            //assert(v[_0].two_norm2() <1e-30);
            //assert(r_w.two_norm2() <1e-30);
            //tmp variables to be updated
            WellVector wSol(wRes.size());
            ResVector resSol(resRes.size());
            ResVector dresSol(resRes.size());
            WellVector dwSol(wSol.size());
            resSol = 0.0;
            wSol = 0.0; 
            //ResVector dresRes(dresRes.size());
            dresSol = 0.0;
            dwSol = 0.0;
            
            //dresRes = 0.0;
            B.mmv(wSol, resRes);// probably does noting as wSol =0
            // Solve reservoir part for "pressure"
            auto tmp_resRes = resRes;// need copy of residual since solver change it
            
            {
                
                Dune::InverseOperatorResult res_result;    
                resSolver_->apply(dresSol, tmp_resRes, res_tol, res_result);
                // update residual
                resSol += dresSol;//update solution
                //update residual
                A.mmv(dresSol, resRes);
                C.mmv(dresSol, wRes); 
            }
            // solve reservoir part again for all secons "second stage of cpr"
            auto tmp_wRes = wRes;
            { // smoother in cpr
                dwSol = 0.0;
                wellSolver_->apply(dwSol, tmp_wRes, well_tol, well_result);
                wSol += dwSol;
                B.mmv(dwSol, resRes);
                D.mmv(dwSol, wRes);
                dresSol = 0.0; 
                tmp_resRes = resRes;
                Dune::InverseOperatorResult res_result;
                //std::cout << "TwoLevelMethodCpr apply input:" << std::endl;
                //Dune::writeMatrixMarket(tmp_resRes,std::cout);

                resSmoother_->apply(dresSol, tmp_resRes, res_tol, res_result);

                resSol += dresSol;//update solution
                //update residual
                A.mmv(dresSol, resRes); 
                C.mmv(dresSol, wRes); 
                //std::cout << "TwoLevelMethodCpr apply result:" << std::endl;
                //Dune::writeMatrixMarket(resSol,std::cout);

             }
            {
                dwSol = 0.0;
                tmp_wRes = wRes;
                wellSolver_->apply(dwSol, tmp_wRes, well_tol, well_result);
                wSol += dwSol;
            }
            v[_0] = resSol;
            v[_1] = wSol;
            if(false){
                std::cout << "WellSol1 output:" << std::endl;
                Dune::writeMatrixMarket(wSol,std::cout);
                WellVector tmp=r_w;//;(wRes.size());
                //wRes = 0.0;
                C.mmv(resSol,tmp);
                auto tmp_tmp = tmp;
                wellSolver_->apply(wSol, tmp_tmp, well_tol, well_result);
                D.mmv(wSol,tmp);
                std::cout << "WellSol2 output:" << std::endl;
                Dune::writeMatrixMarket(wSol,std::cout);
                v[_1] = wSol;
                std::cout << "SystemPreconditioner apply output:" << std::endl;
                Dune::writeMatrixMarket(tmp,std::cout);
                std::cout << "Norms: wRes " << tmp.two_norm() << "rhs" << r_w.two_norm() << std::endl;
                std::cout << "-----------D----------" << std::endl;
                Dune::writeMatrixMarket(D,std::cout);
                std::cout << "-----------C----------" << std::endl;
                Dune::writeMatrixMarket(C,std::cout);
            }

    }
    
    SystemPreconditionerParallel::SystemPreconditionerParallel(const SystemMatrix& S,const std::function<ResVector()> &weightsCalculator,int pressureIndex, 
                                               const Opm::PropertyTree& prm,const SystemComm& syscomm) :
      S_(S), prm_(prm), syscomm_(syscomm) {
            auto rop = std::make_unique<ResOperator>(S[_0][_0],syscomm_[_0]);
            auto wop = std::make_unique<WellOperator>(S[_1][_1]);
            auto resprm = prm_.get_child("reservoir_solver");
            auto wellprm = prm_.get_child("well_solver");
            auto resprmsmoother = prm_.get_child("reservoir_smoother");
            //std::function<ResVector()> weightsCalculatorRes;
            auto rsol = std::make_unique<ResFlexibleSolverType>(*rop,
                                                                syscomm_[_0],   
                                                                resprm,
                                                                weightsCalculator,
                                                                pressureIndex);

            auto rsmoother = std::make_unique<ResFlexibleSolverType>(*rop,
                                                                syscomm_[_0],   
                                                                resprmsmoother,
                                                                weightsCalculator,
                                                                pressureIndex);

            std::function<WellVector()> weightsCalculatorWell;
            auto wsol = std::make_unique<WellFlexibleSolverType>(*wop, wellprm,
                                                             weightsCalculatorWell,
                                                             pressureIndex); 
            this->rop_ = std::move(rop);
            this->wop_ = std::move(wop);
            this->resSolver_ = std::move(rsol);
            this->wellSolver_ = std::move(wsol);
            this->resSmoother_ = std::move(rsmoother);
        }
    
    void 
    SystemPreconditionerParallel::apply(SystemVector& v, const SystemVector& d) {
            // change order?
            double well_tol = prm_.get<double>("well_solver.tol");
            double res_tol = prm_.get<double>("reservoir_solver.tol");
            const auto& A = S_[_0][_0];
            const auto& C = S_[_1][_0];
            const auto& B = S_[_0][_1];
            const auto& D = S_[_1][_1];
            const auto& r_r = d[_0];
            const auto& r_w = d[_1];

            auto wRes = r_w;
            auto resRes = r_r;
            //assert(v[_0].two_norm2() <1e-30);
            //assert(r_w.two_norm2() <1e-30);
            //C.mv(v[_0], wRes); 

            

            //wellSolver_->apply(wSol, wRes, well_tol, well_result); // seems to do things not good
            //v[_1] = wSol;
            // Use reservoir solver
            
            WellVector wSol(wRes.size());
            ResVector resSol(resRes.size());
            ResVector dresSol(resRes.size());
            WellVector dwSol(wSol.size());
            //ResVector dresRes(dresRes.size());
            resSol = 0.0;
            wSol = 0.0;
            dresSol = 0.0;
            dwSol = 0.0;

            //dresRes = 0.0;
            //B.mmv(wSol, resRes);
            
            auto tmp_resRes = resRes;
            {   
                Dune::InverseOperatorResult res_result;
                syscomm_[_0].copyOwnerToAll(tmp_resRes,tmp_resRes);
                resSolver_->apply(dresSol, tmp_resRes, res_tol, res_result);
                // update residual
                resSol += dresSol;//update solution
                //update residual
                A.mmv(dresSol, resRes);
                C.mmv(dresSol, wRes); 
            }
            
             { // smoother in cpr
                dwSol = 0.0;
                auto tmp_wRes = wRes;
                Dune::InverseOperatorResult well_result;
                wellSolver_->apply(dwSol, tmp_wRes, well_tol, well_result);
                wSol += dwSol;
                B.mmv(dwSol, resRes);
                D.mmv(dwSol, wRes);
              
                dresSol = 0.0;
                tmp_resRes = resRes;
                syscomm_[_0].copyOwnerToAll(tmp_resRes,tmp_resRes);
                Dune::InverseOperatorResult res_result;
                //std::cout << "TwoLevelMethodCpr apply rhs:" << std::endl;
                //Dune::writeMatrixMarket(tmp_resRes,std::cout);
                resSmoother_->apply(dresSol, tmp_resRes, res_tol, res_result);
                resSol += dresSol;//update solution
                //update residual
                A.mmv(dresSol, resRes); 
                C.mmv(dresSol, wRes);
                //std::cout << "TwoLevelMethodCpr apply result:" << std::endl;
                //Dune::writeMatrixMarket(dresSol,std::cout);
    
 
             }


            {
                dwSol = 0.0;
                auto tmp_wRes = wRes;
                Dune::InverseOperatorResult well_result;
                wellSolver_->apply(dwSol, tmp_wRes, well_tol, well_result);
                wSol += dwSol;
            }
            syscomm_[_0].copyOwnerToAll(resSol,resSol);
            v[_0] = resSol;
            v[_1] = wSol;
            
            // update well rhs

            // update residual
            // S_[1][_0].mv(vPart, wellPart, -1.0, 1.0);
        
    }

}
const int numResDofs = 3;
const int numWellDofs = 4;
using CommSeq = Dune::Amg::SequentialInformation;
using CommPar = Dune::OwnerOverlapCopyCommunication<int, int>;
using WWMatrix = Dune::BCRSMatrix<Dune::FieldMatrix<double, numWellDofs, numWellDofs>>;
using WellVector = Dune::BlockVector<Dune::FieldVector<double, numWellDofs>>;
using WellOp = Dune::MatrixAdapter<WWMatrix, WellVector, WellVector>;
using WellFlexibleSolverType = Dune::FlexibleSolver<WellOp>;  
template class Dune::FlexibleSolver<WellOp>;
//template class Dune::FlexibleSolver<WellOp>;
/* template Dune::FlexibleSolver<WellOp>::FlexibleSolver(WellOp& op,
                     const CommPar& comm,
                     const Opm::PropertyTree& prm,
                     const std::function<WellVector()>& weightsCalculator,
                     std::size_t pressureIndex); */
//template class Opm::PreconditionerFactory<WellOp, CommSeq>;
//template class Opm::PreconditionerFactory<WellOp, CommPar>;                                          
// Define matrix and vector typ
