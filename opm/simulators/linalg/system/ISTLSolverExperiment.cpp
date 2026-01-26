#include "config.h"
#include <iostream>
#include <dune/istl/bvector.hh>
#include <opm/simulators/linalg/matrixblock.hh>
#include <dune/common/fmatrix.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/multitypeblockmatrix.hh>
#include <dune/istl/multitypeblockvector.hh>
#include "WellMatrixMerger.hpp"
#include "SystemPreconditioner.hpp"
#include <opm/common/ErrorMacros.hpp>

namespace Opm
{
        namespace SystemSolver {
            const int numResDofs = 3;
            const int numWellDofs = 4;

                // Define matrix and vector types
             //using RRMatrix = Dune::BCRSMatrix<Dune::FieldMatrix<double, numResDofs, numResDofs>>;
             using RRMatrix = Dune::BCRSMatrix<Opm::MatrixBlock<double, numResDofs, numResDofs>>;
                //using RWtype = Dune::FieldMatrix<double, numResDofs, numWellDofs>;
                using RWMatrix = Dune::BCRSMatrix<Dune::FieldMatrix<double, numResDofs, numWellDofs>>;
                using WRMatrix = Dune::BCRSMatrix<Dune::FieldMatrix<double, numWellDofs, numResDofs>>;
                using WWMatrix = Dune::BCRSMatrix<Dune::FieldMatrix<double, numWellDofs, numWellDofs>>;
                using RVector = Dune::BlockVector<Dune::FieldVector<double, numResDofs>>;
                using WVector = Dune::BlockVector<Dune::FieldVector<double, numWellDofs>>;
                  using SystemMatrix = Dune::MultiTypeBlockMatrix<
                Dune::MultiTypeBlockVector<RRMatrix, RWMatrix>,
                Dune::MultiTypeBlockVector<WRMatrix, WWMatrix>
        >;

        using SystemMatrix = Dune::MultiTypeBlockMatrix<
                Dune::MultiTypeBlockVector<RRMatrix, RWMatrix>,
                Dune::MultiTypeBlockVector<WRMatrix, WWMatrix>
        >; 
        using SystemVector = Dune::MultiTypeBlockVector<RVector, WVector>;

    Dune::InverseOperatorResult solveSystem(const SystemMatrix& S, SystemVector& x, const SystemVector& b, 
        const std::function<RVector()> &weightCalculator,int pressureIndex, const Opm::PropertyTree& prm)
    {
        // Here we would implement the solver logic for the system S * x = b
        // This is a placeholder implementation
        int verbosity = prm.get<int>("verbosity");           // Reduce output verbosity
        if(verbosity){
            std::cout << "Solving system with merged matrices..." << std::endl;
        }
         const Dune::MatrixAdapter<SystemMatrix, SystemVector, SystemVector> S_linop(S);
         //const Dune::OverlappingSchwarzOperator<SystemVector, SystemVector, Dune::OwnerOverlapCopyCommunication<int, int> > S_linop(S_linop, comm);
    //TailoredPrecondDiag precond(S,prm);
    Opm::PropertyTree precond_prm = prm.get_child("preconditioner");
    SystemPreconditioner precond(S,weightCalculator, pressureIndex, precond_prm);
    
    // Set solver parameters
    double linsolve_tol = prm.get<double>("tol");  // Less strict tolerance
    int max_iter = prm.get<int>("maxiter");           // Limit iterations
   
    Dune::InverseOperatorResult result;
    // Create and run the solver with error handling
    try {
        if(verbosity > 0){
        std::cout << "Solving system with BiCGSTAB solver..." << std::endl;
        }
        const std::string solver_type = prm.get<std::string>("solver");
        using AbstractSolverType = Dune::InverseOperator<SystemVector, SystemVector>;
        std::shared_ptr<AbstractSolverType> linsolver;
        if( solver_type == "bicgstab"){
            linsolver = std::make_shared<Dune::BiCGSTABSolver<SystemVector>>(
                                                                             S_linop,
                                                                             precond,
                                                                             linsolve_tol,
                                                                             max_iter,
                                                                             verbosity
                                                                             );
        } else if ( solver_type == "fgmres"){
          int restart = prm.get<int>("restart", 15);
          linsolver = std::make_shared<Dune::RestartedGMResSolver<SystemVector>>(
                                                                               S_linop,
                                                                               precond,
                                                                               linsolve_tol,
                                                                               restart,
                                                                               max_iter,
                                                                               verbosity
                                                                             );
        }else {
          OPM_THROW(std::invalid_argument,
                      "Properties: Solver " + solver_type + " not known.");
        }
        auto residual(b);
        linsolver->apply(x, residual, result);
        //assert(false);//debug
        // Print results
        if(verbosity > 10){
        std::cout << "\nSolver results:" << std::endl;
        std::cout << "  Converged: " << (result.converged ? "Yes" : "No") << std::endl;
        std::cout << "  Iterations: " << result.iterations << std::endl;
        std::cout << "  Reduction: " << result.reduction << std::endl;
        std::cout << "  Elapsed time: " << result.elapsed << " seconds" << std::endl;
        std::cout << "\nMatrix merger example completed successfully!" << std::endl;
        }
        //printVector("Solution x_r", x[_0]);
        //printVector("Solution x_w", x[_1]);
    }
    catch (const std::exception& e) {
        std::cerr << "Error solving system: " << e.what() << std::endl;
        //return result;
    }
    
    
    return result;
    }

    Dune::InverseOperatorResult solveSystem(const SystemMatrix& S, SystemVector& x, const SystemVector& b, 
        const std::function<RVector()> &weightCalculator,int pressureIndex, const Opm::PropertyTree& prm,
        const Dune::OwnerOverlapCopyCommunication<int, int>& comm)
    {
        // Here we would implement the solver logic for the system S * x = b
        // This is a placeholder implementation
          const bool is_iorank = comm.communicator().rank() == 0;
        const int verbosity = is_iorank ? prm.get<int>("verbosity", 0) : 0;
        if(verbosity){
            std::cout << "Solving system with merged matrices..." << std::endl;
        }
         //const Dune::MatrixAdapter<SystemMatrix, SystemVector, SystemVector> S_linop(S);
        using WellComm = Dune::JacComm;
        using SystemComm = Dune::MultiCommunicator<const Dune::OwnerOverlapCopyCommunication<int, int>&,const WellComm&>;
         WellComm seqComm;
         SystemComm systemComm(comm,seqComm);
         
         const Dune::OverlappingSchwarzOperator<SystemMatrix, SystemVector, SystemVector, SystemComm > S_linop(S, systemComm);
         std::shared_ptr< Dune::ScalarProduct<SystemVector> > scalarproduct = Dune::createScalarProduct<SystemVector, SystemComm>(systemComm, S_linop.category());
         //
    //TailoredPrecondDiag precond(S,prm);
    Opm::PropertyTree precond_prm = prm.get_child("preconditioner");
    SystemPreconditionerParallel precond(S,weightCalculator, pressureIndex, precond_prm,systemComm);
    Dune::BlockPreconditioner<SystemVector, SystemVector, SystemComm,SystemPreconditionerParallel> sysprecond(precond, systemComm);
    
    // Set solver parameters
    double linsolve_tol = prm.get<double>("tol");  // Less strict tolerance
    int max_iter = prm.get<int>("maxiter");           // Limit iterations
   
    Dune::InverseOperatorResult result;
    // Create and run the solver with error handling
    try {
        if(verbosity > 0){
        std::cout << "Solving system with BiCGSTAB solver parallel...rank.." << comm.communicator().rank() << std::endl;
        }
        using AbstractSolverType = Dune::InverseOperator<SystemVector, SystemVector>;
        std::shared_ptr<AbstractSolverType> linsolver;
        const std::string solver_type = prm.get<std::string>("solver");
        if( solver_type == "bicgstab"){
            linsolver = std::make_shared<Dune::BiCGSTABSolver<SystemVector>>(
                                                                             S_linop,
                                                                             *scalarproduct,
                                                                             sysprecond,
                                                                             linsolve_tol,
                                                                             max_iter,
                                                                             verbosity
                                                                             );
        }else if ( solver_type == "fgmres"){
          int restart = prm.get<int>("restart", 15);
          linsolver = std::make_shared<Dune::RestartedGMResSolver<SystemVector>>(
                                                                             S_linop,
                                                                             *scalarproduct,
                                                                             sysprecond,
                                                                             linsolve_tol,
                                                                             restart,
                                                                             max_iter,
                                                                             verbosity
                                                                             );
        }else {
          OPM_THROW(std::invalid_argument,
                      "Properties: Solver " + solver_type + " not known.");
        }
        auto residual(b);
        linsolver->apply(x, residual, result);
        //assert(false);
                      
        // Print results
        if(verbosity > 10){
        std::cout << "\nSolver results:" << std::endl;
        std::cout << "  Converged: " << (result.converged ? "Yes" : "No") << std::endl;
        std::cout << "  Iterations: " << result.iterations << std::endl;
        std::cout << "  Reduction: " << result.reduction << std::endl;
        std::cout << "  Elapsed time: " << result.elapsed << " seconds" << std::endl;
        std::cout << "\nMatrix merger example completed successfully!" << std::endl;
        }
        //printVector("Solution x_r", x[_0]);
        //printVector("Solution x_w", x[_1]);
    }
    catch (const std::exception& e) {
        std::cerr << "Error solving system: " << e.what() << std::endl;
        //return result;
    }
    
    
    return result;
    }
}
    }
