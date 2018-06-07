/*
  Copyright 2018 SINTEF Digital, Mathematics and Cybernetics.
  Copyright 2018 Statoil Petroleum AS.

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

#ifndef OPM_LINEARSOLVERAMGCLNEW_HEADER_INCLUDED
#define OPM_LINEARSOLVERAMGCLNEW_HEADER_INCLUDED

#include <vector>
#include <memory>
#include <opm/common/ErrorMacros.hpp>
#include <opm/autodiff/DebugTimeReport.hpp>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/optional.hpp>
#include <iostream>

#include <amgcl/make_solver.hpp>
#include <amgcl/solver/bicgstab.hpp>
#include <amgcl/amg.hpp>
#include <amgcl/coarsening/smoothed_aggregation.hpp>
#include <amgcl/relaxation/spai0.hpp>
#include <amgcl/adapter/crs_tuple.hpp>
#include <amgcl/runtime.hpp>
//#include <amgcl/solver/runtime.hpp>
#include <amgcl/relaxation/runtime.hpp>
#include <amgcl/preconditioner/dummy.hpp>
#include <amgcl/adapter/crs_tuple.hpp>
#include <amgcl/adapter/zero_copy.hpp>
#include <amgcl/backend/builtin.hpp>
#include <amgcl/preconditioner/runtime.hpp>
#include <amgcl/preconditioner/cpr.hpp>
#include <amgcl/preconditioner/cpr_drs.hpp>

namespace Opm
{
struct AmgOptions
{
    int coarsen_id = 3;
    int coarse_enough = -1;
    bool direct_coarse = false;
    int max_levels = -1;
    int ncycle = 1;
    int npre = 1;// should be one according to stein
    int npost =1;
    int pre_cycles = 1;
};
    class LinearSolverAmgclNew
    {
    public:
        LinearSolverAmgclNew(int block_size, bool use_cpr_drs=true):
            block_size_(block_size),
            use_cpr_drs_(use_cpr_drs),
            solve_time_(1e99),
            setup_time_(0.0)
        {

        }
        void init(const int sz,
                          const std::vector<int>& ptr,
                          const std::vector<int>& col,
                          const std::vector<double>& val,
                          const std::vector<double>& rhs,
                          const double tolerance,
                          const int maxiter){
            /*if(block_size_>2){
                hackScalingFactors(sz, ptr, const_cast<std::vector<double>&>(val), const_cast<std::vector<double>&>(rhs));
            }*/


            auto x = new DebugTimeReport("setup all");
            time::StopWatch clock;
            clock.start();
            // setupPreconditioner(prm);
            //const int block_size = 2;
            const double dd = 0.3;
            const double ps = 0.03;
            prm_.put("precond.block_size", block_size_);
            prm_.put("precond.active_rows", sz);
            AmgOptions opts;
            setCoarseningAMGCL("precond.pprecond.", opts, prm_);
            setRelaxationAMGCL("precond.pprecond.relax.", 3, prm_);
            setRelaxationAMGCL("precond.sprecond.", 3, prm_);

            prm_.put("solver.type", amgcl::runtime::solver::bicgstab);
            prm_.put("solver.tol", tolerance);
            prm_.put("solver.maxiter", maxiter);
            //setupSolver(tolerance, maxiter, prm);
            //if(use_cpr_drs_){
            prm_.put("precond.eps_dd", dd);
            prm_.put("precond.eps_ps", ps);

            solver_ = std::make_shared<solver_t>(boost::tie(sz, ptr, col, val), prm_);
            std::ofstream file("amg_setup_cpr_drs.json");
            // solveRegular(sz, ptr, col, val, rhs, tolerance, maxiter, sol, iters, error);
            //solveCPR(sz, ptr, col, val, rhs, tolerance, maxiter, sol, iters, error);
            setup_time_ = clock.secsSinceStart();
            delete x;
            boost::property_tree::json_parser::write_json(file, prm_);
        }

        void updatePre(const int sz,
                        const std::vector<int>& ptr,
                        const std::vector<int>& col,
                        const std::vector<double>& val,
                        const double tolerance,
                        const int maxiter){
            auto x = new DebugTimeReport("update smoother");
            //if(block_size_>2){
            //    hackScalingFactors(sz, ptr, const_cast<std::vector<double>&>(val), const_cast<std::vector<double>&>(rhs));
            //}
            //solver_.reset();
            //solver_ = std::make_shared<solver_t>(boost::tie(sz, ptr, col, val), prm_);
            //solver_->precond().initUpdate(boost::tie(sz, ptr, col, val));
            solver_->precond().updateSmoother(boost::tie(sz, ptr, col, val));
            //solver_->precond()
            delete x;
        }

        void solve( const std::vector<double>& rhs,
                    std::vector<double>& sol,
                    int& iters,
                    double& error)
        {
                auto y = new DebugTimeReport("solution");
                time::StopWatch clock;
                clock.start();
                boost::tie(iters, error) = (*solver_)(rhs, sol);
                std::cout << (*solver_) << std::endl;
                std::cout << "Iterations" << iters << std::endl;
                std::cout << "Errors " << error << std::endl;
                solve_time_ = clock.secsSinceStart();
                delete y;
            /*}else{
                using Backend = amgcl::backend::builtin<double>;
                using PPrecond = amgcl::runtime::amg<Backend>;
                using SPrecond = amgcl::runtime::relaxation::as_preconditioner<Backend>;
                using Precond = amgcl::preconditioner::cpr<PPrecond, SPrecond>;
                using IterativeSolver = amgcl::runtime::iterative_solver<Backend>;
                using Solver = amgcl::make_solver<Precond, IterativeSolver>;
                auto x = new DebugTimeReport("setup");
                Solver solve(boost::tie(sz, ptr, col, val), prm);
                std::ofstream file("amg_setup_cpr.json");
                boost::property_tree::json_parser::write_json(file, prm);
                delete x;
                auto y = new DebugTimeReport("solution");
                boost::tie(iters, error) = solve(rhs, sol);
                std::cout << solve << std::endl;
                delete y;
            }*/
        }
        void hackScalingFactors(const int sz,
                                const std::vector<int>& ptr,
                                std::vector<double>& val,
                                std::vector<double>& rhs){
            //const int block_size = 2;
            const double b_gas = 100.0;
            const int n = sz / block_size_;
            assert(n*block_size_ == sz);
            for (int blockrow = 0; blockrow < n; ++blockrow) {
                const int gasrow = block_size_*blockrow + 2;
                for (int jj = ptr[gasrow]; jj < ptr[gasrow + 1]; ++jj) {
                    val[jj] /= b_gas;
                }
                rhs[gasrow] /= b_gas;
            }
        }
    private:
        void setCoarseningAMGCL(const std::string& prefix, const AmgOptions& options, boost::property_tree::ptree& prm)
        {
            std::string coarsetype = prefix + "coarsening.type";
            switch(options.coarsen_id) {
            case 1:
                prm.put(coarsetype,  amgcl::runtime::coarsening::smoothed_aggregation);
                break;
            case 2:
                prm.put(coarsetype,  amgcl::runtime::coarsening::ruge_stuben);
                break;
            case 3:
                prm.put(coarsetype,  amgcl::runtime::coarsening::aggregation);
                break;
            case 4:
                prm.put(coarsetype,  amgcl::runtime::coarsening::smoothed_aggr_emin);
                break;
            default:
                OPM_THROW(std::runtime_error, "Unknown coarsening id: " << options.coarsen_id);
            }
            // When is a level coarse enough
            if (options.coarse_enough >= 0){
                prm.put(prefix + "coarse_enough", options.coarse_enough);
            }
            // Use direct solver for coarse sys
            prm.put(prefix + "direct_coarse", options.direct_coarse);
            // Max levels
            if (options.max_levels >= 0){
                prm.put(prefix + "max_levels", options.max_levels);
            }
            // Number of cycles
            if (options.ncycle >= 0){
                prm.put(prefix + "ncycle", options.ncycle);
            }
            // Pre cycles
            if (options.npre >= 0){
                prm.put(prefix + "npre", options.npre);
            }
            // Post cycles
            if (options.npost >= 0){
                prm.put(prefix + "npost", options.npost);
            }
            // Pre cycles (precond)
            if (options.pre_cycles >= 0){
                prm.put(prefix + "pre_cycles", options.pre_cycles);
            }
        }

        void setRelaxationAMGCL(const std::string& relaxParam, const int relax_id, boost::property_tree::ptree& prm)
        {
            std::string relaxType = relaxParam + "type";
            switch(relax_id) {
            case 1:
                prm.put(relaxType,  amgcl::runtime::relaxation::spai0);
                break;
            case 2:
                prm.put(relaxType,  amgcl::runtime::relaxation::gauss_seidel);
                break;
            case 3:
                prm.put(relaxType,  amgcl::runtime::relaxation::ilu0);
                break;
            case 4:
                prm.put(relaxType,  amgcl::runtime::relaxation::iluk);
                break;
            case 5:
                prm.put(relaxType,  amgcl::runtime::relaxation::ilut);
                break;
            case 6:
                prm.put(relaxType,  amgcl::runtime::relaxation::damped_jacobi);
                break;
            case 7:
                prm.put(relaxType,  amgcl::runtime::relaxation::spai1);
                break;
            case 8:
                prm.put(relaxType,  amgcl::runtime::relaxation::chebyshev);
                break;
            default:
                OPM_THROW(std::runtime_error, "Unknown relaxation type: " << relax_id);
            }
        }

        void setupPreconditioner(boost::property_tree::ptree& prm)
        {
            // Use AMG.
            prm.put("precond.class", amgcl::runtime::precond_class::amg);
            AmgOptions opts;
            opts.coarsen_id = 3; // aggregation
            setCoarseningAMGCL("precond.", opts, prm);
            setRelaxationAMGCL("precond.relax.", 1, prm);
            setRelaxationAMGCL("precond.relax.", 3, prm);

            // using Precond = amgcl::runtime::relaxation::as_preconditioner<Backend>;
            // relaxParam = "precond.";
            // prm.put("precond.class", amgcl::runtime::precond_class::relaxation);
            // relaxParam = "precond.relax.";
            // prm.put("precond.class", amgcl::runtime::precond_class::dummy);
        }

        void setupSolver(const double tolerance,
                         const int maxiter,
                         boost::property_tree::ptree& prm)
        {
            prm.put("solver.type", amgcl::runtime::solver::bicgstab);
            prm.put("solver.tol", tolerance);
            prm.put("solver.maxiter", maxiter);
        }


        int block_size_;
        bool use_cpr_drs_;
        using Backend = amgcl::backend::builtin<double>;
        using PPrecond = amgcl::runtime::amg<Backend>;
        using SPrecond = amgcl::runtime::relaxation::as_preconditioner<Backend>;
        using Precond = amgcl::preconditioner::cpr_drs<PPrecond, SPrecond>;
        using IterativeSolver = amgcl::runtime::iterative_solver<Backend>;
        using solver_t = amgcl::make_solver<Precond, IterativeSolver>;
        //using solver_t = amgcl:<Precond, IterativeSolver>;
        std::shared_ptr<solver_t> solver_;
        boost::property_tree::ptree prm_;
        double solve_time_;
        double setup_time_;
    };

} // namespace Opm

#endif // OPM_LINEARSOLVERAMGCL_HEADER_INCLUDED
