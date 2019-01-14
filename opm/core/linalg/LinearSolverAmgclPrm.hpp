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
    class LinearSolverAmgclNew
    {
    public:
      LinearSolverAmgclNew(const boost::property_tree::ptree& prm):
	prm_(prm),
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
            auto x = new DebugTimeReport("setup all");
            time::StopWatch clock;
            clock.start();
	    int block_size = prm_.get<int>("block_size"):
            prm_.put("precond.block_size", block_size);
            prm_.put("precond.active_rows", sz);
            //prm_.put("solver.tol", tolerance);
            //prm_.put("solver.maxiter", maxiter);
            solver_ = std::make_shared<solver_t>(std::tie(sz, ptr, col, val), prm_);
            std::ofstream file("amg_setup_cpr_drs_debug.json");
            setup_time_ = clock.secsSinceStart();
            delete x;
            std::property_tree::json_parser::write_json(file, prm_);
        }

        void updatePre(const int sz,
                        const std::vector<int>& ptr,
                        const std::vector<int>& col,
                        const std::vector<double>& val){
            auto x = new DebugTimeReport("update smoother");
            solver_->precond().updateSmoother(std::tie(sz, ptr, col, val));
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
                std::tie(iters, error) = (*solver_)(rhs, sol);
                std::cout << (*solver_) << std::endl;
                std::cout << "Iterations" << iters << std::endl;
                std::cout << "Errors " << error << std::endl;
                solve_time_ = clock.secsSinceStart();
                delete y;            
        }
        void hackScalingFactors(const int sz,
                                const std::vector<int>& ptr,
                                std::vector<double>& val,
                                std::vector<double>& rhs){
            //const int block_size = 2;
	    int block_size = prm_.get<double>("block_size"):
            const double b_gas = 100.0;
            const int n = sz / block_size;
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
 
      using Backend = amgcl::backend::builtin<double>;
      using PPrecond = amgcl::amg<Backend, amgcl::runtime::coarsening::wrapper, amgcl::runtime::relaxation::wrapper>;
      using SPrecond = amgcl::relaxation::as_preconditioner<Backend, amgcl::runtime::relaxation::wrapper>;
      using Precond = amgcl::preconditioner::cpr_drs<PPrecond, SPrecond>;
      using IterativeSolver = amgcl::runtime::solver::wrapper<Backend>;
      using solver_t = amgcl::make_solver<Precond, IterativeSolver>;
      //using solver_t = amgcl:<Precond, IterativeSolver>;
      std::shared_ptr<solver_t> solver_;
      boost::property_tree::ptree prm_;
      double solve_time_;
      double setup_time_;
    };

} // namespace Opm

#endif // OPM_LINEARSOLVERAMGCL_HEADER_INCLUDED
