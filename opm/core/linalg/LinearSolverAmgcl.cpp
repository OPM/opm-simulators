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

#include <opm/core/linalg/LinearSolverAmgcl.hpp>

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
//#include <amgcl/runtime.hpp>
#include <amgcl/solver/runtime.hpp>
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




  void LinearSolverAmgcl::solveRegular(const int sz,
				       const std::vector<int>& ptr,
				       const std::vector<int>& col,
				       const std::vector<double>& val,
				       const std::vector<double>& rhs,
				       std::vector<double>& sol,
				       int& iters,
				       double& error)
  {
    // Using dynamic choice of preconditioners and solvers, all
    // setup information is therefore put into a property_tree.
    const boost::property_tree::ptree& prm = prm_;
    using Backend = amgcl::backend::builtin<double>;
    using Precond = amgcl::runtime::preconditioner<Backend>;
    using IterativeSolver = amgcl::runtime::solver::wrapper<Backend>;
    using Solver = amgcl::make_solver<Precond, IterativeSolver>;
    Solver solve(std::tie(sz, ptr, col, val), prm);
            std::tie(iters, error) = solve(rhs, sol);
  }


    void LinearSolverAmgcl::hackScalingFactors(const int sz,
                            const std::vector<int>& ptr,
                            std::vector<double>& val,
                            std::vector<double>& rhs)
    {
        //const int block_size = 2;
        const double b_gas = 100.0;
    	int block_size = prm_.get<int>("block_size");
        const int n = sz / block_size;
        assert(n*block_size == sz);
        for (int blockrow = 0; blockrow < n; ++blockrow) {
            const int gasrow = block_size*blockrow + 2;
            for (int jj = ptr[gasrow]; jj < ptr[gasrow + 1]; ++jj) {
                val[jj] /= b_gas;
            }
            rhs[gasrow] /= b_gas;
        }
    }

    void LinearSolverAmgcl::solveCPR(const int sz,
                  const std::vector<int>& ptr,
                  const std::vector<int>& col,
                  const std::vector<double>& val,
                  const std::vector<double>& rhs,
                  std::vector<double>& sol,
                  int& iters,
                  double& error)
    {
        boost::property_tree::ptree prm = prm_;
	bool use_cpr_drs = (prm.get<std::string>("solver_type") == "cpr_drs");
	prm.erase("solver_type");
	prm.erase("block_size");
	prm.erase("verbose");
        if(use_cpr_drs){
	  //prm.put("precond.eps_dd", dd);
          //  prm.put("precond.eps_ps", ps);
            using Backend = amgcl::backend::builtin<double>;
	    using PPrecond = amgcl::amg<Backend, amgcl::runtime::coarsening::wrapper, amgcl::runtime::relaxation::wrapper>;
            using SPrecond = amgcl::relaxation::as_preconditioner<Backend, amgcl::runtime::relaxation::wrapper>;
            using Precond = amgcl::preconditioner::cpr_drs<PPrecond, SPrecond>;
	    using IterativeSolver = amgcl::runtime::solver::wrapper<Backend>;
            using Solver = amgcl::make_solver<Precond, IterativeSolver>;
            auto x = new DebugTimeReport("setup");
            Solver solve(std::tie(sz, ptr, col, val), prm);
            std::ofstream file("amg_setup_cpr_drs.json");
            boost::property_tree::json_parser::write_json(file, prm);
            delete x;
            auto y = new DebugTimeReport("solution");
            std::tie(iters, error) = solve(rhs, sol);
	    if(prm_.get<bool>("verbose")){
	      std::cout << solve << std::endl;
	    }
            delete y;
        }else{
            using Backend = amgcl::backend::builtin<double>;
	    using PPrecond = amgcl::amg<Backend, amgcl::runtime::coarsening::wrapper, amgcl::runtime::relaxation::wrapper>;
            using SPrecond = amgcl::relaxation::as_preconditioner<Backend, amgcl::runtime::relaxation::wrapper>;
            using Precond = amgcl::preconditioner::cpr<PPrecond, SPrecond>;
	    using IterativeSolver = amgcl::runtime::solver::wrapper<Backend>;
            using Solver = amgcl::make_solver<Precond, IterativeSolver>;
            auto x = new DebugTimeReport("setup");
            Solver solve(std::tie(sz, ptr, col, val), prm);
            std::ofstream file("amg_setup_cpr.json");
            boost::property_tree::json_parser::write_json(file, prm);
            delete x;
            auto y = new DebugTimeReport("solution");
            std::tie(iters, error) = solve(rhs, sol);
	    if(prm_.get<bool>("verbose")){
	      std::cout << solve << std::endl;
	    }
            delete y;
        }
    }


    void LinearSolverAmgcl::solve(const int sz,
                                  const std::vector<int>& ptr,
                                  const std::vector<int>& col,
                                  const std::vector<double>& val,
                                  const std::vector<double>& rhs,
                                  std::vector<double>& sol,
                                  int& iters,
                                  double& error)
    {
      int block_size = prm_.get<int>("block_size");
      if(block_size>2){
	hackScalingFactors(sz, ptr, const_cast<std::vector<double>&>(val), const_cast<std::vector<double>&>(rhs));
        }
        DebugTimeReport rep("amgcl-timer");
        // solveRegular(sz, ptr, col, val, rhs, tolerance, maxiter, sol, iters, error);
        solveCPR(sz, ptr, col, val, rhs, sol, iters, error);
    }

} // namespace Opm
