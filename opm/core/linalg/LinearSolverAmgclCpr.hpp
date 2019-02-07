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

#ifndef OPM_LINEARSOLVERAMGCLCPR_HEADER_INCLUDED
#define OPM_LINEARSOLVERAMGCLCPR_HEADER_INCLUDED

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
  class LinearSolverAmgclCpr
  {
  public:
    LinearSolverAmgclCpr(const boost::property_tree::ptree& prm):
      prm_(prm),
      solve_time_(1e99),
      setup_time_(0.0)
    {
      
    }
      void init(const int sz,
		const std::vector<int>& ptr,
		const std::vector<int>& col,
		const std::vector<double>& val){
	auto x = new DebugTimeReport("setup all");
	time::StopWatch clock;
	boost::property_tree::ptree prm = prm_;
	use_drs_ = prm.get<bool>("use_drs");
	prm.erase("solver_type");
	prm.erase("use_drs");
	prm.erase("block_size");
	prm.erase("verbose");
	if(use_drs_){	    
	  solver_cpr_drs_.reset(new solver_cpr_drs_t(std::tie(sz, ptr, col, val), prm));
	  if(prm_.get<bool>("verbose")){
	    std::ofstream file("amg_setup_cpr_drs.json");
	    boost::property_tree::json_parser::write_json(file, prm);
	  }
	}else{
	  solver_cpr_.reset(new solver_cpr_t(std::tie(sz, ptr, col, val), prm));
	  if(prm_.get<bool>("verbose")){
	    std::ofstream file("amg_setup_cpr.json");
	    boost::property_tree::json_parser::write_json(file, prm);
	  }
	}
	delete x;
      }

      // void updatePre(const int sz,
      // 		     const std::vector<int>& ptr,
      // 		     const std::vector<int>& col,
      // 		     const std::vector<double>& val){
      // 	auto x = new DebugTimeReport("update smoother");
      // 	solver_->precond().updateSmoother(std::tie(sz, ptr, col, val));
      // 	//solver_->precond()
      // 	delete x;
      // }

        void solve( const std::vector<double>& rhs,
                    std::vector<double>& sol,
                    int& iters,
                    double& error)
        {
                auto y = new DebugTimeReport("solution");
                time::StopWatch clock;
                clock.start();
		if(use_drs_){
		  std::tie(iters, error) = (*solver_cpr_drs_)(rhs, sol);
		  if(prm_.get<bool>("verbose")){
		    std::cout << (*solver_cpr_drs_) << std::endl;
		  }
		}else{
		  std::tie(iters, error) = (*solver_cpr_)(rhs, sol);
		  if(prm_.get<bool>("verbose")){
		    std::cout << (*solver_cpr_) << std::endl;
		  }
		}
		if(prm_.get<bool>("verbose")){
		  std::cout << "Iterations" << iters << std::endl;
		  std::cout << "Errors " << error << std::endl;
		}
                solve_time_ = clock.secsSinceStart();
                delete y;            
        }
      
    private:
      using Backend = amgcl::backend::builtin<double>;
      using PPrecond = amgcl::amg<Backend, amgcl::runtime::coarsening::wrapper, amgcl::runtime::relaxation::wrapper>;
      using SPrecond = amgcl::relaxation::as_preconditioner<Backend, amgcl::runtime::relaxation::wrapper>;
      using Precond_cpr_drs = amgcl::preconditioner::cpr_drs<PPrecond, SPrecond>;
      using IterativeSolver = amgcl::runtime::solver::wrapper<Backend>;
      using solver_cpr_drs_t = amgcl::make_solver<Precond_cpr_drs, IterativeSolver>;
      using Precond_cpr = amgcl::preconditioner::cpr<PPrecond, SPrecond>;
      //using solver_t = amgcl:<Precond, IterativeSolver>;
      using solver_cpr_t = amgcl::make_solver<Precond_cpr, IterativeSolver>;
      std::shared_ptr<solver_cpr_drs_t> solver_cpr_drs_;
      std::shared_ptr<solver_cpr_t> solver_cpr_;
      boost::property_tree::ptree prm_;
      double solve_time_;
      double setup_time_;
      bool use_drs_;
    };

} // namespace Opm

#endif // OPM_LINEARSOLVERAMGCL_HEADER_INCLUDED
