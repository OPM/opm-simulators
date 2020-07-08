/*
  Copyright 2019 SINTEF Digital, Mathematics and Cybernetics.

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

#include <config.h>

#include <opm/simulators/linalg/setupPropertyTree.hpp>

#include <boost/version.hpp>
#include <boost/property_tree/json_parser.hpp>

namespace Opm
{

/// Set up a property tree intended for FlexibleSolver by either reading
/// the tree from a JSON file or creating a tree giving the default solver
/// and preconditioner. If the latter, the parameters --linear-solver-reduction,
/// --linear-solver-maxiter and --linear-solver-verbosity are used, but if reading
/// from file the data in the JSON file will override any other options.
template<class TypeTag>
boost::property_tree::ptree
setupPropertyTree(const FlowLinearSolverParameters& p)
{
    boost::property_tree::ptree prm;
    if (p.linear_solver_configuration_ == "file") {
#if BOOST_VERSION / 100 % 1000 > 48
        if (p.linear_solver_configuration_json_file_ == "none"){
            OPM_THROW(std::invalid_argument,
                      "--linear-solver-configuration=file requires that a filename "
                      << "is passed with "
                      << "--linear-solver-configuration-json-file=filename.");
        }else{
            boost::property_tree::read_json(p.linear_solver_configuration_json_file_, prm);
        }
#else
       OPM_THROW(std::invalid_argument,
                 "--linear-solver-configuration=file not supported with "
                 << "boost version. Needs version > 1.48.");
#endif
    }
    else
    {
        std::string conf =  p.linear_solver_configuration_;
        // Support old UseCpr if not configuration was set
        if((conf == "cpr_trueimpes") || (conf == "cpr_quasiimpes")){
            prm.put("tol", p.linear_solver_reduction_);
            if (EWOMS_PARAM_IS_SET(TypeTag, int, LinearSolverMaxIter))
                prm.put("maxiter", p.linear_solver_maxiter_);// Trust that the user knows what he does
            else
                prm.put("maxiter", 20); // Use our own default.
            prm.put("verbosity", p.linear_solver_verbosity_);
            prm.put("solver", "bicgstab");
            prm.put("preconditioner.type", "cpr");
            if (conf == "cpr_quasiimpes" )
            {
                prm.put("preconditioner.weight_type","quasiimpes");
            }
            else
            {
                prm.put("preconditioner.weight_type","trueimpes");
            }
            prm.put("preconditioner.finesmoother.type", "ParOverILU0");
            prm.put("preconditioner.finesmoother.relaxation", 1.0);
            prm.put("preconditioner.pressure_var_index",1);
            prm.put("preconditioner.verbosity",0);
            prm.put("preconditioner.coarsesolver.maxiter",1);
            prm.put("preconditioner.coarsesolver.tol",1e-1);
            prm.put("preconditioner.coarsesolver.solver","loopsolver");
            prm.put("preconditioner.coarsesolver.verbosity",0);
            prm.put("preconditioner.coarsesolver.preconditioner.type","amg");
            prm.put("preconditioner.coarsesolver.preconditioner.alpha",0.333333333333);
            prm.put("preconditioner.coarsesolver.preconditioner.relaxation",1.0);
            if (EWOMS_PARAM_IS_SET(TypeTag, int, CprMaxEllIter))
                prm.put("preconditioner.coarsesolver.preconditioner.iterations", p.cpr_max_ell_iter_);
            else
                prm.put("preconditioner.coarsesolver.preconditioner.iterations",1);
            prm.put("preconditioner.coarsesolver.preconditioner.coarsenTarget",1200);
            prm.put("preconditioner.coarsesolver.preconditioner.pre_smooth",1);
            prm.put("preconditioner.coarsesolver.preconditioner.post_smooth",1);
            prm.put("preconditioner.coarsesolver.preconditioner.beta",1e-5);
            prm.put("preconditioner.coarsesolver.preconditioner.smoother","ILU0");
            prm.put("preconditioner.coarsesolver.preconditioner.verbosity",0);
            prm.put("preconditioner.coarsesolver.preconditioner.maxlevel",15);
            prm.put("preconditioner.coarsesolver.preconditioner.skip_isolated",0);
        } else {
            if(conf != "ilu0"){
                OPM_THROW(std::invalid_argument, conf  << "is not a valid setting for --linear-solver-configuration."
                          << " Please use ilu0, cpr_trueimpes, or cpr_quasiimpes");
            }
            prm.put("tol", p.linear_solver_reduction_);
            prm.put("maxiter", p.linear_solver_maxiter_);
            prm.put("verbosity", p.linear_solver_verbosity_);
            prm.put("solver", "bicgstab");
            prm.put("preconditioner.type", "ParOverILU0");
            prm.put("preconditioner.relaxation", p.ilu_relaxation_);
            prm.put("preconditioner.ilulevel", p.ilu_fillin_level_);
        }
    }
    return prm;
}

} // namespace Opm
