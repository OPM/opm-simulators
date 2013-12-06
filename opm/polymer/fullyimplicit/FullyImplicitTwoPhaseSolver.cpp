/**/

#include <opm/autodiff/polymer/FullyImplicitTwoPhaseSolver.hpp>

#include <opm/core/pressure/tpfa/trans_tpfa.h>
#include <opm/autodiff/AutoDiffBlock.hpp>
#include <opm/autodiff/AutoDiffHelpers.hpp>
#include <opm/autodiff/polymer/IncompPropsAdInterface.hpp>

#include <opm/core/grid.h>
#include <opm/core/linalg/LinearSolverInterface.hpp>
#include <opm/core/props/rock/RockCompressibility.hpp>
#include <opm/core/simulator/TwophaseState.hpp>
#include <opm/core/simulator/WellState.hpp>
#include <opm/core/utility/ErrorMacros.hpp>

#include <cassert>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <Eigen/Eigen>
#include <algorithm>
namespace Opm {

namespace {
    
    std::vector<int>
    buildAllCells(const int nc)
    {
        std::vector<int> all_cells(nc);
        for (int c = 0; c < nc; ++c) { all_cells[c] = c; }

        return all_cells;
    }
    struct Chop01 {
        double operator()(double x) const { return std::max(std::min(x, 1.0), 0.0); }
    };

}//anonymous namespace





typedef AutoDiffBlock<double> ADB;
typedef ADB::V V;
typedef ADB::M M;
typedef Eigen::Array<double,
                     Eigen::Dynamic,
                     Eigen::Dynamic,
                     Eigen::RowMajor> DataBlock;





    FullyImplicitTwoPhaseSolver::
    FullyImplicitTwoPhaseSolver(const UnstructuredGrid&         grid,
                                const IncompPropsAdInterface&   fluid,
                                const LinearSolverInterface&    linsolver)
        : grid_ (grid)
        , fluid_(fluid)
        , linsolver_(linsolver)
        , cells_ (buildAllCells(grid.number_of_cells))
        , ops_(grid)
        , residual_(std::vector<ADB>(fluid.numPhases(), ADB::null()))
     {
     }





    void
    FullyImplicitTwoPhaseSolver::
    step(const double   dt,
         TwophaseState& x,
         const std::vector<double>& src)
    {
        
        V pvol(grid_.number_of_cells);
        // Pore volume
        const typename V::Index nc = grid_.number_of_cells;
        V rho = V::Constant(pvol.size(), 1, *fluid_.porosity());
        std::transform(grid_.cell_volumes, grid_.cell_volumes + nc,
                       rho.data(), pvol.data(),
                       std::multiplies<double>());

        const V pvdt = pvol / dt;

        const SolutionState old_state = constantState(x);
        const double atol  = 1.0e-12;
        const double rtol  = 5.0e-8;
        const int    maxit = 15;

        assemble(pvdt, old_state, x, src);

        const double r0  = residualNorm();
        int          it  = 0;
        std::cout << "\nIteration         Residual\n"
                  << std::setw(9) << it << std::setprecision(9)
                  << std::setw(18) << r0 << std::endl;
        bool resTooLarge = r0 > atol;
        while (resTooLarge && (it < maxit)) {
            const V dx = solveJacobianSystem();
            updateState(dx, x);

            assemble(pvdt, old_state, x, src);

            const double r = residualNorm();

            resTooLarge = (r > atol) && (r > rtol*r0);

            it += 1;
            std::cout << std::setw(9) << it << std::setprecision(9)
                      << std::setw(18) << r << std::endl;
        }

        if (resTooLarge) {
            std::cerr << "Failed to compute converged solution in " << it << " iterations. Ignoring!\n";
            // OPM_THROW(std::runtime_error, "Failed to compute converged solution in " << it << " iterations.");
        }
    }





    FullyImplicitTwoPhaseSolver::SolutionState::SolutionState(const int np)
        : pressure   (    ADB::null())
        , saturation (np, ADB::null())
    {
    }





    FullyImplicitTwoPhaseSolver::SolutionState
    FullyImplicitTwoPhaseSolver::constantState(const TwophaseState& x)
    {
        const int nc = grid_.number_of_cells;
        const int np = x.numPhases();

        SolutionState state(np);

        // Pressure.
        assert (not x.pressure().empty());
        const V p = Eigen::Map<const V>(& x.pressure()[0], nc);
        state.pressure = ADB::constant(p);

        // Saturation.
        assert (not x.saturation().empty());
        const DataBlock s_all = Eigen::Map<const DataBlock>(& x.saturation()[0], nc, np);
        for (int phase = 0; phase < np; ++phase) {
            state.saturation[phase] = ADB::constant(s_all.col(phase));
   //     state.saturation[1] = ADB::constant(s_all.col(1));
        }
        return state;
    }





    FullyImplicitTwoPhaseSolver::SolutionState
    FullyImplicitTwoPhaseSolver::variableState(const TwophaseState& x)
    {
        const int nc = grid_.number_of_cells;
        const int np = x.numPhases();

        std::vector<V> vars0;
        vars0.reserve(np); 

        // Initial pressure.
        assert (not x.pressure().empty());
        const V p = Eigen::Map<const V>(& x.pressure()[0], nc);
        vars0.push_back(p);

        // Initial saturation.
        assert (not x.saturation().empty());
        const DataBlock s_all = Eigen::Map<const DataBlock>(& x.saturation()[0], nc, np);
        const V sw = s_all.col(0);
        vars0.push_back(sw);


        std::vector<ADB> vars = ADB::variables(vars0);

        SolutionState state(np);

        // Pressure.
        int nextvar = 0;
        state.pressure = vars[ nextvar++ ];

        // Saturation.
        const std::vector<int>& bpat = vars[0].blockPattern();
        {
            ADB so  = ADB::constant(V::Ones(nc, 1), bpat);
            ADB& sw = vars[ nextvar++ ];
            state.saturation[0] = sw;
            so = so - sw;
            state.saturation[1] = so;
        }

        assert(nextvar == int(vars.size()));

        return state;
    }
  




    void
    FullyImplicitTwoPhaseSolver::
    assemble(const V&             pvdt,
             const SolutionState& old_state,
             const TwophaseState& x   ,
             const std::vector<double>&  src)
    {
        // Create the primary variables.
        const SolutionState state = variableState(x);

        // -------- Mass balance equations --------
        const V trans = subset(transmissibility(), ops_.internal_faces);
        const std::vector<ADB> kr = computeRelPerm(state);
        for (int phase = 0; phase < fluid_.numPhases(); ++phase) {
            const ADB mflux = computeMassFlux(phase, trans, kr, state);
            ADB source = accumSource(phase, kr, src);
            residual_[phase] =
                pvdt*(state.saturation[phase] - old_state.saturation[phase])
                + ops_.div*mflux - source;
        }

    }
   





    ADB 
    FullyImplicitTwoPhaseSolver::accumSource(const int phase,
                                                 const std::vector<ADB>& kr,
                                                 const std::vector<double>& src) const
    {
        //extract the source to out and in source.
        std::vector<double> outsrc;
        std::vector<double> insrc;
        std::vector<double>::const_iterator it;
        for (it = src.begin(); it != src.end(); ++it) {
            if (*it < 0) {
                outsrc.push_back(*it);
                insrc.push_back(0.0);
            } else if (*it > 0) {
                insrc.push_back(*it);
                outsrc.push_back(0.0);
            } else {
                outsrc.emplace_back(0);
                insrc.emplace_back(0);
            }
        }
        const V source = Eigen::Map<const V>(& src[0], grid_.number_of_cells);
        const V outSrc = Eigen::Map<const V>(& outsrc[0], grid_.number_of_cells);
        const V inSrc = Eigen::Map<const V>(& insrc[0], grid_.number_of_cells);
        
        // compute the out-fracflow.
        ADB f_out = computeFracFlow(phase, kr);
        // compute the in-fracflow.
        V f_in;
        if (phase == 1) {
            f_in = V::Zero(grid_.number_of_cells);
        } else if (phase == 0) {
            f_in = V::Ones(grid_.number_of_cells);
        }
        return f_out * outSrc + f_in * inSrc;
     }





    ADB
    FullyImplicitTwoPhaseSolver::computeFracFlow(int                     phase,
                                                 const std::vector<ADB>& kr) const
    {
        const double* mus = fluid_.viscosity();
        ADB  mob_phase = kr[phase] / V::Constant(kr[phase].size(), 1, mus[phase]);
        ADB  mob_wat = kr[0] / V::Constant(kr[0].size(), 1, mus[0]);
        ADB  mob_oil= kr[1] / V::Constant(kr[1].size(), 1, mus[1]);
        ADB  total_mob = mob_wat + mob_oil;
        ADB f = mob_phase / total_mob;

        return f;
    }





    V 
    FullyImplicitTwoPhaseSolver::solveJacobianSystem() const
    {
        const int np = fluid_.numPhases();
    	if (np != 2) {
	        OPM_THROW(std::logic_error, "Only two-phase ok in FullyImplicitTwoPhaseSolver.");
	    }
	    ADB mass_res = collapseJacs(vertcat(residual_[0], residual_[1]));
        const Eigen::SparseMatrix<double, Eigen::RowMajor> matr = mass_res.derivative()[0];
        V dx(V::Zero(mass_res.size()));
        Opm::LinearSolverInterface::LinearSolverReport rep
            = linsolver_.solve(matr.rows(), matr.nonZeros(),
                               matr.outerIndexPtr(), matr.innerIndexPtr(), matr.valuePtr(),
                               mass_res.value().data(), dx.data());
        if (!rep.converged) {
            OPM_THROW(std::runtime_error,
                      "FullyImplicitBlackoilSolver::solveJacobianSystem(): "
                      "Linear solver convergence failure.");
        }
        return dx;
    }





    void FullyImplicitTwoPhaseSolver::updateState(const V& dx,
                                                  TwophaseState& state) const
    {
        const int np = fluid_.numPhases();
        const int nc = grid_.number_of_cells;
        const V null;
        assert(null.size() == 0);
        const V zero = V::Zero(nc);
        const V one = V::Constant(nc, 1.0);

        // Extract parts of dx corresponding to each part.
        const V dp = subset(dx, Span(nc));
        int varstart = nc;
        const V dsw = subset(dx, Span(nc, 1, varstart));
        varstart += dsw.size();
        assert(varstart == dx.size());

        // Pressure update.
        const V p_old = Eigen::Map<const V>(&state.pressure()[0], nc);
        const V p = p_old - dp;
        std::copy(&p[0], &p[0] + nc, state.pressure().begin());


        // Saturation updates.
        const double dsmax = 0.3;
        const DataBlock s_old = Eigen::Map<const DataBlock>(& state.saturation()[0], nc, np);
        V so = one;
        const V sw_old = s_old.col(0);
        const V dsw_limited = sign(dsw) * dsw.abs().min(dsmax);
        const V sw = (sw_old - dsw_limited).unaryExpr(Chop01());
        so -= sw;
        for (int c = 0; c < nc; ++c) {
            state.saturation()[c*np] = sw[c];
        }
        for (int c = 0; c < nc; ++c) {
            state.saturation()[c*np + 1] = so[c];
        }

    }
   




    std::vector<ADB>
    FullyImplicitTwoPhaseSolver::computeRelPerm(const SolutionState& state) const
    {

        const ADB sw = state.saturation[0];
        const ADB so = state.saturation[1];

        return fluid_.relperm(sw, so, cells_);
    }
    
    
    
    
    
    
    
    
    
    
    ADB
    FullyImplicitTwoPhaseSolver::computeMassFlux(const int               phase ,
                                                 const V&                trans,
                                                 const std::vector<ADB>& kr    ,
                                                 const SolutionState&    state ) const
    {
//        const ADB tr_mult = transMult(state.pressure);
        const double* mus = fluid_.viscosity();
        ADB  mob = kr[phase] / V::Constant(kr[phase].size(), 1, mus[phase]);

        const ADB dp = ops_.ngrad * state.pressure;
        const ADB head = trans * dp;

        UpwindSelector<double> upwind(grid_, ops_, head.value());

        return upwind.select(mob) * head;
    }
  

   
   
   
    double
    FullyImplicitTwoPhaseSolver::residualNorm() const
    {
        double r = 0;
        for (std::vector<ADB>::const_iterator
                 b = residual_.begin(),
                 e = residual_.end();
             b != e; ++b)
        {
            r = std::max(r, (*b).value().matrix().norm());
        }

        return r;
    }
   
   
   
   
   
    V
    FullyImplicitTwoPhaseSolver::transmissibility() const
    {
        const V::Index nc = grid_.number_of_cells;
        V htrans(grid_.cell_facepos[nc]);
        V trans(grid_.cell_facepos[nc]);
        UnstructuredGrid* ug = const_cast<UnstructuredGrid*>(& grid_);
        tpfa_htrans_compute(ug, fluid_.permeability(), htrans.data());
        tpfa_trans_compute (ug, htrans.data()     , trans.data());
        
        return trans;
    }






}//namespace Opm
