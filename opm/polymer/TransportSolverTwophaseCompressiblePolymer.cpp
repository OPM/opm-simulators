/*
  Copyright 2012 SINTEF ICT, Applied Mathematics.

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

#include <opm/polymer/TransportSolverTwophaseCompressiblePolymer.hpp>
#include <opm/core/props/BlackoilPropertiesInterface.hpp>
#include <opm/core/grid.h>
#include <opm/core/transport/reorder/reordersequence.h>
#include <opm/common/utility/numeric/RootFinders.hpp>
#include <opm/core/utility/miscUtilities.hpp>
#include <opm/core/utility/miscUtilitiesBlackoil.hpp>
#include <opm/core/pressure/tpfa/trans_tpfa.h>
#include <opm/common/ErrorMacros.hpp>
#include <cmath>
#include <list>
#include <iostream>
// Choose error policy for scalar solves here.
typedef Opm::RegulaFalsi<Opm::WarnAndContinueOnError> RootFinder;


class Opm::TransportSolverTwophaseCompressiblePolymer::ResidualEquation
{
public:
    GradientMethod gradient_method;
    int cell;
    double s0;
    double c0;
    double cmax0;
    double influx;  // B_i sum_j b_j min(v_ij, 0)*f(s_j) - B_i q_w
    double influx_polymer;
    double outflux; // sum_j max(v_ij, 0) - B_i q
    double porevolume0;
    double porevolume;
    double porosity0;
    double porosity;
    double B_cell0;
    double B_cell;
    double dtpv;    // dt/pv(i)
    double dps;
    double rhor;
    double ads0;

    TransportSolverTwophaseCompressiblePolymer& tm;

    ResidualEquation(TransportSolverTwophaseCompressiblePolymer& tmodel, int cell_index);
    void computeResidual(const double* x, double* res) const;
    void computeResidual(const double* x, double* res, double& mc, double& ff) const;
    double computeResidualS(const double* x) const;
    double computeResidualC(const double* x) const;
    void computeGradientResS(const double* x, double* res, double* gradient) const;
    void computeGradientResC(const double* x, double* res, double* gradient) const;
    void computeJacobiRes(const double* x, double* dres_s_dsdc, double* dres_c_dsdc) const;



private:
    void computeResAndJacobi(const double* x, const bool if_res_s, const bool if_res_c,
                             const bool if_dres_s_dsdc, const bool if_dres_c_dsdc,
                             double* res, double* dres_s_dsdc,
                             double* dres_c_dsdc, double& mc, double& ff) const;
};

class Opm::TransportSolverTwophaseCompressiblePolymer::ResidualCGrav {
public:
    const TransportSolverTwophaseCompressiblePolymer& tm;
    const int cell;
    const double s0;
    const double c0;
    const double cmax0;
    const double porevolume;
    const double porosity;
    const double dtpv;    // dt/pv(i)
    const double dps;
    const double rhor;
    double c_ads0;
    double gf[2];
    int nbcell[2];
    mutable double last_s;

    ResidualCGrav(const TransportSolverTwophaseCompressiblePolymer& tmodel,
                  const std::vector<int>& cells,
                  const int pos,
                  const double* gravflux);

    double operator()(double c) const;
    double computeGravResidualS(double s, double c) const;
    double computeGravResidualC(double s, double c) const;
    double lastSaturation() const;
};

class Opm::TransportSolverTwophaseCompressiblePolymer::ResidualSGrav {
public:
    const ResidualCGrav& res_c_eq_;
    double c;

    ResidualSGrav(const ResidualCGrav& res_c_eq, const double c_init = 0.0);
    double operator()(double s) const;
};


namespace
{
    bool check_interval(const double* xmin, const double* xmax, double* x);

    double norm(double* res)
    {
        return std::max(std::abs(res[0]), std::abs(res[1]));
    }

    // Define a piecewise linear curve along which we will look for zero of the "s" or "r" residual.
    // The curve starts at "x", goes along the direction "direction" until it hits the boundary of the box of
    // admissible values for "s" and "x" (which is given by "[x_min[0], x_max[0]]x[x_min[1], x_max[1]]").
    // Then it joins in a straight line the point "end_point".
    class CurveInSCPlane{
    public:
        CurveInSCPlane();
        void setup(const double* x, const double* direction,
                   const double* end_point, const double* x_min,
                   const double* x_max, const double tol,
                   double& t_max_out, double& t_out_out);
        void computeXOfT(double*, const double) const;

    private:
        double direction_[2];
        double end_point_[2];
        double x_max_[2];
        double x_min_[2];
        double t_out_;
        double t_max_; // t_max = t_out + 1
        double x_out_[2];
        double x_[2];
    };

}


namespace Opm
{
    TransportSolverTwophaseCompressiblePolymer::TransportSolverTwophaseCompressiblePolymer(const UnstructuredGrid& grid,
                                                                         const BlackoilPropertiesInterface& props,
                                                                         const PolymerProperties& polyprops,
                                                                         const SingleCellMethod method,
                                                                         const double tol,
                                                                         const int maxit)
        : grid_(grid),
          props_(props),
          polyprops_(polyprops),
          darcyflux_(0),
          porevolume0_(0),
          porevolume_(0),
          source_(0),
          polymer_inflow_c_(0),
          dt_(0.0),
          tol_(tol),
          maxit_(maxit),
          method_(method),
          adhoc_safety_(1.1),
          concentration_(0),
          cmax_(0),
          fractionalflow_(grid.number_of_cells, -1.0),
          mc_(grid.number_of_cells, -1.0),
          gravity_(0),
          mob_(2*grid.number_of_cells, -1.0),
          ia_upw_(grid.number_of_cells + 1, -1),
          ja_upw_(grid.number_of_faces, -1),
          ia_downw_(grid.number_of_cells + 1, -1),
          ja_downw_(grid.number_of_faces, -1)

    {
        const int np = props.numPhases();
        const int num_cells = grid.number_of_cells;
        if (props.numPhases() != 2) {
            OPM_THROW(std::runtime_error, "Property object must have 2 phases");
        }
        visc_.resize(np*num_cells);
        A_.resize(np*np*num_cells);
        A0_.resize(np*np*num_cells);
        smin_.resize(np*num_cells);
        smax_.resize(np*num_cells);
        allcells_.resize(num_cells);
        for (int i = 0; i < num_cells; ++i) {
            allcells_[i] = i;
        }
        props.satRange(num_cells, &allcells_[0], &smin_[0], &smax_[0]);
    }




    void TransportSolverTwophaseCompressiblePolymer::setPreferredMethod(SingleCellMethod method)
    {
        method_ = method;
    }




    void TransportSolverTwophaseCompressiblePolymer::solve(const double* darcyflux,
                                                  const std::vector<double>& initial_pressure,
                                                  const std::vector<double>& pressure,
                                                  const std::vector<double>& temperature,
                                                  const double* porevolume0,
                                                  const double* porevolume,
                                                  const double* source,
                                                  const double* polymer_inflow_c,
                                                  const double dt,
                                                  std::vector<double>& saturation,
                                                  std::vector<double>& surfacevol,
                                                  std::vector<double>& concentration,
                                                  std::vector<double>& cmax)
    {
        darcyflux_ = darcyflux;
        porevolume0_ = porevolume0;
        porevolume_ = porevolume;
        source_ = source;
        dt_ = dt;
        polymer_inflow_c_ = polymer_inflow_c;
        toWaterSat(saturation, saturation_);
        concentration_ = &concentration[0];
        cmax_ = &cmax[0];

#if PROFILING
        res_counts.clear();
#endif

        props_.viscosity(grid_.number_of_cells, &pressure[0], &temperature[0], NULL, &allcells_[0], &visc_[0], NULL);
        props_.matrix(grid_.number_of_cells, &initial_pressure[0], &temperature[0], NULL, &allcells_[0], &A0_[0], NULL);
        props_.matrix(grid_.number_of_cells, &pressure[0], &temperature[0], NULL, &allcells_[0], &A_[0], NULL);

        // Check immiscibility requirement (only done for first cell).
        if (A_[1] != 0.0 || A_[2] != 0.0) {
            OPM_THROW(std::runtime_error, "TransportCompressibleSolverTwophaseCompressibleTwophase requires a property object without miscibility.");
        }
        std::vector<int> seq(grid_.number_of_cells);
        std::vector<int> comp(grid_.number_of_cells + 1);
        int ncomp;
        compute_sequence_graph(&grid_, darcyflux_,
                               &seq[0], &comp[0], &ncomp,
                               &ia_upw_[0], &ja_upw_[0]);
        const int nf = grid_.number_of_faces;
        std::vector<double> neg_darcyflux(nf);
        std::transform(darcyflux, darcyflux + nf, neg_darcyflux.begin(), std::negate<double>());
        compute_sequence_graph(&grid_, &neg_darcyflux[0],
                               &seq[0], &comp[0], &ncomp,
                               &ia_downw_[0], &ja_downw_[0]);
        reorderAndTransport(grid_, darcyflux);
        toBothSat(saturation_, saturation);

        // Compute surface volume as a postprocessing step from saturation and A_
        computeSurfacevol(grid_.number_of_cells, props_.numPhases(), &A_[0], &saturation[0], &surfacevol[0]);
    }




    // Residual for saturation equation, single-cell implicit Euler transport
    //
    //     r(s) = s - s0 + dt/pv*( influx + outflux*f(s) )
    //
    // where influx is water influx, outflux is total outflux.
    // Influxes are negative, outfluxes positive.
    struct TransportSolverTwophaseCompressiblePolymer::ResidualS
    {
        TransportSolverTwophaseCompressiblePolymer::ResidualEquation& res_eq_;
        const double c_;
        explicit ResidualS(TransportSolverTwophaseCompressiblePolymer::ResidualEquation& res_eq,
                           const double c)
            : res_eq_(res_eq),
              c_(c)
        {
        }

        double operator()(double s) const
        {
            double x[2];
            x[0] = s;
            x[1] = c_;
            return res_eq_.computeResidualS(x);
        }
    };

    // Residual for concentration equation, single-cell implicit Euler transport
    //
    //  \TODO doc me
    // where ...
    // Influxes are negative, outfluxes positive.
    struct TransportSolverTwophaseCompressiblePolymer::ResidualC
    {
        mutable double s; // Mutable in order to change it with every operator() call to be the last computed s value.
        TransportSolverTwophaseCompressiblePolymer::ResidualEquation& res_eq_;
        explicit ResidualC(TransportSolverTwophaseCompressiblePolymer::ResidualEquation& res_eq)
            : res_eq_(res_eq)
        {}

        void computeBothResiduals(const double s_arg, const double c_arg, double& res_s, double& res_c, double& mc, double& ff) const
        {
            double x[2];
            double res[2];
            x[0] = s_arg;
            x[1] = c_arg;
            res_eq_.computeResidual(x, res, mc, ff);
            res_s = res[0];
            res_c = res[1];
        }

        double operator()(double c) const
        {
            ResidualS res_s(res_eq_, c);
            int iters_used;
            // Solve for s first.
            // s = modifiedRegulaFalsi(res_s, std::max(tm.smin_[2*cell], dps), tm.smax_[2*cell],
            //                      tm.maxit_, tm.tol_, iters_used);
            s = RootFinder::solve(res_s, res_eq_.s0, 0.0, 1.0,
                                  res_eq_.tm.maxit_, res_eq_.tm.tol_, iters_used);
            double x[2];
            x[0] = s;
            x[1] = c;
            double res = res_eq_.computeResidualC(x);
#ifdef EXTRA_DEBUG_OUTPUT
            std::cout << "c = " << c << "    s = " << s << "    c-residual = " << res << std::endl;
#endif
            return res;
        }

        double lastSaturation() const
        {
            return s;
        }
    };


    // ResidualEquation gathers parameters to construct the residual, computes its
    // value and the values of its derivatives.

    TransportSolverTwophaseCompressiblePolymer::ResidualEquation::ResidualEquation(TransportSolverTwophaseCompressiblePolymer& tmodel, int cell_index)
        : tm(tmodel)
    {
        gradient_method = Analytic;
        cell    = cell_index;
        const int np = tm.props_.numPhases();
        s0      = tm.saturation_[cell];
        c0      = tm.concentration_[cell];
        cmax0   = tm.cmax_[cell];
        double src_flux  = -tm.source_[cell];
        bool src_is_inflow = src_flux < 0.0;
        B_cell0 = 1.0/tm.A0_[np*np*cell + 0];
        B_cell = 1.0/tm.A_[np*np*cell + 0];
        influx  =  src_is_inflow ? B_cell*src_flux : 0.0;
        outflux = !src_is_inflow ? src_flux : 0.0;
        porevolume0 = tm.porevolume0_[cell];
        porevolume  = tm.porevolume_[cell];
        const double vol_cell = tm.grid_.cell_volumes[cell];
        porosity0 = porevolume0/vol_cell;
        porosity  = porevolume/vol_cell;
        dtpv  = tm.dt_/porevolume;
        dps = tm.polyprops_.deadPoreVol();
        rhor = tm.polyprops_.rockDensity();
        tm.polyprops_.adsorption(c0, cmax0, ads0);
        double mc;
        tm.computeMc(tm.polymer_inflow_c_[cell_index], mc);
        influx_polymer = src_is_inflow ? src_flux*mc : 0.0;
        for (int i = tm.grid_.cell_facepos[cell]; i < tm.grid_.cell_facepos[cell+1]; ++i) {
            int f = tm.grid_.cell_faces[i];
            double flux;
            int other;
            // Compute cell flux
            if (cell == tm.grid_.face_cells[2*f]) {
                flux  = tm.darcyflux_[f];
                other = tm.grid_.face_cells[2*f+1];
            } else {
                flux  =-tm.darcyflux_[f];
                other = tm.grid_.face_cells[2*f];
            }
            // Add flux to influx or outflux, if interior.
            if (other != -1) {
                if (flux < 0.0) {
                    const double b_face =tm.A_[np*np*other+ 0];
                    influx  += B_cell*b_face*flux*tm.fractionalflow_[other];
                    influx_polymer += flux*tm.fractionalflow_[other]*tm.mc_[other];
                } else {
                    outflux += flux; // Because B_cell*b_face = 1 for outflow faces
                }
            }
        }
    }


    void TransportSolverTwophaseCompressiblePolymer::ResidualEquation::computeResidual(const double* x, double* res) const
    {
        double dres_s_dsdc[2];
        double dres_c_dsdc[2];
        double mc;
        double ff;
        computeResAndJacobi(x, true, true, false, false, res, dres_s_dsdc, dres_c_dsdc, mc, ff);
    }

    void TransportSolverTwophaseCompressiblePolymer::ResidualEquation::computeResidual(const double* x, double* res, double& mc, double& ff) const
    {
        double dres_s_dsdc[2];
        double dres_c_dsdc[2];
        computeResAndJacobi(x, true, true, false, false, res, dres_s_dsdc, dres_c_dsdc, mc, ff);
    }


    double TransportSolverTwophaseCompressiblePolymer::ResidualEquation::computeResidualS(const double* x) const
    {
        double res[2];
        double dres_s_dsdc[2];
        double dres_c_dsdc[2];
        double mc;
        double ff;
        computeResAndJacobi(x, true, false, false, false, res, dres_s_dsdc, dres_c_dsdc, mc, ff);
        return res[0];
    }

    double TransportSolverTwophaseCompressiblePolymer::ResidualEquation::computeResidualC(const double* x) const
    {
        double res[2];
        double dres_s_dsdc[2];
        double dres_c_dsdc[2];
        double mc;
        double ff;
        computeResAndJacobi(x, false, true, false, false, res, dres_s_dsdc, dres_c_dsdc, mc, ff);
        return res[1];
    }

    void TransportSolverTwophaseCompressiblePolymer::ResidualEquation::computeGradientResS(const double* x, double* res, double* gradient) const
    // If gradient_method == FinDif, use finite difference
    // If gradient_method == Analytic, use analytic expresions
    {
        double dres_c_dsdc[2];
        double mc;
        double ff;
        computeResAndJacobi(x, true, true, true, false, res, gradient, dres_c_dsdc, mc, ff);
    }

    void TransportSolverTwophaseCompressiblePolymer::ResidualEquation::computeGradientResC(const double* x, double* res, double* gradient) const
    // If gradient_method == FinDif, use finite difference
    // If gradient_method == Analytic, use analytic expresions
    {
        double dres_s_dsdc[2];
        double mc;
        double ff;
        computeResAndJacobi(x, true, true, false, true, res, dres_s_dsdc, gradient, mc, ff);
    }

    // Compute the Jacobian of the residual equations.
    void TransportSolverTwophaseCompressiblePolymer::ResidualEquation::computeJacobiRes(const double* x, double* dres_s_dsdc, double* dres_c_dsdc) const
    {
        double res[2];
        double mc;
        double ff;
        computeResAndJacobi(x, false, false, true, true, res, dres_s_dsdc, dres_c_dsdc, mc, ff);
    }

    void TransportSolverTwophaseCompressiblePolymer::ResidualEquation::computeResAndJacobi(const double* x, const bool if_res_s, const bool if_res_c,
                                                                                  const bool if_dres_s_dsdc, const bool if_dres_c_dsdc,
                                                                                  double* res, double* dres_s_dsdc,
                                                                                  double* dres_c_dsdc, double& mc, double& ff) const
    {
        if ((if_dres_s_dsdc || if_dres_c_dsdc) && gradient_method == Analytic) {
            double s = x[0];
            double c = x[1];
            double  dff_dsdc[2];
            double mc_dc;
            double ads_dc;
            double ads;
            tm.fracFlowWithDer(s, c, cmax0, cell, ff, dff_dsdc);
            if (if_dres_c_dsdc) {
                tm.polyprops_.adsorptionWithDer(c, cmax0, ads, ads_dc);
                tm.computeMcWithDer(c, mc, mc_dc);
            } else {
                tm.polyprops_.adsorption(c, cmax0, ads);
                tm.computeMc(c, mc);
            }
            if (if_res_s) {
                res[0] = s - B_cell/B_cell0*porosity0/porosity*s0 + dtpv*(outflux*ff + influx);
#if PROFILING
                tm.res_counts.push_back(Newton_Iter(true, cell, x[0], x[1]));
#endif
            }
            if (if_res_c) {
                // Not clear if the rock compressibility should be
                // considered as a constant in the adsorption term.
                res[1] = (1 - dps)*s*c - (1 - dps)*B_cell/B_cell0*porosity0/porosity*s0*c0
                    + rhor*B_cell/porosity*((1.0 - porosity)*ads - (1.0 - porosity0)*ads0)
                    + dtpv*(outflux*ff*mc + influx_polymer);
#if PROFILING
                tm.res_counts.push_back(Newton_Iter(false, cell, x[0], x[1]));
#endif
            }
            if (if_dres_s_dsdc) {
                dres_s_dsdc[0] = 1 + dtpv*outflux*dff_dsdc[0];
                dres_s_dsdc[1] = dtpv*outflux*dff_dsdc[1];
            }
            if (if_dres_c_dsdc) {
                dres_c_dsdc[0] = (1.0 - dps)*c + dtpv*outflux*dff_dsdc[0]*mc;
                dres_c_dsdc[1] = (1 - dps)*s + rhor*B_cell/porosity*(1.0 - porosity)*ads_dc
                    + dtpv*outflux*(dff_dsdc[1]*mc + ff*mc_dc);
            }

        } else if (if_res_c || if_res_s) {
            double s = x[0];
            double c = x[1];
            tm.fracFlow(s, c, cmax0, cell, ff);
            if (if_res_s) {
                res[0] = s - B_cell/B_cell0*porosity0/porosity*s0 + dtpv*(outflux*ff + influx);
#if PROFILING
                tm.res_counts.push_back(Newton_Iter(true, cell, x[0], x[1]));
#endif
            }
            if (if_res_c) {
                tm.computeMc(c, mc);
                double ads;
                tm.polyprops_.adsorption(c, cmax0, ads);
                res[1] = (1 - dps)*s*c - (1 - dps)*B_cell/B_cell0*porosity0/porosity*s0*c0
                    + rhor*B_cell/porosity*((1.0 - porosity)*ads - (1.0 - porosity0)*ads0)
                    + dtpv*(outflux*ff*mc + influx_polymer);
#if PROFILING
                tm.res_counts.push_back(Newton_Iter(false, cell, x[0], x[1]));
#endif
            }
        }

        if ((if_dres_c_dsdc || if_dres_s_dsdc) && gradient_method == FinDif) {
            double epsi = 1e-8;
            double res_epsi[2];
            double res_0[2];
            double x_epsi[2];
            computeResidual(x, res_0);
            if (if_dres_s_dsdc) {
                x_epsi[0] = x[0] + epsi;
                x_epsi[1] = x[1];
                computeResidual(x_epsi, res_epsi);
                dres_s_dsdc[0] =  (res_epsi[0] - res_0[0])/epsi;
                x_epsi[0] = x[0];
                x_epsi[1] = x[1] + epsi;
                computeResidual(x_epsi, res_epsi);
                dres_s_dsdc[1] = (res_epsi[0] - res_0[0])/epsi;
            }
            if (if_dres_c_dsdc) {
                x_epsi[0] = x[0] + epsi;
                x_epsi[1] = x[1];
                computeResidual(x_epsi, res_epsi);
                dres_c_dsdc[0] =  (res_epsi[1] - res_0[1])/epsi;
                x_epsi[0] = x[0];
                x_epsi[1] = x[1] + epsi;
                computeResidual(x_epsi, res_epsi);
                dres_c_dsdc[1] = (res_epsi[1] - res_0[1])/epsi;
            }
        }
    }

    // Compute the "s" residual along the curve "curve" for a given residual equation "res_eq".
    // The operator() is sent to a root solver.
    class TransportSolverTwophaseCompressiblePolymer::ResSOnCurve
    {
    public:
        ResSOnCurve(const TransportSolverTwophaseCompressiblePolymer::ResidualEquation& res_eq);
        double operator()(const double t) const;
        CurveInSCPlane curve;
    private:
        const TransportSolverTwophaseCompressiblePolymer::ResidualEquation& res_eq_;
    };

    // Compute the "c" residual along the curve "curve" for a given residual equation "res_eq".
    // The operator() is sent to a root solver.
    class TransportSolverTwophaseCompressiblePolymer::ResCOnCurve
    {
    public:
        ResCOnCurve(const TransportSolverTwophaseCompressiblePolymer::ResidualEquation& res_eq);
        double operator()(const double t) const;
        CurveInSCPlane curve;
    private:
        const TransportSolverTwophaseCompressiblePolymer::ResidualEquation& res_eq_;
    };

    TransportSolverTwophaseCompressiblePolymer::ResSOnCurve::ResSOnCurve(const TransportSolverTwophaseCompressiblePolymer::ResidualEquation& res_eq)
        : res_eq_(res_eq)
    {
    }

    double TransportSolverTwophaseCompressiblePolymer::ResSOnCurve::operator()(const double t) const
    {
        double x_of_t[2];
        double x_c[2];
        curve.computeXOfT(x_of_t, t);
        res_eq_.tm.scToc(x_of_t, x_c);
        return res_eq_.computeResidualS(x_c);
    }

    TransportSolverTwophaseCompressiblePolymer::ResCOnCurve::ResCOnCurve(const TransportSolverTwophaseCompressiblePolymer::ResidualEquation& res_eq)
        : res_eq_(res_eq)
    {
    }

    double TransportSolverTwophaseCompressiblePolymer::ResCOnCurve::operator()(const double t) const
    {
        double x_of_t[2];
        double x_c[2];
        curve.computeXOfT(x_of_t, t);
        res_eq_.tm.scToc(x_of_t, x_c);
        return res_eq_.computeResidualC(x_c);
    }



    void TransportSolverTwophaseCompressiblePolymer::solveSingleCell(const int cell)
    {
        switch (method_) {
        case Bracketing:
            solveSingleCellBracketing(cell);
            break;
        case Newton:
            solveSingleCellNewton(cell, true);
            break;
        case NewtonC:
            solveSingleCellNewton(cell, false);
            break;
        case Gradient:
            solveSingleCellGradient(cell);
            break;
        default:
            OPM_THROW(std::runtime_error, "Unknown method " << method_);
        }
    }


    void TransportSolverTwophaseCompressiblePolymer::solveSingleCellBracketing(int cell)
    {

        ResidualEquation res_eq(*this, cell);
        ResidualC res(res_eq);
        const double a = 0.0;
        const double b = polyprops_.cMax()*adhoc_safety_; // Add 10% to account for possible non-monotonicity of hyperbolic system.
        int iters_used;

        // Check if current state is an acceptable solution.
        double res_sc[2];
        double mc, ff;
        res.computeBothResiduals(saturation_[cell], concentration_[cell], res_sc[0], res_sc[1], mc, ff);
        if (norm(res_sc) < tol_) {
            fractionalflow_[cell] = ff;
            mc_[cell] = mc;
            return;
        }

        concentration_[cell] = RootFinder::solve(res, a, b, maxit_, tol_, iters_used);
        cmax_[cell] = std::max(cmax_[cell], concentration_[cell]);
        saturation_[cell] = res.lastSaturation();
        fracFlow(saturation_[cell], concentration_[cell], cmax_[cell], cell,
                 fractionalflow_[cell]);
        computeMc(concentration_[cell], mc_[cell]);
    }



    // Newton method, where we first try a Newton step. Then, if it does not work well, we look for
    // the zero of either the residual in s or the residual in c along a specified piecewise linear
    // curve. In these cases, we can use a robust 1d solver.
    void TransportSolverTwophaseCompressiblePolymer::solveSingleCellGradient(int cell)
    {
        int iters_used_falsi = 0;
        const int max_iters_split = maxit_;
        int iters_used_split = 0;

        // Check if current state is an acceptable solution.
        ResidualEquation res_eq(*this, cell);
        double x[2] = {saturation_[cell], saturation_[cell]*concentration_[cell]};
        double res[2];
        double mc;
        double ff;
        double x_c[2];
        scToc(x, x_c);
        res_eq.computeResidual(x_c, res, mc, ff);
        if (norm(res) <= tol_) {
            cmax_[cell] = std::max(cmax_[cell], concentration_[cell]);
            fractionalflow_[cell] = ff;
            mc_[cell] = mc;
            return;
        }

        double x_min[2] = { 0.0, 0.0 };
        double x_max[2] = { 1.0, polyprops_.cMax()*adhoc_safety_ };
        double x_min_res_s[2] = { x_min[0], x_min[1] };
        double x_max_res_s[2] = { x_max[0], x_min[0] };
        double x_min_res_sc[2] = { x_min[0], x_min[1] };
        double x_max_res_sc[2] = { x_max[0], x_max[1] };
        double t;
        double t_max;
        double t_out;
        double direction[2];
        double end_point[2];
        double gradient[2];
        ResSOnCurve res_s_on_curve(res_eq);
        ResCOnCurve res_c_on_curve(res_eq);
        bool if_res_s;

        while ((norm(res) > tol_) && (iters_used_split < max_iters_split)) {
            if (std::abs(res[0]) < std::abs(res[1])) {
                if (res[0] < -tol_) {
                    direction[0] = x_max_res_s[0] - x[0];
                    direction[1] = x_max_res_s[1] - x[1];
                    if_res_s = true;
                } else if (res[0] > tol_) {
                    direction[0] = x_min_res_s[0] - x[0];
                    direction[1] = x_min_res_s[1] - x[1];
                    if_res_s = true;
                } else {
                    scToc(x, x_c);
                    res_eq.computeGradientResS(x_c, res, gradient);
                    // dResS/d(s_) = dResS/ds - c/s*dResS/ds
                    // dResS/d(sc_) = -1/s*dResS/dc
                    if (x[0] > 1e-2*tol_) {
                        // With s,c variables, we should have
                        // direction[0] = -gradient[1];
                        // direction[1] = gradient[0];
                        // With s, sc variables, we get:
                        scToc(x, x_c);
                        direction[0] = 1.0/x[0]*gradient[1];
                        direction[1] = gradient[0] - x_c[1]/x[0]*gradient[1];
                    } else {
                        // acceptable approximation for nonlinear relative permeability.
                        direction[0] = 0.0;
                        direction[1] = gradient[0];
                    }
                    if_res_s = false;
                }
            } else {
                if (res[1] < -tol_) {
                    direction[0] = x_max_res_sc[0] - x[0];
                    direction[1] = x_max_res_sc[1] - x[1];
                    if_res_s = false;
                } else if (res[1] > tol_) {
                    direction[0] = x_min_res_sc[0] - x[0];
                    direction[1] = x_min_res_sc[1] - x[1];
                    if_res_s = false;
                } else {
                    res_eq.computeGradientResC(x, res, gradient);
                    // dResC/d(s_) = dResC/ds - c/s*dResC/ds
                    // dResC/d(sc_) = -1/s*dResC/dc
                    if (x[0] > 1e-2*tol_) {
                        // With s,c variables, we should have
                        // direction[0] = -gradient[1];
                        // direction[1] = gradient[0];
                        // With s, sc variables, we get:
                        scToc(x, x_c);
                        direction[0] = 1.0/x[0]*gradient[1];
                        direction[1] = gradient[0] - x_c[1]/x[0]*gradient[1];
                    } else {
                        // We take 1.0/s*gradient[1]: wrong for linear permeability,
                        // acceptable for nonlinear relative permeability.
                        direction[0] = 1.0 - res_eq.dps;
                        direction[1] = gradient[0];
                    }
                    if_res_s = true;
                }
            }
            if (if_res_s) {
                if (res[0] < 0) {
                    end_point[0] = x_max_res_s[0];
                    end_point[1] = x_max_res_s[1];
                    res_s_on_curve.curve.setup(x, direction, end_point, x_min, x_max, tol_, t_max, t_out);
                    if (res_s_on_curve(t_out) >= 0) {
                        t_max = t_out;
                    }
                } else {
                    end_point[0] = x_min_res_s[0];
                    end_point[1] = x_min_res_s[1];
                    res_s_on_curve.curve.setup(x, direction, end_point, x_min, x_max, tol_, t_max, t_out);
                    if (res_s_on_curve(t_out) <= 0) {
                        t_max = t_out;
                    }
                }
                // Note: In some experiments modifiedRegularFalsi does not yield a result under the given tolerance.
                t = RootFinder::solve(res_s_on_curve, 0., t_max, maxit_, tol_, iters_used_falsi);
                res_s_on_curve.curve.computeXOfT(x, t);
            } else {
                if (res[1] < 0) {
                    end_point[0] = x_max_res_sc[0];
                    end_point[1] = x_max_res_sc[1];
                    res_c_on_curve.curve.setup(x, direction, end_point, x_min, x_max, tol_, t_max, t_out);
                    if (res_c_on_curve(t_out) >= 0) {
                        t_max = t_out;
                    }
                } else {
                    end_point[0] = x_min_res_sc[0];
                    end_point[1] = x_min_res_sc[1];
                    res_c_on_curve.curve.setup(x, direction, end_point, x_min, x_max, tol_, t_max, t_out);
                    if (res_c_on_curve(t_out) <= 0) {
                        t_max = t_out;
                    }
                }
                t = RootFinder::solve(res_c_on_curve, 0., t_max, maxit_, tol_, iters_used_falsi);
                res_c_on_curve.curve.computeXOfT(x, t);

            }
            scToc(x, x_c);
            res_eq.computeResidual(x_c, res, mc, ff);
            iters_used_split += 1;
        }



        if ((iters_used_split >=  max_iters_split) && (norm(res) > tol_)) {
            OPM_MESSAGE("Newton for single cell did not work in cell number " << cell);
            solveSingleCellBracketing(cell);
        } else {
            scToc(x, x_c);
            concentration_[cell] = x_c[1];
            cmax_[cell] = std::max(cmax_[cell], concentration_[cell]);
            saturation_[cell] = x[0];
            fractionalflow_[cell] = ff;
            mc_[cell] = mc;
        }
    }

    void TransportSolverTwophaseCompressiblePolymer::solveSingleCellNewton(int cell, bool use_sc,
                                                                  bool use_explicit_step)
    {
        const int max_iters_split = maxit_;
        int iters_used_split = 0;

        // Check if current state is an acceptable solution.
        ResidualEquation res_eq(*this, cell);
        double x[2] = {saturation_[cell], concentration_[cell]};
        double res[2];
        double mc;
        double ff;
        res_eq.computeResidual(x, res, mc, ff);
        if (norm(res) <= tol_) {
            cmax_[cell] = std::max(cmax_[cell], concentration_[cell]);
            fractionalflow_[cell] = ff;
            mc_[cell] = mc;
            return;
        }

        if (use_explicit_step) {
            // x is updated to an explicit step.
            x[0] = saturation_[cell]-res[0];
            if ((x[0]>1) || (x[0]<0)) {
                // If we are outside the allowed domain for s, we
                // reset s to 0.5, which should not far from the
                // inflexion point of the residual, that is, the point
                // where Newton's method performs best.
                x[0] = 0.5;
                x[1] = x[1];
            }
            if (x[0]>0) {
                x[1] =  concentration_[cell]*saturation_[cell]-res[1];
                x[1] = x[1]/x[0];
                if(x[1]> polyprops_.cMax()){
                    x[1]= polyprops_.cMax()/2.0;
                }
                if(x[1]<0){
                    x[1]=0;
                }
            } else {
                x[1]=0;
            }
            res_eq.computeResidual(x, res, mc, ff);
        }

        const double x_min[2] = { 0.0, 0.0 };
        const double x_max[2] = { 1.0, polyprops_.cMax()*adhoc_safety_ };
        bool successfull_newton_step = true;

        // initialize x_new to avoid warning
        double x_new[2] = {0.0, 0.0};
        double res_new[2];

        if (use_sc) {
            // We switch to variables x[0] = s, x[1] = sc.
            x[1] = x[0]*x[1];
        }

        // x_c will contain the s-c  variable when use_sc = true
        double x_c[2];

        // Variables to store the Jacobian.
        double dFx_dx;
        double dFx_dy;
        double dFy_dx;
        double dFy_dy;

        while ((norm(res) > tol_) &&
               (iters_used_split < max_iters_split)  &&
               successfull_newton_step) {
            double dres_s_dsdc[2];
            double dres_c_dsdc[2];
            if (use_sc) {
                // Convert from (s, c) to (s, sc) variables.
                scToc(x, x_c);
                double x_c_app[2];
                // The computation of the Jacobi fails for s=0 (we have an undetermined fraction 0/0).
                // When s is close to zero we replace x_c with x_c_app as defined now.
                x_c_app[1] = x_c[1];
                if (x_c[0] < 1e-2*tol_) {
                    x_c_app[0] = 1e-2*tol_;
                } else {
                    x_c_app[0] = x_c[0];
                }
                res_eq.computeJacobiRes(x_c_app, dres_s_dsdc, dres_c_dsdc);
                dFx_dx = (dres_s_dsdc[0]-x_c_app[1]*dres_s_dsdc[1]);
                dFx_dy = (dres_s_dsdc[1]/x_c_app[0]);
                dFy_dx = (dres_c_dsdc[0]-x_c_app[1]*dres_c_dsdc[1]);
                dFy_dy = (dres_c_dsdc[1]/x_c_app[0]);
            } else {
                res_eq.computeJacobiRes(x, dres_s_dsdc, dres_c_dsdc);
                dFx_dx= dres_s_dsdc[0];
                dFx_dy= dres_s_dsdc[1];
                dFy_dx= dres_c_dsdc[0];
                dFy_dy= dres_c_dsdc[1];
            }
            double det = dFx_dx*dFy_dy - dFy_dx*dFx_dy;
            double alpha = 1.0;
            int max_lin_it = 100;
            int lin_it = 0;
            res_new[0] = res[0]*2;
            res_new[1] = res[1]*2;
            while((norm(res_new)>norm(res)) && (lin_it<max_lin_it)) {
                x_new[0] = x[0] - alpha*(res[0]*dFy_dy - res[1]*dFx_dy)/det;
                x_new[1] = x[1] - alpha*(res[1]*dFx_dx - res[0]*dFy_dx)/det;
                if (use_sc) {
                    scToc(x_new, x_c);
                    check_interval(x_min, x_max, x_c);
                    res_eq.computeResidual(x_c, res_new, mc, ff);
                } else {
                    check_interval(x_min, x_max, x);
                    res_eq.computeResidual(x, res_new, mc, ff);
                }
                alpha = alpha/2.0;
                lin_it = lin_it + 1;
            }
            if (lin_it>=max_lin_it) {
                successfull_newton_step = false;
            } else  {
                if (use_sc) {
                    scToc(x_new, x);
                } else {
                    x[0] = x_new[0];
                    x[1] = x_new[1];
                }
                res[0] = res_new[0];
                res[1] = res_new[1];
                iters_used_split += 1;
                successfull_newton_step = true;;
            }
        }

        if ((iters_used_split >=  max_iters_split) && (norm(res) > tol_)) {
            OPM_MESSAGE("Newton for single cell did not work in cell number " << cell);
            solveSingleCellBracketing(cell);
        } else {
            concentration_[cell] = x[1];
            cmax_[cell] = std::max(cmax_[cell], concentration_[cell]);
            saturation_[cell] = x[0];
            fractionalflow_[cell] = ff;
            mc_[cell] = mc;
        }
    }


    void TransportSolverTwophaseCompressiblePolymer::solveMultiCell(const int num_cells, const int* cells)
    {
        double max_s_change = 0.0;
        double max_c_change = 0.0;
        int num_iters = 0;
        // Must store state variables before we start.
        std::vector<double> s0(num_cells);
        std::vector<double> c0(num_cells);
        std::vector<double> cmax0(num_cells);
        // Must set initial fractional flows etc. before we start.
        for (int i = 0; i < num_cells; ++i) {
            const int cell = cells[i];
            fracFlow(saturation_[cell], concentration_[cell], cmax_[cell],
                     cell, fractionalflow_[cell]);
            computeMc(concentration_[cell], mc_[cell]);
            s0[i] = saturation_[cell];
            c0[i] = concentration_[cell];
            cmax0[i] = cmax_[i];
        }
        do {
            // int max_s_change_cell = -1;
            // int max_c_change_cell = -1;
            max_s_change = 0.0;
            max_c_change = 0.0;
            for (int i = 0; i < num_cells; ++i) {
                const int cell = cells[i];
                const double old_s = saturation_[cell];
                const double old_c = concentration_[cell];
                saturation_[cell] = s0[i];
                concentration_[cell] = c0[i];
                cmax_[cell] = cmax0[i];
                solveSingleCell(cell);
                // std::cout << "cell = " << cell << "    delta s = " << saturation_[cell] - old_s << std::endl;
                // if (max_s_change < std::fabs(saturation_[cell] - old_s)) {
                //     max_s_change_cell = cell;
                // }
                // if (max_c_change < std::fabs(concentration_[cell] - old_c)) {
                //     max_c_change_cell = cell;
                // }
                max_s_change = std::max(max_s_change, std::fabs(saturation_[cell] - old_s));
                max_c_change = std::max(max_c_change, std::fabs(concentration_[cell] - old_c));
            }
            // std::cout << "Iter = " << num_iters << "    max_s_change = " << max_s_change
            //        << "    in cell " << max_change_cell << std::endl;
        } while (((max_s_change > tol_) || (max_c_change > tol_)) && ++num_iters < maxit_);
        if (max_s_change > tol_) {
            OPM_THROW(std::runtime_error, "In solveMultiCell(), we did not converge after "
                  << num_iters << " iterations. Delta s = " << max_s_change);
        }
        if (max_c_change > tol_) {
            OPM_THROW(std::runtime_error, "In solveMultiCell(), we did not converge after "
                  << num_iters << " iterations. Delta c = " << max_c_change);
        }
        // std::cout << "Solved " << num_cells << " cell multicell problem in "
        //           << num_iters << " iterations." << std::endl;
    }

    void TransportSolverTwophaseCompressiblePolymer::fracFlow(double s, double c, double cmax,
                                                     int cell, double& ff) const
    {
        double dummy[2];
        fracFlowBoth(s, c, cmax, cell, ff,  dummy, false);
    }

    void TransportSolverTwophaseCompressiblePolymer::fracFlowWithDer(double s, double c, double cmax,
                                                            int cell, double& ff,
                                                            double* dff_dsdc) const
    {
        fracFlowBoth(s, c, cmax, cell, ff, dff_dsdc, true);
    }

    void TransportSolverTwophaseCompressiblePolymer::fracFlowBoth(double s, double c, double cmax, int cell,
                                                         double& ff, double* dff_dsdc,
                                                         bool if_with_der) const
    {
        double relperm[2];
        double drelperm_ds[4];
        double sat[2] = {s, 1 - s};
        if (if_with_der) {
            props_.relperm(1, sat, &cell, relperm, drelperm_ds);
        } else {
            props_.relperm(1, sat, &cell, relperm, 0);
        }
        double mob[2];
        double dmob_ds[4];
        double dmob_dc[2];
        double dmobwat_dc;
        const int np = props_.numPhases();
        polyprops_.effectiveMobilitiesBoth(c, cmax, &visc_[np*cell], relperm, drelperm_ds,
                                           mob, dmob_ds, dmobwat_dc, if_with_der);

        ff = mob[0]/(mob[0] + mob[1]);
        if (if_with_der) {
            dmob_dc[0] = dmobwat_dc;
            dmob_dc[1] = 0.;
            //dff_dsdc[0] = (dmob_ds[0]*mob[1] + dmob_ds[3]*mob[0])/((mob[0] + mob[1])*(mob[0] + mob[1])); // derivative with respect to s
            // at the moment the dmob_ds only have diagonal elements since the saturation is derivated out in effectiveMobilitiesBoth
            dff_dsdc[0] = ((dmob_ds[0]-dmob_ds[2])*mob[1] - (dmob_ds[1]-dmob_ds[3])*mob[0])/((mob[0] + mob[1])*(mob[0] + mob[1])); // derivative with respect to s
            dff_dsdc[1] = (dmob_dc[0]*mob[1] - dmob_dc[1]*mob[0])/((mob[0] + mob[1])*(mob[0] + mob[1])); // derivative with respect to c
        }
    }

    void TransportSolverTwophaseCompressiblePolymer::computeMc(double c, double& mc) const
    {
        polyprops_.computeMc(c, mc);
    }

    void TransportSolverTwophaseCompressiblePolymer::computeMcWithDer(double c, double& mc,
                                                             double &dmc_dc) const
    {
        polyprops_.computeMcWithDer(c, mc, dmc_dc);
    }



    TransportSolverTwophaseCompressiblePolymer::ResidualSGrav::ResidualSGrav(const ResidualCGrav& res_c_eq,
                                                                    const double c_init)
        : res_c_eq_(res_c_eq),
          c(c_init)
    {
    }

    double TransportSolverTwophaseCompressiblePolymer::ResidualSGrav::operator()(double s) const
    {
        return res_c_eq_.computeGravResidualS(s, c);
    }




    // Residual for concentration equation for gravity segregation
    //
    // res_c = s*(1 - dps)*c - s0*( - dps)*c0
    //     +  dtpv*sum_{j adj i}( mc * gravmod_ij * gf_ij ).
    //  \TODO doc me
    // where ...
    // Influxes are negative, outfluxes positive.

    TransportSolverTwophaseCompressiblePolymer::ResidualCGrav::ResidualCGrav(const TransportSolverTwophaseCompressiblePolymer& tmodel,
                                                                    const std::vector<int>& cells,
                                                                    const int pos,
                                                                    const double* gravflux) // Always oriented towards next in column. Size = colsize - 1.
        : tm(tmodel),
          cell(cells[pos]),
          s0(tm.saturation_[cell]),
          c0(tm.concentration_[cell]),
          cmax0(tm.cmax0_[cell]),
          porevolume(tm.porevolume_[cell]),
          porosity(porevolume/tm.grid_.cell_volumes[cell]),
          dtpv(tm.dt_/porevolume),
          dps(tm.polyprops_.deadPoreVol()),
          rhor(tm.polyprops_.rockDensity())

    {
        last_s = s0;
        nbcell[0] = -1;
        gf[0] = 0.0;
        if (pos > 0) {
            nbcell[0] = cells[pos - 1];
            gf[0] = -gravflux[pos - 1];
        }
        nbcell[1] = -1;
        gf[1] = 0.0;
        if (pos < int(cells.size() - 1)) {
            nbcell[1] = cells[pos + 1];
            gf[1] = gravflux[pos];
        }

        tm.polyprops_.adsorption(c0, cmax0, c_ads0);
    }

    double TransportSolverTwophaseCompressiblePolymer::ResidualCGrav::operator()(double c) const
    {

        ResidualSGrav res_s(*this);
        res_s.c = c;
        int iters_used;
        last_s =  RootFinder::solve(res_s, last_s, 0.0, 1.0,
                                    tm.maxit_, tm.tol_,
                                    iters_used);
        return computeGravResidualC(last_s, c);

    }

    double TransportSolverTwophaseCompressiblePolymer::ResidualCGrav::computeGravResidualS(double s, double c) const
    {

        double mobcell[2];
        tm.mobility(s, c, cell, mobcell);

        double res = s - s0;

        for (int nb = 0; nb < 2; ++nb) {
            if (nbcell[nb] != -1) {
                double m[2];
                if (gf[nb] < 0.0) {
                    m[0] = mobcell[0];
                    m[1] = tm.mob_[2*nbcell[nb] + 1];
                } else {
                    m[0] = tm.mob_[2*nbcell[nb]];
                    m[1] = mobcell[1];
                }
                if (m[0] + m[1] > 0.0) {
                    res += -dtpv*gf[nb]*m[0]*m[1]/(m[0] + m[1]);
                }
            }
        }
        return res;
    }

    double TransportSolverTwophaseCompressiblePolymer::ResidualCGrav::computeGravResidualC(double s, double c) const
    {

        double mobcell[2];
        tm.mobility(s, c, cell,  mobcell);
        double c_ads;
        tm.polyprops_.adsorption(c, cmax0, c_ads);

        double res = (1 - dps)*s*c - (1 - dps)*s0*c0
            + rhor*((1.0 - porosity)/porosity)*(c_ads - c_ads0);
        for (int nb = 0; nb < 2; ++nb) {
            if (nbcell[nb] != -1) {
                double m[2];
                double mc;
                if (gf[nb] < 0.0) {
                    m[0] = mobcell[0];
                    tm.computeMc(c, mc);
                    m[1] = tm.mob_[2*nbcell[nb] + 1];
                } else {
                    m[0] = tm.mob_[2*nbcell[nb]];
                    mc = tm.mc_[nbcell[nb]];
                    m[1] = mobcell[1];
                }
                if (m[0] + m[1] > 0.0) {
                    res += -dtpv*gf[nb]*mc*m[0]*m[1]/(m[0] + m[1]);
                }
            }
        }
        return res;
    }

    double TransportSolverTwophaseCompressiblePolymer::ResidualCGrav::lastSaturation() const
    {
        return last_s;
    }


    void TransportSolverTwophaseCompressiblePolymer::mobility(double s, double c, int cell, double* mob) const
    {
        double sat[2] = { s, 1.0 - s };
        double relperm[2];
        const int np = props_.numPhases();
        props_.relperm(1, sat, &cell, relperm, 0);
        polyprops_.effectiveMobilities(c, cmax0_[cell], &visc_[np*cell], relperm, mob);
    }

    void TransportSolverTwophaseCompressiblePolymer::initGravity(const double* grav)
    {
        // Set up transmissibilities.
        std::vector<double> htrans(grid_.cell_facepos[grid_.number_of_cells]);
        const int nf = grid_.number_of_faces;
        trans_.resize(nf);
        gravflux_.resize(nf);
        tpfa_htrans_compute(const_cast<UnstructuredGrid*>(&grid_), props_.permeability(), &htrans[0]);
        tpfa_trans_compute(const_cast<UnstructuredGrid*>(&grid_), &htrans[0], &trans_[0]);

        // Remember gravity vector.
        gravity_ = grav;
    }

    void TransportSolverTwophaseCompressiblePolymer::initGravityDynamic()
    {
        // Set up gravflux_ = T_ij g [   (b_w,i rho_w,S - b_o,i rho_o,S) (z_i - z_f)
        //                             + (b_w,j rho_w,S - b_o,j rho_o,S) (z_f - z_j) ]
        // But b_w,i * rho_w,S = rho_w,i, which we compute with a call to props_.density().
        // We assume that we already have stored T_ij in trans_.
        // We also assume that the A_ matrices are updated from an earlier call to solve().
        const int nc = grid_.number_of_cells;
        const int nf = grid_.number_of_faces;
        const int np = props_.numPhases();
        assert(np == 2);
        const int dim = grid_.dimensions;
        density_.resize(nc*np);
        props_.density(grid_.number_of_cells, &A_[0], grid_.global_cell, &density_[0]);
        std::fill(gravflux_.begin(), gravflux_.end(), 0.0);
        for (int f = 0; f < nf; ++f) {
            const int* c = &grid_.face_cells[2*f];
            const double signs[2] = { 1.0, -1.0 };
            if (c[0] != -1 && c[1] != -1) {
                for (int ci = 0; ci < 2; ++ci) {
                    double gdz = 0.0;
                    for (int d = 0; d < dim; ++d) {
                        gdz += gravity_[d]*(grid_.cell_centroids[dim*c[ci] + d] - grid_.face_centroids[dim*f + d]);
                    }
                    gravflux_[f] += signs[ci]*trans_[f]*gdz*(density_[2*c[ci]] - density_[2*c[ci] + 1]);
                }
            }
        }
    }


    void TransportSolverTwophaseCompressiblePolymer::solveSingleCellGravity(const std::vector<int>& cells,
                                                                   const int pos,
                                                                   const double* gravflux)
    {
        const int cell = cells[pos];
        ResidualCGrav res_c(*this, cells, pos, gravflux);

        // Check if current state is an acceptable solution.
        double res_sc[2];
        res_sc[0]=res_c.computeGravResidualS(saturation_[cell], concentration_[cell]);
        res_sc[1]=res_c.computeGravResidualC(saturation_[cell], concentration_[cell]);

        if (norm(res_sc) < tol_) {
            return;
        }

        const double a = 0.0;
        const double b = polyprops_.cMax()*adhoc_safety_; // Add 10% to account for possible non-monotonicity of hyperbolic system.
        int iters_used;
        concentration_[cell] = RootFinder::solve(res_c, concentration_[cell],
                                                 a, b, maxit_, tol_, iters_used);
        saturation_[cell] = res_c.lastSaturation();
        cmax_[cell] = std::max(cmax0_[cell], concentration_[cell]);
        computeMc(concentration_[cell], mc_[cell]);
        mobility(saturation_[cell], concentration_[cell], cell, &mob_[2*cell]);
    }

    int TransportSolverTwophaseCompressiblePolymer::solveGravityColumn(const std::vector<int>& cells)
    {
        // Set up column gravflux.
        const int nc = cells.size();
        std::vector<double> col_gravflux(nc - 1);
        for (int ci = 0; ci < nc - 1; ++ci) {
            const int cell = cells[ci];
            const int next_cell = cells[ci + 1];
            for (int j = grid_.cell_facepos[cell]; j < grid_.cell_facepos[cell+1]; ++j) {
                const int face = grid_.cell_faces[j];
                const int c1 = grid_.face_cells[2*face + 0];
                const int c2 = grid_.face_cells[2*face + 1];
                if (c1 == next_cell || c2 == next_cell) {
                    const double gf = gravflux_[face];
                    col_gravflux[ci] = (c1 == cell) ? gf : -gf;
                }
            }
        }

        // Store initial saturation s0
        s0_.resize(nc);
        c0_.resize(nc);
        for (int ci = 0; ci < nc; ++ci) {
            s0_[ci] = saturation_[cells[ci]];
            c0_[ci] = concentration_[cells[ci]];
        }

        // Solve single cell problems, repeating if necessary.
        double max_sc_change = 0.0;
        int num_iters = 0;
        do {
            max_sc_change = 0.0;
            for (int ci = 0; ci < nc; ++ci) {
                const int ci2 = nc - ci - 1;
                double old_s[2] = { saturation_[cells[ci]],
                                    saturation_[cells[ci2]] };
                double old_c[2] = { concentration_[cells[ci]],
                                    concentration_[cells[ci2]] };
                saturation_[cells[ci]] = s0_[ci];
                concentration_[cells[ci]] = c0_[ci];
                solveSingleCellGravity(cells, ci, &col_gravflux[0]);
                saturation_[cells[ci2]] = s0_[ci2];
                concentration_[cells[ci2]] = c0_[ci2];
                solveSingleCellGravity(cells, ci2, &col_gravflux[0]);
                max_sc_change = std::max(max_sc_change, 0.25*(std::fabs(saturation_[cells[ci]] - old_s[0]) +
                                                              std::fabs(concentration_[cells[ci]] - old_c[0]) +
                                                              std::fabs(saturation_[cells[ci2]] - old_s[1]) +
                                                              std::fabs(concentration_[cells[ci2]] - old_c[1])));
            }
            // std::cout << "Iter = " << num_iters << "    max_s_change = " << max_s_change << std::endl;
        } while (max_sc_change > tol_ && ++num_iters < maxit_);

        if (max_sc_change > tol_) {
            OPM_THROW(std::runtime_error, "In solveGravityColumn(), we did not converge after "
                  << num_iters << " iterations. Delta s = " << max_sc_change);
        }
        return num_iters + 1;
    }


    void TransportSolverTwophaseCompressiblePolymer::solveGravity(const std::vector<std::vector<int> >& columns,
                                                         const double dt,
                                                         std::vector<double>& saturation,
                                                         std::vector<double>& surfacevol,
                                                         std::vector<double>& concentration,
                                                         std::vector<double>& cmax)
    {

        // Assume that solve() has already been called, so that A_ and
        // porosity_ are current.
        initGravityDynamic();

        // initialize variables.
        dt_ = dt;
        toWaterSat(saturation, saturation_);
        concentration_ = &concentration[0];
        cmax_ = &cmax[0];
        const int nc = grid_.number_of_cells;
        cmax0_.resize(nc);
        std::copy(cmax.begin(), cmax.end(), &cmax0_[0]);

        // Initialize mobilities.
        const int np = props_.numPhases();
        mob_.resize(np*nc);

        for (int cell = 0; cell < nc; ++cell) {
            mobility(saturation_[cell], concentration_[cell], cell, &mob_[np*cell]);
        }


        // Solve on all columns.
        int num_iters = 0;
        // std::cout << "Gauss-Seidel column solver # columns: " << columns.size() << std::endl;
        for (std::vector<std::vector<int> >::size_type i = 0; i < columns.size(); i++) {
            // std::cout << "==== new column" << std::endl;
            num_iters += solveGravityColumn(columns[i]);
        }
        std::cout << "Gauss-Seidel column solver average iterations: "
                  << double(num_iters)/double(columns.size()) << std::endl;

        toBothSat(saturation_, saturation);
        // Compute surface volume as a postprocessing step from saturation and A_
        computeSurfacevol(grid_.number_of_cells, props_.numPhases(), &A_[0], &saturation[0], &surfacevol[0]);

    }

    void TransportSolverTwophaseCompressiblePolymer::scToc(const double* x, double* x_c) const {
        x_c[0] = x[0];
        if (x[0] < 1e-2*tol_) {
            x_c[1] = 0.5*polyprops_.cMax();
        } else {
            x_c[1] = x[1]/x[0];
        }
    }


} // namespace Opm


namespace
{
    bool check_interval(const double* xmin, const double* xmax, double* x) {
        bool test = false;
        if (x[0] < xmin[0]) {
            test = true;
            x[0] = xmin[0];
        } else if (x[0] > xmax[0]) {
            test = true;
            x[0] = xmax[0];
        }
        if (x[1] < xmin[1]) {
            test = true;
            x[1] = xmin[1];
        } else if (x[1] > xmax[1]) {
            test = true;
            x[1] = xmax[1];
        }
        return test;
    }


    CurveInSCPlane::CurveInSCPlane()
    {
    }

    // Setup the curve (see comment above).
    // The curve is parametrized by t in [0, t_max], t_out is equal to t when the curve hits the bounding
    // rectangle. x_out=(s_out, c_out) denotes the values of s and c at that point.
    void CurveInSCPlane::setup(const double* x, const double* direction,
                               const double* end_point, const double* x_min,
                               const double* x_max, const double tol,
                               double& t_max_out, double& t_out_out)
    {
        x_[0] = x[0];
        x_[1] = x[1];
        x_max_[0] = x_max[0];
        x_max_[1] = x_max[1];
        x_min_[0] = x_min[0];
        x_min_[1] = x_min[1];
        direction_[0] = direction[0];
        direction_[1] = direction[1];
        end_point_[0] = end_point[0];
        end_point_[1] = end_point[1];
        const double size_direction = std::abs(direction_[0]) + std::abs(direction_[1]);
        if (size_direction < tol) {
            direction_[0] = end_point_[0]-x_[0];
            direction_[1] = end_point_[1]-x_[1];
        } else if ((end_point_[0]-x_[0])*direction_[0] + (end_point_[1]-x_[1])*direction_[1] < 0) {
            direction_[0] *= -1.0;
            direction_[1] *= -1.0;
        }
        bool t0_exists = true;
        double t0 = 0; // dummy default value (so that compiler does not complain).
        if (direction_[0] > 0) {
            t0 = (x_max_[0] - x_[0])/direction_[0];
        } else if (direction_[0] < 0) {
            t0 = (x_min_[0] - x_[0])/direction_[0];
        } else {
            t0_exists = false;
        }
        bool t1_exists = true;
        double t1 = 0; // dummy default value.
        if (direction_[1] > 0) {
            t1 = (x_max_[1] - x_[1])/direction_[1];
        } else if (direction_[1] < 0) {
            t1 = (x_min_[1] - x_[1])/direction_[1];
        } else {
            t1_exists = false;
        }
        if (t0_exists) {
            if (t1_exists) {
                t_out_ = std::min(t0, t1);
            } else {
                t_out_ = t0;
            }
        } else if (t1_exists) {
            t_out_ = t1;
        } else {
            OPM_THROW(std::runtime_error, "Direction illegal: is a zero vector.");
        }
        x_out_[0] = x_[0] + t_out_*direction_[0];
        x_out_[1] = x_[1] + t_out_*direction_[1];
        t_max_ = t_out_ + 1;
        t_max_out = t_max_;
        t_out_out = t_out_;
    }


    // Compute x=(s,c) for a given t (t is the parameter for the piecewise linear curve)
    void CurveInSCPlane::computeXOfT(double* x_of_t, const double t) const {
        if (t <= t_out_) {
            x_of_t[0] = x_[0] + t*direction_[0];
            x_of_t[1] = x_[1] + t*direction_[1];
        } else {
            x_of_t[0] = 1/(t_max_-t_out_)*((t_max_ - t)*x_out_[0] + end_point_[0]*(t - t_out_));
            x_of_t[1] = 1/(t_max_-t_out_)*((t_max_ - t)*x_out_[1] + end_point_[1]*(t - t_out_));
        }
    }

} // Anonymous namespace


/* Local Variables:    */
/* c-basic-offset:4    */
/* End:                */
