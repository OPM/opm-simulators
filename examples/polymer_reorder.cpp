/*
  Copyright 2012 SINTEF ICT, Applied Mathematics.
  Copyright 2012 Statoil ASA.

  This file is part of the Open Porous Media Project (OPM).

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

#if HAVE_CONFIG_H
#include "config.h"
#endif // HAVE_CONFIG_H

#include <opm/core/pressure/IncompTpfa.hpp>
#include <opm/core/pressure/FlowBCManager.hpp>

#include <opm/core/grid.h>
#include <opm/core/GridManager.hpp>
#include <opm/core/newwells.h>
#include <opm/core/WellsManager.hpp>
#include <opm/core/utility/ErrorMacros.hpp>
#include <opm/core/utility/SimulatorTimer.hpp>
#include <opm/core/utility/StopWatch.hpp>
#include <opm/core/utility/Units.hpp>
#include <opm/core/utility/writeVtkData.hpp>
#include <opm/core/utility/miscUtilities.hpp>
#include <opm/core/utility/parameters/ParameterGroup.hpp>

#include <opm/core/fluid/IncompPropertiesBasic.hpp>
#include <opm/core/fluid/IncompPropertiesFromDeck.hpp>

#include <opm/core/linalg/LinearSolverUmfpack.hpp>
// #define EXPERIMENT_ISTL
#ifdef EXPERIMENT_ISTL
#include <opm/core/linalg/LinearSolverIstl.hpp>
#endif

#include <opm/polymer/TransportModelPolymer.hpp>
#include <opm/polymer/PolymerProperties.hpp>
#include <opm/polymer/polymerUtilities.hpp>

#include <boost/filesystem/convenience.hpp>
#include <boost/scoped_ptr.hpp>
#include <boost/lexical_cast.hpp>

#include <cassert>
#include <cstddef>

#include <algorithm>
#include <tr1/array>
#include <functional>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <iterator>
#include <vector>
#include <numeric>



class AdHocProps : public Opm::IncompPropertiesBasic
{
public:
    AdHocProps(const Opm::parameter::ParameterGroup& param, int dim, int num_cells)
        : Opm::IncompPropertiesBasic(param, dim, num_cells)
    {
        ASSERT(numPhases() == 2);
        sw_.resize(3);
        sw_[0] = 0.2;
        sw_[1] = 0.7;
        sw_[2] = 1.0;
        krw_.resize(3);
        krw_[0] = 0.0;
        krw_[1] = 0.7;
        krw_[2] = 1.0;
        so_.resize(2);
        so_[0] = 0.3;
        so_[1] = 0.8;
        kro_.resize(2);
        kro_[0] = 0.0;
        kro_[1] = 1.0;
    }

    virtual void relperm(const int n,
                         const double* s,
                         const int* /*cells*/,
                         double* kr,
                         double* dkrds) const
    {
        // ASSERT(dkrds == 0);
        // We assume two phases flow
        for (int i = 0; i < n; ++i) {
            kr[2*i] = krw(s[2*i]);
            kr[2*i+1] = kro(s[2*i+1]);
            if (dkrds != 0) {
                dkrds[2*i] = krw_dsw(s[2*i]);
                dkrds[2*i+3] = kro_dso(s[2*i+1]);
                dkrds[2*i+1] = -dkrds[2*i+3];
                dkrds[2*i+2] = -dkrds[2*i];
            }
        }
    }

    virtual void satRange(const int n,
                          const int* /*cells*/,
                          double* smin,
                          double* smax) const
    {
        const int np = 2;
        for (int i = 0; i < n; ++i) {
            smin[np*i + 0] = sw_[0];
            smax[np*i + 0] = sw_.back();
            smin[np*i + 1] = 1.0 - sw_[0];
            smax[np*i + 1] = 1.0 - sw_.back();
        }
    }

private:
    double krw(double s) const
    {
        return Opm::linearInterpolation(sw_, krw_, s);
    }

    double krw_dsw(double s) const
    {
        return Opm::linearInterpolationDerivative(sw_, krw_, s);
    }


    double kro(double s) const
    {
        return Opm::linearInterpolation(so_, kro_, s);
    }

    double kro_dso(double s) const
    {
        return Opm::linearInterpolationDerivative(so_, kro_, s);
    }

    std::vector<double> sw_;
    std::vector<double> krw_;
    std::vector<double> so_;
    std::vector<double> kro_;
};


class ReservoirState
{
public:
    ReservoirState(const UnstructuredGrid* g, const double init_sat = 0.0)
        : press_ (g->number_of_cells, 0.0),
          fpress_(g->number_of_faces, 0.0),
          flux_  (g->number_of_faces, 0.0),
          sat_   (2 * g->number_of_cells, 0.0),
          concentration_(g->number_of_cells, 0.0),
          cmax_(g->number_of_cells, 0.0)
    {
        for (int cell = 0; cell < g->number_of_cells; ++cell) {
            sat_[2*cell] = init_sat;
            sat_[2*cell + 1] = 1.0 - init_sat;
        }
    }

    enum ExtremalSat { MinSat, MaxSat };

    void setToMinimumWaterSat(const Opm::IncompPropertiesInterface& props)
    {
        const int n = props.numCells();
        std::vector<int> cells(n);
        for (int i = 0; i < n; ++i) {
            cells[i] = i;
        }
        setWaterSat(cells, props, MinSat);
    }

    void setWaterSat(const std::vector<int>& cells,
                     const Opm::IncompPropertiesInterface& props,
                     ExtremalSat es)
    {
        const int n = cells.size();
        std::vector<double> smin(2*n);
        std::vector<double> smax(2*n);
        props.satRange(n, &cells[0], &smin[0], &smax[0]);
        const double* svals = (es == MinSat) ? &smin[0] : &smax[0];
        for (int ci = 0; ci < n; ++ci) {
            const int cell = cells[ci];
            sat_[2*cell] = svals[2*ci];
            sat_[2*cell + 1] = 1.0 - sat_[2*cell];
        }
    }

    // Initialize saturations so that there is water below woc,
    // and oil above.
    // TODO: add 'anitialiasing', obtaining a more precise woc
    //       by f. ex. subdividing cells cut by the woc.
    void initWaterOilContact(const UnstructuredGrid& grid,
                             const Opm::IncompPropertiesInterface& props,
                             const double woc)
    {
        // Find out which cells should have water and which should have oil.
        std::vector<int> oil;
        std::vector<int> water;
        const int num_cells = grid.number_of_cells;
        oil.reserve(num_cells);
        water.reserve(num_cells);
        const int dim = grid.dimensions;
        for (int c = 0; c < num_cells; ++c) {
            const double z = grid.cell_centroids[dim*c + dim - 1];
            if (z > woc) {
                // Z is depth, we put water in the deepest parts
                // (even if oil is heavier...).
                water.push_back(c);
            } else {
                oil.push_back(c);
            }
        }

        // Set saturations.
        setWaterSat(oil, props, MinSat);
        setWaterSat(water, props, MaxSat);
    }

    int numPhases() const { return sat_.size()/press_.size(); }

    std::vector<double>& pressure    () { return press_ ; }
    std::vector<double>& facepressure() { return fpress_; }
    std::vector<double>& faceflux    () { return flux_  ; }
    std::vector<double>& saturation  () { return sat_   ; }
    std::vector<double>& concentration() { return concentration_; }
    std::vector<double>& cmax() { return cmax_; }

    const std::vector<double>& pressure    () const { return press_ ; }
    const std::vector<double>& facepressure() const { return fpress_; }
    const std::vector<double>& faceflux    () const { return flux_  ; }
    const std::vector<double>& saturation  () const { return sat_   ; }
    const std::vector<double>& concentration() const { return concentration_; }
    const std::vector<double>& cmax() const { return cmax_; }

private:
    std::vector<double> press_ ;
    std::vector<double> fpress_;
    std::vector<double> flux_  ;
    std::vector<double> sat_   ;
    std::vector<double> concentration_;
    std::vector<double> cmax_;
};




static void outputState(const UnstructuredGrid& grid,
                        const ReservoirState& state,
                        const int step,
                        const std::string& output_dir)
{
    // Write data in VTK format.
    std::ostringstream vtkfilename;
    vtkfilename << output_dir << "/output-" << std::setw(3) << std::setfill('0') << step << ".vtu";
    std::ofstream vtkfile(vtkfilename.str().c_str());
    if (!vtkfile) {
        THROW("Failed to open " << vtkfilename.str());
    }
    Opm::DataMap dm;
    dm["saturation"] = &state.saturation();
    dm["pressure"] = &state.pressure();
    dm["concentration"] = &state.concentration();
    dm["cmax"] = &state.cmax();
    std::vector<double> cell_velocity;
    Opm::estimateCellVelocity(grid, state.faceflux(), cell_velocity);
    dm["velocity"] = &cell_velocity;
    Opm::writeVtkData(grid, dm, vtkfile);

    // Write data (not grid) in Matlab format
    for (Opm::DataMap::const_iterator it = dm.begin(); it != dm.end(); ++it) {
        std::ostringstream fname;
        fname << output_dir << "/" << it->first << "-" << std::setw(3) << std::setfill('0') << step << ".dat";
        std::ofstream file(fname.str().c_str());
        if (!file) {
            THROW("Failed to open " << fname.str());
        }
        const std::vector<double>& d = *(it->second);
        std::copy(d.begin(), d.end(), std::ostream_iterator<double>(file, "\n"));
    }
}




class PolymerInflow
{
public:
    PolymerInflow(const double starttime,
                  const double endtime,
                  const double amount)
        : stime_(starttime), etime_(endtime), amount_(amount)
    {
    }
    double operator()(double time)
    {
        if (time >= stime_ && time < etime_) {
            return amount_;
        } else {
            return 0.0;
        }
    }
private:
    double stime_;
    double etime_;
    double amount_;
};



static void outputWaterCut(const Opm::Watercut& watercut,
                           const std::string& output_dir)
{
    // Write water cut curve.
    std::string fname = output_dir  + "/watercut.txt";
    std::ofstream os(fname.c_str());
    if (!os) {
        THROW("Failed to open " << fname);
    }
    watercut.write(os);
}


// ----------------- Main program -----------------
int
main(int argc, char** argv)
{
    std::cout << "\n================    Test program for incompressible two-phase flow with polymer    ===============\n\n";
    Opm::parameter::ParameterGroup param(argc, argv, false);
    std::cout << "---------------    Reading parameters     ---------------" << std::endl;

    // Reading various control parameters.
    const bool output = param.getDefault("output", true);
    std::string output_dir;
    if (output) {
        output_dir = param.getDefault("output_dir", std::string("output"));
        // Ensure that output dir exists
        boost::filesystem::path fpath(output_dir);
        create_directories(fpath);
    }

    // If we have a "deck_filename", grid and props will be read from that.
    bool use_deck = param.has("deck_filename");
    boost::scoped_ptr<Opm::GridManager> grid;
    boost::scoped_ptr<Opm::IncompPropertiesInterface> props;
    boost::scoped_ptr<Opm::WellsManager> wells;
    Opm::SimulatorTimer simtimer;
    double water_oil_contact = 0.0;
    bool woc_set = false;
    Opm::PolymerProperties polydata;
    if (use_deck) {
        std::string deck_filename = param.get<std::string>("deck_filename");
        Opm::EclipseGridParser deck(deck_filename);
        // Grid init
        grid.reset(new Opm::GridManager(deck));
        // Rock and fluid init
        const int* gc = grid->c_grid()->global_cell;
        std::vector<int> global_cell(gc, gc + grid->c_grid()->number_of_cells);
        props.reset(new Opm::IncompPropertiesFromDeck(deck, global_cell));
        // props.reset(new AdHocProps(param, grid->c_grid()->dimensions, grid->c_grid()->number_of_cells));
        // Wells init.
        wells.reset(new Opm::WellsManager(deck, *grid->c_grid(), props->permeability()));
        // Timer init.
        if (deck.hasField("TSTEP")) {
            simtimer.init(deck);
        } else {
            simtimer.init(param);
        }
        // Water-oil contact.
        if (deck.hasField("EQUIL")) {
            water_oil_contact = deck.getEQUIL().equil[0].water_oil_contact_depth_;
            woc_set = true;
        } else if (param.has("water_oil_contact")) {
            water_oil_contact = param.get<double>("water_oil_contact");
            woc_set = true;
        }
        polydata.readFromDeck(deck);
    } else {
        // Grid init.
        const int nx = param.getDefault("nx", 100);
        const int ny = param.getDefault("ny", 100);
        const int nz = param.getDefault("nz", 1);
        const double dx = param.getDefault("dx", 1.0);
        const double dy = param.getDefault("dy", 1.0);
        const double dz = param.getDefault("dz", 1.0);
        grid.reset(new Opm::GridManager(nx, ny, nz, dx, dy, dz));
        // Rock and fluid init.
        // props.reset(new Opm::IncompPropertiesBasic(param, grid->c_grid()->dimensions, grid->c_grid()->number_of_cells));
        props.reset(new AdHocProps(param, grid->c_grid()->dimensions, grid->c_grid()->number_of_cells));
        // Wells init.
        wells.reset(new Opm::WellsManager());
        // Timer init.
        simtimer.init(param);
        if (param.has("water_oil_contact")) {
            water_oil_contact = param.get<double>("water_oil_contact");
            woc_set = true;
        }
        // Setting polydata defaults to mimic a simple example case.
        double c_max = param.getDefault("c_max_limit", 5.0);
        double mix_param = param.getDefault("mix_param", 1.0);
        double rock_density = param.getDefault("rock_density", 1000.0);
        double dead_pore_vol = param.getDefault("dead_pore_vol", 0.15);
        std::vector<double> c_vals_visc(2, -1e100);
        c_vals_visc[0] = 0.0;
        c_vals_visc[1] = 7.0;
        std::vector<double> visc_mult_vals(2, -1e100);
        visc_mult_vals[0] = 1.0;
        // polydata.visc_mult_vals[1] = param.getDefault("c_max_viscmult", 30.0);
        visc_mult_vals[1] = 20.0;
        std::vector<double> c_vals_ads(3, -1e100);
        c_vals_ads[0] = 0.0;
        c_vals_ads[1] = 2.0;
        c_vals_ads[2] = 8.0;
        std::vector<double> ads_vals(3, -1e100);
        ads_vals[0] = 0.0;
        // polydata.ads_vals[1] = param.getDefault("c_max_ads", 0.0025);
        ads_vals[1] = 0.0015;
        ads_vals[2] = 0.0025;
        polydata.set(c_max, mix_param, rock_density, dead_pore_vol, c_vals_visc, visc_mult_vals, c_vals_ads, ads_vals);
    }



    double poly_start = param.getDefault("poly_start_days", 300.0)*Opm::unit::day;
    double poly_end = param.getDefault("poly_end_days", 800.0)*Opm::unit::day;
    double poly_amount = param.getDefault("poly_amount", polydata.cMax());
    PolymerInflow poly_inflow(poly_start, poly_end, poly_amount);

    // Extra rock init.
    std::vector<double> porevol;
    computePorevolume(*grid->c_grid(), *props, porevol);
    double tot_porevol = std::accumulate(porevol.begin(), porevol.end(), 0.0);

    // Gravity init.
    double gravity[3] = { 0.0 };
    double g = param.getDefault("gravity", 0.0);
    bool use_gravity = g != 0.0;
    if (use_gravity) {
        gravity[grid->c_grid()->dimensions - 1] = g;
        if (props->density()[0] == props->density()[1]) {
            std::cout << "**** Warning: nonzero gravity, but zero density difference." << std::endl;
        }
    }

    // Solvers init.
    // Pressure solver.
#ifdef EXPERIMENT_ISTL
    Opm::LinearSolverIstl linsolver(param);
#else
    Opm::LinearSolverUmfpack linsolver;
#endif // EXPERIMENT_ISTL
    const double *grav = use_gravity ? &gravity[0] : 0;
    Opm::IncompTpfa psolver(*grid->c_grid(), props->permeability(), grav, linsolver);
    // Reordering solver.
    const double nltol = param.getDefault("nl_tolerance", 1e-9);
    const int maxit = param.getDefault("nl_maxiter", 30);
    Opm::TransportModelPolymer::SingleCellMethod method;
    std::string method_string = param.getDefault("single_cell_method", std::string("Bracketing"));
    if (method_string == "Bracketing") {
        method = Opm::TransportModelPolymer::Bracketing;
    } else if (method_string == "Newton") {
        method = Opm::TransportModelPolymer::Newton;
    } else {
        THROW("Unknown method: " << method_string);
    }
    Opm::TransportModelPolymer tmodel(*grid->c_grid(), props->porosity(), &porevol[0], *props, polydata,
                                      method, nltol, maxit);

    // Boundary conditions.
    Opm::FlowBCManager bcs;

    // State-related and source-related variables init.
    int num_cells = grid->c_grid()->number_of_cells;
    std::vector<double> totmob;
    std::vector<double> omega; // Will remain empty if no gravity.
    double init_sat = param.getDefault("init_sat", 0.0);
    ReservoirState state(grid->c_grid(), init_sat);
    if (!param.has("init_sat")) {
        state.setToMinimumWaterSat(*props);
    }
    // We need a separate reorder_sat, because the reorder
    // code expects a scalar sw, not both sw and so.
    std::vector<double> reorder_sat(num_cells);
    std::vector<double> src(num_cells, 0.0);
    int scenario = param.getDefault("scenario", woc_set ? 4 : 0);
    switch (scenario) {
    case 0:
        {
            std::cout << "==== Scenario 0: simple wells or single-cell source and sink.\n";
            if (wells->c_wells()) {
                Opm::wellsToSrc(*wells->c_wells(), num_cells, src);
            } else {
                double flow_per_sec = 0.1*tot_porevol/Opm::unit::day;
                if (param.has("injection_rate_per_day")) {
                    flow_per_sec = param.get<double>("injection_rate_per_day")/Opm::unit::day;
                }
                src[0] = flow_per_sec;
                src[num_cells - 1] = -flow_per_sec;
            }
            break;
        }
    case 1:
        {
            std::cout << "==== Scenario 1: half source, half sink.\n";
            double flow_per_sec = 0.1*porevol[0]/Opm::unit::day;
            std::fill(src.begin(), src.begin() + src.size()/2, flow_per_sec);
            std::fill(src.begin() + src.size()/2, src.end(), -flow_per_sec);
            break;
        }
    case 2:
        {
            std::cout << "==== Scenario 2: gravity convection.\n";
            if (!use_gravity) {
                std::cout << "**** Warning: running gravity convection scenario, but gravity is zero." << std::endl;
            }
            if (use_deck) {
                std::cout << "**** Warning: running gravity convection scenario, which expects a cartesian grid."
                          << std::endl;
            }
            if (grid->c_grid()->cartdims[2] <= 1) {
                std::cout << "**** Warning: running gravity convection scenario, which expects nz > 1." << std::endl;
            }
            std::vector<int> left_cells;
            left_cells.reserve(num_cells/2);
            const int *glob_cell = grid->c_grid()->global_cell;
            for (int cell = 0; cell < num_cells; ++cell) {
                const int* cd = grid->c_grid()->cartdims;
                const int gc = glob_cell == 0 ? cell : glob_cell[cell];
                bool left = (gc % cd[0]) < cd[0]/2;
                if (left) {
                    left_cells.push_back(cell);
                }
            }
            state.setWaterSat(left_cells, *props, ReservoirState::MaxSat);
            break;
        }
    case 3:
        {
            std::cout << "==== Scenario 3: gravity segregation.\n";
            if (!use_gravity) {
                std::cout << "**** Warning: running gravity segregation scenario, but gravity is zero." << std::endl;
            }
            if (use_deck) {
                std::cout << "**** Warning: running gravity segregation scenario, which expects a cartesian grid."
                          << std::endl;
            }
            if (grid->c_grid()->cartdims[2] <= 1) {
                std::cout << "**** Warning: running gravity segregation scenario, which expects nz > 1." << std::endl;
            }
            std::vector<double>& sat = state.saturation();
            const int *glob_cell = grid->c_grid()->global_cell;
            // Water on top
            for (int cell = 0; cell < num_cells; ++cell) {
                const int* cd = grid->c_grid()->cartdims;
                const int gc = glob_cell == 0 ? cell : glob_cell[cell];
                bool top = (gc / cd[0] / cd[1]) < cd[2]/2;
                sat[2*cell] = top ? 1.0 : 0.0;
                sat[2*cell + 1 ] = 1.0 - sat[2*cell];
            }
            break;
        }
    case 4:
        {
            std::cout << "==== Scenario 4: water-oil contact and simple wells or sources\n";
            if (!use_gravity) {
                std::cout << "**** Warning: initializing segregated water and oil zones, but gravity is zero." << std::endl;
            }
            state.initWaterOilContact(*grid->c_grid(), *props, water_oil_contact);
            if (wells->c_wells()) {
                Opm::wellsToSrc(*wells->c_wells(), num_cells, src);
            } else {
                double flow_per_sec = 0.01*tot_porevol/Opm::unit::day;
                src[0] = flow_per_sec;
                src[grid->c_grid()->number_of_cells - 1] = -flow_per_sec;
            }
            break;
        }
    default:
        {
            THROW("==== Scenario " << scenario << " is unknown.");
        }
    }
    std::vector<double> reorder_src = src;

    // Dirichlet boundary conditions.
    if (param.getDefault("use_pside", false)) {
        int pside = param.get<int>("pside");
        double pside_pressure = param.get<double>("pside_pressure");
        bcs.pressureSide(*grid->c_grid(), Opm::FlowBCManager::Side(pside), pside_pressure);
    }

    // The allcells vector is used in calls to computeTotalMobility()
    // and computeTotalMobilityOmega().
    std::vector<int> allcells(num_cells);
    for (int cell = 0; cell < num_cells; ++cell) {
        allcells[cell] = cell;
    }

    // Warn if any parameters are unused.
    if (param.anyUnused()) {
        std::cout << "--------------------   Unused parameters:   --------------------\n";
        param.displayUsage();
        std::cout << "----------------------------------------------------------------" << std::endl;
    }

    // Write parameters used for later reference.
    if (output) {
        param.writeParam(output_dir + "/spu_2p.param");
    }

    // Main simulation loop.
    Opm::time::StopWatch pressure_timer;
    double ptime = 0.0;
    Opm::time::StopWatch transport_timer;
    double ttime = 0.0;
    Opm::time::StopWatch total_timer;
    total_timer.start();
    std::cout << "\n\n================    Starting main simulation loop     ===============" << std::endl;
    double init_satvol[2] = { 0.0 };
    double satvol[2] = { 0.0 };
    double injected[2] = { 0.0 };
    double produced[2] = { 0.0 };
    double tot_injected[2] = { 0.0 };
    double tot_produced[2] = { 0.0 };
    Opm::computeSaturatedVol(porevol, state.saturation(), init_satvol);
    std::cout << "\nInitial saturations are    " << init_satvol[0]/tot_porevol
              << "    " << init_satvol[1]/tot_porevol << std::endl;
    Opm::Watercut watercut;
    watercut.push(0.0, 0.0, 0.0);
    for (; !simtimer.done(); ++simtimer) {
        // Report timestep and (optionally) write state to disk.
        simtimer.report(std::cout);
        if (output) {
            outputState(*grid->c_grid(), state, simtimer.currentStepNum(), output_dir);
        }

        // Solve pressure.
        if (use_gravity) {
            computeTotalMobilityOmega(*props, polydata, allcells, state.saturation(), state.concentration(),
                                      totmob, omega);
        } else {
            computeTotalMobility(*props, polydata, allcells, state.saturation(), state.concentration(),
                                 totmob);
        }
        pressure_timer.start();
        psolver.solve(totmob, omega, src, bcs.c_bcs(), state.pressure(), state.faceflux());
        pressure_timer.stop();
        double pt = pressure_timer.secsSinceStart();
        std::cout << "Pressure solver took:  " << pt << " seconds." << std::endl;
        ptime += pt;

        // Process transport sources (to include bdy terms).
        Opm::computeTransportSource(*grid->c_grid(), src, state.faceflux(), 1.0, reorder_src);

        // Find inflow rate.
        const double current_time = simtimer.currentTime();
        const double stepsize = simtimer.currentStepLength();
        const double inflowc0 = poly_inflow(current_time + 1e-5*stepsize);
        const double inflowc1 = poly_inflow(current_time + (1.0 - 1e-5)*stepsize);
        if (inflowc0 != inflowc1) {
            std::cout << "**** Warning: polymer inflow rate changes during timestep. Using rate near start of step.";
        }
        const double inflow_c = inflowc0;

        // Solve transport.
        transport_timer.start();
        Opm::toWaterSat(state.saturation(), reorder_sat);
        tmodel.solve(&state.faceflux()[0], &reorder_src[0], stepsize, inflow_c,
                     &reorder_sat[0], &state.concentration()[0], &state.cmax()[0]);
        Opm::toBothSat(reorder_sat, state.saturation());
        transport_timer.stop();
        double tt = transport_timer.secsSinceStart();
        std::cout << "Transport solver took: " << tt << " seconds." << std::endl;
        ttime += tt;

        // Report volume balances.
        Opm::computeSaturatedVol(porevol, state.saturation(), satvol);
        Opm::computeInjectedProduced(*props, state.saturation(), src, simtimer.currentStepLength(), injected, produced);
        tot_injected[0] += injected[0];
        tot_injected[1] += injected[1];
        tot_produced[0] += produced[0];
        tot_produced[1] += produced[1];
        std::cout.precision(5);
        const int width = 18;
        std::cout << "\nVolume balance report (all numbers relative to total pore volume).\n";
        std::cout << "    Saturated volumes:     "
                  << std::setw(width) << satvol[0]/tot_porevol
                  << std::setw(width) << satvol[1]/tot_porevol << std::endl;
        std::cout << "    Injected volumes:      "
                  << std::setw(width) << injected[0]/tot_porevol
                  << std::setw(width) << injected[1]/tot_porevol << std::endl;
        std::cout << "    Produced volumes:      "
                  << std::setw(width) << produced[0]/tot_porevol
                  << std::setw(width) << produced[1]/tot_porevol << std::endl;
        std::cout << "    Total inj volumes:     "
                  << std::setw(width) << tot_injected[0]/tot_porevol
                  << std::setw(width) << tot_injected[1]/tot_porevol << std::endl;
        std::cout << "    Total prod volumes:    "
                  << std::setw(width) << tot_produced[0]/tot_porevol
                  << std::setw(width) << tot_produced[1]/tot_porevol << std::endl;
        std::cout << "    In-place + prod - inj: "
                  << std::setw(width) << (satvol[0] + tot_produced[0] - tot_injected[0])/tot_porevol
                  << std::setw(width) << (satvol[1] + tot_produced[1] - tot_injected[1])/tot_porevol << std::endl;
        std::cout << "    Init - now - pr + inj: "
                  << std::setw(width) << (init_satvol[0] - satvol[0] - tot_produced[0] + tot_injected[0])/tot_porevol
                  << std::setw(width) << (init_satvol[1] - satvol[1] - tot_produced[1] + tot_injected[1])/tot_porevol
                  << std::endl;
        std::cout.precision(8);

        watercut.push(simtimer.currentTime() + simtimer.currentStepLength(),
                      produced[0]/(produced[0] + produced[1]),
                      tot_produced[0]/tot_porevol);
    }
    total_timer.stop();

    std::cout << "\n\n================    End of simulation     ===============\n"
              << "Total time taken: " << total_timer.secsSinceStart()
              << "\n  Pressure time:  " << ptime
              << "\n  Transport time: " << ttime << std::endl;

    if (output) {
        outputState(*grid->c_grid(), state, simtimer.currentStepNum(), output_dir);
        outputWaterCut(watercut, output_dir);
    }
}
