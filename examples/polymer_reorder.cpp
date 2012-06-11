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
#include <opm/core/wells/WellsManager.hpp>
#include <opm/core/utility/ErrorMacros.hpp>
#include <opm/core/utility/initState.hpp>
#include <opm/core/utility/SimulatorTimer.hpp>
#include <opm/core/utility/StopWatch.hpp>
#include <opm/core/utility/Units.hpp>
#include <opm/core/utility/writeVtkData.hpp>
#include <opm/core/utility/miscUtilities.hpp>
#include <opm/core/utility/parameters/ParameterGroup.hpp>

#include <opm/core/fluid/IncompPropertiesBasic.hpp>
#include <opm/core/fluid/IncompPropertiesFromDeck.hpp>
#include <opm/core/fluid/RockCompressibility.hpp>

#include <opm/core/linalg/LinearSolverFactory.hpp>

#include <opm/core/transport/transport_source.h>
#include <opm/core/transport/CSRMatrixUmfpackSolver.hpp>
#include <opm/core/transport/NormSupport.hpp>
#include <opm/core/transport/ImplicitAssembly.hpp>
#include <opm/core/transport/ImplicitTransport.hpp>
#include <opm/core/transport/JacobianSystem.hpp>
#include <opm/core/transport/CSRMatrixBlockAssembler.hpp>
#include <opm/core/transport/SinglePointUpwindTwoPhase.hpp>

#include <opm/core/utility/ColumnExtract.hpp>

#include <opm/polymer/PolymerState.hpp>
#include <opm/polymer/SinglePointUpwindTwoPhasePolymer.hpp>
#include <opm/polymer/GravityColumnSolverPolymer.hpp>
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




static void outputState(const UnstructuredGrid& grid,
                        const Opm::PolymerState& state,
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
    dm["cmax"] = &state.maxconcentration();
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


static void outputWellReport(const Opm::WellReport& wellreport,
                             const std::string& output_dir)
{
    // Write well report.
    std::string fname = output_dir  + "/wellreport.txt";
    std::ofstream os(fname.c_str());
    if (!os) {
        THROW("Failed to open " << fname);
    }
    wellreport.write(os);
}



// --------------- Types needed to define transport solver ---------------

class PolymerFluid2pWrappingProps
{
public:
    PolymerFluid2pWrappingProps(const Opm::IncompPropertiesInterface& props, const Opm::PolymerProperties& polyprops)
        : props_(props),
          polyprops_(polyprops),
          smin_(props.numCells()*props.numPhases()),
          smax_(props.numCells()*props.numPhases())
    {
        if (props.numPhases() != 2) {
            THROW("PolymerFluid2pWrapper requires 2 phases.");
        }
        const int num_cells = props.numCells();
        std::vector<int> cells(num_cells);
        for (int c = 0; c < num_cells; ++c) {
            cells[c] = c;
        }
        props.satRange(num_cells, &cells[0], &smin_[0], &smax_[0]);
    }

    double density(int phase) const
    {
        return props_.density()[phase];
    }

    template <class PolyC,
              class CAds,
              class DCAdsDc>

    void adsorption(const PolyC& c, const PolyC& cmax, CAds& cads, DCAdsDc& dcadsdc)
    {
        polyprops_.adsorptionWithDer(c, cmax, cads, dcadsdc);
    }

    const double* porosity() const
    {
        return props_.porosity();
    }

    double deadporespace() const
    {
        return polyprops_.deadPoreVol();
    }

    double rockdensity() const
    {
        return polyprops_.rockDensity();
    }

    template <class Sat,
              class PolyC,
              class Mob,
              class DMobDs,
              class DMobWatDc>
    void mobility(int cell, const Sat& s, const PolyC& c, const PolyC& cmax,
                  Mob& mob, DMobDs& dmobds, DMobWatDc& dmobwatdc) const
    {
        const double* visc = props_.viscosity();
        double relperm[2];
        double drelpermds[4];
        props_.relperm(1, &s[0], &cell, relperm, drelpermds);
        polyprops_.effectiveMobilitiesWithDer(c, cmax, visc, relperm, drelpermds, mob, dmobds, dmobwatdc);
    }

    template <class Sat,
              class Pcap,
              class DPcap>
    void pc(int c, const Sat& s, Pcap& pcap, DPcap& dpcap) const
    {
        double pcow[2];
        double dpcow[4];
        props_.capPress(1, &s[0], &c, pcow, dpcow);
        pcap = pcow[0];
        ASSERT(pcow[1] == 0.0);
        dpcap = dpcow[0];
        ASSERT(dpcow[1] == 0.0);
        ASSERT(dpcow[2] == 0.0);
        ASSERT(dpcow[3] == 0.0);
    }

    double s_min(int c) const
    {
        return smin_[2*c + 0];
    }

    double s_max(int c) const
    {
        return smax_[2*c + 0];
    }

    double cMax() const
    {
        return polyprops_.cMax();
    }

    template <class PolyC,
              class Mc,
              class DMcDc>
    void computeMc(const PolyC& c, Mc& mc,
                   DMcDc& dmcdc) const
    {
        polyprops_.computeMcWithDer(c, mc, dmcdc);
    }

private:
    const Opm::IncompPropertiesInterface& props_;
    const Opm::PolymerProperties& polyprops_;
    std::vector<double> smin_;
    std::vector<double> smax_;
};

typedef PolymerFluid2pWrappingProps TwophaseFluidPolymer;
typedef Opm::SinglePointUpwindTwoPhasePolymer<TwophaseFluidPolymer> FluxModel;

using namespace Opm::ImplicitTransportDefault;

typedef NewtonVectorCollection< ::std::vector<double> >      NVecColl;
typedef JacobianSystem        < struct CSRMatrix, NVecColl > JacSys;

template <class Vector>
class MaxNorm {
public:
    static double
    norm(const Vector& v) {
        return AccumulationNorm <Vector, MaxAbs>::norm(v);
    }
};

typedef Opm::ImplicitTransport<FluxModel,
                               JacSys        ,
                               MaxNorm       ,
                               VectorNegater ,
                               VectorZero    ,
                               MatrixZero    ,
                               VectorAssign  > TransportSolver;





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



// ----------------- Main program -----------------
int
main(int argc, char** argv)
{
    using namespace Opm;

    std::cout << "\n================    Test program for incompressible two-phase flow with polymer    ===============\n\n";
    Opm::parameter::ParameterGroup param(argc, argv, false);
    std::cout << "---------------    Reading parameters     ---------------" << std::endl;

    // Reading various control parameters.
    const bool guess_old_solution = param.getDefault("guess_old_solution", false);
    const bool use_reorder = param.getDefault("use_reorder", true);
    const bool output = param.getDefault("output", true);
    std::string output_dir;
    int output_interval = 1;
    if (output) {
        output_dir = param.getDefault("output_dir", std::string("output"));
        // Ensure that output dir exists
        boost::filesystem::path fpath(output_dir);
        try {
            create_directories(fpath);
        }
        catch (...) {
            THROW("Creating directories failed: " << fpath);
        }
        output_interval = param.getDefault("output_interval", output_interval);
    }
    const int num_transport_substeps = param.getDefault("num_transport_substeps", 1);

    // If we have a "deck_filename", grid and props will be read from that.
    bool use_deck = param.has("deck_filename");
    boost::scoped_ptr<Opm::GridManager> grid;
    boost::scoped_ptr<Opm::IncompPropertiesInterface> props;
    boost::scoped_ptr<Opm::WellsManager> wells;
    boost::scoped_ptr<Opm::RockCompressibility> rock_comp;
    Opm::SimulatorTimer simtimer;
    Opm::PolymerState state;
    Opm::PolymerProperties polyprop;
    bool check_well_controls = false;
    int max_well_control_iterations = 0;
    double gravity[3] = { 0.0 };
    if (use_deck) {
        std::string deck_filename = param.get<std::string>("deck_filename");
        Opm::EclipseGridParser deck(deck_filename);
        // Grid init
        grid.reset(new Opm::GridManager(deck));
        // Rock and fluid init
        const int* gc = grid->c_grid()->global_cell;
        std::vector<int> global_cell(gc, gc + grid->c_grid()->number_of_cells);
        props.reset(new Opm::IncompPropertiesFromDeck(deck, global_cell));
        // Wells init.
        wells.reset(new Opm::WellsManager(deck, *grid->c_grid(), props->permeability()));
        check_well_controls = param.getDefault("check_well_controls", false);
        max_well_control_iterations = param.getDefault("max_well_control_iterations", 10);
        // Timer init.
        if (deck.hasField("TSTEP")) {
            simtimer.init(deck);
        } else {
            simtimer.init(param);
        }
        // Rock compressibility.
        rock_comp.reset(new Opm::RockCompressibility(deck));
        // Gravity.
        gravity[2] = deck.hasField("NOGRAV") ? 0.0 : Opm::unit::gravity;
        // Init state variables (saturation and pressure).
        initStateFromDeck(*grid->c_grid(), *props, deck, gravity[2], state);
        // Init polymer properties.
        polyprop.readFromDeck(deck);
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
        props.reset(new Opm::IncompPropertiesBasic(param, grid->c_grid()->dimensions, grid->c_grid()->number_of_cells));
        // Wells init.
        wells.reset(new Opm::WellsManager());
        // Timer init.
        simtimer.init(param);
        // Rock compressibility.
        rock_comp.reset(new Opm::RockCompressibility(param));
        // Gravity.
        gravity[2] = param.getDefault("gravity", 0.0);
        // Init state variables (saturation and pressure).
        initStateBasic(*grid->c_grid(), *props, param, gravity[2], state);
        // Init Polymer state
        if (param.has("poly_init")) {
            double poly_init = param.getDefault("poly_init", 0.0);
            for (int cell = 0; cell < grid->c_grid()->number_of_cells; ++cell) {
                double smin[2], smax[2];
                props->satRange(1, &cell, smin, smax);
                if (state.saturation()[2*cell] > 0.5*(smin[0] + smax[0])) {
                    state.concentration()[cell] = poly_init;
                    state.maxconcentration()[cell] = poly_init;
                } else {
                    state.saturation()[2*cell + 0] = 0.;
                    state.saturation()[2*cell + 1] = 1.;
                    state.concentration()[cell] = 0.;
                    state.maxconcentration()[cell] = 0.;
                }
            }
        }
        // Init polymer properties.
        // Setting defaults to provide a simple example case.
        double c_max = param.getDefault("c_max_limit", 5.0);
        double mix_param = param.getDefault("mix_param", 1.0);
        double rock_density = param.getDefault("rock_density", 1000.0);
        double dead_pore_vol = param.getDefault("dead_pore_vol", 0.15);
        double res_factor = param.getDefault("res_factor", 1.) ; // res_factor = 1 gives no change in permeability
        double c_max_ads = param.getDefault("c_max_ads", 1.);
        int ads_index = param.getDefault<int>("ads_index", Opm::PolymerProperties::NoDesorption);
        std::vector<double> c_vals_visc(2, -1e100);
        c_vals_visc[0] = 0.0;
        c_vals_visc[1] = 7.0;
        std::vector<double> visc_mult_vals(2, -1e100);
        visc_mult_vals[0] = 1.0;
        // polyprop.visc_mult_vals[1] = param.getDefault("c_max_viscmult", 30.0);
        visc_mult_vals[1] = 20.0;
        std::vector<double> c_vals_ads(3, -1e100);
        c_vals_ads[0] = 0.0;
        c_vals_ads[1] = 2.0;
        c_vals_ads[2] = 8.0;
        std::vector<double> ads_vals(3, -1e100);
        ads_vals[0] = 0.0;
        ads_vals[1] = 0.0015;
        ads_vals[2] = 0.0025;
        // ads_vals[1] = 0.0;
        // ads_vals[2] = 0.0;
        polyprop.set(c_max, mix_param, rock_density, dead_pore_vol, res_factor, c_max_ads,
                     static_cast<Opm::PolymerProperties::AdsorptionBehaviour>(ads_index),
                     c_vals_visc,  visc_mult_vals, c_vals_ads, ads_vals);
    }

    // Initialize polymer inflow function.
    double poly_start = param.getDefault("poly_start_days", 300.0)*Opm::unit::day;
    double poly_end = param.getDefault("poly_end_days", 800.0)*Opm::unit::day;
    double poly_amount = param.getDefault("poly_amount", polyprop.cMax());
    PolymerInflow poly_inflow(poly_start, poly_end, poly_amount);

    // Extra fluid init for transport solver.
    TwophaseFluidPolymer fluid(*props, polyprop);

    // Warn if gravity but no density difference.
    bool use_gravity = (gravity[0] != 0.0 || gravity[1] != 0.0 || gravity[2] != 0.0);
    if (use_gravity) {
        if (props->density()[0] == props->density()[1]) {
            std::cout << "**** Warning: nonzero gravity, but zero density difference." << std::endl;
        }
    }
    bool use_segregation_split = false;
    bool use_column_solver = false;
    bool use_gauss_seidel_gravity = false;
    if (use_gravity && use_reorder) {
        use_segregation_split = param.getDefault("use_segregation_split", use_segregation_split);
        if (use_segregation_split) {
            use_column_solver = param.getDefault("use_column_solver", use_column_solver);
            if (use_column_solver) {
                use_gauss_seidel_gravity = param.getDefault("use_gauss_seidel_gravity", use_gauss_seidel_gravity);
            }
        }
    }

    // Check that rock compressibility is not used with solvers that do not handle it.
    int nl_pressure_maxiter = 0;
    double nl_pressure_tolerance = 0.0;
    if (rock_comp->isActive()) {
        if (!use_reorder) {
            THROW("Cannot run implicit (non-reordering) transport solver with rock compressibility yet.");
        }
        nl_pressure_maxiter = param.getDefault("nl_pressure_maxiter", 10);
        nl_pressure_tolerance = param.getDefault("nl_pressure_tolerance", 1.0); // in Pascal
    }

    // Source-related variables init.
    int num_cells = grid->c_grid()->number_of_cells;
    std::vector<double> totmob;
    std::vector<double> omega; // Will remain empty if no gravity.
    std::vector<double> rc; // Will remain empty if no rock compressibility.

    // Extra rock init.
    std::vector<double> porevol;
    if (rock_comp->isActive()) {
        computePorevolume(*grid->c_grid(), props->porosity(), *rock_comp, state.pressure(), porevol);
    } else {
        computePorevolume(*grid->c_grid(), props->porosity(), porevol);
    }
    double tot_porevol_init = std::accumulate(porevol.begin(), porevol.end(), 0.0);

    // We need a separate reorder_sat, because the reorder
    // code expects a scalar sw, not both sw and so.
    std::vector<double> reorder_sat(num_cells);
    std::vector<double> src(num_cells, 0.0);

    // Initialising src
    if (wells->c_wells()) {
        // Do nothing, wells will be the driving force, not source terms.
        // Opm::wellsToSrc(*wells->c_wells(), num_cells, src);
    } else {
        const double default_injection = use_gravity ? 0.0 : 0.1;
        const double flow_per_sec = param.getDefault<double>("injected_volume_per_day", default_injection)/Opm::unit::day;
        src[0] = flow_per_sec;
        src[num_cells - 1] = -flow_per_sec;
    }

    std::vector<double> reorder_src = src;

    // Boundary conditions.
    Opm::FlowBCManager bcs;
    if (param.getDefault("use_pside", false)) {
        int pside = param.get<int>("pside");
        double pside_pressure = param.get<double>("pside_pressure");
        bcs.pressureSide(*grid->c_grid(), Opm::FlowBCManager::Side(pside), pside_pressure);
    }

    // Solvers init.
    // Linear solver.
    Opm::LinearSolverFactory linsolver(param);
    // Pressure solver.
    const double *grav = use_gravity ? &gravity[0] : 0;
    Opm::IncompTpfa psolver(*grid->c_grid(), props->permeability(), grav, linsolver, wells->c_wells());
    // Reordering solver.
    const double nl_tolerance = param.getDefault("nl_tolerance", 1e-9);
    const int nl_maxiter = param.getDefault("nl_maxiter", 30);
    Opm::TransportModelPolymer::SingleCellMethod method;
    std::string method_string = param.getDefault("single_cell_method", std::string("Bracketing"));
    if (method_string == "Bracketing") {
        method = Opm::TransportModelPolymer::Bracketing;
    } else if (method_string == "Newton") {
        method = Opm::TransportModelPolymer::Newton;
    } else {
        THROW("Unknown method: " << method_string);
    }

    Opm::TransportModelPolymer reorder_model(*grid->c_grid(), props->porosity(), &porevol[0], *props, polyprop,
                                             method, nl_tolerance, nl_maxiter);

    if (use_gauss_seidel_gravity) {
        reorder_model.initGravity(grav);
    }
    // Non-reordering solver.
    FluxModel  fmodel(fluid, *grid->c_grid(), porevol, grav, guess_old_solution);
    if (use_gravity) {
        fmodel.initGravityTrans(*grid->c_grid(), psolver.getHalfTrans());
    }
    TransportSolver tsolver(fmodel);
    // Column-based gravity segregation solver.
    std::vector<std::vector<int> > columns;
    if (use_column_solver) {
        Opm::extractColumn(*grid->c_grid(), columns);
    }
    Opm::GravityColumnSolverPolymer<FluxModel, TwophaseFluidPolymer> colsolver(fmodel, fluid, *grid->c_grid(), nl_tolerance, nl_maxiter);

    // // // Not implemented for polymer.
    // // Control init.
    // Opm::ImplicitTransportDetails::NRReport  rpt;
    // Opm::ImplicitTransportDetails::NRControl ctrl;
    // if (!use_reorder || use_segregation_split) {
    //     ctrl.max_it = param.getDefault("max_it", 20);
    //     ctrl.verbosity = param.getDefault("verbosity", 0);
    //     ctrl.max_it_ls = param.getDefault("max_it_ls", 5);
    // }
    // // Linear solver init.
    // using Opm::ImplicitTransportLinAlgSupport::CSRMatrixUmfpackSolver;
    // CSRMatrixUmfpackSolver linsolve;

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
    double init_polymass = 0.0;
    double satvol[2] = { 0.0 };
    double polymass = 0.0;
    double polymass_adsorbed = 0.0;
    double injected[2] = { 0.0 };
    double produced[2] = { 0.0 };
    double polyinj = 0.0;
    double polyprod = 0.0;
    double tot_injected[2] = { 0.0 };
    double tot_produced[2] = { 0.0 };
    double tot_polyinj = 0.0;
    double tot_polyprod = 0.0;
    Opm::computeSaturatedVol(porevol, state.saturation(), init_satvol);
    std::cout << "\nInitial saturations are    " << init_satvol[0]/tot_porevol_init
              << "    " << init_satvol[1]/tot_porevol_init << std::endl;
    Opm::Watercut watercut;
    watercut.push(0.0, 0.0, 0.0);
    Opm::WellReport wellreport;
    std::vector<double> well_bhp;
    std::vector<double> well_perfrates;
    std::vector<double> fractional_flows;
    std::vector<double> well_resflows_phase;
    int num_wells = 0;
    if (wells->c_wells()) {
        num_wells = wells->c_wells()->number_of_wells;
        well_bhp.resize(num_wells, 0.0);
        well_perfrates.resize(wells->c_wells()->well_connpos[num_wells], 0.0);
        well_resflows_phase.resize((wells->c_wells()->number_of_phases)*(wells->c_wells()->number_of_wells), 0.0);
        wellreport.push(*props, *wells->c_wells(), state.saturation(), 0.0, well_bhp, well_perfrates);
    }
    for (; !simtimer.done(); ++simtimer) {
        // Report timestep and (optionally) write state to disk.
        simtimer.report(std::cout);
        if (output && (simtimer.currentStepNum() % output_interval == 0)) {
            outputState(*grid->c_grid(), state, simtimer.currentStepNum(), output_dir);
        }

        // Solve pressure.
        if (use_gravity) {
            computeTotalMobilityOmega(*props, polyprop, allcells, state.saturation(), state.concentration(), state.maxconcentration(),
                                      totmob, omega);
        } else {
            computeTotalMobility(*props, polyprop, allcells, state.saturation(), state.concentration(), state.maxconcentration(),
                                 totmob);
        }
        std::vector<double> wdp;
        if (wells->c_wells()) {
            Opm::computeWDP(*wells->c_wells(), *grid->c_grid(), state.saturation(), props->density(), gravity[2], true, wdp);
        }
        if (check_well_controls) {
            computeFractionalFlow(*props, allcells, state.saturation(), fractional_flows);
        }
        if (check_well_controls) {
            wells->applyExplicitReinjectionControls(well_resflows_phase, well_resflows_phase);
        }
        bool well_control_passed = !check_well_controls;
        int well_control_iteration = 0;
        do {
            pressure_timer.start();
            if (rock_comp->isActive()) {
                rc.resize(num_cells);
                std::vector<double> initial_pressure = state.pressure();
                std::vector<double> initial_porevolume(num_cells);
                computePorevolume(*grid->c_grid(), props->porosity(), *rock_comp, initial_pressure, initial_porevolume);
                std::vector<double> pressure_increment(num_cells + num_wells);
                std::vector<double> prev_pressure(num_cells + num_wells);
                for (int iter = 0; iter < nl_pressure_maxiter; ++iter) {

                    for (int cell = 0; cell < num_cells; ++cell) {
                        rc[cell] = rock_comp->rockComp(state.pressure()[cell]);
                    }
                    computePorevolume(*grid->c_grid(), props->porosity(), *rock_comp, state.pressure(), porevol);
                    std::copy(state.pressure().begin(), state.pressure().end(), prev_pressure.begin());
                    std::copy(well_bhp.begin(), well_bhp.end(), prev_pressure.begin() + num_cells);
                    // prev_pressure = state.pressure();

                    // compute pressure increment
                    psolver.solveIncrement(totmob, omega, src, wdp, bcs.c_bcs(), porevol, rc,
                                           prev_pressure, initial_porevolume, simtimer.currentStepLength(),
                                           pressure_increment);

                    double max_change = 0.0;
                    for (int cell = 0; cell < num_cells; ++cell) {
                        state.pressure()[cell] += pressure_increment[cell];
                        max_change = std::max(max_change, std::fabs(pressure_increment[cell]));
                    }
                    for (int well = 0; well < num_wells; ++well) {
                        well_bhp[well] += pressure_increment[num_cells + well];
                        max_change = std::max(max_change, std::fabs(pressure_increment[num_cells + well]));
                    }

                    std::cout << "Pressure iter " << iter << "   max change = " << max_change << std::endl;
                    if (max_change < nl_pressure_tolerance) {
                        break;
                    }
                }
                psolver.computeFaceFlux(totmob, omega, src, wdp, bcs.c_bcs(), state.pressure(), state.faceflux(),
                                        well_bhp, well_perfrates);
            } else {
                psolver.solve(totmob, omega, src, wdp, bcs.c_bcs(), state.pressure(), state.faceflux(),
                              well_bhp, well_perfrates);
            }
            pressure_timer.stop();
            double pt = pressure_timer.secsSinceStart();
            std::cout << "Pressure solver took:  " << pt << " seconds." << std::endl;
            ptime += pt;


            if (check_well_controls) {
                Opm::computePhaseFlowRatesPerWell(*wells->c_wells(),
                                                  fractional_flows,
                                                  well_perfrates,
                                                  well_resflows_phase);
                std::cout << "Checking well conditions." << std::endl;
                // For testing we set surface := reservoir
                well_control_passed = wells->conditionsMet(well_bhp, well_resflows_phase, well_resflows_phase);
                ++well_control_iteration;
                if (!well_control_passed && well_control_iteration > max_well_control_iterations) {
                    THROW("Could not satisfy well conditions in " << max_well_control_iterations << " tries.");
                }
                if (!well_control_passed) {
                    std::cout << "Well controls not passed, solving again." << std::endl;
                } else {
                    std::cout << "Well conditions met." << std::endl;
                }
            }
        } while (!well_control_passed);

        // Process transport sources (to include bdy terms and well flows).
        Opm::computeTransportSource(*grid->c_grid(), src, state.faceflux(), 1.0,
                                    wells->c_wells(), well_perfrates, reorder_src);


        // Find inflow rate.
        const double current_time = simtimer.currentTime();
        double stepsize = simtimer.currentStepLength();
        const double inflowc0 = poly_inflow(current_time + 1e-5*stepsize);
        const double inflowc1 = poly_inflow(current_time + (1.0 - 1e-5)*stepsize);
        if (inflowc0 != inflowc1) {
            std::cout << "**** Warning: polymer inflow rate changes during timestep. Using rate near start of step.";
        }
        const double inflow_c = inflowc0;


        // Solve transport.
        transport_timer.start();
        if (num_transport_substeps != 1) {
            stepsize /= double(num_transport_substeps);
            std::cout << "Making " << num_transport_substeps << " transport substeps." << std::endl;
        }
        for (int tr_substep = 0; tr_substep < num_transport_substeps; ++tr_substep) {
            if (use_reorder) {
                Opm::toWaterSat(state.saturation(), reorder_sat);
                reorder_model.solve(&state.faceflux()[0], &porevol[0], &reorder_src[0], stepsize, inflow_c,
                                    &reorder_sat[0], &state.concentration()[0], &state.maxconcentration()[0]);
                Opm::toBothSat(reorder_sat, state.saturation());
                Opm::computeInjectedProduced(*props, polyprop, state.saturation(), state.concentration(), state.maxconcentration(),
                                             reorder_src, simtimer.currentStepLength(), inflow_c,
                                             injected, produced, polyinj, polyprod);
                if (use_segregation_split) {
                    if (use_column_solver) {
                        if (use_gauss_seidel_gravity) {
                            reorder_model.solveGravity(columns, &porevol[0], stepsize, reorder_sat,
                                                       state.concentration(), state.maxconcentration());
                            Opm::toBothSat(reorder_sat, state.saturation());
                        } else {
                            colsolver.solve(columns, stepsize, state.saturation(), state.concentration(),
                                            state.maxconcentration());
                        }
                    } else {
                        THROW("use_segregation_split option for polymer is only implemented in the use_column_solver case.");
                    }
                }
            } else {
                THROW("Implicit transport solver not implemented for polymer.");
            }
        }
        transport_timer.stop();
        double tt = transport_timer.secsSinceStart();
        std::cout << "Transport solver took: " << tt << " seconds." << std::endl;
        ttime += tt;

        // Report volume balances.
        Opm::computeSaturatedVol(porevol, state.saturation(), satvol);
        polymass = Opm::computePolymerMass(porevol, state.saturation(), state.concentration(), polyprop.deadPoreVol());
        polymass_adsorbed = Opm::computePolymerAdsorbed(*props, polyprop, porevol, state.maxconcentration());
        tot_injected[0] += injected[0];
        tot_injected[1] += injected[1];
        tot_produced[0] += produced[0];
        tot_produced[1] += produced[1];
        tot_polyinj += polyinj;
        tot_polyprod += polyprod;
        std::cout.precision(5);
        const int width = 18;
        std::cout << "\nVolume and polymer mass balance: "
            "   water(pv)           oil(pv)       polymer(kg)\n";
        std::cout << "    Saturated volumes:     "
                  << std::setw(width) << satvol[0]/tot_porevol_init
                  << std::setw(width) << satvol[1]/tot_porevol_init
                  << std::setw(width) << polymass << std::endl;
        std::cout << "    Adsorbed volumes:      "
                  << std::setw(width) << 0.0
                  << std::setw(width) << 0.0
                  << std::setw(width) << polymass_adsorbed << std::endl;
        std::cout << "    Injected volumes:      "
                  << std::setw(width) << injected[0]/tot_porevol_init
                  << std::setw(width) << injected[1]/tot_porevol_init
                  << std::setw(width) << polyinj << std::endl;
        std::cout << "    Produced volumes:      "
                  << std::setw(width) << produced[0]/tot_porevol_init
                  << std::setw(width) << produced[1]/tot_porevol_init
                  << std::setw(width) << polyprod << std::endl;
        std::cout << "    Total inj volumes:     "
                  << std::setw(width) << tot_injected[0]/tot_porevol_init
                  << std::setw(width) << tot_injected[1]/tot_porevol_init
                  << std::setw(width) << tot_polyinj << std::endl;
        std::cout << "    Total prod volumes:    "
                  << std::setw(width) << tot_produced[0]/tot_porevol_init
                  << std::setw(width) << tot_produced[1]/tot_porevol_init
                  << std::setw(width) << tot_polyprod << std::endl;
        std::cout << "    In-place + prod - inj: "
                  << std::setw(width) << (satvol[0] + tot_produced[0] - tot_injected[0])/tot_porevol_init
                  << std::setw(width) << (satvol[1] + tot_produced[1] - tot_injected[1])/tot_porevol_init
                  << std::setw(width) << (polymass + tot_polyprod - tot_polyinj + polymass_adsorbed) << std::endl;
        std::cout << "    Init - now - pr + inj: "
                  << std::setw(width) << (init_satvol[0] - satvol[0] - tot_produced[0] + tot_injected[0])/tot_porevol_init
                  << std::setw(width) << (init_satvol[1] - satvol[1] - tot_produced[1] + tot_injected[1])/tot_porevol_init
                  << std::setw(width) << (init_polymass - polymass - tot_polyprod + tot_polyinj - polymass_adsorbed)
                  << std::endl;
        std::cout.precision(8);

        watercut.push(simtimer.currentTime() + simtimer.currentStepLength(),
                      produced[0]/(produced[0] + produced[1]),
                      tot_produced[0]/tot_porevol_init);
        if (wells->c_wells()) {
            wellreport.push(*props, *wells->c_wells(), state.saturation(),
                            simtimer.currentTime() + simtimer.currentStepLength(),
                            well_bhp, well_perfrates);
        }
    }
    total_timer.stop();

    std::cout << "\n\n================    End of simulation     ===============\n"
              << "Total time taken: " << total_timer.secsSinceStart()
              << "\n  Pressure time:  " << ptime
              << "\n  Transport time: " << ttime << std::endl;

    if (output) {
        outputState(*grid->c_grid(), state, simtimer.currentStepNum(), output_dir);
        outputWaterCut(watercut, output_dir);
        if (wells->c_wells()) {
            outputWellReport(wellreport, output_dir);
        }
    }
}
