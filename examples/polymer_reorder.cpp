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

#include <opm/polymer/IncompTpfaPolymer.hpp>
#include <opm/core/pressure/FlowBCManager.hpp>

#include <opm/core/grid.h>
#include <opm/core/GridManager.hpp>
#include <opm/core/newwells.h>
#include <opm/core/wells/WellsManager.hpp>
#include <opm/core/utility/ErrorMacros.hpp>
#include <opm/core/utility/initState.hpp>
#include <opm/core/simulator/SimulatorTimer.hpp>
#include <opm/core/utility/StopWatch.hpp>
#include <opm/core/utility/Units.hpp>
#include <opm/core/utility/writeVtkData.hpp>
#include <opm/core/utility/miscUtilities.hpp>
#include <opm/core/utility/parameters/ParameterGroup.hpp>

#include <opm/core/fluid/IncompPropertiesBasic.hpp>
#include <opm/core/fluid/IncompPropertiesFromDeck.hpp>
#include <opm/core/fluid/RockCompressibility.hpp>

#include <opm/core/linalg/LinearSolverFactory.hpp>
//#include <opm/core/linalg/LinearSolverAGMG.hpp>

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
#include <opm/core/simulator/WellState.hpp>
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
#include <list>



static void outputState(const UnstructuredGrid& grid,
                        const Opm::PolymerState& state,
                        const int step,
                        const std::string& output_dir,
                        const Opm::TransportModelPolymer& reorder_model)
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
    dm["faceflux"] = &state.faceflux();
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

    std::ostringstream fname;
    fname << output_dir << "/" << "residualcounts" << "-" << std::setw(3) << std::setfill('0') << step << ".dat";
    std::ofstream file(fname.str().c_str());
    if (!file) {
        THROW("Failed to open " << fname.str());
    }

    typedef std::list<Opm::TransportModelPolymer::Newton_Iter> ListRes;

    const ListRes& res_counts = reorder_model.res_counts;
    for (ListRes::const_iterator it = res_counts.begin(); it != res_counts.end(); ++it) {
        file << it->res_s << "," << it->cell << "," << std::setprecision(15) << it->s << "," << std::setprecision(15) << it->c << "\n";
    }
    file.close();

    
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

class IncompPropertiesCorey : public Opm::IncompPropertiesBasic {

private:
    std::vector<double> exponents_;
    int np_;

    double corey_kr(double s, int p) const  {
        return std::pow(s, exponents_[p]);
    }

    double corey_dkrds(double s, int p) const {
        return exponents_[p]*std::pow(s, exponents_[p] - 1.0);
    }

public:
    IncompPropertiesCorey(const Opm::parameter::ParameterGroup& param,
                          const int dim,
                          const int num_cells,
                          const std::vector<double> exponents
                          ) : IncompPropertiesBasic(param, dim, num_cells) {
        exponents_ = exponents;
        np_ = numPhases();
    }
    /// \param[in]  n      Number of data points.
    /// \param[in]  s      Array of nP saturation values.
    /// \param[in]  cells  Array of n cell indices to be associated with the s values.
    /// \param[out] kr     Array of nP relperm values, array must be valid before calling.
    /// \param[out] dkrds  If non-null: array of nP^2 relperm derivative values,
    ///                    array must be valid before calling.
    ///                    The P^2 derivative matrix is
    ///                           m_{ij} = \frac{dkr_i}{ds^j},
    ///                    and is output in Fortran order (m_00 m_10 m_20 m_01 ...)
    virtual void relperm(const int n,
                         const double* s,
                         const int* /*cells*/,
                         double* kr,
                         double* dkrds) const {

        if (dkrds == 0) {
            // #pragma omp parallel for
            for (int i = 0; i < n; ++i) {
                for (int p = 0; p < np_; ++p) {
                    kr[i*np_ + p] = corey_kr(s[i*np_ + p], p);
                }
            }
            return;
        }
        // #pragma omp parallel for
        for (int i = 0; i < n; ++i) {
            std::fill(dkrds + i*np_*np_, dkrds + (i+1)*np_*np_, 0.0);
            for (int p = 0; p < np_; ++p) {
                kr[i*np_ + p] = corey_kr(s[i*np_ + p], p);
                // Only diagonal elements in derivative.
                dkrds[i*np_*np_ + p*np_ + p] = corey_dkrds(s[i*np_ + p], p);
            }
        }
    }
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
    bool use_deck = param.getDefault("use_deck", true);
    use_deck = param.has("deck_filename") && use_deck;
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
        props.reset(new Opm::IncompPropertiesFromDeck(deck, *grid->c_grid()));
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
        bool use_corey = false;
        use_corey = param.getDefault("use_corey", false);
        if (use_corey) {
            std::vector<double> exponents(2, 1.0);
            exponents[0] = param.getDefault("n1", 1.0);
            exponents[1] = param.getDefault("n2", 1.0);
            props.reset(new IncompPropertiesCorey(param, grid->c_grid()->dimensions, grid->c_grid()->number_of_cells, exponents));
        } else {
            props.reset(new IncompPropertiesBasic(param, grid->c_grid()->dimensions, grid->c_grid()->number_of_cells));
        }
        // Wells init.
        wells.reset(new Opm::WellsManager());
        // Timer init.
        simtimer.init(param);
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
	bool use_deck_fluid = param.getDefault("use_deck_fluid", false);
	if(!use_deck_fluid){
	    // Rock compressibility.
	    rock_comp.reset(new Opm::RockCompressibility(param));
	    double c_max = param.getDefault("c_max_limit", 5.0);
	    double mix_param = param.getDefault("mix_param", 1.0);
	    double rock_density = param.getDefault("rock_density", 1000.0);
	    double dead_pore_vol = param.getDefault("dead_pore_vol", 0.1);
	    double res_factor = param.getDefault("res_factor", 1.) ; // res_factor = 1 gives no change in permeability
	    double c_max_ads = param.getDefault("c_max_ads", 1.);
	    int ads_index = param.getDefault<int>("ads_index", Opm::PolymerProperties::NoDesorption);
	    std::vector<double> c_vals_visc(2, -1e100);
	    c_vals_visc[0] = 0.0;
	    c_vals_visc[1] = c_max;
	    std::vector<double> visc_mult_vals(2, -1e100);
	    visc_mult_vals[0] = 1.0;
	    visc_mult_vals[1] = param.getDefault("c_max_viscmult", 30.0);
	    std::vector<double> c_vals_ads(2, -1e100);
	    c_vals_ads[0] = 0.0;
	    c_vals_ads[1] = 8.0;
	    // Here we set up adsorption equal to zero.
	    std::vector<double> ads_vals(2, -1e100);
	    ads_vals[0] = 0.0;
	    ads_vals[1] = 0.0;
	    polyprop.set(c_max, mix_param, rock_density, dead_pore_vol, res_factor, c_max_ads,
			 static_cast<Opm::PolymerProperties::AdsorptionBehaviour>(ads_index),
			 c_vals_visc,  visc_mult_vals, c_vals_ads, ads_vals);
	}else{
	    std::string deck_filename = param.get<std::string>("deck_filename");
	    Opm::EclipseGridParser deck(deck_filename);
	    rock_comp.reset(new Opm::RockCompressibility(deck));
	    polyprop.readFromDeck(deck);
	}
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
    double nl_pressure_residual_tolerance = 0.0;
    double nl_pressure_change_tolerance = 0.0;
    if (rock_comp->isActive()) {
        if (!use_reorder) {
            THROW("Cannot run implicit (non-reordering) transport solver with rock compressibility yet.");
        }
        nl_pressure_residual_tolerance = param.getDefault("nl_pressure_residual_tolerance", 0.0);
        nl_pressure_change_tolerance = param.getDefault("nl_pressure_change_tolerance", 1.0); // In Pascal.
        nl_pressure_maxiter = param.getDefault("nl_pressure_maxiter", 10);
    }

    // Source-related variables init.
    int num_cells = grid->c_grid()->number_of_cells;

    // Extra rock init.
    std::vector<double> porevol;
    if (rock_comp->isActive()) {
        computePorevolume(*grid->c_grid(), props->porosity(), *rock_comp, state.pressure(), porevol);
    } else {
        computePorevolume(*grid->c_grid(), props->porosity(), porevol);
    }
    double tot_porevol_init = std::accumulate(porevol.begin(), porevol.end(), 0.0);

    // Initialising src
    std::vector<double> src(num_cells, 0.0);
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
    //Opm::LinearSolverAGMG linsolver;
    // Pressure solver.
    const double *grav = use_gravity ? &gravity[0] : 0;
    Opm::IncompTpfaPolymer psolver(*grid->c_grid(), *props, rock_comp.get(), polyprop, linsolver,
                                   nl_pressure_residual_tolerance, nl_pressure_change_tolerance,
                                   nl_pressure_maxiter,
                                   grav, wells->c_wells(), src, bcs.c_bcs());
    // Reordering solver.
    const double nl_tolerance = param.getDefault("nl_tolerance", 1e-9);
    const int nl_maxiter = param.getDefault("nl_maxiter", 30);
    Opm::TransportModelPolymer::SingleCellMethod method;
    std::string method_string = param.getDefault("single_cell_method", std::string("Bracketing"));
    if (method_string == "Bracketing") {
        method = Opm::TransportModelPolymer::Bracketing;
    } else if (method_string == "Newton") {
        method = Opm::TransportModelPolymer::Newton;
    } else if (method_string == "Gradient") {
        method = Opm::TransportModelPolymer::Gradient;
    } else if (method_string == "NewtonSimpleSC") {
        method = Opm::TransportModelPolymer::NewtonSimpleSC;
    } else if (method_string == "NewtonSimpleC") {
        method = Opm::TransportModelPolymer::NewtonSimpleC;
    } else {
        THROW("Unknown method: " << method_string);
    }

    Opm::TransportModelPolymer reorder_model(*grid->c_grid(), *props, polyprop,
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
    Opm::WellState well_state;
    well_state.init(wells->c_wells(), state);
    std::vector<double> fractional_flows;
    std::vector<double> well_resflows_phase;
    int num_wells = 0;
    if (wells->c_wells()) {
        num_wells = wells->c_wells()->number_of_wells;
        well_resflows_phase.resize((wells->c_wells()->number_of_phases)*(wells->c_wells()->number_of_wells), 0.0);
        wellreport.push(*props, *wells->c_wells(), state.saturation(), 0.0, well_state.bhp(), well_state.perfRates());
    }
    for (; !simtimer.done(); ++simtimer) {
        // Report timestep and (optionally) write state to disk.
        simtimer.report(std::cout);
        if (output && (simtimer.currentStepNum() % output_interval == 0)) {
            outputState(*grid->c_grid(), state, simtimer.currentStepNum(), output_dir, reorder_model);
        }

        // Solve pressure.
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
            std::vector<double> initial_pressure = state.pressure();
            psolver.solve(simtimer.currentStepLength(), state, well_state);
            if (!rock_comp->isActive()) {
                // Compute average pressures of previous and last
                // step, and total volume.
                double av_prev_press = 0.;
                double av_press = 0.;
                double tot_vol = 0.;
                for (int cell = 0; cell < num_cells; ++cell) {
                    av_prev_press += initial_pressure[cell]*grid->c_grid()->cell_volumes[cell];
                    av_press      += state.pressure()[cell]*grid->c_grid()->cell_volumes[cell];
                    tot_vol       += grid->c_grid()->cell_volumes[cell];
                }
                // Renormalization constant
                const double ren_const = (av_prev_press - av_press)/tot_vol;
                for (int cell = 0; cell < num_cells; ++cell) {
                    state.pressure()[cell] += ren_const;
                }
                for (int well = 0; well < num_wells; ++well) {
                    well_state.bhp()[well] += ren_const;
                }
            }
            pressure_timer.stop();
            double pt = pressure_timer.secsSinceStart();
            std::cout << "Pressure solver took:  " << pt << " seconds." << std::endl;
            ptime += pt;


            if (check_well_controls) {
                Opm::computePhaseFlowRatesPerWell(*wells->c_wells(),
                                                  fractional_flows,
                                                  well_state.perfRates(),
                                                  well_resflows_phase);
                std::cout << "Checking well conditions." << std::endl;
                // For testing we set surface := reservoir
                well_control_passed = wells->conditionsMet(well_state.bhp(), well_resflows_phase, well_resflows_phase);
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

        // Update pore volumes if rock is compressible.
        if (rock_comp->isActive()) {
            computePorevolume(*grid->c_grid(), props->porosity(), *rock_comp, state.pressure(), porevol);
        }

        // Process transport sources (to include bdy terms and well flows).
        Opm::computeTransportSource(*grid->c_grid(), src, state.faceflux(), 1.0,
                                    wells->c_wells(), well_state.perfRates(), reorder_src);


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
                reorder_model.solve(&state.faceflux()[0], &porevol[0], &reorder_src[0], stepsize, inflow_c,
                                    state.saturation(), state.concentration(), state.maxconcentration());
                Opm::computeInjectedProduced(*props, polyprop, state.saturation(), state.concentration(), state.maxconcentration(),
                                             reorder_src, simtimer.currentStepLength(), inflow_c,
                                             injected, produced, polyinj, polyprod);
                if (use_segregation_split) {
                    if (use_column_solver) {
                        if (use_gauss_seidel_gravity) {
                            reorder_model.solveGravity(columns, &porevol[0], stepsize, state.saturation(),
                                                       state.concentration(), state.maxconcentration());
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
                            well_state.bhp(), well_state.perfRates());
        }
    }
    total_timer.stop();

    std::cout << "\n\n================    End of simulation     ===============\n"
              << "Total time taken: " << total_timer.secsSinceStart()
              << "\n  Pressure time:  " << ptime
              << "\n  Transport time: " << ttime << std::endl;

    if (output) {
        outputState(*grid->c_grid(), state, simtimer.currentStepNum(), output_dir, reorder_model);
        outputWaterCut(watercut, output_dir);
        if (wells->c_wells()) {
            outputWellReport(wellreport, output_dir);
        }
    }
}
