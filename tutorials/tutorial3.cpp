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


#if HAVE_CONFIG_H
#include "config.h"
#endif // HAVE_CONFIG_H

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <cassert>
#include <opm/core/grid.h>
#include <opm/core/GridManager.hpp>
#include <opm/core/utility/writeVtkData.hpp>
#include <opm/core/linalg/LinearSolverUmfpack.hpp>
#include <opm/core/pressure/IncompTpfa.hpp>
#include <opm/core/pressure/FlowBCManager.hpp>
#include <opm/core/fluid/IncompPropertiesBasic.hpp>

#include <opm/core/transport/transport_source.h>
#include <opm/core/transport/CSRMatrixUmfpackSolver.hpp>
#include <opm/core/transport/NormSupport.hpp>
#include <opm/core/transport/ImplicitAssembly.hpp>
#include <opm/core/transport/ImplicitTransport.hpp>
#include <opm/core/transport/JacobianSystem.hpp>
#include <opm/core/transport/CSRMatrixBlockAssembler.hpp>
#include <opm/core/transport/SinglePointUpwindTwoPhase.hpp>

#include <opm/core/TwophaseState.hpp>

#include <opm/core/utility/miscUtilities.hpp>
#include <opm/core/utility/Units.hpp>
#include <opm/core/utility/parameters/ParameterGroup.hpp>

/// \page tutorial3 Multiphase flow
/// Multiphase flow

class TwophaseFluid
{
public:
    TwophaseFluid(const Opm::IncompPropertiesInterface& props)
    : props_(props),
      smin_(props.numCells()*props.numPhases()),
      smax_(props.numCells()*props.numPhases())
    {
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

    template <class Sat,
    class Mob,
    class DMob>
    void mobility(int c, const Sat& s, Mob& mob, DMob& dmob) const
    {
        props_.relperm(1, &s[0], &c, &mob[0], &dmob[0]);
        const double* mu = props_.viscosity();
        mob[0] /= mu[0];
        mob[1] /= mu[1];
        // Recall that we use Fortran ordering for kr derivatives,
        // therefore dmob[i*2 + j] is row j and column i of the
        // matrix.
        // Each row corresponds to a kr function, so which mu to
        // divide by also depends on the row, j.
        dmob[0*2 + 0] /= mu[0];
        dmob[0*2 + 1] /= mu[1];
        dmob[1*2 + 0] /= mu[0];
        dmob[1*2 + 1] /= mu[1];
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
        dpcap = dpcow[0];
    }

    double s_min(int c) const
    {
        return smin_[2*c + 0];
    }

    double s_max(int c) const
    {
        return smax_[2*c + 0];
    }

private:
    const Opm::IncompPropertiesInterface& props_;
    std::vector<double> smin_;
    std::vector<double> smax_;
};

typedef Opm::SinglePointUpwindTwoPhase<TwophaseFluid> TransportModel;

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

typedef Opm::ImplicitTransport<TransportModel,
        JacSys        ,
        MaxNorm       ,
        VectorNegater ,
        VectorZero    ,
        MatrixZero    ,
        VectorAssign  > TransportSolver;
 


int main ()
{

    int dim = 3;
    int nx = 20;
    int ny = 20;
    int nz = 1;
    double dx = 10.;
    double dy = 10.;
    double dz = 10.;
    using namespace Opm;
    GridManager grid(nx, ny, nz, dx, dy, dz);
    int num_cells = grid.c_grid()->number_of_cells;
    int num_faces = grid.c_grid()->number_of_faces;

    int num_phases = 2;
    using namespace unit;
    using namespace prefix;
    std::vector<double> rho(2, 1000.);
    std::vector<double> mu(2, 1.*centi*Poise);
    double porosity = 0.5;
    double k = 10*milli*darcy;
        
    SaturationPropsBasic::RelPermFunc rel_perm_func;
    rel_perm_func = SaturationPropsBasic::Linear;

    IncompPropertiesBasic props(num_phases, rel_perm_func, rho, mu,
                                     porosity, k, dim, num_cells);



    const double *grav = 0;
    std::vector<double> omega; 

    LinearSolverUmfpack linsolver;
    IncompTpfa psolver(*grid.c_grid(), props.permeability(), grav, linsolver);


    TwophaseFluid fluid(props);

    std::vector<double> src(num_cells, 0.0);
    src[0] = 1.;
    src[num_cells-1] = -1.;


    std::vector<double> empty_wdp; 
    std::vector<double> empty_well_bhp;
    std::vector<double> empty_well_flux;

    TransportSource* tsrc = create_transport_source(2, 2);
    double ssrc[]   = { 1.0, 0.0 };
    double ssink[]  = { 0.0, 1.0 };
    double zdummy[] = { 0.0, 0.0 };
    for (int cell = 0; cell < num_cells; ++cell) {
        if (src[cell] > 0.0) {
            append_transport_source(cell, 2, 0, src[cell], ssrc, zdummy, tsrc);
        } else if (src[cell] < 0.0) {
            append_transport_source(cell, 2, 0, src[cell], ssink, zdummy, tsrc);
        }
    }
    
    std::vector<double> porevol;
    computePorevolume(*grid.c_grid(), props, porevol);

    const bool guess_old_solution = true;
    TransportModel  model  (fluid, *grid.c_grid(), porevol, grav, guess_old_solution);
    TransportSolver tsolver(model);

    double dt = 0.1*day;
    int num_time_steps = 20;

    TwophaseState state;
    ImplicitTransportDetails::NRReport  rpt;
    ImplicitTransportDetails::NRControl ctrl;
    std::vector<double> totmob;

    std::vector<int> allcells(num_cells);
    for (int cell = 0; cell < num_cells; ++cell) {
        allcells[cell] = cell;
    }

    FlowBCManager bcs;

    // Linear solver init.
    using ImplicitTransportLinAlgSupport::CSRMatrixUmfpackSolver;
    CSRMatrixUmfpackSolver linsolve;

    state.init(*grid.c_grid());
    // By default: initialise water saturation to minimum everywhere.
    state.setWaterSat(allcells, props, TwophaseState::MinSat);
    

    std::ostringstream vtkfilename;

    for (int i = 0; i < num_time_steps; ++i) {

        computeTotalMobility(props, allcells, state.saturation(), totmob);
        psolver.solve(totmob, omega, src, empty_wdp, bcs.c_bcs(), 
                      state.pressure(), state.faceflux(), empty_well_bhp, 
                      empty_well_flux);
        tsolver.solve(*grid.c_grid(), tsrc, dt, ctrl, state, linsolve, rpt);

        vtkfilename.str(""); 
        vtkfilename << "tutorial3-" << std::setw(3) << std::setfill('0') << i << ".vtu";
        std::ofstream vtkfile(vtkfilename.str().c_str());
        Opm::DataMap dm;
        dm["saturation"] = &state.saturation();
        dm["pressure"] = &state.pressure();
        Opm::writeVtkData(*grid.c_grid(), dm, vtkfile);

    }
}


/// \page tutorial3
/// <TABLE>
/// <TR>
/// <TD> \image html tutorial3-000.png </TD>
/// <TD> \image html tutorial3-005.png </TD>
/// </TR>
/// <TR>
/// <TD> \image html tutorial3-010.png </TD>
/// <TD> \image html tutorial3-015.png </TD>
/// </TR>
/// <TR>
/// <TD> \image html tutorial3-019.png </TD>
/// <TD> </TD>
/// </TR>
/// </TABLE>

/// \page tutorial3
/// \section sourcecode Complete source code.
/// \include tutorial3.cpp 
