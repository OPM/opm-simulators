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
/// The Darcy law gives 
///   \f[u_\alpha= -\frac1{\mu_\alpha} K_\alpha\nabla p_\alpha\f]
/// where \f$\mu_\alpha\f$ and \f$K_\alpha\f$ represent the viscosity
/// and the permeability tensor for each phase \f$\alpha\f$. In the two phase
/// case, we have either \f$\alpha=w\f$ or \f$\alpha=o\f$. 
/// In this tutorial, we do not take into account capillary pressure so that 
/// \f$p=p_w=p_o\f$ and gravity
/// effects. We  denote by \f$K\f$ the absolute permeability tensor and each phase
/// permeability is defined through its relative permeability by the expression
/// \f[K_\alpha=k_{r\alpha}K.\f]
/// The phase mobility are defined as
///  \f[\lambda_\alpha=\frac{k_{r\alpha}}{\mu_\alpha}\f]
/// so that the Darcy law may be rewritten as
///  \f[u_\alpha= -\lambda_\alpha K\nabla p.\f]
/// The conservation of mass for each phase writes:
/// \f[\frac{\partial}{\partial t}(\phi\rho_\alpha s_\alpha)+\nabla\cdot(\rho_\alpha u_\alpha)=q_\alpha\f]
/// where \f$s_\alpha\f$ denotes the saturation of the phase \f$\alpha\f$ and \f$q_\alpha\f$ is a source term. Let
/// us consider a two phase flow with oil and water. We assume that the phases are incompressible. Since
/// \f$s_w+s_o=1\f$, we get
///  \f[ \nabla\cdot u=\frac{q_w}{\rho_w}+\frac{q_o}{\rho_o}.\f]
/// where we define
///   \f[u=u_w+u_o.\f]
/// Let the total mobility be equal to
/// \f[\lambda=\lambda_w+\lambda_o\f]
/// Then, we have 
/// \f[u=-\lambda K\nabla p.\f]
/// The set of equations
/// \f[\nabla\cdot u=\frac{q_w}{\rho_w}+\frac{q_o}{\rho_o},\quad u=-\lambda K\nabla p.\f]
/// is referred to as the <strong>pressure equation</strong>. We introduce 
/// the fractional flow \f$f_w\f$
/// as
/// \f[f_w=\frac{\lambda_w}{\lambda_w+\lambda_o}\f]
/// and obtain
/// \f[\phi\frac{\partial s}{\partial t}+\nabla\cdot(f_w u)=\frac{q_w}{\rho_w}\f]
/// which is referred to as the <strong>transport equation</strong>. The pressure and 
/// transport equation are coupled. In this tutorial, we implement a splitting scheme, 
/// where, at each time step, we decouple the two equations. We solve first
/// the pressure equation and then update the water saturation by solving
/// the transport equation assuming that \f$u\f$ is constant in time in the time step 
/// interval we are considering.

/// \page tutorial3
/// \section commentedcode3 Commented code:
/// \page tutorial3
/// \details
/// We define a class which computes mobility, capillary pressure and
/// the minimum and maximum saturation value for each cell. 
/// \code
class TwophaseFluid
{
public:
    /// \endcode
    /// \page tutorial3
    /// \details Constructor operator. Takes in the fluid properties defined  
    /// \c props
    /// \code
    TwophaseFluid(const Opm::IncompPropertiesInterface& props);
    /// \endcode

    /// \page tutorial3
    /// \details Density for each phase.
    /// \code
    double density(int phase) const;
    /// \endcode

    /// \page tutorial3
    /// \details Computes the mobility and its derivative with respect to saturation
    /// for a given cell \c c and saturation \c sat. The template classes \c Sat, 
    /// \c Mob, \c DMob are typically arrays. By using templates, we do not have to 
    /// investigate how these array objects are implemented 
    /// (as long as they have an \c operator[] method).
    /// \code
    template <class Sat, class Mob, class DMob>  
    void mobility(int c, const Sat& s, Mob& mob, DMob& dmob) const;
    /// \endcode

    /// \page tutorial3
    /// \details Computes the capillary pressure  and its derivative with respect 
    ///  to saturation
    /// for a given cell \c c and saturation \c sat.
    /// \code
    template <class Sat, class Pcap, class DPcap>  
    void pc(int c, const Sat& s, Pcap& pcap, DPcap& dpcap) const;
    /// \endcode

    /// \page tutorial3
    /// \details Returns the minimum permitted saturation.
    /// \code
    double s_min(int c) const;
    /// \endcode

    /// \page tutorial3
    /// \details Returns the maximum permitted saturation
    /// \code
    double s_max(int c) const;
    /// \endcode

    /// \page tutorial3
    /// \details Private variables
    /// \code
private:
    const Opm::IncompPropertiesInterface& props_;
    std::vector<double> smin_;
    std::vector<double> smax_;
};
/// \endcode

/// \page tutorial3
/// \details We set up the transport model.
/// \code
typedef Opm::SinglePointUpwindTwoPhase<TwophaseFluid> TransportModel;
/// \endcode

/// \page tutorial3
/// \details
/// The transport equation is nonlinear. We use an implicit transport solver
/// which implements a Newton-Raphson solver. 
/// We define the format of the objects 
/// which will be used by the solver.
/// \code
using namespace Opm::ImplicitTransportDefault;
typedef NewtonVectorCollection< ::std::vector<double> >  NVecColl;
typedef JacobianSystem< struct CSRMatrix, NVecColl > JacSys;

template <class Vector>
class MaxNorm {
public:
    static double
    norm(const Vector& v) {
        return AccumulationNorm <Vector, MaxAbs>::norm(v);
    }
};
/// \endcode

/// \page tutorial3
/// \details
/// We set up the solver.
/// \code
typedef Opm::ImplicitTransport<TransportModel,
        JacSys        ,
        MaxNorm       ,
        VectorNegater ,
        VectorZero    ,
        MatrixZero    ,
        VectorAssign  > TransportSolver;
/// \endcode


/// \page tutorial3
/// \details
/// Main function
/// \code
int main ()
{
    /// \endcode
    /// \page tutorial3
    /// \details
    /// We define the grid. A cartesian grid with 1200 cells.
    /// \code
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
    /// \endcode

    /// \page tutorial3
    /// \details
    /// We define the properties of the fluid.\n 
    /// Number of phases.
    /// \code
    int num_phases = 2;
    using namespace unit;
    using namespace prefix;
    /// \endcode
    /// \page tutorial3
    /// \details density vector (one component per phase).
    /// \code
    std::vector<double> rho(2, 1000.);
    /// \endcode
    /// \page tutorial3
    /// \details viscosity vector (one component per phase).
    /// \code
    std::vector<double> mu(2, 1.*centi*Poise);
    /// \endcode
    /// \page tutorial3
    /// \details porosity and permeability of the rock.
    /// \code
    double porosity = 0.5;
    double k = 10*milli*darcy;
    /// \endcode

    /// \page tutorial3
    /// \details We define the relative permeability function. We use a basic fluid
    /// description and set this function to be linear. For more realistic fluid, the 
    /// saturation function is given by the data.
    /// \code
    SaturationPropsBasic::RelPermFunc rel_perm_func;
    rel_perm_func = SaturationPropsBasic::Linear;
    /// \endcode

    /// \page tutorial3
    /// \details We construct a basic fluid with the properties we have defined above.
    /// Each property is constant and hold for all cells.
    /// \code
    IncompPropertiesBasic props(num_phases, rel_perm_func, rho, mu,
                                     porosity, k, dim, num_cells);
    TwophaseFluid fluid(props);
    /// \endcode

    
    /// \page tutorial3
    /// \details Gravity parameters. Here, we set zero gravity.
    /// \code
    const double *grav = 0;
    std::vector<double> omega; 
    /// \endcode

    /// \page tutorial3
    /// \details We set up the pressure solver.
    /// \code
    LinearSolverUmfpack linsolver;
    IncompTpfa psolver(*grid.c_grid(), props.permeability(), grav, linsolver);
    /// \endcode

    

    /// \page tutorial3
    /// \details We set up the source term
    /// \code
    std::vector<double> src(num_cells, 0.0);
    src[0] = 1.;
    src[num_cells-1] = -1.;
    /// \endcode

    /// \page tutorial3
    /// \details We set up the wells. Here, there are no well and we let them empty.
    /// \code
    std::vector<double> empty_wdp; 
    std::vector<double> empty_well_bhp;
    std::vector<double> empty_well_flux;
    /// \endcode

    /// \page tutorial3
    /// \details We set up the source term for the transport solver.
    /// \code
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
    /// \endcode

    /// \page tutorial3
    /// \details We compute the pore volume
    /// \code
    std::vector<double> porevol;
    computePorevolume(*grid.c_grid(), props, porevol);
    /// \endcode

    /// \page tutorial3
    /// \details We set up the transport solver.
    /// \code
    const bool guess_old_solution = true;
    TransportModel  model  (fluid, *grid.c_grid(), porevol, grav, guess_old_solution);
    TransportSolver tsolver(model);
    /// \endcode

    /// \page tutorial3
    /// \details Time integration parameters
    /// \code
    double dt = 0.1*day;
    int num_time_steps = 20;
    /// \endcode

    /// \page tutorial3
    /// \details Control paramaters for the implicit solver.
    /// \code 
    ImplicitTransportDetails::NRReport  rpt;
    ImplicitTransportDetails::NRControl ctrl;
    /// \endcode



    /// \page tutorial3
    /// \details We define a vector which contains all cell indexes. We use this 
    /// vector to set up parameters on the whole domains.
    std::vector<int> allcells(num_cells);
    for (int cell = 0; cell < num_cells; ++cell) {
        allcells[cell] = cell;
    }


    /// \page tutorial3
    /// \details We set up the boundary conditions. Letting bcs empty is equivalent 
    /// to no flow boundary conditions.
    /// \code
    FlowBCManager bcs;
    /// \endcode

    /// \page tutorial3
    /// \details 
    /// Linear solver init.
    /// \code
    using ImplicitTransportLinAlgSupport::CSRMatrixUmfpackSolver;
    CSRMatrixUmfpackSolver linsolve;
    /// \endcode

    /// \page tutorial3
    /// \details 
    /// We initialise water saturation to minimum everywhere.
    /// \code
    TwophaseState state;
    state.init(*grid.c_grid());
    state.setWaterSat(allcells, props, TwophaseState::MinSat);
    /// \endcode
    
    /// \page tutorial3
    /// \details We introduce a vector which contains the total mobility 
    /// on all cells.
    /// \code
    std::vector<double> totmob;
    /// \endcode

    /// \page tutorial3
    /// \details This string will contain the name of a VTK output vector.
    /// \code
    std::ostringstream vtkfilename;
    /// \endcode


    /// \page tutorial3
    /// \details Loop over the time steps.
    /// \code
    for (int i = 0; i < num_time_steps; ++i) {
        /// \endcode
        /// \page tutorial3
        /// \details Compute the total mobility. It is needed by the pressure solver
        /// \code
        computeTotalMobility(props, allcells, state.saturation(), totmob);
        /// \endcode
        /// \page tutorial3
        /// \details Solve the pressure equation
        /// \code
        psolver.solve(totmob, omega, src, empty_wdp, bcs.c_bcs(), 
                      state.pressure(), state.faceflux(), empty_well_bhp, 
                      empty_well_flux);
        /// \endcode
        /// \page tutorial3
        /// \details  Transport solver
        /// \code
        tsolver.solve(*grid.c_grid(), tsrc, dt, ctrl, state, linsolve, rpt);
        /// \endcode

        /// \page tutorial3
        /// \details Write the output to file.
        /// \code
        vtkfilename.str(""); 
        vtkfilename << "tutorial3-" << std::setw(3) << std::setfill('0') << i << ".vtu";
        std::ofstream vtkfile(vtkfilename.str().c_str());
        Opm::DataMap dm;
        dm["saturation"] = &state.saturation();
        dm["pressure"] = &state.pressure();
        Opm::writeVtkData(*grid.c_grid(), dm, vtkfile);
    }
}
/// \endcode

/// \page tutorial3
/// \details Implementation of the TwophaseFluid class.
/// \code
TwophaseFluid::TwophaseFluid(const Opm::IncompPropertiesInterface& props)
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

double TwophaseFluid::density(int phase) const
{
    return props_.density()[phase];
}

template <class Sat,
          class Mob,
          class DMob>
void TwophaseFluid::mobility(int c, const Sat& s, Mob& mob, DMob& dmob) const
{
    props_.relperm(1, &s[0], &c, &mob[0], &dmob[0]);
    const double* mu = props_.viscosity();
    mob[0] /= mu[0];
    mob[1] /= mu[1];
    /// \endcode
    /// \page tutorial3
    /// \details We
    /// recall that we use Fortran ordering for kr derivatives,
    /// therefore dmob[i*2 + j] is row j and column i of the
    /// matrix.
    /// Each row corresponds to a kr function, so which mu to
    /// divide by also depends on the row, j.
    /// \code
    dmob[0*2 + 0] /= mu[0];
    dmob[0*2 + 1] /= mu[1];
    dmob[1*2 + 0] /= mu[0];
    dmob[1*2 + 1] /= mu[1];
}

template <class Sat,
          class Pcap,
          class DPcap>
void TwophaseFluid::pc(int c, const Sat& s, Pcap& pcap, DPcap& dpcap) const
{
    pcap = 0.;
    dpcap = 0.;
}

double TwophaseFluid::s_min(int c) const
{
    return smin_[2*c + 0];
}

double TwophaseFluid::s_max(int c) const
{
    return smax_[2*c + 0];
}
/// \endcode


/// \page tutorial3
/// \section results3 Results.
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
/// \details
/// \section completecode3 Complete source code:
/// \include tutorial3.cpp 
