/*===========================================================================
//
// File: spu_2p.cpp
//
// Created: 2011-10-05 10:29:01+0200
//
// Authors: Ingeborg S. Ligaarden <Ingeborg.Ligaarden@sintef.no>
//          Jostein R. Natvig     <Jostein.R.Natvig@sintef.no>
//          Halvor M. Nilsen      <HalvorMoll.Nilsen@sintef.no>
//          Atgeirr F. Rasmussen  <atgeirr@sintef.no>
//          Bård Skaflestad       <Bard.Skaflestad@sintef.no>
//
//==========================================================================*/


/*
  Copyright 2011 SINTEF ICT, Applied Mathematics.
  Copyright 2011 Statoil ASA.

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
#include "config.h"

#include <opm/core/linalg/sparse_sys.h>

#include <opm/core/pressure/tpfa/ifs_tpfa.h>
#include <opm/core/pressure/tpfa/trans_tpfa.h>

#include <opm/core/utility/cart_grid.h>
#include <opm/core/utility/ErrorMacros.hpp>
#include <opm/core/utility/Units.hpp>
#include <opm/core/utility/cpgpreprocess/cgridinterface.h>
#include <opm/core/utility/parameters/ParameterGroup.hpp>

#include <opm/core/fluid/SimpleFluid2p.hpp>
#include <opm/core/fluid/IncompPropertiesBasic.hpp>
#include <opm/core/fluid/IncompPropertiesFromDeck.hpp>

#include <opm/core/transport/transport_source.h>
#include <opm/core/transport/CSRMatrixUmfpackSolver.hpp>
#include <opm/core/transport/NormSupport.hpp>
#include <opm/core/transport/ImplicitAssembly.hpp>
#include <opm/core/transport/ImplicitTransport.hpp>
#include <opm/core/transport/JacobianSystem.hpp>
#include <opm/core/transport/CSRMatrixBlockAssembler.hpp>
#include <opm/core/transport/SinglePointUpwindTwoPhase.hpp>

#include <opm/core/transport/reorder/twophasetransport.hpp>

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




namespace Opm
{


    /// Concrete grid class constructing a
    /// corner point grid from a deck,
    /// or a cartesian grid.
    class Grid
    {
    public:
	Grid(const Opm::EclipseGridParser& deck)
	{
	    // Extract data from deck.
	    const std::vector<double>& zcorn = deck.getFloatingPointValue("ZCORN");
	    const std::vector<double>& coord = deck.getFloatingPointValue("COORD");
	    const std::vector<int>& actnum = deck.getIntegerValue("ACTNUM");
	    std::vector<int> dims;
	    if (deck.hasField("DIMENS")) {
		dims = deck.getIntegerValue("DIMENS");
	    } else if (deck.hasField("SPECGRID")) {
		dims = deck.getSPECGRID().dimensions;
	    } else {
		THROW("Deck must have either DIMENS or SPECGRID.");
	    }

	    // Collect in input struct for preprocessing.
	    struct grdecl grdecl;
	    grdecl.zcorn = &zcorn[0];
	    grdecl.coord = &coord[0];
	    grdecl.actnum = &actnum[0];
	    grdecl.dims[0] = dims[0];
	    grdecl.dims[1] = dims[1];
	    grdecl.dims[2] = dims[2];

	    // Process and compute.
	    ug_ = preprocess(&grdecl, 0.0);
	    compute_geometry(ug_);
	}

	Grid(int nx, int ny)
	{
	    ug_ = create_cart_grid_2d(nx, ny);
	}

	Grid(int nx, int ny, int nz)
	{
	    ug_ = create_cart_grid_3d(nx, ny, nz);
	}

	~Grid()
	{
	    free_grid(ug_);
	}

	virtual const UnstructuredGrid* c_grid() const
	{
	    return ug_;
	}

    private:
	// Disable copying and assignment.
	Grid(const Grid& other);
	Grid& operator=(const Grid& other);
	struct UnstructuredGrid* ug_;
    };

} // namespace Opm




class ReservoirState {
public:
    ReservoirState(const UnstructuredGrid* g, const int num_phases = 2)
        : press_ (g->number_of_cells, 0.0),
          fpress_(g->number_of_faces, 0.0),
          flux_  (g->number_of_faces, 0.0),
          sat_   (num_phases * g->number_of_cells, 0.0)
    {
	for (int cell = 0; cell < g->number_of_cells; ++cell) {
	    sat_[num_phases*cell + num_phases - 1] = 1.0;
	}
    }

    int numPhases() const { return sat_.size()/press_.size(); }

    ::std::vector<double>& pressure    () { return press_ ; }
    ::std::vector<double>& facepressure() { return fpress_; }
    ::std::vector<double>& faceflux    () { return flux_  ; }
    ::std::vector<double>& saturation  () { return sat_   ; }

    const ::std::vector<double>& pressure    () const { return press_ ; }
    const ::std::vector<double>& facepressure() const { return fpress_; }
    const ::std::vector<double>& faceflux    () const { return flux_  ; }
    const ::std::vector<double>& saturation  () const { return sat_   ; }

private:
    ::std::vector<double> press_ ;
    ::std::vector<double> fpress_;
    ::std::vector<double> flux_  ;
    ::std::vector<double> sat_   ;
};




class PressureSolver {
public:
    PressureSolver(const UnstructuredGrid* g,
		   const Opm::IncompPropertiesInterface& props)
        : htrans_(g->cell_facepos[ g->number_of_cells ]),
          trans_ (g->number_of_faces),
          gpress_(g->cell_facepos[ g->number_of_cells ])
    {
	UnstructuredGrid* gg = const_cast<UnstructuredGrid*>(g);
        tpfa_htrans_compute(gg, props.permeability(), &htrans_[0]);

        h_ = ifs_tpfa_construct(gg);
    }

    ~PressureSolver()
    {
        ifs_tpfa_destroy(h_);
    }

    template <class State>
    void
    solve(const UnstructuredGrid*      g     ,
          const ::std::vector<double>& totmob,
          const ::std::vector<double>& src   ,
          State&                       state )
    {
	UnstructuredGrid* gg = const_cast<UnstructuredGrid*>(g);
        tpfa_eff_trans_compute(gg, &totmob[0], &htrans_[0], &trans_[0]);

        // No gravity
        ::std::fill(gpress_.begin(), gpress_.end(), double(0.0));

        ifs_tpfa_assemble(gg, &trans_[0], &src[0], &gpress_[0], h_);

        using Opm::ImplicitTransportLinAlgSupport::CSRMatrixUmfpackSolver;

        CSRMatrixUmfpackSolver linsolve;
        linsolve.solve(h_->A, h_->b, h_->x);

        ifs_tpfa_press_flux(gg, &trans_[0], h_,
                            &state.pressure()[0],
                            &state.faceflux()[0]);
    }

private:
    ::std::vector<double> htrans_;
    ::std::vector<double> trans_ ;
    ::std::vector<double> gpress_;

    struct ifs_tpfa_data* h_;
};




static void
compute_porevolume(const UnstructuredGrid* g,
                   const Opm::IncompPropertiesInterface& props,
                   std::vector<double>& porevol)
{
    int num_cells = g->number_of_cells;
    porevol.resize(num_cells);
    const double* poro = props.porosity();
    ::std::transform(poro, poro + num_cells,
                     g->cell_volumes,
                     porevol.begin(),
                     ::std::multiplies<double>());
}


static void
compute_totmob(const Opm::IncompPropertiesInterface& props,
	       const std::vector<double>& s,
	       std::vector<double>& totmob)
{
    int num_cells = props.numCells();
    int num_phases = props.numPhases();
    totmob.resize(num_cells);
    ASSERT(int(s.size()) == num_cells*num_phases);
    std::vector<int> cells(num_cells);
    for (int cell = 0; cell < num_cells; ++cell) {
	cells[cell] = cell;
    }
    std::vector<double> kr(num_cells*num_phases);
    props.relperm(num_cells, &s[0], &cells[0], &kr[0], 0);
    const double* mu = props.viscosity();
    for (int cell = 0; cell < num_cells; ++cell) {
	totmob[cell] = 0;
	for (int phase = 0; phase < num_phases; ++phase) {	
	    totmob[cell] += kr[2*cell + phase]/mu[phase];
	}
    }
}


template <class State>
void outputState(const UnstructuredGrid* grid,
		 const State& state,
		 const int step,
		 const std::string& output_dir)
{
    std::ostringstream vtkfilename;
    vtkfilename << output_dir << "/output-" << std::setw(3) << std::setfill('0') << step << ".vtu";
    std::ofstream vtkfile(vtkfilename.str().c_str());
    if (!vtkfile) {
	THROW("Failed to open " << vtkfilename.str());
    }
    writeVtkDataGeneralGrid(grid, state, vtkfile);
}




template <class State>
void writeVtkDataAllCartesian(const std::tr1::array<int, 3>& dims,
			      const std::tr1::array<double, 3>& cell_size,
			      const State& state,
			      std::ostream& vtk_file)
{
    // Dimension is hardcoded in the prototype and the next two lines,
    // but the rest is flexible (allows dimension == 2 or 3).
    int dimension = 3;
    int num_cells = dims[0]*dims[1]*dims[2];

    ASSERT(dimension == 2 || dimension == 3);
    ASSERT(num_cells = dims[0]*dims[1]* (dimension == 2 ? 1 : dims[2]));

    vtk_file << "# vtk DataFile Version 2.0\n";
    vtk_file << "Structured Grid\n \n";
    vtk_file << "ASCII \n";
    vtk_file << "DATASET STRUCTURED_POINTS\n";

    vtk_file << "DIMENSIONS "
	     << dims[0] + 1 << " "
	     << dims[1] + 1 << " ";
    if (dimension == 3) {
	vtk_file << dims[2] + 1;
    } else {
	vtk_file << 1;
    }
    vtk_file << "\n";
	
    vtk_file << "ORIGIN " << 0.0 << " " << 0.0 << " " << 0.0 << "\n";

    vtk_file << "SPACING " << cell_size[0] << " " << cell_size[1];
    if (dimension == 3) {
	vtk_file << " " << cell_size[2];
    } else {
	vtk_file << " " << 0.0;
    }
    vtk_file << "\n";

    vtk_file << "CELL_DATA " << num_cells << '\n';
    vtk_file << "SCALARS pressure float" << '\n';
    vtk_file << "LOOKUP_TABLE pressure_table " << '\n';
    for (int i = 0; i < num_cells; ++i) {
	vtk_file << state.pressure()[i] << '\n';
    }

    ASSERT(state.numPhases() == 2);
    vtk_file << "SCALARS saturation float" << '\n';
    vtk_file << "LOOKUP_TABLE saturation_table " << '\n';
    for (int i = 0; i < num_cells; ++i) {
	double s = state.saturation()[2*i];
    	if (s > 1e-10) {
    	    vtk_file << s << '\n';
    	} else {
    	    vtk_file << 0.0 << '\n';
    	}
    }
}

typedef std::map<std::string, std::string> PMap;


struct Tag
{
    Tag(const std::string& tag, const PMap& props, std::ostream& os)
	: name_(tag), os_(os)
    {
	indent(os);
	os << "<" << tag;
	for (PMap::const_iterator it = props.begin(); it != props.end(); ++it) {
	    os << " " << it->first << "=\"" << it->second << "\""; 
	}
	os << ">\n";
	++indent_;
    }
    Tag(const std::string& tag, std::ostream& os)
	: name_(tag), os_(os)
    {
	indent(os);
	os << "<" << tag << ">\n";
	++indent_;
    }
    ~Tag()
    {
	--indent_;
	indent(os_);
	os_ << "</" << name_ << ">\n";
    }
    static void indent(std::ostream& os)
    {
	for (int i = 0; i < indent_; ++i) {
	    os << "  ";
	}
    }
private:
    static int indent_;
    std::string name_;
    std::ostream& os_;
};

int Tag::indent_ = 0;


template <class State>
void writeVtkDataGeneralGrid(const UnstructuredGrid* grid,
			     const State& state,
			     std::ostream& os)
{
    if (grid->dimensions != 3) {
	THROW("Vtk output for 3d grids only");
    }
    os.precision(12);
    os << "<?xml version=\"1.0\"?>\n";
    PMap pm;
    pm["type"] = "UnstructuredGrid";
    Tag vtkfiletag("VTKFile", pm, os);
    Tag ugtag("UnstructuredGrid", os);
    int num_pts = grid->number_of_nodes;
    int num_cells = grid->number_of_cells;
    pm.clear();
    pm["NumberOfPoints"] = boost::lexical_cast<std::string>(num_pts);
    pm["NumberOfCells"] = boost::lexical_cast<std::string>(num_cells);
    Tag piecetag("Piece", pm, os);
    {
	Tag pointstag("Points", os);
	pm.clear();
	pm["type"] = "Float64";
	pm["Name"] = "Coordinates";
	pm["NumberOfComponents"] = "3";
	pm["format"] = "ascii";
	Tag datag("DataArray", pm, os);
	for (int i = 0; i < num_pts; ++i) {
	    Tag::indent(os);
	    os << grid->node_coordinates[3*i + 0] << ' '
	       << grid->node_coordinates[3*i + 1] << ' '
	       << grid->node_coordinates[3*i + 2] << '\n';
	}
    }
    {
	Tag cellstag("Cells", os);
	pm.clear();
	pm["type"] = "Int32";
	pm["NumberOfComponents"] = "1";
	pm["format"] = "ascii";
	std::vector<int> cell_numpts;
	cell_numpts.reserve(num_cells);
	{
	    pm["Name"] = "connectivity";
	    Tag t("DataArray", pm, os);
	    int hf = 0;
	    for (int c = 0; c < num_cells; ++c) {
		std::set<int> cell_pts;
		for (; hf < grid->cell_facepos[c+1]; ++hf) {
		    int f = grid->cell_faces[hf];
		    const int* fnbeg = grid->face_nodes + grid->face_nodepos[f];
		    const int* fnend = grid->face_nodes + grid->face_nodepos[f+1];
		    cell_pts.insert(fnbeg, fnend);
		}
		cell_numpts.push_back(cell_pts.size());
		Tag::indent(os);
		std::copy(cell_pts.begin(), cell_pts.end(),
			  std::ostream_iterator<int>(os, " "));
		os << '\n';
	    }
	}
	{
	    pm["Name"] = "offsets";
	    Tag t("DataArray", pm, os);
	    int offset = 0;
	    const int num_per_line = 10;
	    for (int c = 0; c < num_cells; ++c) {
		if (c % num_per_line == 0) {
		    Tag::indent(os);
		}
		offset += cell_numpts[c];
		os << offset << ' ';
		if (c % num_per_line == num_per_line - 1
		    || c == num_cells - 1) {
		    os << '\n';
		}
	    }
	}
	std::vector<int> cell_foffsets;
	cell_foffsets.reserve(num_cells);
	{
	    pm["Name"] = "faces";
	    Tag t("DataArray", pm, os);
	    const int* fp = grid->cell_facepos;
	    int offset = 0;
	    for (int c = 0; c < num_cells; ++c) {
		Tag::indent(os);
		os << fp[c+1] - fp[c] << '\n';
		++offset;
		for (int hf = fp[c]; hf < fp[c+1]; ++hf) {
		    int f = grid->cell_faces[hf];
		    const int* np = grid->face_nodepos;
		    int f_num_pts = np[f+1] - np[f];
		    Tag::indent(os);
		    os << f_num_pts << ' ';
		    ++offset;
		    std::copy(grid->face_nodes + np[f],
			      grid->face_nodes + np[f+1],
			      std::ostream_iterator<int>(os, " "));
		    os << '\n';
		    offset += f_num_pts;
		}
		cell_foffsets.push_back(offset);
	    }
	}
	{
	    pm["Name"] = "faceoffsets";
	    Tag t("DataArray", pm, os);
	    const int num_per_line = 10;
	    for (int c = 0; c < num_cells; ++c) {
		if (c % num_per_line == 0) {
		    Tag::indent(os);
		}
		os << cell_foffsets[c] << ' ';
		if (c % num_per_line == num_per_line - 1
		    || c == num_cells - 1) {
		    os << '\n';
		}
	    }
	}
	{
	    pm["type"] = "UInt8";
	    pm["Name"] = "types";
	    Tag t("DataArray", pm, os);
	    const int num_per_line = 10;
	    for (int c = 0; c < num_cells; ++c) {
		if (c % num_per_line == 0) {
		    Tag::indent(os);
		}
		os << "42 ";
		if (c % num_per_line == num_per_line - 1
		    || c == num_cells - 1) {
		    os << '\n';
		}
	    }
	}
    }
    {
	pm.clear();
	pm["Scalars"] = "saturation";
	Tag celldatatag("CellData", pm, os);
	pm.clear();
	pm["type"] = "Int32";
	pm["NumberOfComponents"] = "1";
	pm["format"] = "ascii";
	pm["type"] = "Float64";
	{
	    pm["Name"] = "pressure";
	    Tag ptag("DataArray", pm, os);
	    const int num_per_line = 5;
	    for (int c = 0; c < num_cells; ++c) {
		if (c % num_per_line == 0) {
		    Tag::indent(os);
		}
		os << state.pressure()[c] << ' ';
		if (c % num_per_line == num_per_line - 1
		    || c == num_cells - 1) {
		    os << '\n';
		}
	    }	    
	}
	{
	    pm["Name"] = "saturation";
	    Tag ptag("DataArray", pm, os);
	    const int num_per_line = 5;
	    for (int c = 0; c < num_cells; ++c) {
		if (c % num_per_line == 0) {
		    Tag::indent(os);
		}
		os << state.saturation()[2*c] << ' ';
		if (c % num_per_line == num_per_line - 1
		    || c == num_cells - 1) {
		    os << '\n';
		}
	    }	    
	}
    }
}




static void toWaterSat(const std::vector<double>& sboth, std::vector<double>& sw)
{
    int num = sboth.size()/2;
    sw.resize(num);
    for (int i = 0; i < num; ++i) {
	sw[i] = sboth[2*i];
    }
}

static void toBothSat(const std::vector<double>& sw, std::vector<double>& sboth)
{
    int num = sw.size();
    sboth.resize(2*num);
    for (int i = 0; i < num; ++i) {
	sboth[2*i] = sw[i];
	sboth[2*i + 1] = 1.0 - sw[i];
    }
}




// --------------- Types needed to define transport solver ---------------

class SimpleFluid2pWrappingProps
{
public:
    SimpleFluid2pWrappingProps(const Opm::IncompPropertiesInterface& props)
	: props_(props)
    {
	if (props.numPhases() != 2) {
	    THROW("SimpleFluid2pWrapper requires 2 phases.");
	}
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
	double pc[2];
	double dpc[4];
	props_.capPress(1, &s[0], &c, pc, dpc);
	pcap = pc[0];
	ASSERT(pc[1] == 0.0);
	dpcap = dpc[0];
	ASSERT(dpc[1] == 0.0);
	ASSERT(dpc[2] == 0.0);
	ASSERT(dpc[3] == 0.0);
    }
 
    /// \todo Properly implement s_min() and s_max().
    ///       We must think about how to do this in
    ///       the *Properties* classes.
    double s_min(int c) const { (void) c; return 0.0; }
    double s_max(int c) const { (void) c; return 1.0; }

private:
    const Opm::IncompPropertiesInterface& props_;
};

typedef SimpleFluid2pWrappingProps TwophaseFluid;
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



// ----------------- Main program -----------------
int
main(int argc, char** argv)
{
    std::cout << "\n================    Test program for incompressible two-phase flow     ===============\n\n";
    Opm::parameter::ParameterGroup param(argc, argv, false);
    std::cout << "---------------    Reading parameters     ---------------" << std::endl;

    // Reading various control parameters.
    const int num_psteps = param.getDefault("num_psteps", 1);
    const double stepsize_days = param.getDefault("stepsize_days", 1.0);
    const double stepsize = Opm::unit::convert::from(stepsize_days, Opm::unit::day);
    const bool guess_old_solution = param.getDefault("guess_old_solution", false);
    const bool use_reorder = param.getDefault("use_reorder", true);
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
    boost::scoped_ptr<Opm::Grid> grid;
    boost::scoped_ptr<Opm::IncompPropertiesInterface> props;
    if (use_deck) {
	std::string deck_filename = param.get<std::string>("deck_filename");
	Opm::EclipseGridParser deck(deck_filename);
	// Grid init
	grid.reset(new Opm::Grid(deck));
	// Rock and fluid init
	const int* gc = grid->c_grid()->global_cell;
	std::vector<int> global_cell(gc, gc + grid->c_grid()->number_of_cells);
	props.reset(new Opm::IncompPropertiesFromDeck(deck, global_cell));
    } else {
	// Grid init.
	const int nx = param.getDefault("nx", 100);
	const int ny = param.getDefault("ny", 100);
	const int nz = param.getDefault("nz", 1);
	grid.reset(new Opm::Grid(nx, ny, nz));
	// Rock and fluid init.
	props.reset(new Opm::IncompPropertiesBasic(param, grid->c_grid()->dimensions, grid->c_grid()->number_of_cells));
    }

    // Extra rock init.
    std::vector<double> porevol;
    compute_porevolume(grid->c_grid(), *props, porevol);
    double tot_porevol = std::accumulate(porevol.begin(), porevol.end(), 0.0);

    // Extra fluid init for transport solver.
    TwophaseFluid fluid(*props);

    // Solvers init.
    PressureSolver psolver(grid->c_grid(), *props);
    TransportModel  model  (fluid, *grid->c_grid(), porevol, 0, guess_old_solution);
    TransportSolver tsolver(model);

    // State-related and source-related variables init.
    std::vector<double> totmob;
    ReservoirState state(grid->c_grid(), props->numPhases());
    // We need a separate reorder_sat, because the reorder
    // code expects a scalar sw, not both sw and so.
    std::vector<double> reorder_sat(grid->c_grid()->number_of_cells);
    double flow_per_sec = 0.1*tot_porevol/Opm::unit::day;
    std::vector<double> src   (grid->c_grid()->number_of_cells, 0.0);
    src[0]                         =  flow_per_sec;
    src[grid->c_grid()->number_of_cells - 1] = -flow_per_sec;
    TransportSource* tsrc = create_transport_source(2, 2);
    double ssrc[]   = { 1.0, 0.0 };
    double ssink[]  = { 0.0, 1.0 };
    double zdummy[] = { 0.0, 0.0 };
    append_transport_source(0, 2, 0, src[0], ssrc, zdummy, tsrc);
    append_transport_source(grid->c_grid()->number_of_cells - 1, 2, 0,
			    src.back(), ssink, zdummy, tsrc);
    std::vector<double> reorder_src = src;

    // Control init.
    Opm::ImplicitTransportDetails::NRReport  rpt;
    Opm::ImplicitTransportDetails::NRControl ctrl;
    double current_time = 0.0;
    double total_time = stepsize*num_psteps;
    ctrl.max_it = param.getDefault("max_it", 20);
    ctrl.verbosity = param.getDefault("verbosity", 0);
    ctrl.max_it_ls = param.getDefault("max_it_ls", 5);

    // Linear solver init.
    using Opm::ImplicitTransportLinAlgSupport::CSRMatrixUmfpackSolver;
    CSRMatrixUmfpackSolver linsolve;

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
    std::cout << "\n\n================    Starting main simulation loop     ===============" << std::endl;
    for (int pstep = 0; pstep < num_psteps; ++pstep) {
        std::cout << "\n\n---------------    Simulation step number " << pstep
                  << "    ---------------"
                  << "\n      Current time (days)     " << Opm::unit::convert::to(current_time, Opm::unit::day)
                  << "\n      Current stepsize (days) " << Opm::unit::convert::to(stepsize, Opm::unit::day)
                  << "\n      Total time (days)       " << Opm::unit::convert::to(total_time, Opm::unit::day)
                  << "\n" << std::endl;

	if (output) {
	    outputState(grid->c_grid(), state, pstep, output_dir);
	}

	compute_totmob(*props, state.saturation(), totmob);
	psolver.solve(grid->c_grid(), totmob, src, state);

	if (use_reorder) {
	    toWaterSat(state.saturation(), reorder_sat);
	    // We must treat reorder_src here,
	    // if we are to handle anything but simple water
	    // injection, since it is expected to be
	    // equal to total outflow (if negative)
	    // and water inflow (if positive).
	    // Also, for anything but noflow boundaries,
	    // boundary flows must be accumulated into
	    // source term following the same convention.
	    twophasetransport(&porevol[0],
			      &reorder_src[0],
			      stepsize,
			      const_cast<UnstructuredGrid*>(grid->c_grid()),
			      props.get(),
			      &state.faceflux()[0],
			      &reorder_sat[0]);
	    toBothSat(reorder_sat, state.saturation());
	} else {
	    tsolver.solve(*grid->c_grid(), tsrc, stepsize, ctrl, state, linsolve, rpt);
	    std::cout << rpt;
	}

	current_time += stepsize;
    }

    if (output) {
	outputState(grid->c_grid(), state, num_psteps, output_dir);
    }

    destroy_transport_source(tsrc);
}
