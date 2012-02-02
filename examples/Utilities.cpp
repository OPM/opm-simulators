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

#include "config.h"

#include "Utilities.hpp"

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

#include <opm/core/transport/CSRMatrixUmfpackSolver.hpp>

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
#include <numeric>




namespace Opm
{






    void
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


    void
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



    void writeVtkDataAllCartesian(const std::tr1::array<int, 3>& dims,
				  const std::tr1::array<double, 3>& cell_size,
				  const std::vector<double>& pressure,
				  const std::vector<double>& saturation,
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
	    vtk_file << pressure[i] << '\n';
	}

	vtk_file << "SCALARS saturation float" << '\n';
	vtk_file << "LOOKUP_TABLE saturation_table " << '\n';
	for (int i = 0; i < num_cells; ++i) {
	    double s = saturation[2*i];
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


    void writeVtkDataGeneralGrid(const UnstructuredGrid* grid,
				 const std::vector<double>& pressure,
				 const std::vector<double>& saturation,
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
		    os << pressure[c] << ' ';
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
		    os << saturation[2*c] << ' ';
		    if (c % num_per_line == num_per_line - 1
			|| c == num_cells - 1) {
			os << '\n';
		    }
		}	    
	    }
	}
    }




    void toWaterSat(const std::vector<double>& sboth, std::vector<double>& sw)
    {
	int num = sboth.size()/2;
	sw.resize(num);
	for (int i = 0; i < num; ++i) {
	    sw[i] = sboth[2*i];
	}
    }

    void toBothSat(const std::vector<double>& sw, std::vector<double>& sboth)
    {
	int num = sw.size();
	sboth.resize(2*num);
	for (int i = 0; i < num; ++i) {
	    sboth[2*i] = sw[i];
	    sboth[2*i + 1] = 1.0 - sw[i];
	}
    }



} // namespace Opm

