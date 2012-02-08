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

#ifndef OPM_UTILITIES_HEADER_INCLUDED
#define OPM_UTILITIES_HEADER_INCLUDED

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




    class PressureSolver
    {
    public:
	PressureSolver(const UnstructuredGrid* g,
		       const IncompPropertiesInterface& props)
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
	    std::fill(gpress_.begin(), gpress_.end(), double(0.0));

	    ifs_tpfa_assemble(gg, &trans_[0], &src[0], &gpress_[0], h_);

	    using ImplicitTransportLinAlgSupport::CSRMatrixUmfpackSolver;

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




    void
    compute_porevolume(const UnstructuredGrid* g,
		       const Opm::IncompPropertiesInterface& props,
		       std::vector<double>& porevol);


    void
    compute_totmob(const Opm::IncompPropertiesInterface& props,
		   const std::vector<double>& s,
		   std::vector<double>& totmob);




    void writeVtkDataAllCartesian(const std::tr1::array<int, 3>& dims,
				  const std::tr1::array<double, 3>& cell_size,
				  const std::vector<double>& pressure,
				  const std::vector<double>& saturation,
				  std::ostream& vtk_file);



    typedef std::map<std::string, const std::vector<double>*> DataMap;

    void writeVtkDataGeneralGrid(const UnstructuredGrid* grid,
				 const DataMap& data,
				 std::ostream& os);




    void toWaterSat(const std::vector<double>& sboth, std::vector<double>& sw);

    void toBothSat(const std::vector<double>& sw, std::vector<double>& sboth);



} // namespace Opm




#endif // OPM_UTILITIES_HEADER_INCLUDED
