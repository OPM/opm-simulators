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


#include <opm/core/WellsManager.hpp>
#include <opm/core/eclipse/EclipseGridParser.hpp>
#include <opm/core/grid.h>
#include <opm/core/newwells.h>
#include <opm/core/utility/ErrorMacros.hpp>
#include <opm/core/utility/Units.hpp>

#include <tr1/array>
#include <cmath>


// Helper structs and functions for the implementation.
namespace
{

    struct WellData
    {
	well_type type;
	control_type control;
	double target;
	double reference_bhp_depth;
	surface_component injected_phase;
    };


    struct PerfData
    {
	int cell;
	double well_index;
    };


    int prod_control_mode(const std::string& control)
    {
	const int num_prod_control_modes = 8;
	static std::string prod_control_modes[num_prod_control_modes] =
	    {std::string("ORAT"), std::string("WRAT"), std::string("GRAT"),
	     std::string("LRAT"), std::string("RESV"), std::string("BHP"),
	     std::string("THP"), std::string("GRUP") };
	int m = -1;
	for (int i=0; i<num_prod_control_modes; ++i) {
	    if (control == prod_control_modes[i]) {
		m = i;
		break;
	    }
	}
	if (m >= 0) {
	    return m;
	} else {
	    THROW("Unknown well control mode = " << control << " in input file");
	}
    }


    int inje_control_mode(const std::string& control)
    {
	const int num_inje_control_modes = 5;
	static std::string inje_control_modes[num_inje_control_modes] =
	    {std::string("RATE"), std::string("RESV"), std::string("BHP"),
	     std::string("THP"), std::string("GRUP") };
	int m = -1;
	for (int i=0; i<num_inje_control_modes; ++i) {
	    if (control == inje_control_modes[i]) {
		m = i;
		break;
	    }
	}
 
	if (m >= 0) {
	    return m;
	} else {
	    THROW("Unknown well control mode = " << control << " in input file");
	}
    }


    std::tr1::array<double, 3> getCubeDim(const UnstructuredGrid& grid, int cell)
    {
	using namespace std;
	tr1::array<double, 3> cube;
	int num_local_faces = grid.cell_facepos[cell + 1] - grid.cell_facepos[cell];
	vector<double> x(num_local_faces);
	vector<double> y(num_local_faces);
	vector<double> z(num_local_faces);
	for (int lf=0; lf<num_local_faces; ++ lf) {
	    int face = grid.cell_faces[grid.cell_facepos[cell] + lf];
	    const double* centroid = &grid.face_centroids[grid.dimensions*face];
	    x[lf] = centroid[0];
	    y[lf] = centroid[1];
	    z[lf] = centroid[2];
	}
	cube[0] = *max_element(x.begin(), x.end()) - *min_element(x.begin(), x.end());
	cube[1] = *max_element(y.begin(), y.end()) - *min_element(y.begin(), y.end());
	cube[2] = *max_element(z.begin(), z.end()) - *min_element(z.begin(), z.end());
	return cube;
    }

    // Use the Peaceman well model to compute well indices.
    // radius is the radius of the well.
    // cubical contains [dx, dy, dz] of the cell.
    // (Note that the well model asumes that each cell is a cuboid).
    // cell_permeability is the permeability tensor of the given cell.
    // returns the well index of the cell.
    double computeWellIndex(const double radius,
			    const std::tr1::array<double, 3>& cubical,
			    const double* cell_permeability,
			    const double skin_factor)
    {
	using namespace std;
	// sse: Using the Peaceman model.
	// NOTE: The formula is valid for cartesian grids, so the result can be a bit
	// (in worst case: there is no upper bound for the error) off the mark.
	const double permx = cell_permeability[0];
	const double permy = cell_permeability[3*1 + 1];
	double effective_perm = sqrt(permx*permy);
	// sse: The formula for r_0 can be found on page 39 of
	// "Well Models for Mimetic Finite Differerence Methods and Improved Representation
	//  of Wells in Multiscale Methods" by Ingeborg Skjelkvåle Ligaarden.
	assert(permx > 0.0);
	assert(permy > 0.0);
	double kxoy = permx / permy;
	double kyox = permy / permx;
	double r0_denominator = pow(kyox, 0.25) + pow(kxoy, 0.25);
	double r0_numerator = sqrt((sqrt(kyox)*cubical[0]*cubical[0]) +
				   (sqrt(kxoy)*cubical[1]*cubical[1]));
	assert(r0_denominator > 0.0);
	double r0 = 0.28 * r0_numerator / r0_denominator;
	assert(radius > 0.0);
	assert(r0 > 0.0);
	if (r0 < radius) {
	    std::cout << "ERROR: Too big well radius detected.";
	    std::cout << "Specified well radius is " << radius
		      << " while r0 is " << r0 << ".\n";
	}
        const long double two_pi = 6.2831853071795864769252867665590057683943387987502116419498;
	double wi_denominator = log(r0 / radius) + skin_factor;
	double wi_numerator = two_pi * cubical[2];
	assert(wi_denominator > 0.0);
	double wi = effective_perm * wi_numerator / wi_denominator;
	assert(wi > 0.0);
	return wi;
    }

} // anonymous namespace





namespace Opm
{


    /// Default constructor.
    WellsManager::WellsManager()
	: w_(0)
    {
    }



    /// Construct wells from deck.
    WellsManager::WellsManager(const Opm::EclipseGridParser& deck,
			       const UnstructuredGrid& grid,
			       const double* permeability)
	: w_(0)
    {
	if (grid.dimensions != 3) {
	    THROW("We cannot initialize wells from a deck unless the corresponding grid is 3-dimensional.");
	}
	// NOTE: Implementation copied and modified from dune-porsol's class BlackoilWells.
	std::vector<std::string> keywords;
	keywords.push_back("WELSPECS");
	keywords.push_back("COMPDAT");
// 	keywords.push_back("WELTARG");
	if (!deck.hasFields(keywords)) {
	    MESSAGE("Missing well keywords in deck, initializing no wells.");
            return;
	}
	if (!(deck.hasField("WCONINJE") || deck.hasField("WCONPROD")) ) {
	    THROW("Needed field is missing in file");
	}

	// These data structures will be filled in this constructor,
	// then used to initialize the Wells struct.
	std::vector<std::string> well_names;
        std::vector<WellData> well_data;
	std::vector<std::vector<PerfData> > wellperf_data;

	// Get WELSPECS data
	const WELSPECS& welspecs = deck.getWELSPECS();
	const int num_wells = welspecs.welspecs.size();
	well_names.reserve(num_wells);
	well_data.reserve(num_wells);
	wellperf_data.resize(num_wells);
	for (int w = 0; w < num_wells; ++w) {
	    well_names.push_back(welspecs.welspecs[w].name_);
	    WellData wd;
	    well_data.push_back(wd);
            well_data.back().reference_bhp_depth = welspecs.welspecs[w].datum_depth_BHP_;
            if (welspecs.welspecs[w].datum_depth_BHP_ < 0.0) {
                // Set refdepth to a marker value, will be changed
                // after getting perforation data to the centroid of
                // the cell of the top well perforation.
                well_data.back().reference_bhp_depth = -1e100;
	    }
        }

	// global_cell is a map from compressed cells to Cartesian grid cells.
	// We must make the inverse lookup.
	const int* global_cell = grid.global_cell;
	const int* cpgdim = grid.cartdims;
	std::map<int,int> cartesian_to_compressed;
	for (int i = 0; i < grid.number_of_cells; ++i) {
	    cartesian_to_compressed.insert(std::make_pair(global_cell[i], i));
	}

	// Get COMPDAT data   
	const COMPDAT& compdat = deck.getCOMPDAT();
	const int num_compdat  = compdat.compdat.size();
	for (int kw = 0; kw < num_compdat; ++kw) {
	    // Extract well name, or the part of the well name that
	    // comes before the '*'.
	    std::string name = compdat.compdat[kw].well_;
	    std::string::size_type len = name.find('*');
	    if (len != std::string::npos) {
		name = name.substr(0, len);
	    }
	    // Look for well with matching name.
	    bool found = false;
	    for (int wix = 0; wix < num_wells; ++wix) {
		if (well_names[wix].compare(0,len, name) == 0) { // equal
		    // We have a matching name.
		    int ix = compdat.compdat[kw].grid_ind_[0] - 1;
		    int jy = compdat.compdat[kw].grid_ind_[1] - 1;
		    int kz1 = compdat.compdat[kw].grid_ind_[2] - 1;
                    int kz2 = compdat.compdat[kw].grid_ind_[3] - 1;
                    for (int kz = kz1; kz <= kz2; ++kz) {
                        int cart_grid_indx = ix + cpgdim[0]*(jy + cpgdim[1]*kz);
                        std::map<int, int>::const_iterator cgit = 
                            cartesian_to_compressed.find(cart_grid_indx);
                        if (cgit == cartesian_to_compressed.end()) {
                            THROW("Cell with i,j,k indices " << ix << ' ' << jy << ' '
                                  << kz << " not found in grid!");
                        }
                        int cell = cgit->second;
                        PerfData pd;
                        pd.cell = cell;
                        if (compdat.compdat[kw].connect_trans_fac_ > 0.0) {
                            pd.well_index = compdat.compdat[kw].connect_trans_fac_;
                        } else {
                            double radius = 0.5*compdat.compdat[kw].diameter_;
                            if (radius <= 0.0) {
                                radius = 0.5*unit::feet;
                                MESSAGE("**** Warning: Well bore internal radius set to " << radius);
                            }
			    std::tr1::array<double, 3> cubical = getCubeDim(grid, cell);
                            const double* cell_perm = &permeability[grid.dimensions*grid.dimensions*cell];  
                            pd.well_index = computeWellIndex(radius, cubical, cell_perm,
                                                             compdat.compdat[kw].skin_factor_);
                        }
                        wellperf_data[wix].push_back(pd);
                    }
                    found = true;
                    break;
                }
            }
	    if (!found) {
		THROW("Undefined well name: " << compdat.compdat[kw].well_
		      << " in COMPDAT");
	    }
	}

	// Set up reference depths that were defaulted. Count perfs.
	int num_perfs = 0;
        for (int w = 0; w < num_wells; ++w) {
	    num_perfs += wellperf_data[w].size();
            if (well_data[w].reference_bhp_depth == -1e100) {
                // It was defaulted. Set reference depth to minimum perforation depth.
                double min_depth = 1e100;
                int num_wperfs = wellperf_data[w].size();
                for (int perf = 0; perf < num_wperfs; ++perf) {
                    double depth = grid.cell_centroids[3*wellperf_data[w][perf].cell + 2];
                    min_depth = std::min(min_depth, depth);
                }
                well_data[w].reference_bhp_depth = min_depth;
            }
        }
 
	// Get WCONINJE data
        if (deck.hasField("WCONINJE")) {
            const WCONINJE& wconinjes = deck.getWCONINJE();
            const int num_wconinjes = wconinjes.wconinje.size();
            for (int kw = 0; kw < num_wconinjes; ++kw) {
		// Extract well name, or the part of the well name that
		// comes before the '*'.
                std::string name = wconinjes.wconinje[kw].well_;
                std::string::size_type len = name.find('*');
                if (len != std::string::npos) {
                    name = name.substr(0, len);
                }
                bool well_found = false;
                for (int wix = 0; wix < num_wells; ++wix) {
                    if (well_names[wix].compare(0,len, name) == 0) { //equal
                        well_found = true;
                        well_data[wix].type = INJECTOR;
                        int m = inje_control_mode(wconinjes.wconinje[kw].control_mode_);
                        switch(m) {
                        case 0:  // RATE
                            well_data[wix].control = RATE;
                            well_data[wix].target = wconinjes.wconinje[kw].surface_flow_max_rate_;
                            break;
                        case 1:  // RESV
                            well_data[wix].control = RATE;
                            well_data[wix].target = wconinjes.wconinje[kw].fluid_volume_max_rate_;
                            break;
                        case 2:  // BHP
                            well_data[wix].control = BHP;
                            well_data[wix].target = wconinjes.wconinje[kw].BHP_limit_;
                            break;
                        case 3:  // THP
                            well_data[wix].control = BHP;
                            well_data[wix].target = wconinjes.wconinje[kw].THP_limit_;
                            break;
                        case 4:
                            break;
                        default:
                            THROW("Unknown well control mode; WCONIJE  = "
                                  << wconinjes.wconinje[kw].control_mode_
                                  << " in input file");
                        }
                        if (wconinjes.wconinje[kw].injector_type_ == "WATER") {
                            well_data[wix].injected_phase = WATER;
                        } else if (wconinjes.wconinje[kw].injector_type_ == "OIL") {
                            well_data[wix].injected_phase = OIL;
                        } else if (wconinjes.wconinje[kw].injector_type_ == "GAS") {
                            well_data[wix].injected_phase = GAS;
                        } else {
			    THROW("Error in injector specification, found no known fluid type.");
			}
		    }
		}
                if (!well_found) {
                    THROW("Undefined well name: " << wconinjes.wconinje[kw].well_
                          << " in WCONINJE");
                }
            }
        }

	// Get WCONPROD data
        if (deck.hasField("WCONPROD")) {
            const WCONPROD& wconprods = deck.getWCONPROD();
            const int num_wconprods   = wconprods.wconprod.size();
            for (int kw = 0; kw < num_wconprods; ++kw) {
                std::string name = wconprods.wconprod[kw].well_;
                std::string::size_type len = name.find('*');
                if (len != std::string::npos) {
                    name = name.substr(0, len);
                }

                bool well_found = false;
                for (int wix = 0; wix < num_wells; ++wix) {
                    if (well_names[wix].compare(0,len, name) == 0) { //equal
                        well_found = true;
                        well_data[wix].type = PRODUCER;
                        int m = prod_control_mode(wconprods.wconprod[kw].control_mode_);
                        switch(m) {
                        case 0:  // ORAT
                            well_data[wix].control = RATE;
                            well_data[wix].target = wconprods.wconprod[kw].oil_max_rate_;
                            break;
                        case 1:  // WRAT
                            well_data[wix].control = RATE;
                            well_data[wix].target = wconprods.wconprod[kw].water_max_rate_;
                            break;
                        case 2:  // GRAT
                            well_data[wix].control = RATE;
                            well_data[wix].target = wconprods.wconprod[kw].gas_max_rate_;
                            break;
                        case 3:  // LRAT
                            well_data[wix].control = RATE;
                            well_data[wix].target = wconprods.wconprod[kw].liquid_max_rate_;
                            break;
                        case 4:  // RESV 
                            well_data[wix].control = RATE;
                            well_data[wix].target = wconprods.wconprod[kw].fluid_volume_max_rate_;
                            break;
                        case 5:  // BHP
                            well_data[wix].control = BHP; 
                            well_data[wix].target = wconprods.wconprod[kw].BHP_limit_;
                            break;
                        case 6:  // THP 
                            well_data[wix].control = BHP;
                            well_data[wix].target = wconprods.wconprod[kw].THP_limit_;
                            break;
                        case 7:
                            // Handle group here.
                            break;
                        default:
                            THROW("Unknown well control mode; WCONPROD  = "
                                  << wconprods.wconprod[kw].control_mode_
                                  << " in input file");
                        }
                    }
                }
                if (!well_found) {
                    THROW("Undefined well name: " << wconprods.wconprod[kw].well_
                          << " in WCONPROD");
                }
            }
        }

	// Get WELTARG data
        if (deck.hasField("WELTARG")) {
            const WELTARG& weltargs = deck.getWELTARG();
            const int num_weltargs  = weltargs.weltarg.size();	
            for (int kw = 0; kw < num_weltargs; ++kw) {
                std::string name = weltargs.weltarg[kw].well_;
                std::string::size_type len = name.find('*');
                if (len != std::string::npos) {
                    name = name.substr(0, len);
                }
                bool well_found = false;
                for (int wix = 0; wix < num_wells; ++wix) {
                    if (well_names[wix].compare(0,len, name) == 0) { //equal
                        well_found = true;
                        well_data[wix].target = weltargs.weltarg[kw].new_value_;
                        break;
                    }
                }
                if (!well_found) {
                    THROW("Undefined well name: " << weltargs.weltarg[kw].well_
                          << " in WELTARG");
                }
            }
        }

        // Debug output.
#define EXTRA_OUTPUT
#ifdef EXTRA_OUTPUT
	std::cout << "\t WELL DATA" << std::endl;
	for(int i = 0; i< num_wells; ++i) {
	    std::cout << i << ": " << well_data[i].type << "  "
		      << well_data[i].control << "  " << well_data[i].target
		      << std::endl;
	}

	std::cout << "\n\t PERF DATA" << std::endl;
	for(int i=0; i< int(wellperf_data.size()); ++i) {
	    for(int j=0; j< int(wellperf_data[i].size()); ++j) {
		std::cout << i << ": " << wellperf_data[i][j].cell << "  "
			  << wellperf_data[i][j].well_index << std::endl;
	    }
	}
#endif

	// Set up the Wells struct.
	w_ = wells_create(num_wells, num_perfs);
	if (!w_) {
	    THROW("Failed creating Wells struct.");
	}
	double fracs[3][3] = { { 1.0, 0.0, 0.0 },
			       { 0.0, 1.0, 0.0 },
			       { 0.0, 0.0, 1.0 } };
	for (int w = 0; w < num_wells; ++w) {
	    int nperf = wellperf_data[w].size();
	    std::vector<int> cells(nperf);
	    std::vector<double> wi(nperf);
	    for (int perf = 0; perf < nperf; ++perf) {
		cells[perf] = wellperf_data[w][perf].cell;
		wi[perf] = wellperf_data[w][perf].well_index;
	    }
	    const double* zfrac = (well_data[w].type == INJECTOR) ? fracs[well_data[w].injected_phase] : 0;
	    int ok = wells_add(well_data[w].type, well_data[w].reference_bhp_depth, nperf,
			       zfrac, &cells[0], &wi[0], w_);
	    if (!ok) {
		THROW("Failed to add a well.");
	    }
	    // We only append a single control at this point.
	    // TODO: Handle multiple controls.
	    ok = well_controls_append(well_data[w].control, well_data[w].target, w_->ctrls[w]);
	    if (!ok) {
		THROW("Failed to add well controls.");
	    }
	}
        

    }



    /// Destructor.
    WellsManager::~WellsManager()
    {
	wells_destroy(w_);
    }




    /// Access the managed Wells.
    /// The method is named similarly to c_str() in std::string,
    /// to make it clear that we are returning a C-compatible struct.
    const Wells* WellsManager::c_wells() const
    {
	return w_;
    }




} // namespace Opm
