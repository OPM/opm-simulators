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
#include <opm/core/WellCollection.hpp>
#include <opm/core/fluid/blackoil/phaseUsageFromDeck.hpp>

#include <tr1/array>
#include <cmath>


// Helper structs and functions for the implementation.
namespace
{

    struct WellData
    {
        WellType type;
        // WellControlType control;
        // double target;
        double reference_bhp_depth;
        // Opm::InjectionSpecification::InjectorType injected_phase;
    };


    struct PerfData
    {
        int cell;
        double well_index;
    };

    namespace ProductionControl
    {
        enum Mode { ORAT, WRAT, GRAT,
                    LRAT, CRAT, RESV,
                    BHP, THP, GRUP };
        Mode mode(const std::string& control)
        {
            const int num_prod_control_modes = 9;
            static std::string prod_control_modes[num_prod_control_modes] =
                {std::string("ORAT"), std::string("WRAT"), std::string("GRAT"),
                 std::string("LRAT"), std::string("CRAT"), std::string("RESV"),
                 std::string("BHP"), std::string("THP"), std::string("GRUP") };
            int m = -1;
            for (int i=0; i<num_prod_control_modes; ++i) {
                if (control == prod_control_modes[i]) {
                    m = i;
                    break;
                }
            }
            if (m >= 0) {
                return static_cast<Mode>(m);
            } else {
                THROW("Unknown well control mode = " << control << " in input file");
            }
        }
    } // namespace ProductionControl


    namespace InjectionControl
    {
        enum Mode { RATE, RESV, BHP,
                    THP, GRUP };
        Mode mode(const std::string& control)
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
                return static_cast<Mode>(m);
            } else {
                THROW("Unknown well control mode = " << control << " in input file");
            }
        }
    } // namespace InjectionControl


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
//      keywords.push_back("WELTARG");
        if (!deck.hasFields(keywords)) {
            MESSAGE("Missing well keywords in deck, initializing no wells.");
            return;
        }
        if (!(deck.hasField("WCONINJE") || deck.hasField("WCONPROD")) ) {
            THROW("Needed field is missing in file");
        }

        // Obtain phase usage data.
        PhaseUsage pu = phaseUsageFromDeck(deck);

        // These data structures will be filled in this constructor,
        // then used to initialize the Wells struct.
        std::vector<std::string> well_names;
        std::vector<WellData> well_data;
        std::vector<std::vector<PerfData> > wellperf_data;

        // For easy lookup:
        std::map<std::string, int> well_names_to_index;

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
            well_names_to_index[welspecs.welspecs[w].name_] = w;
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

        if (global_cell) {
            for (int i = 0; i < grid.number_of_cells; ++i) {
                cartesian_to_compressed.insert(std::make_pair(global_cell[i], i));
            }
        }
        else {
            for (int i = 0; i < grid.number_of_cells; ++i) {
                cartesian_to_compressed.insert(std::make_pair(i, i));
            }
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
                    // Extract corresponding WELSPECS defintion for
                    // purpose of default location specification.
                    const WelspecsLine& wspec = welspecs.welspecs[wix];

                    // We have a matching name.
                    int ix = compdat.compdat[kw].grid_ind_[0] - 1;
                    int jy = compdat.compdat[kw].grid_ind_[1] - 1;
                    int kz1 = compdat.compdat[kw].grid_ind_[2] - 1;
                    int kz2 = compdat.compdat[kw].grid_ind_[3] - 1;

                    if (ix < 0) {
                        // Defaulted I location.  Extract from WELSPECS.
                        ix = wspec.I_ - 1;
                    }
                    if (jy < 0) {
                        // Defaulted J location.  Extract from WELSPECS.
                        jy = wspec.J_ - 1;
                    }
                    if (kz1 < 0) {
                        // Defaulted KZ1.  Use top layer.
                        kz1 = 0;
                    }
                    if (kz2 < 0) {
                        // Defaulted KZ2.  Use bottom layer.
                        kz2 = cpgdim[2] - 1;
                    }

                    for (int kz = kz1; kz <= kz2; ++kz) {
                        int cart_grid_indx = ix + cpgdim[0]*(jy + cpgdim[1]*kz);
                        std::map<int, int>::const_iterator cgit =
                            cartesian_to_compressed.find(cart_grid_indx);
                        if (cgit == cartesian_to_compressed.end()) {
                            THROW("Cell with i,j,k indices " << ix << ' ' << jy << ' '
                                  << kz << " not found in grid (well = " << name << ')');
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
        ASSERT(grid.dimensions == 3);
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

        // Create the well data structures.
        w_ = create_wells(pu.num_phases, num_wells, num_perfs);
        if (!w_) {
            THROW("Failed creating Wells struct.");
        }

        // Classify wells
        if (deck.hasField("WCONINJE")) {
            const std::vector<WconinjeLine>& lines = deck.getWCONINJE().wconinje;
            for (size_t i = 0 ; i < lines.size(); ++i) {
                const std::map<std::string, int>::const_iterator it = well_names_to_index.find(lines[i].well_);
                if (it != well_names_to_index.end()) {
                    const int well_index = it->second;
                    well_data[well_index].type = INJECTOR;
                }
                else {
                    THROW("Unseen well name: " << lines[i].well_ << " first seen in WCONINJE");
                }
            }
        }
        if (deck.hasField("WCONPROD")) {
            const std::vector<WconprodLine>& lines = deck.getWCONPROD().wconprod;
            for (size_t i = 0; i < lines.size(); ++i) {
                const std::map<std::string, int>::const_iterator it = well_names_to_index.find(lines[i].well_);
                if (it != well_names_to_index.end()) {
                    const int well_index = it->second;
                    well_data[well_index].type = PRODUCER;
                } else {
                    THROW("Unseen well name: " << lines[i].well_ << " first seen in WCONPROD");
                }
                
            }
        }
        
        // Add wells.
        for (int w = 0; w < num_wells; ++w) {
            const int w_num_perf = wellperf_data[w].size();
            std::vector<int> perf_cells(w_num_perf);
            std::vector<double> perf_prodind(w_num_perf);
            for (int perf = 0; perf < w_num_perf; ++perf) {
                perf_cells[perf] = wellperf_data[w][perf].cell;
                perf_prodind[perf] = wellperf_data[w][perf].well_index;
            }
            const double* comp_frac = NULL;
            // We initialize all wells with a null component fraction,
            // and must (for injection wells) overwrite it later.
            int ok = add_well(well_data[w].type, well_data[w].reference_bhp_depth, w_num_perf,
                              comp_frac, &perf_cells[0], &perf_prodind[0], w_);
            if (!ok) {
                THROW("Failed adding well " << well_names[w] << " to Wells data structure.");
            }
        }

        // Get WCONINJE data, add injection controls to wells.
        if (deck.hasField("WCONINJE")) {
            const WCONINJE& wconinjes = deck.getWCONINJE();
            const int num_wconinjes = wconinjes.wconinje.size();
            for (int kw = 0; kw < num_wconinjes; ++kw) {
                const WconinjeLine& wci_line = wconinjes.wconinje[kw];
                // Extract well name, or the part of the well name that
                // comes before the '*'.
                std::string name = wci_line.well_;
                std::string::size_type len = name.find('*');
                if (len != std::string::npos) {
                    name = name.substr(0, len);
                }
                bool well_found = false;
                for (int wix = 0; wix < num_wells; ++wix) {
                    if (well_names[wix].compare(0,len, name) == 0) { //equal
                        well_found = true;
                        ASSERT(well_data[wix].type == w_->type[wix]);
                        if (well_data[wix].type != INJECTOR) {
                            THROW("Found WCONINJE entry for a non-injector well: " << well_names[wix]);
                        }

                        // Add all controls that are present in well.
                        int ok = 1;
                        int control_pos[5] = { -1, -1, -1, -1, -1 };
                        if (ok && wci_line.surface_flow_max_rate_ >= 0.0) {
                            control_pos[InjectionControl::RATE] = w_->ctrls[wix]->num;
                            const double distr[3] = { 1.0, 1.0, 1.0 };
                            ok = append_well_controls(SURFACE_RATE, wci_line.surface_flow_max_rate_,
                                                      distr, wix, w_);
                        }
                        if (ok && wci_line.reservoir_flow_max_rate_ >= 0.0) {
                            control_pos[InjectionControl::RESV] = w_->ctrls[wix]->num;
                            const double distr[3] = { 1.0, 1.0, 1.0 };
                            ok = append_well_controls(RESERVOIR_RATE, wci_line.reservoir_flow_max_rate_,
                                                      distr, wix, w_);
                        }
                        if (ok && wci_line.BHP_limit_ > 0.0) {
                            control_pos[InjectionControl::BHP] = w_->ctrls[wix]->num;
                            ok = append_well_controls(BHP, wci_line.BHP_limit_,
                                                      NULL, wix, w_);
                        }
                        if (ok && wci_line.THP_limit_ > 0.0) {
                            THROW("We cannot handle THP limit for well " << well_names[wix]);
                        }
                        if (!ok) {
                            THROW("Failure occured appending controls for well " << well_names[wix]);
                        }
                        InjectionControl::Mode mode = InjectionControl::mode(wci_line.control_mode_);
                        const int cpos = control_pos[mode];
                        if (cpos == -1 && mode != InjectionControl::GRUP) {
                            THROW("Control for " << wci_line.control_mode_ << " not specified in well " << well_names[wix]);
                        }
                        set_current_control(wix, cpos, w_);

                        // Set well component fraction.
                        double cf[3] = { 0.0, 0.0, 0.0 };
                        if (wci_line.injector_type_[0] == 'W') {
                            if (!pu.phase_used[BlackoilPhases::Aqua]) {
                                THROW("Water phase not used, yet found water-injecting well.");
                            }
                            cf[pu.phase_pos[BlackoilPhases::Aqua]] = 1.0;
                        } else if (wci_line.injector_type_[0] == 'O') {
                            if (!pu.phase_used[BlackoilPhases::Liquid]) {
                                THROW("Oil phase not used, yet found oil-injecting well.");
                            }
                            cf[pu.phase_pos[BlackoilPhases::Liquid]] = 1.0;
                        } else if (wci_line.injector_type_[0] == 'G') {
                            if (!pu.phase_used[BlackoilPhases::Vapour]) {
                                THROW("Water phase not used, yet found water-injecting well.");
                            }
                            cf[pu.phase_pos[BlackoilPhases::Vapour]] = 1.0;
                        }
                        std::copy(cf, cf + pu.num_phases, w_->comp_frac + wix*pu.num_phases);
                    }
                }
                if (!well_found) {
                    THROW("Undefined well name: " << wci_line.well_
                          << " in WCONINJE");
                }
            }
        }

        // Get WCONPROD data
        if (deck.hasField("WCONPROD")) {
            const WCONPROD& wconprods = deck.getWCONPROD();
            const int num_wconprods   = wconprods.wconprod.size();
            for (int kw = 0; kw < num_wconprods; ++kw) {
                const WconprodLine& wcp_line = wconprods.wconprod[kw];
                std::string name = wcp_line.well_;
                std::string::size_type len = name.find('*');
                if (len != std::string::npos) {
                    name = name.substr(0, len);
                }
                bool well_found = false;
                for (int wix = 0; wix < num_wells; ++wix) {
                    if (well_names[wix].compare(0,len, name) == 0) { //equal
                        well_found = true;
                        ASSERT(well_data[wix].type == w_->type[wix]);
                        if (well_data[wix].type != PRODUCER) {
                            THROW("Found WCONPROD entry for a non-producer well: " << well_names[wix]);
                        }
                        // Add all controls that are present in well.
                        int control_pos[9] = { -1, -1, -1, -1, -1, -1, -1, -1, -1 };
                        int ok = 1;
                        if (ok && wcp_line.oil_max_rate_ >= 0.0) {
                            if (!pu.phase_used[BlackoilPhases::Liquid]) {
                                THROW("Oil phase not active and ORAT control specified.");
                            }
                            control_pos[ProductionControl::ORAT] = w_->ctrls[wix]->num;
                            double distr[3] = { 0.0, 0.0, 0.0 };
                            distr[pu.phase_pos[BlackoilPhases::Liquid]] = 1.0;
                            ok = append_well_controls(SURFACE_RATE, -wcp_line.oil_max_rate_,
                                                      distr, wix, w_);
                        }
                        if (ok && wcp_line.water_max_rate_ >= 0.0) {
                            if (!pu.phase_used[BlackoilPhases::Aqua]) {
                                THROW("Water phase not active and WRAT control specified.");
                            }
                            control_pos[ProductionControl::WRAT] = w_->ctrls[wix]->num;
                            double distr[3] = { 0.0, 0.0, 0.0 };
                            distr[pu.phase_pos[BlackoilPhases::Aqua]] = 1.0;
                            ok = append_well_controls(SURFACE_RATE, -wcp_line.water_max_rate_,
                                                      distr, wix, w_);
                        }
                        if (ok && wcp_line.gas_max_rate_ >= 0.0) {
                            if (!pu.phase_used[BlackoilPhases::Vapour]) {
                                THROW("Gas phase not active and GRAT control specified.");
                            }
                            control_pos[ProductionControl::GRAT] = w_->ctrls[wix]->num;
                            double distr[3] = { 0.0, 0.0, 0.0 };
                            distr[pu.phase_pos[BlackoilPhases::Vapour]] = 1.0;
                            ok = append_well_controls(SURFACE_RATE, -wcp_line.gas_max_rate_,
                                                      distr, wix, w_);
                        }
                        if (ok && wcp_line.liquid_max_rate_ >= 0.0) {
                            if (!pu.phase_used[BlackoilPhases::Aqua]) {
                                THROW("Water phase not active and LRAT control specified.");
                            }
                            if (!pu.phase_used[BlackoilPhases::Liquid]) {
                                THROW("Oil phase not active and LRAT control specified.");
                            }
                            control_pos[ProductionControl::LRAT] = w_->ctrls[wix]->num;
                            double distr[3] = { 0.0, 0.0, 0.0 };
                            distr[pu.phase_pos[BlackoilPhases::Aqua]] = 1.0;
                            distr[pu.phase_pos[BlackoilPhases::Liquid]] = 1.0;
                            ok = append_well_controls(SURFACE_RATE, -wcp_line.liquid_max_rate_,
                                                      distr, wix, w_);
                        }
                        if (ok && wcp_line.reservoir_flow_max_rate_ >= 0.0) {
                            control_pos[ProductionControl::RESV] = w_->ctrls[wix]->num;
                            double distr[3] = { 1.0, 1.0, 1.0 };
                            ok = append_well_controls(RESERVOIR_RATE, -wcp_line.reservoir_flow_max_rate_,
                                                 distr, wix, w_);
                        }
                        if (ok && wcp_line.BHP_limit_ > 0.0) {
                            control_pos[ProductionControl::BHP] = w_->ctrls[wix]->num;
                            ok = append_well_controls(BHP, wcp_line.BHP_limit_,
                                                 NULL, wix, w_);
                        }
                        if (ok && wcp_line.THP_limit_ > 0.0) {
                            THROW("We cannot handle THP limit for well " << well_names[wix]);
                        }
                        if (!ok) {
                            THROW("Failure occured appending controls for well " << well_names[wix]);
                        }
                        ProductionControl::Mode mode = ProductionControl::mode(wcp_line.control_mode_);
                        const int cpos = control_pos[mode];
                        if (cpos == -1 && mode != ProductionControl::GRUP) {
                            THROW("Control mode type " << mode << " not present in well " << well_names[wix]);
                        }
                        set_current_control(wix, cpos, w_);
                    }
                }
                if (!well_found) {
                    THROW("Undefined well name: " << wcp_line.well_
                          << " in WCONPROD");
                }
            }
        }

        // Get WELTARG data
        if (deck.hasField("WELTARG")) {
            THROW("We currently do not handle WELTARG.");
            /*
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
            */
        }

        // Debug output.
#define EXTRA_OUTPUT
#ifdef EXTRA_OUTPUT
        /*
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
        */
#endif


        // Build the well_collection_ well group hierarchy.
        if (deck.hasField("GRUPTREE")) {
            std::cout << "Found gruptree" << std::endl;
            const GRUPTREE& gruptree = deck.getGRUPTREE();
            std::map<std::string, std::string>::const_iterator it = gruptree.tree.begin();
            for( ; it != gruptree.tree.end(); ++it) {
                well_collection_.addChild(it->first, it->second, deck);
            }
        }
        for (size_t i = 0; i < welspecs.welspecs.size(); ++i) {
            WelspecsLine line = welspecs.welspecs[i];
            well_collection_.addChild(line.name_, line.group_, deck);
        }



        // Set the guide rates:
        if (deck.hasField("WGRUPCON")) {
            std::cout << "Found Wgrupcon" << std::endl;
            WGRUPCON wgrupcon = deck.getWGRUPCON();
            const std::vector<WgrupconLine>& lines = wgrupcon.wgrupcon;
            std::cout << well_collection_.getLeafNodes().size() << std::endl;
            for (size_t i = 0; i < lines.size(); i++) {
                std::string name = lines[i].well_;
                const int wix = well_names_to_index[name];
                WellNode& wellnode = *well_collection_.getLeafNodes()[wix];
                ASSERT(wellnode.name() == name);
                if (well_data[wix].type == PRODUCER) {
                    wellnode.prodSpec().guide_rate_ = lines[i].guide_rate_;
                    if (lines[i].phase_ == "OIL") {
                        wellnode.prodSpec().guide_rate_type_ = ProductionSpecification::OIL;
                    } else {
                        THROW("Guide rate type " << lines[i].phase_ << " specified for producer "
                              << name << " in WGRUPCON, cannot handle.");
                    }
                } else if (well_data[wix].type == INJECTOR) {
                    wellnode.injSpec().guide_rate_ = lines[i].guide_rate_;
                    if (lines[i].phase_ == "RAT") {
                        wellnode.injSpec().guide_rate_type_ = InjectionSpecification::RAT;
                    } else {
                        THROW("Guide rate type " << lines[i].phase_ << " specified for injector "
                              << name << " in WGRUPCON, cannot handle.");
                    }
                } else {
                    THROW("Unknown well type " << well_data[wix].type << " for well " << name);
                }
            }
        }
        well_collection_.setWellsPointer(w_);
        well_collection_.applyGroupControls();
    }



    /// Destructor.
    WellsManager::~WellsManager()
    {
        destroy_wells(w_);
    }


    /// Does the "deck" define any wells?
    bool WellsManager::empty() const
    {
        return (w_ == 0) || (w_->number_of_wells == 0);
    }



    /// Access the managed Wells.
    /// The method is named similarly to c_str() in std::string,
    /// to make it clear that we are returning a C-compatible struct.
    const Wells* WellsManager::c_wells() const
    {
        return w_;
    }

    const WellCollection& WellsManager::wellCollection() const
    {
        return well_collection_;
    }

    bool WellsManager::conditionsMet(const std::vector<double>& well_bhp,
                                     const std::vector<double>& well_reservoirrates_phase,
                                     const std::vector<double>& well_surfacerates_phase)
    {
        return well_collection_.conditionsMet(well_bhp,
                                              well_reservoirrates_phase,
                                              well_surfacerates_phase);
    }

    /// Applies explicit reinjection controls. This must be called at each timestep to be correct.
    /// \param[in]    well_reservoirrates_phase
    ///                         A vector containing reservoir rates by phase for each well.
    ///                         Is assumed to be ordered the same way as the related Wells-struct,
    ///                         with all phase rates of a single well adjacent in the array.
    /// \param[in]    well_surfacerates_phase
    ///                         A vector containing surface rates by phase for each well.
    ///                         Is assumed to be ordered the same way as the related Wells-struct,
    ///                         with all phase rates of a single well adjacent in the array.

    void WellsManager::applyExplicitReinjectionControls(const std::vector<double>& well_reservoirrates_phase,
                                                        const std::vector<double>& well_surfacerates_phase)
    {
        well_collection_.applyExplicitReinjectionControls(well_reservoirrates_phase, well_surfacerates_phase);
    }

} // namespace Opm
