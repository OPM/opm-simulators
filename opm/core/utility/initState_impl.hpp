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

#ifndef OPM_INITSTATE_IMPL_HEADER_INCLUDED
#define OPM_INITSTATE_IMPL_HEADER_INCLUDED


#include <opm/core/utility/parameters/ParameterGroup.hpp>
#include <opm/core/grid.h>
#include <opm/core/eclipse/EclipseGridParser.hpp>
#include <opm/core/utility/ErrorMacros.hpp>
#include <opm/core/utility/MonotCubicInterpolator.hpp>
#include <opm/core/utility/Units.hpp>
#include <opm/core/fluid/IncompPropertiesInterface.hpp>
#include <opm/core/fluid/BlackoilPropertiesInterface.hpp>
#include <cmath>

namespace Opm
{

    namespace
    {

        // Find the cells that are below and above a depth.
        // TODO: add 'anitialiasing', obtaining a more precise split
        //       by f. ex. subdividing cells cut by the split depths.
        void cellsBelowAbove(const UnstructuredGrid& grid,
                             const double depth,
                             std::vector<int>& below,
                             std::vector<int>& above)
        {
            const int num_cells = grid.number_of_cells;
            below.reserve(num_cells);
            above.reserve(num_cells);
            const int dim = grid.dimensions;
            for (int c = 0; c < num_cells; ++c) {
                const double z = grid.cell_centroids[dim*c + dim - 1];
                if (z > depth) {
                    below.push_back(c);
                } else {
                    above.push_back(c);
                }
            }
        }

        enum WaterInit { WaterBelow, WaterAbove };

        // Initialize saturations so that there is water below woc,
        // and oil above.
        // If invert is true, water is instead above, oil below.
        template <class Props, class State>
        void initWaterOilContact(const UnstructuredGrid& grid,
                                 const Props& props,
                                 const double woc,
                                 const WaterInit waterinit,
                                 State& state)
        {
            // Find out which cells should have water and which should have oil.
            std::vector<int> water;
            std::vector<int> oil;
            // Assuming that water should go below the woc, but warning if oil is heavier.
            // if (props.density()[0] < props.density()[1]) {
            //     std::cout << "*** warning: water density is less than oil density, "
            //         "but initialising water below woc." << std::endl;
            // }
            switch (waterinit) {
            case WaterBelow:
                cellsBelowAbove(grid, woc, water, oil);
                break;
            case WaterAbove:
                cellsBelowAbove(grid, woc, oil, water);
            }
            // Set saturations.
            state.setFirstSat(oil, props, State::MinSat);
            state.setFirstSat(water, props, State::MaxSat);
        }


        // Initialize hydrostatic pressures depending only on gravity,
        // (constant) phase densities and a water-oil contact depth.
        // The pressure difference between points is equal to
        //     delta_p = delta_z * gravity * rho
        // where rho is the (assumed constant) oil density above the
        // woc, water density below woc.
        // Note that even if there is (immobile) water in the oil zone,
        // it does not contribute to the pressure difference.
        // Note that by manipulating the densities parameter,
        // it is possible to initialise with water on top instead of
        // on the bottom etc.
        template <class State>
        void initHydrostaticPressure(const UnstructuredGrid& grid,
                                     const double* densities,
                                     const double woc,
                                     const double gravity,
                                     const double datum_z,
                                     const double datum_p,
                                     State& state)
        {
            std::vector<double>& p = state.pressure();
            const int num_cells = grid.number_of_cells;
            const int dim = grid.dimensions;
            // Compute pressure at woc
            const double rho_datum = datum_z > woc ? densities[0] : densities[1];
            const double woc_p = datum_p + (woc - datum_z)*gravity*rho_datum;
            for (int c = 0; c < num_cells; ++c) {
                // Compute pressure as delta from woc pressure.
                const double z = grid.cell_centroids[dim*c + dim - 1];
                const double rho = z > woc ? densities[0] : densities[1];
                p[c] = woc_p + (z - woc)*gravity*rho;
            }
        }


        // Facade to initHydrostaticPressure() taking a property object,
        // for similarity to initHydrostaticPressure() for compressible fluids.
        template <class State>
        void initHydrostaticPressure(const UnstructuredGrid& grid,
                                     const IncompPropertiesInterface& props,
                                     const double woc,
                                     const double gravity,
                                     const double datum_z,
                                     const double datum_p,
                                     State& state)
        {
            const double* densities = props.density();
            initHydrostaticPressure(grid, densities, woc, gravity, datum_z, datum_p, state);
        }



        // Helper functor for initHydrostaticPressure() for compressible fluids.
        struct Density
        {
            const BlackoilPropertiesInterface& props_;
            Density(const BlackoilPropertiesInterface& props) : props_(props) {}
            double operator()(const double pressure, const int phase)
            {
                ASSERT(props_.numPhases() == 2);
                const double surfvol[2][2] = { { 1.0, 0.0 },
                                               { 0.0, 1.0 } };
                // We do not handle multi-region PVT/EQUIL at this point.
                const int* cells = 0;
                double A[4] = { 0.0 };
                props_.matrix(1, &pressure, surfvol[phase], cells, A, 0);
                double rho[2] = { 0.0 };
                props_.density(1, A, rho);
                return rho[phase];
            }
        };

        // Initialize hydrostatic pressures depending only on gravity,
        // phase densities that may vary with pressure and a water-oil
        // contact depth. The pressure ODE is given as
        //     \grad p = \rho g \grad z
        // where rho is the oil density above the woc, water density
        // below woc. Note that even if there is (immobile) water in
        // the oil zone, it does not contribute to the pressure there.
        template <class State>
        void initHydrostaticPressure(const UnstructuredGrid& grid,
                                     const BlackoilPropertiesInterface& props,
                                     const double woc,
                                     const double gravity,
                                     const double datum_z,
                                     const double datum_p,
                                     State& state)
        {
            ASSERT(props.numPhases() == 2);

            // Obtain max and min z for which we will need to compute p.
            const int num_cells = grid.number_of_cells;
            const int dim = grid.dimensions;
            double zlim[2] = { 1e100, -1e100 };
            for (int c = 0; c < num_cells; ++c) {
                const double z = grid.cell_centroids[dim*c + dim - 1];
                zlim[0] = std::min(zlim[0], z);
                zlim[1] = std::max(zlim[1], z);
            }

            // We'll use a minimum stepsize of 1/100 of the z range.
            const double hmin = (zlim[1] - zlim[0])/100.0;

            // Store p(z) results in an associative array.
            std::map<double, double> press_by_z;
            press_by_z[datum_z] = datum_p;

            // Set up density evaluator.
            Density rho(props);

            // Solve the ODE from datum_z to woc.
            int phase = (datum_z > woc) ? 0 : 1;
            int num_steps = int(std::ceil(std::fabs(woc - datum_z)/hmin));
            double pval = datum_p;
            double zval = datum_z;
            double h = (woc - datum_z)/double(num_steps);
            for (int i = 0; i < num_steps; ++i) {
                zval += h;
                const double dp = rho(pval, phase)*gravity;
                pval += h*dp;
                press_by_z[zval] = pval;
            }
            double woc_p = pval;

            // Solve the ODE from datum_z to the end of the interval.
            double z_end = (datum_z > woc) ? zlim[1] : zlim[0];
            num_steps = int(std::ceil(std::fabs(z_end - datum_z)/hmin));
            pval = datum_p;
            zval = datum_z;
            h = (z_end - datum_z)/double(num_steps);
            for (int i = 0; i < num_steps; ++i) {
                zval += h;
                const double dp = rho(pval, phase)*gravity;
                pval += h*dp;
                press_by_z[zval] = pval;
            }

            // Solve the ODE from woc to the other end of the interval.
            // Switching phase and z_end.
            phase = (phase + 1) % 2;
            z_end = (datum_z > woc) ? zlim[0] : zlim[1];
            pval = woc_p;
            zval = woc;
            h = (z_end - datum_z)/double(num_steps);
            for (int i = 0; i < num_steps; ++i) {
                zval += h;
                const double dp = rho(pval, phase)*gravity;
                pval += h*dp;
                press_by_z[zval] = pval;
            }

            // Create monotone spline for interpolating solution.
            std::vector<double> zv, pv;
            zv.reserve(press_by_z.size());
            pv.reserve(press_by_z.size());
            for (std::map<double, double>::const_iterator it = press_by_z.begin();
                 it != press_by_z.end(); ++it) {
                zv.push_back(it->first);
                pv.push_back(it->second);
            }
            MonotCubicInterpolator press(zv, pv);

            // Evaluate pressure at each cell centroid.
            std::vector<double>& p = state.pressure();
            for (int c = 0; c < num_cells; ++c) {
                const double z = grid.cell_centroids[dim*c + dim - 1];
                p[c] = press(z);
            }
        }

    } // anonymous namespace


    /// Initialize a twophase state from parameters.
    /// The following parameters are accepted (defaults):
    ///    convection_testcase   (false)    Water in the 'left' part of the grid.
    ///    ref_pressure          (100)      Initial pressure in bar for all cells
    ///                                     (if convection_testcase is true),
    ///                                     or pressure at woc depth.
    ///    segregation_testcase  (false)    Water above the woc instead of below.
    ///    water_oil_contact     (none)     Depth of water-oil contact (woc).
    ///    init_saturation       (none)     Initial water saturation for all cells.
    /// If convection_testcase is true, the saturation is initialised
    /// as indicated, and pressure is initialised to a constant value
    /// ('ref_pressure').
    /// If segregation_testcase is true, the saturation is initialised
    /// as indicated, and pressure is initialised hydrostatically.
    /// Otherwise we have 3 cases:
    ///   1) If 'water_oil_contact' is given, saturation is initialised
    ///      accordingly.
    ///   2) If 'water_oil_contact' is not given, but 'init_saturation'
    ///      is given, water saturation is set to that value everywhere.
    ///   3) If neither are given, water saturation is set to minimum.
    /// In all three cases, pressure is initialised hydrostatically.
    /// In case 2) and 3), the depth of the first cell is used as reference depth.
    template <class State>
    void initStateBasic(const UnstructuredGrid& grid,
                        const IncompPropertiesInterface& props,
                        const parameter::ParameterGroup& param,
                        const double gravity,
                        State& state)
    {
        const int num_phases = props.numPhases();
        if (num_phases != 2) {
            THROW("initStateTwophaseBasic(): currently handling only two-phase scenarios.");
        }
        state.init(grid, num_phases);
        const int num_cells = props.numCells();
        // By default: initialise water saturation to minimum everywhere.
        std::vector<int> all_cells(num_cells);
        for (int i = 0; i < num_cells; ++i) {
            all_cells[i] = i;
        }
        state.setFirstSat(all_cells, props, State::MinSat);
        const bool convection_testcase = param.getDefault("convection_testcase", false);
        const bool segregation_testcase = param.getDefault("segregation_testcase", false);
        if (convection_testcase) {
            // Initialise water saturation to max in the 'left' part.
            std::vector<int> left_cells;
            left_cells.reserve(num_cells/2);
            const int *glob_cell = grid.global_cell;
            const int* cd = grid.cartdims;
            for (int cell = 0; cell < num_cells; ++cell) {
                const int gc = glob_cell == 0 ? cell : glob_cell[cell];
                bool left = (gc % cd[0]) < cd[0]/2;
                if (left) {
                    left_cells.push_back(cell);
                }
            }
            state.setFirstSat(left_cells, props, State::MaxSat);
            const double init_p = param.getDefault("ref_pressure", 100)*unit::barsa;
            std::fill(state.pressure().begin(), state.pressure().end(), init_p);
        } else if (segregation_testcase) {
            // Warn against error-prone usage.
            if (gravity == 0.0) {
                std::cout << "**** Warning: running gravity segregation scenario, but gravity is zero." << std::endl;
            }
            if (grid.cartdims[2] <= 1) {
                std::cout << "**** Warning: running gravity segregation scenario, which expects nz > 1." << std::endl;
            }
            // Initialise water saturation to max *above* water-oil contact.
            const double woc = param.get<double>("water_oil_contact");
            initWaterOilContact(grid, props, woc, WaterAbove, state);
            // Initialise pressure to hydrostatic state.
            const double ref_p = param.getDefault("ref_pressure", 100)*unit::barsa;
            double dens[2] = { props.density()[1], props.density()[0] };
            initHydrostaticPressure(grid, dens, woc, gravity, woc, ref_p, state);
        } else if (param.has("water_oil_contact")) {
            // Warn against error-prone usage.
            if (gravity == 0.0) {
                std::cout << "**** Warning: running gravity convection scenario, but gravity is zero." << std::endl;
            }
            if (grid.cartdims[2] <= 1) {
                std::cout << "**** Warning: running gravity convection scenario, which expects nz > 1." << std::endl;
            }
            // Initialise water saturation to max below water-oil contact.
            const double woc = param.get<double>("water_oil_contact");
            initWaterOilContact(grid, props, woc, WaterBelow, state);
            // Initialise pressure to hydrostatic state.
            const double ref_p = param.getDefault("ref_pressure", 100)*unit::barsa;
            initHydrostaticPressure(grid, props.density(), woc, gravity, woc, ref_p, state);
        } else if (param.has("init_saturation")) {
            // Initialise water saturation to init_saturation parameter.
            const double init_saturation = param.get<double>("init_saturation");
            for (int cell = 0; cell < num_cells; ++cell) {
                state.saturation()[2*cell] = init_saturation;
                state.saturation()[2*cell + 1] = 1.0 - init_saturation;
            }
            // Initialise pressure to hydrostatic state.
            const double ref_p = param.getDefault("ref_pressure", 100)*unit::barsa;
            const double rho =  props.density()[0]*init_saturation + props.density()[1]*(1.0 - init_saturation);
            const double dens[2] = { rho, rho };
            const double ref_z = grid.cell_centroids[0 + grid.dimensions - 1];
            initHydrostaticPressure(grid, dens, ref_z, gravity, ref_z, ref_p, state);
        } else {
            // Use default: water saturation is minimum everywhere.
            // Initialise pressure to hydrostatic state.
            const double ref_p = param.getDefault("ref_pressure", 100)*unit::barsa;
            const double rho =  props.density()[1];
            const double dens[2] = { rho, rho };
            const double ref_z = grid.cell_centroids[0 + grid.dimensions - 1];
            initHydrostaticPressure(grid, dens, ref_z, gravity, ref_z, ref_p, state);
        }
    }


    /// Initialize a blackoil state from parameters.
    /// The following parameters are accepted (defaults):
    ///    convection_testcase   (false)    Water in the 'left' part of the grid.
    ///    ref_pressure          (100)      Initial pressure in bar for all cells
    ///                                     (if convection_testcase is true),
    ///                                     or pressure at woc depth.
    ///    water_oil_contact     (none)     Depth of water-oil contact (woc).
    /// If convection_testcase is true, the saturation is initialised
    /// as indicated, and pressure is initialised to a constant value
    /// ('ref_pressure').
    /// Otherwise we have 2 cases:
    ///   1) If 'water_oil_contact' is given, saturation is initialised
    ///      accordingly.
    ///   2) Water saturation is set to minimum.
    /// In both cases, pressure is initialised hydrostatically.
    /// In case 2), the depth of the first cell is used as reference depth.
    template <class State>
    void initStateBasic(const UnstructuredGrid& grid,
                        const BlackoilPropertiesInterface& props,
                        const parameter::ParameterGroup& param,
                        const double gravity,
                        State& state)
    {
        // TODO: Refactor to exploit similarity with IncompProp* case.
        const int num_phases = props.numPhases();
        if (num_phases != 2) {
            THROW("initStateTwophaseBasic(): currently handling only two-phase scenarios.");
        }
        state.init(grid, num_phases);
        const int num_cells = props.numCells();
        // By default: initialise water saturation to minimum everywhere.
        std::vector<int> all_cells(num_cells);
        for (int i = 0; i < num_cells; ++i) {
            all_cells[i] = i;
        }
        state.setFirstSat(all_cells, props, State::MinSat);
        const bool convection_testcase = param.getDefault("convection_testcase", false);
        if (convection_testcase) {
            // Initialise water saturation to max in the 'left' part.
            std::vector<int> left_cells;
            left_cells.reserve(num_cells/2);
            const int *glob_cell = grid.global_cell;
            const int* cd = grid.cartdims;
            for (int cell = 0; cell < num_cells; ++cell) {
                const int gc = glob_cell == 0 ? cell : glob_cell[cell];
                bool left = (gc % cd[0]) < cd[0]/2;
                if (left) {
                    left_cells.push_back(cell);
                }
            }
            state.setFirstSat(left_cells, props, State::MaxSat);
            const double init_p = param.getDefault("ref_pressure", 100)*unit::barsa;
            std::fill(state.pressure().begin(), state.pressure().end(), init_p);
        } else if (param.has("water_oil_contact")) {
            // Warn against error-prone usage.
            if (gravity == 0.0) {
                std::cout << "**** Warning: running gravity convection scenario, but gravity is zero." << std::endl;
            }
            if (grid.cartdims[2] <= 1) {
                std::cout << "**** Warning: running gravity convection scenario, which expects nz > 1." << std::endl;
            }
            // Initialise water saturation to max below water-oil contact.
            const double woc = param.get<double>("water_oil_contact");
            initWaterOilContact(grid, props, woc, WaterBelow, state);
            // Initialise pressure to hydrostatic state.
            const double ref_p = param.getDefault("ref_pressure", 100)*unit::barsa;
            initHydrostaticPressure(grid, props, woc, gravity, woc, ref_p, state);
        } else {
            // Use default: water saturation is minimum everywhere.
            // Initialise pressure to hydrostatic state.
            const double ref_p = param.getDefault("ref_pressure", 100)*unit::barsa;
            const double ref_z = grid.cell_centroids[0 + grid.dimensions - 1];
            const double woc = -1e100;
            initHydrostaticPressure(grid, props, woc, gravity, ref_z, ref_p, state);
        }
    }


    /// Initialize a state from input deck.
    /// If EQUIL is present:
    ///   - saturation is set according to the water-oil contact,
    ///   - pressure is set to hydrostatic equilibrium.
    /// Otherwise:
    ///   - saturation is set according to SWAT,
    ///   - pressure is set according to PRESSURE.
    template <class Props, class State>
    void initStateFromDeck(const UnstructuredGrid& grid,
                           const Props& props,
                           const EclipseGridParser& deck,
                           const double gravity,
                           State& state)
    {
        const int num_phases = props.numPhases();
        if (num_phases != 2) {
            THROW("initStateFromDeck(): currently handling only two-phase scenarios.");
        }
        state.init(grid, num_phases);
        if (deck.hasField("EQUIL")) {
            // Set saturations depending on oil-water contact.
            const EQUIL& equil= deck.getEQUIL();
            if (equil.equil.size() != 1) {
                THROW("No region support yet.");
            }
            const double woc = equil.equil[0].water_oil_contact_depth_;
            initWaterOilContact(grid, props, woc, WaterBelow, state);
            // Set pressure depending on densities and depths.
            const double datum_z = equil.equil[0].datum_depth_;
            const double datum_p = equil.equil[0].datum_depth_pressure_;
            initHydrostaticPressure(grid, props, woc, gravity, datum_z, datum_p, state);
        } else if (deck.hasField("SWAT") && deck.hasField("PRESSURE")) {
            // Set saturations from SWAT, pressure from PRESSURE.
            std::vector<double>& s = state.saturation();
            std::vector<double>& p = state.pressure();
            const std::vector<double>& sw_deck = deck.getFloatingPointValue("SWAT");
            const std::vector<double>& p_deck = deck.getFloatingPointValue("PRESSURE");
            const int num_cells = grid.number_of_cells;
            for (int c = 0; c < num_cells; ++c) {
                int c_deck = (grid.global_cell == NULL) ? c : grid.global_cell[c];
                s[2*c] = sw_deck[c_deck];
                s[2*c + 1] = 1.0 - s[2*c];
                p[c] = p_deck[c_deck];
            }
        } else {
            THROW("initStateFromDeck(): we must either have EQUIL, or both SWAT and PRESSURE.");
        }
    }


} // namespace Opm


#endif // OPM_INITSTATE_IMPL_HEADER_INCLUDED
