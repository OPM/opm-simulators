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
#include <opm/core/grid/GridHelpers.hpp>
#include <opm/core/io/eclipse/EclipseGridParser.hpp>
#include <opm/core/utility/ErrorMacros.hpp>
#include <opm/core/utility/MonotCubicInterpolator.hpp>
#include <opm/core/utility/Units.hpp>
#include <opm/core/props/IncompPropertiesInterface.hpp>
#include <opm/core/props/BlackoilPropertiesInterface.hpp>
#include <opm/core/props/phaseUsageFromDeck.hpp>
#include <opm/core/utility/miscUtilitiesBlackoil.hpp>

#include <opm/parser/eclipse/Utility/EquilWrapper.hpp>

#include <iostream>
#include <cmath>

namespace Opm
{

    namespace
    {
#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunneeded-internal-declaration"
#endif /* __clang__ */
        // Find the cells that are below and above a depth.
        // TODO: add 'anitialiasing', obtaining a more precise split
        //       by f. ex. subdividing cells cut by the split depths.
       template<class T>
        void cellsBelowAbove(int number_of_cells,
                             T begin_cell_centroids,
                             int dimension,
                             const double depth,
                             std::vector<int>& below,
                             std::vector<int>& above)
        {
            below.reserve(number_of_cells);
            above.reserve(number_of_cells);
            for (int c = 0; c < number_of_cells; ++c) {
                const double z = UgGridHelpers
                    ::getCoordinate(UgGridHelpers::increment(begin_cell_centroids, c,
                                                             dimension),
                                    dimension-1);
                if (z > depth) {
                    below.push_back(c);
                } else {
                    above.push_back(c);
                }
            }
        }
#ifdef __clang__
#pragma clang diagnostic pop
#endif /* __clang__ */

        enum WaterInit { WaterBelow, WaterAbove };

        // Initialize saturations so that there is water below woc,
        // and oil above.
        // If invert is true, water is instead above, oil below.
        template <class T, class Props, class State>
        void initWaterOilContact(int number_of_cells,
                                 T begin_cell_centroids,
                                 int dimensions,
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
                cellsBelowAbove(number_of_cells, begin_cell_centroids, dimensions,
                                woc, water, oil);
                break;
            case WaterAbove:
                cellsBelowAbove(number_of_cells, begin_cell_centroids, dimensions,
                                woc, oil, water);
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
        template <class T, class State>
        void initHydrostaticPressure(int number_of_cells,
                                     T begin_cell_centroids,
                                     int dimensions,
                                     const double* densities,
                                     const double woc,
                                     const double gravity,
                                     const double datum_z,
                                     const double datum_p,
                                     State& state)
        {
            std::vector<double>& p = state.pressure();
            // Compute pressure at woc
            const double rho_datum = datum_z > woc ? densities[0] : densities[1];
            const double woc_p = datum_p + (woc - datum_z)*gravity*rho_datum;
            for (int c = 0; c < number_of_cells; ++c) {
                // Compute pressure as delta from woc pressure.
                const double z = UgGridHelpers::
                    getCoordinate(UgGridHelpers::increment(begin_cell_centroids, c, dimensions),
                                  dimensions-1);
                const double rho = z > woc ? densities[0] : densities[1];
                p[c] = woc_p + (z - woc)*gravity*rho;
            }
        }


        // Facade to initHydrostaticPressure() taking a property object,
        // for similarity to initHydrostaticPressure() for compressible fluids.
        template <class T, class State>
        void initHydrostaticPressure(int number_of_cells,
                                     T begin_cell_centroids,
                                     int dimensions,
                                     const IncompPropertiesInterface& props,
                                     const double woc,
                                     const double gravity,
                                     const double datum_z,
                                     const double datum_p,
                                     State& state)
        {
            const double* densities = props.density();
            initHydrostaticPressure(number_of_cells, begin_cell_centroids, dimensions,
                                    densities, woc, gravity, datum_z, datum_p, state);
        }



        // Helper functor for initHydrostaticPressure() for compressible fluids.
        struct Density
        {
            const BlackoilPropertiesInterface& props_;
            Density(const BlackoilPropertiesInterface& props) : props_(props) {}
            double operator()(const double pressure, const int phase)
            {
                assert(props_.numPhases() == 2);
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
       template <class T, class State>
        void initHydrostaticPressure(int number_of_cells,
                                     T begin_cell_centroids,
                                     int dimensions,
                                     const BlackoilPropertiesInterface& props,
                                     const double woc,
                                     const double gravity,
                                     const double datum_z,
                                     const double datum_p,
                                     State& state)
        {
            assert(props.numPhases() == 2);

            // Obtain max and min z for which we will need to compute p.
            double zlim[2] = { 1e100, -1e100 };
            for (int c = 0; c < number_of_cells; ++c) {
                const double z = UgGridHelpers::
                    getCoordinate(UgGridHelpers::increment(begin_cell_centroids, c, dimensions),
                                  dimensions-1);;
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
            for (int c = 0; c < number_of_cells; ++c) {
                const double z = UgGridHelpers::
                    getCoordinate(UgGridHelpers::increment(begin_cell_centroids, c, dimensions),
                                  dimensions-1);
                p[c] = press(z);
            }
        }

        // Initialize face pressures to distance-weighted average of adjacent cell pressures.
    template <class State, class FaceCells, class FCI, class CCI>
        void initFacePressure(int dimensions,
                              int number_of_faces,
                              FaceCells face_cells,
                              FCI begin_face_centroids,
                              CCI begin_cell_centroids,
                              State& state)
        {
            const std::vector<double>& cp = state.pressure();
            std::vector<double>& fp = state.facepressure();
            for (int f = 0; f < number_of_faces; ++f) {
                double dist[2] = { 0.0, 0.0 };
                double press[2] = { 0.0, 0.0 };
                int bdy_idx = -1;
                for (int j = 0; j < 2; ++j) {
                    const int c = face_cells(f, j);
                    if (c >= 0) {
                        dist[j] = 0.0;
                        for (int dd = 0; dd < dimensions; ++dd) {
                            double diff = UgGridHelpers::
                                getCoordinate(UgGridHelpers::increment(begin_face_centroids, f,
                                                                       dimensions),
                                              dd)
                                - UgGridHelpers::
                                getCoordinate(UgGridHelpers::increment(begin_cell_centroids, c,
                                                                       dimensions),
                                              dd);
                            dist[j] += diff*diff;
                        }
                        dist[j] = std::sqrt(dist[j]);
                        press[j] = cp[c];
                    } else {
                        bdy_idx = j;
                    }
                }
                if (bdy_idx == -1) {
                    fp[f] = press[0]*(dist[1]/(dist[0] + dist[1])) + press[1]*(dist[0]/(dist[0] + dist[1]));
                } else {
                    fp[f] = press[(bdy_idx + 1) % 2];
                }
            }
        }


    } // anonymous namespace


    /// Initialize a twophase state from parameters.
    template <class State>
    void initStateBasic(const UnstructuredGrid& grid,
                        const IncompPropertiesInterface& props,
                        const parameter::ParameterGroup& param,
                        const double gravity,
                        State& state)
    {
        initStateBasic(grid.number_of_cells, grid.global_cell, grid.cartdims,
                       grid.number_of_faces, UgGridHelpers::faceCells(grid),
                       grid.face_centroids, grid.cell_centroids,
                       grid.dimensions, props, param,
                       gravity, state);
    }

    template <class FaceCells, class CCI, class FCI, class State>
    void initStateBasic(int number_of_cells,
                        const int* global_cell,
                        const int* cartdims,
                        int number_of_faces,
                        FaceCells face_cells,
                        FCI begin_face_centroids,
                        CCI begin_cell_centroids,
                        int dimensions,
                        const IncompPropertiesInterface& props,
                        const parameter::ParameterGroup& param,
                        const double gravity,
                        State& state)
{
        const int num_phases = props.numPhases();
        if (num_phases != 2) {
            OPM_THROW(std::runtime_error, "initStateTwophaseBasic(): currently handling only two-phase scenarios.");
        }
        state.init(number_of_cells, number_of_faces, num_phases);
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
            for (int cell = 0; cell < num_cells; ++cell) {
                const int gc = global_cell == 0 ? cell : global_cell[cell];
                bool left = (gc % cartdims[0]) < cartdims[0]/2;
                if (left) {
                    left_cells.push_back(cell);
                }
            }
            state.setFirstSat(left_cells, props, State::MaxSat);
            const double init_p = param.getDefault("ref_pressure", 100.0)*unit::barsa;
            std::fill(state.pressure().begin(), state.pressure().end(), init_p);
        } else if (segregation_testcase) {
            // Warn against error-prone usage.
            if (gravity == 0.0) {
                std::cout << "**** Warning: running gravity segregation scenario, but gravity is zero." << std::endl;
            }
            if (cartdims[2] <= 1) {
                std::cout << "**** Warning: running gravity segregation scenario, which expects nz > 1." << std::endl;
            }
            // Initialise water saturation to max *above* water-oil contact.
            const double woc = param.get<double>("water_oil_contact");
            initWaterOilContact(number_of_cells, begin_cell_centroids, dimensions,
                                props, woc, WaterAbove, state);
            // Initialise pressure to hydrostatic state.
            const double ref_p = param.getDefault("ref_pressure", 100.0)*unit::barsa;
            double dens[2] = { props.density()[1], props.density()[0] };
            initHydrostaticPressure(number_of_cells, begin_cell_centroids, dimensions,
                                    dens, woc, gravity, woc, ref_p, state);
        } else if (param.has("water_oil_contact")) {
            // Warn against error-prone usage.
            if (gravity == 0.0) {
                std::cout << "**** Warning: running gravity convection scenario, but gravity is zero." << std::endl;
            }
            if (cartdims[2] <= 1) {
                std::cout << "**** Warning: running gravity convection scenario, which expects nz > 1." << std::endl;
            }
            // Initialise water saturation to max below water-oil contact.
            const double woc = param.get<double>("water_oil_contact");
            initWaterOilContact(number_of_cells, begin_cell_centroids, dimensions,
                                props, woc, WaterBelow, state);
            // Initialise pressure to hydrostatic state.
            const double ref_p = param.getDefault("ref_pressure", 100.0)*unit::barsa;
            initHydrostaticPressure(number_of_cells, begin_cell_centroids, dimensions,
                                    props.density(), woc, gravity, woc, ref_p, state);
        } else if (param.has("init_saturation")) {
            // Initialise water saturation to init_saturation parameter.
            const double init_saturation = param.get<double>("init_saturation");
            for (int cell = 0; cell < num_cells; ++cell) {
                state.saturation()[2*cell] = init_saturation;
                state.saturation()[2*cell + 1] = 1.0 - init_saturation;
            }
            // Initialise pressure to hydrostatic state.
            const double ref_p = param.getDefault("ref_pressure", 100.0)*unit::barsa;
            const double rho =  props.density()[0]*init_saturation + props.density()[1]*(1.0 - init_saturation);
            const double dens[2] = { rho, rho };
            const double ref_z = getCoordinate(begin_cell_centroids, dimensions - 1);
            initHydrostaticPressure(number_of_cells, begin_cell_centroids, dimensions,
                                    dens, ref_z, gravity, ref_z, ref_p, state);
        } else {
            // Use default: water saturation is minimum everywhere.
            // Initialise pressure to hydrostatic state.
            const double ref_p = param.getDefault("ref_pressure", 100.0)*unit::barsa;
            const double rho =  props.density()[1];
            const double dens[2] = { rho, rho };
            const double ref_z = getCoordinate(begin_cell_centroids, dimensions - 1);
            initHydrostaticPressure(number_of_cells, begin_cell_centroids, dimensions,
                                    dens, ref_z, gravity, ref_z, ref_p, state);
        }

        // Finally, init face pressures.
        initFacePressure(dimensions, number_of_faces, face_cells, begin_face_centroids,
                         begin_cell_centroids, state);
    }


    /// Initialize a blackoil state from parameters.
    template <class State>
    void initStateBasic(const UnstructuredGrid& grid,
                        const BlackoilPropertiesInterface& props,
                        const parameter::ParameterGroup& param,
                        const double gravity,
                        State& state)
    {
        initStateBasic(grid.number_of_cells, grid.global_cell, grid.cartdims,
                       grid.number_of_faces, UgGridHelpers::faceCells(grid), 
                       grid.face_centroids, grid.cell_centroids, grid.dimensions,
                       props, param, gravity, state);
    }

    template <class FaceCells, class FCI, class CCI, class State>
    void initStateBasic(int number_of_cells,
                        const int* global_cell,
                        const int* cartdims,
                        int number_of_faces,
                        FaceCells face_cells,
                        FCI begin_face_centroids,
                        CCI begin_cell_centroids,
                        int dimensions,
                        const BlackoilPropertiesInterface& props,
                        const parameter::ParameterGroup& param,
                        const double gravity,
                        State& state)
    {
        // TODO: Refactor to exploit similarity with IncompProp* case.
        const int num_phases = props.numPhases();
        if (num_phases != 2) {
            OPM_THROW(std::runtime_error, "initStateTwophaseBasic(): currently handling only two-phase scenarios.");
        }
        state.init(number_of_cells, number_of_faces, num_phases);
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
            for (int cell = 0; cell < num_cells; ++cell) {
                const int gc = global_cell == 0 ? cell : global_cell[cell];
                bool left = (gc % cartdims[0]) < cartdims[0]/2;
                if (left) {
                    left_cells.push_back(cell);
                }
            }
            state.setFirstSat(left_cells, props, State::MaxSat);
            const double init_p = param.getDefault("ref_pressure", 100.0)*unit::barsa;
            std::fill(state.pressure().begin(), state.pressure().end(), init_p);
        } else if (param.has("water_oil_contact")) {
            // Warn against error-prone usage.
            if (gravity == 0.0) {
                std::cout << "**** Warning: running gravity convection scenario, but gravity is zero." << std::endl;
            }
            if (cartdims[2] <= 1) {
                std::cout << "**** Warning: running gravity convection scenario, which expects nz > 1." << std::endl;
            }
            // Initialise water saturation to max below water-oil contact.
            const double woc = param.get<double>("water_oil_contact");
            initWaterOilContact(number_of_cells, begin_cell_centroids, dimensions,
                                props, woc, WaterBelow, state);
            // Initialise pressure to hydrostatic state.
            const double ref_p = param.getDefault("ref_pressure", 100.0)*unit::barsa;
            initHydrostaticPressure(number_of_cells, begin_cell_centroids, dimensions,
                                    props, woc, gravity, woc, ref_p, state);
        } else {
            // Use default: water saturation is minimum everywhere.
            // Initialise pressure to hydrostatic state.
            const double ref_p = param.getDefault("ref_pressure", 100.0)*unit::barsa;
            const double ref_z = UgGridHelpers::getCoordinate(begin_cell_centroids, dimensions - 1);
            const double woc = -1e100;
            initHydrostaticPressure(number_of_cells, begin_cell_centroids, dimensions,
                                    props, woc, gravity, ref_z, ref_p, state);
        }

        // Finally, init face pressures.
        initFacePressure(dimensions, number_of_faces, face_cells, begin_face_centroids,
                         begin_cell_centroids, state);
    }


    /// Initialize a state from input deck.
    template <class Props, class State>
    void initStateFromDeck(const UnstructuredGrid& grid,
                           const Props& props,
                           Opm::DeckConstPtr newParserDeck,
                           const double gravity,
                           State& state)
    {
        const int num_phases = props.numPhases();
        const PhaseUsage pu = phaseUsageFromDeck(newParserDeck);
        if (num_phases != pu.num_phases) {
            OPM_THROW(std::runtime_error, "initStateFromDeck():  user specified property object with " << num_phases << " phases, "
                  "found " << pu.num_phases << " phases in deck.");
        }
        state.init(grid, num_phases);
        if (newParserDeck->hasKeyword("EQUIL")) {
            if (num_phases != 2) {
                OPM_THROW(std::runtime_error, "initStateFromDeck(): EQUIL-based init currently handling only two-phase scenarios.");
            }
            if (pu.phase_used[BlackoilPhases::Vapour]) {
                OPM_THROW(std::runtime_error, "initStateFromDeck(): EQUIL-based init currently handling only oil-water scenario (no gas).");
            }
            // Set saturations depending on oil-water contact.
            EquilWrapper equil(newParserDeck->getKeyword("EQUIL"));
            if (equil.numRegions() != 1) {
                OPM_THROW(std::runtime_error, "initStateFromDeck(): No region support yet.");
            }
            const double woc = equil.waterOilContactDepth(0);
            initWaterOilContact(grid, props, woc, WaterBelow, state);
            // Set pressure depending on densities and depths.
            const double datum_z = equil.datumDepth(0);
            const double datum_p = equil.datumDepthPressure(0);
            initHydrostaticPressure(grid, props, woc, gravity, datum_z, datum_p, state);
        } else if (newParserDeck->hasKeyword("PRESSURE")) {
            // Set saturations from SWAT/SGAS, pressure from PRESSURE.
            std::vector<double>& s = state.saturation();
            std::vector<double>& p = state.pressure();
            const std::vector<double>& p_deck = newParserDeck->getKeyword("PRESSURE")->getSIDoubleData();
            const int num_cells = grid.number_of_cells;
            if (num_phases == 2) {
                if (!pu.phase_used[BlackoilPhases::Aqua]) {
                    // oil-gas: we require SGAS
                    if (!newParserDeck->hasKeyword("SGAS")) {
                        OPM_THROW(std::runtime_error, "initStateFromDeck(): missing SGAS keyword in 2-phase init");
                    }
                    const std::vector<double>& sg_deck = newParserDeck->getKeyword("SGAS")->getSIDoubleData();
                    const int gpos = pu.phase_pos[BlackoilPhases::Vapour];
                    const int opos = pu.phase_pos[BlackoilPhases::Liquid];
                    for (int c = 0; c < num_cells; ++c) {
                        int c_deck = (grid.global_cell == NULL) ? c : grid.global_cell[c];
                        s[2*c + gpos] = sg_deck[c_deck];
                        s[2*c + opos] = 1.0 - sg_deck[c_deck];
                        p[c] = p_deck[c_deck];
                    }
                } else {
                    // water-oil or water-gas: we require SWAT
                    if (!newParserDeck->hasKeyword("SWAT")) {
                        OPM_THROW(std::runtime_error, "initStateFromDeck(): missing SWAT keyword in 2-phase init");
                    }
                    const std::vector<double>& sw_deck = newParserDeck->getKeyword("SWAT")->getSIDoubleData();
                    const int wpos = pu.phase_pos[BlackoilPhases::Aqua];
                    const int nwpos = (wpos + 1) % 2;
                    for (int c = 0; c < num_cells; ++c) {
                        int c_deck = (grid.global_cell == NULL) ? c : grid.global_cell[c];
                        s[2*c + wpos] = sw_deck[c_deck];
                        s[2*c + nwpos] = 1.0 - sw_deck[c_deck];
                        p[c] = p_deck[c_deck];
                    }
                }
            } else if (num_phases == 3) {
                const bool has_swat_sgas = newParserDeck->hasKeyword("SWAT") && newParserDeck->hasKeyword("SGAS");
                if (!has_swat_sgas) {
                    OPM_THROW(std::runtime_error, "initStateFromDeck(): missing SGAS or SWAT keyword in 3-phase init.");
                }
                const int wpos = pu.phase_pos[BlackoilPhases::Aqua];
                const int gpos = pu.phase_pos[BlackoilPhases::Vapour];
                const int opos = pu.phase_pos[BlackoilPhases::Liquid];
                const std::vector<double>& sw_deck = newParserDeck->getKeyword("SWAT")->getSIDoubleData();
                const std::vector<double>& sg_deck = newParserDeck->getKeyword("SGAS")->getSIDoubleData();
                for (int c = 0; c < num_cells; ++c) {
                    int c_deck = (grid.global_cell == NULL) ? c : grid.global_cell[c];
                    s[3*c + wpos] = sw_deck[c_deck];
                    s[3*c + opos] = 1.0 - (sw_deck[c_deck] + sg_deck[c_deck]);
                    s[3*c + gpos] = sg_deck[c_deck];
                    p[c] = p_deck[c_deck];
                 }
            } else {
                OPM_THROW(std::runtime_error, "initStateFromDeck(): init with SWAT etc. only available with 2 or 3 phases.");
            }
        } else {
            OPM_THROW(std::runtime_error, "initStateFromDeck(): we must either have EQUIL, or PRESSURE and SWAT/SOIL/SGAS.");
        }

        // Finally, init face pressures.
        initFacePressure(grid, state);
    }

    /// Initialize a state from input deck.
    template <class Props, class State>
    void initStateFromDeck(const UnstructuredGrid& grid,
                           const Props& props,
                           const EclipseGridParser& deck,
                           const double gravity,
                           State& state)
    {
        initStateFromDeck(grid.number_of_cells, grid.global_cell,
                          grid.number_of_faces, UgGridHelpers::faceCells(grid),
                          grid.face_centroids, grid.cell_centroids,
                          grid.dimensions, props, deck,
                          gravity, state);
    }

    template <class FaceCells, class FCI, class CCI, class Props, class State>
    void initStateFromDeck(int number_of_cells,
                           const int* global_cell,
                           int number_of_faces,
                           FaceCells face_cells,
                           FCI begin_face_centroids,
                           CCI begin_cell_centroids,
                           int dimensions,
                           const Props& props,
                           const EclipseGridParser& deck,
                           const double gravity,
                           State& state)
    {
        const int num_phases = props.numPhases();
        const PhaseUsage pu = phaseUsageFromDeck(deck);
        if (num_phases != pu.num_phases) {
            OPM_THROW(std::runtime_error, "initStateFromDeck():  user specified property object with " << num_phases << " phases, "
                  "found " << pu.num_phases << " phases in deck.");
        }
        state.init(number_of_cells, number_of_faces, num_phases);
        if (deck.hasField("EQUIL")) {
            if (num_phases != 2) {
                OPM_THROW(std::runtime_error, "initStateFromDeck(): EQUIL-based init currently handling only two-phase scenarios.");
            }
            if (pu.phase_used[BlackoilPhases::Vapour]) {
                OPM_THROW(std::runtime_error, "initStateFromDeck(): EQUIL-based init currently handling only oil-water scenario (no gas).");
            }
            // Set saturations depending on oil-water contact.
            const EQUIL& equil= deck.getEQUIL();
            if (equil.equil.size() != 1) {
                OPM_THROW(std::runtime_error, "initStateFromDeck(): No region support yet.");
            }
            const double woc = equil.equil[0].water_oil_contact_depth_;
            initWaterOilContact(number_of_cells, begin_cell_centroids, dimensions,
                                props, woc, WaterBelow, state);
            // Set pressure depending on densities and depths.
            const double datum_z = equil.equil[0].datum_depth_;
            const double datum_p = equil.equil[0].datum_depth_pressure_;
            initHydrostaticPressure(number_of_cells, begin_cell_centroids, dimensions,
                                    props, woc, gravity, datum_z, datum_p, state);
        } else if (deck.hasField("PRESSURE")) {
            // Set saturations from SWAT/SGAS, pressure from PRESSURE.
            std::vector<double>& s = state.saturation();
            std::vector<double>& p = state.pressure();
            const std::vector<double>& p_deck = deck.getFloatingPointValue("PRESSURE");
            if (num_phases == 2) {
                if (!pu.phase_used[BlackoilPhases::Aqua]) {
                    // oil-gas: we require SGAS
                    if (!deck.hasField("SGAS")) {
                        OPM_THROW(std::runtime_error, "initStateFromDeck(): missing SGAS keyword in 2-phase init");
                    }
                    const std::vector<double>& sg_deck = deck.getFloatingPointValue("SGAS");
                    const int gpos = pu.phase_pos[BlackoilPhases::Vapour];
                    const int opos = pu.phase_pos[BlackoilPhases::Liquid];
                    for (int c = 0; c < number_of_cells; ++c) {
                        int c_deck = (global_cell == NULL) ? c : global_cell[c];
                        s[2*c + gpos] = sg_deck[c_deck];
                        s[2*c + opos] = 1.0 - sg_deck[c_deck];
                        p[c] = p_deck[c_deck];
                    }
                } else {
                    // water-oil or water-gas: we require SWAT
                    if (!deck.hasField("SWAT")) {
                        OPM_THROW(std::runtime_error, "initStateFromDeck(): missing SWAT keyword in 2-phase init");
                    }
                    const std::vector<double>& sw_deck = deck.getFloatingPointValue("SWAT");
                    const int wpos = pu.phase_pos[BlackoilPhases::Aqua];
                    const int nwpos = (wpos + 1) % 2;
                    for (int c = 0; c < number_of_cells; ++c) {
                        int c_deck = (global_cell == NULL) ? c : global_cell[c];
                        s[2*c + wpos] = sw_deck[c_deck];
                        s[2*c + nwpos] = 1.0 - sw_deck[c_deck];
                        p[c] = p_deck[c_deck];
                    }
                }
            } else if (num_phases == 3) {
                const bool has_swat_sgas = deck.hasField("SWAT") && deck.hasField("SGAS");
                if (!has_swat_sgas) {
                    OPM_THROW(std::runtime_error, "initStateFromDeck(): missing SGAS or SWAT keyword in 3-phase init.");
                }
                const int wpos = pu.phase_pos[BlackoilPhases::Aqua];
                const int gpos = pu.phase_pos[BlackoilPhases::Vapour];
                const int opos = pu.phase_pos[BlackoilPhases::Liquid];
                const std::vector<double>& sw_deck = deck.getFloatingPointValue("SWAT");
                const std::vector<double>& sg_deck = deck.getFloatingPointValue("SGAS");
                for (int c = 0; c < number_of_cells; ++c) {
                    int c_deck = (global_cell == NULL) ? c : global_cell[c];
                    s[3*c + wpos] = sw_deck[c_deck];
                    s[3*c + opos] = 1.0 - (sw_deck[c_deck] + sg_deck[c_deck]);
                    s[3*c + gpos] = sg_deck[c_deck];
                    p[c] = p_deck[c_deck];
                 }
            } else {
                OPM_THROW(std::runtime_error, "initStateFromDeck(): init with SWAT etc. only available with 2 or 3 phases.");
            }
        } else {
            OPM_THROW(std::runtime_error, "initStateFromDeck(): we must either have EQUIL, or PRESSURE and SWAT/SOIL/SGAS.");
        }

        // Finally, init face pressures.
        initFacePressure(dimensions, number_of_faces, face_cells, begin_face_centroids,
                         begin_cell_centroids, state);
    }

    /// Initialize surface volume from pressure and saturation by z = As.
    /// Here  saturation is used as an intial guess for z in the
    /// computation of A.
    template <class Props, class State>
    void initBlackoilSurfvol(const UnstructuredGrid& grid,
                             const Props& props,
                             State& state)
    {
        initBlackoilSurfvol(grid.number_of_cells, props, state);
    }

    template <class Props, class State>
    void initBlackoilSurfvol(int number_of_cells,
                             const Props& props,
                             State& state)
    {
        state.surfacevol() = state.saturation();
        const int np = props.numPhases();
        const int nc = number_of_cells;
        std::vector<double> allA(nc*np*np);
        std::vector<int> allcells(nc);
        for (int c = 0; c < nc; ++c) {
            allcells[c] = c;
        }
        // Assuming that using the saturation as z argument here does not change
        // the outcome. This is not guaranteed unless we have only a single phase
        // per cell.
        props.matrix(nc, &state.pressure()[0], &state.surfacevol()[0], &allcells[0], &allA[0], 0);
        for (int c = 0; c < nc; ++c) {
            // Using z = As
            double* z = &state.surfacevol()[c*np];
            const double* s = &state.saturation()[c*np];
            const double* A = &allA[c*np*np];

            for (int row = 0; row < np; ++row) { z[row] = 0.0; }

            for (int col = 0; col < np; ++col) {
                for (int row = 0; row < np; ++row) {
                    // Recall: A has column-major ordering.
                    z[row] += A[row + np*col]*s[col];
                }
            }
        }
    }   
    /// Initialize surface volume from pressure and saturation by z = As.
    /// Here the solution gas/oil ratio or vapor oil/gas ratio is used to
    /// compute an intial z for the computation of A.
    template <class Props, class State>
    void initBlackoilSurfvolUsingRSorRV(const UnstructuredGrid& grid,
                             const Props& props,
                             State& state)
    {
        initBlackoilSurfvolUsingRSorRV(grid.number_of_cells, props, state);
    }

    template <class Props, class State>
    void initBlackoilSurfvolUsingRSorRV(int number_of_cells,
                                        const Props& props,
                                        State& state)
    {
        
        if (props.numPhases() != 3) {
            OPM_THROW(std::runtime_error, "initBlackoilSurfvol() is only supported in three-phase simulations.");
        }
        const std::vector<double>& rs = state.gasoilratio();
        const std::vector<double>& rv = state.rv();

        //make input for computation of the A matrix
        state.surfacevol() = state.saturation();
        const PhaseUsage pu = props.phaseUsage();

        const int np = props.numPhases();

        std::vector<double> allA_a(number_of_cells*np*np);
        std::vector<double> allA_l(number_of_cells*np*np);
        std::vector<double> allA_v(number_of_cells*np*np);

        std::vector<int> allcells(number_of_cells);
        std::vector<double> z_init(number_of_cells*np);
        std::fill(z_init.begin(),z_init.end(),0.0);

        for (int c = 0; c < number_of_cells; ++c) {
            allcells[c] = c;
        }

        std::vector<double> capPressures(number_of_cells*np);
        props.capPress(number_of_cells,&state.saturation()[0],&allcells[0],&capPressures[0],NULL);

        std::vector<double> Pw(number_of_cells);
        std::vector<double> Pg(number_of_cells);

        for (int c = 0; c < number_of_cells; ++c){
            Pw[c] = state.pressure()[c] + capPressures[c*np + BlackoilPhases::Aqua];
            Pg[c] = state.pressure()[c] + capPressures[c*np + BlackoilPhases::Vapour];
        }


        double z_tmp;

        // Water phase
        if(pu.phase_used[BlackoilPhases::Aqua])
           for (int c = 0; c <  number_of_cells ; ++c){
               for (int p = 0; p < np ; ++p){
                   if (p == BlackoilPhases::Aqua)
                       z_tmp = 1;
                   else
                       z_tmp = 0;

                   z_init[c*np + p] = z_tmp;
               }
           }
        props.matrix(number_of_cells, &state.pressure()[0], &z_init[0], &allcells[0], &allA_a[0], 0);

        // Liquid phase
        if(pu.phase_used[BlackoilPhases::Liquid]){
            for (int c = 0; c <  number_of_cells ; ++c){
                for (int p = 0; p < np ; ++p){
                     if(p == BlackoilPhases::Vapour){
                         if(state.saturation()[np*c + p] > 0)
                             z_tmp = 1e10;
                         else
                             z_tmp = rs[c];
                     }
                     else if(p == BlackoilPhases::Liquid)
                         z_tmp = 1;
                     else
                         z_tmp = 0;

                     z_init[c*np + p] = z_tmp;

                }
            }
        }
        props.matrix(number_of_cells, &state.pressure()[0], &z_init[0], &allcells[0], &allA_l[0], 0);

        if(pu.phase_used[BlackoilPhases::Vapour]){
            for (int c = 0; c <  number_of_cells ; ++c){
                for (int p = 0; p < np ; ++p){
                     if(p == BlackoilPhases::Liquid){
                         if(state.saturation()[np*c + p] > 0)
                             z_tmp = 1e10;
                         else
                             z_tmp = rv[c];
                     }
                     else if(p == BlackoilPhases::Vapour)
                         z_tmp = 1;
                     else
                         z_tmp = 0;

                     z_init[c*np + p] = z_tmp;

                }
            }
        }
        props.matrix(number_of_cells, &state.pressure()[0], &z_init[0], &allcells[0], &allA_v[0], 0);

        for (int c = 0; c < number_of_cells; ++c) {
            // Using z = As
            double* z = &state.surfacevol()[c*np];
            const double* s = &state.saturation()[c*np];
            const double* A_a = &allA_a[c*np*np];
            const double* A_l = &allA_l[c*np*np];
            const double* A_v = &allA_v[c*np*np];

            for (int row = 0; row < np; ++row) { z[row] = 0.0; }

            for (int col = 0; col < np; ++col) {
                z[0] += A_a[0 + np*col]*s[col];
                z[1] += A_l[1 + np*col]*s[col];
                z[2] += A_v[2 + np*col]*s[col];

            }
            double ztmp = z[2];
            z[2] += z[1]*rs[c];
            z[1] += ztmp*rv[c];

        }
    }


    /// Initialize a blackoil state from input deck.
    template <class Props, class State>
    void initBlackoilStateFromDeck(const UnstructuredGrid& grid,
                                   const Props& props,
                                   const EclipseGridParser& deck,
                                   const double gravity,
                                   State& state)
    {
        initBlackoilStateFromDeck(grid.number_of_cells, grid.global_cell,
                                  grid.number_of_faces, UgGridHelpers::faceCells(grid),
                                  grid.face_centroids, grid.cell_centroids,grid.dimensions,
                                  props, deck, gravity, state);
    }

    template <class FaceCells, class FCI, class CCI, class Props, class State>
    void initBlackoilStateFromDeck(int number_of_cells,
                                   const int* global_cell,
                                   int number_of_faces,
                                   FaceCells face_cells,
                                   FCI begin_face_centroids,
                                   CCI begin_cell_centroids,
                                   int dimensions,
                                   const Props& props,
                                   const EclipseGridParser& deck,
                                   const double gravity,
                                   State& state)
    {
        initStateFromDeck(number_of_cells, global_cell, number_of_faces,
                          face_cells, begin_face_centroids, begin_cell_centroids,
                          dimensions, props, deck, gravity, state);
        if (deck.hasField("RS")) {
            const std::vector<double>& rs_deck = deck.getFloatingPointValue("RS");
            for (int c = 0; c < number_of_cells; ++c) {
                int c_deck = (global_cell == NULL) ? c : global_cell[c];
                state.gasoilratio()[c] = rs_deck[c_deck];
            }
            initBlackoilSurfvolUsingRSorRV(number_of_cells, props, state);
            computeSaturation(props,state);
        } else if (deck.hasField("RV")){
            const std::vector<double>& rv_deck = deck.getFloatingPointValue("RV");
            for (int c = 0; c < number_of_cells; ++c) {
                int c_deck = (global_cell == NULL) ? c : global_cell[c];
                state.rv()[c] = rv_deck[c_deck];
            }
            initBlackoilSurfvolUsingRSorRV(number_of_cells, props, state);
            computeSaturation(props,state);
        }

        else {
            OPM_THROW(std::runtime_error, "Temporarily, we require the RS or the RV field.");
        }
    }

    /// Initialize a blackoil state from input deck.
    template <class Props, class State>
    void initBlackoilStateFromDeck(const UnstructuredGrid& grid,
                                   const Props& props,
                                   Opm::DeckConstPtr newParserDeck,
                                   const double gravity,
                                   State& state)
    {
        initStateFromDeck(grid, props, newParserDeck, gravity, state);
        if (newParserDeck->hasKeyword("RS")) {
            const std::vector<double>& rs_deck = newParserDeck->getKeyword("RS")->getSIDoubleData();
            const int num_cells = grid.number_of_cells;
            for (int c = 0; c < num_cells; ++c) {
                int c_deck = (grid.global_cell == NULL) ? c : grid.global_cell[c];
                state.gasoilratio()[c] = rs_deck[c_deck];
            }
            initBlackoilSurfvolUsingRSorRV(grid, props, state);
            computeSaturation(props,state);
        } else if (newParserDeck->hasKeyword("RV")){
            const std::vector<double>& rv_deck = newParserDeck->getKeyword("RV")->getSIDoubleData();
            const int num_cells = grid.number_of_cells;
            for (int c = 0; c < num_cells; ++c) {
                int c_deck = (grid.global_cell == NULL) ? c : grid.global_cell[c];
                state.rv()[c] = rv_deck[c_deck];
            }
            initBlackoilSurfvolUsingRSorRV(grid, props, state);
            computeSaturation(props,state);
        }
        else {
            OPM_THROW(std::runtime_error, "Temporarily, we require the RS or the RV field.");
        }
    }


} // namespace Opm


#endif // OPM_INITSTATE_IMPL_HEADER_INCLUDED
