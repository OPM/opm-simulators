/*
  Copyright 2014 SINTEF ICT, Applied Mathematics.
  Copyright 2017 IRIS

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


#ifdef HAVE_CONFIG_H
#include "config.h"
#endif // HAVE_CONFIG_H

#include <opm/grid/UnstructuredGrid.h>
#include <opm/grid/GridManager.hpp>
#include <opm/common/ErrorMacros.hpp>
#include <opm/core/simulator/initStateEquil.hpp>
#include <opm/common/utility/parameters/ParameterGroup.hpp>
#include <opm/core/props/BlackoilPhases.hpp>
#include <opm/core/props/phaseUsageFromDeck.hpp>
#include <opm/core/simulator/BlackoilState.hpp>
#include <opm/grid/utility/compressedToCartesian.hpp>

#include <opm/parser/eclipse/Parser/ParseContext.hpp>
#include <opm/parser/eclipse/Parser/Parser.hpp>
#include <opm/parser/eclipse/Deck/Deck.hpp>

#include <opm/material/fluidmatrixinteractions/EclMaterialLawManager.hpp>

#include <boost/filesystem.hpp>

#include <fstream>

namespace
{
    void warnIfUnusedParams(const Opm::ParameterGroup& param)
    {
        if (param.anyUnused()) {
            std::cout << "--------------------   Unused parameters:   --------------------\n";
            param.displayUsage();
            std::cout << "----------------------------------------------------------------" << std::endl;
        }
    }

    void outputData(const std::string& output_dir,
                    const std::string& name,
                    const std::vector<double>& data)
    {
        std::ostringstream fname;
        fname << output_dir << "/" << name;
        boost::filesystem::path fpath = fname.str();
        try {
            create_directories(fpath);
        }
        catch (...) {
            OPM_THROW(std::runtime_error, "Creating directories failed: " << fpath);
        }
        fname << "/" << "initial.txt";
        std::ofstream file(fname.str().c_str());
        if (!file) {
            OPM_THROW(std::runtime_error, "Failed to open " << fname.str());
        }
        std::copy(data.begin(), data.end(), std::ostream_iterator<double>(file, "\n"));
    }

    /// Convert saturations from a vector of individual phase saturation vectors
    /// to an interleaved format where all values for a given cell come before all
    /// values for the next cell, all in a single vector.
    template <class FluidSystem>
    void convertSats(std::vector<double>& sat_interleaved, const std::vector< std::vector<double> >& sat, const Opm::PhaseUsage& pu)
    {
        assert(sat.size() == 3);
        const auto nc = sat[0].size();
        const auto np = sat_interleaved.size() / nc;
        for (size_t c = 0; c < nc; ++c) {
            if ( FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx)) {
                const int opos = pu.phase_pos[Opm::BlackoilPhases::Liquid];
                const std::vector<double>& sat_p = sat[ FluidSystem::oilPhaseIdx];
                sat_interleaved[np*c + opos] = sat_p[c];
            }
            if ( FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)) {
                const int wpos = pu.phase_pos[Opm::BlackoilPhases::Aqua];
                const std::vector<double>& sat_p = sat[ FluidSystem::waterPhaseIdx];
                sat_interleaved[np*c + wpos] = sat_p[c];
            }
            if ( FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)) {
                const int gpos = pu.phase_pos[Opm::BlackoilPhases::Vapour];
                const std::vector<double>& sat_p = sat[ FluidSystem::gasPhaseIdx];
                sat_interleaved[np*c + gpos] = sat_p[c];
            }
        }
    }

} // anon namespace

// ----------------- Main program -----------------
int
main(int argc, char** argv)
try
{
    using namespace Opm;

    // Setup.
    ParameterGroup param(argc, argv);
    std::cout << "---------------    Reading parameters     ---------------" << std::endl;
    const std::string deck_filename = param.get<std::string>("deck_filename");
    Opm::ParseContext parseContext;
    Opm::Parser parser;
    const Opm::Deck& deck = parser.parseFile(deck_filename , parseContext);
    const Opm::EclipseState eclipseState(deck, parseContext);
    const double grav = param.getDefault("gravity", unit::gravity);
    GridManager gm(eclipseState.getInputGrid());
    const UnstructuredGrid& grid = *gm.c_grid();
    warnIfUnusedParams(param);

    // Create material law manager.
    std::vector<int> compressedToCartesianIdx
        = Opm::compressedToCartesian(grid.number_of_cells, grid.global_cell);

    typedef BlackOilFluidSystem<double> FluidSystem;

    // Forward declaring the MaterialLawManager template.
    typedef Opm::ThreePhaseMaterialTraits<double,
    /*wettingPhaseIdx=*/FluidSystem::waterPhaseIdx,
    /*nonWettingPhaseIdx=*/FluidSystem::oilPhaseIdx,
    /*gasPhaseIdx=*/FluidSystem::gasPhaseIdx> MaterialTraits;
    typedef Opm::EclMaterialLawManager<MaterialTraits> MaterialLawManager;

    MaterialLawManager materialLawManager = MaterialLawManager();
    materialLawManager.initFromDeck(deck, eclipseState, compressedToCartesianIdx);

    // Initialisation.
    //initBlackoilSurfvolUsingRSorRV(UgGridHelpers::numCells(grid), props, state);
    BlackoilState state( UgGridHelpers::numCells(grid) , UgGridHelpers::numFaces(grid), 3);
    FluidSystem::initFromDeck(deck, eclipseState);
    PhaseUsage pu = phaseUsageFromDeck(deck);

    typedef EQUIL::DeckDependent::InitialStateComputer<FluidSystem> ISC;

    ISC isc(materialLawManager, eclipseState, grid, grav);

    const bool oil = FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx);
    const int oilpos = FluidSystem::oilPhaseIdx;
    const int waterpos = FluidSystem::waterPhaseIdx;
    const int ref_phase = oil ? oilpos : waterpos;

    state.pressure() = isc.press()[ref_phase];
    convertSats<FluidSystem>(state.saturation(), isc.saturation(), pu);
    state.gasoilratio() = isc.rs();
    state.rv() = isc.rv();

    // Output.
    const std::string output_dir = param.getDefault<std::string>("output_dir", "output");
    outputData(output_dir, "pressure", state.pressure());
    outputData(output_dir, "saturation", state.saturation());
    outputData(output_dir, "rs", state.gasoilratio());
    outputData(output_dir, "rv", state.rv());
}
catch (const std::exception& e) {
    std::cerr << "Program threw an exception: " << e.what() << "\n";
    throw;
}
