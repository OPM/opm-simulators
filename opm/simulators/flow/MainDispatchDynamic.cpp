/*
  Copyright 2013, 2014, 2015 SINTEF ICT, Applied Mathematics.
  Copyright 2014 Dr. Blatt - HPC-Simulation-Software & Services
  Copyright 2015 IRIS AS
  Copyright 2014 STATOIL ASA.

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

#include <config.h>

#include <opm/simulators/flow/Main.hpp>

#include <flow/flow_biofilm.hpp>
#include <flow/flow_blackoil.hpp>
#include <flow/flow_blackoil_legacyassembly.hpp>
#include <flow/flow_brine.hpp>
#include <flow/flow_brine_precsalt_vapwat.hpp>
#include <flow/flow_brine_saltprecipitation.hpp>
#include <flow/flow_energy.hpp>
#include <flow/flow_extbo.hpp>
#include <flow/flow_foam.hpp>
#include <flow/flow_gasoil.hpp>
#include <flow/flow_gasoil_energy.hpp>
#include <flow/flow_gasoildiffuse.hpp>
#include <flow/flow_gaswater.hpp>
#include <flow/flow_gaswater_brine.hpp>
#include <flow/flow_gaswater_dissolution.hpp>
#include <flow/flow_gaswater_dissolution_diffuse.hpp>
#include <flow/flow_gaswater_energy.hpp>
#include <flow/flow_gaswater_saltprec_energy.hpp>
#include <flow/flow_gaswater_saltprec_vapwat.hpp>
#include <flow/flow_gaswater_solvent.hpp>
#include <flow/flow_micp.hpp>
#include <flow/flow_oilwater.hpp>
#include <flow/flow_oilwater_brine.hpp>
#include <flow/flow_oilwater_polymer.hpp>
#include <flow/flow_oilwater_polymer_injectivity.hpp>
#include <flow/flow_onephase.hpp>
#include <flow/flow_onephase_energy.hpp>
#include <flow/flow_polymer.hpp>
#include <flow/flow_solvent.hpp>
#include <flow/flow_solvent_foam.hpp>

#include <cstdlib>
#include <iostream>

// ---------------------------------------------------------------------------
// Implementation of dispatchDynamic_()
// ---------------------------------------------------------------------------

int Opm::Main::dispatchDynamic_()
{
    const auto& rspec = this->eclipseState_->runspec();
    const auto& phases = rspec.phases();

    this->setupVanguard();

    // run the actual simulator
    //
    // TODO: make sure that no illegal combinations like thermal and
    //       twophase are requested.
    const bool thermal = eclipseState_->getSimulationConfig().isThermal();

    // Single-phase case
    if (rspec.micp()) {
        return this->runMICP(phases);
    }

    // water-only case
    else if (phases.size() == 1 && phases.active(Phase::WATER) && !thermal) {
        return this->runWaterOnly(phases);
    }

    // water-only case with energy
    else if (phases.size() == 2 && phases.active(Phase::WATER) && thermal) {
        return this->runWaterOnlyEnergy(phases);
    }

    // Biofilm case
    else if (rspec.biof()) {
        return this->runBiofilm(phases);
    }

    // Twophase cases
    else if (phases.size() == 2 && !thermal) {
        return this->runTwoPhase(phases);
    }

    // Polymer case
    else if (phases.active(Phase::POLYMER)) {
        return this->runPolymer(phases);
    }

    // Foam case
    else if (phases.active(Phase::FOAM) && !phases.active(Phase::SOLVENT)) {
        return this->runFoam();
    }

    // Solvent case
    else if (phases.active(Phase::SOLVENT)) {
        return this->runSolvent(phases);
    }

    // Brine case
    else if (phases.active(Phase::BRINE) && !thermal) {
        return this->runBrine(phases);
    }

    // Extended BO case
    else if (phases.active(Phase::ZFRACTION)) {
        return this->runExtendedBlackOil();
    }

    // Energy case
    else if (thermal) {
        return this->runThermal(phases);
    }

    // Blackoil case
    else if (phases.size() == 3) {
        return this->runBlackOil();
    }

    else {
        if (outputCout_) {
            std::cerr << "No suitable configuration found, valid are "
                      << "Twophase, polymer, foam, brine, solvent, "
                      << "energy, and blackoil.\n";
        }

        return EXIT_FAILURE;
    }
}

int Opm::Main::runMICP(const Phases& phases)
{
    if (!phases.active(Phase::WATER) || (phases.size() > 2)) {
        if (outputCout_) {
            std::cerr << "No valid configuration is found for MICP simulation, "
                      << "the only valid option is water + MICP\n";
        }

        return EXIT_FAILURE;
    }

    return flowMICPMain(this->argc_,
                        this->argv_,
                        this->outputCout_,
                        this->outputFiles_);
}

int Opm::Main::runTwoPhase(const Phases& phases)
{
    const bool diffusive = eclipseState_->getSimulationConfig().isDiffusive();
    const bool disgasw = eclipseState_->getSimulationConfig().hasDISGASW();
    const bool vapwat = eclipseState_->getSimulationConfig().hasVAPWAT();

    // oil-gas
    if (phases.active(Phase::OIL) && phases.active(Phase::GAS)) {
        if (diffusive) {
            return flowGasOilDiffuseMain(argc_, argv_, outputCout_, outputFiles_);
        }

        return flowGasOilMain(argc_, argv_, outputCout_, outputFiles_);
    }

    // oil-water
    else if (phases.active(Phase::OIL) && phases.active(Phase::WATER)) {
        if (diffusive) {
            if (outputCout_) {
                std::cerr << "The DIFFUSE option is not available for "
                             "the two-phase water/oil model.\n";
            }

            return EXIT_FAILURE;
        }

        return flowOilWaterMain(argc_, argv_, outputCout_, outputFiles_);
    }

    // gas-water
    else if (phases.active(Phase::GAS) && phases.active(Phase::WATER)) {
        if (disgasw || vapwat) {
            if (diffusive) {
                return flowGasWaterDissolutionDiffuseMain(argc_, argv_,
                                                          outputCout_,
                                                          outputFiles_);
            }

            return flowGasWaterDissolutionMain(argc_, argv_, outputCout_, outputFiles_);
        }

        if (diffusive) {
            if (outputCout_) {
                std::cerr << "The DIFFUSE option is not available for "
                             "the two-phase gas/water model without "
                             "disgasw or vapwat.\n";
            }
            return EXIT_FAILURE;
        }

        return flowGasWaterMain(argc_, argv_, outputCout_, outputFiles_);
    }
    else {
        if (outputCout_) {
            std::cerr << "No suitable configuration found, valid "
                         "are Twophase (oilwater, oilgas and gaswater), "
                         "polymer, solvent, or blackoil.\n";
        }

        return EXIT_FAILURE;
    }
}


int Opm::Main::runBiofilm(const Phases& phases)
    {
        if (!(phases.active(Phase::WATER) && phases.active(Phase::GAS)) || (phases.size() != 2)) {
            if (outputCout_) {
                std::cerr << "Biofilm option can only be used for two-phase water/gas "
                          << "model (i.e. in combination with WATER and GAS)." << std::endl;
            }

            return EXIT_FAILURE;
        }
        return flowBiofilmMain(this->argc_,
                        this->argv_,
                        this->outputCout_,
                        this->outputFiles_);
    }

int Opm::Main::runPolymer(const Phases& phases)
{
    if (! phases.active(Phase::WATER)) {
        if (outputCout_) {
            std::cerr << "No valid configuration is found for polymer "
                         "simulation, valid options include "
                         "oilwater + polymer and blackoil + polymer\n";
        }

        return EXIT_FAILURE;
    }

    // Need to track the polymer molecular weight
    // for the injectivity study
    if (phases.active(Phase::POLYMW)) {
        // only oil water two phase for now
        assert (phases.size() == 4);
        return flowOilWaterPolymerInjectivityMain(argc_, argv_, outputCout_, outputFiles_);
    }

    if (phases.size() == 3) { // oil water polymer case
        return flowOilWaterPolymerMain(argc_, argv_, outputCout_, outputFiles_);
    }

    return flowPolymerMain(argc_, argv_, outputCout_, outputFiles_);
}

int Opm::Main::runFoam()
{
    return flowFoamMain(argc_, argv_, outputCout_, outputFiles_);
}

int Opm::Main::runWaterOnly(const Phases& phases)
{
    if (!phases.active(Phase::WATER) || phases.size() != 1) {
        if (outputCout_) {
            std::cerr << "No valid configuration is found for "
                         "water-only simulation, valid options include "
                         "water, water + thermal\n";
        }

        return EXIT_FAILURE;
    }

    return flowWaterOnlyMain(argc_, argv_, outputCout_, outputFiles_);
}

int Opm::Main::runWaterOnlyEnergy(const Phases& phases)
{
    if (!phases.active(Phase::WATER) || phases.size() != 2) {
        if (outputCout_) {
            std::cerr << "No valid configuration is found for water-only "
                         "simulation, valid options include "
                         "water, water + thermal\n";
        }

        return EXIT_FAILURE;
    }

    return flowWaterOnlyEnergyMain(argc_, argv_, outputCout_, outputFiles_);
}

int Opm::Main::runBrine(const Phases& phases)
{
    if (! phases.active(Phase::WATER) || phases.size() == 2) {
        if (outputCout_) {
            std::cerr << "No valid configuration is found for brine "
                         "simulation, valid options include "
                         "oilwater + brine, gaswater + brine "
                         "and blackoil + brine\n";
        }

        return EXIT_FAILURE;
    }

    if (phases.size() == 3) {
        if (phases.active(Phase::OIL)) {
            // oil water brine case
            return flowOilWaterBrineMain(argc_, argv_, outputCout_, outputFiles_);
        }

        if (phases.active(Phase::GAS)) {
            // gas water brine case
            if (eclipseState_->getSimulationConfig().hasPRECSALT() &&
                eclipseState_->getSimulationConfig().hasVAPWAT())
            {
                // Case with water vaporization into gas phase and salt precipitation
                return flowGasWaterSaltprecVapwatMain(argc_, argv_,
                                                      outputCout_,
                                                      outputFiles_);
            }
            else {
                return flowGasWaterBrineMain(argc_, argv_, outputCout_, outputFiles_);
            }
        }
    }
    else if (eclipseState_->getSimulationConfig().hasPRECSALT()) {
        if (eclipseState_->getSimulationConfig().hasVAPWAT()) {
            //case with water vaporization into gas phase and salt precipitation
            return flowBrinePrecsaltVapwatMain(argc_, argv_, outputCout_, outputFiles_);
        }
        else {
            return flowBrineSaltPrecipitationMain(argc_, argv_, outputCout_, outputFiles_);
        }
    }
    else {
        return flowBrineMain(argc_, argv_, outputCout_, outputFiles_);
    }

    return EXIT_FAILURE;
}

int Opm::Main::runSolvent(const Phases& phases)
{
    if (phases.active(Phase::FOAM)) {
        return flowSolventFoamMain(argc_, argv_, outputCout_, outputFiles_);
    }

    // solvent + gas + water
    if (!phases.active(Phase::OIL) &&
        phases.active(Phase::WATER) &&
        phases.active(Phase::GAS))
    {
        return flowGasWaterSolventMain(argc_, argv_, outputCout_, outputFiles_);
    }

    // solvent + gas + water + oil
    if (phases.active(Phase::OIL) &&
        phases.active(Phase::WATER) &&
        phases.active(Phase::GAS))
    {
        return flowSolventMain(argc_, argv_, outputCout_, outputFiles_);
    }

    if (outputCout_) {
        std::cerr << "No valid configuration is found for solvent "
                     "simulation, valid options include "
                     "gas + water + solvent and gas + oil + water + solvent\n";
    }

    return EXIT_FAILURE;
}

int Opm::Main::runExtendedBlackOil()
{
    return flowExtboMain(argc_, argv_, outputCout_, outputFiles_);
}

int Opm::Main::runThermal(const Phases& phases)
{
    // oil-gas-thermal
    if (!phases.active(Phase::WATER) &&
        phases.active(Phase::OIL) &&
        phases.active( Phase::GAS))
    {
        return flowGasOilEnergyMain(argc_, argv_, outputCout_, outputFiles_);
    }

    // water-gas-thermal
    if (!phases.active(Phase::OIL) &&
        phases.active(Phase::WATER) &&
        phases.active(Phase::GAS))
    {
        if (phases.active(Phase::BRINE)) {
            return flowGasWaterSaltprecEnergyMain(argc_, argv_, outputCout_, outputFiles_);
        }

        return flowGasWaterEnergyMain(argc_, argv_, outputCout_, outputFiles_);
    }

    return flowEnergyMain(argc_, argv_, outputCout_, outputFiles_);
}

int Opm::Main::runBlackOil()
{
    if (this->eclipseState_->getSimulationConfig().isDiffusive()) {
        // Use the traditional linearizer, as the TpfaLinearizer does not
        // support the diffusion module yet.
        return flowBlackoilMain(argc_, argv_, outputCout_, outputFiles_);
    }

    return flowBlackoilTpfaMain(argc_, argv_, outputCout_, outputFiles_);
}
