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


#include <opm/core/wells/WellsManager.hpp>
#include <opm/core/io/eclipse/EclipseGridParser.hpp>
#include <opm/core/grid.h>
#include <opm/core/wells.h>
#include <opm/core/well_controls.h>
#include <opm/core/utility/ErrorMacros.hpp>
#include <opm/core/wells/WellCollection.hpp>
#include <opm/core/props/phaseUsageFromDeck.hpp>

#include <opm/parser/eclipse/EclipseState/Schedule/ScheduleEnums.hpp>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <map>
#include <string>
#include <utility>
#include <iostream>


// Helper structs and functions for the implementation.
namespace WellsManagerDetail
{



    namespace ProductionControl
    {
        namespace Details {
            std::map<std::string, Mode>
            init_mode_map() {
                std::map<std::string, Mode> m;

                m.insert(std::make_pair("ORAT", ORAT));
                m.insert(std::make_pair("WRAT", WRAT));
                m.insert(std::make_pair("GRAT", GRAT));
                m.insert(std::make_pair("LRAT", LRAT));
                m.insert(std::make_pair("CRAT", CRAT));
                m.insert(std::make_pair("RESV", RESV));
                m.insert(std::make_pair("BHP" , BHP ));
                m.insert(std::make_pair("THP" , THP ));
                m.insert(std::make_pair("GRUP", GRUP));

                return m;
            }
        } // namespace Details

        Mode mode(const std::string& control)
        {
            static std::map<std::string, Mode>
                mode_map = Details::init_mode_map();

            std::map<std::string, Mode>::iterator
                p = mode_map.find(control);

            if (p != mode_map.end()) {
                return p->second;
            }
            else {
                OPM_THROW(std::runtime_error, "Unknown well control mode = "
                      << control << " in input file");
            }
        }


        Mode mode(Opm::WellProducer::ControlModeEnum controlMode)
        {
            switch( controlMode ) {
            case Opm::WellProducer::ORAT:
                return ORAT;
            case Opm::WellProducer::WRAT:
                return WRAT;
            case Opm::WellProducer::GRAT:
                return GRAT;
            case Opm::WellProducer::LRAT:
                return LRAT;
            case Opm::WellProducer::CRAT:
                return CRAT;
            case Opm::WellProducer::RESV:
                return RESV;
            case Opm::WellProducer::BHP:
                return BHP;
            case Opm::WellProducer::THP:
                return THP;
            case Opm::WellProducer::GRUP:
                return GRUP;
            default:
                throw std::invalid_argument("unhandled enum value");
            }
        }
    } // namespace ProductionControl


    namespace InjectionControl
    {

        namespace Details {
            std::map<std::string, Mode>
            init_mode_map() {
                std::map<std::string, Mode> m;

                m.insert(std::make_pair("RATE", RATE));
                m.insert(std::make_pair("RESV", RESV));
                m.insert(std::make_pair("BHP" , BHP ));
                m.insert(std::make_pair("THP" , THP ));
                m.insert(std::make_pair("GRUP", GRUP));

                return m;
            }
        } // namespace Details

        Mode mode(const std::string& control)
        {
            static std::map<std::string, Mode>
                mode_map = Details::init_mode_map();

            std::map<std::string, Mode>::iterator
                p = mode_map.find(control);

            if (p != mode_map.end()) {
                return p->second;
            }
            else {
                OPM_THROW(std::runtime_error, "Unknown well control mode = "
                      << control << " in input file");
            }
        }

        Mode mode(Opm::WellInjector::ControlModeEnum controlMode)
        {
            switch ( controlMode  ) {
            case Opm::WellInjector::GRUP:
                return GRUP;
            case Opm::WellInjector::RESV:
                return RESV;
            case Opm::WellInjector::RATE:
                return RATE;
            case Opm::WellInjector::THP:
                return THP;
            case Opm::WellInjector::BHP:
                return BHP;
            default:
                throw std::invalid_argument("unhandled enum value");
            }
        }

    } // namespace InjectionControl

    // Use the Peaceman well model to compute well indices.
    // radius is the radius of the well.
    // cubical contains [dx, dy, dz] of the cell.
    // (Note that the well model asumes that each cell is a cuboid).
    // cell_permeability is the permeability tensor of the given cell.
    // returns the well index of the cell.
    double computeWellIndex(const double radius,
                            const std::array<double, 3>& cubical,
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


    /// Construct from existing wells object.
WellsManager::WellsManager(struct Wells* W, bool checkCellExistence)
    : w_(clone_wells(W)), checkCellExistence_(checkCellExistence)
    {
    }

    /// Construct wells from deck.
    WellsManager::WellsManager(const Opm::EclipseStateConstPtr eclipseState,
                               const size_t timeStep,
                               const UnstructuredGrid& grid,
                               const double* permeability,
                               bool checkCellExistence)
        : w_(0), checkCellExistence_(checkCellExistence)
    {
        if (UgGridHelpers::dimensions(grid) != 3) {
            OPM_THROW(std::runtime_error, "We cannot initialize wells from a deck unless the corresponding grid is 3-dimensional.");
        }

        if (eclipseState->getSchedule()->numWells() == 0) {
            OPM_MESSAGE("No wells specified in Schedule section, initializing no wells");
            return;
        }

        std::map<int,int> cartesian_to_compressed;
        setupCompressedToCartesian(UgGridHelpers::globalCell(grid), 
                                   UgGridHelpers::numCells(grid), cartesian_to_compressed);

        // Obtain phase usage data.
        PhaseUsage pu = phaseUsageFromDeck(eclipseState);

        // These data structures will be filled in this constructor,
        // then used to initialize the Wells struct.
        std::vector<std::string> well_names;
        std::vector<WellData> well_data;


        // For easy lookup:
        std::map<std::string, int> well_names_to_index;

        ScheduleConstPtr schedule = eclipseState->getSchedule();
        std::vector<WellConstPtr> wells = schedule->getWells(timeStep);

        well_names.reserve(wells.size());
        well_data.reserve(wells.size());

        createWellsFromSpecs(wells, timeStep, UgGridHelpers::cell2Faces(grid),
                             UgGridHelpers::cartDims(grid),
                             UgGridHelpers::beginFaceCentroids(grid),
                             UgGridHelpers::beginCellCentroids(grid),
                             UgGridHelpers::dimensions(grid),
                             well_names, well_data, well_names_to_index, pu, cartesian_to_compressed, permeability);
        setupWellControls(wells, timeStep, well_names, pu);

        {
            GroupTreeNodeConstPtr fieldNode = eclipseState->getSchedule()->getGroupTree(timeStep)->getNode("FIELD");
            GroupConstPtr fieldGroup = eclipseState->getSchedule()->getGroup(fieldNode->name());
            well_collection_.addField(fieldGroup, timeStep, pu);
            addChildGroups(fieldNode, eclipseState->getSchedule(), timeStep, pu);
        }

        for (auto wellIter = wells.begin(); wellIter != wells.end(); ++wellIter ) {
            well_collection_.addWell((*wellIter), timeStep, pu);
        }

        well_collection_.setWellsPointer(w_);
        well_collection_.applyGroupControls();


        setupGuideRates(wells, timeStep, well_data, well_names_to_index);

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
    }



    /// Construct wells from deck.
    WellsManager::WellsManager(const Opm::EclipseGridParser& deck,
                               const UnstructuredGrid& grid,
                               const double* permeability,
                               bool checkCellExistence)
        : w_(0), checkCellExistence_(checkCellExistence)
    {
        init(deck, grid.number_of_cells, grid.global_cell, grid.cartdims, grid.dimensions,
             grid.cell_centroids, UgGridHelpers::cell2Faces(grid), grid.face_centroids,
             permeability);
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

    void WellsManager::setupCompressedToCartesian(const int* global_cell, int number_of_cells,
                                                  std::map<int,int>& cartesian_to_compressed ) {
        // global_cell is a map from compressed cells to Cartesian grid cells.
        // We must make the inverse lookup.

        if (global_cell) {
            for (int i = 0; i < number_of_cells; ++i) {
                cartesian_to_compressed.insert(std::make_pair(global_cell[i], i));
            }
        }
        else {
            for (int i = 0; i < number_of_cells; ++i) {
                cartesian_to_compressed.insert(std::make_pair(i, i));
            }
        }

    }



    void WellsManager::setupWellControls(std::vector<WellConstPtr>& wells, size_t timeStep,
                                         std::vector<std::string>& well_names, const PhaseUsage& phaseUsage) {
        int well_index = 0;
        for (auto wellIter= wells.begin(); wellIter != wells.end(); ++wellIter) {
            WellConstPtr well = (*wellIter);

            if ( !( well->getStatus( timeStep ) == WellCommon::SHUT || well->getStatus( timeStep ) == WellCommon::OPEN) ) {
                OPM_THROW(std::runtime_error, "Currently we do not support well status " << WellCommon::Status2String(well->getStatus( timeStep )));
            }

            if (well->isInjector(timeStep)) {
                clear_well_controls(well_index, w_);
                int ok = 1;
                int control_pos[5] = { -1, -1, -1, -1, -1 };

                if (well->hasInjectionControl(timeStep , WellInjector::RATE)) {
                    control_pos[WellsManagerDetail::InjectionControl::RATE] = well_controls_get_num(w_->ctrls[well_index]);
                    double distr[3] = { 0.0, 0.0, 0.0 };
                    WellInjector::TypeEnum injectorType = well->getInjectorType(timeStep);

                    if (injectorType == WellInjector::TypeEnum::WATER) {
                        distr[phaseUsage.phase_pos[BlackoilPhases::Aqua]] = 1.0;
                    } else if (injectorType == WellInjector::TypeEnum::OIL) {
                        distr[phaseUsage.phase_pos[BlackoilPhases::Liquid]] = 1.0;
                    } else if (injectorType == WellInjector::TypeEnum::GAS) {
                        distr[phaseUsage.phase_pos[BlackoilPhases::Vapour]] = 1.0;
                    }

                    ok = append_well_controls(SURFACE_RATE,
                                              well->getSurfaceInjectionRate( timeStep ) ,
                                              distr,
                                              well_index,
                                              w_);
                }

                if (ok && well->hasInjectionControl(timeStep , WellInjector::RESV)) {
                    control_pos[WellsManagerDetail::InjectionControl::RESV] = well_controls_get_num(w_->ctrls[well_index]);
                    double distr[3] = { 0.0, 0.0, 0.0 };
                    WellInjector::TypeEnum injectorType = well->getInjectorType(timeStep);

                    if (injectorType == WellInjector::TypeEnum::WATER) {
                        distr[phaseUsage.phase_pos[BlackoilPhases::Aqua]] = 1.0;
                    } else if (injectorType == WellInjector::TypeEnum::OIL) {
                        distr[phaseUsage.phase_pos[BlackoilPhases::Liquid]] = 1.0;
                    } else if (injectorType == WellInjector::TypeEnum::GAS) {
                        distr[phaseUsage.phase_pos[BlackoilPhases::Vapour]] = 1.0;
                    }

                    ok = append_well_controls(RESERVOIR_RATE,
                                              well->getReservoirInjectionRate( timeStep ),
                                              distr,
                                              well_index,
                                              w_);
                }

                if (ok && well->hasInjectionControl(timeStep , WellInjector::BHP)) {
                    control_pos[WellsManagerDetail::InjectionControl::BHP] = well_controls_get_num(w_->ctrls[well_index]);
                    ok = append_well_controls(BHP,
                                              well->getBHPLimit(timeStep),
                                              NULL,
                                              well_index,
                                              w_);
                }

                if (ok && well->hasInjectionControl(timeStep , WellInjector::THP)) {
                    OPM_THROW(std::runtime_error, "We cannot handle THP limit for well " << well_names[well_index]);
                }


                if (!ok) {
                    OPM_THROW(std::runtime_error, "Failure occured appending controls for well " << well_names[well_index]);
                }


                {
                    WellsManagerDetail::InjectionControl::Mode mode = WellsManagerDetail::InjectionControl::mode( well->getInjectorControlMode(timeStep) );
                    int cpos = control_pos[mode];
                    if (cpos == -1 && mode != WellsManagerDetail::InjectionControl::GRUP) {
                        OPM_THROW(std::runtime_error, "Control not specified in well " << well_names[well_index]);
                    }

                    // We need to check if the well is shut or not
                    if (well->getStatus( timeStep ) == WellCommon::SHUT) {
                        cpos = ~cpos;
                    }
                    set_current_control(well_index, cpos, w_);
                }


                // Set well component fraction.
                double cf[3] = { 0.0, 0.0, 0.0 };
                if (well->getInjectorType(timeStep) == WellInjector::WATER) {
                    if (!phaseUsage.phase_used[BlackoilPhases::Aqua]) {
                        OPM_THROW(std::runtime_error, "Water phase not used, yet found water-injecting well.");
                    }
                    cf[phaseUsage.phase_pos[BlackoilPhases::Aqua]] = 1.0;
                } else if (well->getInjectorType(timeStep) == WellInjector::OIL) {
                    if (!phaseUsage.phase_used[BlackoilPhases::Liquid]) {
                        OPM_THROW(std::runtime_error, "Oil phase not used, yet found oil-injecting well.");
                    }
                    cf[phaseUsage.phase_pos[BlackoilPhases::Liquid]] = 1.0;
                } else if (well->getInjectorType(timeStep) == WellInjector::GAS) {
                    if (!phaseUsage.phase_used[BlackoilPhases::Vapour]) {
                        OPM_THROW(std::runtime_error, "Gas phase not used, yet found gas-injecting well.");
                    }
                    cf[phaseUsage.phase_pos[BlackoilPhases::Vapour]] = 1.0;
                }
                std::copy(cf, cf + phaseUsage.num_phases, w_->comp_frac + well_index*phaseUsage.num_phases);

            }

            if (well->isProducer(timeStep)) {
                // Add all controls that are present in well.
                // First we must clear existing controls, in case the
                // current WCONPROD line is modifying earlier controls.
                clear_well_controls(well_index, w_);
                int control_pos[9] = { -1, -1, -1, -1, -1, -1, -1, -1, -1 };
                int ok = 1;
                if (ok && well->hasProductionControl(timeStep , WellProducer::ORAT)) {
                    if (!phaseUsage.phase_used[BlackoilPhases::Liquid]) {
                        OPM_THROW(std::runtime_error, "Oil phase not active and ORAT control specified.");
                    }

                    control_pos[WellsManagerDetail::ProductionControl::ORAT] = well_controls_get_num(w_->ctrls[well_index]);
                    double distr[3] = { 0.0, 0.0, 0.0 };
                    distr[phaseUsage.phase_pos[BlackoilPhases::Liquid]] = 1.0;
                    ok = append_well_controls(SURFACE_RATE,
                                              -well->getOilRate( timeStep ),
                                              distr,
                                               well_index,
                                              w_);
                }

                if (ok && well->hasProductionControl(timeStep , WellProducer::WRAT)) {
                    if (!phaseUsage.phase_used[BlackoilPhases::Aqua]) {
                        OPM_THROW(std::runtime_error, "Water phase not active and WRAT control specified.");
                    }
                    control_pos[WellsManagerDetail::ProductionControl::WRAT] = well_controls_get_num(w_->ctrls[well_index]);
                    double distr[3] = { 0.0, 0.0, 0.0 };
                    distr[phaseUsage.phase_pos[BlackoilPhases::Aqua]] = 1.0;
                    ok = append_well_controls(SURFACE_RATE,
                                              -well->getWaterRate(timeStep),
                                              distr,
                                              well_index,
                                              w_);
                }

                if (ok && well->hasProductionControl(timeStep , WellProducer::GRAT)) {
                    if (!phaseUsage.phase_used[BlackoilPhases::Vapour]) {
                        OPM_THROW(std::runtime_error, "Gas phase not active and GRAT control specified.");
                    }
                    control_pos[WellsManagerDetail::ProductionControl::GRAT] = well_controls_get_num(w_->ctrls[well_index]);
                    double distr[3] = { 0.0, 0.0, 0.0 };
                    distr[phaseUsage.phase_pos[BlackoilPhases::Vapour]] = 1.0;
                    ok = append_well_controls(SURFACE_RATE,
                                              -well->getGasRate( timeStep ),
                                              distr,
                                              well_index,
                                              w_);
                }

                if (ok && well->hasProductionControl(timeStep , WellProducer::LRAT)) {
                    if (!phaseUsage.phase_used[BlackoilPhases::Aqua]) {
                        OPM_THROW(std::runtime_error, "Water phase not active and LRAT control specified.");
                    }
                    if (!phaseUsage.phase_used[BlackoilPhases::Liquid]) {
                        OPM_THROW(std::runtime_error, "Oil phase not active and LRAT control specified.");
                    }
                    control_pos[WellsManagerDetail::ProductionControl::LRAT] = well_controls_get_num(w_->ctrls[well_index]);
                    double distr[3] = { 0.0, 0.0, 0.0 };
                    distr[phaseUsage.phase_pos[BlackoilPhases::Aqua]] = 1.0;
                    distr[phaseUsage.phase_pos[BlackoilPhases::Liquid]] = 1.0;
                    ok = append_well_controls(SURFACE_RATE,
                                              -well->getLiquidRate(timeStep),
                                              distr,
                                              well_index,
                                              w_);
                }

                if (ok && well->hasProductionControl(timeStep , WellProducer::RESV)) {
                    control_pos[WellsManagerDetail::ProductionControl::RESV] = well_controls_get_num(w_->ctrls[well_index]);
                    double distr[3] = { 1.0, 1.0, 1.0 };
                    ok = append_well_controls(RESERVOIR_RATE,
                                              -well->getResVRate(timeStep),
                                              distr,
                                              well_index,
                                              w_);
                }

                if (ok && well->hasProductionControl(timeStep , WellProducer::BHP)) {
                    control_pos[WellsManagerDetail::ProductionControl::BHP] = well_controls_get_num(w_->ctrls[well_index]);
                    ok = append_well_controls(BHP,
                                              well->getBHPLimit( timeStep ) ,
                                              NULL,
                                              well_index,
                                              w_);
                }

                if (ok && well->hasProductionControl(timeStep , WellProducer::THP)) {
                    OPM_THROW(std::runtime_error, "We cannot handle THP limit for well " << well_names[well_index]);
                }

                if (!ok) {
                    OPM_THROW(std::runtime_error, "Failure occured appending controls for well " << well_names[well_index]);
                }

                WellsManagerDetail::ProductionControl::Mode mode = WellsManagerDetail::ProductionControl::mode(well->getProducerControlMode(timeStep));
                int cpos = control_pos[mode];
                if (cpos == -1 && mode != WellsManagerDetail::ProductionControl::GRUP) {
                    OPM_THROW(std::runtime_error, "Control mode type " << mode << " not present in well " << well_names[well_index]);
                }
                // If it's shut, we complement the cpos
                if (well->getStatus(timeStep) == WellCommon::SHUT) {
                    cpos = ~cpos; // So we can easily retrieve the cpos later
                }
                set_current_control(well_index, cpos, w_);
            }
            well_index++;
        }

    }

    void WellsManager::addChildGroups(GroupTreeNodeConstPtr parentNode, ScheduleConstPtr schedule, size_t timeStep, const PhaseUsage& phaseUsage) {
        for (auto childIter = parentNode->begin(); childIter != parentNode->end(); ++childIter) {
            GroupTreeNodeConstPtr childNode = (*childIter).second;
            well_collection_.addGroup(schedule->getGroup(childNode->name()), parentNode->name(), timeStep, phaseUsage);
            addChildGroups(childNode, schedule, timeStep, phaseUsage);
        }
    }




    void WellsManager::setupGuideRates(std::vector<WellConstPtr>& wells, const size_t timeStep, std::vector<WellData>& well_data, std::map<std::string, int>& well_names_to_index)
    {
        for (auto wellIter = wells.begin(); wellIter != wells.end(); ++wellIter ) {
            WellConstPtr well = *wellIter;
            const int wix = well_names_to_index[well->name()];
            WellNode& wellnode = *well_collection_.getLeafNodes()[wix];

            if (well->getGuideRatePhase(timeStep) != GuideRate::UNDEFINED) {
                if (well_data[wix].type == PRODUCER) {
                    wellnode.prodSpec().guide_rate_ = well->getGuideRate(timeStep);
                    if (well->getGuideRatePhase(timeStep) == GuideRate::OIL) {
                        wellnode.prodSpec().guide_rate_type_ = ProductionSpecification::OIL;
                    } else {
                        OPM_THROW(std::runtime_error, "Guide rate type " << GuideRate::GuideRatePhaseEnum2String(well->getGuideRatePhase(timeStep)) << " specified for producer "
                                  << well->name() << " in WGRUPCON, cannot handle.");
                    }
                } else if (well_data[wix].type == INJECTOR) {
                    wellnode.injSpec().guide_rate_ = well->getGuideRate(timeStep);
                    if (well->getGuideRatePhase(timeStep) == GuideRate::RAT) {
                        wellnode.injSpec().guide_rate_type_ = InjectionSpecification::RAT;
                    } else {
                        OPM_THROW(std::runtime_error, "Guide rate type " << GuideRate::GuideRatePhaseEnum2String(well->getGuideRatePhase(timeStep)) << " specified for injector "
                                  << well->name() << " in WGRUPCON, cannot handle.");
                    }
                } else {
                    OPM_THROW(std::runtime_error, "Unknown well type " << well_data[wix].type << " for well " << well->name());
                }
            }
        }
    }

} // namespace Opm
