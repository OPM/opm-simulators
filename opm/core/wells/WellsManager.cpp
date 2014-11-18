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

    // Compute direction permutation corresponding to completion's
    // direction.  First two elements of return value are directions
    // perpendicular to completion while last element is direction
    // along completion.
    inline std::array< std::array<double,3>::size_type, 3 >
    directionIndices(const Opm::CompletionDirection::DirectionEnum direction)
    {
        typedef std::array<double,3>::size_type idx_t;
        typedef std::array<idx_t,3>             permutation;

        switch (direction) {
        case Opm::CompletionDirection::DirectionEnum::X:
            return permutation {{ idx_t(1), idx_t(2), idx_t(0) }};

        case Opm::CompletionDirection::DirectionEnum::Y:
            return permutation {{ idx_t(2), idx_t(0), idx_t(1) }};

        case Opm::CompletionDirection::DirectionEnum::Z:
            return permutation {{ idx_t(0), idx_t(1), idx_t(2) }};
        }
        // All enum values should be handled above. Therefore
        // we should never reach this one. Anyway for the sake
        // of reduced warnings we throw an exception.
        throw std::invalid_argument("unhandled enum value");

    }

    // Permute (diagonal) permeability components according to
    // completion's direction.
    inline std::array<double,3>
    permComponents(const Opm::CompletionDirection::DirectionEnum direction,
                   const double*                                 perm)
    {
        const auto p = directionIndices(direction);
        const std::array<double,3>::size_type d = 3;

        std::array<double,3>
            K = {{ perm[ p[0]*(d + 1) ] ,
                   perm[ p[1]*(d + 1) ] ,
                   perm[ p[2]*(d + 1) ] }};

        return K;
    }

    // Permute cell's geometric extent according to completion's
    // direction.  Honour net-to-gross ratio.
    //
    // Note: 'extent' is intentionally accepted by modifiable value
    // rather than reference-to-const to support NTG manipulation.
    inline std::array<double,3>
    effectiveExtent(const Opm::CompletionDirection::DirectionEnum direction,
                    const double                                  ntg,
                    std::array<double,3>                          extent)
    {
        // Vertical extent affected by net-to-gross ratio.
        extent[2] *= ntg;

        const auto p = directionIndices(direction);

        std::array<double,3>
            D = {{ extent[ p[0] ] ,
                   extent[ p[1] ] ,
                   extent[ p[2] ] }};

        return D;
    }

    // Compute Peaceman's effective radius of single completion.
    inline double
    effectiveRadius(const std::array<double,3>& K,
                    const std::array<double,3>& D)
    {
        const double K01   = K[0] / K[1];
        const double K10   = K[1] / K[0];

        const double D0_sq = D[0] * D[0];
        const double D1_sq = D[1] * D[1];

        const double num = std::sqrt((std::sqrt(K10) * D0_sq) +
                                     (std::sqrt(K01) * D1_sq));
        const double den = std::pow(K01, 0.25) + std::pow(K10, 0.25);

        // Note: Analytic constant 0.28 derived for infintely sized
        // formation with repeating well placement.
        return 0.28 * (num / den);
    }

    // Use the Peaceman well model to compute well indices.
    // radius is the radius of the well.
    // cubical contains [dx, dy, dz] of the cell.
    // (Note that the well model asumes that each cell is a cuboid).
    // cell_permeability is the permeability tensor of the given cell.
    // returns the well index of the cell.
    double
    computeWellIndex(const double                                  radius,
                     const std::array<double, 3>&                  cubical,
                     const double*                                 cell_permeability,
                     const double                                  skin_factor,
                     const Opm::CompletionDirection::DirectionEnum direction,
                     const double                                  ntg)
    {
        const std::array<double,3>& K =
            permComponents(direction, cell_permeability);

        const std::array<double,3>& D =
            effectiveExtent(direction, ntg, cubical);

        const double r0 = effectiveRadius(K, D);
        const double Kh = std::sqrt(K[0] * K[1]) * D[2];

        // Angle of completion exposed to flow.  We assume centre
        // placement so there's complete exposure (= 2\pi).
        const double angle =
            6.2831853071795864769252867665590057683943387987502116419498;

        double rw = radius;
        if (r0 < rw) {
            std::cerr << "Completion radius exceeds effective radius\n"
                      << "Reset to effective\n";

            rw = r0;
        }

        // NOTE: The formula is originally derived and valid for
        // Cartesian grids only.
        return (angle * Kh) / (std::log(r0 / rw) + skin_factor);
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
    WellsManager::WellsManager(struct Wells* W)
        : w_(clone_wells(W))
    {
    }

    /// Construct wells from deck.
    WellsManager::WellsManager(const Opm::EclipseStateConstPtr eclipseState,
                               const size_t timeStep,
                               const UnstructuredGrid& grid,
                               const double* permeability)
        : w_(0)
    {
        init(eclipseState, timeStep, UgGridHelpers::numCells(grid),
             UgGridHelpers::globalCell(grid), UgGridHelpers::cartDims(grid), 
             UgGridHelpers::dimensions(grid), UgGridHelpers::beginCellCentroids(grid),
             UgGridHelpers::cell2Faces(grid), UgGridHelpers::beginFaceCentroids(grid),
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

            if (well->getStatus(timeStep) == WellCommon::STOP) {
                // STOPed wells are added to the well list with the given controll and closed.
                well_controls_stop_well(w_->ctrls[well_index]);
            }

            if (well->getStatus(timeStep) == WellCommon::SHUT) {
                //SHUT wells are not added to the well list
                continue;
            }

            if (well->isInjector(timeStep)) {
                const WellInjectionProperties& injectionProperties = well->getInjectionProperties(timeStep);
                int ok = 1;
                int control_pos[5] = { -1, -1, -1, -1, -1 };
                                
                clear_well_controls(well_index, w_);
                if (injectionProperties.hasInjectionControl(WellInjector::RATE)) {
                    control_pos[WellsManagerDetail::InjectionControl::RATE] = well_controls_get_num(w_->ctrls[well_index]);
                    double distr[3] = { 0.0, 0.0, 0.0 };
                    WellInjector::TypeEnum injectorType = injectionProperties.injectorType;

                    if (injectorType == WellInjector::TypeEnum::WATER) {
                        distr[phaseUsage.phase_pos[BlackoilPhases::Aqua]] = 1.0;
                    } else if (injectorType == WellInjector::TypeEnum::OIL) {
                        distr[phaseUsage.phase_pos[BlackoilPhases::Liquid]] = 1.0;
                    } else if (injectorType == WellInjector::TypeEnum::GAS) {
                        distr[phaseUsage.phase_pos[BlackoilPhases::Vapour]] = 1.0;
                    }

                    ok = append_well_controls(SURFACE_RATE,
                                              injectionProperties.surfaceInjectionRate, 
                                              distr,
                                              well_index,
                                              w_);
                }

                if (ok && injectionProperties.hasInjectionControl(WellInjector::RESV)) {
                    control_pos[WellsManagerDetail::InjectionControl::RESV] = well_controls_get_num(w_->ctrls[well_index]);
                    double distr[3] = { 0.0, 0.0, 0.0 };
                    WellInjector::TypeEnum injectorType = injectionProperties.injectorType;

                    if (injectorType == WellInjector::TypeEnum::WATER) {
                        distr[phaseUsage.phase_pos[BlackoilPhases::Aqua]] = 1.0;
                    } else if (injectorType == WellInjector::TypeEnum::OIL) {
                        distr[phaseUsage.phase_pos[BlackoilPhases::Liquid]] = 1.0;
                    } else if (injectorType == WellInjector::TypeEnum::GAS) {
                        distr[phaseUsage.phase_pos[BlackoilPhases::Vapour]] = 1.0;
                    }

                    ok = append_well_controls(RESERVOIR_RATE,
                                              injectionProperties.reservoirInjectionRate, 
                                              distr,
                                              well_index,
                                              w_);
                }

                if (ok && injectionProperties.hasInjectionControl(WellInjector::BHP)) {
                    control_pos[WellsManagerDetail::InjectionControl::BHP] = well_controls_get_num(w_->ctrls[well_index]);
                    control_pos[WellsManagerDetail::InjectionControl::BHP] = well_controls_get_num(w_->ctrls[well_index]);
                    ok = append_well_controls(BHP,
                                              injectionProperties.BHPLimit, 
                                              NULL,
                                              well_index,
                                              w_);
                }

                if (ok && injectionProperties.hasInjectionControl(WellInjector::THP)) {
                    OPM_THROW(std::runtime_error, "We cannot handle THP limit for well " << well_names[well_index]);
                }

                if (!ok) {
                    OPM_THROW(std::runtime_error, "Failure occured appending controls for well " << well_names[well_index]);
                }

                if (injectionProperties.controlMode != WellInjector::CMODE_UNDEFINED) {
                    WellsManagerDetail::InjectionControl::Mode mode = WellsManagerDetail::InjectionControl::mode( injectionProperties.controlMode );
                    int cpos = control_pos[mode];
                    if (cpos == -1 && mode != WellsManagerDetail::InjectionControl::GRUP) {
                        OPM_THROW(std::runtime_error, "Control not specified in well " << well_names[well_index]);
                    }

                    set_current_control(well_index, cpos, w_);
                }

                // Set well component fraction.
                double cf[3] = { 0.0, 0.0, 0.0 };
                {
                    WellInjector::TypeEnum injectorType = injectionProperties.injectorType;

                    if (injectorType == WellInjector::WATER) {
                        if (!phaseUsage.phase_used[BlackoilPhases::Aqua]) {
                            OPM_THROW(std::runtime_error, "Water phase not used, yet found water-injecting well.");
                        }
                        cf[phaseUsage.phase_pos[BlackoilPhases::Aqua]] = 1.0;
                    } else if (injectorType == WellInjector::OIL) {
                        if (!phaseUsage.phase_used[BlackoilPhases::Liquid]) {
                            OPM_THROW(std::runtime_error, "Oil phase not used, yet found oil-injecting well.");
                        }
                        cf[phaseUsage.phase_pos[BlackoilPhases::Liquid]] = 1.0;
                    } else if (injectorType == WellInjector::GAS) {
                        if (!phaseUsage.phase_used[BlackoilPhases::Vapour]) {
                            OPM_THROW(std::runtime_error, "Gas phase not used, yet found gas-injecting well.");
                        }
                        cf[phaseUsage.phase_pos[BlackoilPhases::Vapour]] = 1.0;
                    }
                    std::copy(cf, cf + phaseUsage.num_phases, w_->comp_frac + well_index*phaseUsage.num_phases);
                }
            }

            if (well->isProducer(timeStep)) {
                // Add all controls that are present in well.
                // First we must clear existing controls, in case the
                // current WCONPROD line is modifying earlier controls.
                const WellProductionProperties& productionProperties = well->getProductionProperties(timeStep);
                int control_pos[9] = { -1, -1, -1, -1, -1, -1, -1, -1, -1 };
                int ok = 1;
                
                clear_well_controls(well_index, w_);
                if (ok && productionProperties.hasProductionControl(WellProducer::ORAT)) {
                    if (!phaseUsage.phase_used[BlackoilPhases::Liquid]) {
                        OPM_THROW(std::runtime_error, "Oil phase not active and ORAT control specified.");
                    }

                    control_pos[WellsManagerDetail::ProductionControl::ORAT] = well_controls_get_num(w_->ctrls[well_index]);
                    double distr[3] = { 0.0, 0.0, 0.0 };
                    distr[phaseUsage.phase_pos[BlackoilPhases::Liquid]] = 1.0;
                    ok = append_well_controls(SURFACE_RATE,
                                              -productionProperties.OilRate, 
                                              distr,
                                               well_index,
                                              w_);
                }

                if (ok && productionProperties.hasProductionControl(WellProducer::WRAT)) {
                    if (!phaseUsage.phase_used[BlackoilPhases::Aqua]) {
                        OPM_THROW(std::runtime_error, "Water phase not active and WRAT control specified.");
                    }
                    control_pos[WellsManagerDetail::ProductionControl::WRAT] = well_controls_get_num(w_->ctrls[well_index]);
                    double distr[3] = { 0.0, 0.0, 0.0 };
                    distr[phaseUsage.phase_pos[BlackoilPhases::Aqua]] = 1.0;
                    ok = append_well_controls(SURFACE_RATE,
                                              -productionProperties.WaterRate, 
                                              distr,
                                              well_index,
                                              w_);
                }

                if (ok && productionProperties.hasProductionControl(WellProducer::GRAT)) {
                    if (!phaseUsage.phase_used[BlackoilPhases::Vapour]) {
                        OPM_THROW(std::runtime_error, "Gas phase not active and GRAT control specified.");
                    }
                    control_pos[WellsManagerDetail::ProductionControl::GRAT] = well_controls_get_num(w_->ctrls[well_index]);
                    double distr[3] = { 0.0, 0.0, 0.0 };
                    distr[phaseUsage.phase_pos[BlackoilPhases::Vapour]] = 1.0;
                    ok = append_well_controls(SURFACE_RATE,
                                              -productionProperties.GasRate,
                                              distr,
                                              well_index,
                                              w_);
                }

                if (ok && productionProperties.hasProductionControl(WellProducer::LRAT)) {
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
                                              -productionProperties.LiquidRate , 
                                              distr,
                                              well_index,
                                              w_);
                }

                if (ok && productionProperties.hasProductionControl(WellProducer::RESV)) {
                    control_pos[WellsManagerDetail::ProductionControl::RESV] = well_controls_get_num(w_->ctrls[well_index]);
                    double distr[3] = { 1.0, 1.0, 1.0 };
                    ok = append_well_controls(RESERVOIR_RATE,
                                              -productionProperties.ResVRate , 
                                              distr,
                                              well_index,
                                              w_);
                }

                if (ok && productionProperties.hasProductionControl(WellProducer::BHP)) {
                    control_pos[WellsManagerDetail::ProductionControl::BHP] = well_controls_get_num(w_->ctrls[well_index]);
                    ok = append_well_controls(BHP,
                                              productionProperties.BHPLimit , 
                                              NULL,
                                              well_index,
                                              w_);
                }

                if (ok && productionProperties.hasProductionControl(WellProducer::THP)) {
                    OPM_THROW(std::runtime_error, "We cannot handle THP limit for well " << well_names[well_index]);
                }

                if (!ok) {
                    OPM_THROW(std::runtime_error, "Failure occured appending controls for well " << well_names[well_index]);
                }

                if (productionProperties.controlMode != WellProducer::CMODE_UNDEFINED) {
                    WellsManagerDetail::ProductionControl::Mode mode = WellsManagerDetail::ProductionControl::mode(productionProperties.controlMode);
                    int cpos = control_pos[mode];
                    if (cpos == -1 && mode != WellsManagerDetail::ProductionControl::GRUP) {
                        OPM_THROW(std::runtime_error, "Control mode type " << mode << " not present in well " << well_names[well_index]);
                    }
                    else {
                        set_current_control(well_index, cpos, w_);
                    }
                }

                // Set well component fraction to match preferred phase for the well.
                double cf[3] = { 0.0, 0.0, 0.0 };
                {
                    switch (well->getPreferredPhase()) {
                    case Phase::WATER:
                        if (!phaseUsage.phase_used[BlackoilPhases::Aqua]) {
                            OPM_THROW(std::runtime_error, "Water phase not used, yet found water-preferring well.");
                        }
                        cf[phaseUsage.phase_pos[BlackoilPhases::Aqua]] = 1.0;
                        break;
                    case Phase::OIL:
                        if (!phaseUsage.phase_used[BlackoilPhases::Liquid]) {
                            OPM_THROW(std::runtime_error, "Oil phase not used, yet found oil-preferring well.");
                        }
                        cf[phaseUsage.phase_pos[BlackoilPhases::Liquid]] = 1.0;
                        break;
                    case Phase::GAS:
                        if (!phaseUsage.phase_used[BlackoilPhases::Vapour]) {
                            OPM_THROW(std::runtime_error, "Gas phase not used, yet found gas-preferring well.");
                        }
                        cf[phaseUsage.phase_pos[BlackoilPhases::Vapour]] = 1.0;
                        break;
                    default:
                        OPM_THROW(std::logic_error, "Unknown preferred phase: " << well->getPreferredPhase());
                    }
                    std::copy(cf, cf + phaseUsage.num_phases, w_->comp_frac + well_index*phaseUsage.num_phases);
                }
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
