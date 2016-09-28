/*
  Copyright 2014 SINTEF ICT, Applied Mathematics.

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

#ifndef OPM_THRESHOLDPRESSURES_HEADER_INCLUDED
#define OPM_THRESHOLDPRESSURES_HEADER_INCLUDED

#include <opm/core/props/BlackoilPropertiesFromDeck.hpp>
#include <opm/core/props/BlackoilPhases.hpp>

#include <vector>
#include <opm/parser/eclipse/EclipseState/SimulationConfig/SimulationConfig.hpp>
#include <opm/parser/eclipse/EclipseState/SimulationConfig/ThresholdPressure.hpp>
#include <opm/parser/eclipse/EclipseState/EclipseState.hpp>
#include <opm/parser/eclipse/Deck/Deck.hpp>
#include <opm/parser/eclipse/EclipseState/Grid/NNC.hpp>


namespace Opm
{
/// \brief Compute the maximum gravity corrected pressure difference of all
///        equilibration regions given a reservoir state.
/// \tparam    Grid           Type of grid object (UnstructuredGrid or CpGrid).
/// \param[out] maxDp         The resulting pressure difference between equilibration regions
/// \param[in] deck           Input deck, EQLOPTS and THPRES are accessed from it.
/// \param[in] eclipseState   Processed eclipse state, EQLNUM is accessed from it.
/// \param[in] grid           The grid to which the thresholds apply.
/// \param[in] initialState   The state of the reservoir
/// \param[in] props          The object which calculates fluid properties
/// \param[in] gravity        The gravity constant
template <class Grid>
void computeMaxDp(std::map<std::pair<int, int>, double>& maxDp,
                  const DeckConstPtr& deck,
                  EclipseStateConstPtr eclipseState,
                  const Grid& grid,
                  const BlackoilState& initialState,
                  const BlackoilPropertiesFromDeck& props,
                  const double gravity)
{

    const PhaseUsage& pu = props.phaseUsage();

    const auto& eqlnum = eclipseState->get3DProperties().getIntGridProperty("EQLNUM");
    const auto& eqlnumData = eqlnum.getData();

    const int numPhases = initialState.numPhases();
    const int numCells = UgGridHelpers::numCells(grid);
    const int numPvtRegions = deck->getKeyword("TABDIMS").getRecord(0).getItem("NTPVT").get< int >(0);

    // retrieve the minimum (residual!?) and the maximum saturations for all cells
    std::vector<double> minSat(numPhases*numCells);
    std::vector<double> maxSat(numPhases*numCells);
    std::vector<int> allCells(numCells);
    for (int cellIdx = 0; cellIdx < numCells; ++cellIdx) {
        allCells[cellIdx] = cellIdx;
    }
    props.satRange(numCells, allCells.data(), minSat.data(), maxSat.data());

    // retrieve the surface densities
    std::vector<std::vector<double> > surfaceDensity(numPvtRegions);
    const auto& densityKw = deck->getKeyword("DENSITY");
    for (int regionIdx = 0; regionIdx < numPvtRegions; ++regionIdx) {
        surfaceDensity[regionIdx].resize(numPhases);

        if (pu.phase_used[BlackoilPhases::Aqua]) {
            const int wpos = pu.phase_pos[BlackoilPhases::Aqua];
            surfaceDensity[regionIdx][wpos] =
                densityKw.getRecord(regionIdx).getItem("WATER").getSIDouble(0);
        }

        if (pu.phase_used[BlackoilPhases::Liquid]) {
            const int opos = pu.phase_pos[BlackoilPhases::Liquid];
            surfaceDensity[regionIdx][opos] =
                densityKw.getRecord(regionIdx).getItem("OIL").getSIDouble(0);
        }

        if (pu.phase_used[BlackoilPhases::Vapour]) {
            const int gpos = pu.phase_pos[BlackoilPhases::Vapour];
            surfaceDensity[regionIdx][gpos] =
                densityKw.getRecord(regionIdx).getItem("GAS").getSIDouble(0);
        }
    }

    // retrieve the PVT region of each cell. note that we need c++ instead of
    // Fortran indices.
    const int* gc = UgGridHelpers::globalCell(grid);
    std::vector<int> pvtRegion(numCells);
    const auto& cartPvtRegion = eclipseState->get3DProperties().getIntGridProperty("PVTNUM").getData();
    for (int cellIdx = 0; cellIdx < numCells; ++cellIdx) {
        const int cartCellIdx = gc ? gc[cellIdx] : cellIdx;
        pvtRegion[cellIdx] = std::max(0, cartPvtRegion[cartCellIdx] - 1);
    }

    // compute the initial "phase presence" of each cell (required to calculate
    // the inverse formation volume factors
    std::vector<PhasePresence> cond(numCells);
    for (int cellIdx = 0; cellIdx < numCells; ++cellIdx) {
        if (pu.phase_used[BlackoilPhases::Aqua]) {
            const double sw = initialState.saturation()[numPhases*cellIdx + pu.phase_pos[BlackoilPhases::Aqua]];
            if (sw > 0.0) {
                cond[cellIdx].setFreeWater();
            }
        }

        if (pu.phase_used[BlackoilPhases::Liquid]) {
            const double so = initialState.saturation()[numPhases*cellIdx + pu.phase_pos[BlackoilPhases::Liquid]];
            if (so > 0.0) {
                cond[cellIdx].setFreeOil();
            }
        }

        if (pu.phase_used[BlackoilPhases::Vapour]) {
            const double sg = initialState.saturation()[numPhases*cellIdx + pu.phase_pos[BlackoilPhases::Vapour]];
            if (sg > 0.0) {
                cond[cellIdx].setFreeGas();
            }
        }
    }

    // calculate the initial fluid densities for the gravity correction.
    std::vector<std::vector<double>> rho(numPhases);
    for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
        rho[phaseIdx].resize(numCells);
    }

    // compute the capillary pressures of the active phases
    std::vector<double> capPress(numCells*numPhases);
    std::vector<int> cellIdxArray(numCells);
    for (int cellIdx = 0; cellIdx < numCells; ++ cellIdx) {
        cellIdxArray[cellIdx] = cellIdx;
    }
    props.capPress(numCells, initialState.saturation().data(), cellIdxArray.data(), capPress.data(), NULL);

    // compute the absolute pressure of each active phase: for some reason, E100
    // defines the capillary pressure for the water phase as p_o - p_w while it
    // uses p_g - p_o for the gas phase. (it would be more consistent to use the
    // oil pressure as reference for both the other phases.) probably this is
    // done to always have a positive number for the capillary pressure (as long
    // as the medium is hydrophilic)
    std::vector<std::vector<double> > phasePressure(numPhases);
    for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
        phasePressure[phaseIdx].resize(numCells);
    }

    for (int cellIdx = 0; cellIdx < numCells; ++ cellIdx) {
        // we currently hard-code the oil phase as the reference phase!
        assert(pu.phase_used[BlackoilPhases::Liquid]);

        const int opos = pu.phase_pos[BlackoilPhases::Liquid];
        phasePressure[opos][cellIdx] = initialState.pressure()[cellIdx];

        if (pu.phase_used[BlackoilPhases::Aqua]) {
            const int wpos = pu.phase_pos[BlackoilPhases::Aqua];
            phasePressure[wpos][cellIdx] =
                initialState.pressure()[cellIdx]
                + (capPress[cellIdx*numPhases + opos] - capPress[cellIdx*numPhases + wpos]);
        }

        if (pu.phase_used[BlackoilPhases::Vapour]) {
            const int gpos = pu.phase_pos[BlackoilPhases::Vapour];
            phasePressure[gpos][cellIdx] =
                initialState.pressure()[cellIdx]
                + (capPress[cellIdx*numPhases + gpos] - capPress[cellIdx*numPhases + opos]);
        }
    }

    // calculate the densities of the active phases for each cell
    if (pu.phase_used[BlackoilPhases::Aqua]) {
        const int wpos = pu.phase_pos[BlackoilPhases::Aqua];
        const auto& pvtw = props.waterPvt();
        for (int cellIdx = 0; cellIdx < numCells; ++ cellIdx) {
            int pvtRegionIdx = pvtRegion[cellIdx];

            double T = initialState.temperature()[cellIdx];
            double p = phasePressure[wpos][cellIdx];
            double b = pvtw.inverseFormationVolumeFactor(pvtRegionIdx, T, p);

            rho[wpos][cellIdx] = surfaceDensity[pvtRegionIdx][wpos]*b;
        }
    }

    if (pu.phase_used[BlackoilPhases::Liquid]) {
        const int opos = pu.phase_pos[BlackoilPhases::Liquid];
        const auto& pvto = props.oilPvt();
        for (int cellIdx = 0; cellIdx < numCells; ++ cellIdx) {
            int pvtRegionIdx = pvtRegion[cellIdx];

            double T = initialState.temperature()[cellIdx];
            double p = phasePressure[opos][cellIdx];
            double Rs = initialState.gasoilratio()[cellIdx];
            double RsSat = pvto.saturatedGasDissolutionFactor(pvtRegionIdx, T, p);

            double b;
            if (Rs >= RsSat) {
                b = pvto.saturatedInverseFormationVolumeFactor(pvtRegionIdx, T, p);
            }
            else {
                b = pvto.inverseFormationVolumeFactor(pvtRegionIdx, T, p, Rs);
            }

            rho[opos][cellIdx] = surfaceDensity[pvtRegionIdx][opos]*b;
            if (pu.phase_used[BlackoilPhases::Vapour]) {
                int gpos = pu.phase_pos[BlackoilPhases::Vapour];
                rho[opos][cellIdx] += surfaceDensity[pvtRegionIdx][gpos]*Rs*b;
            }
        }
    }

    if (pu.phase_used[BlackoilPhases::Vapour]) {
        const int gpos = pu.phase_pos[BlackoilPhases::Vapour];
        const auto& pvtg = props.gasPvt();
        for (int cellIdx = 0; cellIdx < numCells; ++ cellIdx) {
            int pvtRegionIdx = pvtRegion[cellIdx];

            double T = initialState.temperature()[cellIdx];
            double p = phasePressure[gpos][cellIdx];
            double Rv = initialState.rv()[cellIdx];
            double RvSat = pvtg.saturatedOilVaporizationFactor(pvtRegionIdx, T, p);

            double b;
            if (Rv >= RvSat) {
                b = pvtg.saturatedInverseFormationVolumeFactor(pvtRegionIdx, T, p);
            }
            else {
                b = pvtg.inverseFormationVolumeFactor(pvtRegionIdx, T, p, Rv);
            }
            rho[gpos][cellIdx] = surfaceDensity[pvtRegionIdx][gpos]*b;
            if (pu.phase_used[BlackoilPhases::Liquid]) {
                int opos = pu.phase_pos[BlackoilPhases::Liquid];
                rho[gpos][cellIdx] += surfaceDensity[pvtRegionIdx][opos]*Rv*b;
            }
        }
    }

    // Calculate the maximum pressure potential difference between all PVT region
    // transitions of the initial solution.
    const int num_faces = UgGridHelpers::numFaces(grid);
    const auto& fc = UgGridHelpers::faceCells(grid);
    for (int face = 0; face < num_faces; ++face) {
        const int c1 = fc(face, 0);
        const int c2 = fc(face, 1);
        if (c1 < 0 || c2 < 0) {
            // Boundary face, skip this.
            continue;
        }
        const int gc1 = (gc == 0) ? c1 : gc[c1];
        const int gc2 = (gc == 0) ? c2 : gc[c2];
        const int eq1 = eqlnumData[gc1];
        const int eq2 = eqlnumData[gc2];

        if (eq1 == eq2) {
            // not an equilibration region boundary. skip this.
            continue;
        }

        // update the maximum pressure potential difference between the two
        // regions
        const auto barrierId = std::make_pair(std::min(eq1, eq2), std::max(eq1, eq2));
        if (maxDp.count(barrierId) == 0) {
            maxDp[barrierId] = 0.0;
        }

        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            const double z1 = UgGridHelpers::cellCenterDepth(grid, c1);
            const double z2 = UgGridHelpers::cellCenterDepth(grid, c2);
            const double zAvg = (z1 + z2)/2; // average depth

            const double rhoAvg = (rho[phaseIdx][c1] + rho[phaseIdx][c2])/2;

            const double s1 = initialState.saturation()[numPhases*c1 + phaseIdx];
            const double s2 = initialState.saturation()[numPhases*c2 + phaseIdx];

            const double sResid1 = minSat[numPhases*c1 + phaseIdx];
            const double sResid2 = minSat[numPhases*c2 + phaseIdx];

            // compute gravity corrected pressure potentials at the average depth
            const double p1 = phasePressure[phaseIdx][c1] + rhoAvg*gravity*(zAvg - z1);
            const double p2 = phasePressure[phaseIdx][c2] + rhoAvg*gravity*(zAvg - z2);

            if ((p1 > p2 && s1 > sResid1) || (p2 > p1 && s2 > sResid2))
                maxDp[barrierId] = std::max(maxDp[barrierId], std::abs(p1 - p2));
        }
    }
}

    /// \brief Get a vector of pressure thresholds from EclipseState.
    /// This function looks at EQLOPTS, THPRES and EQLNUM to determine
    /// pressure thresholds.  It does not consider the case where the
    /// threshold values are defaulted, in which case they should be
    /// determined from the initial, equilibrated simulation state.
    /// \tparam    Grid           Type of grid object (UnstructuredGrid or CpGrid).
    /// \param[in] deck           Input deck, EQLOPTS and THPRES are accessed from it.
    /// \param[in] eclipseState   Processed eclipse state, EQLNUM is accessed from it.
    /// \param[in] maxDp          The maximum gravity corrected pressure differences between
    ///                           the equilibration regions.
    /// \param[in] grid           The grid to which the thresholds apply.
    /// \return                   A vector of pressure thresholds, one
    ///                           for each face in the grid. A value
    ///                           of zero means no threshold for that
    ///                           particular face. An empty vector is
    ///                           returned if there is no THPRES
    ///                           feature used in the deck.



    template <class Grid>
    std::vector<double> thresholdPressures(const DeckConstPtr& /* deck */,
                                           EclipseStateConstPtr eclipseState,
                                           const Grid& grid,
                                           const std::map<std::pair<int, int>, double>& maxDp)
    {
        const SimulationConfig& simulationConfig = eclipseState->getSimulationConfig();
        std::vector<double> thpres_vals;
        if (simulationConfig.hasThresholdPressure()) {
            std::shared_ptr<const ThresholdPressure> thresholdPressure = simulationConfig.getThresholdPressure();
            const auto& eqlnum = eclipseState->get3DProperties().getIntGridProperty("EQLNUM");
            const auto& eqlnumData = eqlnum.getData();

            // Set threshold pressure values for each cell face.
            const int num_faces = UgGridHelpers::numFaces(grid);
            const auto& fc = UgGridHelpers::faceCells(grid);
            const int* gc = UgGridHelpers::globalCell(grid);
            thpres_vals.resize(num_faces, 0.0);
            for (int face = 0; face < num_faces; ++face) {
                const int c1 = fc(face, 0);
                const int c2 = fc(face, 1);
                if (c1 < 0 || c2 < 0) {
                    // Boundary face, skip it.
                    continue;
                }
                const int gc1 = (gc == 0) ? c1 : gc[c1];
                const int gc2 = (gc == 0) ? c2 : gc[c2];
                const int eq1 = eqlnumData[gc1];
                const int eq2 = eqlnumData[gc2];

                if (thresholdPressure->hasRegionBarrier(eq1,eq2)) {
                    if (thresholdPressure->hasThresholdPressure(eq1,eq2)) {
                        thpres_vals[face] = thresholdPressure->getThresholdPressure(eq1,eq2);
                    }
                    else {
                        // set the threshold pressure for faces of PVT regions where the third item
                        // has been defaulted to the maximum pressure potential difference between
                        // these regions
                        const auto barrierId = std::make_pair(std::min(eq1, eq2), std::max(eq1, eq2));
                        if (maxDp.count(barrierId) > 0)
                            thpres_vals[face] = maxDp.at(barrierId);
                        else
                            thpres_vals[face] = 0.0;
                    }
                }

            }
        }
        return thpres_vals;
    }

    /// \brief Get a vector of pressure thresholds from either EclipseState
    /// or maxDp (for defaulted values) for all Non-neighbour connections (NNCs).
    /// \param[in] nnc            The NNCs,
    /// \param[in] eclipseState   Processed eclipse state, EQLNUM is accessed from it.
    /// \param[in] maxDp          The maximum gravity corrected pressure differences between
    ///                           the equilibration regions.
    /// \return                   A vector of pressure thresholds, one
    ///                           for each NNC in the grid. A value
    ///                           of zero means no threshold for that
    ///                           particular connection. An empty vector is
    ///                           returned if there is no THPRES
    ///                           feature used in the deck.
     std::vector<double> thresholdPressuresNNC(EclipseStateConstPtr eclipseState,
                                               const NNC& nnc,
                                               const std::map<std::pair<int, int>, double>& maxDp)
    {
        const SimulationConfig& simulationConfig = eclipseState->getSimulationConfig();
        std::vector<double> thpres_vals;
        if (simulationConfig.hasThresholdPressure()) {
            std::shared_ptr<const ThresholdPressure> thresholdPressure = simulationConfig.getThresholdPressure();
            const auto& eqlnum = eclipseState->get3DProperties().getIntGridProperty("EQLNUM");
            const auto& eqlnumData = eqlnum.getData();

            // Set values for each NNC

            thpres_vals.resize(nnc.numNNC(), 0.0);
            for (size_t i = 0 ; i < nnc.numNNC(); ++i) {
                const int gc1 = nnc.nncdata()[i].cell1;
                const int gc2 = nnc.nncdata()[i].cell2;
                const int eq1 = eqlnumData[gc1];
                const int eq2 = eqlnumData[gc2];

                if (thresholdPressure->hasRegionBarrier(eq1,eq2)) {
                    if (thresholdPressure->hasThresholdPressure(eq1,eq2)) {
                        thpres_vals[i] = thresholdPressure->getThresholdPressure(eq1,eq2);
                    } else {
                        // set the threshold pressure for NNC of PVT regions where the third item
                        // has been defaulted to the maximum pressure potential difference between
                        // these regions
                        const auto barrierId = std::make_pair(eq1, eq2);
                        thpres_vals[i] = maxDp.at(barrierId);
                    }
                }
            }
        }
        return thpres_vals;
    }
}

#endif // OPM_THRESHOLDPRESSURES_HEADER_INCLUDED
