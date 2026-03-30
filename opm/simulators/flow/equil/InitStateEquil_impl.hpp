// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
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

  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the list of
  copyright holders.
*/
#ifndef OPM_INIT_STATE_EQUIL_IMPL_HPP
#define OPM_INIT_STATE_EQUIL_IMPL_HPP

#include <opm/simulators/flow/equil/InitStateEquilBase_impl.hpp>
#include <opm/simulators/flow/equil/InitStateEquil.hpp>

#include <opm/input/eclipse/EclipseState/Tables/RsvdTable.hpp>
#include <opm/input/eclipse/EclipseState/Tables/RvvdTable.hpp>
#include <opm/input/eclipse/EclipseState/Tables/RvwvdTable.hpp>
#include <opm/input/eclipse/EclipseState/Tables/PbvdTable.hpp>
#include <opm/input/eclipse/EclipseState/Tables/PdvdTable.hpp>

namespace Opm {
namespace EQUIL {

namespace DeckDependent {

template<class FluidSystem,
         class Grid,
         class GridView,
         class ElementMapper,
         class CartesianIndexMapper>
template<class MaterialLawManager>
InitialStateComputer<FluidSystem,
                     Grid,
                     GridView,
                     ElementMapper,
                     CartesianIndexMapper>::
InitialStateComputer(MaterialLawManager& materialLawManager,
                     const EclipseState& eclipseState,
                     const Grid& grid,
                     const GridView& gridView,
                     const CartesianIndexMapper& cartMapper,
                     const Scalar grav,
                     const int num_pressure_points,
                     const bool applySwatInit)
    : Base(grid.size(/*codim=*/0),
           eclipseState.getTableManager().rtemp(),
           cartMapper,
           num_pressure_points)
    , rs_(grid.size(/*codim=*/0))
    , rv_(grid.size(/*codim=*/0))
    , rvw_(grid.size(/*codim=*/0))
{
    //Check for presence of kw SWATINIT
    if (applySwatInit) {
        if (eclipseState.fieldProps().has_double("SWATINIT")) {
            if constexpr (std::is_same_v<Scalar,double>) {
                this->swatInit_ = eclipseState.fieldProps().get_double("SWATINIT");
            } else {
                const auto& input = eclipseState.fieldProps().get_double("SWATINIT");
                this->swatInit_.resize(input.size());
                std::ranges::copy(input, this->swatInit_.begin());
            }
        }
    }

    // Querry cell depth, cell top-bottom.
    // numerical aquifer cells might be specified with different depths.
    const auto& num_aquifers = eclipseState.aquifer().numericalAquifers();
    this->updateCellProps_(gridView, num_aquifers);

    // Get the equilibration records.
    const std::vector<EquilRecord> rec = getEquil(eclipseState);
    const auto& tables = eclipseState.getTableManager();
    // Create (inverse) region mapping.
    const RegionMapping<> eqlmap(equilnum(eclipseState, grid));
    const int invalidRegion = -1;
    this->regionPvtIdx_.resize(rec.size(), invalidRegion);
    this->setRegionPvtIdx(eclipseState, eqlmap);

    // Create Rs functions.
    rsFunc_.reserve(rec.size());

    auto getArray = [](const std::vector<double>& input)
    {
        if constexpr (std::is_same_v<Scalar,double>) {
            return input;
        } else {
            std::vector<Scalar> output;
            output.resize(input.size());
            std::ranges::copy(input, output.begin());
            return output;
        }
    };

    if (FluidSystem::enableDissolvedGas()) {
        for (std::size_t i = 0; i < rec.size(); ++i) {
            if (eqlmap.cells(i).empty()) {
                rsFunc_.push_back(std::shared_ptr<Miscibility::RsVD<FluidSystem>>());
                continue;
            }
            const int pvtIdx = this->regionPvtIdx_[i];
            if (!rec[i].liveOilInitConstantRs()) {
                const TableContainer& rsvdTables = tables.getRsvdTables();
                const TableContainer& pbvdTables = tables.getPbvdTables();
                if (rsvdTables.size() > 0) {
                    const RsvdTable& rsvdTable = rsvdTables.getTable<RsvdTable>(i);
                    auto depthColumn = getArray(rsvdTable.getColumn("DEPTH").vectorCopy());
                    auto rsColumn = getArray(rsvdTable.getColumn("RS").vectorCopy());
                    rsFunc_.push_back(std::make_shared<Miscibility::RsVD<FluidSystem>>(pvtIdx,
                                                                                       depthColumn, rsColumn));
                } else if (pbvdTables.size() > 0) {
                    const PbvdTable& pbvdTable = pbvdTables.getTable<PbvdTable>(i);
                    auto depthColumn = getArray(pbvdTable.getColumn("DEPTH").vectorCopy());
                    auto pbubColumn = getArray(pbvdTable.getColumn("PBUB").vectorCopy());
                    rsFunc_.push_back(std::make_shared<Miscibility::PBVD<FluidSystem>>(pvtIdx,
                                                                                       depthColumn, pbubColumn));

                } else {
                    throw std::runtime_error("Cannot initialise: RSVD or PBVD table not available.");
                }

            }
            else {
                if (rec[i].gasOilContactDepth() != rec[i].datumDepth()) {
                    throw std::runtime_error("Cannot initialise: when no explicit RSVD table is given, \n"
                                             "datum depth must be at the gas-oil-contact. "
                                             "In EQUIL region "+std::to_string(i + 1)+"  (counting from 1), this does not hold.");
                }
                const Scalar pContact = rec[i].datumDepthPressure();
                const Scalar TContact = 273.15 + 20; // standard temperature for now
                rsFunc_.push_back(std::make_shared<Miscibility::RsSatAtContact<FluidSystem>>(pvtIdx, pContact, TContact));
            }
        }
    }
    else {
        for (std::size_t i = 0; i < rec.size(); ++i) {
            rsFunc_.push_back(std::make_shared<Miscibility::NoMixing<Scalar>>());
        }
    }

    rvFunc_.reserve(rec.size());
    if (FluidSystem::enableVaporizedOil()) {
        for (std::size_t i = 0; i < rec.size(); ++i) {
            if (eqlmap.cells(i).empty()) {
                rvFunc_.push_back(std::shared_ptr<Miscibility::RvVD<FluidSystem>>());
                continue;
            }
            const int pvtIdx = this->regionPvtIdx_[i];
            if (!rec[i].wetGasInitConstantRv()) {
                const TableContainer& rvvdTables = tables.getRvvdTables();
                const TableContainer& pdvdTables = tables.getPdvdTables();

                if (rvvdTables.size() > 0) {
                    const RvvdTable& rvvdTable = rvvdTables.getTable<RvvdTable>(i);
                    auto depthColumn = getArray(rvvdTable.getColumn("DEPTH").vectorCopy());
                    auto rvColumn = getArray(rvvdTable.getColumn("RV").vectorCopy());
                    rvFunc_.push_back(std::make_shared<Miscibility::RvVD<FluidSystem>>(pvtIdx,
                                                                                       depthColumn, rvColumn));
                } else if (pdvdTables.size() > 0) {
                    const PdvdTable& pdvdTable = pdvdTables.getTable<PdvdTable>(i);
                    auto depthColumn = getArray(pdvdTable.getColumn("DEPTH").vectorCopy());
                    auto pdewColumn = getArray(pdvdTable.getColumn("PDEW").vectorCopy());
                    rvFunc_.push_back(std::make_shared<Miscibility::PDVD<FluidSystem>>(pvtIdx,
                                                                                       depthColumn, pdewColumn));
                } else {
                    throw std::runtime_error("Cannot initialise: RVVD or PDCD table not available.");
                }
            }
            else {
                if (rec[i].gasOilContactDepth() != rec[i].datumDepth()) {
                    throw std::runtime_error(
                              "Cannot initialise: when no explicit RVVD table is given, \n"
                              "datum depth must be at the gas-oil-contact. "
                              "In EQUIL region "+std::to_string(i + 1)+" (counting from 1), this does not hold.");
                }
                const Scalar pContact = rec[i].datumDepthPressure() + rec[i].gasOilContactCapillaryPressure();
                const Scalar TContact = 273.15 + 20; // standard temperature for now
                rvFunc_.push_back(std::make_shared<Miscibility::RvSatAtContact<FluidSystem>>(pvtIdx,pContact, TContact));
            }
        }
    }
    else {
        for (std::size_t i = 0; i < rec.size(); ++i) {
            rvFunc_.push_back(std::make_shared<Miscibility::NoMixing<Scalar>>());
        }
    }

    rvwFunc_.reserve(rec.size());
    if (FluidSystem::enableVaporizedWater()) {
        for (std::size_t i = 0; i < rec.size(); ++i) {
            if (eqlmap.cells(i).empty()) {
                rvwFunc_.push_back(std::shared_ptr<Miscibility::RvwVD<FluidSystem>>());
                continue;
            }
            const int pvtIdx = this->regionPvtIdx_[i];
            if (!rec[i].humidGasInitConstantRvw()) {
                const TableContainer& rvwvdTables = tables.getRvwvdTables();

                if (rvwvdTables.size() > 0) {
                    const RvwvdTable& rvwvdTable = rvwvdTables.getTable<RvwvdTable>(i);
                    auto depthColumn = getArray(rvwvdTable.getColumn("DEPTH").vectorCopy());
                    auto rvwvdColumn = getArray(rvwvdTable.getColumn("RVWVD").vectorCopy());
                    rvwFunc_.push_back(std::make_shared<Miscibility::RvwVD<FluidSystem>>(pvtIdx,
                                                                                       depthColumn, rvwvdColumn));
                } else {
                    throw std::runtime_error("Cannot initialise: RVWVD table not available.");
                }
            }
            else {
                const auto oilActive = FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx);
                if (oilActive) {
                    if (rec[i].gasOilContactDepth() != rec[i].datumDepth()) {
                        rvwFunc_.push_back(std::make_shared<Miscibility::NoMixing<Scalar>>());
                        const auto msg = "No explicit RVWVD table is given for EQUIL region " + std::to_string(i + 1) +". \n"
                                        "and datum depth is not at the gas-oil-contact. \n"
                                        "Rvw is set to 0.0 in all cells. \n";
                        OpmLog::warning(msg);
                    } else {
                        // pg = po + Pcgo = po + (pg - po)
                        // for gas-condensate with initial no oil zone: water-oil contact depth (OWC) equal gas-oil contact depth (GOC)
                        const Scalar pContact = rec[i].datumDepthPressure() + rec[i].gasOilContactCapillaryPressure();
                        const Scalar TContact = 273.15 + 20; // standard temperature for now
                        rvwFunc_.push_back(std::make_shared<Miscibility::RvwSatAtContact<FluidSystem>>(pvtIdx,pContact, TContact));
                    }
                }
                else {
                     // two-phase gas-water sytem:  water-oil contact depth is taken equal to gas-water contact depth (GWC)
                     // and water-oil capillary pressure (Pcwo) is taken equal to gas-water capillary pressure (Pcgw) at GWC
                     if (rec[i].waterOilContactDepth() != rec[i].datumDepth()) {
                        rvwFunc_.push_back(std::make_shared<Miscibility::NoMixing<Scalar>>());
                        const auto msg = "No explicit RVWVD table is given for EQUIL region " + std::to_string(i + 1) +". \n"
                                         "and datum depth is not at the gas-water-contact. \n"
                                         "Rvw is set to 0.0 in all cells. \n";
                        OpmLog::warning(msg);
                    } else {
                        // pg = pw + Pcgw = pw + (pg - pw)
                        const Scalar pContact = rec[i].datumDepthPressure() + rec[i].waterOilContactCapillaryPressure();
                        const Scalar TContact = 273.15 + 20; // standard temperature for now
                        rvwFunc_.push_back(std::make_shared<Miscibility::RvwSatAtContact<FluidSystem>>(pvtIdx,pContact, TContact));
                    }
                }
            }
        }
    }
    else {
        for (std::size_t i = 0; i < rec.size(); ++i) {
            rvwFunc_.push_back(std::make_shared<Miscibility::NoMixing<Scalar>>());
        }
    }


    // EXTRACT the initial temperature
    this->updateInitialTemperature_(eclipseState, eqlmap);

    // EXTRACT the initial salt concentration
    this->updateInitialSaltConcentration_(eclipseState, eqlmap);

    // EXTRACT the initial salt saturation
    this->updateInitialSaltSaturation_(eclipseState, eqlmap);

    // Compute pressures, saturations, rs and rv factors.
    const auto& comm = grid.comm();
    calcPressSatRsRv(eqlmap, rec, materialLawManager, gridView, comm, grav);

    // modify the pressure and saturation for numerical aquifer cells
    this->applyNumericalAquifers_(gridView, num_aquifers, eclipseState.runspec().co2Storage() || eclipseState.runspec().h2Storage());

    // Modify oil pressure in no-oil regions so that the pressures of present phases can
    // be recovered from the oil pressure and capillary relations.
}

template<class FluidSystem,
         class Grid,
         class GridView,
         class ElementMapper,
         class CartesianIndexMapper>
template<class RMap, class MaterialLawManager, class Comm>
void InitialStateComputer<FluidSystem,
                          Grid,
                          GridView,
                          ElementMapper,
                          CartesianIndexMapper>::
calcPressSatRsRv(const RMap& reg,
                 const std::vector<EquilRecord>& rec,
                 MaterialLawManager& materialLawManager,
                 const GridView& gridView,
                 const Comm& comm,
                 const Scalar grav)
{
    using PhaseSat = Details::PhaseSaturations<
        MaterialLawManager, FluidSystem, EquilReg<Scalar>, typename RMap::CellId
    >;

    auto ptable = Details::PressureTable<FluidSystem, EquilReg<Scalar>>{ grav, this->num_pressure_points_ };
    auto psat   = PhaseSat { materialLawManager, this->swatInit_ };
    auto vspan  = std::array<Scalar, 2>{};

    std::vector<int> regionIsEmpty(rec.size(), 0);
    for (std::size_t r = 0; r < rec.size(); ++r) {
        const auto& cells = reg.cells(r);

        Details::verticalExtent(cells, this->cellZMinMax_, comm, vspan);

        const auto acc = rec[r].initializationTargetAccuracy();
        if (acc > 0) {
            // The grid blocks are treated as being tilted
            // First check if the region has cells
            if (cells.empty()) {
                regionIsEmpty[r] = 1;
                continue;
            }
            const auto eqreg = EquilReg {
                rec[r], this->rsFunc_[r], this->rvFunc_[r], this->rvwFunc_[r],
                this->tempVdTable_[r], this->saltVdTable_[r], this->regionPvtIdx_[r]
            };
            // Ensure contacts are within the span
            vspan[0] = std::min(vspan[0], std::min(eqreg.zgoc(), eqreg.zwoc()));
            vspan[1] = std::max(vspan[1], std::max(eqreg.zgoc(), eqreg.zwoc()));
            ptable.equilibrate(eqreg, vspan);
            // For titled blocks, we can use a simple weightening based on title of the grid
            // this->equilibrateTiltedFaultBlockSimple(cells, eqreg, gridView, acc, ptable, psat);
            this->equilibrateTiltedFaultBlock(cells, eqreg, gridView, acc, ptable, psat);
        }
        else if (acc == 0) {
            if (cells.empty()) {
                regionIsEmpty[r] = 1;
                continue;
            }
            const auto eqreg = EquilReg {
                rec[r], this->rsFunc_[r], this->rvFunc_[r], this->rvwFunc_[r],
                this->tempVdTable_[r], this->saltVdTable_[r], this->regionPvtIdx_[r]
            };
            vspan[0] = std::min(vspan[0], std::min(eqreg.zgoc(), eqreg.zwoc()));
            vspan[1] = std::max(vspan[1], std::max(eqreg.zgoc(), eqreg.zwoc()));
            ptable.equilibrate(eqreg, vspan);
            // Centre-point method
            this->equilibrateCellCentres(cells, eqreg, ptable, psat);
        }
        else if (acc < 0) {
            if (cells.empty()) {
                regionIsEmpty[r] = 1;
                continue;
            }
            const auto eqreg = EquilReg {
                rec[r], this->rsFunc_[r], this->rvFunc_[r], this->rvwFunc_[r],
                this->tempVdTable_[r], this->saltVdTable_[r], this->regionPvtIdx_[r]
            };
            vspan[0] = std::min(vspan[0], std::min(eqreg.zgoc(), eqreg.zwoc()));
            vspan[1] = std::max(vspan[1], std::max(eqreg.zgoc(), eqreg.zwoc()));
            ptable.equilibrate(eqreg, vspan);
            // Horizontal subdivision
            this->equilibrateHorizontal(cells, eqreg, -acc, ptable, psat);
        }
    }
    comm.min(regionIsEmpty.data(),regionIsEmpty.size());
    if (comm.rank() == 0) {
        for (std::size_t r = 0; r < rec.size(); ++r) {
            if (regionIsEmpty[r]) //region is empty on all partitions
                OpmLog::warning("Equilibration region " + std::to_string(r + 1)
                                 + " has no active cells");
        }
    }
}

template<class FluidSystem,
         class Grid,
         class GridView,
         class ElementMapper,
         class CartesianIndexMapper>
template<class CellRange, class EquilibrationMethod>
void InitialStateComputer<FluidSystem,
                          Grid,
                          GridView,
                          ElementMapper,
                          CartesianIndexMapper>::
cellLoop(const CellRange&      cells,
         EquilibrationMethod&& eqmethod)
{
    const auto oilPos = FluidSystem::oilPhaseIdx;
    const auto gasPos = FluidSystem::gasPhaseIdx;
    const auto watPos = FluidSystem::waterPhaseIdx;

    const auto oilActive = FluidSystem::phaseIsActive(oilPos);
    const auto gasActive = FluidSystem::phaseIsActive(gasPos);
    const auto watActive = FluidSystem::phaseIsActive(watPos);

    auto pressures   = Details::PhaseQuantityValue<Scalar>{};
    auto saturations = Details::PhaseQuantityValue<Scalar>{};
    Scalar Rs          = 0.0;
    Scalar Rv          = 0.0;
    Scalar Rvw         = 0.0;

    for (const auto& cell : cells) {
        eqmethod(cell, pressures, saturations, Rs, Rv, Rvw);

        if (oilActive) {
            this->pp_ [oilPos][cell] = pressures.oil;
            this->sat_[oilPos][cell] = saturations.oil;
        }

        if (gasActive) {
            this->pp_ [gasPos][cell] = pressures.gas;
            this->sat_[gasPos][cell] = saturations.gas;
        }

        if (watActive) {
            this->pp_ [watPos][cell] = pressures.water;
            this->sat_[watPos][cell] = saturations.water;
        }

        if (oilActive && gasActive) {
            this->rs_[cell] = Rs;
            this->rv_[cell] = Rv;
        }

        if (watActive && gasActive) {
            this->rvw_[cell] = Rvw;
        }
    }
}

template<class FluidSystem,
         class Grid,
         class GridView,
         class ElementMapper,
         class CartesianIndexMapper>
template<class CellRange, class PressTable, class PhaseSat>
void InitialStateComputer<FluidSystem,
                          Grid,
                          GridView,
                          ElementMapper,
                          CartesianIndexMapper>::
equilibrateCellCentres(const CellRange&         cells,
                       const EquilReg<Scalar>&  eqreg,
                       const PressTable&        ptable,
                       PhaseSat&                psat)
{
    using CellPos = typename PhaseSat::Position;
    using CellID  = std::remove_cv_t<std::remove_reference_t<
        decltype(std::declval<CellPos>().cell)>>;
    this->cellLoop(cells, [this, &eqreg,  &ptable, &psat]
        (const CellID                 cell,
         Details::PhaseQuantityValue<Scalar>& pressures,
         Details::PhaseQuantityValue<Scalar>& saturations,
         Scalar&                      Rs,
         Scalar&                      Rv,
         Scalar&                      Rvw) -> void
    {
        const auto pos = CellPos {
            cell, this->cellCenterDepth_[cell]
        };

        saturations = psat.deriveSaturations(pos, eqreg, ptable);
        pressures   = psat.correctedPhasePressures();

        const auto temp = this->temperature_[cell];

        Rs = eqreg.dissolutionCalculator()
            (pos.depth, pressures.oil, temp, saturations.gas);

        Rv = eqreg.evaporationCalculator()
            (pos.depth, pressures.gas, temp, saturations.oil);

        Rvw = eqreg.waterEvaporationCalculator()
            (pos.depth, pressures.gas, temp, saturations.water);
    });
}

template<class FluidSystem,
         class Grid,
         class GridView,
         class ElementMapper,
         class CartesianIndexMapper>
template<class CellRange, class PressTable, class PhaseSat>
void InitialStateComputer<FluidSystem,
                          Grid,
                          GridView,
                          ElementMapper,
                          CartesianIndexMapper>::
equilibrateHorizontal(const CellRange&        cells,
                      const EquilReg<Scalar>& eqreg,
                      const int               acc,
                      const PressTable&       ptable,
                      PhaseSat&               psat)
{
    using CellPos = typename PhaseSat::Position;
    using CellID  = std::remove_cv_t<std::remove_reference_t<
        decltype(std::declval<CellPos>().cell)>>;

    this->cellLoop(cells, [this, acc, &eqreg, &ptable, &psat]
        (const CellID                 cell,
         Details::PhaseQuantityValue<Scalar>& pressures,
         Details::PhaseQuantityValue<Scalar>& saturations,
         Scalar&                      Rs,
         Scalar&                      Rv,
         Scalar&                      Rvw) -> void
    {
        pressures  .reset();
        saturations.reset();

        Scalar totfrac = 0.0;
        for (const auto& [depth, frac] : Details::horizontalSubdivision(cell, this->cellZSpan_[cell], acc)) {
            const auto pos = CellPos { cell, depth };

            saturations.axpy(psat.deriveSaturations(pos, eqreg, ptable), frac);
            pressures  .axpy(psat.correctedPhasePressures(), frac);

            totfrac += frac;
        }

        if (totfrac > 0.) {
            saturations /= totfrac;
            pressures /= totfrac;
        } else {
            // Fall back to centre point method for zero-thickness cells.
            const auto pos = CellPos {
                    cell, this->cellCenterDepth_[cell]
            };

            saturations = psat.deriveSaturations(pos, eqreg, ptable);
            pressures   = psat.correctedPhasePressures();
        }

        const auto temp = this->temperature_[cell];
        const auto cz   = this->cellCenterDepth_[cell];

        Rs = eqreg.dissolutionCalculator()
            (cz, pressures.oil, temp, saturations.gas);

        Rv = eqreg.evaporationCalculator()
            (cz, pressures.gas, temp, saturations.oil);

        Rvw = eqreg.waterEvaporationCalculator()
            (cz, pressures.gas, temp, saturations.water);
    });
}

template<class FluidSystem, class Grid, class GridView, class ElementMapper, class CartesianIndexMapper>
template<class CellRange, class PressTable, class PhaseSat>
void InitialStateComputer<FluidSystem, Grid, GridView, ElementMapper, CartesianIndexMapper>::
equilibrateTiltedFaultBlockSimple(const CellRange& cells,
                             const EquilReg<Scalar>& eqreg,
                             const GridView&         gridView,
                             const int               acc,
                             const PressTable&       ptable,
                             PhaseSat&               psat)
{
    using CellPos = typename PhaseSat::Position;
    using CellID  = std::remove_cv_t<std::remove_reference_t<
        decltype(std::declval<CellPos>().cell)>>;

    this->cellLoop(cells, [this, acc, &eqreg, &ptable, &psat, &gridView]
        (const CellID                 cell,
         Details::PhaseQuantityValue<Scalar>& pressures,
         Details::PhaseQuantityValue<Scalar>& saturations,
         Scalar&                      Rs,
         Scalar&                      Rv,
         Scalar&                      Rvw) -> void
    {
        pressures.reset();
        saturations.reset();
        Scalar totalWeight = 0.0;

        // We assume grid blocks are treated as being tilted
        const auto& [zmin, zmax] = this->cellZMinMax_[cell];
        const Scalar cellThickness = zmax - zmin;
        const Scalar halfThickness = cellThickness / 2.0;

        // Calculate dip parameters from corner point geometry
        Scalar dipAngle, dipAzimuth;
        Details::computeBlockDip(this->cellCorners_[cell], dipAngle, dipAzimuth);

        // Reference point for TVD calculations
                std::array<Scalar, 3> referencePoint = {
            this->cellCenterXY_[cell].first,
            this->cellCenterXY_[cell].second,
            this->cellCenterDepth_[cell]
        };

        //  We have acc levels within each half (upper and lower) of the block
        const int numLevelsPerHalf = std::min(20, acc);

        // Create subdivisions for upper and lower halves with cross-section weighting
        std::vector<std::pair<Scalar, Scalar>> levels;

        // Subdivide upper and lower halves separately
        for (int side = 0; side < 2; ++side) {
            Scalar halfStart = (side == 0) ? zmin : zmin + halfThickness;

            for (int i = 0; i < numLevelsPerHalf; ++i) {
                // Calculate depth at the center of this subdivision
                Scalar depth = halfStart + (i + 0.5) * (halfThickness / numLevelsPerHalf);

                // A simple way: we can estimate cross-section weight based on dip angle
                // For horizontal cells: weight = 1.0, for tilted cells: weight decreases with dip
                Scalar crossSectionWeight = (halfThickness / numLevelsPerHalf);

                // Apply dip correction to weight (cross-section area decreases with dip)
                if (std::abs(dipAngle) > 1e-10) {
                    crossSectionWeight /= std::cos(dipAngle);
                }

                levels.emplace_back(depth, crossSectionWeight);
            }
        }

        for (const auto& [depth, weight] : levels) {
            // Convert measured depth to True Vertical Depth for tilted blocks
            const auto& [x, y] = this->cellCenterXY_[cell];
            Scalar tvd = Details::calculateTrueVerticalDepth(
                depth, x, y, dipAngle, dipAzimuth, referencePoint);

            const auto pos = CellPos{cell, tvd};

            auto localSaturations = psat.deriveSaturations(pos, eqreg, ptable);
            auto localPressures = psat.correctedPhasePressures();

            // Apply cross-section weighted averaging
            saturations.axpy(localSaturations, weight);
            pressures.axpy(localPressures, weight);
            totalWeight += weight;
        }

        // Normalize results
        if (totalWeight > 1e-10) {
            saturations /= totalWeight;
            pressures /= totalWeight;
        } else {
            // Fallback to center point method using TVD
            const auto& [x, y] = this->cellCenterXY_[cell];
            Scalar tvdCenter = Details::calculateTrueVerticalDepth(
                this->cellCenterDepth_[cell], x, y, dipAngle, dipAzimuth, referencePoint);
            const auto pos = CellPos{cell, tvdCenter};
            saturations = psat.deriveSaturations(pos, eqreg, ptable);
            pressures = psat.correctedPhasePressures();
        }

        // Compute solution ratios at cell center TVD
        const auto temp = this->temperature_[cell];
        const auto& [x, y] = this->cellCenterXY_[cell];
        Scalar tvdCenter = Details::calculateTrueVerticalDepth(
            this->cellCenterDepth_[cell], x, y, dipAngle, dipAzimuth, referencePoint);

        Rs = eqreg.dissolutionCalculator()(tvdCenter, pressures.oil, temp, saturations.gas);
        Rv = eqreg.evaporationCalculator()(tvdCenter, pressures.gas, temp, saturations.oil);
        Rvw = eqreg.waterEvaporationCalculator()(tvdCenter, pressures.gas, temp, saturations.water);
    });
}

template<class FluidSystem, class Grid, class GridView, class ElementMapper, class CartesianIndexMapper>
template<class CellRange, class PressTable, class PhaseSat>
void InitialStateComputer<FluidSystem, Grid, GridView, ElementMapper, CartesianIndexMapper>::
equilibrateTiltedFaultBlock(const CellRange&        cells,
                             const EquilReg<Scalar>& eqreg,
                             const GridView&         gridView,
                             const int               acc,
                             const PressTable&       ptable,
                             PhaseSat&               psat)
{
    using CellPos = typename PhaseSat::Position;
    using CellID  = std::remove_cv_t<std::remove_reference_t<
        decltype(std::declval<CellPos>().cell)>>;

    std::vector<typename GridView::template Codim<0>::Entity> entityMap(gridView.size(0));
    for (const auto& entity : entities(gridView, Dune::Codim<0>())) {
        CellID idx = gridView.indexSet().index(entity);
        entityMap[idx] = entity;
    }

    // Face Area Calculation
    auto polygonArea = [](const std::vector<std::array<Scalar, 2>>& pts) {
        if (pts.size() < 3) return Scalar(0);
        Scalar area = 0;
        for (size_t i = 0; i < pts.size(); ++i) {
            size_t j = (i + 1) % pts.size();
            area += pts[i][0] * pts[j][1] - pts[j][0] * pts[i][1];
        }
        return std::abs(area) * Scalar(0.5);
    };

    // Compute horizontal cross-section at given depth
    auto computeCrossSectionArea = [&](const CellID cell, Scalar depth) -> Scalar {
        try {
            const auto& entity = entityMap[cell];
            const auto& geometry = entity.geometry();
            const int numCorners = geometry.corners();

            std::vector<std::array<Scalar, 3>> corners(numCorners);
            for (int i = 0; i < numCorners; ++i) {
                const auto& corner = geometry.corner(i);
                corners[i] = {static_cast<Scalar>(corner[0]), static_cast<Scalar>(corner[1]), static_cast<Scalar>(corner[2])};
            }

            // Find all intersections between horizontal plane and cell edges
            std::vector<std::array<Scalar, 2>> intersectionPoints;
            const Scalar tol = 1e-10;

            // Check all edges between corners (could be optimized further)
            for (size_t i = 0; i < corners.size(); ++i) {
                for (size_t j = i + 1; j < corners.size(); ++j) {
                    Scalar za = corners[i][2];
                    Scalar zb = corners[j][2];

                    if ((za - depth) * (zb - depth) <= 0.0 && std::abs(za - zb) > tol) {
                        // Edge crosses the horizontal plane
                        Scalar t = (depth - za) / (zb - za);
                        Scalar x = corners[i][0] + t * (corners[j][0] - corners[i][0]);
                        Scalar y = corners[i][1] + t * (corners[j][1] - corners[i][1]);
                        intersectionPoints.push_back({x, y});
                    }
                }
            }

            // Remove duplicates
            if (intersectionPoints.size() > 1) {
                auto pointsEqual = [tol](const std::array<Scalar, 2>& a, const std::array<Scalar, 2>& b) {
                    return std::abs(a[0] - b[0]) < tol && std::abs(a[1] - b[1]) < tol;
                };

                intersectionPoints.erase(
                    std::unique(intersectionPoints.begin(), intersectionPoints.end(), pointsEqual),
                    intersectionPoints.end()
                );
            }

            if (intersectionPoints.size() < 3) {
                // No valid grid found, use fallback
                return 0.0;
            }

            // Order points counter-clockwise around centroid
            Scalar cx = 0, cy = 0;
            for (const auto& p : intersectionPoints) {
                cx += p[0]; cy += p[1];
            }
            cx /= intersectionPoints.size();
            cy /= intersectionPoints.size();

            // Sorting
            auto angleCompare = [cx, cy](const std::array<Scalar, 2>& a, const std::array<Scalar, 2>& b) {
                return std::atan2(a[1] - cy, a[0] - cx) < std::atan2(b[1] - cy, b[0] - cx);
            };

            std::ranges::sort(intersectionPoints, angleCompare);

            return polygonArea(intersectionPoints);

        } catch (const std::exception& e) {
            return 0.0;
        }
    };

    auto cellProcessor = [this, acc, &eqreg, &ptable, &psat, &computeCrossSectionArea]
        (const CellID                 cell,
         Details::PhaseQuantityValue<Scalar>& pressures,
         Details::PhaseQuantityValue<Scalar>& saturations,
         Scalar&                      Rs,
         Scalar&                      Rv,
         Scalar&                      Rvw) -> void
    {
        pressures.reset();
        saturations.reset();
        Scalar totalWeight = 0.0;

        const auto& zmin = this->cellZMinMax_[cell].first;
        const auto& zmax = this->cellZMinMax_[cell].second;
        const Scalar cellThickness = zmax - zmin;
        const Scalar halfThickness = cellThickness / 2.0;

        // Calculate dip parameters from corner point geometry
        Scalar dipAngle, dipAzimuth;
        Details::computeBlockDip(this->cellCorners_[cell], dipAngle, dipAzimuth);

        // Reference point for TVD calculations
        std::array<Scalar, 3> referencePoint = {
            this->cellCenterXY_[cell].first,
            this->cellCenterXY_[cell].second,
            this->cellCenterDepth_[cell]
        };

        // We have acc levels within each half (upper and lower) of the block
        const int numLevelsPerHalf = std::min(20, acc);

        // Create subdivisions for upper and lower halves with cross-section weighting
        std::vector<std::pair<Scalar, Scalar>> levels;

        // Subdivide upper and lower halves separately
        for (int side = 0; side < 2; ++side) {
            Scalar halfStart = (side == 0) ? zmin : zmin + halfThickness;

            for (int i = 0; i < numLevelsPerHalf; ++i) {
                // Calculate depth at the center of this subdivision
                Scalar depth = halfStart + (i + 0.5) * (halfThickness / numLevelsPerHalf);

                // Compute cross-section area at this depth
                Scalar crossSectionArea = computeCrossSectionArea(cell, depth);

                // Weight is proportional to: area × Δz (volume element)
                Scalar weight = crossSectionArea * (halfThickness / numLevelsPerHalf);

                levels.emplace_back(depth, weight);
            }
        }

        // Maybe not necessary (for debug)
        bool hasValidAreas = false;
        for (const auto& level : levels) {
            if (level.second > 1e-10) {
                hasValidAreas = true;
                break;
            }
        }

        if (!hasValidAreas) {
            // Fallback to dip-based weighting as used in equilibrateTiltedFaultBlockSimple
            levels.clear();
            for (int side = 0; side < 2; ++side) {
                Scalar halfStart = (side == 0) ? zmin : zmin + halfThickness;
                for (int i = 0; i < numLevelsPerHalf; ++i) {
                    Scalar depth = halfStart + (i + 0.5) * (halfThickness / numLevelsPerHalf);
                    Scalar weight = (halfThickness / numLevelsPerHalf);
                    if (std::abs(dipAngle) > 1e-10) {
                        weight /= std::cos(dipAngle);
                    }
                    levels.emplace_back(depth, weight);
                }
            }
        }

        for (const auto& level : levels) {
            Scalar depth = level.first;
            Scalar weight = level.second;

            // Convert measured depth to True Vertical Depth for tilted blocks
            const auto& xy = this->cellCenterXY_[cell];
            Scalar tvd = Details::calculateTrueVerticalDepth(
                depth, xy.first, xy.second, dipAngle, dipAzimuth, referencePoint);

            const auto pos = CellPos{cell, tvd};

            auto localSaturations = psat.deriveSaturations(pos, eqreg, ptable);
            auto localPressures = psat.correctedPhasePressures();

            // Apply cross-section weighted averaging
            saturations.axpy(localSaturations, weight);
            pressures.axpy(localPressures, weight);
            totalWeight += weight;
        }

        if (totalWeight > 1e-10) {
            saturations /= totalWeight;
            pressures /= totalWeight;
        } else {
            // Fallback to center point method using TVD
            const auto& xy = this->cellCenterXY_[cell];
            Scalar tvdCenter = Details::calculateTrueVerticalDepth(
                this->cellCenterDepth_[cell], xy.first, xy.second, dipAngle, dipAzimuth, referencePoint);
            const auto pos = CellPos{cell, tvdCenter};
            saturations = psat.deriveSaturations(pos, eqreg, ptable);
            pressures = psat.correctedPhasePressures();
        }

        // Compute solution ratios at cell center TVD
        const auto temp = this->temperature_[cell];
        const auto& xy = this->cellCenterXY_[cell];
        Scalar tvdCenter = Details::calculateTrueVerticalDepth(
            this->cellCenterDepth_[cell], xy.first, xy.second, dipAngle, dipAzimuth, referencePoint);

        Rs = eqreg.dissolutionCalculator()(tvdCenter, pressures.oil, temp, saturations.gas);
        Rv = eqreg.evaporationCalculator()(tvdCenter, pressures.gas, temp, saturations.oil);
        Rvw = eqreg.waterEvaporationCalculator()(tvdCenter, pressures.gas, temp, saturations.water);
    };

    this->cellLoop(cells, cellProcessor);
}
}
} // namespace EQUIL
} // namespace Opm

#endif // OPM_INIT_STATE_EQUIL_IMPL_HPP
